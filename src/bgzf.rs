//! BGZF reading and writing with (de)compression on the global rayon pool.
//!
//! Scoring runs on the same pool, so `-@ N` bounds all of it together (main.rs `init_threads`).
//! Framing and compression follow noodles-bgzf 0.52 (zlib-rs, level 6, blocks of at most
//! `MAX_BUF_SIZE` bytes), so the bytes written are the same as before and do not depend on the
//! thread count. Blocks are handed out and written in input order.

use std::collections::VecDeque;
use std::io::{self, Read, Write};
use std::sync::mpsc;

const HEADER_SIZE: usize = 18;
const TRAILER_SIZE: usize = 8;
const MAX_ISIZE: usize = 1 << 16;
/// Worst-case DEFLATE overhead at level 0 (as in noodles-bgzf).
const LEVEL_0_OVERHEAD: usize = 15;
/// Largest uncompressed block; its compressed frame always fits in 64 KiB.
const MAX_BUF_SIZE: usize = MAX_ISIZE - HEADER_SIZE - TRAILER_SIZE - LEVEL_0_OVERHEAD;
const LEVEL: i32 = 6;
const WINDOW_BITS: u8 = 15;
/// SAM spec §4.1.2 end-of-file marker.
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

/// A block being (de)compressed; its result arrives through the receiver.
type Pending = mpsc::Receiver<io::Result<Vec<u8>>>;

/// Run `f` on the pool when `parallel`, otherwise right away.
fn submit(parallel: bool, f: impl FnOnce() -> io::Result<Vec<u8>> + Send + 'static) -> Pending {
    let (tx, rx) = mpsc::sync_channel(1);
    if parallel {
        rayon::spawn(move || {
            let _ = tx.send(f());
        });
    } else {
        let _ = tx.send(f());
    }
    rx
}

fn wait(rx: Pending) -> io::Result<Vec<u8>> {
    wait_for(&rx).unwrap_or_else(|| Err(io::Error::other("BGZF worker stopped")))
}

/// Receive from `rx`; `None` if its sender is gone. On a pool thread, run other pool jobs while
/// waiting instead of idling, so a thread that waits still counts as one doing work.
pub fn wait_for<T>(rx: &mpsc::Receiver<T>) -> Option<T> {
    loop {
        match rx.try_recv() {
            Ok(v) => return Some(v),
            Err(mpsc::TryRecvError::Disconnected) => return None,
            Err(mpsc::TryRecvError::Empty) => match rayon::yield_now() {
                Some(rayon::Yield::Executed) => {}
                // Nothing queued: wait briefly, then look for new jobs again.
                Some(rayon::Yield::Idle) => {
                    match rx.recv_timeout(std::time::Duration::from_micros(100)) {
                        Ok(v) => return Some(v),
                        Err(mpsc::RecvTimeoutError::Disconnected) => return None,
                        Err(mpsc::RecvTimeoutError::Timeout) => {}
                    }
                }
                None => return rx.recv().ok(),
            },
        }
    }
}

/// Blocks kept in flight: enough to keep every pool thread busy, and bounded.
fn depth(parallel: bool) -> usize {
    if parallel {
        2 * rayon::current_num_threads() + 2
    } else {
        1
    }
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

pub struct Reader<R> {
    inner: R,
    parallel: bool,
    pending: VecDeque<Pending>,
    block: Vec<u8>,
    pos: usize,
    done: bool,
}

impl<R: Read> Reader<R> {
    pub fn new(inner: R, parallel: bool) -> Self {
        Reader {
            inner,
            parallel,
            pending: VecDeque::new(),
            block: Vec::new(),
            pos: 0,
            done: false,
        }
    }

    fn fill(&mut self) {
        while !self.done && self.pending.len() < depth(self.parallel) {
            match read_frame(&mut self.inner) {
                Ok(Some(frame)) => self
                    .pending
                    .push_back(submit(self.parallel, move || inflate(&frame))),
                Ok(None) => self.done = true,
                // Queued behind the blocks before it, so earlier data is read first.
                Err(e) => {
                    self.pending.push_back(submit(false, || Err(e)));
                    self.done = true;
                }
            }
        }
    }
}

impl<R: Read> Read for Reader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos == self.block.len() {
            self.fill();
            let Some(rx) = self.pending.pop_front() else {
                return Ok(0);
            };
            self.block = wait(rx)?;
            self.pos = 0;
        }
        let n = buf.len().min(self.block.len() - self.pos);
        buf[..n].copy_from_slice(&self.block[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Read one whole BGZF frame; `None` at the end of the input.
fn read_frame(r: &mut impl Read) -> io::Result<Option<Vec<u8>>> {
    let mut frame = vec![0; HEADER_SIZE];
    let mut got = 0;
    while got < HEADER_SIZE {
        match r.read(&mut frame[got..]) {
            Ok(0) if got == 0 => return Ok(None),
            Ok(0) => return Err(invalid("truncated BGZF block")),
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    let valid = frame[0..4] == [0x1f, 0x8b, 0x08, 0x04]
        && frame[10..12] == [6, 0]
        && frame[12..14] == *b"BC"
        && frame[14..16] == [2, 0];
    if !valid {
        return Err(invalid("invalid BGZF header"));
    }
    let size = usize::from(u16::from_le_bytes([frame[16], frame[17]])) + 1;
    if size < HEADER_SIZE + TRAILER_SIZE {
        return Err(invalid("invalid BGZF block size"));
    }
    frame.resize(size, 0);
    r.read_exact(&mut frame[HEADER_SIZE..]).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            invalid("truncated BGZF block")
        } else {
            e
        }
    })?;
    Ok(Some(frame))
}

fn inflate(frame: &[u8]) -> io::Result<Vec<u8>> {
    let (body, trailer) = frame.split_at(frame.len() - TRAILER_SIZE);
    let crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
    let isize = u32::from_le_bytes(trailer[4..].try_into().unwrap()) as usize;
    if isize > MAX_ISIZE {
        return Err(invalid("invalid BGZF ISIZE"));
    }
    let mut data = vec![0; isize];
    let status = zlib_rs::Inflate::new(false, WINDOW_BITS)
        .decompress(
            &body[HEADER_SIZE..],
            &mut data,
            zlib_rs::InflateFlush::Finish,
        )
        .map_err(|_| invalid("invalid BGZF block data"))?;
    if status != zlib_rs::Status::StreamEnd {
        return Err(invalid("invalid BGZF block data"));
    }
    if zlib_rs::crc32::crc32(0, &data) != crc {
        return Err(invalid("block data checksum mismatch"));
    }
    Ok(data)
}

pub struct Writer<W: Write> {
    inner: W,
    parallel: bool,
    buf: Vec<u8>,
    pending: VecDeque<Pending>,
}

impl<W: Write> Writer<W> {
    pub fn new(inner: W, parallel: bool) -> Self {
        Writer {
            inner,
            parallel,
            buf: Vec::with_capacity(MAX_BUF_SIZE),
            pending: VecDeque::new(),
        }
    }

    /// Write the finished blocks at the front; wait for the oldest while too many are in flight.
    fn drain(&mut self, keep: usize) -> io::Result<()> {
        while let Some(rx) = self.pending.front() {
            let frame = if self.pending.len() > keep {
                wait(self.pending.pop_front().unwrap())?
            } else {
                match rx.try_recv() {
                    Ok(result) => {
                        self.pending.pop_front();
                        result?
                    }
                    Err(mpsc::TryRecvError::Empty) => return Ok(()),
                    Err(mpsc::TryRecvError::Disconnected) => {
                        return Err(io::Error::other("BGZF worker stopped"));
                    }
                }
            };
            self.inner.write_all(&frame)?;
        }
        Ok(())
    }

    /// Write every remaining block and the end-of-file marker.
    pub fn finish(&mut self) -> io::Result<()> {
        self.flush()?;
        self.drain(0)?;
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = buf.len().min(MAX_BUF_SIZE - self.buf.len());
        self.buf.extend_from_slice(&buf[..n]);
        if self.buf.len() == MAX_BUF_SIZE {
            self.flush()?;
        }
        Ok(n)
    }

    /// Ends the current block, as noodles-bgzf does; the block is written once compressed.
    fn flush(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            let data = std::mem::replace(&mut self.buf, Vec::with_capacity(MAX_BUF_SIZE));
            self.pending
                .push_back(submit(self.parallel, move || deflate(&data)));
        }
        self.drain(depth(self.parallel))
    }
}

/// One BGZF frame holding `data` (at most `MAX_BUF_SIZE` bytes).
fn deflate(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut cdata = vec![0; zlib_rs::compress_bound(data.len())];
    let mut len = None;
    for level in [LEVEL, 0] {
        let mut encoder = zlib_rs::Deflate::new(level, false, WINDOW_BITS);
        let status = encoder
            .compress(data, &mut cdata, zlib_rs::DeflateFlush::Finish)
            .map_err(|_| io::Error::other("BGZF compression failed"))?;
        if status != zlib_rs::Status::StreamEnd {
            return Err(io::Error::other("BGZF compression failed"));
        }
        let n = encoder.total_out() as usize;
        if n <= MAX_BUF_SIZE + LEVEL_0_OVERHEAD {
            len = Some(n);
            break;
        }
    }
    let len = len.ok_or_else(|| io::Error::other("BGZF block does not fit"))?;
    let size = HEADER_SIZE + len + TRAILER_SIZE;
    let mut frame = Vec::with_capacity(size);
    frame.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ]);
    frame.extend_from_slice(&((size - 1) as u16).to_le_bytes());
    frame.extend_from_slice(&cdata[..len]);
    frame.extend_from_slice(&zlib_rs::crc32::crc32(0, data).to_le_bytes());
    frame.extend_from_slice(&(data.len() as u32).to_le_bytes());
    Ok(frame)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Same bytes as noodles-bgzf's writer, which wrote pmdshodh's BAMs before this module.
    #[test]
    fn output_matches_noodles_bgzf_for_any_thread_setting() {
        let data: Vec<u8> = (0..300_000u32)
            .flat_map(|i| (i * 7 % 251).to_le_bytes())
            .collect();
        let mut expected = noodles::bgzf::io::Writer::new(Vec::new());
        expected.write_all(&data[..1000]).unwrap();
        expected.flush().unwrap();
        expected.write_all(&data[1000..]).unwrap();
        let expected = expected.finish().unwrap();
        for parallel in [false, true] {
            let mut w = Writer::new(Vec::new(), parallel);
            w.write_all(&data[..1000]).unwrap();
            w.flush().unwrap();
            w.write_all(&data[1000..]).unwrap();
            w.finish().unwrap();
            assert!(w.inner == expected, "parallel = {parallel}");
            let mut back = Vec::new();
            Reader::new(&w.inner[..], parallel)
                .read_to_end(&mut back)
                .unwrap();
            assert!(back == data, "round trip, parallel = {parallel}");
        }
    }

    #[test]
    fn truncated_and_corrupt_input_is_an_error() {
        let mut w = Writer::new(Vec::new(), false);
        w.write_all(b"pmdshodh").unwrap();
        w.finish().unwrap();
        let good = w.inner;
        let read = |bytes: &[u8]| {
            let mut out = Vec::new();
            Reader::new(bytes, false)
                .read_to_end(&mut out)
                .map(|_| out)
                .map_err(|e| e.to_string())
        };
        assert_eq!(read(&good).unwrap(), b"pmdshodh");
        assert_eq!(read(&good[..10]).unwrap_err(), "truncated BGZF block");
        assert_eq!(read(&good[..30]).unwrap_err(), "truncated BGZF block");
        let mut bad = good.clone();
        bad[good.len() - 28 - 8] ^= 1; // CRC of the data block
        assert_eq!(read(&bad).unwrap_err(), "block data checksum mismatch");
        let mut bad = good.clone();
        bad[12] = b'X';
        assert_eq!(read(&bad).unwrap_err(), "invalid BGZF header");
    }
}
