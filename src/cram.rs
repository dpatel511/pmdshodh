//! CRAM input decoded container by container on the thread pool (docs/design.md §11).
//!
//! noodles' record reader decodes one container at a time on the reading thread. Containers are
//! independent, so the reading thread reads several raw containers (cheap) and the pool decodes
//! them at once; records keep their input order.

use std::io::{self, Read};
use std::ops::Range;
use std::sync::Arc;

use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use rayon::prelude::*;

pub struct Reader {
    inner: cram::io::Reader<Box<dyn Read + Send>>,
    /// The end-of-file container was read; nothing after it is.
    eof: bool,
}

impl Reader {
    pub fn new(inner: Box<dyn Read + Send>) -> Self {
        Reader {
            inner: cram::io::Reader::new(inner),
            eof: false,
        }
    }

    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_header()
    }

    /// Read up to `n` containers; fewer (none at all) at the end of the input.
    pub fn read_containers(&mut self, n: usize) -> io::Result<Vec<cram::io::reader::Container>> {
        let mut containers = Vec::with_capacity(n);
        while !self.eof && containers.len() < n {
            let bytes = read_container_bytes(self.inner.get_mut())?;
            let mut container = cram::io::reader::Container::default();
            if cram::io::Reader::new(&bytes[..]).read_container(&mut container)? == 0 {
                self.eof = true;
            } else {
                containers.push(container);
            }
        }
        Ok(containers)
    }
}

/// One container's bytes, header and body, for noodles to parse. noodles' `read_container`
/// allocates the body length the header claims before reading it, so a corrupt header could cost
/// up to 2 GB (audit finding B10); here the body is read as it arrives, and a container longer than
/// the input is an error that costs only the bytes the input has.
fn read_container_bytes(input: &mut impl Read) -> io::Result<Vec<u8>> {
    // CRAM 3 §7.3: length, reference sequence ID, alignment start, span, records (ITF8), record
    // counter, bases (LTF8), blocks, landmarks (ITF8 array), CRC-32.
    let mut head = Recorded {
        inner: &mut *input,
        bytes: Vec::new(),
    };
    let mut len = [0; 4];
    head.read_exact(&mut len)?;
    for _ in 0..4 {
        read_itf8(&mut head)?;
    }
    for _ in 0..2 {
        read_ltf8(&mut head)?;
    }
    read_itf8(&mut head)?;
    for _ in 0..read_itf8(&mut head)?.max(0) {
        read_itf8(&mut head)?;
    }
    head.read_exact(&mut [0; 4])?;
    let mut bytes = head.bytes;
    // A negative length is left for noodles to reject.
    let len = u64::try_from(i32::from_le_bytes(len)).unwrap_or(0);
    let read = input.take(len).read_to_end(&mut bytes)?;
    if (read as u64) < len {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("container of {len} bytes is cut short: the input ends after {read}"),
        ));
    }
    Ok(bytes)
}

/// A reader that keeps a copy of the bytes read through it.
struct Recorded<'a, R> {
    inner: &'a mut R,
    bytes: Vec<u8>,
}

impl<R: Read> Read for Recorded<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

fn read_u8(reader: &mut impl Read) -> io::Result<u8> {
    let mut b = [0];
    reader.read_exact(&mut b)?;
    Ok(b[0])
}

/// An ITF8 integer (CRAM 3 §2.3): the leading 1 bits of the first byte count the bytes after it,
/// at most 4; the fifth byte gives only its low 4 bits.
fn read_itf8(reader: &mut impl Read) -> io::Result<i32> {
    let b0 = read_u8(reader)?;
    let extra = b0.leading_ones().min(4);
    let mut n = if extra == 4 {
        u32::from(b0 & 0x0f)
    } else {
        u32::from(b0 & (0x7f >> extra))
    };
    for i in 0..extra {
        let b = u32::from(read_u8(reader)?);
        n = if i == 3 {
            n << 4 | b & 0x0f
        } else {
            n << 8 | b
        };
    }
    Ok(n as i32)
}

/// An LTF8 integer (CRAM 3 §2.3): as ITF8, up to 8 bytes after the first. Only its length matters.
fn read_ltf8(reader: &mut impl Read) -> io::Result<()> {
    for _ in 0..read_u8(reader)?.leading_ones() {
        read_u8(reader)?;
    }
    Ok(())
}

/// The reference sequence a container's records are all aligned to, if it has exactly one.
///
/// noodles 0.119 keeps the container header's reference sequence context type private, so the ID
/// is read from its `Debug` output (`Some(Context { reference_sequence_id: 3, .. })`). noodles is
/// pinned, and a test checks this; if the format changes this returns `None`, which only loses the
/// up-front reference loading in `Reference::cram_batch`, not correctness.
pub fn reference_sequence_id(container: &cram::io::reader::Container) -> Option<usize> {
    let context = format!("{:?}", container.header().reference_sequence_context());
    let digits = context.strip_prefix("Some(Context { reference_sequence_id: ")?;
    let end = digits.find(|c: char| !c.is_ascii_digit())?;
    digits[..end].parse().ok()
}

/// Decoded records in input order, with each mapped record's reference window (alignment start,
/// CIGAR reference span; uppercase) back to back in `bases`.
#[derive(Default)]
pub struct Decoded {
    pub records: Vec<RecordBuf>,
    pub bases: Vec<u8>,
    pub windows: Vec<Option<Range<usize>>>,
}

/// Decode containers and cut their records' reference windows on the pool, in input order.
pub fn decode(
    containers: &[cram::io::reader::Container],
    repository: &fasta::Repository,
    header: &sam::Header,
) -> io::Result<Decoded> {
    let parts: Vec<io::Result<Decoded>> = containers
        .par_iter()
        .map(|c| decode_container(c, repository, header))
        .collect();
    let mut all = Decoded::default();
    for part in parts {
        let part = part?;
        let offset = all.bases.len();
        all.bases.extend_from_slice(&part.bases);
        all.windows.extend(
            part.windows
                .into_iter()
                .map(|w| w.map(|r| r.start + offset..r.end + offset)),
        );
        all.records.extend(part.records);
    }
    Ok(all)
}

/// As noodles' `cram::io::BufReader` decodes a container, plus the windows.
fn decode_container(
    container: &cram::io::reader::Container,
    repository: &fasta::Repository,
    header: &sam::Header,
) -> io::Result<Decoded> {
    let compression_header = container.compression_header()?;
    let mut out = Decoded {
        records: Vec::with_capacity(container.header().record_count()),
        ..Decoded::default()
    };
    let mut sequence: Option<(usize, Arc<fasta::record::Sequence>)> = None;
    for slice in container.slices() {
        let slice = slice?;
        let (core, external) = slice.decode_blocks()?;
        for record in slice.records(
            repository.clone(),
            header,
            &compression_header,
            &core,
            &external,
        )? {
            // The window comes first: it checks that a mapped read ends within its reference
            // sequence, and converting one that does not panics in noodles (audit finding B9).
            let window = window(&record, header, repository, &mut sequence, &mut out.bases)?;
            let record = RecordBuf::try_from_alignment_record(header, &record)?;
            out.windows.push(window);
            out.records.push(record);
        }
    }
    Ok(out)
}

/// Append a mapped record's reference bases to `bases` and return their range (`None` if
/// unmapped). `sequence` keeps the last sequence used.
fn window(
    record: &impl sam::alignment::Record,
    header: &sam::Header,
    repository: &fasta::Repository,
    sequence: &mut Option<(usize, Arc<fasta::record::Sequence>)>,
    bases: &mut Vec<u8>,
) -> io::Result<Option<Range<usize>>> {
    use sam::alignment::Record;
    if Record::flags(record)?.is_unmapped() {
        return Ok(None);
    }
    let invalid = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg);
    let id = Record::reference_sequence_id(record, header)
        .transpose()?
        .ok_or_else(|| invalid("mapped read has no reference sequence".into()))?;
    let start = Record::alignment_start(record)
        .transpose()?
        .ok_or_else(|| invalid("mapped read has no position".into()))?;
    let span = Record::alignment_span(record).transpose()?.unwrap_or(0);
    let (name, _) = header
        .reference_sequences()
        .get_index(id)
        .ok_or_else(|| invalid("reference sequence ID is not in the header".into()))?;
    if sequence.as_ref().is_none_or(|(i, _)| *i != id) {
        let s = repository
            .get(name)
            .transpose()?
            .ok_or_else(|| invalid(format!("reference has no sequence named {name}")))?;
        *sequence = Some((id, s));
    }
    let bases_of = sequence.as_ref().unwrap().1.as_ref().as_ref();
    let (from, to) = (usize::from(start) - 1, usize::from(start) - 1 + span);
    if to > bases_of.len() {
        let read = Record::name(record).map_or_else(|| "*".into(), |n| n.to_string());
        return Err(invalid(format!(
            "read {read}: alignment ends at {to} but reference sequence {name} has {} bases",
            bases_of.len()
        )));
    }
    let at = bases.len();
    bases.extend_from_slice(&bases_of[from..to]);
    Ok(Some(at..bases.len()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::io::Write as _;

    /// A CRAM written by noodles with one container of reads on the second of two sequences.
    fn cram() -> (Vec<u8>, sam::Header, Vec<RecordBuf>, fasta::Repository) {
        let c1 = b"ACGTACGTACGTACGTACGT".to_vec();
        let c2 = b"TTGCAACGTTGCAACGTTGCAACGT".to_vec();
        let def = |n: &str| fasta::record::Definition::new(n, None);
        let repository = fasta::Repository::new(vec![
            fasta::Record::new(def("c1"), c1.into()),
            fasta::Record::new(def("c2"), c2.into()),
        ]);
        let text = "@SQ\tSN:c1\tLN:20\n@SQ\tSN:c2\tLN:25\n\
                    r1\t0\tc2\t2\t60\t10M\t*\t0\t0\tTGCAACGTTG\tIIIIIIIIII\n\
                    r2\t16\tc2\t5\t60\t4M1D6M\t*\t0\t0\tAACGTGCAAC\tIIIIIIIIII\n";
        let mut r = sam::io::Reader::new(text.as_bytes());
        let header = r.read_header().unwrap();
        let records: Vec<RecordBuf> = r.record_bufs(&header).map(Result::unwrap).collect();
        let mut w = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        w.write_header(&header).unwrap();
        for record in &records {
            w.write_alignment_record(&header, record).unwrap();
        }
        w.try_finish(&header).unwrap();
        (w.into_inner(), header, records, repository)
    }

    #[test]
    fn containers_decode_in_order_and_name_their_reference() {
        let (bytes, header, records, repository) = cram();
        let mut reader = Reader::new(Box::new(io::Cursor::new(bytes)));
        let read = reader.read_header().unwrap();
        // The writer adds M5 checksums; the sequences are the same.
        assert!(
            read.reference_sequences()
                .keys()
                .eq(header.reference_sequences().keys())
        );
        let containers = reader.read_containers(4).unwrap();
        assert_eq!(containers.len(), 1);
        assert_eq!(reference_sequence_id(&containers[0]), Some(1));
        let decoded = decode(&containers, &repository, &header).unwrap();
        let names = |r: &[RecordBuf]| {
            r.iter()
                .map(|r| r.name().map(|n| n.to_string()))
                .collect::<Vec<_>>()
        };
        assert_eq!(names(&decoded.records), names(&records));
        assert_eq!(decoded.records[1].sequence(), records[1].sequence());
        // Windows: c2 from position 2 (10M, 10 bases) and from position 5 (4M1D6M, 11 bases).
        let window = |i: usize| &decoded.bases[decoded.windows[i].clone().unwrap()];
        assert_eq!(window(0), b"TGCAACGTTG");
        assert_eq!(window(1), b"AACGTTGCAAC");
        // After the end-of-file container, nothing more is read (a second read would fail).
        assert!(reader.read_containers(4).unwrap().is_empty());
        assert!(reader.read_containers(4).unwrap().is_empty());
    }

    /// CRC-32 (IEEE), as CRAM 3 container headers carry.
    fn crc32(bytes: &[u8]) -> u32 {
        let mut crc = !0u32;
        for &b in bytes {
            crc ^= u32::from(b);
            for _ in 0..8 {
                crc = if crc & 1 == 1 {
                    (crc >> 1) ^ 0xedb8_8320
                } else {
                    crc >> 1
                };
            }
        }
        !crc
    }

    /// A container header with a correct checksum claims 2 GB, and 100 bytes follow it: an error,
    /// read as the bytes arrive. noodles allocates the claimed length before reading it, so a
    /// corrupt input of a few kilobytes cost 1.6 GB (audit finding B10).
    #[test]
    fn container_longer_than_the_input_is_cut_short() {
        let (_, header, _, repository) = cram();
        let mut w = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_writer(Vec::new());
        w.write_header(&header).unwrap();
        let mut bytes = w.into_inner();
        // Length, then reference sequence 0, start 1, span 10, 1 record, record counter 0,
        // 10 bases, 1 block, no landmarks (ITF8 and LTF8 values below 128 are one byte).
        let mut container = i32::MAX.to_le_bytes().to_vec();
        container.extend([0, 1, 10, 1, 0, 10, 1, 0]);
        let crc = crc32(&container);
        container.extend(crc.to_le_bytes());
        bytes.extend(container);
        bytes.extend([0; 100]);
        let mut reader = Reader::new(Box::new(io::Cursor::new(bytes)));
        reader.read_header().unwrap();
        let err = reader.read_containers(4).err().unwrap();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(
            err.to_string(),
            "container of 2147483647 bytes is cut short: the input ends after 100"
        );
    }
}
