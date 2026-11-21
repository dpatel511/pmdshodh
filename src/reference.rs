//! Reference FASTA access for `-r` (docs/design.md §11).
//!
//! Sequences are read through the FASTA index (`.fai`, plus `.gzi` for bgzip) in fixed-size
//! blocks kept in a small least-recently-used cache. Memory stays bounded for sorted and unsorted
//! input, and no whole chromosome or genome is loaded.
//!
//! CRAM decoding is the exception: noodles needs each reference sequence whole, so
//! [`Reference::cram_batch`] loads whole sequences and clears them as the input moves on.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use noodles::core::Region;
use noodles::fasta;

/// Bases per cached block.
const BLOCK: usize = 1 << 20;
/// Blocks kept (64 MiB of sequence).
const CACHED_BLOCKS: usize = 64;

pub struct Reference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    /// Sequence name and length, in index order.
    sequences: Vec<(String, usize)>,
    ids: HashMap<String, usize>,
    blocks: HashMap<(usize, usize), Vec<u8>>,
    /// Cached block keys, least recently used first.
    order: VecDeque<(usize, usize)>,
    /// Set when decoding CRAM.
    cram: Option<CramCache>,
}

/// Whole sequences cached for CRAM decoding before the cache is cleared at the next change of
/// reference sequence (256 Mi bases; the longest human chromosome has 249 Mi).
const CRAM_CACHED_BASES: usize = 1 << 28;

struct CramCache {
    repository: fasta::Repository,
    /// Bases loaded since the cache was last cleared.
    loaded: Arc<AtomicUsize>,
    last: Option<usize>,
    /// Loads sequences for [`Reference::cram_batch`] outside the repository's lock.
    loader: Loader,
    /// A sequence loaded by `loader`, handed to the repository by [`UppercaseAdapter`].
    staged: Staged,
    /// Header IDs of the sequences in the repository.
    cached: HashSet<usize>,
}

type Staged = Arc<Mutex<Option<(Vec<u8>, Vec<u8>)>>>;

/// Reads whole sequences, uppercased as htslib does when it decodes CRAM. Otherwise read bases
/// copied from a soft-masked reference would come out in lowercase.
struct Loader {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    /// The FASTA itself when it is not compressed: sequences are then read in one piece at their
    /// index offsets, which is several times faster than noodles' line-by-line reader.
    plain: Option<File>,
    /// The header's `M5` checksums (lowercase hex) by sequence name. noodles checks a slice
    /// against its reference only when the slice has one reference sequence; multi-reference
    /// slices carry no checksum, so each whole sequence is checked here when it is loaded.
    md5: HashMap<Vec<u8>, String>,
}

impl Loader {
    fn open(path: &Path, header: &noodles::sam::Header) -> io::Result<Self> {
        use noodles::sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
        use std::io::Read;
        let mut plain = File::open(path)?;
        let mut magic = [0u8; 2];
        let compressed = plain.read(&mut magic)? == 2 && magic == [0x1f, 0x8b];
        let md5 = header
            .reference_sequences()
            .iter()
            .filter_map(|(name, seq)| {
                let m5 = seq.other_fields().get(&MD5_CHECKSUM)?;
                Some((name.to_vec(), m5.to_string().to_ascii_lowercase()))
            })
            .collect();
        Ok(Loader {
            reader: fasta::io::indexed_reader::Builder::default().build_from_path(path)?,
            plain: (!compressed).then_some(plain),
            md5,
        })
    }

    /// Read and check a whole sequence (see [`Loader::read`] and [`check_md5`]).
    fn load(&mut self, name: &[u8]) -> io::Result<Vec<u8>> {
        let bases = self.read(name)?;
        if let Some(expected) = self.md5.get(name) {
            check_md5(name, &bases, expected)?;
        }
        Ok(bases)
    }

    /// Read a whole sequence, uppercased, without checking it.
    fn read(&mut self, name: &[u8]) -> io::Result<Vec<u8>> {
        let index: &[fasta::fai::Record] = self.reader.index().as_ref();
        let record = index
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == name)
            .cloned();
        let Some(record) = record else {
            // Not in the index: noodles' error names the sequence.
            return self
                .reader
                .query(&Region::new(name, ..))
                .map(|record| record.sequence().as_ref().to_ascii_uppercase());
        };
        let len = record.length() as usize;
        let mut bases = match &mut self.plain {
            Some(file) => read_checked(file, &record, 0, len)?,
            None => read_checked(self.reader.get_mut(), &record, 0, len)?,
        };
        bases.make_ascii_uppercase();
        Ok(bases)
    }
}

/// Check an uppercased sequence against the CRAM header's `M5` (lowercase hex).
fn check_md5(name: &[u8], bases: &[u8], expected: &str) -> io::Result<()> {
    use md5::Digest;
    let actual: String = md5::Md5::digest(bases)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    if actual != expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sequence {} does not match the CRAM header's M5 checksum \
                 (is -r the reference the CRAM was written with?)",
                String::from_utf8_lossy(name)
            ),
        ));
    }
    Ok(())
}

/// The result of an `M5` check, which may run on the thread pool.
pub type Md5Check = std::sync::mpsc::Receiver<io::Result<()>>;

/// The repository's adapter. noodles holds the repository's write lock while the adapter runs,
/// which stops every decoding thread, so [`Reference::cram_batch`] loads sequences beforehand and
/// stages them here; the adapter loads only what was not staged (containers with several
/// reference sequences).
struct UppercaseAdapter {
    loader: Loader,
    staged: Staged,
    loaded: Arc<AtomicUsize>,
}

impl fasta::repository::Adapter for UppercaseAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let staged = self
            .staged
            .lock()
            .unwrap()
            .take_if(|(n, _)| n.as_slice() == name);
        let bases = match staged {
            Some((_, bases)) => Ok(bases),
            None => self.loader.load(name),
        };
        Some(bases.map(|bases| {
            self.loaded.fetch_add(bases.len(), Ordering::Relaxed);
            fasta::Record::new(fasta::record::Definition::new(name, None), bases.into())
        }))
    }
}

/// Read `len` bases of a sequence from 0-based `start`, at the byte offsets its index gives, with
/// case kept. Every line break must be where the index puts it, and a sequence's first base must
/// follow its header line; otherwise the FASTA has changed since it was indexed (e.g. rewrapped),
/// and the bases at those offsets would be the wrong ones.
fn read_checked<R: io::Read + io::Seek>(
    reader: &mut R,
    record: &fasta::fai::Record,
    start: usize,
    len: usize,
) -> io::Result<Vec<u8>> {
    use std::io::SeekFrom;
    if len == 0 {
        return Ok(Vec::new());
    }
    let (line_bases, line_width) = (
        record.line_base_count().get() as usize,
        record.line_width().get() as usize,
    );
    let mismatch = || {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "sequence {} does not have the line lengths its FASTA index gives (index it again \
                 with samtools faidx)",
                record.name()
            ),
        )
    };
    let offset =
        |i: usize| record.position() + (i / line_bases * line_width + i % line_bases) as u64;
    let after_header = start == 0 && record.position() > 0;
    let from = offset(start) - u64::from(after_header);
    let to = offset(start + len - 1) + 1;
    reader.seek(SeekFrom::Start(from))?;
    let mut raw = vec![0; (to - from) as usize];
    reader.read_exact(&mut raw).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            mismatch()
        } else {
            e
        }
    })?;
    let mut raw = &raw[..];
    if after_header {
        if raw[0] != b'\n' {
            return Err(mismatch());
        }
        raw = &raw[1..];
    }
    // A line at a time: copy its bases, and check that its line ending is all line breaks. If, in
    // addition, the span has no other line breaks (counted in one pass), there are none among the
    // bases.
    let line_break = |b: &u8| *b == b'\n' || *b == b'\r';
    let breaks = raw.iter().filter(|b| line_break(b)).count();
    let mut endings = 0;
    let mut bases = Vec::with_capacity(len);
    let mut column = start % line_bases;
    while !raw.is_empty() {
        let in_bases = column < line_bases;
        let end = if in_bases { line_bases } else { line_width };
        let (part, rest) = raw.split_at((end - column).min(raw.len()));
        if in_bases {
            bases.extend_from_slice(part);
        } else if part.iter().all(line_break) {
            endings += part.len();
        } else {
            return Err(mismatch());
        }
        column += part.len();
        if column == line_width {
            column = 0;
        }
        raw = rest;
    }
    if breaks != endings {
        return Err(mismatch());
    }
    Ok(bases)
}

impl Reference {
    /// Open an indexed FASTA. Fails if the index is missing or unreadable.
    pub fn open(path: &Path) -> io::Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default().build_from_path(path)?;
        let sequences: Vec<(String, usize)> = reader
            .index()
            .as_ref()
            .iter()
            .map(|r| (r.name().to_string(), r.length() as usize))
            .collect();
        let ids = sequences
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), i))
            .collect();
        Ok(Reference {
            reader,
            sequences,
            ids,
            blocks: HashMap::new(),
            order: VecDeque::new(),
            cram: None,
        })
    }

    /// Prepare for decoding CRAM from this FASTA (opened again from `path`). Sequences are
    /// checked against the `M5` checksums in the CRAM's `header` as they are loaded.
    pub fn enable_cram(&mut self, path: &Path, header: &noodles::sam::Header) -> io::Result<()> {
        let loaded = Arc::new(AtomicUsize::new(0));
        let staged: Staged = Arc::default();
        let adapter = UppercaseAdapter {
            loader: Loader::open(path, header)?,
            staged: Arc::clone(&staged),
            loaded: Arc::clone(&loaded),
        };
        self.cram = Some(CramCache {
            repository: fasta::Repository::new(adapter),
            loaded,
            last: None,
            loader: Loader::open(path, header)?,
            staged,
            cached: HashSet::new(),
        });
        Ok(())
    }

    /// The sequence repository for decoding a batch of CRAM containers whose single reference
    /// sequences are `ids` (in input order; containers with none or several are left out), with
    /// those sequences loaded. Loading here, before the containers are decoded in parallel, loads
    /// each sequence once (noodles' repository would load it once per thread that misses it), and
    /// outside the repository's lock, so threads decoding earlier batches keep going.
    ///
    /// When the batch moves to another reference sequence, the cache is cleared once it holds more
    /// than [`CRAM_CACHED_BASES`]; sorted input then loads each sequence about once (the one the
    /// batch starts on is loaded again after a clear).
    ///
    /// The sequences loaded here are checked against the header's `M5` on the thread pool (with
    /// one thread, right away), so the reading thread is not held up; the batch must not be used
    /// until every returned check has passed.
    pub fn cram_batch(
        &mut self,
        header: &noodles::sam::Header,
        ids: &[usize],
    ) -> io::Result<(fasta::Repository, Vec<Md5Check>)> {
        let c = self
            .cram
            .as_mut()
            .expect("CRAM input is opened with enable_cram");
        // A batch with no single-reference container (multi-reference slices, as samtools writes
        // for many small contigs, or unmapped reads) also clears it: those sequences are loaded
        // while decoding, so otherwise the cache would grow to the whole reference.
        if (ids.is_empty() || ids.iter().any(|&id| c.last != Some(id)))
            && c.loaded.load(Ordering::Relaxed) > CRAM_CACHED_BASES
        {
            c.repository.clear();
            c.cached.clear();
            c.loaded.store(0, Ordering::Relaxed);
        }
        if let Some(&last) = ids.last() {
            c.last = Some(last);
        }
        let mut checks = Vec::new();
        for &id in ids {
            if c.cached.contains(&id) {
                continue;
            }
            let (name, _) = header.reference_sequences().get_index(id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "container reference sequence ID is not in the header",
                )
            })?;
            let bases = c.loader.read(name)?;
            *c.staged.lock().unwrap() = Some((name.to_vec(), bases));
            let sequence = c.repository.get(name).transpose()?;
            c.cached.insert(id);
            if let (Some(sequence), Some(expected)) = (
                sequence,
                c.loader.md5.get(AsRef::<[u8]>::as_ref(name)).cloned(),
            ) {
                let (tx, rx) = std::sync::mpsc::sync_channel(1);
                let name = name.to_vec();
                let check = move || {
                    let _ = tx.send(check_md5(&name, sequence.as_ref().as_ref(), &expected));
                };
                if rayon::current_num_threads() > 1 {
                    rayon::spawn(check);
                } else {
                    check();
                }
                checks.push(rx);
            }
        }
        Ok((c.repository.clone(), checks))
    }

    /// Length of the named sequence, if the reference has it.
    pub fn length(&self, name: &str) -> Option<usize> {
        self.ids.get(name).map(|&i| self.sequences[i].1)
    }

    /// Append `len` bases of `name` starting at 0-based `start` to `out`, as stored in the FASTA
    /// (case kept).
    pub fn fetch(
        &mut self,
        name: &str,
        start: usize,
        len: usize,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        let invalid = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg);
        let &id = self
            .ids
            .get(name)
            .ok_or_else(|| invalid(format!("reference has no sequence named {name}")))?;
        let seq_len = self.sequences[id].1;
        let end = start.checked_add(len).filter(|&end| end <= seq_len);
        let Some(end) = end else {
            return Err(invalid(format!(
                "alignment ends at {} but reference sequence {name} has {seq_len} bases",
                start.saturating_add(len)
            )));
        };
        let mut pos = start;
        while pos < end {
            let block = pos / BLOCK;
            let from = pos - block * BLOCK;
            let to = (end - block * BLOCK).min(BLOCK);
            let bases = self.block(id, block)?;
            out.extend_from_slice(&bases[from..to]);
            pos = (block + 1) * BLOCK;
        }
        Ok(())
    }

    fn block(&mut self, id: usize, block: usize) -> io::Result<&[u8]> {
        let key = (id, block);
        if self.blocks.contains_key(&key) {
            if let Some(i) = self.order.iter().position(|k| *k == key) {
                self.order.remove(i);
            }
        } else {
            if self.order.len() == CACHED_BLOCKS {
                let oldest = self.order.pop_front().unwrap();
                self.blocks.remove(&oldest);
            }
            let seq_len = self.sequences[id].1;
            let first = block * BLOCK;
            let last = ((block + 1) * BLOCK).min(seq_len);
            let record = self.reader.index().as_ref()[id].clone();
            let bases = read_checked(self.reader.get_mut(), &record, first, last - first)?;
            self.blocks.insert(key, bases);
        }
        self.order.push_back(key);
        Ok(&self.blocks[&key])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A FASTA with one sequence longer than two blocks (wrapped at 60 bases) and its index.
    fn fasta(dir: &Path) -> (std::path::PathBuf, Vec<u8>) {
        let seq: Vec<u8> = (0..2 * BLOCK + 123).map(|i| b"ACGTacgtN"[i % 9]).collect();
        let mut text = b">c1 description\n".to_vec();
        for line in seq.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        text.extend_from_slice(b">c2\nAC\n");
        let path = dir.join("ref.fa");
        std::fs::write(&path, &text).unwrap();
        let index = fasta::fs::index(&path).unwrap();
        let mut w = fasta::fai::io::Writer::new(File::create(dir.join("ref.fa.fai")).unwrap());
        w.write_index(&index).unwrap();
        (path, seq)
    }

    #[test]
    fn fetch_spans_blocks_keeps_case_and_checks_bounds() {
        let dir = std::env::temp_dir().join("pmdshodh_reference_test");
        std::fs::create_dir_all(&dir).unwrap();
        let (path, seq) = fasta(&dir);
        let mut r = Reference::open(&path).unwrap();
        assert_eq!(r.length("c1"), Some(seq.len()));
        assert_eq!(r.length("c2"), Some(2));
        assert_eq!(r.length("c3"), None);
        for (start, len) in [
            (0, 10),
            (BLOCK - 5, 10),
            (2 * BLOCK + 100, 23),
            (7, 2 * BLOCK),
        ] {
            let mut out = Vec::new();
            r.fetch("c1", start, len, &mut out).unwrap();
            assert_eq!(out, &seq[start..start + len], "{start}+{len}");
        }
        let mut out = Vec::new();
        r.fetch("c2", 0, 2, &mut out).unwrap();
        assert_eq!(out, b"AC");
        assert!(r.fetch("c2", 1, 2, &mut Vec::new()).is_err());
        assert!(r.fetch("c3", 0, 1, &mut Vec::new()).is_err());
        // B3: a window length near usize::MAX is out of bounds, not an overflow.
        let e = r.fetch("c2", 1, usize::MAX, &mut Vec::new()).unwrap_err();
        assert!(
            e.to_string().contains("reference sequence c2 has 2 bases"),
            "{e}"
        );
    }

    /// Audit finding B1 (docs/audit-0.3.md §10): a FASTA that no longer matches its index was read
    /// with wrong bases from the second block on. Rewrapped lines, a changed header line and CRLF
    /// line ends under an LF index are now errors; a correctly indexed CRLF FASTA still works.
    #[test]
    fn fasta_that_does_not_match_its_index_is_an_error() {
        let dir = std::env::temp_dir().join("pmdshodh_reference_stale");
        std::fs::create_dir_all(&dir).unwrap();
        let (path, seq) = fasta(&dir);
        let fai = std::fs::read(dir.join("ref.fa.fai")).unwrap();
        let stale = |name: &str, text: &[u8]| {
            let p = dir.join(name);
            std::fs::write(&p, text).unwrap();
            std::fs::write(dir.join(format!("{name}.fai")), &fai).unwrap();
            p
        };
        let wrapped = |width: usize, end: &[u8], header: &[u8]| {
            let mut text = header.to_vec();
            text.extend_from_slice(end);
            for line in seq.chunks(width) {
                text.extend_from_slice(line);
                text.extend_from_slice(end);
            }
            text.extend_from_slice(b">c2");
            text.extend_from_slice(end);
            text.extend_from_slice(b"AC");
            text.extend_from_slice(end);
            text
        };
        for (name, text) in [
            ("rewrapped.fa", wrapped(70, b"\n", b">c1 description")),
            ("header.fa", wrapped(60, b"\n", b">c1 descr")),
            ("crlf.fa", wrapped(60, b"\r\n", b">c1 description")),
        ] {
            let p = stale(name, &text);
            for (start, len) in [(0, 10), (2 * BLOCK + 100, 23)] {
                let mut r = Reference::open(&p).unwrap();
                let e = r.fetch("c1", start, len, &mut Vec::new()).unwrap_err();
                assert!(
                    e.to_string()
                        .contains("does not have the line lengths its FASTA index gives"),
                    "{name} {start}: {e}"
                );
            }
        }
        // The original, and a CRLF FASTA indexed as it is, read correctly.
        let crlf = dir.join("crlf_indexed.fa");
        std::fs::write(&crlf, wrapped(60, b"\r\n", b">c1 description")).unwrap();
        let index = fasta::fs::index(&crlf).unwrap();
        fasta::fai::io::Writer::new(File::create(dir.join("crlf_indexed.fa.fai")).unwrap())
            .write_index(&index)
            .unwrap();
        for p in [&path, &crlf] {
            let mut r = Reference::open(p).unwrap();
            let mut out = Vec::new();
            r.fetch("c1", 2 * BLOCK + 100, 23, &mut out).unwrap();
            assert_eq!(out, &seq[2 * BLOCK + 100..2 * BLOCK + 123]);
        }
    }

    /// Whole sequences for CRAM come out uppercase and unwrapped, from a plain FASTA (fast path,
    /// with and without a final newline) and from a bgzip one (noodles' reader).
    #[test]
    fn cram_adapter_loads_whole_sequences_uppercase() {
        let dir = std::env::temp_dir().join("pmdshodh_reference_cram");
        std::fs::create_dir_all(&dir).unwrap();
        let (path, seq) = fasta(&dir);
        let load = |path: &Path, name: &[u8]| {
            let mut r = Reference::open(path).unwrap();
            r.enable_cram(path, &Default::default()).unwrap();
            let c = r.cram.as_ref().unwrap();
            let s = c.repository.get(name).unwrap().unwrap();
            (
                s.as_ref().as_ref().to_vec(),
                c.loaded.load(Ordering::Relaxed),
            )
        };
        let (c1, loaded) = load(&path, b"c1");
        assert_eq!(c1, seq.to_ascii_uppercase());
        assert_eq!(loaded, seq.len());
        assert_eq!(load(&path, b"c2").0, b"AC");

        // No final newline.
        let text = std::fs::read(&path).unwrap();
        let bare = dir.join("bare.fa");
        std::fs::write(&bare, &text[..text.len() - 1]).unwrap();
        std::fs::copy(dir.join("ref.fa.fai"), dir.join("bare.fa.fai")).unwrap();
        assert_eq!(load(&bare, b"c2").0, b"AC");

        // bgzip, through noodles' reader: one BGZF block, so the .gzi is empty, and the .fai is
        // the uncompressed file's.
        let small = b">s1
acgtNNacgt
AC
";
        std::fs::write(dir.join("small.fa"), small).unwrap();
        let index = fasta::fs::index(dir.join("small.fa")).unwrap();
        let gz = dir.join("small.fa.gz");
        let mut w = noodles::bgzf::io::Writer::new(File::create(&gz).unwrap());
        std::io::Write::write_all(&mut w, small).unwrap();
        w.finish().unwrap();
        let mut iw =
            fasta::fai::io::Writer::new(File::create(dir.join("small.fa.gz.fai")).unwrap());
        iw.write_index(&index).unwrap();
        noodles::bgzf::gzi::fs::write(dir.join("small.fa.gz.gzi"), &Default::default()).unwrap();
        assert_eq!(load(&gz, b"s1").0, b"ACGTNNACGTAC");
    }

    /// The CRAM cache is cleared past its limit also for a batch whose containers all have
    /// several reference sequences (no single-reference IDs), so it cannot grow to the whole
    /// reference.
    #[test]
    fn cram_cache_is_cleared_for_multi_reference_batches() {
        let dir = std::env::temp_dir().join("pmdshodh_reference_cram_clear");
        std::fs::create_dir_all(&dir).unwrap();
        let (path, _) = fasta(&dir);
        let mut r = Reference::open(&path).unwrap();
        let header = noodles::sam::Header::default();
        r.enable_cram(&path, &header).unwrap();
        let c = r.cram.as_ref().unwrap();
        c.repository.get(b"c2").unwrap().unwrap();
        c.loaded.store(CRAM_CACHED_BASES + 1, Ordering::Relaxed);
        r.cram_batch(&header, &[]).unwrap();
        let c = r.cram.as_ref().unwrap();
        assert_eq!(c.loaded.load(Ordering::Relaxed), 0);
        assert!(c.repository.is_empty());
    }

    #[test]
    fn missing_index_is_an_error() {
        let dir = std::env::temp_dir().join("pmdshodh_reference_noindex");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("ref.fa");
        std::fs::write(&path, ">c1\nACGT\n").unwrap();
        let _ = std::fs::remove_file(dir.join("ref.fa.fai"));
        assert!(Reference::open(&path).is_err());
    }
}
