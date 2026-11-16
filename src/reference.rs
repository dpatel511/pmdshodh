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

use noodles::core::{Position, Region};
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

    fn load(&mut self, name: &[u8]) -> io::Result<Vec<u8>> {
        let index: &[fasta::fai::Record] = self.reader.index().as_ref();
        let record = index
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == name)
            .cloned();
        let bases = match (&mut self.plain, record) {
            (Some(file), Some(record)) => read_plain(file, &record),
            _ => self
                .reader
                .query(&Region::new(name, ..))
                .map(|record| record.sequence().as_ref().to_ascii_uppercase()),
        }?;
        if let Some(expected) = self.md5.get(name) {
            use md5::Digest;
            let actual: String = md5::Md5::digest(&bases)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            if actual != *expected {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference sequence {} does not match the CRAM header's M5 checksum \
                         (is -r the reference the CRAM was written with?)",
                        String::from_utf8_lossy(name)
                    ),
                ));
            }
        }
        Ok(bases)
    }
}

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

/// Read a sequence of an uncompressed FASTA at its index offset in one read, dropping line breaks
/// and uppercasing in place.
fn read_plain(file: &mut File, record: &fasta::fai::Record) -> io::Result<Vec<u8>> {
    use std::io::{Read, Seek, SeekFrom};
    let len = record.length() as usize;
    let (line_bases, line_width) = (
        record.line_base_count().get() as usize,
        record.line_width().get() as usize,
    );
    let span = if len == 0 {
        0
    } else {
        let last_line = (len - 1) / line_bases;
        last_line * line_width + (len - last_line * line_bases)
    };
    let mut bases = vec![0; span];
    file.seek(SeekFrom::Start(record.position()))?;
    file.read_exact(&mut bases)?;
    let mut kept = 0;
    for i in 0..span {
        let b = bases[i];
        if b != b'\n' && b != b'\r' {
            bases[kept] = b.to_ascii_uppercase();
            kept += 1;
        }
    }
    if kept != len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "sequence {} does not have the line lengths its FASTA index gives",
                record.name()
            ),
        ));
    }
    bases.truncate(len);
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
    pub fn cram_batch(
        &mut self,
        header: &noodles::sam::Header,
        ids: &[usize],
    ) -> io::Result<fasta::Repository> {
        let c = self
            .cram
            .as_mut()
            .expect("CRAM input is opened with enable_cram");
        if ids.iter().any(|&id| c.last != Some(id))
            && c.loaded.load(Ordering::Relaxed) > CRAM_CACHED_BASES
        {
            c.repository.clear();
            c.cached.clear();
            c.loaded.store(0, Ordering::Relaxed);
        }
        if let Some(&last) = ids.last() {
            c.last = Some(last);
        }
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
            let bases = c.loader.load(name)?;
            *c.staged.lock().unwrap() = Some((name.to_vec(), bases));
            c.repository.get(name).transpose()?;
            c.cached.insert(id);
        }
        Ok(c.repository.clone())
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
        if start + len > seq_len {
            return Err(invalid(format!(
                "alignment ends at {} but reference sequence {name} has {seq_len} bases",
                start + len
            )));
        }
        let mut pos = start;
        while pos < start + len {
            let block = pos / BLOCK;
            let from = pos - block * BLOCK;
            let to = (start + len - block * BLOCK).min(BLOCK);
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
            let (name, seq_len) = &self.sequences[id];
            let first = block * BLOCK + 1;
            let last = ((block + 1) * BLOCK).min(*seq_len);
            let position = |p: usize| Position::try_from(p).map_err(io::Error::other);
            let region = Region::new(name.as_str(), position(first)?..=position(last)?);
            let record = self.reader.query(&region)?;
            self.blocks.insert(key, record.sequence().as_ref().to_vec());
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
