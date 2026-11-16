//! Reference FASTA access for `-r` (docs/design.md §11).
//!
//! Sequences are read through the FASTA index (`.fai`, plus `.gzi` for bgzip) in fixed-size
//! blocks kept in a small least-recently-used cache. Memory stays bounded for sorted and unsorted
//! input, and no whole chromosome or genome is loaded.

use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io;
use std::path::Path;

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
        })
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
