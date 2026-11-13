//! Turn one alignment into informative sites and sum their score terms.
//!
//! Behaviour follows PMDtools 0.60 as measured in docs/blackbox.md (sections 1.2, 1.4, 1.5, 1.6b),
//! except for the departures listed in docs/design.md section 5.

use crate::model::{Params, site_term_ds, site_term_ss};

/// Library preparation, which decides which mismatches count as damage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Library {
    /// C→T from the 5′ end and G→A from the 3′ end (PMDtools default).
    DoubleStranded,
    /// C→T in molecule orientation, measured from both ends (PMDtools `--ss`).
    SingleStranded,
}

/// UDG treatment, which decides which C/G sites keep damage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Treatment {
    /// All sites (PMDtools default, `--UDGminus`).
    None,
    /// CpG sites plus the terminal base (PMDtools `--UDGhalf`).
    Half,
    /// CpG sites only (PMDtools `--CpG` / `--UDGplus`).
    Full,
}

/// Scoring settings for one run.
#[derive(Debug, Clone, Copy)]
pub struct Settings {
    pub library: Library,
    pub treatment: Treatment,
    pub params: Params,
    /// Bases with a lower phred quality are ignored (PMDtools `-q`).
    pub min_baseq: u8,
    /// Rebuild the reference the way PMDtools 0.60 does, including its bug with combined CIGAR
    /// operations (`--pmdtools-compat`; docs/blackbox.md 1.6c). Off by default.
    pub pmdtools_compat: bool,
}

/// One CIGAR operation with its length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// `M`.
    Match(usize),
    /// `=`.
    SequenceMatch(usize),
    /// `X`.
    SequenceMismatch(usize),
    Insertion(usize),
    Deletion(usize),
    SoftClip(usize),
    HardClip(usize),
    /// `N` or `P`: not supported, as in PMDtools.
    Unsupported,
}

/// Score of one read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadScore {
    pub pmds: f64,
    /// Number of informative sites that contributed (0 gives PMDS 0, as in PMDtools).
    pub n_sites: u32,
}

/// Reference base at each SEQ position, rebuilt from the CIGAR and the MD tag.
/// Inserted and soft-clipped positions have no reference base (`None`).
///
/// MD deletions take exactly as many letters as the CIGAR deletion has, so MD strings without the
/// `0` separator after a deletion (written by MALT; PMDtools issue #5) are read correctly.
/// Any other disagreement between CIGAR, SEQ and MD is an error.
pub fn reference_along_read(
    seq: &[u8],
    cigar: &[CigarOp],
    md: &[u8],
) -> Result<Vec<Option<u8>>, String> {
    let read_len: usize = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::SequenceMatch(n)
            | CigarOp::SequenceMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n) => *n,
            _ => 0,
        })
        .sum();
    if read_len != seq.len() {
        return Err(format!(
            "CIGAR covers {read_len} read bases but SEQ has {}",
            seq.len()
        ));
    }

    let mut md = MdCursor { md, pos: 0, run: 0 };
    let mut refs = Vec::with_capacity(seq.len());
    for op in cigar {
        match *op {
            CigarOp::Match(n) | CigarOp::SequenceMatch(n) | CigarOp::SequenceMismatch(n) => {
                for _ in 0..n {
                    let read_base = seq[refs.len()].to_ascii_uppercase();
                    refs.push(Some(md.next_aligned()?.unwrap_or(read_base)));
                }
            }
            CigarOp::Insertion(n) | CigarOp::SoftClip(n) => refs.extend((0..n).map(|_| None)),
            CigarOp::Deletion(n) => md.deletion(n)?,
            CigarOp::HardClip(_) => {}
            CigarOp::Unsupported => return Err("unsupported CIGAR operation (N or P)".into()),
        }
    }
    md.finish()?;
    Ok(refs)
}

/// Reference base at each SEQ position as PMDtools 0.60 rebuilds it, bug included
/// (docs/blackbox.md 1.6c). Only for `--pmdtools-compat`; call [`reference_along_read`] first so
/// that invalid records are still rejected.
///
/// - No `I` and no `S` in the CIGAR: the MD reference without deleted bases is laid along SEQ.
/// - Otherwise: gaps are placed by a counter that `M` and `X` advance, `I` and `S` mark and
///   advance, `D` also advances and `=` ignores; the MD reference *with* deleted bases fills the
///   remaining positions in order.
///
/// MD is read as PMDtools reads it: after `^`, every letter up to the next digit is deleted, and
/// lowercase letters are kept as they are, so they never count as C or G (sites, CpG context or
/// the damage table). The default reconstruction reads them as uppercase.
///
/// For a reverse-strand read, PMDtools' reverse complement keeps only `ACGTNRY` and `acgt` and
/// drops every other reference character (e.g. lowercase `n`, most IUPAC codes), so the rest of
/// the reference shifts towards the molecule's 5′ end (the SAM end) and the positions left over at
/// the SAM start have no reference base (docs/blackbox.md, 2026-11-13).
pub fn reference_along_read_pmdtools(
    seq: &[u8],
    cigar: &[CigarOp],
    md: &[u8],
    reverse: bool,
) -> Vec<Option<u8>> {
    let refs = pmdtools_reference_sam_order(seq, cigar, md);
    let dropped = |r: &Option<u8>| r.is_some_and(|b| !b"ACGTNRYacgt".contains(&b));
    if !reverse || !refs.iter().any(dropped) {
        return refs;
    }
    let kept: Vec<Option<u8>> = refs.iter().copied().filter(|r| !dropped(r)).collect();
    std::iter::repeat_n(None, refs.len() - kept.len())
        .chain(kept)
        .collect()
}

fn pmdtools_reference_sam_order(seq: &[u8], cigar: &[CigarOp], md: &[u8]) -> Vec<Option<u8>> {
    // MD reference: None = matches the read, Some(base) = mismatch or deleted base.
    let mut expansion: Vec<(Option<u8>, bool)> = Vec::new(); // (base, is_deleted)
    let mut i = 0;
    while i < md.len() {
        if md[i].is_ascii_digit() {
            let start = i;
            while i < md.len() && md[i].is_ascii_digit() {
                i += 1;
            }
            let n: usize = std::str::from_utf8(&md[start..i])
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
            expansion.extend(std::iter::repeat_n((None, false), n));
        } else if md[i] == b'^' {
            i += 1;
            while i < md.len() && md[i].is_ascii_alphabetic() {
                expansion.push((Some(md[i]), true));
                i += 1;
            }
        } else {
            expansion.push((Some(md[i]), false));
            i += 1;
        }
    }

    let n = seq.len();
    let fill = |e: Option<&(Option<u8>, bool)>, i: usize| -> Option<u8> {
        e.map(|(base, _)| base.unwrap_or(seq[i].to_ascii_uppercase()))
    };
    let simple = !cigar
        .iter()
        .any(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::SoftClip(_)));
    if simple {
        let aligned: Vec<_> = expansion.iter().filter(|(_, del)| !del).collect();
        return (0..n).map(|i| fill(aligned.get(i).copied(), i)).collect();
    }

    let mut gap = vec![false; n];
    let mut pos = 0usize;
    for op in cigar {
        match *op {
            CigarOp::Match(k) | CigarOp::SequenceMismatch(k) | CigarOp::Deletion(k) => pos += k,
            CigarOp::Insertion(k) | CigarOp::SoftClip(k) => {
                for g in gap.iter_mut().skip(pos).take(k) {
                    *g = true;
                }
                pos += k;
            }
            CigarOp::SequenceMatch(_) | CigarOp::HardClip(_) | CigarOp::Unsupported => {}
        }
    }
    let mut next = expansion.iter();
    (0..n)
        .map(|i| if gap[i] { None } else { fill(next.next(), i) })
        .collect()
}

/// Reads the MD tag in step with the CIGAR.
struct MdCursor<'a> {
    md: &'a [u8],
    pos: usize,
    /// Matching bases left in the current MD number.
    run: usize,
}

impl MdCursor<'_> {
    fn read_number(&mut self) -> Option<usize> {
        let start = self.pos;
        while self.pos < self.md.len() && self.md[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos == start {
            return None;
        }
        std::str::from_utf8(&self.md[start..self.pos])
            .ok()?
            .parse()
            .ok()
    }

    /// Next aligned base: `None` if it matches the read, or the reference base if it is a mismatch.
    fn next_aligned(&mut self) -> Result<Option<u8>, String> {
        loop {
            if self.run > 0 {
                self.run -= 1;
                return Ok(None);
            }
            if let Some(n) = self.read_number() {
                self.run = n;
                continue;
            }
            return match self.md.get(self.pos) {
                Some(b) if b.is_ascii_alphabetic() => {
                    self.pos += 1;
                    Ok(Some(b.to_ascii_uppercase()))
                }
                Some(b'^') => Err("MD has a deletion where the CIGAR has an aligned base".into()),
                Some(_) => Err("MD contains an invalid character".into()),
                None => Err("MD is shorter than the CIGAR".into()),
            };
        }
    }

    /// Consume a deletion of `n` bases: `^` followed by exactly `n` letters.
    fn deletion(&mut self, n: usize) -> Result<(), String> {
        if self.run > 0 {
            return Err("CIGAR has a deletion inside an MD match run".into());
        }
        while let Some(k) = self.read_number() {
            if k != 0 {
                return Err("CIGAR has a deletion where MD has matching bases".into());
            }
        }
        if self.md.get(self.pos) != Some(&b'^') {
            return Err("CIGAR has a deletion that MD does not have".into());
        }
        self.pos += 1;
        let letters = &self.md[self.pos..];
        if letters.len() < n || !letters[..n].iter().all(u8::is_ascii_alphabetic) {
            return Err("MD deletion is shorter than the CIGAR deletion".into());
        }
        self.pos += n;
        Ok(())
    }

    /// Everything left in MD must be zero-length matches.
    fn finish(&mut self) -> Result<(), String> {
        while self.run == 0 {
            match self.read_number() {
                Some(k) => self.run = k,
                None => break,
            }
        }
        if self.run > 0 || self.pos < self.md.len() {
            return Err("MD is longer than the CIGAR".into());
        }
        Ok(())
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// True if the reference base at the next SEQ position is G.
/// With `wrap`, a C at the last position takes the first position's reference base as its next base.
fn g_follows(refs: &[Option<u8>], i: usize, wrap: bool) -> bool {
    if i + 1 < refs.len() {
        refs[i + 1] == Some(b'G')
    } else {
        wrap && refs.first() == Some(&Some(b'G'))
    }
}

/// True if the reference base at the previous SEQ position is C.
/// With `wrap`, a G at the first position takes the last position's reference base as its
/// previous base.
fn c_precedes(refs: &[Option<u8>], i: usize, wrap: bool) -> bool {
    if i > 0 {
        refs[i - 1] == Some(b'C')
    } else {
        wrap && refs.last() == Some(&Some(b'C'))
    }
}

/// Score one read. `seq`, `qual` and `refs` are in SAM orientation and have the same length;
/// `qual` holds raw phred values.
pub fn score_read(
    seq: &[u8],
    qual: &[u8],
    refs: &[Option<u8>],
    reverse: bool,
    s: &Settings,
) -> ReadScore {
    let n = seq.len();
    let mut score = ReadScore {
        pmds: 0.0,
        n_sites: 0,
    };
    match s.library {
        Library::DoubleStranded => {
            // Scored in SAM orientation; the ds model is symmetric under reverse complement.
            // PMDtools 0.60 works in molecule orientation, where a G at the molecule's first base
            // takes the molecule's last base as its previous base. In SAM orientation that is the
            // first G of a forward read, or the last C of a reverse read (docs/blackbox.md 1.6c).
            let wrap_g = s.pmdtools_compat && !reverse;
            let wrap_c = s.pmdtools_compat && reverse;
            for i in 0..n {
                let (Some(r), b, q) = (refs[i], seq[i].to_ascii_uppercase(), qual[i]) else {
                    continue;
                };
                if q < s.min_baseq {
                    continue;
                }
                if r == b'C' && (b == b'C' || b == b'T') {
                    let keep = match s.treatment {
                        Treatment::None => true,
                        Treatment::Half => i == 0 || g_follows(refs, i, wrap_c),
                        Treatment::Full => g_follows(refs, i, wrap_c),
                    };
                    if keep {
                        score.pmds += site_term_ds(b == b'C', (i + 1) as u32, q, &s.params);
                        score.n_sites += 1;
                    }
                } else if r == b'G' && (b == b'G' || b == b'A') {
                    let keep = match s.treatment {
                        Treatment::None => true,
                        Treatment::Half => i == n - 1 || c_precedes(refs, i, wrap_g),
                        Treatment::Full => c_precedes(refs, i, wrap_g),
                    };
                    if keep {
                        score.pmds += site_term_ds(b == b'G', (n - i) as u32, q, &s.params);
                        score.n_sites += 1;
                    }
                }
            }
        }
        Library::SingleStranded => {
            // Scored in molecule orientation: reverse-strand reads are reverse-complemented.
            let (mrefs, mseq, mqual): (Vec<Option<u8>>, Vec<u8>, Vec<u8>) = if reverse {
                (
                    refs.iter().rev().map(|r| r.map(complement)).collect(),
                    seq.iter()
                        .rev()
                        .map(|b| complement(b.to_ascii_uppercase()))
                        .collect(),
                    qual.iter().rev().copied().collect(),
                )
            } else {
                (
                    refs.to_vec(),
                    seq.iter().map(u8::to_ascii_uppercase).collect(),
                    qual.to_vec(),
                )
            };
            for i in 0..n {
                let (Some(r), b, q) = (mrefs[i], mseq[i], mqual[i]) else {
                    continue;
                };
                if q < s.min_baseq || r != b'C' || !(b == b'C' || b == b'T') {
                    continue;
                }
                let keep = match s.treatment {
                    Treatment::None => true,
                    // Only the molecule's 5′ first base is exempt (docs/blackbox.md 1.6b).
                    Treatment::Half => i == 0 || g_follows(&mrefs, i, false),
                    Treatment::Full => g_follows(&mrefs, i, false),
                };
                if keep {
                    let (z, y) = ((i + 1) as u32, (n - i) as u32);
                    score.pmds += site_term_ss(b == b'C', z, y, q, &s.params);
                    score.n_sites += 1;
                }
            }
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs(seq: &str, cigar: &[CigarOp], md: &str) -> Result<String, String> {
        let r = reference_along_read(seq.as_bytes(), cigar, md.as_bytes())?;
        Ok(r.iter().map(|b| b.map_or('-', char::from)).collect())
    }

    #[test]
    fn mismatches_and_matches_from_md() {
        assert_eq!(refs("TAAC", &[CigarOp::Match(4)], "0C3"), Ok("CAAC".into()));
    }

    #[test]
    fn clips_and_insertions_have_no_reference_base() {
        let cigar = [
            CigarOp::SoftClip(2),
            CigarOp::Match(2),
            CigarOp::Insertion(1),
            CigarOp::Match(1),
            CigarOp::HardClip(5),
        ];
        assert_eq!(refs("GGTAAC", &cigar, "0C2"), Ok("--CA-C".into()));
    }

    #[test]
    fn deletion_with_and_without_zero_separator() {
        // PMDtools issue #5: MALT writes ^CCTTT0C instead of ^CCTT0T0C.
        let cigar = [CigarOp::Match(4), CigarOp::Deletion(4), CigarOp::Match(3)];
        let spec = refs("AAAACTA", &cigar, "4^CCTT0T0C1");
        let malt = refs("AAAACTA", &cigar, "4^CCTTT0C1");
        assert_eq!(spec, Ok("AAAATCA".into()));
        assert_eq!(malt, spec);
    }

    #[test]
    fn md_that_contradicts_cigar_is_an_error() {
        assert!(refs("TAAC", &[CigarOp::Match(4)], "0C2").is_err()); // MD too short
        assert!(refs("TAAC", &[CigarOp::Match(4)], "0C4").is_err()); // MD too long
        assert!(refs("TAAC", &[CigarOp::Match(3)], "0C2").is_err()); // CIGAR vs SEQ length
        assert!(refs("TAAC", &[CigarOp::Match(4)], "0^C3").is_err()); // deletion not in CIGAR
        assert!(refs("TAAC", &[CigarOp::Match(2), CigarOp::Unsupported], "4").is_err());
    }

    fn compat(seq: &str, cigar: &[CigarOp], md: &str) -> String {
        reference_along_read_pmdtools(seq.as_bytes(), cigar, md.as_bytes(), false)
            .iter()
            .map(|b| b.map_or('-', char::from))
            .collect()
    }

    // Expected strings are PMDtools 0.60 --printalignments reference lines (docs/blackbox.md probes 15–16).
    #[test]
    fn compat_reproduces_pmdtools_reference_lines() {
        use CigarOp::*;
        let read = "ACGTACGTACTTCCAAGGTT";
        assert_eq!(
            compat(read, &[Match(10), Deletion(2), Match(10)], "10^GG10"),
            "ACGTACGTACTTCCAAGGTT"
        );
        assert_eq!(
            compat(
                &format!("NN{read}"),
                &[SoftClip(2), Match(10), Deletion(2), Match(10)],
                "10^GG10"
            ),
            "--ACGTACGTACGGCCAAGGTT"
        );
        assert_eq!(
            compat(
                &format!("{read}NNN"),
                &[Match(10), Deletion(2), Match(10), SoftClip(3)],
                "10^GG10"
            ),
            "ACGTACGTACGGCCAAGGTTNN-"
        );
        assert_eq!(
            compat(
                "ACGTACGTACNNTTCCAAGGTT",
                &[SequenceMatch(10), Insertion(2), SequenceMatch(10)],
                "20"
            ),
            "--GTACGTACNNTTCCAAGGTT"
        );
    }

    // Departure from PMDtools (docs/blackbox.md 1.6c): the default reconstruction is correct.
    #[test]
    fn default_reconstruction_fixes_clip_plus_deletion() {
        let cigar = [
            CigarOp::SoftClip(2),
            CigarOp::Match(10),
            CigarOp::Deletion(2),
            CigarOp::Match(10),
        ];
        let read = "NNACGTACGTACTTCCAAGGTT";
        assert_eq!(
            refs(read, &cigar, "10^GG10"),
            Ok("--ACGTACGTACTTCCAAGGTT".into())
        );
        assert_ne!(compat(read, &cigar, "10^GG10"), "--ACGTACGTACTTCCAAGGTT");
    }

    fn default_settings(library: Library, treatment: Treatment) -> Settings {
        Settings {
            library,
            treatment,
            params: Params::default(),
            min_baseq: 0,
            pmdtools_compat: false,
        }
    }

    fn score(seq: &str, ref_seq: &str, reverse: bool, s: &Settings) -> ReadScore {
        let refs: Vec<Option<u8>> = ref_seq.bytes().map(Some).collect();
        score_read(seq.as_bytes(), &vec![40; seq.len()], &refs, reverse, s)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    // Expected values are PMDtools 0.60 outputs from docs/blackbox.md probes 1, 6, 9 and 12.
    #[test]
    fn terminal_damage_scores_like_pmdtools() {
        let ds = default_settings(Library::DoubleStranded, Treatment::None);
        let a = "A".repeat(29);
        assert!(close(
            score(&format!("T{a}"), &format!("C{a}"), false, &ds).pmds,
            5.029717285061266
        ));
        assert!(close(
            score(&format!("{a}A"), &format!("{a}G"), false, &ds).pmds,
            5.029717285061266
        ));
    }

    #[test]
    fn no_informative_sites_scores_zero() {
        let ds = default_settings(Library::DoubleStranded, Treatment::None);
        let r = score("AAAA", "AAAA", false, &ds);
        assert_eq!(
            r,
            ReadScore {
                pmds: 0.0,
                n_sites: 0
            }
        );
    }

    #[test]
    fn single_stranded_uses_molecule_orientation() {
        let ss = default_settings(Library::SingleStranded, Treatment::None);
        let a = "A".repeat(29);
        // Forward: C→T at the first base. Reverse: SAM G→A at the last base is the molecule's 5′ C→T.
        let fwd = score(&format!("T{a}"), &format!("C{a}"), false, &ss).pmds;
        let rev = score(&format!("{a}A"), &format!("{a}G"), true, &ss).pmds;
        assert!(close(fwd, 9.681909363978177));
        assert!(close(rev, 9.681909363978177));
        // G→A on a forward read is not damage in ss mode.
        assert_eq!(
            score(&format!("{a}A"), &format!("{a}G"), false, &ss).n_sites,
            0
        );
    }

    #[test]
    fn cpg_modes_select_sites_like_pmdtools() {
        let full = default_settings(Library::DoubleStranded, Treatment::Full);
        let half = default_settings(Library::DoubleStranded, Treatment::Half);
        let a = "A".repeat(28);
        // Terminal C→T, not CpG: dropped by full, kept by half.
        assert_eq!(
            score(&format!("TA{a}"), &format!("CA{a}"), false, &full).n_sites,
            0
        );
        assert!(close(
            score(&format!("TA{a}"), &format!("CA{a}"), false, &half).pmds,
            5.029717285061266
        ));
        // Interior CpG C→T (z = 5) plus its G (match, z = 25): kept by both.
        let read = format!("AAAATG{}", "A".repeat(24));
        let reference = format!("AAAACG{}", "A".repeat(24));
        assert!(close(
            score(&read, &reference, false, &full).pmds,
            3.7003221261207053
        ));
        assert!(close(
            score(&read, &reference, false, &half).pmds,
            3.7003221261207053
        ));
    }

    // Departure from PMDtools (docs/blackbox.md 1.6c): CpG context does not wrap around the read
    // ends by default. PMDtools (and --pmdtools-compat) treats a forward read's first G, or a
    // reverse read's last C, as CpG when the base at the other end completes the pair.
    #[test]
    fn cpg_context_wraps_only_in_compat_mode() {
        let mut s = default_settings(Library::DoubleStranded, Treatment::Full);
        let (read, reference) = ("GAAAC", "GAAAC");
        assert_eq!(score(read, reference, false, &s).n_sites, 0);
        assert_eq!(score(read, reference, true, &s).n_sites, 0);
        s.pmdtools_compat = true;
        let fwd = score(read, reference, false, &s);
        let rev = score(read, reference, true, &s);
        assert_eq!((fwd.n_sites, rev.n_sites), (1, 1));
        // Forward: the first G (3′ distance 5). Reverse: the last C (5′ distance 5). Same term.
        assert!(close(fwd.pmds, site_term_ds(true, 5, 40, &s.params)));
        assert!(close(rev.pmds, fwd.pmds));
    }

    #[test]
    fn low_quality_bases_are_ignored() {
        let mut s = default_settings(Library::DoubleStranded, Treatment::None);
        s.min_baseq = 41;
        assert_eq!(score("TAAA", "CAAA", false, &s).n_sites, 0);
    }
}
