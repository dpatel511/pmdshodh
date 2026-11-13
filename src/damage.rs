//! Damage tables: read bases (A, C, G, T) at reference C and G, by distance from each end of the
//! molecule. Reverse-strand reads are reverse-complemented first (molecule orientation), as in
//! PMDtools `--deamination` (docs/blackbox.md). Index 0 is the terminal base.

use crate::sites::{Settings, Treatment};

/// Read-base order in every count array.
pub const BASES: [u8; 4] = *b"ACGT";

/// Counts per position for one molecule end: `[ref C, ref G][position][read base A, C, G, T]`.
pub type EndCounts = [Vec<[u64; 4]>; 2];

#[derive(Debug, Clone, PartialEq)]
pub struct DamageCounts {
    pub range: usize,
    pub five_prime: EndCounts,
    pub three_prime: EndCounts,
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

/// Index of `b` in [`BASES`].
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl DamageCounts {
    pub fn new(range: usize) -> Self {
        let end = || [vec![[0; 4]; range], vec![[0; 4]; range]];
        DamageCounts {
            range,
            five_prime: end(),
            three_prime: end(),
        }
    }

    /// Add one read (SAM orientation in; `refs` from the site module).
    /// Sites follow the treatment and base-quality settings used for scoring. CpG context is
    /// taken from the reference; with `pmdtools_compat` it is taken as PMDtools does for its
    /// table (C: next *read* base is G; G: previous reference base is C, wrapping at the start).
    pub fn add_read(
        &mut self,
        seq: &[u8],
        qual: &[u8],
        refs: &[Option<u8>],
        reverse: bool,
        s: &Settings,
    ) {
        let n = seq.len();
        // Molecule orientation without copying: molecule position j is SAM position at(j),
        // complemented for reverse-strand reads.
        let at = |j: usize| if reverse { n - 1 - j } else { j };
        let mref = |j: usize| {
            let r = refs[at(j)];
            if reverse { r.map(complement) } else { r }
        };
        let mbase = |j: usize| {
            let b = seq[at(j)].to_ascii_uppercase();
            if reverse { complement(b) } else { b }
        };
        for i in 0..n {
            let (Some(r), Some(b)) = (mref(i), base_index(mbase(i))) else {
                continue;
            };
            if qual[at(i)] < s.min_baseq {
                continue;
            }
            let ref_index = match r {
                b'C' => 0,
                b'G' => 1,
                _ => continue,
            };
            let cpg = if ref_index == 0 {
                if s.pmdtools_compat {
                    i + 1 < n && mbase(i + 1) == b'G'
                } else {
                    i + 1 < n && mref(i + 1) == Some(b'G')
                }
            } else if i > 0 {
                mref(i - 1) == Some(b'C')
            } else {
                s.pmdtools_compat && mref(n - 1) == Some(b'C')
            };
            let terminal = if ref_index == 0 { i == 0 } else { i == n - 1 };
            let keep = match s.treatment {
                Treatment::None => true,
                Treatment::Half => terminal || cpg,
                Treatment::Full => cpg,
            };
            if !keep {
                continue;
            }
            if i < self.range {
                self.five_prime[ref_index][i][b] += 1;
            }
            if n - 1 - i < self.range {
                self.three_prime[ref_index][n - 1 - i][b] += 1;
            }
        }
    }

    /// C→T rate at the 5′ terminal base, its binomial standard error, and the number of sites.
    pub fn first_position(&self) -> (Option<f64>, Option<f64>, u64) {
        let Some(c) = self.five_prime[0].first() else {
            return (None, None, 0);
        };
        let n: u64 = c.iter().sum();
        if n == 0 {
            return (None, None, 0);
        }
        let p = c[3] as f64 / n as f64;
        (Some(p), Some((p * (1.0 - p) / n as f64).sqrt()), n)
    }
}

impl DamageCounts {
    /// Add another set of counts (same range). Integer sums, so the order does not matter.
    pub fn merge(&mut self, other: &DamageCounts) {
        for (mine, theirs) in [
            (&mut self.five_prime, &other.five_prime),
            (&mut self.three_prime, &other.three_prime),
        ] {
            for (a, b) in mine.iter_mut().zip(theirs) {
                for (x, y) in a.iter_mut().zip(b) {
                    for (u, v) in x.iter_mut().zip(y) {
                        *u += v;
                    }
                }
            }
        }
    }
}
