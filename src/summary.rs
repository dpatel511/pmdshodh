//! Run summary written as JSON for MultiQC (docs/design.md section 6). Changing any field requires
//! bumping `FORMAT_VERSION` (CLAUDE.md §7). No timestamps: identical runs give identical files.

use std::collections::BTreeMap;

use noodles::sam::alignment::record::Flags;
use serde_json::{Value, json};

use crate::damage::{BASES, DamageCounts, EndCounts};
use crate::model::NULL_BACKGROUND;
use crate::reads::{ReadFilters, SkipReason};
use crate::sites::{Library, ReadScore, Settings, Treatment};

pub const FORMAT_VERSION: u32 = 1;
const HIST_MIN: f64 = -20.0;
const HIST_MAX: f64 = 40.0;
const HIST_BIN: f64 = 0.5;

/// Facts about the run that are not counts.
pub struct RunInfo<'a> {
    pub subcommand: &'a str,
    pub command: String,
    pub input: String,
    pub settings: Settings,
    pub filters: ReadFilters,
    pub threshold: Option<f64>,
    pub upper_threshold: Option<f64>,
}

pub struct Summary {
    input: u64,
    scored: u64,
    kept: u64,
    removed_by_threshold: u64,
    no_informative_sites: u64,
    skipped: [u64; SkipReason::ALL.len()],
    flags: BTreeMap<&'static str, u64>,
    histogram: Vec<u64>,
    underflow: u64,
    overflow: u64,
    length_scored: BTreeMap<usize, u64>,
    length_kept: BTreeMap<usize, u64>,
    pub damage: DamageCounts,
}

impl Summary {
    pub fn new(range: usize) -> Self {
        let bins = ((HIST_MAX - HIST_MIN) / HIST_BIN).round() as usize;
        Summary {
            input: 0,
            scored: 0,
            kept: 0,
            removed_by_threshold: 0,
            no_informative_sites: 0,
            skipped: [0; SkipReason::ALL.len()],
            flags: [
                "secondary",
                "supplementary",
                "duplicate",
                "qc_fail",
                "paired",
                "reverse",
            ]
            .into_iter()
            .map(|k| (k, 0))
            .collect(),
            histogram: vec![0; bins],
            underflow: 0,
            overflow: 0,
            length_scored: BTreeMap::new(),
            length_kept: BTreeMap::new(),
            damage: DamageCounts::new(range),
        }
    }

    /// Count an input record by its flags (before any filtering).
    pub fn add_input(&mut self, flags: Flags) {
        self.input += 1;
        for (key, set) in [
            ("secondary", flags.is_secondary()),
            ("supplementary", flags.is_supplementary()),
            ("duplicate", flags.is_duplicate()),
            ("qc_fail", flags.is_qc_fail()),
            ("paired", flags.is_segmented()),
            ("reverse", flags.is_reverse_complemented()),
        ] {
            if set {
                *self.flags.get_mut(key).unwrap() += 1;
            }
        }
    }

    pub fn add_skipped(&mut self, reason: SkipReason) {
        self.skipped[reason as usize] += 1;
    }

    /// Count a scored read. `kept` is `None` when there is no filtering (score, damage).
    pub fn add_scored(&mut self, score: &ReadScore, length: usize, kept: Option<bool>) {
        self.scored += 1;
        if score.n_sites == 0 {
            self.no_informative_sites += 1;
        }
        if score.pmds < HIST_MIN {
            self.underflow += 1;
        } else if score.pmds >= HIST_MAX {
            self.overflow += 1;
        } else {
            let last = self.histogram.len() - 1;
            let i = ((score.pmds - HIST_MIN) / HIST_BIN).floor() as usize;
            self.histogram[i.min(last)] += 1;
        }
        *self.length_scored.entry(length).or_default() += 1;
        match kept {
            Some(true) => {
                self.kept += 1;
                *self.length_kept.entry(length).or_default() += 1;
            }
            Some(false) => self.removed_by_threshold += 1,
            None => {}
        }
    }

    pub fn skipped(&self) -> &[u64] {
        &self.skipped
    }

    pub fn to_json(&self, info: &RunInfo) -> Value {
        let s = &info.settings;
        let filtering = info.subcommand == "filter";
        let end = |counts: &EndCounts| {
            json!({
                "C": base_columns(&counts[0]),
                "G": base_columns(&counts[1]),
            })
        };
        let (rate, se, n) = self.damage.first_position();
        let lengths = |m: &BTreeMap<usize, u64>| -> Value {
            m.iter().map(|(k, v)| (k.to_string(), json!(v))).collect()
        };
        json!({
            "format_version": FORMAT_VERSION,
            "tool": {"name": "pmdshodh", "version": env!("CARGO_PKG_VERSION")},
            "subcommand": info.subcommand,
            "command": info.command,
            "input": info.input,
            "mode": {
                "library": match s.library { Library::DoubleStranded => "ds", Library::SingleStranded => "ss" },
                "treatment": match s.treatment { Treatment::None => "none", Treatment::Half => "half", Treatment::Full => "full" },
                "pmdtools_compat": s.pmdtools_compat,
            },
            "model": {
                "pmd_p": s.params.p,
                "pmd_c": s.params.c,
                "pi_ancient": s.params.pi_ancient,
                "pi_contamination": s.params.pi_contamination,
                "null_background": NULL_BACKGROUND,
            },
            "filters": {
                "threshold": info.threshold,
                "upper_threshold": info.upper_threshold,
                "min_mapq": info.filters.min_mapq,
                "min_baseq": s.min_baseq,
                "no_clips": info.filters.no_clips,
                "no_indels": info.filters.no_indels,
            },
            "reads": {
                "input": self.input,
                "scored": self.scored,
                "kept": if filtering { json!(self.kept) } else { Value::Null },
                "removed_by_threshold": if filtering { json!(self.removed_by_threshold) } else { Value::Null },
                "no_informative_sites": self.no_informative_sites,
                "skipped": SkipReason::ALL.iter().map(|r| (r.name().to_string(), json!(self.skipped[*r as usize]))).collect::<serde_json::Map<_, _>>(),
                "flags": self.flags,
            },
            "score_histogram": {
                "bin_width": HIST_BIN,
                "min": HIST_MIN,
                "max": HIST_MAX,
                "underflow": self.underflow,
                "overflow": self.overflow,
                "counts": self.histogram,
            },
            "length_histogram": {
                "scored": lengths(&self.length_scored),
                "kept": if filtering { lengths(&self.length_kept) } else { Value::Null },
            },
            "damage": {
                "range": self.damage.range,
                "orientation": "molecule",
                "reads_used": "scored",
                "five_prime": end(&self.damage.five_prime),
                "three_prime": end(&self.damage.three_prime),
                "first_position": {"c_to_t": rate, "se": se, "n": n},
            },
        })
    }
}

/// `{ "A": [...], "C": [...], "G": [...], "T": [...] }`, position 0 = terminal base.
fn base_columns(rows: &[[u64; 4]]) -> Value {
    BASES
        .iter()
        .enumerate()
        .map(|(j, b)| {
            (
                (*b as char).to_string(),
                json!(rows.iter().map(|r| r[j]).collect::<Vec<_>>()),
            )
        })
        .collect::<serde_json::Map<_, _>>()
        .into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Params;

    fn settings() -> Settings {
        Settings {
            library: Library::DoubleStranded,
            treatment: Treatment::None,
            params: Params::default(),
            min_baseq: 0,
            pmdtools_compat: false,
        }
    }

    #[test]
    fn histogram_bins_and_overflow() {
        let mut s = Summary::new(5);
        for pmds in [-25.0, -20.0, -0.25, 0.0, 0.49, 39.99, 40.0] {
            s.add_scored(&ReadScore { pmds, n_sites: 1 }, 30, None);
        }
        assert_eq!((s.underflow, s.overflow), (1, 1));
        assert_eq!(s.histogram[0], 1); // -20.0
        assert_eq!(s.histogram[39], 1); // -0.25 in [-0.5, 0)
        assert_eq!(s.histogram[40], 2); // 0.0 and 0.49 in [0, 0.5)
        assert_eq!(s.histogram[119], 1); // 39.99
    }

    #[test]
    fn damage_counts_use_molecule_orientation() {
        let mut d = DamageCounts::new(3);
        let refs = |s: &str| s.bytes().map(Some).collect::<Vec<_>>();
        // Forward: 5′ C→T at position 0. Reverse: SAM 3′ G→A is the molecule's 5′ C→T.
        d.add_read(b"TAAA", &[40; 4], &refs("CAAA"), false, &settings());
        d.add_read(b"AAAA", &[40; 4], &refs("AAAG"), true, &settings());
        assert_eq!(d.five_prime[0][0], [0, 0, 0, 2]); // two ref C at 5′ terminal, both read T
        let (rate, se, n) = d.first_position();
        assert_eq!((rate, se, n), (Some(1.0), Some(0.0), 2));
    }
}
