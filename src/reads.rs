//! Classify one alignment record: skip it (with a reason), reject it (invalid data) or score it.

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::damage::DamageCounts;
use crate::sites::{
    CigarOp, ReadScore, Settings, reference_along_read, reference_along_read_pmdtools, score_read,
};

/// Read-level filters (PMDtools `-m`, `--noclips`, `--noindels`).
#[derive(Debug, Clone, Copy, Default)]
pub struct ReadFilters {
    pub min_mapq: u8,
    pub no_clips: bool,
    pub no_indels: bool,
}

/// Why a read was not scored. Every skipped read is counted by reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SkipReason {
    Unmapped,
    LowMapq,
    Clipped,
    Indel,
    UnsupportedCigar,
    NoSequence,
    NoQualities,
}

impl SkipReason {
    pub const ALL: [SkipReason; 7] = [
        SkipReason::Unmapped,
        SkipReason::LowMapq,
        SkipReason::Clipped,
        SkipReason::Indel,
        SkipReason::UnsupportedCigar,
        SkipReason::NoSequence,
        SkipReason::NoQualities,
    ];

    pub fn name(self) -> &'static str {
        match self {
            SkipReason::Unmapped => "unmapped",
            SkipReason::LowMapq => "min_mapq",
            SkipReason::Clipped => "clips",
            SkipReason::Indel => "indels",
            SkipReason::UnsupportedCigar => "unsupported_cigar",
            SkipReason::NoSequence => "no_sequence",
            SkipReason::NoQualities => "no_qualities",
        }
    }
}

/// Result of looking at one read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Outcome {
    /// Scored; `length` is the SEQ length.
    Scored {
        score: ReadScore,
        length: usize,
    },
    Skipped(SkipReason),
}

/// Look at one record. `Err` means the record is invalid input (e.g. a mapped read without MD).
/// Scored reads are also added to `damage`, if given.
pub fn process(
    record: &dyn Record,
    settings: &Settings,
    filters: &ReadFilters,
    damage: Option<&mut DamageCounts>,
) -> Result<Outcome, String> {
    let fields = Fields {
        flags: &|| record.flags(),
        mapping_quality: &|| record.mapping_quality().transpose(),
        cigar: &|| record.cigar().iter().collect(),
        sequence: &|| record.sequence().iter().collect(),
        quality_scores: &|| record.quality_scores().iter().collect(),
        md: &|| record.data().get(&Tag::MISMATCHED_POSITIONS),
    };
    process_fields(&fields, settings, filters, damage)
}

/// [`process`] for a BAM record, reading its fields directly instead of through the generic
/// record trait (no boxed per-base iterators). Same outcomes and messages.
pub fn process_bam(
    record: &noodles::bam::Record,
    settings: &Settings,
    filters: &ReadFilters,
    damage: Option<&mut DamageCounts>,
) -> Result<Outcome, String> {
    let fields = Fields {
        flags: &|| Ok(record.flags()),
        mapping_quality: &|| Ok(record.mapping_quality()),
        cigar: &|| record.cigar().iter().collect(),
        sequence: &|| record.sequence().iter().collect(),
        quality_scores: &|| Ok(record.quality_scores().as_ref().to_vec()),
        md: &|| record.data().get(&Tag::MISMATCHED_POSITIONS),
    };
    process_fields(&fields, settings, filters, damage)
}

type Field<'a, T> = &'a dyn Fn() -> T;

/// Field accessors, called only when needed, so a skipped read is not decoded further.
struct Fields<'a> {
    flags: Field<'a, std::io::Result<noodles::sam::alignment::record::Flags>>,
    mapping_quality:
        Field<'a, std::io::Result<Option<noodles::sam::alignment::record::MappingQuality>>>,
    cigar: Field<'a, std::io::Result<Vec<noodles::sam::alignment::record::cigar::Op>>>,
    sequence: Field<'a, Vec<u8>>,
    quality_scores: Field<'a, std::io::Result<Vec<u8>>>,
    md: Field<'a, Option<std::io::Result<Value<'a>>>>,
}

fn process_fields(
    record: &Fields<'_>,
    settings: &Settings,
    filters: &ReadFilters,
    damage: Option<&mut DamageCounts>,
) -> Result<Outcome, String> {
    let flags = (record.flags)().map_err(|e| format!("invalid FLAG: {e}"))?;
    if flags.is_unmapped() {
        return Ok(Outcome::Skipped(SkipReason::Unmapped));
    }
    // A missing mapping quality (255) is kept by any threshold, as in PMDtools.
    let mapq = match (record.mapping_quality)().map_err(|e| format!("invalid MAPQ: {e}"))? {
        Some(m) => m.get(),
        None => 255,
    };
    if mapq < filters.min_mapq {
        return Ok(Outcome::Skipped(SkipReason::LowMapq));
    }

    let mut cigar = Vec::new();
    for op in (record.cigar)().map_err(|e| format!("invalid CIGAR: {e}"))? {
        let n = op.len();
        cigar.push(match op.kind() {
            Kind::Match => CigarOp::Match(n),
            Kind::SequenceMatch => CigarOp::SequenceMatch(n),
            Kind::SequenceMismatch => CigarOp::SequenceMismatch(n),
            Kind::Insertion => CigarOp::Insertion(n),
            Kind::Deletion => CigarOp::Deletion(n),
            Kind::SoftClip => CigarOp::SoftClip(n),
            Kind::HardClip => CigarOp::HardClip(n),
            Kind::Skip | Kind::Pad => CigarOp::Unsupported,
        });
    }
    if cigar.contains(&CigarOp::Unsupported) {
        return Ok(Outcome::Skipped(SkipReason::UnsupportedCigar));
    }
    let has = |f: fn(&CigarOp) -> bool| cigar.iter().any(f);
    if filters.no_clips && has(|op| matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_))) {
        return Ok(Outcome::Skipped(SkipReason::Clipped));
    }
    if filters.no_indels && has(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_))) {
        return Ok(Outcome::Skipped(SkipReason::Indel));
    }

    let seq = (record.sequence)();
    if seq.is_empty() {
        return Ok(Outcome::Skipped(SkipReason::NoSequence));
    }
    let qual = (record.quality_scores)().map_err(|e| format!("invalid QUAL: {e}"))?;
    if qual.is_empty() || qual.iter().all(|&q| q == 255) {
        return Ok(Outcome::Skipped(SkipReason::NoQualities));
    }
    if qual.len() != seq.len() {
        return Err(format!(
            "QUAL has {} values but SEQ has {} bases",
            qual.len(),
            seq.len()
        ));
    }

    let md = match (record.md)() {
        None => return Err("mapped read has no MD tag (add it with samtools calmd)".into()),
        Some(Ok(Value::String(md))) => md.to_vec(),
        Some(Ok(_)) => return Err("MD tag is not a string".into()),
        Some(Err(e)) => return Err(format!("invalid MD tag: {e}")),
    };
    // Always validate with the correct reconstruction, so invalid records are rejected in both modes.
    let mut refs = reference_along_read(&seq, &cigar, &md)?;
    let reverse = flags.is_reverse_complemented();
    if settings.pmdtools_compat {
        refs = reference_along_read_pmdtools(&seq, &cigar, &md, reverse);
    }
    if let Some(d) = damage {
        d.add_read(&seq, &qual, &refs, reverse, settings);
    }
    Ok(Outcome::Scored {
        score: score_read(&seq, &qual, &refs, reverse, settings),
        length: seq.len(),
    })
}
