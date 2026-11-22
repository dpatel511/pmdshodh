//! Classify one alignment record: skip it (with a reason), reject it (invalid data) or score it.

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::damage::DamageCounts;
use crate::sites::{
    CigarOp, ReadScore, Settings, md_from_reference, reference_along_read,
    reference_along_read_from_fasta, reference_along_read_pmdtools, score_read,
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
    /// `--pmdtools-compat` only: the rebuilt reference has no A, C, G or T (audit finding B5).
    NoReferenceBases,
}

impl SkipReason {
    pub const ALL: [SkipReason; 8] = [
        SkipReason::Unmapped,
        SkipReason::LowMapq,
        SkipReason::Clipped,
        SkipReason::Indel,
        SkipReason::UnsupportedCigar,
        SkipReason::NoSequence,
        SkipReason::NoQualities,
        SkipReason::NoReferenceBases,
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
            SkipReason::NoReferenceBases => "no_reference_bases",
        }
    }
}

/// How a scored read's MD tag compared with the reference given with `-r`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MdCheck {
    /// No `-r`: reference bases came from MD.
    NotChecked,
    /// `-r` given; the read has no MD tag.
    Absent,
    /// `-r` given; MD gives the same reference bases.
    Agrees,
    /// `-r` given; MD is invalid or gives other reference bases (written before realignment, or
    /// against reference N bases that the aligner replaced, as bwa does).
    Disagrees,
}

/// Result of looking at one read.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Outcome {
    /// Scored; `length` is the SEQ length.
    Scored {
        score: ReadScore,
        length: usize,
        md: MdCheck,
    },
    Skipped(SkipReason),
}

/// Look at one record. `Err` means the record is invalid input (e.g. a mapped read without MD).
/// Scored reads are also added to `damage`, if given.
pub fn process(
    record: &dyn Record,
    settings: &Settings,
    filters: &ReadFilters,
    reference: Option<&[u8]>,
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
    process_fields(&fields, settings, filters, reference, damage)
}

/// [`process`] for a BAM record, reading its fields directly instead of through the generic
/// record trait (no boxed per-base iterators). Same outcomes and messages.
pub fn process_bam(
    record: &noodles::bam::Record,
    settings: &Settings,
    filters: &ReadFilters,
    reference: Option<&[u8]>,
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
    process_fields(&fields, settings, filters, reference, damage)
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

/// `reference`: with `-r`, the reference bases from the alignment start, `reference_span` long.
///
/// Unmapped reads (flag 4) are skipped. PMDtools ignores flags, so with `--pmdtools-compat` a read
/// flagged unmapped that still has an alignment (bwa writes CIGAR and MD for reads overhanging a
/// reference end) is scored from its MD as PMDtools scores it; one with no alignment, or one that
/// cannot be scored, still counts as unmapped (audit finding B4). No reference window is fetched
/// for unmapped reads, since they may overhang the reference.
fn process_fields(
    record: &Fields<'_>,
    settings: &Settings,
    filters: &ReadFilters,
    reference: Option<&[u8]>,
    damage: Option<&mut DamageCounts>,
) -> Result<Outcome, String> {
    let flags = (record.flags)().map_err(|e| format!("invalid FLAG: {e}"))?;
    if !flags.is_unmapped() {
        return score_fields(record, settings, filters, reference, damage);
    }
    let placed = settings.pmdtools_compat && (record.cigar)().is_ok_and(|c| !c.is_empty());
    if !placed {
        return Ok(Outcome::Skipped(SkipReason::Unmapped));
    }
    Ok(score_fields(record, settings, filters, None, damage)
        .unwrap_or(Outcome::Skipped(SkipReason::Unmapped)))
}

fn score_fields(
    record: &Fields<'_>,
    settings: &Settings,
    filters: &ReadFilters,
    reference: Option<&[u8]>,
    damage: Option<&mut DamageCounts>,
) -> Result<Outcome, String> {
    let flags = (record.flags)().map_err(|e| format!("invalid FLAG: {e}"))?;
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
    // PMDtools drops every read with a hard clip, whatever its flags, with the message it gives
    // for `N` and `P` (docs/blackbox.md 1.6f, audit finding B6). Default mode scores them (D3).
    if settings.pmdtools_compat && has(|op| matches!(op, CigarOp::HardClip(_))) {
        return Ok(Outcome::Skipped(SkipReason::UnsupportedCigar));
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

    let reverse = flags.is_reverse_complemented();
    let (refs, md_check) = match reference {
        // `-r`: the reference is the truth; MD, if present, is only compared with it.
        Some(window) => {
            let refs = reference_along_read_from_fasta(&seq, &cigar, window)?;
            let md_check = match (record.md)() {
                None => MdCheck::Absent,
                Some(Ok(Value::String(md))) => {
                    match reference_along_read(&seq, &cigar, md.as_ref()) {
                        // MD can only give a `=` read base back as `=`.
                        Ok(from_md)
                            if from_md
                                .iter()
                                .zip(&refs)
                                .zip(&seq)
                                .all(|((m, r), &b)| m == r || b == b'=') =>
                        {
                            MdCheck::Agrees
                        }
                        _ => MdCheck::Disagrees,
                    }
                }
                Some(_) => MdCheck::Disagrees,
            };
            let refs = if settings.pmdtools_compat {
                // As `samtools calmd | pmdtools`: PMDtools' reconstruction of a calmd MD tag.
                let md = md_from_reference(&seq, &cigar, window);
                reference_along_read_pmdtools(&seq, &cigar, &md, reverse)
            } else {
                refs
            };
            (refs, md_check)
        }
        None => {
            let md = match (record.md)() {
                None => {
                    return Err(
                        "mapped read has no MD tag (pass the reference with -r, or add MD with samtools calmd)"
                            .into(),
                    );
                }
                Some(Ok(Value::String(md))) => md.to_vec(),
                Some(Ok(_)) => return Err("MD tag is not a string".into()),
                Some(Err(e)) => return Err(format!("invalid MD tag: {e}")),
            };
            // Always validate with the correct reconstruction, so invalid records are rejected in both modes.
            let mut refs = reference_along_read(&seq, &cigar, &md)?;
            if settings.pmdtools_compat {
                refs = reference_along_read_pmdtools(&seq, &cigar, &md, reverse);
            }
            (refs, MdCheck::NotChecked)
        }
    };
    // PMDtools drops a read whose rebuilt reference has no uppercase A, C, G or T ("bad reference
    // sequence reconstruction"): all N, all `=` read bases, or lowercase MD letters (black-box
    // probes, docs/blackbox.md 1.6e). Default mode scores such a read: it has no sites.
    if settings.pmdtools_compat
        && !refs
            .iter()
            .any(|r| matches!(r, Some(b'A' | b'C' | b'G' | b'T')))
    {
        return Ok(Outcome::Skipped(SkipReason::NoReferenceBases));
    }
    // `=` read bases mean "same as the reference" (`samtools calmd -e` writes them). They are
    // resolved from `-r` for scoring; from MD alone their reference base is unknown. PMDtools
    // never counts them as C or G, and `--pmdtools-compat` keeps that. Kept records are written
    // unchanged.
    let seq = if settings.pmdtools_compat || !seq.contains(&b'=') {
        seq
    } else if reference.is_some() {
        seq.iter()
            .zip(&refs)
            .map(|(&b, r)| match (b, r) {
                (b'=', Some(r)) => *r,
                _ => b,
            })
            .collect()
    } else {
        return Err(
            "read has '=' bases, whose reference bases MD does not give (pass the reference with -r)"
                .into(),
        );
    };
    if let Some(d) = damage {
        d.add_read(&seq, &qual, &refs, reverse, settings);
    }
    Ok(Outcome::Scored {
        score: score_read(&seq, &qual, &refs, reverse, settings),
        length: seq.len(),
        md: md_check,
    })
}
