//! MD and CIGAR: arbitrary read bases, CIGAR operations and MD tags through the reconstruction
//! (an error or a result, never a panic); and for valid alignments, the reference rebuilt from
//! the MD that `md_from_reference` writes equals the reference taken from the window directly.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pmdshodh::sites::{
    CigarOp, md_from_reference, reference_along_read, reference_along_read_from_fasta,
    reference_span,
};

fn op(kind: u8, n: usize) -> CigarOp {
    match kind % 8 {
        0 => CigarOp::Match(n),
        1 => CigarOp::SequenceMatch(n),
        2 => CigarOp::SequenceMismatch(n),
        3 => CigarOp::Insertion(n),
        4 => CigarOp::Deletion(n),
        5 => CigarOp::SoftClip(n),
        6 => CigarOp::HardClip(n),
        _ => CigarOp::Unsupported,
    }
}

fuzz_target!(|data: &[u8]| {
    // Layout: number of operations, (kind, length) per operation, a split byte, then read bases
    // and the MD tag or reference window.
    let Some((&n_ops, rest)) = data.split_first() else {
        return;
    };
    let n_ops = (n_ops % 8) as usize;
    if rest.len() < 2 * n_ops + 1 {
        return;
    }
    let (ops, rest) = rest.split_at(2 * n_ops);
    let cigar: Vec<CigarOp> = ops
        .chunks(2)
        .map(|c| op(c[0], (c[1] % 16) as usize))
        .collect();
    let (&split, rest) = rest.split_first().unwrap();
    let (seq, tail) = rest.split_at((split as usize).min(rest.len()));

    // Arbitrary MD.
    let _ = reference_along_read(seq, &cigar, tail);

    // Valid alignments: read length matches the CIGAR, bases from ACGTN, reference from the window.
    let read_len: usize = cigar
        .iter()
        .map(|o| match o {
            CigarOp::Match(n)
            | CigarOp::SequenceMatch(n)
            | CigarOp::SequenceMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n) => *n,
            _ => 0,
        })
        .sum();
    if cigar.contains(&CigarOp::Unsupported) || seq.len() != read_len {
        return;
    }
    let seq: Vec<u8> = seq.iter().map(|b| b"ACGTN"[(*b % 5) as usize]).collect();
    let span = reference_span(&cigar);
    if tail.len() < span {
        return;
    }
    let window: Vec<u8> = tail[..span]
        .iter()
        .map(|b| b"ACGTNacgtn"[(*b % 10) as usize])
        .collect();
    let md = md_from_reference(&seq, &cigar, &window);
    let from_md = reference_along_read(&seq, &cigar, &md);
    let from_fasta = reference_along_read_from_fasta(&seq, &cigar, &window);
    match (&from_md, &from_fasta) {
        (Ok(a), Ok(b)) => assert_eq!(a, b, "MD {}", String::from_utf8_lossy(&md)),
        _ => assert_eq!(
            from_md.is_ok(),
            from_fasta.is_ok(),
            "{from_md:?} {from_fasta:?}"
        ),
    }
});
