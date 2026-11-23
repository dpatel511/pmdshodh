//! Arbitrary bytes as a CRAM file, read as `score -r` reads CRAM: the header is checked against the
//! reference (every sequence present, same length), containers are read in batches, their
//! reference sequences are staged with `Reference::cram_batch`, and the pool decoder cuts each
//! record's reference window. Each record is then scored, in default and compat mode, with its
//! window. The reference is tests/data/pairs.fa and past_end.fa in one indexed FASTA, the
//! references of the seed CRAMs. Expect an outcome with a finite score, or an error; never a
//! panic. Added after audit finding B9 (docs/audit-0.3.md, O4), which no other target could reach.
#![no_main]

use std::io::Cursor;
use std::path::PathBuf;
use std::sync::OnceLock;

use libfuzzer_sys::fuzz_target;
use pmdshodh::cram::{Reader, decode, reference_sequence_id};
use pmdshodh::model::Params;
use pmdshodh::reads::{Outcome, ReadFilters, process};
use pmdshodh::reference::Reference;
use pmdshodh::sites::{Library, Settings, Treatment};

const FASTAS: [&str; 2] = [
    include_str!("../../tests/data/pairs.fa"),
    include_str!("../../tests/data/past_end.fa"),
];

/// The fixture sequences in one FASTA, one line each, with its `.fai`, written once.
fn fasta() -> &'static PathBuf {
    static PATH: OnceLock<PathBuf> = OnceLock::new();
    PATH.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("pmdshodh-fuzz-cram-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (mut fa, mut fai) = (String::new(), String::new());
        for text in FASTAS {
            for entry in text.split('>').filter(|e| !e.trim().is_empty()) {
                let (name_line, seq) = entry.split_once('\n').unwrap();
                let name = name_line.split_whitespace().next().unwrap();
                let seq: String = seq.split_whitespace().collect();
                fa.push_str(&format!(">{name}\n"));
                let offset = fa.len();
                fa.push_str(&format!("{seq}\n"));
                let n = seq.len();
                fai.push_str(&format!("{name}\t{n}\t{offset}\t{n}\t{}\n", n + 1));
            }
        }
        let path = dir.join("ref.fa");
        std::fs::write(&path, fa).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), fai).unwrap();
        path
    })
}

fn check(outcome: Result<Outcome, String>) {
    if let Ok(Outcome::Scored { score, .. }) = outcome {
        assert!(score.pmds.is_finite(), "non-finite score {}", score.pmds);
    }
}

fuzz_target!(|data: &[u8]| {
    let mut reader = Reader::new(Box::new(Cursor::new(data.to_vec())));
    let Ok(header) = reader.read_header() else {
        return;
    };
    let path = fasta();
    let mut reference = Reference::open(path).unwrap();
    // main.rs `check_reference`: a header sequence missing from the reference, or of another
    // length, stops the run before any record is read.
    if !header
        .reference_sequences()
        .iter()
        .all(|(name, seq)| reference.length(&name.to_string()) == Some(usize::from(seq.length())))
    {
        return;
    }
    if reference.enable_cram(path, &header).is_err() {
        return;
    }
    let filters = ReadFilters {
        min_mapq: 0,
        no_clips: false,
        no_indels: false,
    };
    for _ in 0..16 {
        let Ok(containers) = reader.read_containers(4) else {
            return;
        };
        if containers.is_empty() {
            return;
        }
        let ids: Vec<usize> = containers
            .iter()
            .filter_map(reference_sequence_id)
            .collect();
        let Ok((repository, _md5_checks)) = reference.cram_batch(&header, &ids) else {
            return;
        };
        let Ok(decoded) = decode(&containers, &repository, &header) else {
            return;
        };
        for (record, window) in decoded.records.iter().zip(&decoded.windows) {
            let window = window.clone().map(|r| &decoded.bases[r]);
            for pmdtools_compat in [false, true] {
                let settings = Settings {
                    library: Library::DoubleStranded,
                    treatment: Treatment::Half,
                    params: Params::default(),
                    min_baseq: 0,
                    pmdtools_compat,
                };
                check(process(record, &settings, &filters, window, None));
            }
        }
    }
});
