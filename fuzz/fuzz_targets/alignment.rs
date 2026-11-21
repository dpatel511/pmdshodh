//! Arbitrary SAM records (after a fixed header) through the per-read pipeline, in every mode, with
//! MD and with a reference window. Expect an outcome with a finite score, or an error; never a
//! panic.
#![no_main]

use std::io::BufReader;

use libfuzzer_sys::fuzz_target;
use noodles::sam;
use pmdshodh::damage::DamageCounts;
use pmdshodh::model::Params;
use pmdshodh::reads::{Outcome, ReadFilters, process};
use pmdshodh::sites::{Library, Settings, Treatment};

const HEADER: &[u8] = b"@HD\tVN:1.6\n@SQ\tSN:c\tLN:1000\n";
const REFERENCE: &[u8] = b"ACGTCGNACGTTGCAcgtacgNNCGCGTATA";

fn check(outcome: Result<Outcome, String>) {
    if let Ok(Outcome::Scored { score, .. }) = outcome {
        assert!(score.pmds.is_finite(), "non-finite score {}", score.pmds);
    }
}

fuzz_target!(|data: &[u8]| {
    let Some((&mode, body)) = data.split_first() else {
        return;
    };
    let settings = Settings {
        library: if mode & 1 == 0 {
            Library::DoubleStranded
        } else {
            Library::SingleStranded
        },
        treatment: match (mode >> 1) % 3 {
            0 => Treatment::None,
            1 => Treatment::Half,
            _ => Treatment::Full,
        },
        params: Params::default(),
        min_baseq: if mode & 8 == 0 { 0 } else { 20 },
        pmdtools_compat: mode & 16 != 0,
    };
    let filters = ReadFilters {
        min_mapq: if mode & 32 == 0 { 0 } else { 30 },
        no_clips: mode & 64 != 0,
        no_indels: mode & 128 != 0,
    };
    let text = [HEADER, body].concat();
    let mut reader = sam::io::Reader::new(BufReader::new(&text[..]));
    if reader.read_header().is_err() {
        return;
    }
    let mut damage = DamageCounts::new(30);
    for record in reader.records().take(64) {
        let Ok(record) = record else {
            break;
        };
        check(process(&record, &settings, &filters, None, Some(&mut damage)));
        // With -r: a reference window as long as the alignment's reference span.
        let span: usize = record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        let window: Vec<u8> = REFERENCE.iter().copied().cycle().take(span).collect();
        check(process(&record, &settings, &filters, Some(&window), Some(&mut damage)));
    }
});
