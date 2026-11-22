//! Edge cases from the 0.3 audit coverage map (docs/audit-0.3.md §3.1, group a): read filters and
//! skip reasons, full UDG from the command line, `=`/`X` CIGARs, malformed MD tags, reference
//! checks, damaged BGZF input and command-line errors. Also the audit's findings, and the tests
//! its edge-case matrix lacked (§4, mapped in §4.5).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_pmdshodh");

fn tmp_dir(test: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("edge_{test}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().unwrap()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

/// A SAM file with one reference `c` of 1,000 bases and the given records (tab-separated lines).
fn sam(dir: &Path, name: &str, records: &[&str]) -> PathBuf {
    let mut text = String::from("@HD\tVN:1.6\n@SQ\tSN:c\tLN:1000\n");
    for r in records {
        text.push_str(r);
        text.push('\n');
    }
    let path = dir.join(name);
    fs::write(&path, text).unwrap();
    path
}

/// The JSON summary of `damage` with extra options.
fn damage_json(input: &Path, extra: &[&str]) -> serde_json::Value {
    let mut args = vec!["damage", "--quiet"];
    args.extend_from_slice(extra);
    args.push(s(input));
    let o = run(&args);
    assert!(o.status.success(), "{extra:?}: {}", stderr(&o));
    serde_json::from_slice(&o.stdout).unwrap()
}

/// An input error (exit 2) whose message contains `says`.
fn input_error(args: &[&str], says: &str) {
    let o = run(args);
    assert_eq!(o.status.code(), Some(2), "{args:?}: {}", stderr(&o));
    assert!(stderr(&o).contains(says), "{args:?}: {}", stderr(&o));
}

const GOOD: &str = "good\t0\tc\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9";

/// One read per skip reason. Each is counted once, under its own reason, and only when its
/// option is given; `N` and `P` CIGARs, `*` SEQ and `*` QUAL are always skipped.
#[test]
fn every_skip_reason_is_counted_once() {
    let dir = tmp_dir("skips");
    let input = sam(
        &dir,
        "skips.sam",
        &[
            GOOD,
            "lowmapq\t0\tc\t11\t10\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            // MAPQ 255 means "not available": it passes any --min-mapq, as in PMDtools.
            "mapq255\t0\tc\t11\t255\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            "softclip\t0\tc\t11\t60\t2S8M\t*\t0\t0\tGGTAAAAAAA\tIIIIIIIIII\tMD:Z:0C7",
            "hardclip\t0\tc\t11\t60\t2H10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            "insertion\t0\tc\t11\t60\t5M1I4M\t*\t0\t0\tTAAAAGAAAA\tIIIIIIIIII\tMD:Z:0C8",
            "deletion\t0\tc\t11\t60\t5M1D5M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C4^A5",
            "intron\t0\tc\t11\t60\t5M100N5M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            "padding\t0\tc\t11\t60\t5M1P5M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            "noseq\t0\tc\t11\t60\t10M\t*\t0\t0\t*\t*\tMD:Z:0C9",
            "noqual\t0\tc\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\t*\tMD:Z:0C9",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII",
        ],
    );
    let skipped = |j: &serde_json::Value, reason: &str| j["reads"]["skipped"][reason].as_u64();
    let with = damage_json(&input, &["--min-mapq", "20", "--no-clips", "--no-indels"]);
    assert_eq!(with["reads"]["input"], 12);
    for (reason, n) in [
        ("min_mapq", 1),
        ("clips", 2),
        ("indels", 2),
        ("unsupported_cigar", 2),
        ("no_sequence", 1),
        ("no_qualities", 1),
        ("unmapped", 1),
    ] {
        assert_eq!(skipped(&with, reason), Some(n), "{reason}");
    }
    assert_eq!(with["reads"]["scored"], 2, "good and mapq255");
    // Without the options, clipped, indel and low-MAPQ reads are scored.
    let without = damage_json(&input, &[]);
    for reason in ["min_mapq", "clips", "indels"] {
        assert_eq!(skipped(&without, reason), Some(0), "{reason}");
    }
    assert_eq!(without["reads"]["scored"], 7);
    // PMDtools' option names do the same.
    let aliases = damage_json(&input, &["--requiremapq", "20", "--noclips", "--noindels"]);
    assert_eq!(aliases["reads"], with["reads"]);
    // The stderr line names each reason with a count (without --quiet).
    let o = run(&[
        "damage",
        "--min-mapq",
        "20",
        s(&input),
        "-o",
        s(&dir.join("d.json")),
    ]);
    assert!(o.status.success());
    assert!(
        stderr(&o).contains(
            "skipped: unmapped=1, min_mapq=1, unsupported_cigar=2, no_sequence=1, no_qualities=1"
        ),
        "{}",
        stderr(&o)
    );
}

/// A QUAL whose length differs from SEQ, and an MD tag that is not a string, are input errors.
#[test]
fn bad_qual_length_and_non_string_md_are_input_errors() {
    let dir = tmp_dir("bad_fields");
    let short_qual = sam(
        &dir,
        "q.sam",
        &["r1\t0\tc\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIII\tMD:Z:0C9"],
    );
    let o = run(&["score", "--quiet", s(&short_qual)]);
    assert_eq!(o.status.code(), Some(2), "{}", stderr(&o));
    assert!(stderr(&o).contains("r1"), "{}", stderr(&o));
    let int_md = sam(
        &dir,
        "md.sam",
        &["r2\t0\tc\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:i:10"],
    );
    input_error(&["score", "--quiet", s(&int_md)], "MD tag is not a string");
}

/// `--treatment full`, `--UDGplus` and `--CpG` are the same mode: only CpG sites count, so the
/// scores and damage tables differ from no treatment, and the summary records "full".
#[test]
fn full_udg_options_are_equivalent_and_count_only_cpg() {
    let input = data("reads_random.sam");
    let full = damage_json(&input, &["--treatment", "full"]);
    assert_eq!(full["mode"]["treatment"], "full");
    let without_command = |mut j: serde_json::Value| {
        j.as_object_mut().unwrap().remove("command");
        j
    };
    for alias in ["--UDGplus", "--CpG"] {
        let j = damage_json(&input, &[alias]);
        assert_eq!(without_command(j), without_command(full.clone()), "{alias}");
    }
    let none = damage_json(&input, &[]);
    let sites = |j: &serde_json::Value| -> u64 {
        ["A", "C", "G", "T"]
            .iter()
            .map(|b| {
                j["damage"]["five_prime"]["C"][b]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_u64().unwrap())
                    .sum::<u64>()
            })
            .sum()
    };
    assert!(
        sites(&full) > 0 && sites(&full) < sites(&none),
        "fewer sites with full UDG"
    );
    let score = |extra: &[&str]| {
        let mut args = vec!["score", "--quiet"];
        args.extend_from_slice(extra);
        args.push(s(&input));
        run(&args).stdout
    };
    assert_eq!(score(&["--treatment", "full"]), score(&["--CpG"]));
    assert_ne!(score(&["--treatment", "full"]), score(&[]));
}

/// `=` and `X` CIGAR operations align bases like `M`: a read written with them scores the same.
#[test]
fn sequence_match_and_mismatch_cigar_ops_score_like_m() {
    let dir = tmp_dir("eqx");
    let m = sam(&dir, "m.sam", &[GOOD]);
    let eqx = sam(
        &dir,
        "eqx.sam",
        &["good\t0\tc\t11\t60\t1X9=\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9"],
    );
    let score = |p: &Path| run(&["score", "--quiet", s(p)]).stdout;
    assert_eq!(score(&eqx), score(&m));
    assert_eq!(String::from_utf8(score(&m)).unwrap().lines().count(), 2);
    // With -r, `=` and `X` count towards the reference window fetched for the read.
    let fasta = dir.join("c.fa");
    fs::write(
        &fasta,
        format!(">c\n{}C{}\n", "A".repeat(10), "A".repeat(989)),
    )
    .unwrap();
    fs::write(dir.join("c.fa.fai"), "c\t1000\t3\t1000\t1001\n").unwrap();
    let with_r = |p: &Path| run(&["score", "--quiet", "-r", s(&fasta), s(p)]).stdout;
    assert_eq!(with_r(&eqx), score(&m));
    assert_eq!(with_r(&m), score(&m));
}

/// MD tags that cannot describe their CIGAR are input errors naming the problem.
#[test]
fn malformed_md_tags_are_input_errors() {
    let dir = tmp_dir("bad_md");
    for (i, (cigar, md, says)) in [
        ("10M", "5!4", "MD contains an invalid character"),
        (
            "5M1D5M",
            "11",
            "CIGAR has a deletion inside an MD match run",
        ),
        (
            "5M1D5M",
            "4A5",
            "CIGAR has a deletion where MD has matching bases",
        ),
        (
            "5M1D5M",
            "5A4",
            "CIGAR has a deletion that MD does not have",
        ),
        (
            "5M2D5M",
            "5^A5",
            "MD deletion is shorter than the CIGAR deletion",
        ),
    ]
    .iter()
    .enumerate()
    {
        let record =
            format!("r{i}\t0\tc\t11\t60\t{cigar}\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:{md}");
        let input = sam(&dir, &format!("md{i}.sam"), &[record.as_str()]);
        input_error(&["score", "--quiet", s(&input)], says);
    }
}

/// With -r, a header sequence missing from the FASTA is an input error. (An empty FASTA sequence
/// cannot occur: samtools faidx leaves it out of the index, and SAM requires LN >= 1. A FASTA whose
/// line lengths differ from its index is audit finding B1, tested with its fix.)
#[test]
fn reference_problems_are_input_errors() {
    let dir = tmp_dir("reference");
    let fasta = data("pairs.fa");
    // c3 is in the header but not in pairs.fa.
    let reads = dir.join("missing.sam");
    fs::write(
        &reads,
        "@SQ\tSN:c1\tLN:300\n@SQ\tSN:c3\tLN:300\nr1\t0\tc1\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\n",
    )
    .unwrap();
    input_error(
        &["score", "--quiet", "-r", s(&fasta), s(&reads)],
        "c3 is not in the reference",
    );
}

/// A BAM written by pmdshodh, for damaging.
fn bam(dir: &Path) -> Vec<u8> {
    let out = dir.join("in.bam");
    let o = run(&[
        "filter",
        "--quiet",
        "-t",
        "-1000000",
        s(&data("reads_random.sam")),
        "-o",
        s(&out),
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    fs::read(out).unwrap()
}

/// Damaged BGZF blocks are input errors that say what is wrong; a BGZF file with only the
/// end-of-file block is empty input.
#[test]
fn damaged_bgzf_input_is_named() {
    let dir = tmp_dir("bgzf");
    let good = bam(&dir);
    let size = usize::from(u16::from_le_bytes([good[16], good[17]])) + 1;
    let mut cases: Vec<(&str, Vec<u8>, &str)> = Vec::new();
    let mut small = good.clone();
    small[16..18].copy_from_slice(&10u16.to_le_bytes());
    cases.push(("small.bam", small, "invalid BGZF block size"));
    let mut isize = good.clone();
    isize[size - 4..size].copy_from_slice(&70_000u32.to_le_bytes());
    cases.push(("isize.bam", isize, "invalid BGZF ISIZE"));
    let mut data_bad = good.clone();
    for b in &mut data_bad[18..40] {
        *b = 0xff;
    }
    cases.push(("deflate.bam", data_bad, "invalid BGZF block data"));
    const EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    cases.push(("eof_only.bam", EOF.to_vec(), "is empty (no data)"));
    for (name, bytes, says) in cases {
        let path = dir.join(name);
        fs::write(&path, bytes).unwrap();
        input_error(&["score", "--quiet", s(&path)], says);
    }
}

/// A 70,000-base read with random bases and qualities: its BAM record is larger than a BGZF
/// block and does not compress, and it is written and read back unchanged.
#[test]
fn incompressible_record_larger_than_a_block_round_trips() {
    let dir = tmp_dir("big_record");
    let mut x: u64 = 0x9E3779B97F4A7C15;
    let mut rand = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let n = 70_000;
    let seq: String = (0..n)
        .map(|_| b"ACGT"[(rand() % 4) as usize] as char)
        .collect();
    let qual: String = (0..n)
        .map(|_| (b'!' + (rand() % 41) as u8) as char)
        .collect();
    let input = dir.join("big.sam");
    fs::write(
        &input,
        format!("@SQ\tSN:c\tLN:100000\nbig\t0\tc\t1\t60\t{n}M\t*\t0\t0\t{seq}\t{qual}\tMD:Z:{n}\n"),
    )
    .unwrap();
    let out = dir.join("big.bam");
    let o = run(&[
        "filter",
        "--quiet",
        "-t",
        "-1000000",
        s(&input),
        "-o",
        s(&out),
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    let back = dir.join("back.sam");
    let o = run(&[
        "filter",
        "--quiet",
        "-t",
        "-1000000",
        s(&out),
        "-o",
        s(&back),
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    let line = fs::read_to_string(&back).unwrap();
    let record = line.lines().find(|l| l.starts_with("big\t")).unwrap();
    let fields: Vec<&str> = record.split('\t').collect();
    assert_eq!((fields[9], fields[10]), (seq.as_str(), qual.as_str()));
}

/// `--threads 0`, an output that is a directory, and an input error in `damage` are reported
/// with exit status 2.
#[test]
fn command_line_errors() {
    let dir = tmp_dir("cli");
    let input = data("reads_random.sam");
    input_error(
        &["score", "--threads", "0", s(&input)],
        "--threads must be at least 1",
    );
    input_error(&["score", s(&input), "-o", s(&dir)], "is a directory");
    let no_md = sam(
        &dir,
        "no_md.sam",
        &["r1\t0\tc\t11\t60\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII"],
    );
    input_error(&["damage", "--quiet", s(&no_md)], "r1");
}

/// A write failure to an existing device (written in place) is a runtime failure, exit 1.
#[cfg(target_os = "linux")]
#[test]
fn write_failure_exits_1() {
    let o = run(&[
        "score",
        "--quiet",
        s(&data("reads_random.sam")),
        "-o",
        "/dev/full",
    ]);
    assert_eq!(o.status.code(), Some(1), "{}", stderr(&o));
    assert!(
        stderr(&o).contains("cannot write /dev/full"),
        "{}",
        stderr(&o)
    );
}

/// Audit finding B3: with -r, a CIGAR covering far more bases than SEQ aborted while allocating
/// (2e16 bases) or used gigabytes before the length check. It is now an input error at once.
#[test]
fn absurd_cigar_lengths_with_reference_are_input_errors() {
    let dir = tmp_dir("absurd_cigar");
    let fasta = dir.join("c.fa");
    fs::write(&fasta, format!(">c\n{}\n", "ACGT".repeat(250))).unwrap();
    fs::write(dir.join("c.fa.fai"), "c\t1000\t3\t1000\t1001\n").unwrap();
    for cigar in ["1M20000000000000000I", "1M4000000000I"] {
        let input = sam(
            &dir,
            "absurd.sam",
            &[&format!("r1\t0\tc\t10\t37\t{cigar}\t*\t0\t0\tAC\tII")],
        );
        for compat in [&[][..], &["--pmdtools-compat"][..]] {
            let args = [
                &["score", "--quiet", "-r", s(&fasta)][..],
                compat,
                &[s(&input)],
            ]
            .concat();
            let start = std::time::Instant::now();
            input_error(&args, "but SEQ has 2");
            assert!(
                start.elapsed().as_secs() < 5,
                "{cigar}: took {:?}",
                start.elapsed()
            );
        }
    }
}

/// Audit finding B4: PMDtools ignores flags, so it scores a read flagged unmapped that still has
/// a CIGAR and MD (bwa writes these for reads overhanging a reference end). `--pmdtools-compat`
/// skipped them. Now compat scores them as if mapped (from MD, also with -r); default still skips
/// them, and so does compat for flag-4 reads with no alignment or an unusable one.
#[test]
fn compat_scores_placed_unmapped_reads_like_pmdtools() {
    let dir = tmp_dir("placed_unmapped");
    let mapped = "r1\t0\tc\t11\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9";
    let placed = mapped.replacen("\t0\tc", "\t4\tc", 1);
    let records = [
        placed.as_str(),
        // No alignment, and an alignment whose MD cannot describe its CIGAR.
        "r2\t4\t*\t0\t0\t*\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII",
        "r3\t4\tc\t11\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:3",
    ];
    let input = sam(&dir, "placed.sam", &records);
    let alone = sam(&dir, "mapped.sam", &[mapped]);
    let fasta = dir.join("c.fa");
    fs::write(
        &fasta,
        format!(">c\n{}C{}\n", "A".repeat(10), "A".repeat(989)),
    )
    .unwrap();
    fs::write(dir.join("c.fa.fai"), "c\t1000\t3\t1000\t1001\n").unwrap();
    let score = |args: &[&str], p: &Path| {
        let o = run(&[&["score", "--quiet"][..], args, &[s(p)]].concat());
        assert!(o.status.success(), "{args:?}: {}", stderr(&o));
        // Read name, score and sites; the flag column keeps each read's own flag.
        String::from_utf8(o.stdout)
            .unwrap()
            .lines()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                format!("{} {} {}", f[0], f[2], f[3])
            })
            .collect::<Vec<_>>()
            .join("\n")
    };
    let want = score(&["--pmdtools-compat"], &alone);
    for extra in [&[][..], &["-r", s(&fasta)][..]] {
        let compat = score(&[&["--pmdtools-compat"][..], extra].concat(), &input);
        assert_eq!(
            compat, want,
            "compat {extra:?}: placed unmapped r1 scored as mapped"
        );
        let j = damage_json(&input, &[&["--pmdtools-compat"][..], extra].concat());
        assert_eq!(j["reads"]["scored"], 1);
        assert_eq!(j["reads"]["skipped"]["unmapped"], 2);
    }
    let default = score(&[], &input);
    assert_eq!(
        default.lines().count(),
        1,
        "default skips all unmapped reads: {default}"
    );
    assert_eq!(damage_json(&input, &[])["reads"]["skipped"]["unmapped"], 3);
}

/// Audit finding B5: PMDtools drops a read whose rebuilt reference has no uppercase A, C, G or T
/// ("bad reference sequence reconstruction"): all N, all `=` read bases, lowercase MD letters. The
/// probes and outcomes below are PMDtools 0.60's (2026-11-22). `--pmdtools-compat` now drops them
/// too, counted as `no_reference_bases`; default mode scores them (0 sites, score 0).
#[test]
fn compat_drops_reads_without_reference_bases_like_pmdtools() {
    let dir = tmp_dir("no_reference_bases");
    let q = "IIIIIIIIII";
    let rec = |name: &str, cigar: &str, seq: &str, md: &str| {
        format!("{name}\t0\tc\t11\t37\t{cigar}\t*\t0\t0\t{seq}\t{q}\tMD:Z:{md}")
    };
    let dropped = [
        rec("all_n", "10M", "ACGTACGTAC", "0N0N0N0N0N0N0N0N0N0N0"),
        rec("n_reads", "10M", "NNNNNNNNNN", "10"),
        rec(
            "n_with_d",
            "5M1D5M",
            "ACGTACGTAC",
            "0N0N0N0N0N^N0N0N0N0N0N0",
        ),
        rec("n_clip", "2S8M", "ACGTACGTAC", "0N0N0N0N0N0N0N0N0"),
        rec("all_eq", "10M", "==========", "10"),
        rec("lower", "10M", "ACGTACGTAC", "0n0n0n0n0n0n0n0n0n0a0"),
    ];
    let kept = [rec("one_a", "10M", "ACGTACGTAC", "0N0N0N0N0N0N0N0N0N1")];
    let all: Vec<&str> = dropped.iter().chain(&kept).map(String::as_str).collect();
    let input = sam(&dir, "probes.sam", &all);
    let names = |args: &[&str]| -> Vec<String> {
        let o = run(&[&["score", "--quiet"][..], args, &[s(&input)]].concat());
        assert!(o.status.success(), "{}", stderr(&o));
        String::from_utf8(o.stdout)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    assert_eq!(names(&["--pmdtools-compat"]), ["one_a"]);
    let j = damage_json(&input, &["--pmdtools-compat"]);
    assert_eq!(j["reads"]["skipped"]["no_reference_bases"], 6);
    assert_eq!(j["reads"]["scored"], 1);
    assert_eq!(j["format_version"], 4);
    // Default mode scores all seven (`=` bases need -r there, so that read is left out).
    let default_input = sam(
        &dir,
        "default.sam",
        &all.iter()
            .copied()
            .filter(|r| !r.starts_with("all_eq"))
            .collect::<Vec<_>>(),
    );
    let o = run(&["score", "--quiet", s(&default_input)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(String::from_utf8(o.stdout).unwrap().lines().count(), 7);
    let j = damage_json(&default_input, &[]);
    assert_eq!(j["reads"]["skipped"]["no_reference_bases"], 0);
}

/// Audit finding B6: PMDtools drops every read with a hard clip (`H`) in its CIGAR, whatever its
/// flags, strand or other operations, printing "only supports cigar operations M, I, S and D"
/// (docs/blackbox.md 1.6f; PMDtools 0.60 probes, 2026-11-22). Default mode scores them (D3,
/// docs/design.md §5), and so did `--pmdtools-compat`. Compat now drops them as
/// `unsupported_cigar`, as PMDtools does for `N` and `P`; with `--no-clips` they count as `clips`.
#[test]
fn compat_drops_hard_clipped_reads_like_pmdtools() {
    let dir = tmp_dir("hard_clips");
    let rec = |name: &str, flag: u16, cigar: &str, seq: &str, md: &str| {
        let q = "I".repeat(seq.len());
        format!("{name}\t{flag}\tc\t11\t37\t{cigar}\t*\t0\t0\t{seq}\t{q}\tMD:Z:{md}")
    };
    let t = "TAAAAAAAAA";
    let dropped = [
        rec("lead", 0, "3H10M", t, "0C9"),
        rec("trail", 0, "10M3H", t, "0C9"),
        rec("both", 16, "1H10M43H", t, "0C9"),
        rec("supp_r1", 2147, "43H10M", t, "0C9"),
        rec("supp_r2", 2195, "10M45H", t, "0C9"),
        rec("with_s", 0, "3H2S8M", "GGTAAAAAAA", "0C7"),
        rec("with_d", 0, "3H5M1D5M", t, "0C4^A5"),
        // Flagged unmapped but placed: compat scores these from MD (B4), but not with `H`.
        rec("placed", 4, "3H10M", t, "0C9"),
    ];
    let kept = [
        rec("plain", 0, "10M", t, "0C9"),
        rec("supp_soft", 2048, "3S10M", "GGGTAAAAAAAAA", "0C9"),
    ];
    let all: Vec<&str> = dropped.iter().chain(&kept).map(String::as_str).collect();
    let input = sam(&dir, "probes.sam", &all);
    let names = |args: &[&str]| -> Vec<String> {
        let o = run(&[&["score", "--quiet"][..], args, &[s(&input)]].concat());
        assert!(o.status.success(), "{}", stderr(&o));
        String::from_utf8(o.stdout)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    assert_eq!(names(&["--pmdtools-compat"]), ["plain", "supp_soft"]);
    let j = damage_json(&input, &["--pmdtools-compat"]);
    assert_eq!(j["reads"]["skipped"]["unsupported_cigar"], 8);
    assert_eq!(j["reads"]["scored"], 2);
    let j = damage_json(&input, &["--pmdtools-compat", "--no-clips"]);
    assert_eq!(j["reads"]["skipped"]["clips"], 9);
    assert_eq!(j["reads"]["skipped"]["unsupported_cigar"], 0);
    // Default mode scores every mapped one (the flag-4 read is skipped as unmapped).
    assert_eq!(names(&[]).len(), 9);
    let j = damage_json(&input, &[]);
    assert_eq!(j["reads"]["skipped"]["unsupported_cigar"], 0);
    assert_eq!(j["reads"]["skipped"]["unmapped"], 1);
}

// ---- Edge-case matrix (docs/audit-0.3.md §4, map in §4.5) ----

/// `score` output as (name, score, sites), one per scored read.
fn scores(args: &[&str]) -> Vec<(String, f64, u64)> {
    let o = run(&[&["score", "--quiet"][..], args].concat());
    assert!(o.status.success(), "{args:?}: {}", stderr(&o));
    String::from_utf8(o.stdout)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[0].to_string(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
            )
        })
        .collect()
}

/// Non-header lines of a SAM file.
fn records_of(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_string)
        .collect()
}

/// §4.1 layout: paired mates (both mapped, one unmapped, mate on another reference), secondary,
/// supplementary, duplicate and QC-fail reads are all scored, as PMDtools scores them, and counted
/// by flag. Filter writes them unchanged, including RNEXT and PNEXT.
#[test]
fn every_read_layout_is_scored_counted_and_kept_unchanged() {
    let dir = tmp_dir("layouts");
    let r = |name: &str, flag: u16, rnext: &str, pnext: u32| {
        format!(
            "{name}\t{flag}\tc\t11\t37\t10M\t{rnext}\t{pnext}\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9"
        )
    };
    let records = [
        r("single", 0, "*", 0),
        r("mate1", 99, "=", 11),
        r("mate2", 147, "=", 11),
        r("mate_unmapped", 73, "=", 11),
        r("mate_elsewhere", 65, "d", 500),
        r("secondary", 256, "*", 0),
        r("supplementary", 2048, "*", 0),
        r("duplicate", 1024, "*", 0),
        r("qc_fail", 512, "*", 0),
    ];
    let mut text = String::from("@HD\tVN:1.6\n@SQ\tSN:c\tLN:1000\n@SQ\tSN:d\tLN:1000\n");
    for rec in &records {
        text.push_str(rec);
        text.push('\n');
    }
    let input = dir.join("layouts.sam");
    fs::write(&input, &text).unwrap();
    for mode in [&[][..], &["--pmdtools-compat"]] {
        let got = scores(&[mode, &[s(&input)]].concat());
        assert_eq!(got.len(), records.len(), "{mode:?}");
        assert!(got.iter().all(|(_, _, n)| *n == 1), "{mode:?}: {got:?}");
    }
    let flags = damage_json(&input, &[])["reads"]["flags"].clone();
    for (flag, n) in [
        ("paired", 4),
        ("reverse", 1),
        ("secondary", 1),
        ("supplementary", 1),
        ("duplicate", 1),
        ("qc_fail", 1),
    ] {
        assert_eq!(flags[flag], n, "{flag}: {flags}");
    }
    let out = dir.join("kept.sam");
    let o = run(&["filter", "--quiet", "-t", "0", s(&input), "-o", s(&out)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(records_of(&out), records);
}

/// §4.1 length: a 1-bp read is scored by default (PMDtools drops it: docs/design.md §5), and a
/// site 1,200 bp from the 5′ end gets a finite, small score (PMDtools crashes on any site at
/// z >= 1000).
#[test]
fn one_base_and_long_reads_are_scored_by_default() {
    let dir = tmp_dir("lengths");
    let input = dir.join("lengths.sam");
    fs::write(
        &input,
        format!(
            "@SQ\tSN:c\tLN:5000\none\t0\tc\t11\t37\t1M\t*\t0\t0\tT\tI\tMD:Z:0C0\n\
             long\t0\tc\t1\t37\t1200M\t*\t0\t0\t{}T\t{}\tMD:Z:1199C0\n",
            "A".repeat(1199),
            "I".repeat(1200)
        ),
    )
    .unwrap();
    let got = scores(&[s(&input)]);
    assert_eq!(got[0], ("one".to_string(), 5.029717285061267, 1));
    let (name, pmds, n) = &got[1];
    assert_eq!((name.as_str(), *n), ("long", 1));
    assert!(pmds.is_finite() && *pmds > 0.0 && *pmds < 5.0, "{pmds}");
}

/// §4.1 SEQ: lowercase SEQ bases (SAM only) are read as uppercase in both modes. PMDtools does not
/// recognise them (`taaa` over MD `0C3` scores 0.0, 2026-11-22); compat does not reproduce that
/// (a documented limitation, docs/design.md §5). IUPAC and N read bases are never damage: PMDtools
/// scores these reads 5.0297 from their one C site, as pmdshodh does.
#[test]
fn lowercase_iupac_and_n_read_bases() {
    let dir = tmp_dir("seq_bases");
    let r = |name: &str, seq: &str, md: &str| {
        format!("{name}\t0\tc\t11\t37\t4M\t*\t0\t0\t{seq}\tIIII\tMD:Z:{md}")
    };
    let input = sam(
        &dir,
        "bases.sam",
        &[
            &r("upper", "TAAA", "0C3"),
            &r("lower", "taaa", "0C3"),
            &r("iupac", "TRYA", "0C0A0A1"),
            &r("n", "TNNA", "0C0A0A1"),
        ],
    );
    for mode in [&[][..], &["--pmdtools-compat"]] {
        let got = scores(&[mode, &[s(&input)]].concat());
        assert_eq!(got.len(), 4, "{mode:?}");
        for (name, pmds, n) in &got {
            assert_eq!((*pmds, *n), (5.029717285061267, 1), "{mode:?} {name}");
        }
    }
}

/// §4.1 CIGAR and MD: a mapped read with CIGAR `*`, and one with an empty MD tag, are input
/// errors in both modes. PMDtools scores the first from its MD (5.0297) and drops the second
/// silently (2026-11-22); both are listed in docs/design.md §5.
#[test]
fn star_cigar_and_empty_md_on_mapped_reads_are_input_errors() {
    let dir = tmp_dir("star_cigar");
    let star = sam(
        &dir,
        "star.sam",
        &["r1\t0\tc\t11\t37\t*\t*\t0\t0\tTAAA\tIIII\tMD:Z:0C3"],
    );
    let empty = sam(
        &dir,
        "empty_md.sam",
        &["r1\t0\tc\t11\t37\t4M\t*\t0\t0\tTAAA\tIIII\tMD:Z:"],
    );
    for mode in [&[][..], &["--pmdtools-compat"]] {
        let star_args = [&["score"][..], mode, &[s(&star)]].concat();
        input_error(&star_args, "CIGAR covers 0 read bases");
        let empty_args = [&["score"][..], mode, &[s(&empty)]].concat();
        input_error(&empty_args, "MD is shorter than the CIGAR");
    }
}

/// §4.1 position: reads at the first and last reference base, and on a reference of length 1,
/// score the same from MD and from `-r`, in default and compat mode.
#[test]
fn reads_at_reference_ends_score_the_same_from_md_and_reference() {
    let dir = tmp_dir("reference_ends");
    let fasta = dir.join("ref.fa");
    fs::write(&fasta, ">c\nCAAAAAAAAC\n>one\nC\n").unwrap();
    fs::write(
        dir.join("ref.fa.fai"),
        "c\t10\t3\t10\t11\none\t1\t19\t1\t2\n",
    )
    .unwrap();
    let input = dir.join("ends.sam");
    fs::write(
        &input,
        "@SQ\tSN:c\tLN:10\n@SQ\tSN:one\tLN:1\n\
         start\t0\tc\t1\t37\t4M\t*\t0\t0\tTAAA\tIIII\tMD:Z:0C3\n\
         end\t16\tc\t7\t37\t4M\t*\t0\t0\tAAAT\tIIII\tMD:Z:3C0\n\
         length_1\t0\tone\t1\t37\t1M\t*\t0\t0\tT\tI\tMD:Z:0C0\n",
    )
    .unwrap();
    let md = scores(&[s(&input)]);
    assert_eq!(md.len(), 3);
    assert!(md.iter().all(|(_, _, n)| *n == 1), "{md:?}");
    assert_eq!(scores(&["-r", s(&fasta), s(&input)]), md);
    // The 1-bp read is left out of compat (PMDtools drops 1-bp reads).
    let compat = |extra: &[&str]| {
        let mut v = scores(&[&["--pmdtools-compat"][..], extra, &[s(&input)]].concat());
        v.retain(|(name, _, _)| name != "length_1");
        v
    };
    assert_eq!(compat(&["-r", s(&fasta)]), compat(&[]));
    assert_eq!(compat(&[]), md[..2]);
}

/// §4.2 model: CpG context is taken along SEQ, skipping deleted reference bases, in both modes,
/// as in PMDtools (docs/blackbox.md 1.6b): a C, a deleted A and then a G count as CpG. PMDtools
/// 0.60 gives 5.0030 (2 sites) without treatment, with --CpG and with --UDGhalf (2026-11-22).
#[test]
fn cpg_context_skips_deleted_reference_bases() {
    let dir = tmp_dir("cpg_deletion");
    let input = sam(
        &dir,
        "cpg.sam",
        &["r1\t0\tc\t11\t37\t1M1D9M\t*\t0\t0\tTGAAAAAAAA\tIIIIIIIIII\tMD:Z:0C0^A9"],
    );
    for treatment in [&[][..], &["--CpG"], &["--UDGhalf"]] {
        for mode in [&[][..], &["--pmdtools-compat"]] {
            let got = scores(&[treatment, mode, &[s(&input)]].concat());
            assert_eq!(got[0].2, 2, "{treatment:?} {mode:?}");
            assert!(
                (got[0].1 - 5.00304487236408).abs() < 1e-12,
                "{treatment:?} {mode:?}: {}",
                got[0].1
            );
        }
    }
}

/// §4.2 parameters: `-q` counts a base at exactly its value and ignores one below; `-q` and `-m`
/// take 0 to 255; the priors take 0 to 1; `--pmd-p 0` and `--pmd-c 0` are accepted.
#[test]
fn quality_and_parameter_boundaries() {
    let dir = tmp_dir("boundaries");
    // Base quality 20 ('5') at the damaged C; MAPQ 37 and 255.
    let input = sam(
        &dir,
        "q.sam",
        &[
            "q20\t0\tc\t11\t37\t4M\t*\t0\t0\tTAAA\t5III\tMD:Z:0C3",
            "mapq255\t0\tc\t11\t255\t4M\t*\t0\t0\tTAAA\tIIII\tMD:Z:0C3",
        ],
    );
    let sites = |args: &[&str]| -> Vec<(String, u64)> {
        scores(&[args, &[s(&input)]].concat())
            .into_iter()
            .map(|(name, _, n)| (name, n))
            .collect()
    };
    let pair = |a: u64, b: u64| vec![("q20".to_string(), a), ("mapq255".to_string(), b)];
    assert_eq!(sites(&["-q", "20"]), pair(1, 1));
    assert_eq!(sites(&["-q", "21"]), pair(0, 1));
    assert_eq!(sites(&["-q", "0"]), pair(1, 1));
    assert_eq!(sites(&["-q", "255"]), pair(0, 0));
    assert_eq!(sites(&["-m", "0"]), pair(1, 1));
    assert_eq!(sites(&["-m", "255"]), vec![("mapq255".to_string(), 1)]);
    for value in ["0", "1"] {
        for prior in ["--pi-ancient", "--pi-contamination"] {
            assert_eq!(sites(&[prior, value]).len(), 2, "{prior} {value}");
        }
    }
    assert_eq!(sites(&["--pmd-p", "0"]).len(), 2);
    assert_eq!(sites(&["--pmd-c", "0"]).len(), 2);
    for (args, says) in [
        (&["--pi-ancient", "1.1"][..], "between 0 and 1"),
        (&["--pi-contamination=-0.1"], "between 0 and 1"),
        (&["-q", "256"], "256"),
        (&["-m", "256"], "256"),
    ] {
        let o = run(&[&["score"][..], args, &[s(&input)]].concat());
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(stderr(&o).contains(says), "{args:?}: {}", stderr(&o));
    }
}

/// §4.2 thresholds: a read scoring exactly `--threshold` is kept, and one scoring exactly
/// `--upper-threshold` is not, as in PMDtools (docs/blackbox.md 1.7).
#[test]
fn reads_exactly_at_a_threshold() {
    let dir = tmp_dir("threshold_exact");
    let input = sam(
        &dir,
        "t.sam",
        &["r1\t0\tc\t11\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9"],
    );
    let score = scores(&[s(&input)])[0].1.to_string();
    let out = dir.join("kept.sam");
    let kept = |bound: &str| {
        let o = run(&["filter", "--quiet", bound, &score, s(&input), "-o", s(&out)]);
        assert!(o.status.success(), "{}", stderr(&o));
        records_of(&out).len()
    };
    assert_eq!(kept("--threshold"), 1);
    assert_eq!(kept("--upper-threshold"), 0);
}

/// §4.2 output: an existing DS:Z tag is left untouched without `--tag` (PMDtools crashes on it
/// unless --writesamfield), and `--tag` replaces an existing tag of its name, leaving one.
#[test]
fn existing_score_tags_are_kept_or_replaced() {
    let dir = tmp_dir("existing_tags");
    let input = sam(
        &dir,
        "tags.sam",
        &["r1\t0\tc\t11\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9\tDS:Z:9.999\tpm:f:1.5"],
    );
    let out = dir.join("kept.sam");
    let kept = |extra: &[&str]| -> String {
        let fixed = ["filter", "--quiet", "-t", "0"];
        let o = run(&[&fixed[..], extra, &[s(&input), "-o", s(&out)]].concat());
        assert!(o.status.success(), "{}", stderr(&o));
        records_of(&out).remove(0)
    };
    let plain = kept(&[]);
    assert!(plain.ends_with("\tDS:Z:9.999\tpm:f:1.5"), "{plain}");
    let ds = kept(&["--tag", "--tag-name", "DS"]);
    assert_eq!(ds.matches("DS:Z:").count(), 1, "{ds}");
    assert!(
        ds.contains("\tDS:Z:5.03") && ds.contains("\tpm:f:1.5"),
        "{ds}"
    );
    let pm = kept(&["--tag"]);
    assert_eq!(pm.matches("pm:f:").count(), 1, "{pm}");
    assert!(
        pm.contains("\tDS:Z:9.999") && !pm.contains("pm:f:1.5"),
        "{pm}"
    );
}

/// §4.3 damage: `-n` 0, 1 and more than the input; `--range 1`; header-only input; and an input
/// whose reads are all skipped, which still gives a summary that `report` renders.
#[test]
fn damage_limits_and_empty_runs() {
    let dir = tmp_dir("damage_limits");
    let input = sam(
        &dir,
        "two.sam",
        &[
            "r1\t0\tc\t11\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
            "r2\t0\tc\t21\t37\t10M\t*\t0\t0\tTAAAAAAAAA\tIIIIIIIIII\tMD:Z:0C9",
        ],
    );
    for (n, want) in [("0", 0), ("1", 1), ("100", 2)] {
        let j = damage_json(&input, &["-n", n]);
        assert_eq!(j["reads"]["input"], want, "-n {n}");
        assert_eq!(j["reads"]["scored"], want, "-n {n}");
    }
    let j = damage_json(&input, &["--range", "1"]);
    assert_eq!(j["damage"]["range"], 1);
    assert_eq!(
        j["damage"]["five_prime"]["C"]["T"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let header_only = sam(&dir, "header.sam", &[]);
    let skipped = sam(
        &dir,
        "skipped.sam",
        &[
            "u1\t4\t*\t0\t0\t*\t*\t0\t0\tTAAA\tIIII",
            "u2\t4\t*\t0\t0\t*\t*\t0\t0\tTAAA\tIIII",
        ],
    );
    for (input, reads) in [(&header_only, 0), (&skipped, 2)] {
        let summary = dir.join("s.json");
        let o = run(&["damage", "--quiet", s(input), "-o", s(&summary)]);
        assert!(o.status.success(), "{}", stderr(&o));
        let j: serde_json::Value = serde_json::from_slice(&fs::read(&summary).unwrap()).unwrap();
        assert_eq!(j["reads"]["input"], reads, "{}", s(input));
        assert_eq!(j["reads"]["scored"], 0, "{}", s(input));
        let page = dir.join("r.html");
        let o = run(&["report", "--quiet", s(&summary), "-o", s(&page)]);
        assert!(o.status.success(), "{}", stderr(&o));
        assert!(!fs::read_to_string(&page).unwrap().contains("NaN"));
    }
}

/// §4.3 report: 50 samples, a Unicode and a 200-character file name, `--range 1000`, a UDG-half
/// compat summary and a CRAM summary in one report; it is byte-identical when made twice.
#[test]
fn report_on_many_and_unusual_summaries() {
    let dir = tmp_dir("report_unusual");
    let fixture = data("reads_random.sam");
    let unicode = dir.join("échantillon_日本.sam");
    let long = dir.join(format!("{}.sam", "x".repeat(200)));
    fs::copy(&fixture, &unicode).unwrap();
    fs::copy(&fixture, &long).unwrap();
    let summary = |name: &str, args: &[&str]| -> PathBuf {
        let out = dir.join(name);
        let o = run(&[&["damage", "--quiet"][..], args, &["-o", s(&out)]].concat());
        assert!(o.status.success(), "{}", stderr(&o));
        out
    };
    let (fa, cram) = (data("pairs.fa"), data("pairs.cram"));
    let mut inputs = vec![
        summary("unicode.json", &[s(&unicode)]),
        summary("long.json", &[s(&long)]),
        summary("range.json", &["--range", "1000", s(&fixture)]),
        summary("cram.json", &["-r", s(&fa), s(&cram)]),
        summary(
            "half.json",
            &["--pmdtools-compat", "--UDGhalf", s(&fixture)],
        ),
    ];
    for i in inputs.len()..50 {
        let p = dir.join(format!("copy{i}.json"));
        fs::copy(&inputs[4], &p).unwrap();
        inputs.push(p);
    }
    let page = |name: &str| {
        let out = dir.join(name);
        let mut args = vec!["report", "--quiet"];
        args.extend(inputs.iter().map(|p| s(p)));
        args.extend(["-o", s(&out)]);
        let o = run(&args);
        assert!(o.status.success(), "{}", stderr(&o));
        fs::read(out).unwrap()
    };
    let (a, b) = (page("a.html"), page("b.html"));
    assert_eq!(a, b, "byte-identical");
    let html = String::from_utf8(a).unwrap();
    assert!(html.contains("échantillon_日本.sam"));
    assert!(html.contains(&"x".repeat(200)));
    assert_eq!(
        html.matches("<section").count(),
        51,
        "50 samples and the overlays"
    );
    assert!(!html.contains("NaN"));
}

/// §4.4 command line: an unknown option, an option given twice, and an option missing its value
/// are usage errors (exit 2) in every subcommand, with nothing on stdout.
#[test]
fn unknown_repeated_and_valueless_options_are_usage_errors() {
    let input = data("reads_random.sam");
    let input = s(&input);
    for cmd in [
        &["score"][..],
        &["damage"],
        &["filter", "-t", "3", "-o", "x.bam"],
    ] {
        for args in [&["--bogus"][..], &["-q", "1", "-q", "2"], &["--range"]] {
            let o = run(&[cmd, args, &[input]].concat());
            assert_eq!(o.status.code(), Some(2), "{cmd:?} {args:?}");
            assert!(o.stdout.is_empty(), "{cmd:?} {args:?}");
        }
    }
    let o = run(&["report", "--bogus", "x.json"]);
    assert_eq!(o.status.code(), Some(2));
}

/// §4.1 CIGAR: an insertion or deletion right after a soft clip. Compat equals PMDtools 0.60
/// (`--printDS`, 2026-11-22: 4.0254, 0.0 and 0.0). Default mode places the deletion correctly
/// (docs/blackbox.md 1.6c), so the C→T scores as the same read without the deletion does.
#[test]
fn indels_next_to_a_soft_clip() {
    let dir = tmp_dir("clip_indel");
    let r = |name: &str, flag: u16, cigar: &str, seq: &str, md: &str| {
        let q = "I".repeat(seq.len());
        format!("{name}\t{flag}\tc\t11\t37\t{cigar}\t*\t0\t0\t{seq}\t{q}\tMD:Z:{md}")
    };
    let input = sam(
        &dir,
        "clip_indel.sam",
        &[
            &r("s_i", 0, "2S1I8M", "GGATAAAAAAA", "0C7"),
            &r("s_d", 0, "2S1D8M", "GGTAAAAAAA", "^A0C7"),
            &r("s_d_rev", 16, "2S1D8M", "GGAAAAAAAT", "^A7C0"),
            &r("no_d", 0, "2S8M", "GGTAAAAAAA", "0C7"),
        ],
    );
    let values = |mode: &[&str]| -> Vec<f64> {
        scores(&[mode, &[s(&input)]].concat())
            .into_iter()
            .map(|(_, p, _)| p)
            .collect()
    };
    let compat = values(&["--pmdtools-compat"]);
    assert_eq!(compat[..3], [4.0254364447734625, 0.0, 0.0]);
    let default = values(&[]);
    assert_eq!(default[0], compat[0]);
    assert_eq!(default[1], default[3], "deletion leaves z unchanged");
    assert!(default[2] > 0.0);
}

/// §4.4 outputs: `-o /dev/stdout` writes in place, as `-o -` does (Unix only).
#[cfg(unix)]
#[test]
fn output_to_dev_stdout() {
    let input = data("reads_random.sam");
    // The @PG line records the command line, so it differs between the two.
    let sam_without_pg = |out: &str| -> String {
        let fixed = ["filter", "--quiet", "-t", "3", "--output-format", "sam"];
        let o = run(&[&fixed[..], &[s(&input), "-o", out]].concat());
        assert!(o.status.success(), "{out}: {}", stderr(&o));
        String::from_utf8(o.stdout)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with("@PG"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let via_dev = sam_without_pg("/dev/stdout");
    assert!(via_dev.lines().any(|l| !l.starts_with('@')), "has records");
    assert_eq!(via_dev, sam_without_pg("-"));
}

/// Audit finding B8: PMDtools 0.60 drops every read whose SEQ is 1 base, silently, whatever its
/// CIGAR op, strand or flags; a 2-base SEQ with one aligned base is scored (docs/blackbox.md
/// 1.6g; probes 2026-11-22). `--pmdtools-compat` scored them. It now drops them, counted as
/// `too_short`; default mode still scores them (docs/design.md §5).
#[test]
fn compat_drops_one_base_reads_like_pmdtools() {
    let dir = tmp_dir("one_base");
    let rec = |name: &str, flag: u16, cigar: &str, seq: &str, md: &str| {
        let q = "I".repeat(seq.len());
        format!("{name}\t{flag}\tc\t11\t37\t{cigar}\t*\t0\t0\t{seq}\t{q}\tMD:Z:{md}")
    };
    let dropped = [
        rec("damaged", 0, "1M", "T", "0C0"),
        rec("reverse", 16, "1M", "T", "0C0"),
        rec("match", 0, "1M", "A", "1"),
        rec("x_op", 0, "1X", "T", "0C0"),
        // Flagged unmapped but placed: compat scores these from MD (B4), but not at 1 base.
        rec("placed", 4, "1M", "T", "0C0"),
    ];
    // PMDtools 0.60 --printDS: 4.6881, 5.0297 and 5.0297.
    let kept = [
        rec("clip", 0, "1S1M", "GT", "0C0"),
        rec("insertion", 0, "1M1I", "TG", "0C0"),
        rec("two", 0, "2M", "TA", "0C1"),
    ];
    let all: Vec<&str> = dropped.iter().chain(&kept).map(String::as_str).collect();
    let input = sam(&dir, "one_base.sam", &all);
    let compat = scores(&["--pmdtools-compat", s(&input)]);
    let names: Vec<&str> = compat.iter().map(|(n, _, _)| n.as_str()).collect();
    assert_eq!(names, ["clip", "insertion", "two"]);
    let values: Vec<f64> = compat.iter().map(|(_, p, _)| *p).collect();
    assert!((values[0] - 4.688132074164015).abs() < 1e-12, "{values:?}");
    assert_eq!(values[1..], [5.029717285061267, 5.029717285061267]);
    let j = damage_json(&input, &["--pmdtools-compat"]);
    assert_eq!(j["reads"]["skipped"]["too_short"], 5);
    assert_eq!(j["reads"]["scored"], 3);
    assert_eq!(j["format_version"], 4);
    // They are not written by filter either, at any threshold.
    let out = dir.join("kept.sam");
    let o = run(&[
        "filter",
        "--quiet",
        "--pmdtools-compat",
        "-t",
        "-5",
        s(&input),
        "-o",
        s(&out),
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(records_of(&out).len(), 3);
    // Default mode scores the seven mapped reads; the flag-4 one is skipped as unmapped.
    assert_eq!(scores(&[s(&input)]).len(), 7);
    let j = damage_json(&input, &[]);
    assert_eq!(j["reads"]["skipped"]["too_short"], 0);
    assert_eq!(j["reads"]["skipped"]["unmapped"], 1);
}
