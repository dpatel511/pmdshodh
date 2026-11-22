//! Edge cases from the 0.3 audit coverage map (docs/audit-0.3.md §3.1, group a): read filters and
//! skip reasons, full UDG from the command line, `=`/`X` CIGARs, malformed MD tags, reference
//! checks, damaged BGZF input and command-line errors.

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
    assert_eq!(j["format_version"], 3);
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
