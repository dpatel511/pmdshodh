//! Command-line behaviour: output content, streams, exit codes (CLAUDE.md §11–12, docs/design.md §4).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles_util::alignment;

const BIN: &str = env!("CARGO_BIN_EXE_pmdshodh");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// A fresh directory per test under Cargo's target temp dir.
fn tmp_dir(test: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(test);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().unwrap()
}

fn read_records(path: &Path) -> (noodles::sam::Header, Vec<RecordBuf>) {
    let mut reader = alignment::io::reader::Builder::default()
        .build_from_path(path)
        .unwrap();
    let header = reader.read_header().unwrap();
    let records = reader
        .records(&header)
        .map(|r| RecordBuf::try_from_alignment_record(&header, &r.unwrap()).unwrap())
        .collect();
    (header, records)
}

fn name(r: &RecordBuf) -> String {
    r.name().unwrap().to_string()
}

/// PMDtools --threshold 3 --writesamfield output, as (read name, DS text).
fn pmdtools_ds3() -> Vec<(String, String)> {
    fs::read_to_string(data("reads_random_pmdtools_ds3.tsv"))
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("read0"))
        .map(|l| {
            let (n, ds) = l.split_once('\t').unwrap();
            (n.to_string(), ds.to_string())
        })
        .collect()
}

#[test]
fn filter_keeps_the_reads_pmdtools_keeps_with_identical_ds_tags() {
    let dir = tmp_dir("filter_ds");
    let out = dir.join("kept.bam");
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "--pmdtools-compat",
        "-t",
        "3",
        "--tag",
        "--tag-name",
        "DS",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let (_, records) = read_records(&out);
    let got: Vec<(String, String)> = records
        .iter()
        .map(|r| {
            let ds = match r.data().get(&Tag::new(b'D', b'S')) {
                Some(Value::String(s)) => s.to_string(),
                other => panic!("{}: DS tag {other:?}", name(r)),
            };
            (name(r), ds)
        })
        .collect();
    let want = pmdtools_ds3();
    assert_eq!(got.len(), 483);
    assert_eq!(got, want);
}

#[test]
fn kept_records_are_unchanged_apart_from_the_score_tag() {
    let dir = tmp_dir("filter_unchanged");
    let out = dir.join("kept.bam");
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        "--tag",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let (_, inputs) = read_records(&input);
    let (_, outputs) = read_records(&out);
    let pm = Tag::new(b'p', b'm');
    let mut it = inputs.iter();
    for o in &outputs {
        let i = it
            .by_ref()
            .find(|i| name(i) == name(o))
            .expect("input order kept");
        assert!(
            matches!(o.data().get(&pm), Some(Value::Float(_))),
            "{}: pm:f tag",
            name(o)
        );
        let mut stripped = o.clone();
        stripped.data_mut().remove(&pm);
        assert_eq!(&stripped, i, "{} changed", name(o));
    }
    assert!(!outputs.is_empty());
}

#[test]
fn output_header_gains_one_pg_line() {
    let dir = tmp_dir("filter_pg");
    let out = dir.join("kept.sam");
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let text = fs::read_to_string(&out).unwrap();
    let pg: Vec<&str> = text.lines().filter(|l| l.starts_with("@PG")).collect();
    assert_eq!(pg.len(), 1);
    assert!(pg[0].contains("ID:pmdshodh") && pg[0].contains("CL:") && pg[0].contains("filter"));
}

#[test]
fn stdout_carries_only_the_bam_and_the_summary_goes_to_stderr() {
    let input = data("reads_random.sam");
    let o = run(&["filter", "-t", "3", input.to_str().unwrap(), "-o", "-"]);
    assert!(o.status.success());
    assert_eq!(&o.stdout[..2], &[0x1f, 0x8b], "BGZF/gzip magic on stdout");
    let stderr = String::from_utf8(o.stderr).unwrap();
    assert!(stderr.contains("1500 reads read") && stderr.contains("kept"));

    let quiet = run(&[
        "filter",
        "-t",
        "3",
        "--quiet",
        input.to_str().unwrap(),
        "-o",
        "-",
    ]);
    assert!(quiet.status.success());
    assert!(quiet.stderr.is_empty());
    // The @PG command line differs (it records --quiet), so compare the decoded records.
    let dir = tmp_dir("filter_stdout");
    let (a, b) = (dir.join("a.bam"), dir.join("b.bam"));
    fs::write(&a, &o.stdout).unwrap();
    fs::write(&b, &quiet.stdout).unwrap();
    let (_, ra) = read_records(&a);
    let (_, rb) = read_records(&b);
    assert!(!ra.is_empty());
    assert_eq!(ra, rb, "--quiet must not change the records");
}

#[test]
fn invalid_input_exits_2_and_leaves_no_output_file() {
    let dir = tmp_dir("filter_invalid");
    let sam = dir.join("no_md.sam");
    fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c\tLN:100\nok\t0\tc\t1\t37\t4M\t*\t0\t0\tTAAC\tIIII\tMD:Z:0C3\nbad\t0\tc\t1\t37\t4M\t*\t0\t0\tTAAC\tIIII\n",
    )
    .unwrap();
    let out = dir.join("out.bam");
    let o = run(&[
        "filter",
        "-t",
        "0",
        sam.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("read bad"));
    assert!(!out.exists(), "no partial output under the final name");
    assert_eq!(
        fs::read_dir(&dir).unwrap().count(),
        1,
        "temporary file removed"
    );
}

#[test]
fn usage_errors_exit_2() {
    let input = data("reads_random.sam");
    let input = input.to_str().unwrap();
    for args in [
        vec!["filter", input, "-o", "x.bam"], // no threshold
        vec![
            "filter",
            "-t",
            "3",
            "--ss",
            "--library",
            "ds",
            input,
            "-o",
            "x.bam",
        ],
        vec!["filter", "-t", "3", input, "-o", "x.cram"], // unsupported output
        vec!["score", "--pmd-p", "1.5", input],           // parameter out of range
        vec!["score", "missing-file.sam"],
    ] {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(o.stdout.is_empty() || args[0] == "score", "{args:?}");
    }
}

#[test]
fn unwritable_output_exits_1() {
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        input.to_str().unwrap(),
        "-o",
        "no/such/dir/out.bam",
    ]);
    assert_eq!(o.status.code(), Some(1));
}

#[test]
fn empty_input_gives_an_empty_valid_bam() {
    let dir = tmp_dir("filter_empty");
    let sam = dir.join("empty.sam");
    fs::write(&sam, "@HD\tVN:1.6\n@SQ\tSN:c\tLN:100\n").unwrap();
    let out = dir.join("out.bam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        sam.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let (header, records) = read_records(&out);
    assert!(records.is_empty());
    assert_eq!(header.reference_sequences().len(), 1);
}

#[test]
fn help_and_version_go_to_stdout_with_status_0() {
    for flag in ["--help", "--version"] {
        let o = run(&[flag]);
        assert!(o.status.success());
        assert!(!o.stdout.is_empty() && o.stderr.is_empty());
    }
}

#[test]
fn ds_tag_text_matches_python_round_formatting() {
    // Cases from PMDtools output: trailing zeros dropped, at least one decimal kept.
    let dir = tmp_dir("filter_ds_text");
    let out = dir.join("kept.sam");
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "--pmdtools-compat",
        "-t",
        "3",
        "--tag",
        "--tag-name",
        "DS",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let text = fs::read_to_string(&out).unwrap();
    let ours: Vec<String> = text
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.rsplit_once("DS:Z:").unwrap().1.to_string())
        .collect();
    let want: Vec<String> = pmdtools_ds3().into_iter().map(|(_, ds)| ds).collect();
    assert_eq!(ours, want);
}

fn damage_json(args: &[&str]) -> serde_json::Value {
    let input = data("reads_random.sam");
    let mut all = vec!["damage", "--quiet"];
    all.extend_from_slice(args);
    all.push(input.to_str().unwrap());
    let o = run(&all);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    serde_json::from_slice(&o.stdout).unwrap()
}

/// With --pmdtools-compat, the damage table equals PMDtools --deamination --range 30
/// (5′ columns at reference C, 3′ columns at reference G; PMDtools prints 5 decimals).
/// UDG-half tables are not compared: PMDtools' rule for them is not yet known (docs/blackbox.md 1.6d).
#[test]
fn damage_table_matches_pmdtools_deamination() {
    let table = fs::read_to_string(data("reads_random_pmdtools_deamination.tsv")).unwrap();
    // UDG-half takes a C's CpG context from the reference, --CpG from the next read base
    // (docs/blackbox.md 1.6d); UDG-half was not compared until that was found (2026-11-20).
    for (treatment, flag) in [
        ("none", "--UDGminus"),
        ("half", "--UDGhalf"),
        ("full", "--CpG"),
    ] {
        let j = damage_json(&["--pmdtools-compat", "--range", "30", flag]);
        let d = &j["damage"];
        let frac = |end: &str, refb: &str, readb: &str, z: usize| {
            let col = |b: &str| d[end][refb][b][z].as_u64().unwrap() as f64;
            col(readb) / ["A", "C", "G", "T"].iter().map(|b| col(b)).sum::<f64>()
        };
        let rows: Vec<Vec<&str>> = table
            .lines()
            .filter(|l| l.starts_with(&format!("{treatment}\t")))
            .map(|l| l.split('\t').collect())
            .collect();
        assert_eq!(rows.len(), 30);
        for row in rows {
            let z: usize = row[1].parse().unwrap();
            let ours = [
                frac("five_prime", "C", "T", z),
                frac("five_prime", "C", "A", z),
                frac("five_prime", "C", "G", z),
                frac("five_prime", "C", "C", z),
                frac("three_prime", "G", "A", z),
                frac("three_prime", "G", "T", z),
                frac("three_prime", "G", "C", z),
                frac("three_prime", "G", "G", z),
            ];
            for (k, (o, want)) in ours.iter().zip(&row[2..]).enumerate() {
                let want: f64 = want.parse().unwrap();
                assert!(
                    (o - want).abs() <= 5e-6,
                    "{treatment} z={z} column {k}: pmdshodh {o:.6}, PMDtools {want}"
                );
            }
        }
    }
}

#[test]
fn summary_json_counts_and_versions() {
    let j = damage_json(&[]);
    assert_eq!(j["format_version"], 2);
    assert_eq!(j["reference"]["bases_from"], "md");
    assert!(j["reference"]["md_tags"].is_null());
    assert_eq!(j["reads"]["input"], 1500);
    assert_eq!(j["reads"]["scored"], 1500);
    assert!(j["reads"]["kept"].is_null());
    let hist: u64 = j["score_histogram"]["counts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap())
        .sum::<u64>()
        + j["score_histogram"]["underflow"].as_u64().unwrap()
        + j["score_histogram"]["overflow"].as_u64().unwrap();
    assert_eq!(hist, 1500, "every scored read is in the histogram");
    assert_eq!(j["mode"]["library"], "ds");
    assert_eq!(j["model"]["null_background"], 0.001);
}

#[test]
fn filter_summary_counts_kept_and_removed() {
    let dir = tmp_dir("filter_summary");
    let (out, js) = (dir.join("kept.bam"), dir.join("s.json"));
    let input = data("reads_random.sam");
    let o = run(&[
        "filter",
        "--pmdtools-compat",
        "-t",
        "3",
        "--summary",
        js.to_str().unwrap(),
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let j: serde_json::Value = serde_json::from_str(&fs::read_to_string(&js).unwrap()).unwrap();
    assert_eq!(j["reads"]["kept"], 483);
    assert_eq!(j["reads"]["removed_by_threshold"], 1500 - 483);
    assert_eq!(j["filters"]["threshold"], 3.0);
}

#[test]
fn damage_max_reads_stops_early_with_status_0() {
    let j = damage_json(&["-n", "100"]);
    assert_eq!(j["reads"]["input"], 100);
}

/// Multi-threaded output must be byte-identical to single-threaded output (CLAUDE.md §12).
/// The input repeats the fixture's records to span several batches and many BGZF blocks.
#[test]
fn output_is_identical_for_any_thread_count() {
    let dir = tmp_dir("threads");
    let text = fs::read_to_string(data("reads_random.sam")).unwrap();
    let (header, body): (Vec<&str>, Vec<&str>) = text.lines().partition(|l| l.starts_with('@'));
    let mut big = header.join("\n") + "\n";
    for _ in 0..17 {
        big += &(body.join("\n") + "\n");
    }
    let sam = dir.join("big.sam");
    fs::write(&sam, big).unwrap();

    let mut outputs = Vec::new();
    for threads in ["1", "4"] {
        let bam = dir.join(format!("out{threads}.bam"));
        let json = dir.join(format!("out{threads}.json"));
        let o = run(&[
            "filter",
            "-t",
            "3",
            "--tag",
            "-@",
            threads,
            "--quiet",
            "--summary",
            json.to_str().unwrap(),
            sam.to_str().unwrap(),
            "-o",
            bam.to_str().unwrap(),
        ]);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        // Score the BAM we just wrote, exercising multi-threaded BGZF reading.
        let tsv = run(&["score", "-@", threads, "--quiet", bam.to_str().unwrap()]);
        assert!(tsv.status.success());
        let mut summary: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&json).unwrap()).unwrap();
        summary["command"] = serde_json::Value::Null; // the command line records -@
        outputs.push((fs::read(&bam).unwrap(), summary, tsv.stdout));
    }
    let (a, b) = (&outputs[0], &outputs[1]);
    assert_eq!(a.1["reads"]["input"], 25_500);
    assert_eq!(a.1, b.1, "summary differs between 1 and 4 threads");
    assert!(a.2 == b.2, "score table differs between 1 and 4 threads");
    // The BAM bytes differ only because @PG CL records the command line (-@ 1 vs -@ 4).
    let decode = |t: &str| {
        let (h, records) = read_records(&dir.join(format!("out{t}.bam")));
        let mut w = noodles::sam::io::Writer::new(Vec::new());
        w.write_header(&h).unwrap();
        let text = String::from_utf8(w.get_ref().clone()).unwrap();
        let header: Vec<String> = text
            .lines()
            .map(|l| {
                l.split('\t')
                    .filter(|f| !f.starts_with("CL:"))
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect();
        (header, records)
    };
    let (ha, ra) = decode("1");
    let (hb, rb) = decode("4");
    assert_eq!(ha, hb, "header differs beyond the @PG command line");
    assert_eq!(ra.len(), rb.len());
    assert!(ra == rb, "records differ between 1 and 4 threads");
}

/// A one-read SAM file (all bases Q40, all-M CIGAR).
fn one_read_sam(path: &Path, flag: u16, seq: &str, md: &str) -> PathBuf {
    let path = path.to_path_buf();
    let n = seq.len();
    let qual = "I".repeat(n);
    let sam = format!(
        "@SQ\tSN:r\tLN:1000\nx\t{flag}\tr\t100\t60\t{n}M\t*\t0\t0\t{seq}\t{qual}\tMD:Z:{md}\n"
    );
    fs::write(&path, sam).unwrap();
    path
}

fn single_score(input: &Path, extra: &[&str]) -> f64 {
    let mut args = vec!["score", "--quiet"];
    args.extend_from_slice(extra);
    args.push(input.to_str().unwrap());
    let o = run(&args);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let out = String::from_utf8(o.stdout).unwrap();
    out.lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .parse()
        .unwrap()
}

// Lowercase MD bases (outside the SAM spec; seen in real data, docs/benchmarks). PMDtools 0.60 never
// counts them as C or G (--printDS, 2026-11-13). Default mode reads them as uppercase (departure,
// docs/design.md §5); --pmdtools-compat reproduces PMDtools.
#[test]
fn lowercase_md_bases_default_reads_uppercase_compat_matches_pmdtools() {
    let dir = tmp_dir("lowercase_md");
    let a = "A".repeat(28);
    // (flag, seq, lowercase MD, uppercase MD, options, PMDtools score with uppercase MD)
    let cases = [
        (
            0,
            format!("T{a}A"),
            "0c29",
            "0C29",
            vec![],
            5.029717285061266,
        ),
        (
            16,
            format!("{a}AA"),
            "29g0",
            "29G0",
            vec![],
            5.029717285061266,
        ),
        (
            0,
            format!("CA{a}"),
            "1g28",
            "1G28",
            vec!["--CpG"],
            1.3220303836508658,
        ),
    ];
    for (flag, seq, lower, upper, opts, pmdtools_upper) in cases {
        let lc = one_read_sam(&dir.join(format!("lower_{lower}.sam")), flag, &seq, lower);
        let uc = one_read_sam(&dir.join(format!("upper_{upper}.sam")), flag, &seq, upper);
        let mut compat = vec!["--pmdtools-compat"];
        compat.extend(&opts);
        let close = |x: f64, y: f64| (x - y).abs() < 1e-12;
        assert!(
            close(single_score(&uc, &compat), pmdtools_upper),
            "{upper} compat"
        );
        assert!(
            close(single_score(&lc, &opts), pmdtools_upper),
            "{lower} default"
        );
        // PMDtools 0.60 --printDS gives 0.0 for each lowercase case.
        assert_eq!(single_score(&lc, &compat), 0.0, "{lower} compat");
    }
}

#[test]
fn score_dash_output_writes_to_stdout() {
    let dir = tmp_dir("score_dash");
    let input = data("reads_random.sam");
    let o = Command::new(BIN)
        .current_dir(&dir)
        .args(["score", "--quiet", input.to_str().unwrap(), "-o", "-"])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(
        String::from_utf8(o.stdout)
            .unwrap()
            .starts_with("read_name\tflag\tpmds\tn_sites\n")
    );
    assert!(!dir.join("-").exists(), "a file named '-' was created");
}

// PMDtools 0.60 drops reference characters other than ACGTNRY/acgt (e.g. lowercase n) when it
// reverse-complements a reverse-strand read, shifting the reference (docs/blackbox.md; seen on
// ENA ERR1892229). --printDS values from 2026-11-13. Default mode does not shift.
#[test]
fn reverse_read_lowercase_n_shift_only_in_pmdtools_compat() {
    let dir = tmp_dir("lowercase_n");
    let seq = "ACGTTGCAACGTTGCAACGTTGCAACGTAA";
    let close = |x: f64, y: f64| (x - y).abs() < 1e-12;
    let rev_n = one_read_sam(&dir.join("rev_n.sam"), 16, seq, "28n0n0");
    let rev_upper = one_read_sam(&dir.join("rev_upper.sam"), 16, seq, "28N0N0");
    let fwd_n = one_read_sam(&dir.join("fwd_n.sam"), 0, seq, "n0n28");
    let compat = ["--pmdtools-compat"];
    assert!(close(single_score(&rev_n, &compat), 20.935623477213277));
    assert!(close(
        single_score(&rev_upper, &compat),
        -0.5514281743729341
    ));
    assert!(close(single_score(&rev_n, &[]), -0.5514281743729341));
    assert!(close(single_score(&fwd_n, &compat), -0.30397668401849764));
}

// Departure from PMDtools (docs/design.md §5): a read with no informative sites scores 0 but
// carries no evidence, so --upper-threshold does not keep it. PMDtools (and --pmdtools-compat)
// keeps it.
#[test]
fn upper_threshold_drops_reads_without_evidence_except_in_compat() {
    let dir = tmp_dir("upper_threshold");
    let input = dir.join("two.sam");
    let (a, q) = ("A".repeat(30), "I".repeat(30));
    let c = format!("C{}", "A".repeat(29));
    fs::write(
        &input,
        format!(
            "@SQ\tSN:r\tLN:1000\nno_sites\t0\tr\t100\t60\t30M\t*\t0\t0\t{a}\t{q}\tMD:Z:30\n\
             one_match\t0\tr\t200\t60\t30M\t*\t0\t0\t{c}\t{q}\tMD:Z:30\n"
        ),
    )
    .unwrap();
    let kept = |extra: &[&str]| -> Vec<String> {
        let out = dir.join("kept.sam");
        let mut args = vec!["filter", "--quiet", "--upper-threshold", "1"];
        args.extend_from_slice(extra);
        args.extend([input.to_str().unwrap(), "-o", out.to_str().unwrap()]);
        let o = run(&args);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        read_records(&out).1.iter().map(name).collect()
    };
    assert_eq!(kept(&[]), ["one_match"]);
    assert_eq!(kept(&["--pmdtools-compat"]), ["no_sites", "one_match"]);
}

// BAM records are read through a BAM-specific path (reads::process_bam); it must give the same
// results as the generic path used for SAM.
#[test]
fn bam_and_sam_input_give_identical_scores_and_damage() {
    let dir = tmp_dir("bam_vs_sam");
    let sam = data("reads_random.sam");
    let bam = dir.join("all.bam");
    // Every scored read passes this threshold, so the BAM holds all reads the SAM scores.
    let o = run(&[
        "filter",
        "--quiet",
        "-t",
        "-1e9",
        sam.to_str().unwrap(),
        "-o",
        bam.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for mode in [&[][..], &["--pmdtools-compat"], &["--ss"], &["--UDGhalf"]] {
        let scores = |input: &Path| {
            let mut args = vec!["score", "--quiet"];
            args.extend_from_slice(mode);
            args.push(input.to_str().unwrap());
            let o = run(&args);
            assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
            o.stdout
        };
        assert_eq!(scores(&sam), scores(&bam), "scores differ for {mode:?}");
        let table = |input: &Path| {
            let mut args = vec!["damage", "--quiet"];
            args.extend_from_slice(mode);
            args.push(input.to_str().unwrap());
            let o = run(&args);
            assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
            let json: serde_json::Value = serde_json::from_slice(&o.stdout).unwrap();
            (json["damage"].clone(), json["reads"]["scored"].clone())
        };
        assert_eq!(
            table(&sam),
            table(&bam),
            "damage table differs for {mode:?}"
        );
    }
}

/// A SAM file of `n` 30-bp reads (several batches of 10,000), with optional broken lines:
/// `no_md` drops the MD tag of that read; `garbage` replaces that line with an unparsable one.
fn many_reads_sam(path: &Path, n: usize, no_md: Option<usize>, garbage: Option<usize>) {
    let (seq, qual) = (format!("T{}", "A".repeat(29)), "I".repeat(30));
    let mut sam = String::from("@HD\tVN:1.6\n@SQ\tSN:c\tLN:100000\n");
    for i in 1..=n {
        let md = if Some(i) == no_md { "" } else { "\tMD:Z:0C29" };
        if Some(i) == garbage {
            sam.push_str("garbage\tline\n");
        } else {
            sam.push_str(&format!(
                "r{i}\t0\tc\t{i}\t37\t30M\t*\t0\t0\t{seq}\t{qual}{md}\n"
            ));
        }
    }
    fs::write(path, sam).unwrap();
}

// With several threads, reading, scoring and output overlap (main.rs scan). Errors must still be
// reported in input order, stop the run with exit 2 and leave no output, as with one thread.
#[test]
fn errors_are_reported_in_input_order_with_any_thread_count() {
    let dir = tmp_dir("pipeline_errors");
    let both = dir.join("both.sam");
    many_reads_sam(&both, 25_000, Some(12_000), Some(21_000));
    let parse_only = dir.join("parse.sam");
    many_reads_sam(&parse_only, 25_000, None, Some(21_000));
    for threads in ["1", "2", "4"] {
        for (input, expect) in [(&both, "read r12000"), (&parse_only, "record 21000")] {
            let out = dir.join("out.bam");
            let o = run(&[
                "filter",
                "-t",
                "0",
                "-@",
                threads,
                input.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]);
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "-@ {threads}: {err}");
            assert!(
                err.contains(expect),
                "-@ {threads}: expected {expect:?} in {err}"
            );
            assert!(
                !out.exists(),
                "-@ {threads}: no output under the final name"
            );
        }
    }
}

#[test]
fn multi_batch_results_do_not_depend_on_thread_count() {
    let dir = tmp_dir("pipeline_batches");
    let input = dir.join("many.sam");
    many_reads_sam(&input, 25_000, None, None);
    let summary = |threads: &str, extra: &[&str]| {
        let mut args = vec!["damage", "--quiet", "-@", threads];
        args.extend_from_slice(extra);
        args.push(input.to_str().unwrap());
        let o = run(&args);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let mut json: serde_json::Value = serde_json::from_slice(&o.stdout).unwrap();
        json.as_object_mut().unwrap().remove("command");
        json
    };
    for extra in [&[][..], &["-n", "15000"]] {
        let one = summary("1", extra);
        assert_eq!(one, summary("2", extra), "{extra:?}");
        assert_eq!(one, summary("4", extra), "{extra:?}");
        let expected = if extra.is_empty() { 25_000 } else { 15_000 };
        assert_eq!(one["reads"]["input"], expected);
    }
    let scores =
        |threads: &str| run(&["score", "--quiet", "-@", threads, input.to_str().unwrap()]).stdout;
    assert_eq!(scores("1"), scores("2"));
    assert_eq!(scores("1"), scores("4"));
}

/// Damage tables are only counted when a summary is written; scores, kept reads and the summary's
/// tables must not depend on that.
#[test]
fn results_do_not_depend_on_summary() {
    let dir = tmp_dir("summary_or_not");
    let input = data("reads_random.sam");
    let input = input.to_str().unwrap();
    let from_damage = damage_json(&[])["damage"].clone();
    assert_ne!(from_damage["first_position"]["n"], 0);
    for threads in ["1", "2"] {
        let js = dir.join(format!("score{threads}.json"));
        let with = run(&[
            "score",
            "--quiet",
            "-@",
            threads,
            "--summary",
            js.to_str().unwrap(),
            input,
        ]);
        let without = run(&["score", "--quiet", "-@", threads, input]);
        assert!(with.status.success() && without.status.success());
        assert!(
            with.stdout == without.stdout,
            "score output depends on --summary (-@ {threads})"
        );
        let j: serde_json::Value = serde_json::from_str(&fs::read_to_string(&js).unwrap()).unwrap();
        assert_eq!(
            j["damage"], from_damage,
            "score --summary damage table (-@ {threads})"
        );

        let (a, b) = (
            dir.join(format!("with{threads}.bam")),
            dir.join(format!("without{threads}.bam")),
        );
        let fjs = dir.join(format!("filter{threads}.json"));
        let fa = ["filter", "-t", "3", "--tag", "--quiet", "-@", threads];
        let o = run(&[
            &fa[..],
            &[
                "--summary",
                fjs.to_str().unwrap(),
                input,
                "-o",
                a.to_str().unwrap(),
            ],
        ]
        .concat());
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let o = run(&[&fa[..], &[input, "-o", b.to_str().unwrap()]].concat());
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let (ra, rb) = (read_records(&a).1, read_records(&b).1);
        assert!(
            !ra.is_empty() && ra == rb,
            "kept reads depend on --summary (-@ {threads})"
        );
        let j: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&fjs).unwrap()).unwrap();
        assert_eq!(
            j["damage"], from_damage,
            "filter --summary damage table (-@ {threads})"
        );
    }
}

/// `--range 0` would give empty damage tables; it is rejected before any input is read.
#[test]
fn range_must_be_at_least_1() {
    let dir = tmp_dir("range_0");
    let input = data("reads_random.sam");
    let input = input.to_str().unwrap();
    let out = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let (s, f, d) = (out("s.tsv"), out("f.bam"), out("d.json"));
    let cases = [
        vec!["score", "--range", "0", input, "-o", &s],
        vec!["filter", "-t", "3", "--range", "0", input, "-o", &f],
        vec!["damage", "--range", "0", input, "-o", &d],
    ];
    for args in cases {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(String::from_utf8_lossy(&o.stderr).contains("--range must be at least 1"));
    }
    assert_eq!(
        fs::read_dir(&dir).unwrap().count(),
        0,
        "no output file is left behind"
    );
}

#[test]
fn subcommand_help_groups_options_and_shows_examples() {
    for (cmd, headings) in [
        (
            "filter",
            &[
                "Output:",
                "Filtering:",
                "Summary:",
                "Scoring model:",
                "Read filters:",
                "Performance:",
            ][..],
        ),
        (
            "score",
            &[
                "Output:",
                "Summary:",
                "Scoring model:",
                "Read filters:",
                "Performance:",
            ][..],
        ),
        (
            "damage",
            &[
                "Input:",
                "Output:",
                "Summary:",
                "Scoring model:",
                "Read filters:",
                "Performance:",
            ][..],
        ),
    ] {
        let o = run(&[cmd, "--help"]);
        assert!(o.status.success() && o.stderr.is_empty());
        let help = String::from_utf8(o.stdout).unwrap();
        assert!(help.contains(&format!("Usage: pmdshodh {cmd} ")), "{help}");
        for h in headings.iter().chain(&["Examples:"]) {
            assert!(help.contains(&format!("\n{h}\n")), "{cmd}: no {h}");
        }
    }
}

/// Every `pmdshodh ...` line in the help examples runs (inputs replaced by a test fixture).
/// The samtools pipe example is skipped: tests do not depend on samtools.
#[test]
fn help_examples_run() {
    let dir = tmp_dir("help_examples");
    let input = data("reads_random.sam");
    let mut ran = 0;
    for cmd in ["filter", "score", "damage"] {
        let help = String::from_utf8(run(&[cmd, "--help"]).stdout).unwrap();
        let examples = help.split("\nExamples:\n").nth(1).unwrap();
        for line in examples.lines().map(str::trim) {
            let Some(rest) = line.strip_prefix("pmdshodh ") else {
                continue;
            };
            let args: Vec<String> = rest
                .split_whitespace()
                .map(|a| match a {
                    "in.bam" => input.to_str().unwrap().to_string(),
                    a if a.ends_with(".bam") || a.ends_with(".json") || a.ends_with(".tsv") => {
                        dir.join(a).to_str().unwrap().to_string()
                    }
                    a => a.to_string(),
                })
                .collect();
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let o = run(&args);
            assert!(
                o.status.success(),
                "{line}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            ran += 1;
        }
    }
    assert_eq!(ran, 6, "two runnable examples per command");
}

/// A synthetic reference (one 2,000-base contig, indexed) and SAM reads aligned to it with damage-like
/// C->T changes near the 5' end, mismatches, deletions, insertions and soft clips. MD tags are
/// written from the reference as samtools calmd writes them. Returns (fasta, sam).
fn reference_and_reads(dir: &Path) -> (PathBuf, PathBuf) {
    let mut x: u64 = 0x2545F4914F6CDD1D;
    let mut rand = |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let refseq: Vec<u8> = (0..2000).map(|_| b"ACGT"[rand(4) as usize]).collect();
    let fasta = dir.join("ref.fa");
    // The second half is soft-masked (lowercase); the MD tags below are uppercase, as calmd
    // writes them.
    let masked: Vec<u8> = refseq
        .iter()
        .enumerate()
        .map(|(i, b)| if i < 1000 { *b } else { b.to_ascii_lowercase() })
        .collect();
    fs::write(&fasta, [b">chrT\n".as_slice(), &masked, b"\n"].concat()).unwrap();
    fs::write(dir.join("ref.fa.fai"), "chrT\t2000\t6\t2000\t2001\n").unwrap();

    let mut sam = String::from("@HD\tVN:1.6\n@SQ\tSN:chrT\tLN:2000\n");
    for i in 0..400 {
        let pos = rand(1900) as usize;
        let shape = rand(4);
        // (CIGAR, read length, reference span)
        let (cigar, ops): (&str, &[(u8, usize)]) = match shape {
            0 => ("40M", &[(b'M', 40)]),
            1 => ("20M2D20M", &[(b'M', 20), (b'D', 2), (b'M', 20)]),
            2 => ("3S37M", &[(b'S', 3), (b'M', 37)]),
            _ => ("15M1I24M", &[(b'M', 15), (b'I', 1), (b'M', 24)]),
        };
        // MD as the SAM spec and samtools calmd write it: the match count before every mismatch
        // and every deletion (0 if none), the reference base(s), and the final count.
        let (mut seq, mut md, mut run, mut r) = (Vec::new(), String::new(), 0, pos);
        for &(op, n) in ops {
            match op {
                b'M' => {
                    for _ in 0..n {
                        let refb = refseq[r];
                        let mut b = refb;
                        if refb == b'C' && seq.len() < 3 && rand(2) == 0 {
                            b = b'T';
                        } else if rand(25) == 0 {
                            b = b"ACGT"[rand(4) as usize];
                        }
                        if b == refb {
                            run += 1;
                        } else {
                            md += &format!("{run}{}", refb as char);
                            run = 0;
                        }
                        seq.push(b);
                        r += 1;
                    }
                }
                b'D' => {
                    md += &format!("{run}^{}", String::from_utf8_lossy(&refseq[r..r + n]));
                    run = 0;
                    r += n;
                }
                _ => seq.extend((0..n).map(|_| b"ACGT"[rand(4) as usize])),
            }
        }
        md += &run.to_string();
        let flag = if i % 2 == 0 { 0 } else { 16 };
        let qual = "I".repeat(seq.len());
        sam += &format!(
            "r{i}\t{flag}\tchrT\t{}\t60\t{cigar}\t*\t0\t0\t{}\t{qual}\tMD:Z:{md}\n",
            pos + 1,
            String::from_utf8(seq).unwrap()
        );
    }
    let path = dir.join("reads.sam");
    fs::write(&path, sam).unwrap();
    (fasta, path)
}

/// `-r` gives the same scores as correct MD tags (default and --pmdtools-compat, 1 and 2
/// threads), scores reads without MD, counts stale MD tags, and rejects a bad reference before
/// reading any record.
#[test]
fn reference_scores_equal_md_scores_and_bad_references_are_rejected() {
    let dir = tmp_dir("reference_r");
    let (fasta, sam) = reference_and_reads(&dir);
    let (fasta, sam_s) = (fasta.to_str().unwrap(), sam.to_str().unwrap());
    let score = |extra: &[&str], input: &str| {
        let mut args = vec!["score", "--quiet"];
        args.extend_from_slice(extra);
        args.push(input);
        let o = run(&args);
        assert!(
            o.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        o.stdout
    };
    for compat in [&[][..], &["--pmdtools-compat"][..]] {
        let from_md = score(compat, sam_s);
        for threads in ["1", "2"] {
            let with_r = score(&[compat, &["-r", fasta, "-@", threads]].concat(), sam_s);
            assert!(
                with_r == from_md,
                "{compat:?} -@ {threads}: -r differs from MD"
            );
        }
    }
    let lines = String::from_utf8(score(&[], sam_s)).unwrap();
    assert_eq!(lines.lines().count(), 401);
    assert!(
        lines
            .lines()
            .skip(1)
            .any(|l| l.split('\t').nth(2) != Some("0")),
        "no read has a score"
    );

    // Without MD tags: -r scores them the same; without -r the error names -r.
    let text = fs::read_to_string(&sam).unwrap();
    let no_md: String = text
        .lines()
        .map(|l| match l.find("\tMD:Z:") {
            Some(i) => &l[..i],
            None => l,
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let no_md_path = dir.join("no_md.sam");
    fs::write(&no_md_path, no_md).unwrap();
    let no_md_s = no_md_path.to_str().unwrap();
    assert!(score(&["-r", fasta], no_md_s) == score(&[], sam_s));
    let o = run(&["score", "--quiet", no_md_s]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("-r"));

    // A stale MD tag (first read's MD shifted by one base) is scored from the reference and counted.
    let stale: String = text.replacen("\tMD:Z:", "\tMD:Z:1A", 1);
    let stale_path = dir.join("stale.sam");
    fs::write(&stale_path, stale).unwrap();
    let js = dir.join("stale.json");
    let o = run(&[
        "score",
        "-r",
        fasta,
        "--summary",
        js.to_str().unwrap(),
        stale_path.to_str().unwrap(),
        "-o",
        dir.join("stale.tsv").to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("1 scored reads have an MD tag that disagrees")
    );
    assert!(fs::read(dir.join("stale.tsv")).unwrap() == score(&[], sam_s));
    let j: serde_json::Value = serde_json::from_str(&fs::read_to_string(&js).unwrap()).unwrap();
    assert_eq!(j["reference"]["bases_from"], "reference");
    assert_eq!(j["reference"]["md_tags"]["disagree"], 1);
    assert_eq!(j["reference"]["md_tags"]["agree"], 399);
    assert_eq!(j["reference"]["md_tags"]["absent"], 0);

    // Bad references: exit 2 before any output is written.
    let other = dir.join("other.fa");
    fs::write(&other, ">chrT\nACGT\n").unwrap();
    fs::write(dir.join("other.fa.fai"), "chrT\t4\t6\t4\t5\n").unwrap();
    let noindex = dir.join("noindex.fa");
    fs::write(&noindex, ">chrT\nACGT\n").unwrap();
    for (reference, message) in [
        (other.to_str().unwrap(), "does not match the input header"),
        (noindex.to_str().unwrap(), "samtools faidx"),
        ("missing.fa", "no such file"),
    ] {
        let out = dir.join("bad.bam");
        let o = run(&[
            "filter",
            "-t",
            "3",
            "-r",
            reference,
            sam_s,
            "-o",
            out.to_str().unwrap(),
        ]);
        assert_eq!(o.status.code(), Some(2), "{reference}");
        assert!(
            String::from_utf8_lossy(&o.stderr).contains(message),
            "{reference}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert!(!out.exists());
    }
}

/// Write `sam` as CRAM, encoded against an uppercase copy of `fasta` (the CRAM checksums do not
/// depend on case).
fn write_cram(dir: &Path, fasta: &Path, sam: &Path) -> PathBuf {
    use noodles::fasta;
    use noodles::sam::alignment::io::Write as _;
    let upper = dir.join("upper.fa");
    fs::write(&upper, fs::read(fasta).unwrap().to_ascii_uppercase()).unwrap();
    let mut fai = fasta.as_os_str().to_owned();
    fai.push(".fai");
    fs::copy(&fai, dir.join("upper.fa.fai")).unwrap();
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(&upper)
        .unwrap();
    let repository =
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader));
    let (header, records) = read_records(sam);
    let path = dir.join("reads.cram");
    let mut w = noodles::cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&path)
        .unwrap();
    w.write_header(&header).unwrap();
    for r in &records {
        w.write_alignment_record(&header, r).unwrap();
    }
    w.try_finish(&header).unwrap();
    path
}

/// CRAM input, decoded with -r against a soft-masked reference, gives the same scores and the same
/// kept records as the SAM it was written from, for 1 and 2 threads. CRAM without -r, or with a
/// reference that differs under the reads (slice MD5), is an input error, and the summary records
/// the input format.
#[test]
fn cram_input_matches_its_sam_and_needs_the_reference() {
    let dir = tmp_dir("cram_input");
    let (fasta, sam) = reference_and_reads(&dir);
    let cram = write_cram(&dir, &fasta, &sam);
    let (fasta, sam, cram) = (
        fasta.to_str().unwrap(),
        sam.to_str().unwrap(),
        cram.to_str().unwrap(),
    );
    let ok = |args: &[&str]| {
        let o = run(args);
        assert!(
            o.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        o.stdout
    };
    for compat in [&[][..], &["--pmdtools-compat"][..]] {
        let from_sam = ok(&[&["score", "--quiet", "-r", fasta], compat, &[sam]].concat());
        for threads in ["1", "2"] {
            let from_cram = ok(&[
                &["score", "--quiet", "-r", fasta, "-@", threads],
                compat,
                &[cram],
            ]
            .concat());
            assert!(
                from_cram == from_sam,
                "{compat:?} -@ {threads}: CRAM scores differ from SAM"
            );
        }
    }
    // Kept records are the SAM records (sequences uppercase although the reference is soft-masked).
    let records = |input: &str| {
        let out = ok(&[
            "filter",
            "--quiet",
            "-r",
            fasta,
            "-t",
            "1",
            "--output-format",
            "sam",
            input,
            "-o",
            "-",
        ]);
        let text = String::from_utf8(out).unwrap();
        text.lines()
            .filter(|l| !l.starts_with('@'))
            .map(str::to_string)
            .collect::<Vec<_>>()
    };
    let kept = records(sam);
    assert!(!kept.is_empty() && kept.len() < 400, "{} kept", kept.len());
    assert_eq!(records(cram), kept);

    // CRAMs written by samtools 1.25 (htslib 1.25) from these same reads.sam and ref.fa, in a
    // directory holding just those files:
    //   samtools view --no-PG -C --output-fmt-option version=3.0 -T ref.fa -o reads_htslib_3.0.cram reads.sam
    // and likewise for 3.1. htslib does not store an MD tag it can regenerate from the reference,
    // and pmdshodh does not regenerate it, so kept records are the SAM records without MD.
    let without_md: Vec<String> = kept
        .iter()
        .map(|l| {
            l.split('\t')
                .filter(|f| !f.starts_with("MD:Z:"))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    for version in ["3.0", "3.1"] {
        let htslib = data(&format!("reads_htslib_{version}.cram"));
        let htslib = htslib.to_str().unwrap();
        for compat in [&[][..], &["--pmdtools-compat"][..]] {
            let from_sam = ok(&[&["score", "--quiet", "-r", fasta], compat, &[sam]].concat());
            let from_htslib = ok(&[&["score", "--quiet", "-r", fasta], compat, &[htslib]].concat());
            assert!(
                from_htslib == from_sam,
                "CRAM {version} {compat:?}: scores differ from SAM"
            );
        }
        assert_eq!(records(htslib), without_md, "CRAM {version}");
    }

    let o = run(&["score", "--quiet", cram]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("needs its reference"));
    assert!(o.stdout.is_empty());

    // A reference with the same names and lengths but one base changed under the reads.
    let wrong = dir.join("wrong.fa");
    let mut text = fs::read(fasta).unwrap();
    let at = ">chrT\n".len() + 500;
    text[at] = if text[at] == b'A' { b'C' } else { b'A' };
    fs::write(&wrong, text).unwrap();
    fs::copy(format!("{fasta}.fai"), dir.join("wrong.fa.fai")).unwrap();
    // Caught by the header M5 check before noodles' slice MD5 check; with several threads the
    // check runs on the pool.
    for threads in ["1", "2"] {
        let o = run(&[
            "score",
            "--quiet",
            "-@",
            threads,
            "-r",
            wrong.to_str().unwrap(),
            cram,
        ]);
        assert_eq!(o.status.code(), Some(2), "-@ {threads}");
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("M5 checksum"),
            "-@ {threads}"
        );
    }

    // --max-reads stops at the same read in CRAM (counted from container headers) as in SAM.
    let mut summaries = Vec::new();
    for (input, format) in [(cram, "cram"), (sam, "sam")] {
        let js = dir.join(format!("{format}.json"));
        ok(&[
            "damage",
            "--quiet",
            "-r",
            fasta,
            "--max-reads",
            "150",
            "-@",
            "2",
            input,
            "-o",
            js.to_str().unwrap(),
        ]);
        let j: serde_json::Value = serde_json::from_str(&fs::read_to_string(&js).unwrap()).unwrap();
        assert_eq!(j["input_format"], format);
        assert_eq!(j["reads"]["input"], 150);
        summaries.push((
            j["reads"].clone(),
            j["damage"].clone(),
            j["score_histogram"].clone(),
        ));
    }
    assert_eq!(summaries[0], summaries[1]);
}

/// A samtools CRAM with one multi-reference slice (two contigs), a proper pair, both strands and
/// an unmapped read. Files: tests/data/pairs.{fa,fa.fai,sam,cram}; the SAM is synthetic (no MD
/// tags), and the CRAM was written by samtools 1.25 in a directory holding just those files:
///   samtools view --no-PG -C --output-fmt-option version=3.1 \
///     --output-fmt-option multi_seq_per_slice=1 -T pairs.fa -o pairs.cram pairs.sam
/// htslib decodes it back to pairs.sam (`--input-fmt-option decode_md=0`). pairs_empty.cram is
/// the same command on the SAM's header lines only.
///
/// Scores and kept records (mate fields included) equal those from the SAM. Multi-reference
/// slices carry no reference MD5, so a reference with one base changed under a read was decoded
/// silently (wrong read bases and scores, exit 0); it must be rejected via the header's M5.
#[test]
fn cram_multi_reference_slice_pairs_and_wrong_reference() {
    let dir = tmp_dir("cram_pairs");
    let (fasta, sam, cram) = (data("pairs.fa"), data("pairs.sam"), data("pairs.cram"));
    let (fasta, sam, cram) = (
        fasta.to_str().unwrap(),
        sam.to_str().unwrap(),
        cram.to_str().unwrap(),
    );
    let ok = |args: &[&str]| {
        let o = run(args);
        assert!(
            o.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        o.stdout
    };
    for threads in ["1", "2"] {
        let score = |input: &str| ok(&["score", "--quiet", "-@", threads, "-r", fasta, input]);
        let from_cram = score(cram);
        assert_eq!(from_cram, score(sam), "-@ {threads}");
        // Five mapped reads scored; the unmapped one is skipped.
        assert_eq!(String::from_utf8(from_cram).unwrap().lines().count(), 6);
        let kept = |input: &str| {
            let out = ok(&[
                "filter",
                "--quiet",
                "-@",
                threads,
                "-r",
                fasta,
                "-t",
                "-1000",
                input,
                "-o",
                "-",
                "--output-format",
                "sam",
            ]);
            String::from_utf8(out)
                .unwrap()
                .lines()
                .filter(|l| !l.starts_with('@'))
                .map(str::to_string)
                .collect::<Vec<_>>()
        };
        let from_cram = kept(cram);
        assert_eq!(from_cram.len(), 5);
        assert_eq!(from_cram, kept(sam), "-@ {threads}");
    }

    // One base changed on c2, inside read s1 (c2:11-50).
    let wrong = dir.join("wrong.fa");
    let mut text = fs::read(fasta).unwrap();
    let at = text.windows(4).position(|w| w == b">c2\n").unwrap() + 4 + 20;
    text[at] = if text[at] == b'A' { b'C' } else { b'A' };
    fs::write(&wrong, text).unwrap();
    fs::copy(format!("{fasta}.fai"), dir.join("wrong.fa.fai")).unwrap();
    let out = dir.join("wrong.bam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        "--quiet",
        "-r",
        wrong.to_str().unwrap(),
        cram,
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("M5 checksum"));
    assert!(!out.exists());

    // An empty CRAM (header only; same samtools command from the SAM header) gives a valid,
    // empty BAM and exit 0.
    let empty = data("pairs_empty.cram");
    let out = dir.join("empty.bam");
    ok(&[
        "filter",
        "--quiet",
        "-r",
        fasta,
        "-t",
        "3",
        empty.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let (header, records) = read_records(&out);
    assert_eq!(header.reference_sequences().len(), 2);
    assert!(records.is_empty());
}

/// `score -o FILE` leaves no file when the run fails part-way, as filter does; a gzip (not BGZF)
/// input is an input error that names bgzip.
#[test]
fn score_output_is_atomic_and_plain_gzip_is_named() {
    let dir = tmp_dir("score_atomic");
    let sam = dir.join("bad.sam");
    many_reads_sam(&sam, 25_000, None, Some(20_001));
    let out = dir.join("scores.tsv");
    let o = run(&["score", sam.to_str().unwrap(), "-o", out.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2));
    assert!(!out.exists());
    assert!(fs::read_dir(&dir).unwrap().all(|e| {
        !e.unwrap()
            .file_name()
            .to_string_lossy()
            .contains("pmdshodh-tmp")
    }));
    // A good run still writes the file.
    let good = dir.join("good.sam");
    many_reads_sam(&good, 10, None, None);
    let o = run(&[
        "score",
        "--quiet",
        good.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success());
    assert_eq!(fs::read_to_string(&out).unwrap().lines().count(), 11);

    // gzip header (FLG without FEXTRA) in front of SAM text.
    let gz = dir.join("reads.sam.gz");
    fs::write(&gz, [0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0x03]).unwrap();
    let o = run(&["score", gz.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("not BGZF: recompress it with bgzip"));
}

/// `=` read bases (same as the reference; `samtools calmd -e` writes them). PMDtools 0.60 never
/// counts them as C, so a read whose matching Cs are `=` scored 5.0297 from 1 site, not 4.2201
/// from 10 like the same read spelled out (both PMDtools values, measured 2026-11-16). Default
/// mode now resolves them from -r and rejects them without -r (MD does not give their reference
/// bases); --pmdtools-compat keeps PMDtools' value. Kept records still carry the `=` bases.
#[test]
fn equals_read_bases_are_resolved_from_the_reference() {
    let dir = tmp_dir("equals_bases");
    let fasta = dir.join("ref.fa");
    fs::write(&fasta, ">c\nCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n").unwrap();
    fs::write(dir.join("ref.fa.fai"), "c\t40\t3\t40\t41\n").unwrap();
    let sam = |name: &str, seq: &str| {
        let path = dir.join(name);
        fs::write(
            &path,
            format!(
                "@SQ\tSN:c\tLN:40\nr\t0\tc\t1\t60\t10M\t*\t0\t0\t{seq}\tIIIIIIIIII\tMD:Z:0C9\n"
            ),
        )
        .unwrap();
        path
    };
    let (eq, spelled) = (
        sam("eq.sam", "T========="),
        sam("spelled.sam", "TCCCCCCCCC"),
    );
    let (eq, spelled, fasta) = (
        eq.to_str().unwrap(),
        spelled.to_str().unwrap(),
        fasta.to_str().unwrap(),
    );
    let score = |args: &[&str]| -> (f64, u32) {
        let o = run(&[&["score", "--quiet"], args].concat());
        assert!(
            o.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        let text = String::from_utf8(o.stdout).unwrap();
        let f: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        (f[2].parse().unwrap(), f[3].parse().unwrap())
    };
    let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
    let (s, n) = score(&["-r", fasta, eq]);
    assert!(close(s, 4.2200864108613265) && n == 10, "{s} {n}");
    assert_eq!(score(&["-r", fasta, eq]), score(&[spelled]));
    for args in [
        &["--pmdtools-compat", eq][..],
        &["--pmdtools-compat", "-r", fasta, eq][..],
    ] {
        let (s, n) = score(args);
        assert!(close(s, 5.029717285061266) && n == 1, "{args:?}: {s} {n}");
    }
    let o = run(&["score", "--quiet", eq]);
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("pass the reference with -r"));

    let o = run(&[
        "filter",
        "--quiet",
        "-r",
        fasta,
        "-t",
        "3",
        eq,
        "-o",
        "-",
        "--output-format",
        "sam",
    ]);
    assert!(o.status.success());
    let text = String::from_utf8(o.stdout).unwrap();
    assert!(
        text.lines()
            .any(|l| l.starts_with("r\t") && l.contains("\tT=========\t"))
    );
}

/// A BAM whose blocks are whole but whose BGZF end-of-file marker is missing (e.g. a copy cut at
/// a block boundary) used to be read silently, losing the reads after the cut. It is now an input
/// error, and filter leaves no output file.
#[test]
fn bam_without_eof_marker_is_an_error() {
    let dir = tmp_dir("bam_no_eof");
    let bam = dir.join("all.bam");
    let o = run(&[
        "filter",
        "--quiet",
        "-t",
        "-1000",
        data("reads_random.sam").to_str().unwrap(),
        "-o",
        bam.to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let bytes = fs::read(&bam).unwrap();
    let cut = dir.join("cut.bam");
    fs::write(&cut, &bytes[..bytes.len() - 28]).unwrap();
    let out = dir.join("out.bam");
    for args in [
        &["score", "--quiet", cut.to_str().unwrap(), "-o", "-"][..],
        &[
            "filter",
            "-t",
            "3",
            cut.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ][..],
    ] {
        let o = run(args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("without the BGZF end-of-file marker"),
            "{args:?}"
        );
    }
    assert!(!out.exists());
}

// Edge cases found by the CLI audit of 2026-11-18 (CHANGELOG.md).

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// The output used to be renamed over the input, replacing it (exit 0): 2,912 reads became 34.
#[test]
fn output_equal_to_input_or_another_output_is_rejected() {
    let dir = tmp_dir("same_file");
    let input = dir.join("in.sam");
    fs::copy(data("reads_random.sam"), &input).unwrap();
    let before = fs::read(&input).unwrap();
    let dotted = dir.join(".").join("in.sam");
    for args in [
        vec!["filter", "-t", "3", s(&input), "-o", s(&input)],
        vec!["filter", "-t", "3", s(&input), "-o", s(&dotted)],
        vec!["score", s(&input), "-o", s(&input)],
        vec!["score", s(&input), "--summary", s(&input)],
        vec!["damage", s(&input), "-o", s(&input)],
    ] {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(
            stderr(&o).contains("is the input file"),
            "{args:?}: {}",
            stderr(&o)
        );
    }
    #[cfg(unix)]
    {
        let link = dir.join("hard.sam");
        fs::hard_link(&input, &link).unwrap();
        let o = run(&["filter", "-t", "3", s(&input), "-o", s(&link)]);
        assert_eq!(o.status.code(), Some(2), "{}", stderr(&o));
    }
    assert_eq!(fs::read(&input).unwrap(), before, "input untouched");
    let out = dir.join("out.bam");
    let o = run(&[
        "filter",
        "-t",
        "3",
        s(&input),
        "-o",
        s(&out),
        "--summary",
        s(&out),
    ]);
    assert_eq!(o.status.code(), Some(2));
    assert!(
        stderr(&o).contains("given for two outputs"),
        "{}",
        stderr(&o)
    );
    assert!(!out.exists());
}

/// Writing under a temporary name and renaming replaced a symlink output with a regular file
/// (for score since 0.2.1). The link is now followed: its target is replaced and the link kept.
#[cfg(unix)]
#[test]
fn symlink_output_replaces_its_target_and_keeps_the_link() {
    let dir = tmp_dir("symlink_output");
    let input = data("reads_random.sam");

    let real = dir.join("real.tsv");
    fs::write(&real, "old\n").unwrap();
    let link = dir.join("link.tsv");
    std::os::unix::fs::symlink(&real, &link).unwrap();
    let o = run(&["score", "--quiet", s(&input), "-o", s(&link)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(
        fs::symlink_metadata(&link)
            .unwrap()
            .file_type()
            .is_symlink()
    );
    assert!(
        fs::read_to_string(&real)
            .unwrap()
            .starts_with("read_name\t")
    );
}

/// Renaming replaced FIFOs (its reader got nothing) and, as root, devices (/dev/full became a
/// regular file). An existing non-regular output is now written in place.
#[cfg(unix)]
#[test]
fn fifo_output_is_written_in_place() {
    use std::os::unix::fs::FileTypeExt;
    let dir = tmp_dir("fifo_output");
    let input = data("reads_random.sam");
    // A FIFO without an extension: written in place, as BAM by default.
    let fifo = dir.join("pipe");
    assert!(
        Command::new("mkfifo")
            .arg(&fifo)
            .status()
            .unwrap()
            .success()
    );
    let reader = {
        let fifo = fifo.clone();
        std::thread::spawn(move || fs::read(fifo).unwrap())
    };
    let o = run(&["filter", "--quiet", "-t", "3", s(&input), "-o", s(&fifo)]);
    assert!(o.status.success(), "{}", stderr(&o));
    // Checked before joining: if the FIFO had been replaced, its reader would wait forever.
    assert!(
        fs::metadata(&fifo).unwrap().file_type().is_fifo(),
        "FIFO kept"
    );
    let got = reader.join().unwrap();
    assert!(got.starts_with(&[0x1f, 0x8b]), "BAM through the FIFO");
    assert!(fs::read_dir(&dir).unwrap().all(|e| {
        !e.unwrap()
            .file_name()
            .to_string_lossy()
            .contains("pmdshodh-tmp")
    }));
}

/// Values in range one by one could still make the damage rate at the read end, p + c, exceed 1,
/// giving NaN scores that filter silently dropped; a NaN threshold kept nothing.
#[test]
fn parameters_that_make_scores_undefined_are_rejected() {
    let input = data("reads_random.sam");
    let input = s(&input);
    for args in [
        vec!["score", "--pmd-p", "0.6", "--pmd-c", "0.5", input],
        vec!["score", "--pmd-p", "1", input],
        vec!["filter", "-t", "nan", input, "-o", "x.bam"],
        vec!["filter", "-t", "-inf", input, "-o", "x.bam"],
        vec!["filter", "--upper-threshold", "inf", input, "-o", "x.bam"],
    ] {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(o.stdout.is_empty(), "{args:?}");
    }
    let o = run(&["score", "--pmd-p", "0.6", "--pmd-c", "0.5", input]);
    assert!(
        stderr(&o).contains("damage rate at the read end"),
        "{}",
        stderr(&o)
    );
    // At most 1 is fine, and a negative value given with a space gets the range message.
    let o = run(&[
        "score", "--quiet", "--pmd-p", "0.5", "--pmd-c", "0.5", input,
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(!String::from_utf8_lossy(&o.stdout).contains("NaN"));
    let o = run(&["score", "--pmd-p", "-0.1", input]);
    assert_eq!(o.status.code(), Some(2));
    assert!(stderr(&o).contains("between 0 and 1"), "{}", stderr(&o));
}

/// --range 100000000 was killed for lack of memory; -@ 100000 tried to start threads for minutes.
#[test]
fn range_and_threads_have_upper_limits() {
    let input = data("reads_random.sam");
    let input = s(&input);
    for (args, msg) in [
        (
            vec!["score", "--range", "1001", input],
            "--range must be at most 1000",
        ),
        (
            vec!["damage", "--range", "100000000", input],
            "--range must be at most 1000",
        ),
        (
            vec!["score", "-@", "1025", input],
            "--threads must be at most 1024",
        ),
    ] {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
        assert!(stderr(&o).contains(msg), "{args:?}: {}", stderr(&o));
    }
    let o = run(&["damage", "--quiet", "--range", "1000", input]);
    assert!(o.status.success(), "{}", stderr(&o));
}

/// A bad --summary path used to fail only after the whole run, with the BAM already written.
#[test]
fn summary_path_is_created_before_any_read() {
    let dir = tmp_dir("summary_early");
    let input = data("reads_random.sam");
    let out = dir.join("out.bam");
    let tsv = dir.join("scores.tsv");
    let bad = Path::new("no").join("such").join("s.json");
    for args in [
        vec![
            "filter",
            "-t",
            "3",
            s(&input),
            "-o",
            s(&out),
            "--summary",
            s(&bad),
        ],
        vec!["score", s(&input), "-o", s(&tsv), "--summary", s(&bad)],
    ] {
        let o = run(&args);
        assert_eq!(o.status.code(), Some(1), "{args:?}");
        assert!(
            stderr(&o).contains(&format!("cannot create {}", bad.display())),
            "{}",
            stderr(&o)
        );
        assert!(
            !stderr(&o).contains("pmdshodh-tmp"),
            "error names the user's path"
        );
    }
    assert_eq!(
        fs::read_dir(&dir).unwrap().count(),
        0,
        "no output, no temporary files"
    );
}

/// 0 bytes is not an alignment file (samtools refuses it); it used to give exit 0 and an empty BAM.
#[test]
fn zero_byte_input_is_an_error() {
    let dir = tmp_dir("zero_bytes");
    let empty = dir.join("empty.bam");
    fs::write(&empty, b"").unwrap();
    let o = run(&["score", s(&empty)]);
    assert_eq!(o.status.code(), Some(2));
    assert!(stderr(&o).contains("is empty"), "{}", stderr(&o));
    let o = Command::new(BIN)
        .args(["filter", "-t", "3", "-", "-o", "-"])
        .stdin(std::process::Stdio::null())
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(
        stderr(&o).contains("standard input is empty"),
        "{}",
        stderr(&o)
    );
    assert!(o.stdout.is_empty());
}

/// filter failed late on a kept read naming a reference the header lacks, as "cannot write output:
/// invalid reference sequence name" with exit 1. It is now an input error naming the read. score
/// and damage still read such SAM (`samtools view in.bam | pmdshodh score -`, as with PMDtools).
#[test]
fn kept_reads_naming_references_missing_from_the_header_are_input_errors() {
    let dir = tmp_dir("unknown_reference");
    let head = "@HD\tVN:1.6\n@SQ\tSN:c\tLN:100\n";
    let read = |rname: &str, rnext: &str| {
        format!("r1\t1\t{rname}\t1\t37\t4M\t{rnext}\t1\t0\tTAAC\tIIII\tMD:Z:0C3\n")
    };
    for (body, msg) in [
        (
            format!("{head}{}", read("x", "*")),
            "read r1: its reference is not in the header",
        ),
        (
            format!("{head}{}", read("c", "y")),
            "read r1: its mate's reference is not in the header",
        ),
        (
            read("c", "*"),
            "read r1: its reference is not in the header",
        ),
    ] {
        let sam = dir.join("in.sam");
        fs::write(&sam, body).unwrap();
        let o = run(&["filter", "-t", "-100", s(&sam), "-o", "-"]);
        assert_eq!(o.status.code(), Some(2));
        assert!(stderr(&o).contains(msg), "{}", stderr(&o));
        // Not kept, not written: no error.
        let o = run(&["filter", "--quiet", "-t", "100", s(&sam), "-o", "-"]);
        assert!(o.status.success(), "{}", stderr(&o));
        for args in [
            vec!["score", "--quiet", s(&sam)],
            vec!["damage", "--quiet", s(&sam)],
        ] {
            let o = run(&args);
            assert!(o.status.success(), "{args:?}: {}", stderr(&o));
        }
    }
}

/// FASTQ, FASTA and binary data used to be parsed as SAM: "invalid header: invalid record" or
/// "record 1: unexpected EOL".
#[test]
fn input_that_is_not_alignments_is_named() {
    let dir = tmp_dir("not_alignments");
    let bytes: Vec<u8> = (0..=255).collect();
    for (content, msg) in [
        (b"@r1\nACGT\n+\nIIII\n".to_vec(), "looks like FASTQ"),
        (b">chr1\nACGT\n".to_vec(), "looks like FASTA"),
        (bytes, "binary data"),
    ] {
        let f = dir.join("in.dat");
        fs::write(&f, content).unwrap();
        let o = run(&["score", s(&f)]);
        assert_eq!(o.status.code(), Some(2));
        assert!(
            stderr(&o).contains("is not SAM, BAM or CRAM"),
            "{}",
            stderr(&o)
        );
        assert!(stderr(&o).contains(msg), "{}", stderr(&o));
    }
}

/// A BAM cut inside a BGZF block now says that the file is probably truncated.
#[test]
fn bam_cut_inside_a_block_says_it_is_truncated() {
    let dir = tmp_dir("bam_cut_block");
    let bam = dir.join("all.bam");
    let input = data("reads_random.sam");
    let o = run(&["filter", "--quiet", "-t", "-1000", s(&input), "-o", s(&bam)]);
    assert!(o.status.success());
    let bytes = fs::read(&bam).unwrap();
    let cut = dir.join("cut.bam");
    fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    let o = run(&["score", "--quiet", s(&cut), "-o", "-"]);
    assert_eq!(o.status.code(), Some(2));
    assert!(stderr(&o).contains("probably truncated"), "{}", stderr(&o));
}

/// `score | head` printed "cannot write output: Broken pipe" and exited 1. Like samtools (killed
/// by SIGPIPE), pmdshodh now stops quietly with status 141.
#[cfg(unix)]
#[test]
fn closed_stdout_ends_the_run_quietly_with_status_141() {
    use std::io::Read;
    let dir = tmp_dir("broken_pipe");
    let sam = dir.join("many.sam");
    many_reads_sam(&sam, 50_000, None, None);
    for args in [
        vec!["score", "--quiet", s(&sam)],
        vec![
            "filter",
            "--quiet",
            "-t",
            "-100",
            "-@",
            "2",
            "--output-format",
            "sam",
            s(&sam),
            "-o",
            "-",
        ],
        vec![
            "filter",
            "--quiet",
            "-t",
            "-100",
            "-@",
            "2",
            s(&sam),
            "-o",
            "-",
        ],
    ] {
        let mut child = Command::new(BIN)
            .args(&args)
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        let mut first = [0u8; 100];
        child.stdout.take().unwrap().read_exact(&mut first).unwrap();
        let o = child.wait_with_output().unwrap();
        assert_eq!(o.status.code(), Some(141), "{args:?}: {}", stderr(&o));
        assert!(o.stderr.is_empty(), "{args:?}: {}", stderr(&o));
    }
}

#[test]
fn output_format_from_extension_in_any_case() {
    let dir = tmp_dir("output_format");
    let input = data("reads_random.sam");
    let upper = dir.join("out.BAM");
    let o = run(&["filter", "--quiet", "-t", "3", s(&input), "-o", s(&upper)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(fs::read(&upper).unwrap().starts_with(&[0x1f, 0x8b]));
    for (args, msg) in [
        (
            vec!["--output-format", "bam", "-o", "x.sam"],
            "contradicts the .sam extension",
        ),
        (vec!["-o", "x.cram"], "CRAM output is not supported"),
        (vec!["-o", "x.txt"], "use a .bam or .sam extension"),
        (
            vec!["--output-format", "sam", "-o", "x.txt"],
            "use a .bam or .sam extension",
        ),
    ] {
        let mut full = vec!["filter", "-t", "3", s(&input)];
        full.extend(args.iter());
        let o = run(&full);
        assert_eq!(o.status.code(), Some(2), "{full:?}");
        assert!(stderr(&o).contains(msg), "{full:?}: {}", stderr(&o));
    }
}

/// On Windows, NUL has no file metadata and NUL.pmdshodh-tmp is the device again, so renaming
/// failed ("cannot move the file to a different disk drive"). Device names are written in place.
#[cfg(windows)]
#[test]
fn windows_device_names_are_written_in_place() {
    let input = data("reads_random.sam");
    for args in [
        vec!["score", "--quiet", s(&input), "-o", "NUL"],
        vec![
            "filter",
            "--quiet",
            "-t",
            "3",
            s(&input),
            "-o",
            "nul",
            "--summary",
            "NUL",
        ],
        vec!["damage", "--quiet", s(&input), "-o", "NUL"],
    ] {
        let o = run(&args);
        assert!(o.status.success(), "{args:?}: {}", stderr(&o));
    }
}

/// BAM stores a CIGAR of more than 65,535 operations in a CG tag. Read from BAM and written as BAM
/// or SAM, such a record gained a second CG tag (264,008 bytes larger for 66,000 operations), so
/// samtools showed CG as an ordinary tag (found 2026-11-18). Kept records must be unchanged.
#[test]
fn long_cigar_records_keep_a_single_cg_tag() {
    use std::io::Read;
    let dir = tmp_dir("long_cigar");
    let ops = 33_000;
    let seq: String = "AC".repeat(ops);
    let sam = dir.join("long.sam");
    fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\n@SQ\tSN:c\tLN:100000\nlong\t0\tc\t1\t37\t{}\t*\t0\t0\t{seq}\t{}\tMD:Z:{ops}\n",
            "1M1I".repeat(ops),
            "I".repeat(2 * ops)
        ),
    )
    .unwrap();
    let cg_tags = |bam: &Path| {
        let mut bytes = Vec::new();
        noodles::bgzf::io::Reader::new(fs::File::open(bam).unwrap())
            .read_to_end(&mut bytes)
            .unwrap();
        bytes.windows(4).filter(|w| w == b"CGBI").count()
    };
    let a = dir.join("a.bam");
    let o = run(&["filter", "--quiet", "-t", "-1000000", s(&sam), "-o", s(&a)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(cg_tags(&a), 1, "SAM to BAM");
    let (b, c) = (dir.join("b.bam"), dir.join("c.sam"));
    for out in [&b, &c] {
        let o = run(&["filter", "--quiet", "-t", "-1000000", s(&a), "-o", s(out)]);
        assert!(o.status.success(), "{}", stderr(&o));
    }
    assert_eq!(cg_tags(&b), 1, "BAM to BAM");
    let line = fs::read_to_string(&c).unwrap();
    let line = line.lines().find(|l| l.starts_with("long\t")).unwrap();
    assert_eq!(line.split('\t').nth(5).unwrap(), "1M1I".repeat(ops));
    assert!(!line.contains("CG:B"), "BAM to SAM: no CG tag in SAM");
}

/// A run ended by SIGINT, SIGTERM or SIGHUP (Ctrl-C, a scheduler's time limit) left its
/// `.pmdshodh-tmp` files behind. They are now removed, and the process still ends by the signal.
#[cfg(unix)]
#[test]
fn a_signal_mid_run_removes_temporary_files() {
    use std::io::Write;
    use std::os::unix::process::ExitStatusExt;
    let mut head = String::from("@HD\tVN:1.6\n@SQ\tSN:c\tLN:100000\n");
    for i in 1..=50 {
        head.push_str(&format!(
            "r{i}\t0\tc\t{i}\t37\t4M\t*\t0\t0\tTAAC\tIIII\tMD:Z:0C3\n"
        ));
    }
    for (signal, name) in [(15, "TERM"), (2, "INT"), (1, "HUP")] {
        let dir = tmp_dir(&format!("signal_{name}"));
        let out = dir.join("out.bam");
        let summary = dir.join("run.json");
        let mut child = Command::new(BIN)
            .args([
                "filter",
                "-t",
                "3",
                "-",
                "-o",
                s(&out),
                "--summary",
                s(&summary),
            ])
            .stdin(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        // The header and a few reads, then stdin stays open: the run waits with its outputs open.
        let mut stdin = child.stdin.take().unwrap();
        stdin.write_all(head.as_bytes()).unwrap();
        stdin.flush().unwrap();
        let temporary = || {
            fs::read_dir(&dir)
                .unwrap()
                .filter(|e| {
                    e.as_ref()
                        .unwrap()
                        .file_name()
                        .to_string_lossy()
                        .ends_with("pmdshodh-tmp")
                })
                .count()
        };
        let start = std::time::Instant::now();
        while temporary() < 2 && start.elapsed().as_secs() < 20 {
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        assert_eq!(temporary(), 2, "both outputs open before the signal");
        let pid = child.id().to_string();
        assert!(
            Command::new("kill")
                .args([&format!("-{name}"), &pid])
                .status()
                .unwrap()
                .success()
        );
        let status = child.wait().unwrap();
        drop(stdin);
        assert_eq!(status.signal(), Some(signal), "SIG{name} ends the process");
        assert_eq!(temporary(), 0, "SIG{name}: temporary files removed");
        assert!(
            !out.exists() && !summary.exists(),
            "SIG{name}: no output under the final names"
        );
    }
}
