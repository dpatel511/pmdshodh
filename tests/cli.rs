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
    // UDG-half tables are not compared: PMDtools' rule for them is not yet known (docs/blackbox.md 1.6d).
    for (treatment, flag) in [("none", "--UDGminus"), ("full", "--CpG")] {
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
    assert_eq!(j["format_version"], 1);
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
