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
