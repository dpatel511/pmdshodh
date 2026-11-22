//! Properties from the 0.3 audit (docs/audit-0.3.md §5), checked on generated inputs: a random
//! reference and up to 25,000 reads (more than two 10,000-record batches) with damage,
//! mismatches, indels, clips, strands, qualities, MAPQ and flags, and MD tags written from the
//! reference.
//!
//! Each case runs the binary several times, so the number of cases is small (PROPTEST_CASES
//! overrides it). The reads per property add up to well over 10,000. PROP_MAX_N caps the reads
//! per case, for runs of many small cases (docs/audit-0.3.md §5.1).

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use proptest::prelude::*;

const BIN: &str = env!("CARGO_BIN_EXE_pmdshodh");
const REF_LEN: usize = 5_000;

/// xorshift64*: the reads come from a seed, so a failing case is reproduced from its seed alone.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, percent: usize) -> bool {
        self.below(100) < percent
    }
}

struct Input {
    dir: PathBuf,
    sam: PathBuf,
    fasta: PathBuf,
    /// SAM record lines, in order.
    records: Vec<String>,
}

/// Write a reference (FASTA and index) and `n` reads with correct MD tags.
fn generate(seed: u64, n: usize, name: &str) -> Input {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("prop_{name}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut rng = Rng(seed | 1);
    // CpG-rich enough for UDG modes, with a few N.
    let reference: Vec<u8> = (0..REF_LEN)
        .map(|_| match rng.below(100) {
            0 => b'N',
            k => b"ACGT"[k % 4],
        })
        .collect();
    let fasta = dir.join("ref.fa");
    let mut text = String::from(">c\n");
    for line in reference.chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    fs::write(&fasta, text).unwrap();
    fs::write(dir.join("ref.fa.fai"), format!("c\t{REF_LEN}\t3\t60\t61\n")).unwrap();

    let mut records = Vec::with_capacity(n);
    for i in 0..n {
        records.push(read(&mut rng, &reference, i));
    }
    let sam = dir.join("in.sam");
    let mut out = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:{REF_LEN}\n");
    for r in &records {
        out.push_str(r);
        out.push('\n');
    }
    fs::write(&sam, out).unwrap();
    Input {
        dir,
        sam,
        fasta,
        records,
    }
}

/// One read as a SAM line: CIGAR of M runs with insertions, deletions and soft clips; C→T damage
/// near the start and G→A near the end; random mismatches; MD written from the reference.
fn read(rng: &mut Rng, reference: &[u8], i: usize) -> String {
    let len = 20 + rng.below(131);
    let mut ops: Vec<(u8, usize)> = Vec::new();
    if rng.chance(10) {
        ops.push((b'S', 1 + rng.below(5)));
    }
    let mut aligned = 0;
    let target = len - ops.iter().map(|o| o.1).sum::<usize>() - if rng.chance(10) { 3 } else { 0 };
    while aligned < target {
        let m = (1 + rng.below(40)).min(target - aligned);
        ops.push((b'M', m));
        aligned += m;
        if aligned < target && rng.chance(15) {
            if rng.chance(50) {
                let k = (1 + rng.below(3)).min(target - aligned);
                ops.push((b'I', k));
                aligned += k;
            } else {
                ops.push((b'D', 1 + rng.below(3)));
            }
        }
    }
    // A deletion or insertion may not end the alignment.
    while matches!(ops.last(), Some((b'D', _)) | Some((b'I', _))) {
        ops.pop();
    }
    let read_len: usize = ops.iter().filter(|o| o.0 != b'D').map(|o| o.1).sum();
    if read_len < len {
        ops.push((b'S', len - read_len));
    }
    let span: usize = ops
        .iter()
        .filter(|o| matches!(o.0, b'M' | b'D'))
        .map(|o| o.1)
        .sum();
    let pos = rng.below(reference.len() - span);

    let mut seq = Vec::new();
    let mut md = String::new();
    let mut run = 0;
    let mut r = pos;
    for &(op, k) in &ops {
        match op {
            b'M' => {
                for _ in 0..k {
                    let refb = reference[r];
                    let mut b = refb;
                    let at = seq.len();
                    if refb == b'C' && at < 10 && rng.chance(40 - 4 * at) {
                        b = b'T';
                    } else if refb == b'G' && at + 10 > len && rng.chance(40) {
                        b = b'A';
                    } else if rng.chance(1) {
                        b = b"ACGT"[rng.below(4)];
                    }
                    if b == refb && refb != b'N' {
                        run += 1;
                    } else {
                        let _ = write!(md, "{run}{}", refb as char);
                        run = 0;
                    }
                    seq.push(b);
                    r += 1;
                }
            }
            b'D' => {
                let _ = write!(md, "{run}^");
                md.push_str(std::str::from_utf8(&reference[r..r + k]).unwrap());
                run = 0;
                r += k;
            }
            _ => seq.extend((0..k).map(|_| b"ACGT"[rng.below(4)])),
        }
    }
    let _ = write!(md, "{run}");
    let cigar: String = ops
        .iter()
        .map(|(op, k)| format!("{k}{}", *op as char))
        .collect();
    let qual: String = (0..seq.len())
        .map(|_| (b'#' + rng.below(39) as u8) as char)
        .collect();
    let mut flag = 0;
    if rng.chance(50) {
        flag |= 16;
    }
    if rng.chance(10) {
        flag |= 1024;
    }
    if rng.chance(5) {
        flag |= 256;
    }
    if rng.chance(20) {
        flag |= 1 | 64;
    }
    let unmapped = rng.chance(2);
    let mapq = rng.below(61);
    if unmapped {
        return format!(
            "r{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{qual}",
            std::str::from_utf8(&seq).unwrap()
        );
    }
    format!(
        "r{i}\t{flag}\tc\t{}\t{mapq}\t{cigar}\t*\t0\t0\t{}\t{qual}\tMD:Z:{md}",
        pos + 1,
        std::str::from_utf8(&seq).unwrap()
    )
}

fn run(args: &[&str]) -> Vec<u8> {
    let o = Command::new(BIN).args(args).output().unwrap();
    assert!(
        o.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    o.stdout
}

/// Output without the command line it records: the summary's `command` field and the SAM `@PG`
/// line both contain `-@ N`.
fn without_command(out: Vec<u8>) -> Vec<u8> {
    if out.starts_with(b"{") {
        let mut j: serde_json::Value = serde_json::from_slice(&out).unwrap();
        j.as_object_mut().unwrap().remove("command");
        return serde_json::to_vec(&j).unwrap();
    }
    String::from_utf8(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .flat_map(|l| [l, "\n"])
        .collect::<String>()
        .into_bytes()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// Per-read (name, score, informative sites) from `score`.
fn scores(out: &[u8]) -> Vec<(String, f64, u64)> {
    String::from_utf8(out.to_vec())
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

/// The exclusive upper bound of reads per case: `default`, or PROP_MAX_N + 1 if that is smaller.
fn reads_below(default: usize) -> usize {
    std::env::var("PROP_MAX_N")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .map_or(default, |m| (m + 1).clamp(2, default))
}

fn config() -> ProptestConfig {
    ProptestConfig {
        cases: std::env::var("PROPTEST_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(12),
        ..ProptestConfig::default()
    }
}

proptest! {
    #![proptest_config(config())]

    /// Properties 1 and 5: `score`, `filter` and `damage` give byte-identical output at -@ 1, 2
    /// and 8, and the same records in SAM or BAM give the same scores and summary.
    #[test]
    fn output_does_not_depend_on_threads_or_input_format(seed in any::<u64>(), n in 1usize..reads_below(25_000)) {
        let input = generate(seed, n, "threads");
        let sam = s(&input.sam);
        for cmd in [
            vec!["score", "--quiet", sam],
            vec!["damage", "--quiet", sam],
            vec!["filter", "--quiet", "-t", "1", "--tag", sam, "-o", "-", "--output-format", "sam"],
        ] {
            let one = without_command(run(&[&cmd[..], &["-@", "1"]].concat()));
            for t in ["2", "8"] {
                let other = without_command(run(&[&cmd[..], &["-@", t]].concat()));
                prop_assert!(other == one, "{cmd:?} -@ {t}");
            }
        }
        // SAM to BAM: all records, unchanged (threshold below any score, but above -10000, below
        // which compat keeps no reads, as PMDtools does: docs/blackbox.md 1.7b).
        let bam = input.dir.join("in.bam");
        run(&["filter", "--quiet", "--pmdtools-compat", "-t", "-9999", sam, "-o", s(&bam)]);
        let from_sam = scores(&run(&["score", "--quiet", "--pmdtools-compat", sam]));
        let from_bam = scores(&run(&["score", "--quiet", "--pmdtools-compat", s(&bam)]));
        prop_assert_eq!(from_bam, from_sam);
    }

    /// Properties 2, 3 and 7: filter keeps exactly the reads its threshold rule selects, as an
    /// in-order subsequence of the input with records unchanged apart from the tag, and the
    /// summary counts add up.
    #[test]
    fn filter_keeps_by_score_in_order_and_counts_add_up(
        seed in any::<u64>(),
        n in 1usize..reads_below(5_000),
        t in -5.0f64..10.0,
        width in 0.5f64..10.0,
    ) {
        let input = generate(seed, n, "filter");
        let sam = s(&input.sam);
        let all = scores(&run(&["score", "--quiet", sam]));
        let upper = t + width;
        let (ts, us) = (t.to_string(), upper.to_string());
        let summary = input.dir.join("f.json");
        let kept_sam = run(&[
            "filter", "--quiet", "-t", &ts, "--upper-threshold", &us, sam,
            "-o", "-", "--output-format", "sam", "--summary", s(&summary),
        ]);
        let kept: Vec<String> = String::from_utf8(kept_sam)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(str::to_string)
            .collect();
        // Rule: score >= t and score < upper; a read with no informative sites is not kept by
        // the upper bound (default mode, docs/design.md §5).
        let want: Vec<&String> = all
            .iter()
            .filter(|(_, p, n_sites)| *p >= t && *p < upper && *n_sites > 0)
            .map(|(name, _, _)| name)
            .collect();
        let got: Vec<String> = kept.iter().map(|l| l.split('\t').next().unwrap().to_string()).collect();
        prop_assert_eq!(got.iter().collect::<Vec<_>>(), want);
        // Kept lines are the input lines, in order.
        let mut it = input.records.iter();
        for line in &kept {
            prop_assert!(it.any(|r| r == line), "not an in-order input record: {line}");
        }
        // Counts.
        let j: serde_json::Value = serde_json::from_slice(&fs::read(&summary).unwrap()).unwrap();
        let r = &j["reads"];
        let skipped: u64 = r["skipped"].as_object().unwrap().values().map(|v| v.as_u64().unwrap()).sum();
        prop_assert_eq!(r["input"].as_u64().unwrap(), n as u64);
        prop_assert_eq!(r["scored"].as_u64().unwrap() + skipped, n as u64);
        prop_assert_eq!(
            r["kept"].as_u64().unwrap() + r["removed_by_threshold"].as_u64().unwrap(),
            r["scored"].as_u64().unwrap()
        );
        prop_assert_eq!(r["kept"].as_u64().unwrap(), kept.len() as u64);
        prop_assert_eq!(r["scored"].as_u64().unwrap(), all.len() as u64);
        // Property 6: report accepts every summary pmdshodh writes, and gives the same page twice.
        let page = |name: &str| {
            let out = input.dir.join(name);
            run(&["report", "--quiet", s(&summary), "-o", s(&out)]);
            fs::read(out).unwrap()
        };
        let first = page("a.html");
        prop_assert!(first.starts_with(b"<!DOCTYPE html>"));
        prop_assert!(first == page("b.html"));
    }

    /// Property 4: scores from the reference (-r) equal scores from correct MD tags, in default
    /// and compat mode, ds and ss.
    #[test]
    fn reference_scores_equal_md_scores(seed in any::<u64>(), n in 1usize..reads_below(5_000)) {
        let input = generate(seed, n, "reference");
        let (sam, fasta) = (s(&input.sam), s(&input.fasta));
        for mode in [&[][..], &["--pmdtools-compat"][..], &["--library", "ss"][..]] {
            let md = run(&[&["score", "--quiet"][..], mode, &[sam]].concat());
            let r = run(&[&["score", "--quiet", "-r", fasta][..], mode, &[sam]].concat());
            prop_assert!(md == r, "{mode:?}");
        }
    }
}

/// Write the input of a failing case for inspection:
/// `PROP_SEED=.. PROP_N=.. cargo test --test properties regenerate -- --ignored`.
#[test]
#[ignore]
fn regenerate() {
    let seed = std::env::var("PROP_SEED").unwrap().parse().unwrap();
    let n = std::env::var("PROP_N").unwrap().parse().unwrap();
    let input = generate(seed, n, "regenerate");
    println!("{}", input.sam.display());
}
