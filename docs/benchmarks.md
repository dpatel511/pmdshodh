# Benchmarks — 2026-11-17

This round benchmarks the 0.2.1 release candidate. It followed [benchmark-plan.md](benchmark-plan.md) on the same laptop and inputs as the earlier rounds. New in this round: `--library ss`, reference-based scoring (`-r`) and CRAM input.

- Each configuration ran **5 times**. All times are medians of 5 runs unless a row says otherwise.
- **Build:** commit `144080e`, the source of 0.2.1. The binary still reported version 0.2.0, because the version number was changed after the round; no source changed.
- **PMDtools was not rerun.** Its timings come from the 2026-11-13 round: same laptop, inputs and PMDtools version, 5 runs each. Its outputs from that round were reused for the parity checks.
- **Background load:** other applications were left running. In the BER11_L_NOuser block the machine sped up and slowed down by about 20% between runs:
  - Run 5 of the alternating `compat` / `score` / metaDMG block was faster for all three tools: 49.2, 40.3 and 51.1 s, against 57.5–60.2, 47.5–50.6 and 53.6–64.5 s for the other runs.
  - `score --library ss -@ 8` ranged from 10.8 to 22.6 s, although its CPU time stayed at 46–52 s.
  - So on BER11_L_NOuser, differences of up to about 20% between rows that ran in different blocks are not meaningful. The Motala3 rows were steadier.

`run.sh` writes the raw logs: commands, `/usr/bin/time -v` output, stderr, and a per-run `times.tsv`. `check.py` writes `summary.tsv` (median and min–max range) and `checks.tsv`. The logs are kept with the benchmark machine and are not committed. The scripts are [benchmarks/run.sh](benchmarks/run.sh), [benchmarks/check.py](benchmarks/check.py) and [benchmarks/classify.py](benchmarks/classify.py).

## Setup

| Item | Value |
|---|---|
| CPU | Intel Core i5-13500H laptop: 12 cores (4 performance + 8 efficiency), 16 threads |
| OS | WSL 2 (Linux 6.18.58.1-microsoft-standard-WSL2), Ubuntu 24.04, on Windows 11; WSL sees 7.6 GB RAM |
| Power | On mains, Windows "Balanced" plan. Other applications were left running |
| Storage | Inputs and outputs on the WSL ext4 disk (`~/bench`), not on `/mnt` |
| Cache | Inputs read once with `cat` before each block of runs (warm cache) |
| Clock | Wall time from `/proc/uptime` (monotonic), because the WSL wall clock that `/usr/bin/time` reads can jump by about 2 s. CPU time and peak memory come from `/usr/bin/time -v` |
| pmdshodh | commit `144080e` (0.2.1 source), `cargo build --release` (Rust 1.99.1, default release profile) |
| PMDtools | 0.60 (Bioconda tag archive, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`), Python 3.12.3; timings and outputs from 2026-11-13 |
| samtools | 1.25 (htslib 1.25), Bioconda |
| metaDMG-cpp | Bioconda package `metadmg` 0.5.3; the binary reports itself as 0.5.0 |
| References | equCab2 (soft-masked) for BER11_L_NOuser and human_g1k_v37 for Motala3, indexed with `samtools faidx`; the same references as the 0.2.0 checks (CHANGELOG.md) |

Timings ran one after another, never in parallel, from 12:15 to 14:30 UTC on 2026-11-17 (PMDtools skipped). One run was kept from an earlier, interrupted start of the same script on the same build: `filter -@ 1` run 1 on BER11_L_NOuser (45.6 s, the fastest of its five).

## Inputs

| Input | Reads | Size | Source |
|---|---|---|---|
| **BER11_L_NOuser** | 14,006,183 | 1.21 GB (md5 `7a73fd9c2895576d8741a6bbca80b10a`) | Ancient horse, untreated library. Derived from ENA [ERR1892229](https://www.ebi.ac.uk/ena/browser/view/ERR1892229) (PRJEB19970, "Ancient genomic changes associated with domestication of the horse"), file `BER11_L.Equus_cab_nucl_wChrUn.realigned.bam` (3,858,832,654 bytes, md5 `6d9b76b98406db5b8c625b67bde93837`). |
| **Motala3** | 23,622,338 | 1.97 GB (md5 `92987d313786940f801c9977fcf515f6`) | Mesolithic human, hg19, ENA [ERR577494](https://www.ebi.ac.uk/ena/browser/view/ERR577494) (PRJEB6272, Lazaridis et al. 2014), file `Motala3.bam` |
| Synthetic 1 M / 10 M / 50 M | 1,000,000 / 10,000,000 / 50,000,000 | — | BER11_L_NOuser records repeated 4 times, then cut to N records. **Synthetic**; used only for the memory test. |
| CRAM copies | same reads | — | The two BAMs written as CRAM by samtools 1.25 with the references above: the files of the 0.2.0 CRAM check |

The ENA files were downloaded on 2026-11-11 and checked against ENA's md5. ENA data are shared openly under INSDC policy. No BAMs are committed.

**Why a subset of ERR1892229:** the file mixes USER-treated (`_user_`) and untreated (`NOuser`) libraries.
- In the first 2 M reads, every `_user_` read (1,514,462) had an MD tag that does not match its alignment: MD covers 5 more reference bases than the CIGAR (4 for one library).
- pmdshodh stops on such reads with an error, as CLAUDE.md requires.
- All 485,538 `NOuser` reads had valid MD tags.

The subset keeps the three `NOuser` read groups:

```bash
samtools view -@ 8 -b -R ber_nouser_rg.txt -o BER11_L_NOuser.bam BER11_L.bam
# ber_nouser_rg.txt: BER11_L_E12.Dab2nd_NOuser_AGCATG, BER11_L_E12.Dab2nd_NOuser_AGCATG.1, BER11_L_E12.Dab2nd_NOuser_AGCATG.2
```

**Damage profiles** (`pmdshodh damage`, 5′ C→T rate):
- BER11_L_NOuser: 18.4% at the first position, 3.0% at position 30.
- Motala3: 1.44% at the first position, 0.17% at position 30. This is consistent with a UDG-treated library; I haven't confirmed it from the paper.

Both are double-stranded libraries and were scored as such, with no UDG option, except in the `--library ss` rows.

## Commands

From the logged `.cmd` files (paths shortened):

```bash
# PMDtools, as nf-core pmdtools/filter runs it (no calmd: the inputs already have MD tags)
samtools view -h in.bam | python3 -W ignore pmdtools.0.60.py --threshold 3 --header | samtools view -b -o out.bam -
# PMDtools per-read scores (for the parity check)
samtools view -h in.bam | python3 -W ignore pmdtools.0.60.py --printDS > out.ds
# pmdshodh
pmdshodh filter --threshold 3 --quiet -@ N in.bam -o out.bam             # N = 1, 2, 4, 8, 12
pmdshodh filter --threshold 3 --quiet --pmdtools-compat -@ N in.bam -o out.bam   # N = 1, 8
pmdshodh score --quiet -@ N in.bam -o out.tsv                             # N = 1, 8
pmdshodh score --quiet --library ss -@ N in.bam -o out.tsv                # N = 1, 8
pmdshodh filter --threshold 3 --quiet -r ref.fa -@ N in.bam -o out.bam    # N = 1, 8
pmdshodh filter --threshold 3 --quiet -r ref.fa -@ N in.cram -o out.bam   # N = 1, 8
# metaDMG
metaDMG-cpp pmd -@ N in.bam > out.txt                                     # N = 1, 8
```

To reproduce:
- Run `RUNS=5 INPUTS="BER11_L_NOuser Motala3" SYNTH="1M 10M 50M" bash docs/benchmarks/run.sh`. The `-r` and CRAM rows run for inputs that have a reference in `REF_MAP` and a CRAM in `CRAMS`.
- Then run `python3 docs/benchmarks/check.py WORK RESULTS BER11_L_NOuser Motala3`.
- For this round, the PMDtools logs, outputs and `times.tsv` rows from the 2026-11-13 round were copied into the results folder first, so `run.sh` skipped PMDtools. It skips any run whose log exists.

## Results

### Filtering at PMD ≥ 3: PMDtools vs pmdshodh

| | BER11_L_NOuser (14.0 M reads) | Motala3 (23.6 M reads) |
|---|---|---|
| PMDtools pipeline (2026-11-13) | 1,579.1 s (26.3 min) | 2,129.0 s (35.5 min) |
| pmdshodh `filter`, 1 thread | 49.5 s (**31.9×**) | 54.8 s (**38.9×**) |
| pmdshodh `filter`, 8 threads | 7.6 s (208×) | 9.8 s (217×) |
| pmdshodh `filter`, 12 threads | 6.6 s (239×) | 9.6 s (222×) |
| pmdshodh `filter --pmdtools-compat`, 1 thread | 57.9 s (27.3×) | 62.8 s (33.9×) |
| pmdshodh `filter --pmdtools-compat`, 8 threads | 8.9 s (177×) | 12.2 s (175×) |

- Speed-ups are relative to the PMDtools pipeline on the same input, measured on 2026-11-13. The two tools were not timed in the same session.
- The 1-thread rows are the like-for-like comparison of the scoring itself. The PMDtools pipeline uses about 1 CPU-second per wall-second (1,566 s CPU for BER11_L_NOuser).
- PMDtools' 5 runs ranged from 1,490 to 1,630 s on BER11_L_NOuser and from 2,077 to 2,207 s on Motala3. pmdshodh's ranges are in `summary.tsv`; for `filter -@ 1` they were 45.6–50.5 s and 53.9–55.7 s.
- The compat rows ran in a later block than the default rows, so the two are not directly comparable on BER11_L_NOuser (see the background-load note above).

### Thread scaling (`pmdshodh filter`)

| `-@` | BER11_L_NOuser | Speed-up | CPU cores used | Motala3 | Speed-up | CPU cores used |
|---|---|---|---|---|---|---|
| 1 | 49.5 s | 1.00× | 1.0 | 54.8 s | 1.00× | 1.0 |
| 2 | 23.6 s | 2.10× | 1.7 | 27.9 s | 1.96× | 1.8 |
| 4 | 10.7 s | 4.63× | 3.7 | 14.3 s | 3.83× | 3.6 |
| 8 | 7.6 s | 6.51× | 7.0 | 9.8 s | 5.59× | 7.1 |
| 12 | 6.6 s | 7.50× | 9.3 | 9.6 s | 5.71× | 9.5 |

- "CPU cores used" is median CPU time (user + system) divided by median wall time. `-@ N` uses at most N cores: all work runs on one pool of N threads (docs/design.md §7).
- **On BER11_L_NOuser, 1 thread used more CPU time in total than 2 or 4:** 47.9 s, against 40.6 s and 39.8 s. That is why the speed-ups at 2 and 4 threads are above 2× and 4×. Motala3 does not show this (52.9, 51.2 and 51.0 s). The 2026-11-14 round saw the same on BER11_L_NOuser. The cause has not been found.
- On Motala3, 12 threads were barely faster than 8 (9.6 s against 9.8 s), although they used more CPU time (90.9 s against 69.6 s).

### Per-read scores: pmdshodh vs metaDMG

| | BER11_L_NOuser | Motala3 |
|---|---|---|
| metaDMG `pmd`, `-@ 1` | 61.7 s | 71.9 s |
| metaDMG `pmd`, `-@ 8` | 51.9 s | 72.4 s |
| pmdshodh `score`, 1 thread | **48.8 s** | **60.6 s** |
| pmdshodh `score`, 8 threads | 9.4 s | 13.8 s |

- **On one thread, pmdshodh is now faster.** The two tools (and `filter --pmdtools-compat`) alternated run by run. `score` was faster in all 5 pairs on both inputs: metaDMG took 13–32% longer per pair on BER11_L_NOuser and 6–20% longer on Motala3.
  - In the 2026-11-14 round (the 0.1.0 candidate), metaDMG was 8–15% faster. Right after that round, before 0.1.0 was tagged, `score` and `filter` stopped counting damage tables that are only needed for `--summary`; that round's profile had put them at about a quarter of the single-thread time.
- metaDMG's `-@ 8` row on BER11_L_NOuser is faster than its `-@ 1` row because of the background-load swings in that block, not threads. On Motala3 the two settings took the same time, and its output was identical for `-@ 1` and `-@ 8` on both inputs.
- pmdshodh with 8 threads is 5.5× faster than metaDMG's faster setting on BER11_L_NOuser and 5.2× faster on Motala3.
- The two tools don't compute the same score, so this compares speed only.
- `score` takes longer than `filter` with several threads because it writes a line for every read (14–24 M lines), whereas `filter` writes only the kept reads.

### Single-stranded scoring, `-r` and CRAM

| | BER11_L_NOuser, 1 thread | 8 threads | Motala3, 1 thread | 8 threads |
|---|---|---|---|---|
| `score` (double-stranded) | 48.8 s | 9.4 s | 60.6 s | 13.8 s |
| `score --library ss` | 40.6 s (39.5–54.4) | 13.1 s (10.8–22.6) | 59.3 s | 14.4 s |
| `filter` (MD tags) | 49.5 s | 7.6 s | 54.8 s | 9.8 s |
| `filter -r`, BAM | 51.3 s | 12.6 s | 65.4 s | 14.8 s |
| `filter -r`, CRAM | 139.6 s | 38.7 s | 201.5 s | 58.5 s |

- **`--library ss` was timed on these double-stranded libraries.** This round has no single-stranded library; these rows measure speed and thread identity only (a real single-stranded library was checked separately: CHANGELOG.md, Unreleased). BER11_L_NOuser's ss rows ran during the background-load swings (ranges in brackets), so the ss and ds rows can't be compared there. On Motala3 they took about the same time.
- **`-r`** looks up every aligned base in the FASTA instead of reading the MD tag. It cost the most with several threads: 7.6 → 12.6 s and 9.8 → 14.8 s at `-@ 8`. The rows ran in different blocks, so the size of the cost is approximate.
- **CRAM** took 2.7–4.0× as long as the same reads as BAM with `-r`. On top of the `-r` work, it decodes the CRAM and loads whole reference sequences, checking each against the header's `M5`. That check now runs on the thread pool (CHANGELOG.md, 0.2.1); this round has no run without it.

### Memory

Peak resident memory (`/usr/bin/time -v`, median of 5; MB = 1024 KiB):

| | 1 thread | 8 threads |
|---|---|---|
| pmdshodh `filter`, BER11_L_NOuser / Motala3 | 11 MB / 10 MB | 38 MB / 33 MB |
| pmdshodh `filter`, synthetic 1 M / 10 M / 50 M reads | 11 / 11 / 11 MB | 36 / 37 / 38 MB |
| pmdshodh `score`, BER11_L_NOuser / Motala3 | 10 MB / 10 MB | 26 MB / 25 MB |
| pmdshodh `score --library ss`, BER11_L_NOuser / Motala3 | 10 MB / 10 MB | 33 MB / 24 MB |
| pmdshodh `filter --pmdtools-compat`, BER11_L_NOuser / Motala3 | 11 MB / 11 MB | 39 MB / 36 MB |
| pmdshodh `filter -r`, BAM, BER11_L_NOuser / Motala3 | 87 MB / 84 MB | 112 MB / 112 MB |
| pmdshodh `filter -r`, CRAM, BER11_L_NOuser / Motala3 | 506 MB / 726 MB | 1,190 MB / 1,248 MB |
| metaDMG `pmd` | 10 MB | 10 MB |
| PMDtools pipeline, largest single process (2026-11-13) | 15 MB | — |

- **pmdshodh's memory is bounded.** With several threads it holds up to four batches of 10,000 records at once, which is why `-@ 8` uses more memory than `-@ 1`.
- On the synthetic inputs, 50× more reads added 2 MB at `-@ 8` and nothing at `-@ 1`.
- Time is linear in the number of reads. On 1 thread, the synthetic inputs took 3.0, 2.9 and 3.0 s per million reads; on 8 threads, 0.60, 0.54 and 0.55 s.
- `-r` added 74–79 MB on BAM input.
- **CRAM** decoding loads whole reference sequences (docs/design.md §11). Their cache is cleared past 256 Mi bases (CHANGELOG.md).
- For the PMDtools pipeline, `/usr/bin/time` reports the largest single process, not the pipeline's total.

## Correctness checks

Run on every benchmark input with `check.py`. All 24 checks passed (12 per input; `check.py` also prints 8 informational lines, which are not checks).

| Check | BER11_L_NOuser | Motala3 |
|---|---|---|
| Per-read scores: PMDtools `--printDS` vs pmdshodh `score --pmdtools-compat` (tolerance 1e-6) | 14,006,183 reads, largest difference 6.0e-14, 0 above tolerance | 23,622,338 reads, largest difference 5.0e-14, 0 above tolerance |
| Kept reads: PMDtools vs `filter --pmdtools-compat` (whole records, in order) | 1,031,193, identical including tags | 301,867, identical including tags |
| `filter` output identical for 1, 2, 4, 8 and 12 threads | 1,029,380 records each | 301,722 records each |
| `filter --pmdtools-compat`, `score` and `score --library ss` identical for 1 and 8 threads | pass | pass |
| `filter -r` (BAM and CRAM) identical for 1 and 8 threads | 1,029,380 records each | 301,706 records each |
| `filter -r` on CRAM vs on BAM, fields 1–11 (CRAM leaves out the MD and NM tags htslib regenerates) | identical | identical |

- **`filter -r` against `filter` with MD tags:** identical on BER11_L_NOuser (1,029,380 records). On Motala3, `-r` keeps 301,706 reads against 301,722, the same counts as the 0.2.0 check. Motala3 has 107 reads whose MD tag disagrees with the reference, and `-r` scores those from the reference (CHANGELOG.md, 0.2.0).

**Default mode vs PMDtools** (kept reads at PMD ≥ 3):
- BER11_L_NOuser: 1,907 reads kept only by PMDtools, 94 only by pmdshodh, 1,029,286 by both.
- Motala3: 242, 97 and 301,625.
- These are the same counts as in the 2026-11-13 and 2026-11-14 rounds. Default scores on these inputs have not changed since then.

Every read whose default score differs from PMDtools' comes from a documented departure (docs/design.md §5); none are unexplained. `classify.py` was not rerun for this round. It was run on the 2026-11-13 round's outputs (2026-11-14), and the kept-read counts above are unchanged.

- **Method:** each read's default score (`score -@ 1`) was compared with PMDtools' `--printDS` score for the same read, in input order. Read names were checked against the BAM, and the three streams had the same length.
- **Tolerance 1e-9:** the counts are the ones in the table.
- **Tolerance 1e-6** (the parity check's): 30 fewer CIGAR-fix reads on BER11_L_NOuser (3,953) and 54 fewer on Motala3 (3,703), which differ by less than that. Nothing else changes, and the reads that cross PMD 3 are the same.
- **Threshold crossings:** the reads that cross PMD 3 add up to the kept-read differences above: 2,001 = 1,907 + 94, and 339 = 242 + 97.

| Cause of the difference (tolerance 1e-9) | BER11_L_NOuser | of which cross PMD 3 | Motala3 | of which cross PMD 3 |
|---|---|---|---|---|
| CIGAR reconstruction fix (`I`/`S` with `D`/`=`) | 3,983 | 1,960 | 3,757 | 339 |
| Lowercase MD bases | 336 | 38 | 0 | 0 |
| Lowercase MD + reverse-read shift | 3 | 3 | 0 | 0 |
| CIGAR fix + lowercase MD | 2 | 0 | 0 | 0 |
| **Unexplained** | **0** | | **0** | |

## Earlier rounds

These rounds were not interleaved with this one, so the differences between rounds are not controlled measurements. The controlled before/after measurements of each speed change are in the CHANGELOG.

- **2026-11-14** (5 runs, 0.1.0 candidate, commits `2fa81df` and `5550981`):
  - `filter` took 68.2 s and 85.5 s on 1 thread (BER11_L_NOuser, Motala3), 10.6 s and 14.1 s on 8 threads, and 9.0 s and 11.3 s on 12.
  - metaDMG `pmd` was 8–15% faster than `score` on one thread (56.6 s against 61.2 s, and 74.4 s against 85.9 s).
  - A `perf` profile of `score -@ 1` put 26.7% of the time in damage counting, which `score` and `filter` only need for `--summary`. The released 0.1.0 skips it without `--summary`.
  - It found a slowdown in `filter --pmdtools-compat -@ 8` (14.5–44.3 s on BER11_L_NOuser): threads waited on glibc's malloc arena lock while the compat reconstruction grew its vectors. Sizing the vectors up front (`5550981`) fixed it (10.8–11.8 s).
- **2026-11-13** (5 runs, commit `b14f4ee`, other applications closed): single-thread `filter` took 51.0 s and 76.7 s, and 8-thread `filter` 10.2 s and 13.4 s. Its `-@ N` used more than N cores. The PMDtools timings in this document come from this round.
- **2026-11-12** (preliminary, 3 runs, commit `a98ae07`):
  - Single-thread `filter` took 84.0 s and 112.7 s, and 8-thread `filter` 18.6 s and 24.7 s.
  - It timed with the WSL wall clock and wasn't interleaved with later rounds, so this is not a controlled comparison.
  - Its PMDtools compat-mode check initially failed on 339 reads; PR #6 fixed the cause before the later rounds.

## Departures from the plan

- **PMDtools was not rerun for this round.** Its timings and outputs are reused from 2026-11-13 (same laptop, inputs and version). This affects only the speed-up ratios, which compare two sessions.
- **No PMDtools run with `samtools calmd`.** All inputs already have MD tags. The `-r` rows are compared with pmdshodh's own MD path; the `calmd | pmdtools` comparison of `-r` is in CHANGELOG.md (0.2.0).
- **`--printDS` and the compat `score` were run once each.** They exist for the parity check.
- **No single-stranded library in this round.** `--library ss` was timed on double-stranded data, and its scores were not compared with PMDtools here. A real single-stranded library was checked separately (CHANGELOG.md, Unreleased).
