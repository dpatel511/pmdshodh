# Benchmarks — 2026-11-14

This round reran pmdshodh and metaDMG on the release candidate. It followed [benchmark-plan.md](benchmark-plan.md) on the same laptop and inputs as the 2026-11-13 round.

- Each configuration ran **5 times**. All times are medians of 5 runs unless a row says otherwise.
- **Builds:**
  - The pmdshodh rows used commit `2fa81df`.
  - The exception is `--pmdtools-compat`, which used `5550981`. That commit fixes a slowdown this round found (see [Compat regression found in this round](#compat-regression-found-in-this-round)).
  - The two commits differ only in compat's reference reconstruction, which the other rows don't run.
- **PMDtools was not rerun.** Its timings come from the 2026-11-13 round: same laptop, inputs and PMDtools version, 5 runs each, 7 h 42 min in total. That round's pmdshodh build was `b14f4ee`, before `-@ N` was changed to mean N threads in total (PR #12).
- **Background load:** other applications were left running during this round, so some timings may include background load. Runs far from their neighbours are pointed out below.

`run.sh` writes the raw logs: commands, `/usr/bin/time -v` output, stderr, and a per-run `times.tsv` with the min–max range. `check.py` writes `summary.tsv` and `checks.tsv`. The logs are kept with the benchmark machine and are not committed. The scripts are [benchmarks/run.sh](benchmarks/run.sh), [benchmarks/check.py](benchmarks/check.py) and [benchmarks/classify.py](benchmarks/classify.py).

## Setup

| Item | Value |
|---|---|
| CPU | Intel Core i5-13500H laptop: 12 cores (4 performance + 8 efficiency), 16 threads |
| OS | WSL 2 (Linux 6.18.58.1-microsoft-standard-WSL2), Ubuntu 24.04, on Windows 11; WSL sees 7.6 GB RAM |
| Power | On mains, Windows "Balanced" plan, sleep off. Other applications were left running (they were closed for the 2026-11-13 round) |
| Storage | Inputs and outputs on the WSL ext4 disk (`~/bench`), not on `/mnt` |
| Cache | Inputs read once with `cat` before each block of runs (warm cache) |
| Clock | Wall time from `/proc/uptime` (monotonic), because the WSL wall clock that `/usr/bin/time` reads can jump by about 2 s. CPU time and peak memory come from `/usr/bin/time -v` |
| pmdshodh | 0.1.0, commits `2fa81df` and `5550981` (compat), `cargo build --release` (Rust 1.99.1, default release profile) |
| PMDtools | 0.60 (Bioconda tag archive, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`), Python 3.12.3; timings from 2026-11-13 |
| samtools | 1.25 (htslib 1.25), Bioconda |
| metaDMG-cpp | Bioconda package `metadmg` 0.5.3; the binary reports itself as 0.5.0 |

Timings ran one after another, never in parallel. The run on `2fa81df` started at 18:57 UTC on 2026-11-14. That run (PMDtools skipped) took 1 h 48 min, ending at 20:44 UTC. The compat rows were rerun on `5550981` from 21:13 to 21:29 UTC.

## Inputs

| Input | Reads | Size | Source |
|---|---|---|---|
| **BER11_L_NOuser** | 14,006,183 | 1.21 GB (md5 `7a73fd9c2895576d8741a6bbca80b10a`) | Ancient horse, untreated library. Derived from ENA [ERR1892229](https://www.ebi.ac.uk/ena/browser/view/ERR1892229) (PRJEB19970, "Ancient genomic changes associated with domestication of the horse"), file `BER11_L.Equus_cab_nucl_wChrUn.realigned.bam` (3,858,832,654 bytes, md5 `6d9b76b98406db5b8c625b67bde93837`). |
| **Motala3** | 23,622,338 | 1.97 GB (md5 `92987d313786940f801c9977fcf515f6`) | Mesolithic human, hg19, ENA [ERR577494](https://www.ebi.ac.uk/ena/browser/view/ERR577494) (PRJEB6272, Lazaridis et al. 2014), file `Motala3.bam` |
| Synthetic 1 M / 10 M / 50 M | 1,000,000 / 10,000,000 / 50,000,000 | — | BER11_L_NOuser records repeated 4 times, then cut to N records. **Synthetic**; used only for the memory test. |

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

Both were scored as double-stranded, with no UDG option.

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
# metaDMG
metaDMG-cpp pmd -@ N in.bam > out.txt                                     # N = 1, 8
```

To reproduce:
- Run `RUNS=5 INPUTS="BER11_L_NOuser Motala3" SYNTH="1M 10M 50M" bash docs/benchmarks/run.sh`.
- Then run `python3 docs/benchmarks/check.py WORK RESULTS BER11_L_NOuser Motala3`.
- For this round, the PMDtools logs and `times.tsv` rows from the 2026-11-13 round were copied into the results folder first, so `run.sh` skipped PMDtools. It skips any run whose log exists.
- The default-vs-PMDtools classification below comes from `python3 docs/benchmarks/classify.py WORK DATA INPUT 1e-9`, run once per input. Without the tolerance argument, it uses 1e-6.

## Results

### Filtering at PMD ≥ 3: PMDtools vs pmdshodh

| | BER11_L_NOuser (14.0 M reads) | Motala3 (23.6 M reads) |
|---|---|---|
| PMDtools pipeline (2026-11-13) | 1,579.1 s (26.3 min) | 2,129.0 s (35.5 min) |
| pmdshodh `filter`, 1 thread | 68.2 s (**23.2×**) | 85.5 s (**24.9×**) |
| pmdshodh `filter`, 8 threads | 10.6 s (149×) | 14.1 s (151×) |
| pmdshodh `filter`, 12 threads | 9.0 s (175×) | 11.3 s (188×) |
| pmdshodh `filter --pmdtools-compat`, 1 thread | 62.8 s (25.1×) | 90.5 s (23.5×) |
| pmdshodh `filter --pmdtools-compat`, 8 threads | 11.8 s (134×) | 16.4 s (130×) |

- Speed-ups are relative to the PMDtools pipeline on the same input.
- The 1-thread rows are the like-for-like comparison of the scoring itself. The PMDtools pipeline uses about 1 CPU-second per wall-second (1,566 s CPU for BER11_L_NOuser).
- PMDtools' 5 runs ranged from 1,490 to 1,630 s on BER11_L_NOuser and from 2,077 to 2,207 s on Motala3. Min–max ranges for every configuration are in `summary.tsv`.
- **The 1-thread rows are slower than in the 2026-11-13 round** (51.0 s and 76.7 s then). `-@ 1` now does reading, BGZF and scoring one after another, where it used to use about 1.2 cores. In alternating runs on BER11_L_NOuser (3 each, same session), `b14f4ee` took 55.3–60.9 s and `5550981` took 62.8–64.9 s, and the new build used no more CPU time (61.4–63.4 s, against 64.5–74.9 s). So `-@ 1` is about 15% slower than before. Most of the rest of the gap from 51.0 s is background load: `b14f4ee` itself took 55–61 s that evening.
- The compat rows ran in a later block (on `5550981`) than the default rows, so the two are not directly comparable. Default and compat were close on 1 thread: compat was faster on BER11_L_NOuser and slower on Motala3.

### Thread scaling (`pmdshodh filter`)

| `-@` | BER11_L_NOuser | Speed-up | CPU cores used | Motala3 | Speed-up | CPU cores used |
|---|---|---|---|---|---|---|
| 1 | 68.2 s | 1.00× | 1.0 | 85.5 s | 1.00× | 1.0 |
| 2 | 31.9 s | 2.14× | 1.8 | 39.4 s | 2.17× | 1.9 |
| 4 | 15.2 s | 4.49× | 3.7 | 20.9 s | 4.09× | 3.5 |
| 8 | 10.6 s | 6.43× | 6.9 | 14.1 s | 6.06× | 6.7 |
| 12 | 9.0 s | 7.58× | 9.6 | 11.3 s | 7.57× | 9.6 |

- "CPU cores used" is median CPU time (user + system) divided by median wall time.
- **`-@ N` now uses at most N cores.** All work runs on one pool of N threads (CHANGELOG.md, docs/design.md §7). This is what pipelines expect when they pass their CPU allocation.
- `-@ 2` is more than 2× faster than `-@ 1` because `-@ 1` used more CPU time in total: 68.1 s on BER11_L_NOuser against 56.7 s for `-@ 2`. That difference was not investigated further.
- **Compared with the 2026-11-13 round** (`b14f4ee`, where `-@ N` used more than N cores):
  - At `-@ 2` and `-@ 4`, runs are slower because they use fewer cores: 20.5 s on 3.1 cores became 31.9 s on 1.8 (BER11_L_NOuser).
  - At `-@ 8` they are about as fast on fewer cores: 10.2 s on 9.2 cores became 10.6 s on 6.9, and 13.4 s on 8.9 became 14.1 s on 6.7.
  - At `-@ 12` they are faster: 10.2 s became 9.0 s, and 12.8 s became 11.3 s.

### Per-read scores: pmdshodh vs metaDMG

| | BER11_L_NOuser | Motala3 |
|---|---|---|
| metaDMG `pmd`, `-@ 1` | **56.6 s** | **74.4 s** |
| metaDMG `pmd`, `-@ 8` | 56.6 s | 75.7 s |
| pmdshodh `score`, 1 thread | 61.2 s | 85.9 s |
| pmdshodh `score`, 8 threads | 11.4 s | 16.8 s |

- **On one thread, metaDMG is faster**: by 8% on BER11_L_NOuser and 15% on Motala3 (8% on both in the 2026-11-13 round).
  - **Changed after this round:** `score` and `filter` no longer count damage tables without `--summary`. On the first 2 M reads of BER11_L_NOuser, `score -@ 1` then took 6.18 s against metaDMG's 8.23 s (CHANGELOG.md). The full inputs have not been rerun.
  - The two tools alternated run by run. metaDMG was faster in 4 of the 5 pairs on BER11_L_NOuser and in all 5 on Motala3.
  - On BER11_L_NOuser, the pair where pmdshodh was faster was metaDMG's slowest run (68.3 s).
- The profile below suggests why: about a quarter of pmdshodh's single-thread time goes to damage tables that `score` and `filter` don't use.
- metaDMG's `-@` option made almost no difference to its time. Its output was identical for `-@ 1` and `-@ 8`.
- pmdshodh with 8 threads is 5.0× faster than metaDMG's faster setting on BER11_L_NOuser and 4.4× faster on Motala3.
- The two tools don't compute the same score, so this compares speed only.
- `score` takes longer than `filter` because it writes a line for every read (14–24 M lines), whereas `filter` writes only the kept reads.

### Memory

Peak resident memory (`/usr/bin/time -v`, median of 5):

| | 1 thread | 8 threads |
|---|---|---|
| pmdshodh `filter`, BER11_L_NOuser / Motala3 | 11 MB / 10 MB | 38 MB / 35 MB |
| pmdshodh `filter`, synthetic 1 M / 10 M / 50 M reads | 11 / 11 / 11 MB | 37 / 38 / 40 MB |
| pmdshodh `score`, BER11_L_NOuser / Motala3 | 10 MB / 10 MB | 27 MB / 24 MB |
| pmdshodh `filter --pmdtools-compat`, BER11_L_NOuser / Motala3 | 11 MB / 11 MB | 41 MB / 34 MB |
| metaDMG `pmd` | 10 MB | 10–11 MB |
| PMDtools pipeline, largest single process (2026-11-13) | 15 MB | — |

- **pmdshodh's memory is bounded.** With several threads it holds up to four batches of 10,000 records at once. That is why `-@ 8` uses more memory than `-@ 1`, and more than in the 2026-11-13 round (27–28 MB).
- On the synthetic inputs, 50× more reads added 3 MB at `-@ 8` and nothing at `-@ 1`.
- Time is linear in the number of reads. On 1 thread, the synthetic inputs took 4.3, 4.2 and 4.2 s per million reads; on 8 threads, 0.8, 0.71 and 0.73 s.
- For the PMDtools pipeline, `/usr/bin/time` reports the largest single process, not the pipeline's total.

## Compat regression found in this round

On `2fa81df`, `filter --pmdtools-compat -@ 8` was slow and erratic.
- **Times:** 14.5–44.3 s on BER11_L_NOuser (median 30.4 s; 13.0 s in the 2026-11-13 round). On Motala3, two of its five runs took 23.6 s and 26.2 s, against 18.4–18.9 s for the other three.
- **Kernel time:** up to 53 s per run, against 4.9 s in the fastest run.
- **Reproducible:** in 4 more back-to-back runs (14.7–42.4 s), voluntary context switches rose from 0.18 M to 3.2 M.

**Cause:**
- `perf record -e context-switches` showed that nearly all of the switches were waits on glibc's malloc arena lock (`__lll_lock_wait_private` under `realloc` and `free`).
- They came from the compat reconstruction, which grows several vectors per read, and from record reading on the serial thread.
- In glibc, `realloc` always takes the arena lock. Under the shared thread pool (PR #12), those locks contended.
- Raising glibc's arena limits or trim threshold did not help.
- Sizing the vectors up front (`5550981`) removed it. Alternating 5 runs of each build gave:

| `filter --pmdtools-compat -@ 8`, BER11_L_NOuser | Wall | System CPU | Voluntary context switches |
|---|---|---|---|
| `2fa81df` | 23.9–66.6 s | 14.9–129.9 s | 0.84–6.12 M |
| `5550981` | 10.8–11.8 s | 1.6–1.9 s | 45–50 k |

**Scope:**
- Scores are unchanged; the parity checks below ran on `5550981`.
- Default `filter` and `score` don't run this code. Default `filter -@ 8` spent 1.8–2.7 s in the kernel in every run (9.4–13.9 s in the 2026-11-13 round).

## Single-thread profile

A `perf record -F 1999 -g` profile of `score -@ 1` on the first 2,000,000 reads of BER11_L_NOuser shows where single-thread time goes. It was a release build with frame pointers and line tables, on the `2fa81df` source; perf 6.8.12.

| Part of `score -@ 1` | Share of samples |
|---|---|
| Damage counting (`DamageCounts::add_read`) | 26.7% |
| Scoring (`score_read`, including `site_term_ds` at 7.6%) | 21.3% |
| BGZF decompression (zlib-rs inflate) | 15.8% |
| Reference reconstruction from MD (`reference_along_read`) | 5.9% |
| BAM sequence decoding | about 6% |

- `filter -@ 1` looks the same: damage counting is 25.6% of its samples. BGZF compression adds about 8% (its two largest zlib-rs functions).
- **Damage counting is the largest single cost.** Its tables are used only by `--summary` and the `damage` command, so `score` and `filter` without `--summary` spend about a quarter of their time on tables nobody reads.
- metaDMG's `pmd` builds no such tables. This is the likely reason it is faster on one thread, but that has not been measured.

Since this round, `score` and `filter` skip damage counting when no `--summary` is asked for. Output is byte-identical (CHANGELOG.md).

## Correctness checks

Run on every benchmark input with `check.py`. All checks passed (compat checks on `5550981`, the rest on `2fa81df`).

| Check | BER11_L_NOuser | Motala3 |
|---|---|---|
| Per-read scores: PMDtools `--printDS` vs pmdshodh `score --pmdtools-compat` (tolerance 1e-6) | 14,006,183 reads, largest difference 6.0e-14, 0 above tolerance | 23,622,338 reads, largest difference 5.0e-14, 0 above tolerance |
| Kept reads: PMDtools vs `filter --pmdtools-compat` (whole records, in order) | 1,031,193, identical including tags | 301,867, identical including tags |
| `filter` output identical for 1, 2, 4, 8 and 12 threads | 1,029,380 records each | 301,722 records each |
| `filter --pmdtools-compat` and `score` identical for 1 and 8 threads | pass | pass |

**Default mode vs PMDtools** (kept reads at PMD ≥ 3):
- BER11_L_NOuser: 1,907 reads kept only by PMDtools, 94 only by pmdshodh, 1,029,286 by both.
- Motala3: 242, 97 and 301,625.
- These are the same counts as in the 2026-11-13 round. Default scoring has not changed since then.

Every read whose default score differs from PMDtools' comes from a documented departure (docs/design.md §5); none are unexplained. This was checked on the 2026-11-13 round's outputs (2026-11-14).

- **Method:** each read's default score (`score -@ 1`) was compared with PMDtools' `--printDS` score for the same read, in input order. Read names were checked against the BAM, and the three streams had the same length.
- **Tolerance 1e-9:** the counts are the ones in the table.
- **Tolerance 1e-6** (the parity check's): 30 fewer CIGAR-fix reads on BER11_L_NOuser (3,953) and 54 fewer on Motala3 (3,703), which differ by less than that. Nothing else changes, and the reads that cross PMD 3 are the same.
- **Threshold crossings:** the reads that cross PMD 3 add up to the kept-read differences above: 2,001 = 1,907 + 94, and 339 = 242 + 97.
- **Preliminary round:** the classification matches it (2026-11-12, commit `a98ae07`).

| Cause of the difference (tolerance 1e-9) | BER11_L_NOuser | of which cross PMD 3 | Motala3 | of which cross PMD 3 |
|---|---|---|---|---|
| CIGAR reconstruction fix (`I`/`S` with `D`/`=`) | 3,983 | 1,960 | 3,757 | 339 |
| Lowercase MD bases | 336 | 38 | 0 | 0 |
| Lowercase MD + reverse-read shift | 3 | 3 | 0 | 0 |
| CIGAR fix + lowercase MD | 2 | 0 | 0 | 0 |
| **Unexplained** | **0** | | **0** | |

## Earlier rounds

- **2026-11-13** (5 runs, commit `b14f4ee`, other applications closed): single-thread `filter` took 51.0 s (BER11_L_NOuser) and 76.7 s (Motala3), and 8-thread `filter` 10.2 s and 13.4 s. Its `-@ N` used more than N cores. The thread-scaling comparison above uses this round.
- **2026-11-12** (preliminary, 3 runs, commit `a98ae07`):
  - Single-thread `filter` took 84.0 s and 112.7 s, and 8-thread `filter` 18.6 s and 24.7 s.
  - It timed with the WSL wall clock and wasn't interleaved with later rounds, so this is not a controlled comparison.
  - Its PMDtools compat-mode check initially failed on 339 reads; PR #6 fixed the cause before the later rounds.

The controlled before/after measurements of each speed change are in the CHANGELOG.

## Departures from the plan

- **PMDtools was not rerun for this round.** Its timings are reused from 2026-11-13 (same laptop, inputs and version). This changes only the speed-up ratios, which compare the two rounds.
- **No PMDtools run with `samtools calmd`.** All inputs already have MD tags, and this run would have needed two reference genomes.
- **`--printDS` and the compat `score` were run once each.** They exist for the parity check.
- **`--ss` is not benchmarked.** No single-stranded library has been tested.
