# Benchmarks — 2026-11-13

Each configuration ran **5 times** on one laptop, following [benchmark-plan.md](benchmark-plan.md). The pmdshodh build was commit `b14f4ee`, which includes the speed work of PRs #9 and #10. All times are medians of 5 runs unless a row says otherwise.

`run.sh` writes the raw logs: commands, `/usr/bin/time -v` output, stderr, and a per-run `times.tsv` with the min–max range. `check.py` writes `summary.tsv` and `checks.tsv`. The logs are kept with the benchmark machine and are not committed. The scripts are [benchmarks/run.sh](benchmarks/run.sh) and [benchmarks/check.py](benchmarks/check.py).

## Setup

| Item | Value |
|---|---|
| CPU | Intel Core i5-13500H laptop: 12 cores (4 performance + 8 efficiency), 16 threads |
| OS | WSL 2 (Linux 6.18.58.1-microsoft-standard-WSL2), Ubuntu 24.04, on Windows 11; WSL sees 7.6 GB RAM |
| Power | On mains, Windows "Balanced" plan, sleep off, other applications closed |
| Storage | Inputs and outputs on the WSL ext4 disk (`~/bench`), not on `/mnt` |
| Cache | Inputs read once with `cat` before each block of runs (warm cache) |
| Clock | Wall time from `/proc/uptime` (monotonic), because the WSL wall clock that `/usr/bin/time` reads can jump by about 2 s. CPU time and peak memory come from `/usr/bin/time -v` |
| pmdshodh | 0.1.0, commit `b14f4ee`, `cargo build --release` (Rust 1.99.1, default release profile) |
| PMDtools | 0.60 (Bioconda tag archive, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`), Python 3.12.3 |
| samtools | 1.25 (htslib 1.25), Bioconda |
| metaDMG-cpp | Bioconda package `metadmg` 0.5.3; the binary reports itself as 0.5.0 |

Timings ran one after another, never in parallel. The whole run took 7 h 42 min (2026-11-13 18:52 to 2026-11-14 02:34 UTC).

## Inputs

| Input | Reads | Size | Source |
|---|---|---|---|
| **BER11_L_NOuser** | 14,006,183 | 1.21 GB (md5 `7a73fd9c2895576d8741a6bbca80b10a`) | Ancient horse, untreated library. Derived from ENA [ERR1892229](https://www.ebi.ac.uk/ena/browser/view/ERR1892229) (PRJEB19970, "Ancient genomic changes associated with domestication of the horse"), file `BER11_L.Equus_cab_nucl_wChrUn.realigned.bam` (3,858,832,654 bytes, md5 `6d9b76b98406db5b8c625b67bde93837`). |
| **Motala3** | 23,622,338 | 1.97 GB (md5 `92987d313786940f801c9977fcf515f6`) | Mesolithic human, hg19, ENA [ERR577494](https://www.ebi.ac.uk/ena/browser/view/ERR577494) (PRJEB6272, Lazaridis et al. 2014), file `Motala3.bam` |
| Synthetic 1 M / 10 M / 50 M | 1,000,000 / 10,000,000 / 50,000,000 | — | BER11_L_NOuser records repeated 4 times, then cut to N records. **Synthetic**; used only for the memory test. |

The ENA files were downloaded on 2026-11-11 and checked against ENA's md5. ENA data are shared openly under INSDC policy. No BAMs are committed.

**Why a subset of ERR1892229:** the file mixes USER-treated (`_user_`) and untreated (`NOuser`) libraries. In the first 2 M reads, every `_user_` read (1,514,462) had an MD tag that does not match its alignment: MD covers 5 more reference bases than the CIGAR (4 for one library). pmdshodh stops on such reads with an error, as CLAUDE.md requires. All 485,538 `NOuser` reads had valid MD tags. The subset keeps the three `NOuser` read groups:

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

Reproduce with `RUNS=5 INPUTS="BER11_L_NOuser Motala3" SYNTH="1M 10M 50M" bash docs/benchmarks/run.sh`, then `python3 docs/benchmarks/check.py WORK RESULTS BER11_L_NOuser Motala3`.

## Results

### Filtering at PMD ≥ 3: PMDtools vs pmdshodh

| | BER11_L_NOuser (14.0 M reads) | Motala3 (23.6 M reads) |
|---|---|---|
| PMDtools pipeline | 1,579.1 s (26.3 min) | 2,129.0 s (35.5 min) |
| pmdshodh `filter`, 1 thread | 51.0 s (**31.0×**) | 76.7 s (**27.8×**) |
| pmdshodh `filter`, 8 threads | 10.2 s (155×) | 13.4 s (159×) |
| pmdshodh `filter`, 12 threads | 10.2 s (155×) | 12.8 s (166×) |
| pmdshodh `filter --pmdtools-compat`, 1 thread | 69.6 s (22.7×) | 98.1 s (21.7×) |
| pmdshodh `filter --pmdtools-compat`, 8 threads | 13.0 s (121×) | 16.9 s (126×) |

- Speed-ups are relative to the PMDtools pipeline on the same input.
- The 1-thread rows are the like-for-like comparison of the scoring itself. The PMDtools pipeline uses about 1 CPU-second per wall-second (1,566 s CPU for BER11_L_NOuser).
- PMDtools' 5 runs ranged from 1,490 to 1,630 s on BER11_L_NOuser and from 2,077 to 2,207 s on Motala3. Min–max ranges for every configuration are in `summary.tsv`. The widest pmdshodh range was `score -@ 8` on BER11_L_NOuser (11.6–18.1 s, one slow run).
- Compat mode is 28–37% slower than default on 1 thread. The speed work of PRs #9 and #10 helped the default path more than compat's reference reconstruction.

### Thread scaling (`pmdshodh filter`)

| `-@` | BER11_L_NOuser | Speed-up | CPU cores used | Motala3 | Speed-up | CPU cores used |
|---|---|---|---|---|---|---|
| 1 | 51.0 s | 1.00× | 1.2 | 76.7 s | 1.00× | 1.1 |
| 2 | 20.5 s | 2.49× | 3.1 | 27.9 s | 2.75× | 2.9 |
| 4 | 14.1 s | 3.62× | 5.7 | 18.1 s | 4.24× | 5.3 |
| 8 | 10.2 s | 5.00× | 9.2 | 13.4 s | 5.72× | 8.9 |
| 12 | 10.2 s | 5.00× | 9.4 | 12.8 s | 5.99× | 9.5 |

- "CPU cores used" is median CPU time (user + system) divided by median wall time.
- **`-@ N` uses more than N cores.** Above 1, reading, scoring and writing run in separate threads, in addition to the N scoring and BGZF threads. This is why going from 1 to 2 is more than 2× faster. Even `-@ 1` uses slightly more than one core.
- Scaling flattens after 8 threads. The CPU has 4 performance cores, and the reader and writer are each a single thread.

### Per-read scores: pmdshodh vs metaDMG

| | BER11_L_NOuser | Motala3 |
|---|---|---|
| metaDMG `pmd`, `-@ 1` | **54.4 s** | **81.2 s** |
| metaDMG `pmd`, `-@ 8` | 51.9 s | 81.6 s |
| pmdshodh `score`, 1 thread | 58.9 s | 87.4 s |
| pmdshodh `score`, 8 threads | 12.5 s | 16.5 s |

- **On one thread, metaDMG is still faster**, by 8% on both inputs (1.08× on each). The two tools alternated run by run, and metaDMG was faster in every one of the 5 pairs on both inputs.
- metaDMG's `-@` option made almost no difference to its time. Its output was identical for `-@ 1` and `-@ 8`.
- pmdshodh with 8 threads is 4.2× faster than metaDMG on BER11_L_NOuser and 4.9× faster on Motala3.
- The two tools don't compute the same score, so this compares speed only.
- `score` takes longer than `filter` because it writes a line for every read (14–24 M lines), whereas `filter` writes only the kept reads.

### Memory

Peak resident memory (`/usr/bin/time -v`, median of 5):

| | 1 thread | 8 threads |
|---|---|---|
| pmdshodh `filter`, BER11_L_NOuser / Motala3 | 12 MB / 11 MB | 28 MB / 25 MB |
| pmdshodh `filter`, synthetic 1 M / 10 M / 50 M reads | 12 / 12 / 12 MB | 27 / 27 / 28 MB |
| pmdshodh `score`, BER11_L_NOuser / Motala3 | 11 MB / 10 MB | 22 MB / 21 MB |
| metaDMG `pmd` | 10–11 MB | 10 MB |
| PMDtools pipeline, largest single process | 15 MB | — |

- **pmdshodh's memory does not grow with input size.** With several threads it holds up to three batches of records at once, so memory is higher than on one thread, but it is bounded.
- On the synthetic inputs, its time is also linear. On 1 thread it took 4.00, 4.06 and 4.04 s per million reads; on 8 threads, 0.70, 0.68 and 0.68 s.
- For the PMDtools pipeline, `/usr/bin/time` reports the largest single process, not the pipeline's total.

## Correctness checks

Run on every benchmark input with `check.py`. All checks passed.

| Check | BER11_L_NOuser | Motala3 |
|---|---|---|
| Per-read scores: PMDtools `--printDS` vs pmdshodh `score --pmdtools-compat` (tolerance 1e-6) | 14,006,183 reads, largest difference 6.0e-14, 0 above tolerance | 23,622,338 reads, largest difference 5.0e-14, 0 above tolerance |
| Kept reads: PMDtools vs `filter --pmdtools-compat` (whole records, in order) | 1,031,193, identical including tags | 301,867, identical including tags |
| `filter` output identical for 1, 2, 4, 8 and 12 threads | 1,029,380 records each | 301,722 records each |
| `filter --pmdtools-compat` and `score` identical for 1 and 8 threads | pass | pass |

**Default mode vs PMDtools** (kept reads at PMD ≥ 3):
- BER11_L_NOuser: 1,907 reads kept only by PMDtools, 94 only by pmdshodh, 1,029,286 by both.
- Motala3: 242, 97 and 301,625.

These counts are the same as in the preliminary round (2026-11-12, commit `a98ae07`). In that round, every read whose default score differed from PMDtools' was traced to a documented departure (docs/design.md §5), with none unexplained. That classification was not rerun on this round's outputs.

| Cause of the difference (preliminary round) | BER11_L_NOuser | of which cross PMD 3 | Motala3 | of which cross PMD 3 |
|---|---|---|---|---|
| CIGAR reconstruction fix (`I`/`S` with `D`/`=`) | 3,983 | 1,960 | 3,757 | 339 |
| Lowercase MD bases | 336 | 38 | 0 | 0 |
| Lowercase MD + reverse-read shift | 3 | 3 | 0 | 0 |
| CIGAR fix + lowercase MD | 2 | 0 | 0 | 0 |
| **Unexplained** | **0** | | **0** | |

## Compared with the preliminary round

A preliminary round ran 3 times per configuration on commit `a98ae07` (2026-11-12). Its single-thread `filter` medians were 84.0 s (BER11_L_NOuser) and 112.7 s (Motala3), against 51.0 s and 76.7 s here. Its 8-thread medians were 18.6 s and 24.7 s, against 10.2 s and 13.4 s.

The two rounds weren't interleaved and the earlier one timed with the WSL wall clock, so this is not a controlled comparison. The controlled before/after measurements of the speed work are in the CHANGELOG. The preliminary round's PMDtools compat-mode check initially failed on 339 reads; PR #6 fixed the cause before this round.

## Departures from the plan

- **No PMDtools run with `samtools calmd`.** All inputs already have MD tags, and this run would have needed two reference genomes.
- **`--printDS` and the compat `score` were run once each.** They exist for the parity check.
- **`--ss` is not benchmarked.** No single-stranded library has been tested.
