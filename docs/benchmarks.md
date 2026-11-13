# Benchmarks — preliminary, 2026-11-12

**Status: preliminary.** Each configuration ran **3 times**, not the 5 that [benchmark-plan.md](benchmark-plan.md) asks for, on one laptop. The pmdshodh build was commit `a98ae07`, which predates the scoring changes of PRs #6 and #7. Those PRs don't change double-stranded default scoring, so the default timings should still apply; that isn't measured. A 5-run benchmark on the final build will replace these numbers.

All numbers are medians of 3 runs. `run.sh` writes the raw logs (commands, `/usr/bin/time -v` output, stderr, and a per-run `times.tsv` with the min–max range), and `check.py` writes `summary.tsv` and `checks.tsv`. The logs for this round are kept with the benchmark machine and are not committed. The scripts are [benchmarks/run.sh](benchmarks/run.sh) and [benchmarks/check.py](benchmarks/check.py).

## Setup

| Item | Value |
|---|---|
| CPU | Intel Core i5-13500H laptop: 12 cores (4 performance + 8 efficiency), 16 threads |
| OS | WSL 2 (Linux 6.18.58.1-microsoft-standard-WSL2), Ubuntu 24.04, on Windows 11; WSL sees 7.6 GB RAM |
| Power | On mains, Windows "Balanced" plan, sleep off |
| Storage | Inputs and outputs on the WSL ext4 disk (`~/bench`), not on `/mnt` |
| Cache | Inputs read once with `cat` before each block of runs (warm cache) |
| pmdshodh | 0.1.0, commit `a98ae074`, `cargo build --release` (Rust 1.99.1, default release profile) |
| PMDtools | 0.60 (Bioconda tag archive, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`), Python 3.12.3 |
| samtools | 1.25 (htslib 1.25), Bioconda |
| metaDMG-cpp | Bioconda package `metadmg` 0.5.3; the binary reports itself as 0.5.0 |

Timings run one after another, never in parallel. The whole run took 5 h 45 min (2026-11-12 20:09 to 2026-11-13 01:54 UTC).

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

Reproduce with `RUNS=3 INPUTS="BER11_L_NOuser Motala3" SYNTH="1M 10M 50M" bash docs/benchmarks/run.sh`, then `python3 docs/benchmarks/check.py WORK RESULTS BER11_L_NOuser Motala3`.

## Results

### Filtering at PMD ≥ 3: PMDtools vs pmdshodh

| | BER11_L_NOuser (14.0 M reads) | Motala3 (23.6 M reads) |
|---|---|---|
| PMDtools pipeline | 1,596.8 s (26.6 min) | 2,161.4 s (36.0 min) |
| pmdshodh `filter`, 1 thread | 84.0 s (**19.0×**) | 112.7 s (**19.2×**) |
| pmdshodh `filter`, 8 threads | 18.6 s (85.8×) | 24.7 s (87.5×) |
| pmdshodh `filter`, 12 threads | 17.5 s (91.2×) | 22.3 s (96.9×) |
| pmdshodh `filter --pmdtools-compat`, 1 thread | 94.5 s | 131.9 s |

Speed-ups are relative to the PMDtools pipeline on the same input. The 1-thread rows are the like-for-like comparison of the scoring itself. The PMDtools pipeline uses about 1 CPU-second per wall-second (1,636 s CPU for BER11_L_NOuser). Compat mode is 12–17% slower than default on 1 thread.

### Thread scaling (`pmdshodh filter`)

| Threads | BER11_L_NOuser | Speed-up | Motala3 | Speed-up |
|---|---|---|---|---|
| 1 | 84.0 s | 1.00× | 112.7 s | 1.00× |
| 2 | 40.0 s | 2.10× | 52.4 s | 2.15× |
| 4 | 27.9 s | 3.01× | 30.1 s | 3.74× |
| 8 | 18.6 s | 4.52× | 24.7 s | 4.56× |
| 12 | 17.5 s | 4.80× | 22.3 s | 5.05× |

- Going from 1 to 2 threads is more than 2× faster. `-@` sets both the scoring threads and the BGZF compression workers, which probably explains it; I haven't measured that.
- Scaling flattens after 8 threads. The CPU has 4 performance cores, and reading and writing records runs on a single thread in the current code.

### Per-read scores: pmdshodh vs metaDMG

| | BER11_L_NOuser | Motala3 |
|---|---|---|
| metaDMG `pmd`, `-@ 1` | **56.2 s** | **84.8 s** |
| metaDMG `pmd`, `-@ 8` | 56.5 s | 81.3 s |
| pmdshodh `score`, 1 thread | 84.3 s | 121.0 s |
| pmdshodh `score`, 8 threads | 25.0 s | 38.3 s |

- **On one thread, metaDMG is faster**: 1.50× on BER11_L_NOuser and 1.43× on Motala3.
- metaDMG's `-@` option made no difference to its time. Its output was identical for `-@ 1` and `-@ 8`.
- pmdshodh with 8 threads is 2.2× faster than metaDMG.
- The two tools don't compute the same score, so this compares speed only.

### Memory

Peak resident memory (`/usr/bin/time -v`, median of 3):

| | 1 thread | 8 threads |
|---|---|---|
| pmdshodh `filter`, BER11_L_NOuser / Motala3 | 11 MB / 11 MB | 16 MB / 14 MB |
| pmdshodh `filter`, synthetic 1 M / 10 M / 50 M reads | 11 / 11 / 11 MB | 14 / 15 / 16 MB |
| metaDMG `pmd` | 10–11 MB | 10–11 MB |
| PMDtools pipeline, largest single process | 15 MB | — |

- **pmdshodh's memory does not grow with input size.** On the synthetic inputs, its time is also linear: 5.8, 5.75 and 5.76 s per million reads on 1 thread, and 1.30, 1.26 and 1.27 s on 8 threads.
- For the PMDtools pipeline, `/usr/bin/time` reports the largest single process, not the pipeline's total.

## Correctness checks

Run on every benchmark input with `check.py`. The tables below are taken from its output for the original run and for the rerun after the fix.

| Check | BER11_L_NOuser | Motala3 |
|---|---|---|
| Per-read scores: PMDtools `--printDS` vs pmdshodh `score --pmdtools-compat` (tolerance 1e-6) | 14,006,183 reads, largest difference 6.0e-14 | 23,622,338 reads, largest difference 5.0e-14 |
| Kept reads: PMDtools vs `filter --pmdtools-compat` (whole records, in order) | 1,031,193, identical including tags | 301,867, identical including tags |
| `filter` output identical for 1, 2, 4, 8 and 12 threads; `score` for 1 and 8 | pass | pass |

**This check failed on the original run** for BER11_L_NOuser:
- 339 reads had different scores (largest difference 15.7), and 35 more reads were kept than PMDtools keeps.
- The cause was two PMDtools behaviours with unusual MD tags that compat mode did not yet reproduce: lowercase reference bases, and characters dropped on reverse-strand reads (docs/blackbox.md).
- PR #6 fixed compat mode. The compat-mode outputs were then regenerated with that code and checked against the same PMDtools outputs, giving the table above. Those regenerated runs were not timed.

**Default mode vs PMDtools** (kept reads at PMD ≥ 3):
- BER11_L_NOuser: 1,907 reads kept only by PMDtools, 94 only by pmdshodh, 1,029,286 by both.
- Motala3: 242, 97 and 301,625.

Every read whose default score differs from PMDtools' comes from a documented departure (docs/design.md §5). None are unexplained:

| Cause of the difference | BER11_L_NOuser | of which cross PMD 3 | Motala3 | of which cross PMD 3 |
|---|---|---|---|---|
| CIGAR reconstruction fix (`I`/`S` with `D`/`=`) | 3,983 | 1,960 | 3,757 | 339 |
| Lowercase MD bases | 336 | 38 | 0 | 0 |
| Lowercase MD + reverse-read shift | 3 | 3 | 0 | 0 |
| CIGAR fix + lowercase MD | 2 | 0 | 0 | 0 |
| **Unexplained** | **0** | | **0** | |

## Departures from the plan

- **3 runs per configuration instead of 5**, for this preliminary round.
- **No PMDtools run with `samtools calmd`.** All inputs already have MD tags, and this run would have needed two reference genomes.
- **`--printDS` and the compat `score` were run once each.** They exist for the parity check.
- **`--ss` is not benchmarked.** No single-stranded library has been tested.

## Speed notes for the next round

These were measured on the first 2 M reads of BER11_L_NOuser, 1 thread, 3 runs each, on an experimental copy. They're not in the released code yet.

| | Time |
|---|---|
| | Median time |
|---|---|
| htslib decoding (`samtools view -c`) | 1.08 s |
| pmdshodh `score` (`a98ae07`) | 12.46 s |
| with precomputed site terms (output byte-identical) | 9.19 s |
| … and damage counting turned off | 6.87 s (7.01 s in a second set of 3 runs) |
| … plus fat LTO, `codegen-units = 1` | 6.85 s: no clear gain; binary 3.95 → 2.94 MB |

Turning damage counting off shows how much it costs at most. It is needed for the summary, so the real gain will be smaller.
