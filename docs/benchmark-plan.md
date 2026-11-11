# Benchmark plan — 2026-11-11

Status: **plan only; nothing measured yet.** No speed or memory claim may be made until this plan has been run and the results are written to `docs/benchmarks.md` with raw logs (CLAUDE.md §7).

## 1. Questions
1. **PMDtools vs pmdshodh.** How long does a PMD-filtered BAM take with PMDtools as nf-core runs it, against `pmdshodh filter` on the same input and threshold?
2. **Scaling.** How does pmdshodh scale with `-@` (1, 2, 4, 8, 12 threads), and where does it stop scaling?
3. **Score-only.** How does `pmdshodh score` compare with `metaDMG-cpp pmd`, which prints per-read scores, single-threaded?
4. **Memory.** What is the peak resident memory, and does it stay flat as input size grows (bounded-memory design)?
5. **Correctness at scale.** On every benchmark input, do the PMDtools and pmdshodh outputs keep the same reads? A speed number for a wrong result is worthless.

## 2. Tools and exact versions to record
| Tool | Version and source | How it is run |
|---|---|---|
| PMDtools | 0.60, Bioconda tag archive (md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`) | `samtools calmd -b in.bam ref.fa \| samtools view -h - \| python3 pmdtools.0.60.py --threshold 3 --header \| samtools view -b -o out.bam -`, the nf-core `pmdtools/filter` pipeline. Also a variant without `calmd` when the BAM already has MD, to separate PMDtools' own cost. |
| samtools | version recorded (`samtools --version`) | as above |
| metaDMG-cpp | release or commit recorded | `metaDMG-cpp pmd in.bam > scores.tsv` (`-@` is ignored; single-threaded) |
| pmdshodh | git commit recorded, **release build** (`cargo build --release`) | `pmdshodh filter -t 3 -@ N in.bam -o out.bam` and `pmdshodh score -@ N in.bam -o scores.tsv` |
| Python | version recorded | for PMDtools |

**Environment:** all tools run inside WSL Ubuntu 24.04 on the same machine, with input and output on the Linux filesystem (`~/bench`, not `/mnt/c` or `/mnt/e`, whose cross-OS file access would dominate the timings).

**Setup this needs** (each needs your go-ahead):
- **samtools and metaDMG-cpp in WSL.** Proposed: micromamba from Bioconda into a dedicated environment.
- **Rust inside WSL**, to build a Linux release binary of pmdshodh.

## 3. Hardware (to be recorded with every result)
- CPU: Intel Core i5-13500H, 12 cores (4 performance + 8 efficiency), 16 threads, laptop part.
- RAM: 15.7 GB physical; WSL currently sees 7 GB.
- Disk: WSL ext4 virtual disk.

**Caveats that must appear next to any number:**
- Performance and efficiency cores differ in speed, so thread-scaling curves will bend.
- Laptops throttle under sustained load. Each run is started on mains power after an idle period, and the power plan is recorded.

## 4. Data
| Input | Purpose | Size | Status |
|---|---|---|---|
| nf-core eager test BAMs (JK2782 mammoth, JK2067 human) | sanity check of the harness; too small to time | 0.5–0.7 MB | already used for parity |
| JK2067 full (nf-core test-datasets) | small real human capture | 6.3 MB | needs approval (small) |
| **A public genome-scale aDNA BAM** with MD tags or a downloadable reference (e.g. an ENA or eager benchmarking sample) | realistic throughput | 1–5 GB | **to be chosen; size and source shown to you before any download** |
| Synthetic scaled inputs: real reads repeated to 1 M, 10 M and 50 M records | memory flatness and thread scaling | generated locally | labelled synthetic in all reports |

## 5. Protocol
1. **Build and check:** `cargo build --release`; record `pmdshodh --version` and the git commit.
2. **Warm-up:** one untimed run per tool and input, so the file sits in the OS cache. Results are reported as warm-cache.
3. **Timing:** `/usr/bin/time -v` around each command, recording wall time, user and system CPU, and maximum resident set size.
   - **At least 5 runs per configuration**, reporting the median and the min–max range. No single-run numbers.
4. **Configurations:**
   - PMDtools pipeline (with and without `calmd`);
   - metaDMG `pmd`;
   - pmdshodh `filter` and `score` at `-@` 1, 2, 4, 8, 12;
   - for pmdshodh: default and `--pmdtools-compat`.
5. **Correctness on every input:** compare the PMDtools filtered read names with pmdshodh's (`--pmdtools-compat`): reads compared, differing reads, largest score difference. Any mismatch stops the benchmark until explained.
6. **Raw data kept:** commands, `/usr/bin/time` output and tool versions go to `docs/benchmarks/raw/` (logs only, no BAMs).

## 6. Reporting rules
- Report what was measured: "pmdshodh filter, 8 threads, took X s (median of 5) vs Y s for the PMDtools pipeline, on input Z, on this hardware". No general "N× faster" headline without the input, hardware and thread count next to it.
- Show single-threaded results too. A multi-threaded tool against a single-threaded one is not a like-for-like comparison of the scoring itself.
- If pmdshodh is slower in any configuration, report it.

## 7. Decisions needed before running
1. Approve installing micromamba (Bioconda samtools and metaDMG-cpp) and Rust inside WSL.
2. Choose the genome-scale BAM. I will propose 2–3 candidates with exact sizes and licences.
