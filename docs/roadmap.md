# Roadmap — agreed 2026-11-11

**Positioning (narrowed 2026-11-15): the maintained per-read damage filter.**
- **Scores:** the same as PMDtools, plus documented fixes.
- **Speed and install:** much faster, and installable with pip and Bioconda.
- **Input:** one step from BAM or CRAM with `-r`, with no separate `samtools calmd`.
- **QC:** damage QC built in, for MultiQC.
- **Drop-in:** it is a drop-in replacement for PMDtools in pipelines such as nf-core/eager.
- **Not a general damage profiler:** pmdshodh does not compete with DamageProfiler, mapDamage2 or AdDeam on that. Its damage tables and the planned HTML report are QC around filtering.
- **What sets it apart from PMDtools:** fitted damage parameters (milestone 6).

The 2026-11-11 plan also aimed to compete directly with metaDMG, DamageProfiler and mapDamage. The landscape checked on 2026-11-15:
- Only PMDtools and ngsBriggs (a preprint, which needs an external contamination estimate) filter individual reads.
- metaDMG `pmd` scores reads but does not filter them.
- PyDamage works on contigs, and bamRefine masks damaged positions.

Each milestone ends with passing tests and, where a claim is made, measured evidence (CLAUDE.md §5, §7).

## Principles
- **The default score stays PMDtools-equivalent**, so published thresholds (e.g. PMD ≥ 3) keep their meaning. Alternative models and fitted parameters are opt-in and labelled in every output.
- **Scope** stays per-read PMD scoring, filtering and damage summaries (CLAUDE.md §7). Taxonomy, contamination estimates, genotyping and quality rescaling remain out of scope unless explicitly agreed.
- **No performance or accuracy claim without measurements** (hardware, versions, commands, repeated runs).

## Milestones

| # | Milestone | Content | Done when |
|---|---|---|---|
| 1 | **CI** | GitHub Actions on Linux and macOS: `fmt`, `clippy`, `test`; `samtools quickcheck` on BAMs written by the tests | Green on every PR |
| 2 | **Benchmarks** | [benchmark-plan.md](benchmark-plan.md): vs the PMDtools pipeline and metaDMG `pmd`, thread scaling, memory | `docs/benchmarks.md` (raw logs kept locally, not committed) |
| 3 | **Easy installation → v0.1.0** | Prebuilt binaries (Linux x86_64/aarch64, macOS, Windows) on GitHub Releases; crates.io; Bioconda recipe; pip (maturin) considered | Users can install without Rust; names re-checked before publishing |
| 4 | **Reference-based scoring and CRAM** | `-r/--reference ref.fa`: score from the reference instead of MD, count stale MD tags, check the header against the reference; CRAM input (CRAM output deferred). Design: docs/design.md §11. **Done 2026-11-16** (v0.2.0): parity results in CHANGELOG.md | Parity: scores with `-r` equal scores after `samtools calmd` on the same data; `--pmdtools-compat -r` equals `calmd \| pmdtools` |
| 5 | **Damage report** | HTML report (`pmdshodh report`): C→T/G→A by position, score histogram, read counts, first-base rate; generated from the JSON summary. PNG deferred (D2, 2026-11-20: it needs a rasterising dependency). **Done 2026-11-20** (v0.3.0): validated on the eager test BAMs, docs/report-proposal.md §9 | Report checked on the eager test BAMs; JSON stays the source of truth |
| 6 | **Fitted damage parameters (opt-in)** | Estimate P and C (and optionally the background) from the sample's damage profile, then score with them; record fitted values and fit quality in the summary | Validated on simulated data with known parameters and on real ds/UDG-half samples; default unchanged |
| 7 | **Integration** | MultiQC module; nf-core module `pmdshodh/filter`; proposal to nf-core/eager (contact needs approval) | Modules accepted upstream |
| later | **Per-reference summaries** | Damage and score summaries per reference sequence for metagenomic BAMs | Agreed again before starting |

## Open items carried forward
- Real single-stranded library for parity testing: done without UDG treatment (2026-11-17, CHANGELOG.md) and USER-treated (I19359, ENA ERR7195621; compat equals PMDtools in ss and ss `--UDGhalf`, docs/audit-0.3.md §8.1, 2026-11-21).
- PMDtools UDG-half damage *table* rule: found 2026-11-20 (docs/blackbox.md 1.6d); `--pmdtools-compat` now matches it.
- eager 2.5.3 `--CpG` for non-UDG libraries: confirm by a run before reporting (your decision).
- Contact with metaDMG developers (your decision).
- Whether `filter` should regenerate MD/NM tags for CRAM input, as samtools does on decoding (your decision; design.md §11 item 8).
