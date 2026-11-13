# Roadmap — agreed 2026-11-11

pmdshodh aims to be both:
- a **drop-in replacement** for PMDtools in pipelines such as nf-core/eager;
- a **standalone tool** that people choose over PMDtools, metaDMG, DamageProfiler or mapDamage when working on their own data.

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
| 4 | **Reference-based scoring and CRAM** | `-r/--reference ref.fa`: score reads without MD tags by reading the reference; CRAM input | Parity: scores with `-r` equal scores from `calmd` + MD on the same data |
| 5 | **Damage report** | HTML (and PNG) report: C→T/G→A by position, score histogram, read counts, first-base rate; generated from the JSON summary | Report checked on the eager test BAMs; JSON stays the source of truth |
| 6 | **Fitted damage parameters (opt-in)** | Estimate P and C (and optionally the background) from the sample's damage profile, then score with them; record fitted values and fit quality in the summary | Validated on simulated data with known parameters and on real ds/UDG-half samples; default unchanged |
| 7 | **Integration** | MultiQC module; nf-core module `pmdshodh/filter`; proposal to nf-core/eager (contact needs approval) | Modules accepted upstream |
| later | **Per-reference summaries** | Damage and score summaries per reference sequence for metagenomic BAMs | Agreed again before starting |

## Open items carried forward
- Real single-stranded library for parity testing.
- PMDtools UDG-half damage *table* rule (scores already match).
- eager 2.5.3 `--CpG` for non-UDG libraries: confirm by a run before reporting (your decision).
- Contact with metaDMG developers (your decision).
