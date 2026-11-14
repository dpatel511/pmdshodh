# pmdshodh

**Fast, maintained post-mortem damage (PMD) scoring and filtering for ancient DNA alignments.**

pmdshodh scores every read in a SAM/BAM file for ancient DNA damage, keeps or removes reads by score, and summarises damage patterns. It uses the scoring model of [PMDtools](https://github.com/pontussk/PMDtools) 0.60, fixes PMDtools' known errors by default (`--pmdtools-compat` reproduces PMDtools exactly), reads and writes alignments directly, and runs on multiple threads.

> **Status:** alpha (v0.1.0, not yet released). Scores have been validated against PMDtools on real and synthetic data. Performance benchmarks are in [docs/benchmarks.md](docs/benchmarks.md); packages and CRAM support are still to come.

## Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Commands](#commands)
- [Options](#options)
- [Output](#output)
- [Compatibility with PMDtools](#compatibility-with-pmdtools)
- [Limitations](#limitations)
- [Development](#development)
- [License](#license)

## Overview

Ancient DNA molecules carry characteristic damage. Cytosine deamination appears as C→T substitutions near the 5′ end of reads, and as G→A near the 3′ end in double-stranded libraries. The PMD score of [Skoglund et al. (2014)](https://doi.org/10.1073/pnas.1318934111) is a per-read log-likelihood ratio comparing a damage model with a no-damage model. Reads scoring above a threshold (commonly 3) are likely to be ancient rather than modern contamination.

pmdshodh implements this score as a single program, with no `samtools calmd | python | samtools view` pipeline. It does not classify taxa, estimate contamination, call genotypes or rescale base qualities.

## Features

- **PMDtools-compatible scores.** It uses the same likelihood model, defaults and options (`--threshold`, `--ss`, `--CpG`, `--UDGhalf`, `-m`, `-q`, …). Single-stranded (`--ss`) scores follow the paper's model and are lower than PMDtools'; see [Differences that change scores](#differences-that-change-scores).
- **Direct SAM/BAM input and output.** Kept records are written unchanged and in input order, with an `@PG` header line.
- **Multi-threaded** scoring and BGZF compression. Results are identical for any thread count.
- **JSON run summary** with read counts, a score histogram and damage tables, designed for MultiQC.
- **Robust input handling.** It never drops reads silently: every skipped read is counted by reason, and invalid input stops with a clear error.
- **Correct by default, exact on request.** Known PMDtools bugs are fixed, and `--pmdtools-compat` reproduces PMDtools exactly for comparisons.

## Installation

pmdshodh is currently installed from source. It needs Rust 1.85 or newer (tested with 1.99.1).

```bash
git clone https://github.com/dpatel511/pmdshodh.git
cd pmdshodh
cargo build --release
./target/release/pmdshodh --version
```

**Input requirement:** alignments must carry an **MD tag**. If your aligner does not write one, add it with `samtools calmd`:

```bash
samtools calmd -b input.bam reference.fa > input.md.bam
```

## Quick start

```bash
# Keep reads with PMD score >= 3 (double-stranded library, no UDG treatment)
pmdshodh filter --threshold 3 sample.bam -o sample.pmd3.bam

# Add the score to each kept read, write a run summary, use 8 threads
pmdshodh filter --threshold 3 --tag --summary sample.json --threads 8 sample.bam -o sample.pmd3.bam

# Single-stranded, UDG-half treated library
pmdshodh filter --threshold 3 --library ss --treatment half sample.bam -o sample.pmd3.bam

# Per-read scores as a table
pmdshodh score sample.bam -o sample.scores.tsv

# Damage summary only
pmdshodh damage sample.bam -o sample.damage.json
```

## Commands

| Command | Purpose | Output |
|---|---|---|
| `filter` | Keep reads by PMD score | SAM or BAM, chosen by the `-o` extension (`-o -` writes BAM to stdout) |
| `score` | Score every read | Tab-separated table: `read_name`, `flag`, `pmds`, `n_sites` |
| `damage` | Summarise damage without filtering | JSON run summary |

`filter` keeps reads with score **≥ `--threshold`** and/or **< `--upper-threshold`**; at least one is required.

Run `pmdshodh <command> --help` for every option.

## Options

PMDtools option names are accepted as aliases (shown in the second column).

### Library and treatment

| Option | PMDtools | Default | Description |
|---|---|---|---|
| `--library ds\|ss` | `--ss` | `ds` | Double- or single-stranded library. Never inferred from the data. |
| `--treatment none\|half\|full` | `--UDGminus`, `--UDGhalf`, `--UDGplus` / `--CpG` | `none` | `half`: CpG sites plus the terminal base; `full`: CpG sites only. |

### Read and base selection

| Option | PMDtools | Default | Description |
|---|---|---|---|
| `-m, --min-mapq N` | `-m`, `--requiremapq` | 0 | Skip reads below this mapping quality. |
| `-q, --min-baseq N` | `-q`, `--requirebaseq` | 0 | Ignore bases below this base quality. |
| `--no-clips` | `--noclips` | off | Skip reads with soft or hard clips. |
| `--no-indels` | `--noindels` | off | Skip reads with insertions or deletions. |

### Model parameters

| Option | PMDtools | Default |
|---|---|---|
| `--pmd-p` | `--PMDpparam` | 0.3 |
| `--pmd-c` | `--PMDconstant` | 0.01 |
| `--pi-ancient` | `--polymorphism_ancient` | 0.001 |
| `--pi-contamination` | `--polymorphism_contamination` | 0.001 |

### Output and performance

| Option | Applies to | Description |
|---|---|---|
| `-t, --threshold X` | `filter` | Keep reads with score ≥ X. |
| `--upper-threshold X` | `filter` | Keep reads with score < X. |
| `--tag`, `--tag-name pm\|DS` | `filter` | Add the score to kept reads as `pm:f` (32-bit float, default) or `DS:Z` (3 decimals, as PMDtools `--writesamfield`). |
| `--summary FILE` | `filter`, `score` | Also write the JSON run summary. |
| `--range N` | all | Positions from each read end in damage tables (default 30). |
| `-n, --max-reads N` | `damage` | Stop after N reads. |
| `-@, --threads N` | all | Worker threads (default 1). |
| `--pmdtools-compat` | all | Reproduce PMDtools 0.60 exactly, including its bugs. |
| `--quiet` | all | Suppress the summary line on stderr. |

## Output

### Alignments (`filter`)

- Kept records are unchanged (identical when decoded) apart from the optional score tag.
- Input order is preserved.
- The input header is copied, and one `@PG` line records the command.
- Output files are written under a temporary name and renamed only on success, so a failed run never leaves a partial file.

### JSON run summary

A single file per run, versioned with `format_version` (currently `1`). It contains:

| Section | Contents |
|---|---|
| `tool`, `command`, `mode`, `model`, `filters` | Version, command line, library and treatment, model parameters, thresholds |
| `reads` | Input, scored, kept and removed counts; skipped reads by reason; flag counts (duplicate, secondary, supplementary, …) |
| `score_histogram`, `length_histogram` | Score distribution (0.5-wide bins from −20 to 40), read lengths |
| `damage` | Counts of read bases A/C/G/T at reference C and G, by position from each molecule end; first-base C→T rate with standard error |

### Exit status and streams

| Exit status | Meaning |
|---|---|
| `0` | Success, including empty input |
| `1` | Runtime failure, e.g. the output cannot be written |
| `2` | Invalid options or invalid input, e.g. a mapped read without an MD tag (the read is named in the message) |

Results go to stdout only when requested (`-o -`, or `score`/`damage` without `-o`). The run summary line and errors go to stderr.

## Compatibility with PMDtools

### Validation

pmdshodh was built from the published model and black-box testing of PMDtools 0.60; PMDtools' source code was not used. Every read in each test set was compared:

| Test set | Reads | Largest score difference |
|---|---|---|
| Single-site score terms (ds and ss, base quality 0–93, 6 parameter sets) | 5,412 | 1.1 × 10⁻¹³ |
| Synthetic reads with clips, indels and `=`/`X` CIGARs, both strands, 8 modes | 1,500 per mode | 8.3 × 10⁻¹⁴ ¹ |
| nf-core/eager test BAMs (mammoth mtDNA, human 1240k capture), ds, `--CpG`, `--UDGhalf`, `--ss` | 13,670 per mode | 3.0 × 10⁻¹⁴ ¹ |
| Damage tables vs `pmdtools --deamination --range 30` (no UDG and full UDG) | 1,500 | identical at all 30 positions |

¹ With `--pmdtools-compat`. Without it, scores differ only where PMDtools has a bug (below).

Methods and full results: [docs/blackbox.md](docs/blackbox.md).

### Differences that change scores

pmdshodh corrects these PMDtools behaviours by default. Use `--pmdtools-compat` to reproduce them.

| PMDtools 0.60 behaviour | Effect | pmdshodh default |
|---|---|---|
| The reference is rebuilt wrongly for reads that combine an insertion or soft clip with a deletion or `=` CIGAR operations | Deleted bases are scored against later read bases, and positions shift | Reference rebuilt correctly from CIGAR and MD |
| CpG context wraps around the read ends in `--CpG` and `--UDGhalf` modes | Changed 6% of read scores (up to 1.62) on the human UDG test BAM | No wrap-around |
| `--ss` counts every C site twice (the double-stranded term plus the paper's single-stranded term, SI Eq. S1) | ss scores are about twice as large as the model gives, and lopsided: a terminal C→T scores 9.68 at the 5′ end and 6.34 at the 3′ end (Q40), so a threshold of 3 is far looser than for ds libraries | One term per C site (Eq. S1): 4.65 at either end. **ss scores are lower than PMDtools'; re-check thresholds tuned on PMDtools ss scores** |
| `--upperthreshold` keeps reads with no informative sites (score 0) | Reads with no C/G evidence at all are output as "non-damaged" | Not kept by `--upper-threshold` (still counted as `no_informative_sites`) |
| Lowercase reference bases in MD (outside the SAM spec) are never counted as C or G | 339 of 14.0 M reads scored differently on a public horse BAM (ENA ERR1892229) | Read as uppercase |
| On reverse-strand reads, MD characters other than `ACGTNRY`/`acgt` (e.g. lowercase `n`) are dropped when the reference is reverse-complemented | The reference shifts against the read; scores up to 15.7 too high | No shift |

**In CpG modes, expect pmdshodh's default results to differ from PMDtools on some reads.**

### Differences that prevent crashes or data loss

| Situation | PMDtools 0.60 | pmdshodh |
|---|---|---|
| Mapped read without an MD tag | Dropped silently | Error (exit 2) |
| MD string without a `0` after a deletion (as written by MALT) | Misparsed silently | Parsed correctly |
| Hard-clipped reads | Dropped | Scored |
| Existing `DS:Z` tag, or `DS:Z:` text in qualities | Crash | Handled |
| Empty input | Crash | Empty output, exit 0 |
| Damage site ≥ 1,000 bp from a read end, or hundreds of mismatches in one read | Crash | Scored |
| `--first` (first-base damage rate with standard error) | Crashes in the 0.60 release | Reported in the JSON summary |

### Behaviour kept for compatibility

These are part of PMDtools' score as computed, and pmdshodh keeps them so that scores mean the same:

- The no-damage model assumes a background damage rate of 0.001.
- Reads are scored regardless of SAM flags (duplicates, secondary and supplementary alignments, paired reads). Remove unwanted reads beforehand; the summary counts them.

The complete list, with reasoning: [docs/design.md](docs/design.md) §5.

## Limitations

- CRAM input is not yet supported.
- Benchmarks cover one laptop and two double-stranded libraries: [docs/benchmarks.md](docs/benchmarks.md). On one thread, metaDMG `pmd` scores reads about 8% faster than `pmdshodh score`.
- No real single-stranded library has been tested; single-stranded parity rests on synthetic reads.
- UDG-half damage *tables* do not yet match PMDtools (UDG-half *scores* do).

## Development

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
```

- **Reference data:** the tests compare against values produced by PMDtools 0.60, stored in `tests/data/`. The scripts in `tests/oracle/` regenerate them; they need PMDtools 0.60 on Linux or WSL.
- **Readable model:** `tests/oracle/pmd_reference.py` is a short Python version of the scoring model.
- **Background:** research, design decisions and black-box findings are in [docs/](docs/).

## License

[MIT](LICENSE).
