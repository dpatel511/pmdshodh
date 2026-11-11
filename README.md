# pmdshodh

Per-read **post-mortem damage (PMD) scoring and filtering** of ancient DNA alignments: a fast, maintained successor to [PMDtools](https://github.com/pontussk/PMDtools), written in Rust.

*pmdshodh* = PMD + Sanskrit *śodh*, "to purify".

> **Status: alpha (v0.1.0, unreleased).** Scores match PMDtools 0.60 on every read tested (see [Parity](#parity-with-pmdtools)). Not yet on Bioconda or crates.io. No speed has been measured yet; see the [benchmark plan](docs/benchmark-plan.md). The command line may still change before the first release.

## What it does

Ancient DNA carries characteristic damage: cytosine deamination shows up as C→T near the 5′ end of reads (and G→A near the 3′ end in double-stranded libraries). PMDtools ([Skoglund et al. 2014, *PNAS*](https://doi.org/10.1073/pnas.1318934111)) gives each read a likelihood-ratio **PMD score**. Ancient DNA studies keep reads above a threshold (often 3) to separate genuinely ancient molecules from modern contamination.

pmdshodh computes the same score and:

- **reads SAM and BAM directly**, with no `samtools calmd | python | samtools view` pipe;
- **writes the filtered BAM or SAM itself**, keeping records unchanged, preserving input order and adding an `@PG` line;
- **runs on several threads** (`-@`), with identical results for any thread count;
- **writes a JSON summary** (read counts, score histogram, damage tables) meant for a MultiQC module;
- **accepts PMDtools' option names** (`--threshold`, `--ss`, `--CpG`, `--UDGhalf`, `-m`, `-q`, …), so existing workflows need few changes;
- **fixes PMDtools' crashes and silent data loss**, and keeps an exact-compatibility switch for comparisons.

It does not classify taxa, estimate contamination, call genotypes or rescale base qualities.

## Installation

From source. This needs Rust 1.85 or newer for the 2024 edition; it has been tested with 1.99.1.

```bash
git clone https://github.com/dpatel511/pmdshodh.git
cd pmdshodh
cargo build --release
./target/release/pmdshodh --version
```

Input reads need an **MD tag**. If your aligner does not write one, add it first, e.g. `samtools calmd -b in.bam ref.fa > in.md.bam`.

## Usage

```bash
# Keep reads with PMD score >= 3 (the common threshold), double-stranded library
pmdshodh filter -t 3 sample.bam -o sample.pmd3.bam

# Same, writing the score on each kept read and a JSON summary, with 8 threads
pmdshodh filter -t 3 --tag --summary sample.pmd3.json -@ 8 sample.bam -o sample.pmd3.bam

# Single-stranded library, UDG-half treated
pmdshodh filter -t 3 --ss --UDGhalf sample.bam -o sample.pmd3.bam

# Per-read score table (TSV: read_name, flag, pmds, n_sites)
pmdshodh score sample.bam -o sample.scores.tsv

# Damage summary only (JSON): C→T / G→A by position, first-base rate with standard error
pmdshodh damage --range 30 sample.bam -o sample.damage.json
```

### Subcommands

| Command | Output |
|---|---|
| `filter` | Reads with PMDS ≥ `--threshold` and/or < `--upper-threshold`, as BAM or SAM (from the `-o` extension; `-o -` writes BAM to stdout). |
| `score` | Tab-separated per-read scores. |
| `damage` | JSON run summary only. |

### Main options (all subcommands)

| Option | PMDtools equivalent | Default |
|---|---|---|
| `--library ds\|ss` | `--ss` | `ds` (never guessed from the data) |
| `--treatment none\|half\|full` | `--UDGminus`, `--UDGhalf`, `--UDGplus`/`--CpG` | `none` |
| `-m, --min-mapq N` | `-m`, `--requiremapq` | 0 |
| `-q, --min-baseq N` | `-q`, `--requirebaseq` | 0 |
| `--no-clips`, `--no-indels` | `--noclips`, `--noindels` | off |
| `--pmd-p`, `--pmd-c`, `--pi-ancient`, `--pi-contamination` | `--PMDpparam`, `--PMDconstant`, `--polymorphism_ancient`, `--polymorphism_contamination` | 0.3, 0.01, 0.001, 0.001 |
| `-@, --threads N` | – | 1 |
| `--pmdtools-compat` | – | off (see below) |
| `--summary FILE`, `--range N` | – (`--range` for tables) | none, 30 |

`filter` also has `--tag` (adds the score to kept reads) and `--tag-name pm|DS`:
- `pm` (default) writes `pm:f`, a 32-bit float;
- `DS` writes `DS:Z` with 3 decimals, as PMDtools `--writesamfield` does.

Run `pmdshodh <command> --help` for the full list.

### Exit status and streams

- **Exit status:**
  - `0`: success, including empty input;
  - `1`: runtime failure (e.g. cannot write output);
  - `2`: invalid options or invalid input, e.g. a mapped read without an MD tag (the message names the read).
- **Streams:** results go to stdout only when requested (`-o -`, or `score`/`damage` without `-o`). The run summary line and errors go to stderr, and `--quiet` silences the summary line.
- **No partial output files:** output files are written under a temporary name and renamed only on success.

## Parity with PMDtools

The scoring model is implemented from the paper and from black-box tests of PMDtools 0.60. PMDtools' code was never read (a clean-room implementation, which keeps pmdshodh's licence independent of PMDtools' GPL-3.0).

Checks so far, each comparing every read:

| Test | Reads | Result |
|---|---|---|
| Single-site score terms (ds and ss, Q 0–93, 6 parameter sets) | 5,412 | match to ≤ 1.1 × 10⁻¹³ |
| Synthetic reads with clips, indels, `=`/`X`, both strands, 8 modes | 1,500 × 8 | match to ≤ 8.3 × 10⁻¹⁴ (with `--pmdtools-compat`) |
| nf-core/eager test BAMs (mammoth mtDNA, human 1240k capture), modes ds, `--CpG`, `--UDGhalf`, `--ss` | 13,670 × 4 | match to ≤ 3.0 × 10⁻¹⁴ (with `--pmdtools-compat`) |
| Damage tables vs `pmdtools --deamination --range 30` (no UDG, full UDG) | 1,500 | match at all 30 positions |

Details and the exact methods are in [docs/blackbox.md](docs/blackbox.md).

### Where pmdshodh deliberately differs

PMDtools 0.60 has bugs that change scores. pmdshodh does the correct thing by default; `--pmdtools-compat` reproduces PMDtools exactly. The two bugs that change scores:

- **Reference rebuilt wrongly for reads that combine an insertion or soft clip with a deletion or `=` CIGAR operations.** Deleted reference bases are scored against the read bases after them, so a deleted G over a read A looks like G→A damage, and later positions shift.
- **CpG context wraps around the read ends** (in `--CpG` and `--UDGhalf` modes). On the human UDG-treated test BAM this changed the score of **6% of reads**, by up to 1.62. So, in CpG modes, expect pmdshodh's default results to differ from PMDtools'.

Other differences, which avoid crashes or silent data loss:

| Situation | PMDtools 0.60 | pmdshodh |
|---|---|---|
| Mapped read without MD tag | dropped silently | error, exit 2 |
| MD without `0` after a deletion (MALT, issue #5) | misparsed silently | read correctly |
| Hard-clipped reads | dropped | scored |
| Input already carrying a `DS:Z` tag; `DS:Z:` text in QUAL (#9) | crash | handled |
| Empty input | crash | empty output, exit 0 |
| Damage site ≥ 1,000 bp from a read end; hundreds of mismatches in one read | crash | scored |
| `--first` (first-base damage with standard error) | crashes in the 0.60 release (#8, #14) | in the JSON summary |

The full list, with the reasoning for each, is in [docs/design.md](docs/design.md) §5.

**What does not change:** the likelihood model itself, including two things that look odd but are what PMDtools computes.
- The null model assumes a background damage rate of 0.001.
- In `--ss` mode each C site contributes two terms, so ss scores are on a larger scale than ds scores.

pmdshodh keeps both so its scores mean the same as PMDtools'.

## JSON summary

`--summary` (or `damage`) writes one JSON file per run (`format_version: 1`). It contains:
- the tool version and command;
- the library mode and model parameters;
- read counts: input, scored, kept, removed, and skipped by reason;
- flag counts (duplicate, secondary, supplementary, …);
- a score histogram and read-length histograms;
- **damage tables**: counts of read bases A/C/G/T at reference C and G, by position from each end of the molecule;
- the first-position C→T rate with its standard error.

Field changes bump `format_version`. A MultiQC module is planned.

## Known limitations

- CRAM input is not supported yet.
- No benchmarks yet ([plan](docs/benchmark-plan.md)).
- No real single-stranded library has been tested; ss parity rests on synthetic reads.
- UDG-half damage *tables* do not yet match PMDtools' (scores do); see [docs/blackbox.md](docs/blackbox.md) §1.6d.
- Reads are scored whatever their SAM flags (duplicates, secondary alignments, unpaired mates of read pairs), as in PMDtools. Remove what you don't want first; the summary counts them.

## Development

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
```

The tests compare pmdshodh against values that PMDtools 0.60 produced, stored in `tests/data/`. The scripts in `tests/oracle/` regenerate them, and need PMDtools 0.60 under Linux or WSL. `tests/oracle/pmd_reference.py` is a short, readable Python version of the scoring model.

Background research is in [docs/](docs/): literature review, feature survey, black-box findings and design decisions.

## Citation

If you use pmdshodh, please cite the method it implements:

> Skoglund P, Northoff BH, Shunkov MV, Derevianko AP, Pääbo S, Krause J, Jakobsson M (2014). Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal. *PNAS* 111(6):2229–2234. [doi:10.1073/pnas.1318934111](https://doi.org/10.1073/pnas.1318934111)

## Licence

[MIT](LICENSE). pmdshodh contains no PMDtools code.
