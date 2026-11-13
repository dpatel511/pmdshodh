# Phase 4: design for pmdshodh v0.1 — draft 2026-11-11

Status: **agreed 2026-11-11** (decisions in section 10). Nothing is implemented. Builds on [literature.md](literature.md), [features.md](features.md) and [blackbox.md](blackbox.md).

## 1. Goal of v0.1

A maintained, fast replacement for `samtools calmd | pmdtools --threshold T | samtools view -b` that:
1. reproduces PMDtools 0.60 scores on every read by default (the spec is blackbox.md §1);
2. reads SAM/BAM/CRAM and writes SAM/BAM directly, multi-threaded, in input order, with bounded memory;
3. writes a versioned JSON summary for a MultiQC module;
4. fixes PMDtools' crashes and silent data loss, listing every departure with a regression test.

Out of scope for v0.1 (CLAUDE.md §7 and features.md §6): base-quality rescaling, masking bases, taxonomy, contamination estimates, genotypes, Briggs-model posteriors.

## 2. Language and libraries

**Recommendation: Rust, using noodles, with a small Python reference model as the test oracle.**

- **Why Rust.** Speed is the reason pmdshodh exists. The scoring is a few arithmetic operations per C/G site, so runtime is dominated by BAM decoding and compression, which needs native code and threads.
- **Why noodles** (MIT, pure Rust; `noodles` 0.119, `noodles-bam` 0.98, `noodles-cram` 0.102, `noodles-bgzf` 0.52, all updated October 2026):
  - builds anywhere Rust builds, including Windows, with no C toolchain;
  - simple Bioconda and wheel builds;
  - multi-threaded BGZF reading and writing.
  - **Risk:** 0.x versions change their API often. Pin exact versions and upgrade deliberately.
- **Alternative: rust-htslib 1.0.2** (bindings to C htslib). More battle-tested CRAM, but it needs a C toolchain, is hard on Windows, and pmd-mask hit a memory leak and a segfault through it (literature.md §5.2). Keep it as the fallback if noodles turns out too slow or incomplete; measure before switching.
- **Other crates:** `clap` (CLI, `-h/--help`, `--version`), `rayon` or plain worker threads (parallel scoring), `serde_json` (summary).
- **The Python oracle** is the fitted model from blackbox.md §1, about 60 lines in `docs/blackbox/harness.py`, moved to `tests/oracle/`. It is written from black-box results, not PMDtools code, and it is the version you can read and check. Rust tests compare against values it produced, and against PMDtools output on the same synthetic reads.
- **Not installed yet:** no Rust toolchain on this machine (CLAUDE.md learnings). Installing `rustup` needs your go-ahead at the start of phase 5.

## 3. Code layout

```
pmdshodh/
  Cargo.toml
  src/
    main.rs        CLI parsing, exit codes, wiring
    model.rs       scoring formulas only (no I/O): per-site terms, ds and ss, null background
    sites.rs       turn one record (SEQ, QUAL, CIGAR, MD, flag) into informative sites with z, y, Q
    modes.rs       library (ds/ss) and treatment (none/half/full) site selection
    filter.rs      threshold logic, read categories, skip reasons
    summary.rs     counts, histograms, damage tables, JSON writer
    io.rs          noodles readers and writers, ordered multi-threaded pipeline
  tests/
    oracle/        Python reference model + fixture generator
    data/          small synthetic SAM fixtures with PMDtools 0.60 scores
    cli.rs         streams, exit codes, output validity, thread determinism
```

`model.rs` and `sites.rs` hold all the science. They should stay short enough for you to review line by line against blackbox.md §1.

## 4. Command-line interface

Subcommands, per CLAUDE.md §11 (no extra executables). PMDtools flag names are accepted as aliases where the meaning is identical.

```
pmdshodh filter  [OPTIONS] <INPUT> -o <OUTPUT>   keep reads by PMD score; write SAM/BAM
pmdshodh score   [OPTIONS] <INPUT> [-o FILE]      per-read score table (TSV)
pmdshodh damage  [OPTIONS] <INPUT> [-o FILE]      damage summary only (JSON)
```

### 4.1 Shared options
| Option | PMDtools alias | Default | Notes |
|---|---|---|---|
| `--library ds\|ss` | `--ss` | `ds` | Never guessed from data (CLAUDE.md §7). Always recorded in the summary. |
| `--treatment none\|half\|full` | `--UDGminus`, `--UDGhalf`, `--UDGplus`, `--CpG` | `none` | `full` = CpG sites only; `half` = CpG plus the first and last base (blackbox.md §1.5). |
| `--pmd-p` | `--PMDpparam` | 0.3 | |
| `--pmd-c` | `--PMDconstant` | 0.01 | |
| `--pi-ancient` | `--polymorphism_ancient` | 0.001 | |
| `--pi-contamination` | `--polymorphism_contamination` | 0.001 | |
| `--min-mapq` | `-m`, `--requiremapq` | 0 | drops reads below |
| `--min-baseq` | `-q`, `--requirebaseq` | 0 | drops bases below (PMDtools behaviour) |
| `--no-clips` | `--noclips` | off | drops reads with soft clips |
| `--no-indels` | `--noindels` | off | drops reads with I or D |
| `--skip-flags <list>` | — | none (D2: score all, like PMDtools) | e.g. `secondary,supplementary,duplicate,qcfail` |
| `-r, --reference <FASTA>` | — | none | needed only to decode CRAM. **Not implemented yet** (step 4 reads SAM/BAM only) |
| `-@, --threads <N>` | — | 1 | output identical for any N |
| `--summary <FILE>` | — | none | JSON summary (section 6) |
| `-n, --max-reads <N>` | `-n` | all | `damage` only; stops cleanly, exit 0 |
| `--quiet` | — | off | suppress routine stderr; errors still printed |

The null-model background damage (0.001) is fixed and documented, not an option, because no PMDtools option changes it.

### 4.2 `filter`
- `--threshold <F>` (PMDtools `-t/--threshold`): keep reads with PMDS **≥ F**.
- `--upper-threshold <F>` (`--upperthreshold`): keep reads with PMDS **< F**. Both may be given.
- `-o <FILE>`: format from extension (`.bam`, `.sam`). `-o -` writes BAM to stdout (`--output-format sam` for SAM). Stdout then carries only the alignment stream (CLAUDE.md §11).
- `--tag`: add the score tag to kept reads, overwriting an existing tag of that name. The default is **`pm:f`** (agreed 2026-11-11), a 32-bit float, so about 7 significant digits. `--tag-name DS` writes `DS:Z` exactly as PMDtools `--writesamfield` does (Python `str(round(x, 3))`; tested against PMDtools on 483 reads) (D1).
- Output to a file is written to `<name>.pmdshodh-tmp` and renamed only when the run succeeds, so a failed run leaves no partial file. Output to stdout cannot be withdrawn: if a reader closes the pipe early, pmdshodh reports the write error and exits 1.
- The header is copied from the input, plus one `@PG` line. Kept records are byte-for-byte unchanged except for the optional tag (CLAUDE.md §7).
- `--header` is accepted and ignored (pmdshodh always writes a header), with a one-line note on stderr unless `--quiet`.

### 4.3 `score`
TSV with a header line: `read_name  flag  pmds  n_sites`. Unlike PMDtools `--printDS`, it includes the read name and flag.
- The likelihood columns were dropped during build step 3: they would need products of probabilities, which underflow on long reads (the cause of a PMDtools crash).
- Only scored reads are listed. Skipped reads are counted on stderr (and in the JSON summary from step 5).

### 4.4 `damage`
Writes only the summary JSON (section 6). It replaces `pmdtools --deamination --range N`, and `--first` (which crashes in the PMDtools release). `--range N` (default 30).

### 4.5 Exit status
- 0: success, including empty input.
- 1: runtime failure (I/O error, write failure, out of memory).
- 2: invalid options, or invalid input data (e.g. a mapped read without MD, MD that contradicts the CIGAR). The message names the read and the reason.

The numbers are fixed once released.

## 5. Read handling and departures from PMDtools

Default = PMDtools behaviour, except where marked. Every departure gets a regression test that shows the PMDtools behaviour and the new one (CLAUDE.md §7).

| Situation | PMDtools 0.60 | pmdshodh v0.1 | Departure? |
|---|---|---|---|
| Unmapped read (flag 4, no CIGAR) | dropped (no MD) | not scored, not written by `filter`, counted as `unmapped` | no (same outcome) |
| Mapped read without MD | **dropped silently**, exit 0 | **error, exit 2** (CLAUDE.md §7) | **yes** |
| MD without the `0` separator after a deletion (MALT, issue #5) | misparsed silently | parsed using the CIGAR deletion length; correct score | **yes**, fixes #5 |
| MD that contradicts the CIGAR (length or base) | scored silently | error, exit 2 | **yes** |
| Hard-clipped read (`H`) | dropped, stderr message | scored (hard clips are not in SEQ, so z is unaffected) | **yes** (D3) |
| 1-bp read | dropped silently | scored | **yes** (trivial) |
| Read with no informative sites | PMDS 0, kept at threshold ≤ 0 and by `--upperthreshold` | PMDS 0, counted as `no_informative_sites`; kept at `--threshold` ≤ 0, but **not kept by `--upper-threshold`** (no evidence of being undamaged); `--pmdtools-compat` keeps it (agreed 2026-11-13) | **yes** (`--upper-threshold` only) |
| Secondary, supplementary, duplicate, QC-fail, paired mates | all scored | all scored by default, counted in summary; `--skip-flags` to drop | no (D2) |
| Existing `DS:Z` tag in input | **crash** unless `--writesamfield` | left untouched; scoring unaffected | **yes**, fixes a crash |
| `DS:Z:` text inside QUAL (issue #9) | crash | not an issue (field-aware parsing) | **yes**, fixes #9 |
| Empty or header-only input | **crash**, exit 1 | empty output with header, zero-count summary, exit 0 | **yes** |
| `--ss` scoring | each C site counted twice (ds term + S1 term): terminal C→T 9.68 at the 5′ end, 6.34 at the 3′ end (Q40) | **one term per site (SI Eq. S1, null background 0.001)**: 4.65 at either end; `--pmdtools-compat` double-counts (agreed 2026-11-13) | **yes** |
| Reference rebuilt for reads combining `I`/`S` with `D` or `=` (blackbox.md 1.6c) | wrong: deleted bases scored against later read bases, positions shifted | **correct by default**; `--pmdtools-compat` reproduces PMDtools (agreed 2026-11-11) | **yes** |
| CpG context at the read ends | wraps around (first G of a forward read, last C of a reverse read) | no wrap by default; `--pmdtools-compat` wraps. Kept as default after real data showed it changes 6% of UDG-data scores (blackbox.md §7; agreed 2026-11-11). Must be stated prominently in user docs. | **yes** |
| Lowercase bases in SAM SEQ | not recognised | treated as uppercase (BAM cannot store lowercase) | **yes** (SAM only) |
| Lowercase reference bases in MD (outside the SAM spec; 483 of 14.0 M reads in ENA ERR1892229, docs/benchmarks) | never counted as C or G: not a site, not CpG context, not in the damage table | **read as uppercase by default**; `--pmdtools-compat` keeps them lowercase and reproduces PMDtools (found 2026-11-13) | **yes** |
| Reverse-strand read whose MD has a reference character other than `ACGTNRY`/`acgt` (e.g. lowercase `n`, most IUPAC codes; 3 reads in ERR1892229) | the character is dropped when the reference is reverse-complemented, so the rest of the reference shifts by one base per dropped character: wrong sites, scores up to +15.7 | **no shift by default**; `--pmdtools-compat` reproduces the shift (found 2026-11-13) | **yes** |
| Site ≥ 1000 bp from an end, or hundreds of mismatches in one read | **crash** (IndexError / ZeroDivisionError), exit 1 | scored (sum of log terms) | **yes** |
| SEQ `*` or QUAL `*` | dropped silently | skipped, counted as `no_sequence` / `no_qualities` | no (same outcome) |
| CIGAR `N` or `P` | dropped with a stderr message | skipped, counted as `unsupported_cigar` | no (same outcome) |
| Null model | background damage 0.001, π_c | same | no |
| Soft clips, insertions in z | counted | same | no |
| `@PG` header line | none | added | allowed by CLAUDE.md §7 |

## 6. Summary JSON (`format_version: 1`)

One file per run. No timestamps, so it is identical for identical inputs and options, whatever the thread count.

```json
{
  "format_version": 1,
  "tool": {"name": "pmdshodh", "version": "0.1.0"},
  "command": "pmdshodh filter --threshold 3 ...",
  "subcommand": "filter",
  "input": "sample.bam",
  "mode": {"library": "ds", "treatment": "none"},
  "model": {"pmd_p": 0.3, "pmd_c": 0.01, "pi_ancient": 0.001,
            "pi_contamination": 0.001, "null_background": 0.001,
            "score": "pmdtools-0.60"},
  "filters": {"threshold": 3.0, "upper_threshold": null, "min_mapq": 0,
              "min_baseq": 0, "no_clips": false, "no_indels": false,
              "skip_flags": []},
  "reads": {
    "input": 0, "scored": 0, "kept": 0, "removed_by_threshold": 0,
    "no_informative_sites": 0,
    "skipped": {"unmapped": 0, "min_mapq": 0, "clips": 0, "indels": 0, "flags": 0},
    "flags": {"secondary": 0, "supplementary": 0, "duplicate": 0,
                    "qc_fail": 0, "paired": 0, "reverse": 0}
  },
  "score_histogram": {"bin_width": 0.5, "min": -20.0, "max": 40.0,
                      "underflow": 0, "overflow": 0, "counts": []},
  "length_histogram": {"scored": {}, "kept": {}},
  "damage": {
    "range": 30, "orientation": "molecule", "reads_used": "scored",
    "five_prime":  {"C": {"A": [], "C": [], "G": [], "T": []}, "G": {"A": [], "C": [], "G": [], "T": []}},
    "three_prime": {"C": {"A": [], "C": [], "G": [], "T": []}, "G": {"A": [], "C": [], "G": [], "T": []}},
    "first_position": {"c_to_t": 0.0, "se": 0.0, "n": 0}
  }
}
```

Notes:
- **Damage counts** (as built in step 5): read bases A/C/G/T at reference C and at reference G, per position (index 0 = terminal base) from each end of the molecule.
  - Reverse reads are reverse-complemented, and `N` read bases are skipped.
  - Raw counts are stored, not rates, so MultiQC can merge and recompute.
  - ds damage is `five_prime.C` (C→T) and `three_prime.G` (G→A); ss damage is C→T in `five_prime.C` and `three_prime.C`.
  - Counts honour `--treatment` and `--min-baseq`.
  - With `--pmdtools-compat`, `five_prime.C` and `three_prime.G` reproduce PMDtools `--deamination` for no and full UDG (blackbox.md 1.6d).
  - JSON keys are sorted alphabetically in the file.
- **Standard error:** `first_position.se` is the binomial standard error √(p(1−p)/n). This is the documented meaning of `--first` ("with a standard error"). PMDtools' own `--first` can't be compared, because it crashes.
- **Histogram:** fixed bins make runs comparable in MultiQC. Exact edges are open to change before release; after release, any field change bumps `format_version`.

## 7. Performance design (no claims yet)

- **Pipeline:** a reader thread decodes BGZF blocks (noodles multi-threaded reader). Batches of records go to N scoring workers. A writer re-orders batches by sequence number and writes with a multi-threaded BGZF writer.
- **Memory:** bounded by (batch size × in-flight batches), independent of file size.
- **Determinism:** the score table and summary JSON are byte-identical for any thread count, and BAM output has identical records and header. The one exception is `@PG CL`, which records the command line, including `-@`. Because the header shares a BGZF block with the first records, compressed bytes then differ from that point. A test checks all of this on 25,500 reads at `-@ 1` and `-@ 4` (CLAUDE.md §12).
- **As built (step 6):**
  - records are read in batches of 10,000;
  - each batch is scored in parallel chunks on a rayon pool of `-@` threads, with per-chunk damage counts summed afterwards;
  - results are applied in input order;
  - with `-@` above 1, the stages overlap (2026-11-13): a reader thread reads batch k+2 and a scoring thread scores batch k+1 while the main thread applies batch k, so at most three batches are in flight. With `-@ 1` the stages run one after another, as before;
  - BGZF input is decompressed, and BAM output compressed, by noodles' multithreaded BGZF reader and writer. BAM always uses the multithreaded writer, so the block layout does not depend on `-@`.
  - Measured speed: docs/benchmarks.md.
- **Benchmark protocol for later speed claims:**
  - same machine, same input;
  - PMDtools 0.60 (Bioconda, with its `samtools calmd`/`view` pipe as in nf-core), metaDMG-cpp `pmd`, pmdshodh at 1, 4 and 8 threads;
  - report hardware, versions, exact commands, wall time and peak memory over at least 3 runs.

## 8. Testing plan

1. **Model unit tests** (`model.rs`): per-site values from blackbox.md probes (z = 1…30, Q = 0…93, ds and ss, parameters changed), checked against PMDtools outputs to ≤ 1e-12.
2. **Site extraction tests** (`sites.rs`): soft clips, insertions, deletions, hard clips, reverse strand in ss, CpG at read edges, MD edge cases (MALT style, `n`/IUPAC, contradictions).
3. **Synthetic parity fixture:** the 500 random reads from blackbox probe3 and probe8 (and new ones with indels and clips), with PMDtools 0.60 scores committed as expected output. These are small text files.
4. **Departure regression tests:** one test per "yes" row in section 5, each showing the PMDtools behaviour (recorded) and pmdshodh's.
5. **CLI tests:** stdout/stderr separation, exit codes 0/1/2, `-h`, `--version`, empty input, `-o -`, BAM validity (`samtools quickcheck` on Linux CI), threads 1 vs 8 byte-identical, input order preserved.
6. **Real-data parity (before any "parity" claim):** compare every read's score with PMDtools 0.60 on at least one public ds and one ss dataset. Report read count, largest difference and tolerance (CLAUDE.md §7).
   - Candidates: nf-core/test-datasets eager BAMs (small; licence to check) and the Bacho Kiro ss library from PMDtools PR #16 (size to check).
   - **No download before you approve the size** (CLAUDE.md §7).

## 9. Delivery order (phase 5 onward)

Each step ends with passing tests before the next starts.
1. Install Rust (your approval); `cargo new`; CI-free local checks (`fmt`, `clippy`, `test`).
2. `model.rs` + unit tests against blackbox values.
3. `sites.rs` + SAM text input only + `score` subcommand → synthetic parity fixture passes.
4. BAM input/output, `filter`, `@PG`, tag, exit codes.
5. Summary JSON + `damage`.
6. Threads + determinism test.
7. Real-data parity run (your approval for data).
8. Benchmarks; then README, licence, packaging (each needs your authorisation).

## 10. Decisions

**Agreed 2026-11-11:**
- **D1: both tags.** By default a local-namespace float tag; `--tag-name DS` writes PMDtools-style `DS:Z` (3 decimals) for compatibility.
  - Local tag name agreed: **`pm:f`** (lowercase, so in the SAM spec's local namespace). Note: `f` is a 32-bit float.
- **D2: score flagged reads like PMDtools**, report per-flag counts, and offer `--skip-flags` to drop them.
- **D3: score hard-clipped reads.** Documented departure with a regression test.
- **D4: Rust + noodles**, with the Python oracle. Pin noodles versions; rust-htslib is the measured fallback.
- **D5: licence** deferred until before the first public commit.

The options considered are kept below for the record.

**[D1] Score tag.**
- Options: `DS:Z` (PMDtools-compatible text, 3 decimals); a local-namespace float tag (e.g. `XP:f`, spec-safe and full precision); or both (local by default, `--tag-name DS` for compatibility).
- Recommendation: **both**.

**[D2] Flagged reads** (secondary, supplementary, duplicate, QC-fail).
- Option 1: score them like PMDtools (parity by default), with `--skip-flags` to drop them.
- Option 2: skip them by default.
- Recommendation: **parity by default, with counts in the summary**, because eager deduplicates beforehand and a default that drops reads changes results for drop-in users.

**[D3] Hard-clipped reads.** Score them (departure; avoids silent loss) or drop like PMDtools. Recommendation: **score them**.

**[D4] Language.** Rust + noodles, with Python oracle (section 2). Recommendation as written. The alternative is Python + pysam: easier for you to own, speed unknown.

**[D5] Licence.** Not needed until publishing; the clean room keeps all options open (MIT/Apache-2.0 is the Rust norm; GPL-3.0 also possible). Recommendation: decide before the first public commit, not now.

**Not design decisions, still yours:**
- whether to confirm and report the eager 2.5.3 `--CpG` finding;
- whether to contact metaDMG.
