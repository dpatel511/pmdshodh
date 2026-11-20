# Audit plan: 0.3.0, before milestone 6

Status: **plan, 2026-11-20.** Nothing here has run yet. The decisions in section 2 were agreed on 2026-11-20. Findings will be recorded in section 10 as they come in.

## 1. Goal, scope and when it is done

Milestone 6 (fitted damage parameters) will change the scoring model. Before that, this audit checks 0.3.0 broadly, in four ways:
- What the tests do not cover (section 3).
- Edge cases, read by read and option by option (section 4).
- Guarantees stated as properties and checked on generated inputs (section 5).
- Every metric on a wide set of real libraries (section 7).

In scope: `filter`, `score`, `damage`, `report`, the JSON summary, input decoding (SAM, BAM, CRAM, `-r`) and the outputs.

Out of scope:
- New features.
- Changing the model or defaults. Any finding there is recorded and decided separately (CLAUDE.md §7).
- Comparison with mapDamage2 or DamageProfiler (decision 2.1).

**Done when** all of the following hold:
- Every uncovered line or branch from section 3 has a test, or a written reason why it needs none.
- Every cell of the matrix in section 4 has a test or is marked not applicable.
- The properties in section 5 hold over at least 10,000 generated cases each.
- The fuzz targets in section 6 have run for 1 hour each with no crash.
- Every check in section 7 passes, or each difference is explained in section 10.

Findings are documented first. Fixes follow in separate PRs, each with a regression test, as 0.3.x patch releases (section 11).

## 2. Decisions (agreed 2026-11-20)

1. **Metrics are compared against PMDtools 0.60 only:** scores, kept reads and `--deamination` tables in `--pmdtools-compat`. Plus internal checks: report against JSON, counts against samtools, threads, `-r` against MD. mapDamage2 and DamageProfiler count differently (all reads, not the scored ones, and in their own orientation), so a comparison would explain differences rather than test correctness. Not now.
2. **`proptest` is added as a dev-dependency.** It is used by tests only and not shipped in the binary or wheels.
3. **Fuzzing covers two parsers:**
   - alignment input (SAM text, BAM records, the MD/CIGAR walk);
   - summary JSON (`report::parse`).

   It uses `cargo-fuzz` on a nightly toolchain, in a `fuzz/` directory outside the published crate and the wheels. It does not change the stable build or CI.
4. **A UDG-treated single-stranded library will be searched for** (section 8). I'll state its source, licence and size before any download.
5. **Local tools:** `cargo-llvm-cov` for coverage and `cargo-fuzz`, installed on the bench machine (WSL) only. They are not project dependencies.

## 3. Part A: coverage map

- Run `cargo llvm-cov --all-targets --html` under WSL on `main`. Record line and branch coverage per file.
- For each uncovered region, record one of:
  - (a) missing test → add it in section 4;
  - (b) unreachable or defensive → note why;
  - (c) platform-only (Windows or Unix branches) → note which CI job covers it.
- Output: a table in section 10 per source file (`main.rs`, `sites.rs`, `model.rs`, `reads.rs`, `damage.rs`, `summary.rs`, `report.rs`, `cram.rs`, `bgzf.rs`, `reference.rs`), with numbers before and after.

## 4. Part B: edge-case matrix

Each row is checked by an existing test (named) or a new one. Where PMDtools can read the input, `--pmdtools-compat` is compared with PMDtools for that case.

### 4.1 Reads and alignments
| Area | Cases |
|---|---|
| Layout | single-end; merged; paired (both mates, one mate unmapped, mate on another reference); forward and reverse; secondary, supplementary, duplicate, QC-fail, unmapped-but-placed |
| Length | 1 bp; shorter than `--range`; exactly `--range` and 2×`--range`; 1,000 bp; 66,000 bp; reads whose two damage windows overlap |
| CIGAR | each op (M, I, D, N, S, H, P, =, X); leading and trailing clips; I or D next to a clip; N introns; `*` CIGAR; more than 65,535 ops (BAM CG tag) |
| SEQ, QUAL | `*` SEQ; `*` QUAL (0xFF); lowercase; N; IUPAC codes; `=` bases with and without `-r`; base qualities at the `--min-baseq` boundary |
| MD | absent; empty; only matches; deletion (`^`); lowercase reference bases; MD that disagrees with CIGAR length; malformed; MD on an unmapped read |
| Position | read at reference start and end; reference of length 1; reads on a reference missing from the header or the FASTA; soft-masked and N reference bases with `-r` |

### 4.2 Scoring and filtering
| Area | Cases |
|---|---|
| Model | ds, ss × UDG none, half, full × default and `--pmdtools-compat`; CpG at read ends and spanning a deletion; reads with 0 informative sites |
| Parameters | boundary values of `--pmd-p`, `--pmd-c`, `--pi-*` (accepted and rejected); p + c at the limit of 1; `--min-mapq`, `--min-baseq` at 0 and 255 |
| Thresholds | score exactly at `--threshold` (kept) and at `--upper-threshold` (dropped); both given; negative thresholds; very large values; NaN and inf rejected |
| Output | record order kept; bytes unchanged except the tag; tag overwrite; `--tag-name DS` formatting; SAM vs BAM; output to stdout; header-only input; `-@ 1` against `-@ 8` |

### 4.3 Summaries and report
| Area | Cases |
|---|---|
| `damage` and `--summary` | `--range` 1 and 1000; `-n` limits (0, 1, more than the input); empty input; all reads skipped; every skip reason counted once |
| JSON | key order stable; integer counts above 2^53 (precision in consumers); `format_version` 1 and 2 |
| `report` | 1 and 50 samples; Unicode and very long file names; labels that still clash after the summary name is added; `range` 1000; histograms with all counts in underflow or overflow; summaries made with `-r`, CRAM, ss, half, full and compat |

### 4.4 Command line and outputs
- Every option with an invalid value, a missing value or a value repeated.
- Unknown options.
- Help and version output, stream use and exit codes.
- Already covered, rechecked only: output paths (symlink, FIFO, `/dev/stdout` as non-root, `>(...)`, devices), signals, closed stdout.

## 5. Part C: properties (`proptest`)

Generated reads, built from random reference windows with random damage, mismatches, indels, clips, strands and qualities. Each property runs on at least 10,000 cases:

1. **Threads:** `filter`, `score` and `damage` output is byte-identical at `-@ 1`, 2 and 8.
2. **Filter:** kept records are a subsequence of the input, in input order, and unchanged except for the optional tag.
3. **Thresholds:** a read is kept if and only if its score from `score` meets the threshold.
4. **`-r` against MD:** scores from the reference equal scores from correct MD tags.
5. **SAM against BAM:** the same records give the same scores and summary.
6. **Determinism:** the same summaries give the same report bytes, and `report::parse` accepts every summary pmdshodh writes.
7. **Counts:** input = scored + skipped by every reason, and kept + removed = scored.
8. **Model:** compat scores equal the Python oracle (`tests/oracle/pmd_reference.py`) on generated reads within 1e-9.

## 6. Part D: fuzzing (`cargo-fuzz`, nightly, outside the crate)

- **Target 1, alignment input.** Arbitrary bytes as SAM, and as a BAM body behind a valid header, through to scoring. Expect a clean error (exit 2) or a result. A panic, hang or memory growth is a finding.
- **Target 2, summary JSON.** Arbitrary JSON, and valid summaries with mutated fields, through `report::parse` and `render`. Expect an error naming the field, or valid HTML with no NaN or inf.
- Run each for 1 hour. Minimise any crash into a regression test in `tests/`.

## 7. Part E: real-data matrix (data already on the bench machine)

| Library | Type | Reads | Source |
|---|---|---|---|
| BER11_L_NOuser | ancient horse, ds, no UDG | 14.0 M | ENA ERR1892229 (docs/benchmarks.md) |
| BER11_L (all) | adds the USER-treated libraries, whose MD tags do not match their CIGARs | 3.9 GB file (count to be recorded) | ENA ERR1892229 |
| Motala3 | Mesolithic human, ds, likely UDG-treated (unconfirmed) | 23.6 M | ENA ERR577494 |
| UKY001.B0102 1240K and MT | ss, no UDG | 14.5 M and 0.28 M | ENA PRJEB37007 |
| JK2782, JK2067 | eager test BAMs: mammoth mtDNA; human 1240k, UDG-treated | 9,892 and 12,473 | nf-core/test-datasets |
| Aligner sweep | bwa aln and mem, bowtie2, samtools outputs, CRAM 3.0 and 3.1, 14 crafted reads | up to 0.2 M each | ERR3991018 on rCRS (CHANGELOG 0.2.3) |

For each library, in each applicable mode (ds or ss; UDG none, half or full; default and compat; BAM and CRAM; with and without `-r`):

| Check | Against | Pass when |
|---|---|---|
| Per-read scores (compat) | PMDtools 0.60 `--printDS` | every read compared; largest difference ≤ 1e-6 |
| Kept reads (compat) | PMDtools `--threshold` | identical read lists |
| Damage tables (compat) | PMDtools `--deamination` | equal to its 5 decimals |
| Default mode against compat | design.md §5 departures | every differing read is explained by a listed departure |
| Read counts | `samtools view -c` and `flagstat` | input, flags and skipped counts agree |
| Threads | `-@ 1` against `-@ 8` | byte-identical outputs |
| `-r` | the same file after `samtools calmd` | identical scores |
| Report | its JSON, using the number-check script (docs/report-proposal.md §9), which this audit commits under `docs/benchmarks/` | 0 FAIL |
| Outputs | `samtools quickcheck`, and `samtools view` record-by-record | valid; kept records byte-identical to the input |
| Speed and memory | docs/benchmarks.md 0.2.1 rows, same commands, median of ≥ 3 now | no regression above run-to-run spread, or explained |

## 8. Additional test data

The libraries in section 7 cover ds without UDG, likely-UDG ds, ss without UDG and a range of aligner outputs. Four kinds are missing. I'll look for public sources and list each with source, licence, size and date before any download:

1. **A UDG-treated single-stranded library** (decision 2.4). It is the only library type never compared with PMDtools.
2. **A double-stranded library with UDG treatment confirmed by its publication** (half and full, if both exist). Motala3 and JK2067 are only presumed treated, so their UDG-mode results cannot be checked against a known truth.
3. **A metagenomic or sedimentary-DNA BAM with many references:** thousands of contigs, short references and reads near contig ends. This tests multi-reference CRAM slices, reference loading and memory with many small sequences.
4. **An nf-core/eager output BAM after mapDamage rescaling:** base qualities lowered at damaged positions. PMDtools users run on such files, so `--min-baseq` and scoring should be checked on them.

Small files (under 100 MB) are preferred. A larger file is only taken if a small one does not exist, and only after you approve its size.

## 9. Order

1. Coverage map (section 3). It sets the test list.
2. Edge-case tests (section 4), as small focused tests. A suspected bug is confirmed by a failing test before it is recorded.
3. Properties (section 5).
4. Real-data matrix (section 7) on existing data, run in the background while sections 2 and 3 proceed.
5. Fuzzing (section 6).
6. New data (section 8), after approval.

Tests that only add coverage, with no behaviour change, go in as `test/` PRs during the audit. Findings go into section 10 as they come in.

## 10. Findings

Recorded as they are found:

| # | Kind | What | Evidence (test, data, command) | Status |
|---|---|---|---|---|

Kinds:
- **B:** wrong output, crash, hang, wrong exit code, or a broken guarantee.
- **D:** documentation does not match behaviour.
- **O:** an observation that needs a decision.

## 11. Fixes

- Each B or D finding gets one `fix/` PR from `main`, with a regression test that fails before the fix. The test is listed by name, and its PR links the finding.
- Scoring changes need the PMDtools parity comparison rerun and reported (CLAUDE.md §12).
- Fixes are released as 0.3.x patch releases, with notes, the same way as 0.2.x.
- O findings wait for your decision.
- Milestone 6 starts after the done-when conditions in section 1 hold.
