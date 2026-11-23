# Audit plan: 0.3.0, before milestone 6

Status: **run 2026-11-20 to 2026-11-22.** The decisions in section 2 were agreed on 2026-11-20. Findings B1–B7 (section 10) are fixed in 0.3.1, and B8 and B9 in 0.3.2. The section 4 map (4.5) found B8 and O1–O3. **Done 2026-11-23:** every done-when condition in section 1 holds (§3.1, §4.5, §5.1, §6.1, §7.1). O1 and O4 were decided on 2026-11-23 (design.md §5; a CRAM fuzz target, §6.1). B10, found on 2026-11-23 by that target, is fixed in 0.3.3.

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
- Output: a table per source file (`main.rs`, `sites.rs`, `model.rs`, `reads.rs`, `damage.rs`, `summary.rs`, `report.rs`, `cram.rs`, `bgzf.rs`, `reference.rs`), with numbers before and after.

### 3.1 Results on `main` at `e136268` (2026-11-20)

The run was `cargo llvm-cov --all-targets` (0.9.2, stable 1.99.1, WSL Ubuntu 24.04) over all 106 tests. The counts are lines and regions; branch coverage needs nightly and was not measured. The child processes in the signal tests end by signal, so they write no coverage data.

| File | Lines | Missed lines | Regions | Missed regions |
|---|---|---|---|---|
| bgzf.rs | 93.4% | 17 | 93.1% | 32 |
| cram.rs | 94.9% | 9 | 90.9% | 33 |
| damage.rs | 99.3% | 1 | 99.6% | 1 |
| main.rs | 94.8% | 59 | 92.0% | 147 |
| model.rs | 100% | 0 | 100% | 0 |
| reads.rs | 85.4% | 24 | 86.6% | 34 |
| reference.rs | 93.6% | 23 | 94.1% | 42 |
| report.rs | 98.4% | 19 | 95.1% | 113 |
| sites.rs | 98.6% | 8 | 98.3% | 16 |
| summary.rs | 99.5% | 1 | 99.6% | 1 |
| **Total** | **96.4%** | **161** | **94.6%** | **419** |

**(a) Missing tests: these become the test list.**
- **Read filters and skip reasons** (`reads.rs` 151–195). No test uses `--min-mapq`, `--no-clips` or `--no-indels`, or a record with an `N` or `P` CIGAR, `*` SEQ or `*` QUAL. No test has a QUAL length that differs from SEQ, or an MD tag that is not a string. Only the unmapped skip reason is tested. The aligner sweep (CHANGELOG 0.2.3) exercised some of these by hand, but the suite does not.
- **Full UDG from the command line** (`main.rs` 369). `--treatment full`, `--UDGplus` and `--CpG` are never run as a command; full UDG is tested only through library unit tests.
- **`=` and `X` CIGAR ops** (`sites.rs` 123–124), where they extend the reference span.
- **Malformed MD tags** (`sites.rs` 351, 364, 368, 373): an invalid character, and deletions that disagree between CIGAR and MD.
- **Reference checks:**
  - a header sequence missing from the FASTA (`main.rs` 962);
  - an empty FASTA sequence (`reference.rs` 179);
  - FASTA line lengths that differ from the `.fai` index (`reference.rs` 196–202);
  - more references than the block cache holds, so blocks are evicted (`reference.rs` 370);
  - a reference sequence reused across CRAM containers (`reference.rs` 292);
  - a CRAM alignment past its reference end (`cram.rs` 164).
- **Damaged BGZF:**
  - a block size below the header (`bgzf.rs` 177);
  - an ISIZE over 64 KiB (195);
  - bad deflate data (206);
  - a BGZF file holding only the EOF block (`main.rs` 519);
  - output of incompressible or very large records (`bgzf.rs` 278, 299).
- **Command line:**
  - `--threads 0` (`main.rs` 1068);
  - an input error in `damage` (1245);
  - the stderr line listing skip reasons (1328);
  - an output path that is a directory (1361);
  - a write failure to a file (1399; `/dev/full` as non-root on Linux).
- **Report validation** (`report.rs` 176, 195, 223, 244, 248):
  - a non-numeric read-length key;
  - `tool.name` other than pmdshodh;
  - `damage.range` 0;
  - a bin width of 0;
  - empty histogram counts.

**(b) Defensive, not reachable from the command line (no test needed; reasons recorded):**
- `sites.rs` 110 and 154: unsupported CIGAR ops, which `reads.rs` already skips.
- `damage.rs` 125: `--range` is at least 1.
- `main.rs` 1020: FLAG is always 16 bits after decoding.
- `report.rs` 398 and 591: the y range is never empty, and an empty length chart returns earlier.
- `bgzf.rs` 53, 60 and 244: worker channels close only if a worker thread panics.
- `bgzf.rs` 164–165 and 184: interrupted and non-EOF read errors.
- `reference.rs` 113 and 295–299: errors that noodles reports earlier.
- `main.rs` 1409–1416: a failed rename after a successful write.

**After the group (a) tests (PR #50):** line coverage on WSL is 97.3%, with 122 lines missed (was 161). These items stay open:
- Reference-cache eviction (more than 64 MiB of reference) and a reference reused across CRAM containers. Both need large inputs, so they move to the real-data matrix (section 7): Motala3 against human_g1k_v37, with `-r` and as CRAM. *Done 2026-11-22:* with the B9 fix, the samtools CRAMs of BER11_L_NOuser (against equCab2) and Motala3 (against human_g1k_v37) give `score -r -@ 8` output byte-identical to their BAMs, in default and compat mode. Both span many containers per chromosome and more than 64 MiB of reference. This repeats the 2026-11-16 check (design.md §11) on current code.
- A CRAM alignment past its reference end: noodles does not write one. It moves to fuzzing (section 6). *Corrected 2026-11-22:* none of the fuzz targets decodes CRAM, and samtools does write such files. Decoding one made noodles panic before this check ran, which is finding B9. The check now runs first, and `cli::cram_read_past_its_reference_end_is_an_input_error` covers it, with a samtools-written fixture.
- The FASTA line-length check (`reference.rs` 196–202): see B1.
- An empty FASTA sequence (`reference.rs` 179) moves to group (b): `samtools faidx` leaves empty sequences out of the index, and SAM requires `LN` ≥ 1.

**(c) Covered on another platform or by unmeasured processes:**
- `main.rs` 1486–1494: Windows device names (the Windows test runs on Windows only).
- `main.rs` 1449–1465: the signal handler, which the Unix signal tests exercise in child processes that end by signal.

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

### 4.5 Map: the test for each cell (2026-11-22, `main` at `c495119` plus PR #63)

Test names are given as `file::name`, where `edge` is `tests/edge_cases.rs`, `cli` is `tests/cli.rs`, `props` is `tests/properties.rs` and `parity` is `tests/read_parity.rs`. Unit tests are named by their module (`sites`, `damage`, `report`, `summary`, `reference`). **New** marks tests added by this map. Where a cell has a PMDtools behaviour, the cited test compares with PMDtools 0.60 output or probes (2026-11-22).

Probing PMDtools for the new tests found three behaviours that compat did not match. **B8**: PMDtools drops 1-bp reads; it gets its own fix PR. **O2**: lowercase SEQ, now a documented limitation. **O3**: CIGAR `*`, now a documented departure. They are in section 10, with **O1** (CpG context across a deletion).

**4.1 Reads and alignments**

| Cell | Test |
|---|---|
| Single-end, forward and reverse | `parity::read_scores_match_pmdtools_0_60_in_all_modes` (1,500 reads, both strands) |
| Merged reads | Not applicable: a merged read is a single-end record (flag 0), the same as above |
| Paired: both mates, one mate unmapped, mate on another reference | **New** `edge::every_read_layout_is_scored_counted_and_kept_unchanged`; `cli::cram_multi_reference_slice_pairs_and_wrong_reference` (a proper pair in CRAM) |
| Secondary, supplementary, duplicate, QC-fail | **New** `edge::every_read_layout_is_scored_counted_and_kept_unchanged` (scored and counted by flag) |
| Unmapped but placed | `edge::compat_scores_placed_unmapped_reads_like_pmdtools` (B4) |
| 1 bp | **New** `edge::one_base_and_long_reads_are_scored_by_default` (default); `edge::compat_drops_one_base_reads_like_pmdtools` (compat, B8) |
| Shorter than `--range`, exactly `--range`, 2 × `--range`, overlapping damage windows | `damage::counting_only_the_ends_matches_counting_every_position` (lengths 1–150 with range 10) |
| 1,000 bp | **New** `edge::one_base_and_long_reads_are_scored_by_default` (a site at z = 1,200) |
| 66,000 bp | `edge::incompressible_record_larger_than_a_block_round_trips` (70,000 bp) |
| CIGAR M, I, D, S | `parity::read_scores_match_pmdtools_0_60_in_all_modes` (fixture with clips and indels) |
| CIGAR N, P | `edge::every_skip_reason_is_counted_once` |
| CIGAR H | `edge::compat_drops_hard_clipped_reads_like_pmdtools` (B6) |
| CIGAR =, X | `edge::sequence_match_and_mismatch_cigar_ops_score_like_m` |
| Leading and trailing clips | `parity::read_scores_match_pmdtools_0_60_in_all_modes`; `sites::clips_and_insertions_have_no_reference_base` |
| I or D next to a clip | **New** `edge::indels_next_to_a_soft_clip`; `sites::default_reconstruction_fixes_clip_plus_deletion` |
| CIGAR `*` on a mapped read | **New** `edge::star_cigar_and_empty_md_on_mapped_reads_are_input_errors` (O3) |
| More than 65,535 operations | `cli::long_cigar_records_keep_a_single_cg_tag` |
| SEQ `*`, QUAL `*` | `edge::every_skip_reason_is_counted_once` |
| Lowercase SEQ | **New** `edge::lowercase_iupac_and_n_read_bases` (O2) |
| N and IUPAC read bases | **New** `edge::lowercase_iupac_and_n_read_bases` |
| `=` read bases, with and without `-r` | `cli::equals_read_bases_are_resolved_from_the_reference` |
| Base quality at the `-q` boundary | **New** `edge::quality_and_parameter_boundaries` |
| MD absent | `cli::invalid_input_exits_2_and_leaves_no_output_file`; `edge::command_line_errors` |
| MD empty | **New** `edge::star_cigar_and_empty_md_on_mapped_reads_are_input_errors` |
| MD only matches; deletion (`^`) | `sites::mismatches_and_matches_from_md`; `sites::deletion_with_and_without_zero_separator` |
| Lowercase MD | `cli::lowercase_md_bases_default_reads_uppercase_compat_matches_pmdtools`; `cli::reverse_read_lowercase_n_shift_only_in_pmdtools_compat` |
| MD that disagrees with the CIGAR; malformed MD | `sites::md_that_contradicts_cigar_is_an_error`; `edge::malformed_md_tags_are_input_errors` |
| MD on an unmapped read | `edge::compat_scores_placed_unmapped_reads_like_pmdtools` |
| Read at the reference start and end; reference of length 1 | **New** `edge::reads_at_reference_ends_score_the_same_from_md_and_reference` |
| Reference missing from the header; missing from the FASTA | `cli::kept_reads_naming_references_missing_from_the_header_are_input_errors`; `edge::reference_problems_are_input_errors` |
| Soft-masked and N reference bases with `-r` | `cli::cram_input_matches_its_sam_and_needs_the_reference` (soft-masked); `props::reference_scores_equal_md_scores` (1% N); `edge::compat_drops_reads_without_reference_bases_like_pmdtools` |

**4.2 Scoring and filtering**

| Cell | Test |
|---|---|
| ds and ss × UDG none, half and full, compat | `parity::read_scores_match_pmdtools_0_60_in_all_modes` |
| The same, default mode | `parity::default_mode_matches_pmdtools_where_its_reference_is_correct`; `edge::full_udg_options_are_equivalent_and_count_only_cpg`; `sites::single_stranded_uses_molecule_orientation` |
| CpG at read ends | `sites::cpg_context_wraps_only_in_compat_mode`; `sites::cpg_modes_select_sites_like_pmdtools` |
| CpG spanning a deletion | **New** `edge::cpg_context_skips_deleted_reference_bases` (O1) |
| 0 informative sites | `sites::no_informative_sites_scores_zero`; `cli::upper_threshold_drops_reads_without_evidence_except_in_compat` |
| `--pmd-p`, `--pmd-c`, `--pi-*` at and beyond their bounds; p + c = 1 | `cli::parameters_that_make_scores_undefined_are_rejected`; `cli::usage_errors_exit_2`; **New** `edge::quality_and_parameter_boundaries` |
| `-m` and `-q` at 0 and 255 | **New** `edge::quality_and_parameter_boundaries`; `edge::every_skip_reason_is_counted_once` (MAPQ 255) |
| Exactly at `--threshold` and at `--upper-threshold` | **New** `edge::reads_exactly_at_a_threshold` |
| Both given; negative thresholds | `props::filter_keeps_by_score_in_order_and_counts_add_up` |
| Very large values; NaN and inf | `cli::compat_filter_ignores_thresholds_outside_pmdtools_range` (B7); `cli::parameters_that_make_scores_undefined_are_rejected` |
| Record order; bytes unchanged except the tag | `props::filter_keeps_by_score_in_order_and_counts_add_up`; `cli::kept_records_are_unchanged_apart_from_the_score_tag` |
| Existing tag overwritten | **New** `edge::existing_score_tags_are_kept_or_replaced` |
| `--tag-name DS` formatting | `cli::ds_tag_text_matches_python_round_formatting` |
| SAM against BAM; stdout | `cli::bam_and_sam_input_give_identical_scores_and_damage`; `cli::output_format_from_extension_in_any_case`; `cli::stdout_carries_only_the_bam_and_the_summary_goes_to_stderr`; `cli::score_dash_output_writes_to_stdout` |
| Header-only input | `cli::empty_input_gives_an_empty_valid_bam`; **New** `edge::damage_limits_and_empty_runs` |
| `-@ 1` against `-@ 8` | `cli::output_is_identical_for_any_thread_count`; `props::output_does_not_depend_on_threads_or_input_format` |

**4.3 Summaries and report**

| Cell | Test |
|---|---|
| `--range` 1 and 1000 | **New** `edge::damage_limits_and_empty_runs`; `cli::range_and_threads_have_upper_limits`; `cli::range_must_be_at_least_1` |
| `-n` 0, 1 and more than the input | **New** `edge::damage_limits_and_empty_runs`; `cli::damage_max_reads_stops_early_with_status_0` |
| Empty input; all reads skipped | **New** `edge::damage_limits_and_empty_runs` (both also through `report`) |
| Every skip reason counted once | `edge::every_skip_reason_is_counted_once` |
| Key order stable | `cli::output_is_identical_for_any_thread_count` (summary compared); `props::filter_keeps_by_score_in_order_and_counts_add_up` (same report twice) |
| Counts above 2^53 | Not applicable: this needs more than 9 × 10^15 reads. Counts are written as JSON integers from `u64`; a consumer that reads them as doubles loses precision only past 2^53 |
| `format_version` 1, 2 and 3 | `report::version_1_lacks_input_format_and_reference`; the `report` unit tests (versions 2 and 3); `cli::summary_json_counts_and_versions` |
| `report` with 1 and 50 samples | `cli::report_from_a_damage_summary_plots_its_rates`; **New** `edge::report_on_many_and_unusual_summaries` |
| Unicode and very long file names | **New** `edge::report_on_many_and_unusual_summaries` |
| Labels that still clash | `report::labels_are_input_file_names_made_unique` |
| `range` 1000 | **New** `edge::report_on_many_and_unusual_summaries` |
| All counts in underflow or overflow | **New** `report::scores_all_outside_the_histogram_range_render`; `summary::histogram_bins_and_overflow` |
| Summaries made with `-r`, CRAM, ss, half, full and compat | `cli::report_shows_md_tag_agreement_with_reference` (`-r`); `cli::report_on_several_summaries_overlays_them_and_is_reproducible` (ss); `report::udg_notes_state_which_sites_are_counted` (half, full); **New** `edge::report_on_many_and_unusual_summaries` (CRAM, compat half) |

**4.4 Command line and outputs**

| Cell | Test |
|---|---|
| Invalid values | `cli::usage_errors_exit_2`; `cli::parameters_that_make_scores_undefined_are_rejected`; `cli::range_and_threads_have_upper_limits`; `edge::command_line_errors`; **New** `edge::quality_and_parameter_boundaries` |
| Missing and repeated values; unknown options | **New** `edge::unknown_repeated_and_valueless_options_are_usage_errors` |
| Help and version; stream use; exit codes | `cli::help_and_version_go_to_stdout_with_status_0`; `cli::subcommand_help_groups_options_and_shows_examples`; `cli::help_examples_run`; `cli::unwritable_output_exits_1`; `edge::write_failure_exits_1` |
| Symlink, FIFO and device outputs | `cli::symlink_output_replaces_its_target_and_keeps_the_link`; `cli::fifo_output_is_written_in_place`; `cli::windows_device_names_are_written_in_place`; `cli::report_fifo_output_and_signal_cleanup` |
| `/dev/stdout` | **New** `edge::output_to_dev_stdout` (Unix). A `>(...)` output is a FIFO made by the shell, so the FIFO test covers it |
| Signals; closed stdout | `cli::a_signal_mid_run_removes_temporary_files`; `cli::closed_stdout_ends_the_run_quietly_with_status_141` |

With these tests, and the B8 fix (PR #64), every cell has a test or is marked not applicable.

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

### 5.1 Results (2026-11-21, `tests/properties.rs`)

- **The generator:** a random 5,000-base reference with 1% `N`, and reads with:
  - M runs with insertions, deletions and soft clips;
  - C→T and G→A damage near the ends, and 1% mismatches;
  - both strands, random qualities and MAPQ;
  - duplicate, secondary and paired flags, and 2% unmapped reads;
  - MD tags identical to `samtools calmd`'s (checked on 592 reads).

  A failing case is reproduced from its seed with the ignored `regenerate` test.
- **Properties 1 and 5** (threads; SAM against BAM): one test, on inputs of up to 25,000 reads, so up to three batches.
- **Properties 2, 3, 6 and 7** (filter subset and order; threshold rule; report; counts): one test.
- **Property 4** (`-r` against MD; default, compat and ss): one test.
- **Property 8** is covered elsewhere, not by proptest: it needs Python, which the test suite does not depend on. It is covered by `model_parity` (precomputed oracle terms) and by the PMDtools comparisons.
- **Results:** 24 cases per property passed on WSL with the B2 fix: about 300,000 reads for properties 1 and 5, and about 60,000 for the others. CI runs 12 cases per property.
- **10,000 small cases per property (2026-11-22, for the done-when condition in section 1).** `PROP_MAX_N=50` caps each case at 50 reads.
  - The three property tests ran 10,000 cases each on WSL, from `main` at `3fc5362` (0.3.1 plus B8), as a release build.
  - They ran in five chunks of 2,000 cases with `PROPTEST_RNG_SEED` 101 to 105, so each chunk can be rerun. Command: `PROP_MAX_N=50 PROPTEST_CASES=2000 PROPTEST_RNG_SEED=101 cargo test --release --test properties`.
  - All passed, taking about 2 minutes per chunk.
  - With at most 50 reads, a case fits in one batch, so properties 1 and 5 across several batches rest on the 24 large cases above.
  - Property 8 is not a proptest property (see above).
- **Findings:**
  - property 4 found B2;
  - the first failure of property 1 was a test bug: the recorded command line contains `-@ N`, so outputs are now compared without it.

## 6. Part D: fuzzing (`cargo-fuzz`, nightly, outside the crate)

- **Target 1, alignment input.** Arbitrary bytes as SAM, and as a BAM body behind a valid header, through to scoring. Expect a clean error (exit 2) or a result. A panic, hang or memory growth is a finding.
- **Target 2, summary JSON.** Arbitrary JSON, and valid summaries with mutated fields, through `report::parse` and `render`. Expect an error naming the field, or valid HTML with no NaN or inf.
- Run each for 1 hour. Minimise any crash into a regression test in `tests/`.

### 6.1 Results (2026-11-21)

Four targets live in `fuzz/`, which is outside the published crate (`cargo package --list` shows no file from it). They were run with cargo-fuzz 0.13.2 and nightly 1.102, as a release build with AddressSanitizer, under WSL Ubuntu 24.04, for 1 hour each, from seeds made from the test data.

| Target | What it runs | Inputs | Result |
|---|---|---|---|
| `alignment` | SAM records through `reads::process` in every mode (ds/ss, UDG none/half/full, compat, filters), with MD and with a `-r` window; scores must be finite | 1.07 M, then 9.54 M after the fix | **B3** (allocation of 4×10¹⁶ bytes); after the fix, clean for 1 hour |
| `bgzf` | arbitrary bytes through the BGZF reader, and a write/read round trip, with and without the thread pool | 21,961 (about 20 per second; the inputs grew to 160 KB) | clean |
| `md` | MD/CIGAR reconstruction on arbitrary input, and MD written from a window giving back that window | 208 M | clean |
| `summary` | arbitrary JSON through `report::parse` and `render`; numeric SVG and data attributes must be finite | 12.0 M | clean |

The first runs produced three failures that were not pmdshodh bugs:
- **A hang in `alignment`:** the harness skipped CIGAR errors with `filter_map(Result::ok)`, and noodles' lazy iterator repeats an error without advancing. pmdshodh's own `collect` stops at the first error, and the binary rejects the record at once.
- **A "timeout" in `bgzf`:** libFuzzer reported it when the 1-hour limit ended a unit (a duration of 2^64 − 1 s). The input runs in 336 ms.
- **A panic in `summary`:** the harness flagged the text "inf<" in a fuzzed skip-reason name (`unsupported_cigarno_inf`), which the page shows escaped. That is not a number.

The harnesses were fixed and the affected targets rerun for a full hour.

**Added 2026-11-23 (O4): `cram`.** Arbitrary bytes as a CRAM file, read as `score -r` reads one. The header must match the reference (every sequence present, with the same length), containers are read four at a time, `Reference::cram_batch` loads their sequences, and `cram::decode` cuts each record's window. Every record is then scored in default and compat mode, and scores must be finite. The reference is `pairs.fa` and `past_end.fa` combined. The 8 seeds are `pairs.cram` and `past_end.cram`, plus six samtools 1.25 CRAMs of one synthetic SAM: CRAM 3.0 and 3.1, small slices in each, an embedded reference, and no reference. Same tools and flags as above, with `-rss_limit_mb=2048`.

| Target | What it runs | Inputs | Result |
|---|---|---|---|
| `cram` | CRAM bytes through header check, container batches, reference loading, pool decoding and `reads::process` in both modes | 0.5 M, then 2.17 M after the fix (7 chunks of 571 s; coverage 5289 to 5350; peak RSS 623 MB) | **B10** (a container claiming 1.65 GB); after the fix, clean for 1 hour |

- **Negative control:** with `cram.rs` from before the B9 fix, the target panics on the `past_end.cram` seed within seconds.
- **Two harness-only panics,** from first versions that built the reference repository by hand: noodles-cram panics (`slice.rs` 350 and 355) on a slice whose reference sequence ID or name is not in the header. The binary rejects both inputs first, with exit 2: "reference … does not match the input header", and "container reference sequence ID is not in the header". The harness now takes the binary's path, through `Reference`.
- **Memory:** before the B10 fix, peak RSS grew to about 1.4 GB over a chunk. A run with `-malloc_limit_mb=256` (195,254 inputs) found no other single large allocation, and ASan reported no leaks. After the fix, peak RSS is about 620 MB in every chunk.

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

### 7.1 Results (2026-11-21)

The run used `docs/benchmarks/audit_matrix.py` on `main` at `ba65da0`, which includes the B1 and B2 fixes, against a v0.3.0 baseline for speed. It ran under WSL Ubuntu 24.04 on an i5-13500H with 16 threads and 7 GB, with samtools 1.25 and PMDtools 0.60 (archive md5 `c059d7ab…`). The BER11_L_NOuser and Motala3 input md5s were checked. Result: **129 PASS, 0 FAIL.**

| Library | Reads | Compat scores vs PMDtools | Compat kept reads (t = 3) | Kept records unchanged, in order | `-r` vs MD |
|---|---|---|---|---|---|
| BER11_L_NOuser | 14,006,183 | 0 above 1e-6 (largest 6e-14) | 1,031,193 = PMDtools | 1,029,380 | 0 differ; 0 MD tags disagree |
| Motala3 | 23,622,338 | 0 (5e-14) | 301,867 = PMDtools | 301,722 | 71 differ, all within the 107 MD tags that disagree with the reference |
| JK2067 | 125,285 | 0 (2.9e-14) | 1,745 = PMDtools | 1,745 | — |
| UKY001 1240K (ss) | 14,461,790 | 0 (8.2e-14) | 7,169,568 = PMDtools | 5,327,200 | — |
| UKY001 MT (ss) | 281,906 | 0 (6.2e-14) | 162,734 = PMDtools | 117,380 | — |
| JK2782 | 1,197 | 0 (3e-14) | 188 = PMDtools | 188 | — |

The kept counts in default mode are lower than compat's, as the documented departures predict (design.md §5). For ss, this is the single count of each C site (SI Eq. S1).

Also checked, all PASS:
- **Read counts** against samtools for every library: input count; duplicate, secondary, supplementary, QC-fail, paired and reverse flags; unmapped skips; input = scored + skipped.
- **`score`, `damage` and `filter`** byte-identical at `-@ 1` and `-@ 8` for every library, apart from the recorded command line.
- **Filter outputs** pass `samtools quickcheck`.
- **Reports** match their JSON: 17 checks per library and 110 for the six-library report.
- **Damage tables**, compat, against PMDtools `--deamination --range 30` on the 200,000-read subsets of BER11_L_NOuser and Motala3, in UDG none, half and full: all within 5e-6, which is PMDtools' 5-decimal rounding.
- **All 17 aligner-sweep BAM and CRAM files:** `-@ 1` equals `-@ 8`, and input counts equal samtools'. The `calmd -e` file needs `-r` (its `=` bases), and passes with it.

Expected input errors (exit 2), as documented:
- `calmd_eq.bam` without `-r`: "read has '=' bases …";
- BER11_L with its USER-treated libraries: "MD is longer than the CIGAR", the broken MD tags noted in docs/benchmarks.md.

Speed and memory, `filter -t 3`, v0.3.0 against `main`, alternating, median of 3:

| | BER11_L_NOuser `-@ 1` | `-@ 8` | Motala3 `-@ 1` | `-@ 8` |
|---|---|---|---|---|
| v0.3.0 | 42.41 s | 8.64 s | 58.99 s | 9.94 s |
| main | 41.69 s | 8.53 s | 56.92 s | 9.95 s |
| Max RSS (both) | 11 MB | 38 MB | 11 MB | 34 MB |

There is no regression; the differences are within run-to-run spread.

**Default mode against compat, read by read (2026-11-22).** This check was not done in the first run, so it was done afterwards with `docs/benchmarks/audit_default_vs_compat.py`.
- **Setup:** a release build whose `src/`, `Cargo.toml` and `Cargo.lock` are identical to `main` at `a87e4e5` (0.3.1 plus B8); samtools 1.25.
- **Inputs:** the six libraries above, the three ENA libraries and the many-reference BAM of section 8.1, and the 13 aligner-sweep BAMs. `calmd_eq.bam` was left out, because default mode needs `-r` for its `=` bases.
- **Modes:** each was run in its own library mode (ds or ss), with UDG none, half (`--UDGhalf`) and full (`--CpG`). The three ss libraries also ran in ds mode, where the ss trigger does not cover every read. In all, 78 runs.
- **What is compared:** every record's default score with its compat score.
- **Departure triggers checked per read:**
  - the reference reconstruction fix: I or S together with D or `=` (blackbox.md 1.6c);
  - lowercase MD;
  - the reverse-read shift;
  - the CpG wrap at the read ends: UDG modes, with the reference along SEQ starting with G and ending with C;
  - ss mode with a C site (one term per site).
- **The rule:** a differing read without any of these triggers counts as unexplained.
- **Reads only one mode scores:** a read scored only by compat must be placed unmapped (B4). A read scored only by default must be hard-clipped (B6), 1-base (B8) or without reference bases (B5), the reads compat drops.
- **Re-scoring:** reads whose only triggers are lowercase MD or the reverse shift are scored again in compat with those removed: the MD uppercased, and the dropped characters written as N. The new score must equal default's.
- **Negative control:** a copy of the script with the CpG-wrap trigger disabled reports JK2782's 26 half- and full-mode reads as unexplained, and FAILs.

**Result: 237 PASS, 0 FAIL. No read is unexplained.**

| Library | ds/ss | differ: none | half | full | Causes (all modes) | Scored by one mode only |
|---|---|---|---|---|---|---|
| BER11_L_NOuser (14.0 M) | ds | 4,324 | 722,512 | 721,421 | CpG wrap; reconstruction 3,983 in none; lowercase MD 336 in none, 3 with the reverse shift | — |
| Motala3 (23.6 M) | ds | 3,757 | 764,537 | 764,373 | CpG wrap; reconstruction | — |
| JK2067 | ds | 2 | 7,762 | 7,762 | CpG wrap; reconstruction | — |
| JK2782 | ds | 0 | 26 | 26 | CpG wrap | — |
| UKY001 1240K (14.5 M) | ss | 14,367,584 | 8,845,859 | 7,396,923 | ss; CpG wrap; reconstruction | — |
| UKY001 MT | ss | 280,385 | 219,538 | 204,124 | ss; CpG wrap | — |
| I19359 MT | ss | 272,575 | 228,307 | 214,851 | ss; CpG wrap; reconstruction | — |
| I8525 MT | ds | 38 | 22,381 | 22,377 | CpG wrap; reconstruction | 20, no reference bases (B5) |
| I0071 MT | ds | 31 | 5,372 | 5,370 | CpG wrap; reconstruction | — |
| Many-reference BAM | ds | 37 | 9,738 | 9,737 | CpG wrap; reconstruction | 2,181, placed unmapped (B4) |
| UKY001 1240K, UKY001 MT, I19359 MT in ds mode | ds | 47, 0, 2 | 513,129, 8,117, 8,357 | 513,127, 8,117, 8,357 | CpG wrap; reconstruction | — |
| 13 sweep BAMs | ds | 0–141 each | 2,013–6,051 each | 2,012–6,051 each | CpG wrap; reconstruction | 18 and 9 hard-clipped (B6), in the two bwa mem files |

- **Re-scoring:** BER11's reads whose only causes are lowercase MD or the reverse shift all equal default after the change: 339 in none (3 of them with both), 91 in half and 71 in full. No other file has such reads.
- **Untreated ds mode:** only the reconstruction fix and lowercase MD change scores.
- **UDG modes:** the CpG wrap adds 2–6% of reads.
- **ss mode:** almost every read differs, as expected from the one-term-per-site model.
- **The reconstruction fix:** default's reconstruction is confirmed separately, because `-r` gives the same scores as MD on BER11 (see above).

**Limits:**
- For the CpG wrap and the ss model, the check is that the trigger is present: the script does not re-score without it.
- For ss, the trigger covers every read with a C site. So the ss-mode runs show mainly that reads without one are equal in both modes. The same libraries in ds mode test the other departures on them.
- Reads that carry a trigger but score the same in both modes are not counted.

**No new findings** from the real-data matrix.

## 8. Additional test data

The libraries in section 7 cover ds without UDG, likely-UDG ds, ss without UDG and a range of aligner outputs. Four kinds are missing. I'll look for public sources and list each with source, licence, size and date before any download:

1. **A UDG-treated single-stranded library** (decision 2.4). It is the only library type never compared with PMDtools.
2. **A double-stranded library with UDG treatment confirmed by its publication** (half and full, if both exist). Motala3 and JK2067 are only presumed treated, so their UDG-mode results cannot be checked against a known truth.
3. **A metagenomic or sedimentary-DNA BAM with many references:** thousands of contigs, short references and reads near contig ends. This tests multi-reference CRAM slices, reference loading and memory with many small sequences.
4. **An nf-core/eager output BAM after mapDamage rescaling:** base qualities lowered at damaged positions. PMDtools users run on such files, so `--min-baseq` and scoring should be checked on them.

Small files (under 100 MB) are preferred. A larger file is only taken if a small one does not exist, and only after you approve its size.

### 8.1 Results (2026-11-21)

**Sources:**
- The AADR v66.p1 annotation table (`v66.p1_1240K.aadr.PUB.anno`, Harvard Dataverse doi:10.7910/DVN/FFIDCW, CC0, 13.4 MB, md5 checked) lists each sample's library type. It has no `ss.half` or `ss.plus`, but 2,149 `ss.USER` samples (single-stranded, USER-treated).
- From it, three mitochondrial BAMs were taken from ENA (INSDC open data), downloaded 2026-11-21, with md5 checked against ENA:

| Gap | ENA run | Sample (publication) | Size | Reads |
|---|---|---|---|---|
| ss, USER-treated (item 1) | ERR7195621 | I19359.MT (Patterson et al. 2021, PRJEB47891) | 8.1 MB | 273,326 |
| ds, UDG-half (item 2) | ERR4589292 | I8525.MT (Narasimhan et al. 2019, PRJEB32466) | 40.4 MB | 755,825 |
| ds, full UDG (item 2) | ERR1995180 | I0071.MT (Lazaridis et al. 2017, PRJEB20914) | 6.2 MB | 141,943 |

- **Item 3, many references:** built locally. Motala3's 310,698 chr21 reads were re-aligned with `bwa aln` to chr21 cut into 3,511 contigs of 10 kb (contigs that are mostly N were dropped), giving a 25 MB BAM and a 9.5 MB CRAM.
- **Item 4, mapDamage rescaling:** already covered. Every read sampled from BER11_L_NOuser (the first 200,000) carries mapDamage's `MR` tag. That file has been in every PMDtools comparison and in the real-data matrix. mapDamage2 2.2.3 refused to rescale it a second time for that reason. The quality changes themselves were not verified.

**Checks.** `docs/benchmarks/audit_matrix.py`-style checks, with compat compared to PMDtools by read name, in the modes each library implies plus the untreated mode. **30 PASS, 0 FAIL** on the three ENA libraries:
- **Compat scores** equal PMDtools' on every read PMDtools scores, in every mode tried: ss and ss with `--UDGhalf`; `--UDGhalf` and untreated; `--CpG` and untreated. Largest difference 4.4e-14.
- **Kept reads** at t = 3 are identical: 18,347; 31,059 and 34,641; 214 and 1,379.
- **Damage tables** match `--deamination` within its 5-decimal rounding.
- **Counts**, `-@ 1` against `-@ 8`, and **report against JSON** all pass.
- **`-r` against MD:** reads differ only where MD tags disagree with rCRS. I19359 and I0071 carry many stale MD tags (29,326 and 23,462 respectively).
- **First-base C→T**, in the order of the table above: 23.35%, 8.83% and 0.79%, as their treatments predict.

**Many-reference BAM and CRAM:**
- `-@ 1` equals `-@ 8`, CRAM `-r` equals BAM MD scores (307,139 reads), and counts equal samtools'.
- CRAM decoding against 3,511 references takes 1.1 s, with 432 MB peak memory.
- This file found **B4**. After the fix, compat equals PMDtools on all 309,320 reads.

**Findings:**
- **B4** (fixed in #57).
- **B5** (fixed in #59): in I8525, 20 reads aligned only to `N` reference bases (MD `0N0N…`) are dropped by PMDtools with "bad reference sequence reconstruction", while compat scores them with 0 sites. The 18-read count difference on the sweep's `bwa_mem_pe.bam` (predating B4) may be the same; not yet checked.
- **B6:** it was not. Those 18 reads are hard-clipped supplementary alignments, which PMDtools drops and compat scored (section 10).

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
| B10 | Bug | A CRAM container header claims its body's length, up to 2^31 − 1 bytes, and noodles-cram allocates that length before reading the body. A corrupt or truncated CRAM of a few kilobytes could therefore cost up to 2 GB before the error, and fail outright where that much memory is not available. With `-@ 1` the corrupt container below happened to be rejected first, at its block checksum, as the first container is decoded before the next is read. | CRAM fuzz target (O4), 2026-11-23: a 7.7 KB input whose second container claims 1,650,128,970 bytes. The binary: exit 2 either way, but 1.6 GB peak RSS with `-@ 2` and `-@ 4` (5 MB with `-@ 1`). | Fix in PR #70: pmdshodh reads the container header itself, reads the body only as far as the input goes, and hands those bytes to noodles. A container longer than the input is an input error (exit 2), "container of N bytes is cut short: the input ends after M". The fuzz input now peaks at 5 MB with every thread count. BER11_L_NOuser.cram output is byte-identical, with no slowdown. Test `cram::tests::container_longer_than_the_input_is_cut_short` (2 GB claimed, 100 bytes given). |
| B9 | Bug | A CRAM whose mapped read runs past the end of its reference sequence made noodles-cram panic while converting the record (exit 101), before pmdshodh's own bounds check in `cram.rs` ran. This happens with `-@ 1` and more threads, and for a read past the end with or without a deletion. samtools 1.25 writes such a CRAM and reads it back. Placed unmapped reads (flag 4) are not affected. | Synthetic, found while writing the reason for the uncovered `cram.rs` check (§3.1), 2026-11-22. `samtools view -C` of an 8M read at position 17 of a 20-base reference: `score -r` panicked at `noodles-cram-0.102.0/src/record/sequence/iter.rs:30`. The same reads as SAM: exit 2. | Fix in PR #67: the reference window, which includes the bounds check, is now taken from the raw CRAM record before it is converted. The read is an input error (exit 2) naming it. Test `cli::cram_read_past_its_reference_end_is_an_input_error` (fixture `tests/data/past_end.cram`). Output is byte-identical before and after on the three sweep CRAMs and `multiref.cram`, in default and compat mode. |
| O4 | Observation | No fuzz target decodes CRAM: decision 2.3 covers alignment text, BAM records, MD/CIGAR and summary JSON. B9 was in that gap. | §6 targets. | Decided 2026-11-23: added the `cram` target (§6.1). Its first hour found B10, and after that fix it was clean for 1 hour. With `past_end.cram` as a seed, it finds B9 in seconds on the code from before that fix. |
| O3 | Observation | A mapped read with CIGAR `*` and an MD tag is scored from its MD by PMDtools (5.0297), and is an input error (exit 2) in pmdshodh, in both modes. Aligners do not write such records. | §4.5 probe (2026-11-22); test `edge::star_cigar_and_empty_md_on_mapped_reads_are_input_errors`. | Decided 2026-11-22: keep the error, in compat too. Listed in design.md §5. |
| O2 | Observation | Lowercase SEQ bases (SAM only; BAM cannot store them) are not recognised by PMDtools (`taaa` over MD `0C3` scores 0.0), and read as uppercase by pmdshodh in both modes (5.0297). | §4.5 probe (2026-11-22); test `edge::lowercase_iupac_and_n_read_bases`. | Decided 2026-11-22: compat does not reproduce it, a documented limitation in design.md §5. |
| O1 | Observation | CpG context is taken along SEQ with deleted reference bases skipped, in default mode as in PMDtools. So a C, a deleted A and then a G count as CpG, although the reference reads CpA there. This changes which sites count with `--CpG` and `--UDGhalf` for reads with a deletion right after a C. | §4.5 probe (2026-11-22): PMDtools and both modes give 5.0030 from 2 sites. Test `edge::cpg_context_skips_deleted_reference_bases`. | Decided 2026-11-23: keep it, in both modes. The read is the molecule, whose sequence carries the methylation, and insertions are already handled the same way (an inserted neighbour is not CpG). Reasoning in design.md §5. |
| B8 | Bug (compat) | PMDtools drops 1-bp reads silently (blackbox.md 1.6). `--pmdtools-compat` scores them, as default mode does by design (design.md §5). | §4.5 probe (2026-11-22): `1M` `T` over MD `0C0`, forward and reverse, and `1M` `A` all match: all dropped; 2 bp is scored. | Fix in PR #64: compat drops reads whose SEQ is 1 base (rule in blackbox.md 1.6g), counted as `too_short`; summary `format_version` 4 (agreed 2026-11-22). Test `edge::compat_drops_one_base_reads_like_pmdtools`. |
| B7 | Bug (compat) | PMDtools filters only when `--threshold` is above −10000 or `--upperthreshold` is below 1,000,000; otherwise it writes no reads, as if neither was given (black-box probes, blackbox.md 1.7b). `filter --pmdtools-compat` applied any bound, so `-t -1000000` kept every read where PMDtools keeps none. | Seen in B6's parity runs: PMDtools `--threshold -1000000` kept 0 of 1,848 scored reads, `-1000` kept all (2026-11-22). First noted as a side finding in blackbox.md §7. | Fix in PR #61: compat keeps no reads when neither bound is in PMDtools' range, with a stderr note; default unchanged. `bwa_mem_pe.bam`: PMDtools and compat both keep 197,046 at `-t -9999` and 0 at `-10000` and `-1000000`. `filter -t 3` records are identical before and after on `bwa_mem_pe`, BER11_L_NOuser and Motala3, in both modes. Test `compat_filter_ignores_thresholds_outside_pmdtools_range`. |
| B6 | Bug (compat) | PMDtools drops every read with a hard clip (`H`) in its CIGAR, whatever its flags (black-box probes, blackbox.md 1.6f). `--pmdtools-compat` scored them, as default mode does by design (D3). bwa mem writes hard clips on supplementary alignments. | Sweep `bwa_mem_pe.bam`: compat 197,064 reads vs PMDtools 197,046; the 18 extra were the file's only supplementary records, all hard-clipped (flags 2131, 2147, 2195, 2211). `bwa_mem_se.bam`: 95,177 vs 95,168, 9 extra (2026-11-22). | Fixed in PR #60: compat drops them as `unsupported_cigar`. Compat then equals PMDtools by read name and flag on both files (≤ 4.2e-14), and only those 18 and 9 lines leave the output. Test `compat_drops_hard_clipped_reads_like_pmdtools`. |
| B5 | Bug (compat) | PMDtools drops reads whose rebuilt reference it calls a "bad reference sequence reconstruction"; seen on reads aligned only to `N` (MD `0N0N…`). `--pmdtools-compat` scores them, with 0 sites and a score of 0. Kept reads are unaffected. | I8525.MT (ERR4589292): 20 of 755,825 reads, in UDG-half and untreated mode (2026-11-21). | Fix in PR #59: compat drops reads whose rebuilt reference has no uppercase A, C, G or T (rule in blackbox.md 1.6e), counted as `no_reference_bases`; summary `format_version` 3. Test `compat_drops_reads_without_reference_bases_like_pmdtools`. |
| B4 | Bug (compat) | `--pmdtools-compat` skipped reads flagged unmapped that still have a CIGAR and MD, which PMDtools scores because it ignores flags. bwa writes these for reads overhanging a reference end. | Many-reference BAM: compat 307,139 reads vs PMDtools 309,320 (2026-11-21). | Fix in PR #57: compat scores them from MD; default skips them, documented in design.md §5. Compat then equals PMDtools on all 309,320. |
| B3 | Bug | With `-r`, a CIGAR covering far more read bases than SEQ had its per-base vector built before the length check. Separately, CIGAR lengths were summed with a plain sum, which could wrap round to SEQ's length, and `fetch` could overflow at `start + len`. | Fuzz target `alignment`, then the release binary with `SEQ` = `AC`: `1M4000000000I` used 7.3 GB for 39 s before exit 2; `1M20000000000000000I` aborted with "memory allocation of 40000000000000002 bytes failed", exit 134 (2026-11-21). Without `-r`, both were exit 2 at once. | Fix in PR #55: the length is checked first with a saturating sum, `reference_span` saturates, and `fetch` uses `checked_add`. Three regression tests. Valid reads give identical output (ber2m, Motala3; default and compat). |
| B2 | Bug | `--pmdtools-compat -r` counted an `N` read base over an `N` reference base as a match when it built MD from the reference. calmd counts it as a mismatch (`5N2` for `ACGTANGT` over `ACGTANGT`). Reads that also hit the 1.6c reconstruction then scored differently from `calmd \| pmdtools`. | Property 4 on seed 48490, 592 reads with calmd-identical MD: compat from MD matched PMDtools 0.60 on all 583 scored reads (≤ 2e-14); compat `-r` differed on 15 reads, by up to 0.085. All 15 had N over N together with the 1.6c trigger (2026-11-21). Real libraries are not affected (byte-identical before and after). | Fix in PR #52: N never matches, as in calmd. Afterwards, 0 of 583 reads differ from PMDtools. |
| B1 | Bug | With `-r` and SAM/BAM input, a FASTA whose `.fai` no longer matches its line lengths (rewrapped after indexing) is read silently, giving wrong reference bases from the second 1 MiB block on. The CRAM path checks line lengths and stops with an error; the SAM/BAM path does not. `samtools faidx` also returns wrong bases here, and `samtools calmd` fails reading the whole sequence. | A synthetic 2.5 Mb reference, indexed at 60 columns and rewrapped to 70, with one read at 2,000,001: `score -r` gave 6.55 from 8 sites instead of −1.17 from 22, with exit 0 (2026-11-21). | Fix in PR #51: every line break is checked against the index, for the `-r` blocks and for CRAM from plain and bgzip FASTA. Regression test `reference::tests::fasta_that_does_not_match_its_index_is_an_error`. Outputs for correctly indexed FASTAs were unchanged on 4 real inputs, and speed was unchanged. |

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
