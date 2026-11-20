# Proposal: HTML damage report (0.3.0, roadmap milestone 5)

Status: **agreed 2026-11-20.** D1–D6 were accepted as recommended (section 6). Merged in PR #44 on 2026-11-20, not yet released. Validated on the eager test BAMs on 2026-11-20 (section 9).

## 1. Goal and scope

A per-run QC report, generated from the JSON run summary, that shows whether the damage filter did what it should. Roadmap milestone 5 (docs/roadmap.md):
- C→T and G→A by position;
- the PMD score histogram;
- read counts;
- the first-base C→T rate.

The milestone is done when the report has been checked on the nf-core/eager test BAMs, with the JSON still the source of truth.

In scope:
- reading one or more summaries;
- one self-contained HTML file;
- plots and tables of what the summary already holds.

Out of scope for 0.3.0:
- a MultiQC module (milestone 7; it will read the same JSON);
- fitted damage parameters (milestone 6);
- new fields in the summary;
- interactive JavaScript;
- any general damage profiling beyond the scored reads (roadmap positioning: QC around filtering, not a DamageProfiler competitor).

## 2. Input: what the summary holds today

`format_version` 2 (src/summary.rs `to_json`):

| Field | Used for |
|---|---|
| `tool`, `subcommand`, `command`, `input`, `input_format` | report header (provenance) |
| `mode` (library, treatment, pmdtools_compat), `model` (pmd_p, pmd_c, pi_ancient, pi_contamination, null_background) | header; which damage curve to plot (section 4) |
| `filters` (threshold, upper_threshold, min_mapq, min_baseq, no_clips, no_indels) | header; threshold lines on the score histogram |
| `reference` (path, bases_from, md_tags agree/disagree/absent) | header; MD agreement row when `-r` was used |
| `reads` (input, scored, kept, removed_by_threshold, no_informative_sites, skipped by reason, flags) | read-count table |
| `score_histogram` (bin_width, min, max, underflow, overflow, counts) | score histogram |
| `length_histogram` (scored, kept) | read-length plot |
| `damage` (range, five_prime / three_prime × reference C/G × read base A/C/G/T counts, first_position c_to_t/se/n) | damage plots and the first-base rate |

The tables hold raw counts, so the report computes rates itself (`design.md` §6). `damage` writes the JSON to `-o`; `filter` and `score` write it with `--summary`.

## 3. Interface (proposed)

```bash
pmdshodh report run.json [more.json ...] -o report.html
```

- **New subcommand `report`** (D1). It reads existing summaries and never re-reads a BAM. You can report on finished runs, and combine several samples, without rerunning anything.
- **Several summaries in one report** (D3): one section per sample, plus overlay plots of all samples' damage curves and score histograms. Each sample is labelled by the file name of its `input`. If two labels clash, the summary file's name is added.
- **Output** follows the same rules as the other outputs:
  - `-o` is required, and `-o -` writes to stdout;
  - the file is written by temporary name and renamed on success, and removed on a signal;
  - an output may not be one of the input summaries.
- **Exit status:** 2 for an unreadable, invalid or unsupported summary, naming the file and field. 1 for write failures.
- `--quiet` as elsewhere; there is no other new option.

## 4. Contents per sample

1. **Header:** pmdshodh version, command line, input and its format, reference, library and treatment, compat mode, model parameters, filters and thresholds.
2. **Read counts:**
   - input, scored, kept, removed by threshold, no informative sites;
   - skipped reads by reason;
   - flag counts;
   - with `-r`, the MD tags that agree, disagree or are absent.
3. **Damage by position.** Following DamageProfiler's MultiQC conventions (docs/features.md §4): percent, positions 1…`range`, with 5′ and 3′ panels.
   - **ds libraries:** 5′ C→T and 3′ G→A.
   - **ss libraries:** 5′ C→T and 3′ **C→T**. Damage shows as C→T at both ends there, and the summary holds 3′ C counts too. It is shown this way on the real ss library UKY001: 25.5% at the 5′ end, 17.7% at the 3′ end.
   - **Other substitutions:** a faint background line for the other substitutions, as a baseline for error.
   - **UDG modes:** a note that the tables count only CpG sites (and the terminal base for half), as in scoring.
4. **First-base C→T rate,** with its standard error and the number of sites.
5. **PMD score histogram:**
   - 0.5-wide bins from −20 to 40, with the counts below and above that range stated;
   - vertical lines at `--threshold` and `--upper-threshold`;
   - for `filter`, the fraction of reads kept.
6. **Read length:** histograms of scored and kept reads (kept only for `filter`).

## 5. How it is built

- **One HTML file with inline SVG**, written by pmdshodh itself.
  - There is no JavaScript, no external CSS or fonts, and no network access, so it opens on a cluster node or as an e-mail attachment.
  - It needs no new dependency.
- **Deterministic:** the same summaries give byte-identical HTML. There is no timestamp, and the version and command line come from the JSON.
- **PNG** (D2): not in 0.3.0. Rasterising SVG needs a large dependency (e.g. resvg), and HTML covers the use. The roadmap's "(and PNG)" would move to later, or be dropped.
- **Size:** a few tens of KB per sample.

## 6. Decisions needed (my recommendation first)

| | Question | Recommendation | Alternative |
|---|---|---|---|
| D1 | How is a report made? | `pmdshodh report` subcommand only | Also `--report FILE` on `filter`, `score` and `damage` (one step, but a second path to keep stable) |
| D2 | PNG output in 0.3.0? | No; HTML with inline SVG only | Add PNG with a rasterising dependency |
| D3 | Several samples per report? | Yes: per-sample sections plus overlays | One report per summary |
| D4 | Damage axis | Percent, positions from 1 (DamageProfiler and MultiQC convention) | Fractions, positions from 0 (as in the JSON and PMDtools `--deamination`) |
| D5 | Summary versions accepted | 2, and 1 where the field exists (1 lacks `input_format` and `reference`; shown as "not recorded") | 2 only, with an error naming the version |
| D6 | ss libraries' 3′ curve | 3′ C→T for ss, 3′ G→A for ds, from `mode.library` | Always G→A (PMDtools' table convention) |

## 7. Plan and checks (after agreement)

1. Branch `feat/report`. Write the CLI, tests and docs to the agreed interface.
2. **Unit tests:**
   - rates and percentages from known counts;
   - SVG coordinates for a small known table;
   - escaping of file names and command lines in HTML;
   - ss/ds curve choice.
3. **CLI tests:**
   - a report from `damage` on `tests/data/reads_random.sam`;
   - several summaries in one report;
   - the same output twice (byte-identical);
   - an invalid JSON or unknown `format_version` gives exit 2;
   - an output equal to an input gives exit 2;
   - `-o -`, a FIFO and the signal cleanup, as for the other outputs.
4. **Validation (the roadmap's "done when"):**
   - reports on the nf-core/eager test BAMs (mammoth mtDNA, human 1240k), checked by eye;
   - every plotted number checked against the JSON by a script.
   - Also on the ds benchmark libraries and the real ss library, to check the ss curve.
5. **Docs:** a README section with an example; `design.md`; CHANGELOG. Then release 0.3.0 the same way as the 0.2.x releases.

## 8. Context for the new session

State on 2026-11-20:
- **Released:** 0.2.4 is on GitHub and PyPI; `main` is at `561bb88`.
- **Bioconda:** the 0.1.0 recipe (bioconda/bioconda-recipes#71666) is still awaiting review; bump it after it merges.
- **Known limitation, documented:** `--pmdtools-compat` does not reproduce PMDtools on `=` read bases of reverse-strand reads, or on reads made only of `=` (design.md §5).
- **Bench machine:** WSL Ubuntu 24.04 at `/root/bench`, with:
  - PMDtools 0.60 (archive md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`);
  - samtools 1.25, bwa 0.7.19 and bowtie2 2.5.5 in the `pmdshodh-bench` conda environment;
  - benchmark BAMs in `data/` and the ss library in `data_ss/`;
  - the aligner sweep in `sweep/`.
- This file is untracked. Commit it with the docs PR once D1–D6 are decided.

## 9. Validation (2026-11-20, branch `feat/report`)

**Data:**
- the two nf-core/eager test BAMs used for parity (docs/blackbox.md §7): nf-core/test-datasets, branch `eager`, MIT licence, re-downloaded 2026-11-20 (SHA-256 `7bf408da…f0eb` and `65f68bde…5d4a`, unchanged);
- the real ss library UKY001.B0102 (MT), to check the ss curve.

These stay on the bench machine (`/root/bench/eager_report`) and are not committed.

**Runs** (WSL Ubuntu 24.04, release build of `451b13c`, as merged in PR #44):
- `filter -t 3 --summary` and `damage` on JK2782 (mammoth mtDNA, ds);
- `filter -t 3 --UDGhalf --summary` and `damage` on JK2067_downsampled_s0.1 (human 1240k, ds, UDG-treated);
- `damage --ss` on UKY001.B0102.MT.

The reports were: JK2782 alone (37 KB), JK2067 alone (28 KB), and all five summaries together (164 KB). Writing the five-sample report twice gave identical files.

**Number check:** a script (kept with the data, not committed) reads every `data-series` in each report. It recomputes from the raw counts in the summary named in that section:
- the 5′ and 3′ damage curves, other substitutions, score-histogram percentages and length bars;
- the first-base text;
- the read-count rows, the kept fraction and the threshold marks;
- and, for the overlays, each sample's three curves.

Results: 17, 17 and 83 PASS lines, 0 FAIL, with a tolerance of 0.0005 percentage points. Changing one plotted value by 0.002 made it report 1 FAIL. A rerun with a CSS-only change (B1 below) gave the same counts.

**By eye:**
- I viewed the overlays and per-sample sections in light and dark mode. One layout bug was found (B1 below).
- Rates for the first three positions (JSON counts):

  | Sample | 5′ C→T | 3′ G→A (ds) or 3′ C→T (ss) |
  |---|---|---|
  | JK2782 | 13.9, 7.0, 6.9% | G→A 15.8, 8.8, 3.7% |
  | JK2067 | 0.9, 0.4, 0.2% | G→A 1.3, 0.8, 0.3% |
  | UKY001 MT (ss) | 28.0, 16.8, 11.4% | C→T 19.5, 17.5, 15.8% (G→A at the 3′ end is 0.6%, as expected for ss) |

  JK2067's near-zero terminal damage fits the design sheet's "full UDG" rather than the README's "half" (docs/blackbox.md §7). The report shows the data as given.
- With `--UDGhalf`, the JK2067 curves beyond position 1 are noisy, because only CpG sites count there. The report says so in its note.

### Bugs and observations from this run

| # | Kind | What | Status |
|---|---|---|---|
| B1 | Bug (layout) | In "All samples", a legend wider than its 440 px chart widened the figure, so the overlay charts stacked into one column instead of two per row. It shows with long labels (5 samples with clashing names). Numbers are unaffected. | Fixed in PR #46: figures are capped at the chart width. |
| O1 | Observation | The read-length chart shares one count axis for scored and kept reads. When few reads are kept (JK2067: 109 of 12,473), the kept bars are hard to see. | Fixed in PR #47: kept reads have their own chart (same x range, own y axis). |
| O2 | Observation (data) | With `--UDGhalf`, curves beyond position 1 count only CpG sites, so they are noisy on small samples (JK2067). The report's note already says which sites are counted. | No change. |
| O3 | Observation (data) | JK2067's near-zero terminal damage fits "full UDG" (eager design sheet) rather than "half" (eager README), as already noted in docs/blackbox.md §7. | No change in pmdshodh. |

None of these change a computed number. The B1 and O1 fixes are in 0.3.0, the report's first release.
