# Phase 2: feature survey — 2026-11-11

Purpose: list what PMDtools and neighbouring tools offer, what pipelines actually use, and what a MultiQC report needs, then propose a v0.1 scope for discussion. Builds on [literature.md](literature.md); nothing here is decided until agreed in phase 4.

Sources:
- PMDtools 0.60 `--help`: captured by running the **Bioconda source archive** (`PMDtools/archive/0.60.tar.gz`, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`, matching the recipe) with Python 3.13.14 on Windows. The script was run, not read. `--version` prints `pmdtools.0.60.py v0.50`. Python 3.13 also printed a `SyntaxWarning` (issue #13) that echoes one source line, the same line quoted in issue #5.
- nf-core/eager 2.5.3 `main.nf` and `nextflow.config` (tag 2.5.3, released 2025-03-18). eager `dev` branch (DSL2 rewrite, unreleased, 1,772 commits ahead of master).
- nf-core/modules `pmdtools/filter`.
- MultiQC `damageprofiler` and `mapdamage` modules.
- READMEs of metaDMG-cpp, ngsBriggs and AuthentiCT.
- SAM tags specification (`hts-specs/SAMtags.tex`).

## 1. PMDtools 0.60 option inventory (from `--help`)

The README documents 8 options; `--help` lists **58**. Each is grouped below with a *proposed* disposition for pmdshodh:
- **v0.1**: needed for parity or for eager.
- **later**.
- **no**: outside scope, or forbidden by CLAUDE.md §7 (never alter kept reads).

### 1.1 Scoring model
| Option | Help text (abridged) | Proposed |
|---|---|---|
| (default) | double-stranded model | v0.1 |
| `--ss` | single stranded | v0.1 |
| `--CpG`, `--UDGplus` | only use Cs and Gs in CpG context | v0.1 |
| `--UDGhalf` | CpG context only, but first and last base used regardless of context | v0.1 |
| `--UDGminus` | use all bases (placeholder) | v0.1 (accept as the default) |
| `--noCpG` | don't use Cs and Gs in CpG context | later |
| `--PMDpparam` | p in the geometric PMD distribution | v0.1 (default 0.3; parity needs it reported) |
| `--PMDconstant` | constant C | v0.1 (default 0.01) |
| `--polymorphism_ancient` | polymorphism between ancient individual and reference | v0.1 (see note) |
| `--polymorphism_contamination` | polymorphism between contaminants and reference | v0.1 (see note) |
| `--adjustss` | strand-aware base-quality penalisation for ss libs | **Open**: what it does to scores |
| `--PMDSprim`, `--PMDSprimthreshold` | no description | **Open**: undocumented variant score |
| `--Ecoli`, `--EcoliCpG`, `--Leipzigsimple` | special terminus rules | no (undocumented, niche) |

**Note on π.** The paper uses one π = 0.001. PMDtools exposes **two** polymorphism rates, one per model. That suggests M_PMD and M_NULL may use different π. This is **Open** and must be resolved by black-box tests, not assumed.

### 1.2 Read and base selection
These change which reads or bases are scored, and so change scores.

| Option | Help text | Proposed |
|---|---|---|
| `-m/--requiremapq` | min mapping quality | v0.1 |
| `-q/--requirebaseq` | min base quality | v0.1 |
| `--noclips`, `--onlyclips` | clip handling | v0.1 for `--noclips`; **Open**: the default includes clipped reads? |
| `--noindels`, `--nodeletions`, `--noinsertions`, `--onlydeletions`, `--onlyinsertions` | indel handling | v0.1 for `--noindels`; others later |
| `--notreverse` | no reverse-complement alignments | later |
| `-c/--chromosome` | only this chromosome | later (users can pass a region) |
| `--readlength`, `--minlength`, `--maxlength` | length filters | later |
| `--minGC`, `--maxGC` | GC of aligned reference | later |
| `--perc_identity` | min percent identity | later |
| `-n/--number` | stop after N reads | v0.1 for damage summaries only, with a clean exit (issue #7) |
| `-r/--refseq` | reference sequence | **Open**: role; eager DSL2 masks via calmd instead (section 2) |

### 1.3 Filtering and output of reads
| Option | Help text | Proposed |
|---|---|---|
| `-t/--threshold` | "only output sequences with PMD score **above** this threshold" | v0.1. **Open: `>` or `≥`?** Reads exactly at the threshold depend on it. |
| `--upperthreshold` | only output sequences with score below this | v0.1 (cheap; gives "keep likely contaminants") |
| `--header` | output the SAM header | n/a: pmdshodh always writes a header |
| `--writesamfield` | add `DS:Z:<PMDS>` field, overwrite if present | v0.1 as an optional score tag (name to agree, section 5) |
| `-p/--printDS` | print PMD scores | v0.1 (per-read score table) |
| `-b/--basic N` | keep reads with C→T within N bp of 5′ | later (simple filter; SI compares it with PMDS) |
| `--terminal` | keep reads with C→T at either end | later |
| `--customterminus` | custom terminus rule (README example `0,-1`) | later |
| `--dry` | print SAM input without filters | no |
| `--samtoolspath` | path to samtools | no (pmdshodh reads BAM directly) |

### 1.4 Damage tables and estimates
| Option | Help text | Proposed |
|---|---|---|
| `-d/--deamination` | base frequencies at reference C/G positions | v0.1, as part of the summary JSON |
| `--range` | positions from terminus (default 30) | v0.1 |
| `--platypus` | big list of base frequencies (CpG split), for `plotPMD.v2.R` | later (CpG split in summary). PR #16 bug noted. |
| `--first` | first-position deamination rate with standard error | v0.1 (eager added it; broken in the tagged release) |
| `--estimate`, `--estimatebase` | "two-terminus estimate of contamination" | later (close to Meyer 2016 conditional damage) |
| `--basecomposition` | base composition | later |
| `--stats` | summary statistics to stderr | replaced by the summary JSON |
| `--printalignments`, `--verbose` | debugging output | later (a debug mode helps parity work) |

### 1.5 Read-modifying options (excluded by CLAUDE.md §7)
`-a/--adjustbaseq`, `--adjustbaseq_all`, `--maskterminaldeaminations`, `--maskterminalbases`, `--maskss`, `--flagss`.

These change qualities or bases, so they are **no** for pmdshodh. Users who need them have mapDamage rescaling, bamUtil trimBam, or pmd-mask.

## 2. What pipelines actually run

### 2.1 nf-core/modules `pmdtools/filter` (MIT)
`samtools calmd <bam> <fasta> [args] | pmdtools --threshold T --header [args2] | samtools view [args3] -Sb -`
- Container `pmdtools:0.60--hdfd78af_5` (Bioconda tag, Python 3.9.7, samtools 1.14).
- `calmd` always recomputes MD against the supplied FASTA.

### 2.2 nf-core/eager 2.5.3 (current release, DSL1)
- **Filter step:** `samtools calmd | pmdtools --threshold ${pmdtools_threshold} ${treatment} --header | samtools view -Sb`.
- **Damage step:** `samtools calmd | pmdtools --deamination [--platypus] --range ${pmdtools_range} ${treatment} -n ${pmdtools_max_reads}`, with a `trap` that turns exit 141 into success (issue #7).
- **Defaults:** threshold 3, range 10, max reads 10,000, platypus off.
- **Reference masking:** an optional BED masks the FASTA (`mask_reference_for_pmdtools`), and `calmd` then runs against the masked FASTA.

**Treatment selection** is `def treatment = udg ? (udg == 'half' ? '--UDGhalf' : '--CpG') : '--UDGminus'`.
- `udg` is always one of the strings `'none'`, `'half'` or `'full'` (validated at input).
- In Groovy any non-empty string is true, so this line picks **`--CpG` for `udg = 'none'`**. Non-UDG libraries would be scored at CpG sites only.
- `--ss` is never passed, although `strandedness` is available in the same tuple.

**Status: inferred from reading the code, not run.** If true, PMD-filtered BAMs from eager 2.5.3 for non-UDG data (and for ss data) are not what users expect. This should be confirmed by a run before anyone is told. Contacting eager is your call (CLAUDE.md §13).

### 2.3 nf-core/eager `dev` (DSL2 rewrite, unreleased)
- **Treatment:** `ext.args2` passes `--UDGminus` / `--UDGhalf` / `--UDGplus` from `meta.damage_treatment`, plus `--ss` when `meta.strandedness == 'single'`. This fixes both 2.5.3 issues.
- **Parameters:**
  - `damage_manipulation_pmdtools_threshold` is an **integer**, default 3;
  - masking via `..._masked_reference` (FASTA) or `..._reference_mask` (BED → `bedtools maskfasta`).
- **PMDtools is only used for filtering.** Damage calculation moved to `damagecalculation_tool` = `damageprofiler` (default) or `mapdamage`.
- **After filtering:** `samtools index` and `samtools flagstat` on the filtered BAM; flagstat feeds MultiQC.
- eager issue #989 (open) plans separate PMDtools roles for estimation (range, max reads, platypus) and manipulation (threshold, mask).

### 2.4 What this means for a drop-in replacement
- **Minimum flags** for an eager DSL2 swap: `--threshold`, `--UDGminus|--UDGhalf|--UDGplus` (and `--CpG`), `--ss`, read BAM, write BAM.
- **Masking works only through MD.** Eager masks by recomputing MD against a masked FASTA. If pmdshodh reads the BAM directly without `calmd`, masking is lost unless pmdshodh takes either:
  - (a) a BED of positions to ignore, or
  - (b) a reference to score against.

  Option (a) leaves kept reads untouched. `calmd` piped into PMDtools rewrites MD and NM in the output BAM, a side effect the current module has silently.
- **Damage tables** are not needed for eager DSL2, but they are for 2.5.3 users and anyone replacing `--deamination`/`--first`.
- A new nf-core module (e.g. `pmdshodh/filter`) would replace the three-command pipe with one command.

## 3. Feature matrix of related tools

| | PMDtools 0.60 | metaDMG-cpp `pmd` | ngsBriggs 0.0.1 | AuthentiCT | DamageProfiler | mapDamage2 |
|---|---|---|---|---|---|---|
| Unit | read | read | read (posterior) + sample | read + sample | sample | sample |
| Input | SAM on stdin, MD required | SAM/BAM/CRAM, MD required | BAM + reference | SAM (stdin or file), MD required | SAM/BAM/CRAM + FASTA | BAM + FASTA |
| Per-read score | `DS:Z:` (`--writesamfield`) or printed | `read<TAB>PMD: score` on stdout | likelihoods; "scores" BAM | score column; `ST:Z:` state tag | — | — |
| Writes filtered BAM | SAM to stdout, via samtools | **no** | `-obam` (filtering undocumented) | no | no | rescaled BAM |
| Threads | no | no (`-@` ignored) | `-nthreads` (recalibration only) | no | no | no |
| ds / ss / UDG modes | ds, `--ss`, CpG, UDGhalf | not documented | ds only (biotin/non-biotin models) | ss only, not UDG | both (plots) | ds, `--single-stranded` |
| Damage summary | tables, `--first`, `--platypus` | via `getdamage` | parameters | parameters | JSON (MultiQC) | text files (MultiQC) |
| Needs | samtools, Python | htslib | htslib, Eigen; build from source | Python, numpy, scipy, pandas | Java | Python, R |
| Bioconda | yes (tag 0.60) | yes | **no** (issue #1, 2024) | yes | yes | yes |
| Licence | GPL-3.0 | none stated | none stated | GPL-3.0 | GPL-3.0 | MIT |
| MultiQC module | no | no | no | no (requested, #4) | yes | yes |

## 4. MultiQC: what existing damage modules show

- **DamageProfiler** (parses one JSON per sample):
  - 5′ C→T and 3′ G→A frequency by position (line plots);
  - forward and reverse read-length distributions;
  - General Stats columns for 5′ C→T and 3′ G→A at the 1st and 2nd base, and read-length median, mean and sd.
  - JSON keys used: `metadata.sample_name`, `metadata.version`, `dmg_5p`, `dmg_3p`, `lendist_fw`, `lendist_rv`, `summary_stats`.
- **mapDamage**: the same 5′/3′ plots and length distributions, parsed from `5pCtoT_freq.txt`, `3pGtoA_freq.txt` and `lgdistribution.txt`.
- **samtools flagstat** is already reported. Eager DSL2 shows read counts after PMD filtering this way.

A pmdshodh MultiQC module would add what no module shows today:
- reads in, kept and skipped by reason;
- the **PMD score distribution** (histogram) and the threshold used;
- library mode and model parameters;
- the 5′/3′ plots for the reads actually scored.

Reusing DamageProfiler's plot conventions (per-position percentages, 1st/2nd base General Stats) keeps reports consistent. A machine-readable JSON with `format_version` (CLAUDE.md §7) is what AuthentiCT users asked for (AuthentiCT #4).

## 5. Score tag name

- **PMDtools writes `DS:Z:<score>`.** `DS` is **not** defined in the SAM tags specification. It is also not in the local namespace (tags starting with X, Y or Z, or containing lowercase letters). So it could be given another meaning in a future spec.
- **Choices for pmdshodh** (CLAUDE.md §7 requires agreement before use):
  - (a) `DS:Z` for drop-in compatibility with tools that read PMDtools output;
  - (b) a local-namespace tag such as `XP:f` or `pm:f`, which is spec-safe and typed as a float;
  - (c) (b) by default, with an option to write `DS:Z`.
- **Open:** do any downstream tools read `DS:Z`? PMDtools itself re-reads it (issue #9); bamRefine was said to "expect PMDtools output" (background doc), unverified.

## 6. Proposed v0.1 scope (for discussion in phase 4)

**Must (v0.1):**
1. Read SAM/BAM/CRAM directly. Write BAM or SAM (`-o -` for stdout). Multi-threaded, preserving input order, bounded memory.
2. Double-stranded score matching PMDtools 0.60 (Bioconda tag) on every read, plus `--ss`, `--CpG`/`--UDGplus`, `--UDGhalf`, and `--UDGminus` as an accepted no-op.
3. `--threshold`, with PMDtools' exact comparison (`>` vs `≥`, to be tested); `--upperthreshold`.
4. Model parameters with PMDtools defaults: P, C and π (one or two values, to be tested). Always recorded in the summary.
5. `--requiremapq`, `--requirebaseq`, `--noclips`, `--noindels`, if tests show they affect scores in common use.
6. Optional per-read score tag (name to agree) and a per-read score table (like `--printDS`).
7. Summary JSON (`format_version`):
   - tool version, command, mode and parameters;
   - reads in, kept, and skipped by reason;
   - score histogram;
   - C→T by 5′ position and G→A by 3′ position (or C→T at both ends for ss);
   - first-position rate with standard error (`--first`).
8. Clear errors: missing MD tag, MD/CIGAR mismatch; no silent zero scores.
9. Accept PMDtools flag names, so the eager config needs minimal change.

**Should (v0.1 if cheap, else v0.2):**
- A position mask from a BED file, so eager's reference masking works without `calmd`.
- A "first N reads" mode for quick damage estimates.
- CpG / non-CpG split in damage tables (with PR #16's 3′ fix).
- A debug option printing per-site contributions (helps parity work and users).

**Later:** `--estimate`-style two-terminus contamination hint, length and GC filters, `--basic`/`--terminal`/`--customterminus`, per-reference summaries for metagenomic BAMs, an optional Briggs-model posterior.

**Won't:** base-quality adjustment, masking or flagging bases, `--dry`, `--samtoolspath`, the undocumented `--Ecoli*`/`--Leipzigsimple` variants.

## 7. Decisions needed from you
1. **v0.1 scope:** accept section 6 as the starting list for phase 4, or trim it?
2. **Score tag:** option (a), (b) or (c) in section 5?
3. **CLI shape:** subcommands (`pmdshodh filter`, `pmdshodh score`, `pmdshodh damage`) or one command with PMDtools-style flags? Section 11 of CLAUDE.md prefers subcommands. PMDtools compatibility suggests accepting its flag names within them.
4. **eager 2.5.3 `--CpG`/`--ss` finding:** confirm by a run during phase 3 before telling anyone? Contacting eager needs your authorisation.

## 8. Open questions carried to phase 3 (black-box tests)
- `--threshold` comparison (`>` vs `≥`); score printed for reads with no informative sites.
- One π or two; the defaults of `--polymorphism_ancient` and `--polymorphism_contamination`.
- Defaults of `--requiremapq` and `--requirebaseq`; whether clipped or indel-containing reads are scored by default.
- What `--UDGhalf`, `--CpG`, `--adjustss`, `--PMDSprim` and `--refseq` change.
- Plus the ten points in literature.md section 2.7.

Running these tests needs SAM input. No samtools is installed here, but PMDtools reads SAM text on stdin, so synthetic SAM files written by hand or by a script are enough. No BAM tooling is needed for phase 3.
