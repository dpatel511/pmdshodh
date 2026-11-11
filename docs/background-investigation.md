# Investigation: fast, maintained aDNA damage-score filter (PMDtools v2) — 2026-11-10

## Question
Is there room for a fast, maintained tool that scores each read for post-mortem damage (PMD) and filters a BAM, with a MultiQC module, for nf-core/eager? Or do pydamage, metaDMG and others already cover it?

## What PMDtools does today
- One Python script, `pmdtools.0.60.py` (2014 paper). It reads SAM from a `samtools view` pipe, needs MD tags, and has no dependencies.
- **Functions:**
  - per-read PMD likelihood score (damage + base quality + polymorphism)
  - `--threshold N` filtering
  - `--ss` (single-stranded libraries)
  - `--CpG` (UDG/USER-treated mammal data)
  - `--deamination` and `--platypus` damage tables (the second feeds an R plot)
  - `--first` (5′ damage rate with standard error, for library screening)
  - `--customterminus`
- **Still used:**
  - nf-core has a `pmdtools/filter` module (BAM + threshold → filtered BAM); eager's default threshold is 3.
  - The SPAAM community book teaches PMD > 3 as "reliably ancient".
  - bamRefine expects PMDtools output.
- **Neglected:** 9 open issues and 2 PRs, none fixed. Several came from eager developers:
  - #6 Python 3 incompatibilities (2018, 13 comments)
  - #13 Python 3 syntax warning (2024)
  - #8 and #14 `--first` errors (2021, 2025)
  - #15 empty deamination table (2025)
  - #7 `-n` fails to close the samtools pipe (2019)
  - #9 faulty `DS:Z:` matching (2021)
  - #5 MD handling (2018)
- **No MultiQC module** for PMDtools.

## The other damage tools
| Tool | Unit | Filters a BAM by read damage? | Notes |
|---|---|---|---|
| PMDtools | read | yes (SAM via pipe) | Python, unmaintained, bugs above |
| metaDMG-cpp `pmd` | read | **no**, prints `read<TAB>PMD: score` to stdout | C++, fast overall; `pmd` is single-threaded (`-@` ignored); no ss/CpG options for `pmd`; main focus is per-taxon damage (`getdamage`, `lca`, `dfit`) |
| ngsBriggs (MER 2025) | read (posterior probability) | writes a "scores" BAM, filtering not documented | C++, v0.0.1, build from source with htslib, no conda; needs an external contamination estimate; claims better accuracy than PMDtools |
| AuthentiCT | sample (+ per-read probability) | no | single-stranded only, not for UDG data; Python |
| pydamage | contig | filters contigs, not reads | for de novo assemblies; no MultiQC module |
| bamdam (2025) | taxon | extracts reads per taxon | post-LCA metagenomics authentication |
| mapDamage2 | sample | no (rescales qualities) | "no longer being actively developed"; has a MultiQC module |
| DamageProfiler | sample | no | Java; damage plots; has a MultiQC module |
| bamRefine / pmd-mask (Rust) | positions | masks bases, does not filter reads | |

**Verdict:** per-read scoring exists in three places:
- PMDtools: unmaintained.
- metaDMG `pmd`: scores only, no output BAM, single-threaded.
- ngsBriggs: early, hard to install.

None is a maintained, multi-threaded, conda-installable **filter** that writes the BAM directly, keeps PMDtools' flags, and reports to MultiQC. The gap is real but narrow. It is a v2 of a utility, not a new method.

## Proposed tool: **pmdshodh** (name chosen 2026-11-10; PMD + Sanskrit śodh "purify")
**Must have (v0.1):**
- Reads BAM/SAM/CRAM directly (no samtools pipe), multi-threaded, streams with constant memory.
- PMD score matching PMDtools' likelihood exactly, as a per-read tag (e.g. `DS:f`), and `--threshold` filtering into a BAM.
- Library modes: double-stranded, `--ss`, UDG-half/`--CpG`, plus trimming of the terminal bases.
- PMDtools-compatible flags (`--threshold`, `--ss`, `--CpG`, `--first`, `--deamination`), so the nf-core module can switch over with one line.
- Summary JSON/TSV: reads in/out, score histogram, C→T/G→A by position, 5′ damage with standard error. This feeds a **MultiQC module**.
- Tests: score parity with PMDtools 0.60 on public eager test data; speed benchmark on a large BAM.

**Later:**
- Optional Briggs-model posterior (ngsBriggs-style) using a supplied contamination estimate.
- Per-reference summary for metagenomic BAMs (several references in one file).
- Length-aware filtering presets.

**Language choice:** Speed is the reason to exist, so the scoring core should be in Rust (rust-htslib or noodles), wrapped for pip/Bioconda with maturin. A pure-Python + pysam version would not be much faster than PMDtools. Cython is a middle path if you want to stay in Python.

## How to stand out and get adopted
1. Fix what users hit: the open PMDtools bugs become regression tests.
2. Show identical scores and the speed-up on the same BAMs.
3. Contribute the nf-core module (`pmdfilter/filter`) and a MultiQC module; offer it to nf-core/eager (#1087 bamRefine shows they are open to damage-tool alternatives).
4. Bioconda recipe from day one.

## Risks and checks before starting
- **Small audience** (aDNA labs) but a loyal one; eager is the main channel.
- **metaDMG could add BAM output to `pmd`.** Ask its developers first; contributing filter output there is the alternative path.
- **Licence:** PMDtools is GPL-3.0. Re-implementing the published model is fine. Copying code means the new tool is GPL too.
- **Ask on nf-core Slack #eager** whether they would switch, and what output they want in MultiQC.
- **Fit:** builds on anvaya (aDNA) and your packaging/test habits; the MultiQC module is a visible second contribution.

## Sources
- PMDtools: https://github.com/pontussk/PMDtools , issues https://github.com/pontussk/PMDtools/issues
- nf-core pmdtools/filter: https://nf-co.re/modules/pmdtools_filter ; eager params: https://nf-co.re/eager/2.3.1/parameters/ ; eager issues: https://github.com/nf-core/eager/issues
- SPAAM authentication chapter: https://www.spaam-community.org/intro-to-ancient-metagenomics-book/authentication.html
- metaDMG-cpp: https://github.com/metaDMG-dev/metaDMG-cpp
- ngsBriggs: https://github.com/RAHenriksen/ngsBriggs , paper https://researchprofiles.ku.dk/en/publications/revisiting-the-briggs-ancient-dna-damage-model-a-fast-maximum-lik/ , preprint https://www.biorxiv.org/content/10.1101/2023.11.06.565746.full.pdf
- AuthentiCT: https://github.com/StephanePeyregne/AuthentiCT
- pydamage: https://github.com/maxibor/pydamage , paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8323603/
- bamdam: https://pypi.org/project/bamdam/
- mapDamage: https://github.com/ginolhac/mapDamage
- DamageProfiler: https://github.com/Integrative-Transcriptomics/DamageProfiler
- bamRefine: https://github.com/etkayapar/bamRefine ; pmd-mask: https://github.com/MaelLefeuvre/pmd-mask
- MultiQC modules: https://docs.seqera.io/multiqc/modules
