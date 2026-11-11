# Phase 1: literature, issues and related tools — 2026-11-11

Purpose:
- Write down the PMD scoring model exactly as published, the biology behind it, and what later papers found about it.
- Record what PMDtools' users, issues and neighbouring tools show about real-world behaviour.

That way pmdshodh can be built from the paper (clean room), and its design can cite evidence.

**Clean-room rule** (agreed 2026-11-11). PMDtools source code is not read. Its behaviour is learned from:
- the paper and its SI;
- the README and `--help` output;
- GitHub issue text;
- running it on test reads.

Exposure so far: issue #5 quotes four lines of the 0.60 parser, and PR #16's description explains one condition in the `--platypus` code. Both were read as bug reports. Neither is used as implementation material.

Status markers:
- **Paper**: stated in a cited paper or its SI (read for this review).
- **Doc/issue**: stated in a README, `--help`, or GitHub issue (not peer reviewed).
- **Derived**: my calculation from the paper's equations, not checked against PMDtools.
- **Open**: not answered by any source read so far. To be settled by black-box tests of PMDtools 0.60 (phase 3) or a design decision (phase 4).

## 1. Damage biology (what the model approximates)

**Briggs et al. 2007** ([doi:10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)):
- Cytosine deamination (C→U, read as T) is concentrated in single-stranded overhangs at molecule ends.
- In double-stranded library prep, T4 polymerase fills in 5′ overhangs. Deamination therefore shows up as C→T near the read's 5′ end, and as G→A near the 3′ end (from the complementary strand).
- Neandertal estimates:
  - single-stranded deamination δss = 0.68;
  - double-stranded deamination δ = 0.0097;
  - nick rate ν = 0.024;
  - overhang length λ = 0.36 (geometric; mean about 1.6–1.8 nt).
- About 21% of reference C are read as T at the first 5′ position. C→T falls quickly over the first ~10 positions and stays above background further in.

**Single-stranded libraries** (Gansauge & Meyer 2013, [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038)):
- Strands are sequenced directly, so deamination appears as C→T at **both** read ends and there is no G→A signal at the 3′ end.
- This is why PMDtools has `--ss` (model in section 2.5).

**UDG treatment** (Briggs et al. 2010, [doi:10.1093/nar/gkp1163](https://doi.org/10.1093/nar/gkp1163)):
- UDG + Endo VIII removes uracils, so most damage disappears.
- Methylated CpG cytosines deaminate to thymine, not uracil, so damage survives at CpG sites.
- This is the basis of `--CpG` for UDG-treated mammalian data.
- In PMDtools issue #7, the author adds that CpG sites "retain almost the entire deamination signal" after UDG.

**Partial UDG / "UDG-half"** (Rohland et al. 2015, [doi:10.1098/rstb.2013.0624](https://doi.org/10.1098/rstb.2013.0624)):
- Damage is restricted to the terminal nucleotides and nearly removed in the interior.
- The default decay (P = 0.3) does not describe these libraries.
- PMDtools has an undocumented-in-README `--UDGhalf` option, which nf-core/eager passes for UDG-half data (section 4.3). What it changes is **Open**.

**Conditional damage as a contamination check** (Meyer et al. 2016, [doi:10.1038/nature17405](https://doi.org/10.1038/nature17405); described in the `condamage` README, paper not read):
- Deamination at the two ends of a molecule is thought to be independent. So the C→T rate at one end, among reads that show C→T at the other end, estimates damage in endogenous reads only.
- Compare that with the overall rate. A large gap suggests undamaged contamination.
- This is cheap to compute in the same pass as scoring. It is a candidate summary metric (section 7).

## 2. The PMD score model (Skoglund et al. 2014)

**Source:** Skoglund P, Northoff BH, Shunkov MV, Derevianko AP, Pääbo S, Krause J, Jakobsson M. *Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal.* PNAS 111(6):2229–34, 2014. [doi:10.1073/pnas.1318934111](https://doi.org/10.1073/pnas.1318934111).
- Main-text equations were read from the PMC version (images of Eqs. 1–5).
- SI equations were read from the SI PDF (`docs/1318934111_pnas.201318934SI.pdf`, downloaded 2026-11-11 from PMC).

### 2.1 Which positions count (Paper)
- Only two kinds of aligned site enter the likelihood:
  1. reference **C**, read base **C or T**;
  2. reference **G**, read base **G or A**.
- All other positions are ignored, except that they count towards the distance z.
- z is measured from the **5′ end** for reference-C sites and from the **3′ end** for reference-G sites. The first position has z = 1.
- When computing z, alignment positions with a gap in the read (deletions relative to the reference) are not counted.

### 2.2 Double-stranded equations (Paper)
Per site, with polymorphism rate π, sequencing error ε and damage rate D_z:

1. P(Match | z) = (1 − π)(1 − ε)(1 − D_z) + (1 − π) ε D_z + π ε (1 − D_z)
2. P(Mismatch | z) = 1 − P(Match | z)
3. D_z = (1 − P)^(z − 1) · P + C
4. L(M | S_i) = P(S_i | M) for the observed state (match, or C→T/G→A mismatch) of site i
5. PMDS = ln( ∏_i L(M_PMD | S_i) / ∏_i L(M_NULL | S_i) )

Notes on the equations:
- The null model M_NULL is the same with D_z = 0 for all z. Then P(Match) = (1 − π)(1 − ε) + π ε.
- The three terms of Eq. 1 are: no event; a sequencing error that reverted a damage-induced T back to C; a sequencing error that reverted a true polymorphism.
- Parameters:
  - π = 0.001;
  - ε = (1/3) · 10^(−Q/10) from the base's phred quality Q, divided by three because only one of the three possible errors gives the modelled mismatch;
  - **P = 0.3**, **C = 0.01**. The authors say P could be raised for older, more damaged samples.
- "log" in Eq. 5 is the natural logarithm.
- A read with no informative sites (I = 0) has an undefined score. The paper discarded such reads.

### 2.3 What the equations imply (Derived; Q = 30, default parameters)

> **Update after phase 3:** PMDtools 0.60 does not compute exactly these numbers. Its null model uses a background damage of 0.001 instead of 0, which lowers every mismatch by about 0.68, and `--ss` counts each site twice. See [blackbox.md](blackbox.md) §1–2. The values below are the paper's equations, not PMDtools output.
Damage rate by position: D_1 = 0.31, D_2 = 0.22, D_3 = 0.157, D_5 = 0.082, D_10 = 0.022, D_20 = 0.0103 (approaching C = 0.01).

| Site contribution to PMDS | z = 1 | z = 2 | z = 3 | z = 20 |
|---|---|---|---|---|
| C→T (or G→A) mismatch | +5.45 | +5.11 | +4.78 | +2.17 |
| match | −0.37 | −0.25 | −0.17 | −0.010 |

**Worked example** (candidate test vector): C→T at 5′ z = 1, reference-C matches at z = 4 and 9, reference-G match at 3′ z = 2 → PMDS = **5.056**.

What follows from this:
- One high-quality terminal C→T alone scores above 5, so it passes the usual threshold of 3.
- C = 0.01 is about 30 times the Q30 error term, so even an interior mismatch adds about +2.2. **Two interior C→T mismatches anywhere in the read pass threshold 3** with no terminal damage. This matches AuthentiCT's finding (section 3) that PMDtools calls more single-internal-substitution reads ancient.
- Matches cost little, so read length barely lowers the score of a read with a mismatch.
- Base quality matters a lot: at low Q, ε grows and a mismatch counts for less. Any quality cap or offset must match PMDtools exactly for parity (**Open**).

### 2.4 Single-stranded model (Paper, SI Eq. S1)
In single-stranded libraries the only observed deamination mismatches are C→T at both ends, and it is impossible to tell 5′ from 3′ overhang damage. The SI therefore uses two damage processes per site: D_z (distance z to the 5′ end) and D_y (distance y to the 3′ end), both defined by Eq. 3:

P(Match | z, y) = (1 − π)(1 − ε)(1 − D_z)(1 − D_y) + (1 − π) ε D_z (1 − D_y) + (1 − π) ε D_y (1 − D_z) + π ε (1 − D_z)(1 − D_y)

Points about this model:
- The SI text concerns C→T only. Whether `--ss` still scores reference-G sites (G→A) is **Open**. The README says `--ss` makes PMDtools "instead look for C->T mismatches also in the 3'-end".
- **Derived**, for a 30-nt read, assuming y = read length − z + 1:
  - C→T at z = 1 or z = 30: **+5.47** each;
  - C→T at z = 15: +2.92;
  - match at either end: −0.38.
- How y is counted (with insertions, deletions and clips) is **Open**, as for z.

### 2.5 Base-quality adjustment (Paper, SI Eq. S2)
The SI also gives a PMD-aware base-quality rescaling for genotyping:

Q_adjusted = −10 · log10( 1 − (1 − D_z)(1 − 10^(−Q/10)) )

- They showed it reduces C/T and G/A genotype ambiguity in Vindija 33.16 mtDNA (SI Fig. S10).
- pmdshodh does not rescale qualities (CLAUDE.md §7), so this is recorded for completeness only.
- Whether PMDtools exposes it (for example as an `--adjustbaseq`-style option) is **Open**; check `--help` in phase 2.

### 2.6 How the paper itself filtered input (Paper, SI)
These are the paper's analysis settings, not necessarily PMDtools defaults:
- MD tags were recomputed with `samtools calmd`.
- For PMDS distributions they **excluded bases with quality < 20**, and **excluded clipped alignments and alignments with gaps** altogether.
- Mapping quality ≥ 30 was required (≥ 90 for ANFO-mapped Neandertal data). The mtDNA contamination experiment applied no mapping-quality filter.
- In the autosomal experiment, bases with Q < 30 were excluded from PMDS computation.
- 454 data without qualities were given Q = 20 everywhere.
- Paired reads were merged with ≥ 11 nt overlap. Overlapping bases had their qualities summed, and the higher-quality base was kept.

Implication: the published numbers come from settings stricter than "score every read". Clipped and gapped reads were never scored in the paper. How PMDtools treats them by default is **Open**.

### 2.7 Points still unspecified (Open; for black-box tests in phase 3)
1. **Reverse-strand reads.** SAM stores them reverse-complemented. Is z measured from the molecule's 5′ end (the right end of the alignment), with C/G and C→T/G→A swapped accordingly? metaDMG fixed exactly this in October 2026 (section 5.1), so it is easy to get wrong.
2. **Insertions** in the read: do they count towards z and y?
3. **Soft-clipped bases**: do they count towards z and y, or is the read skipped (as in the paper's analysis)?
4. **Unmerged paired-end mates.** The 3′ end of a mate is not the molecule's 3′ end. metaDMG's `pmd` keeps only single-end or collapsed reads by default (`-a`), a hint about expected practice.
5. **N bases, IUPAC codes in the reference (e.g. eager #324 shows `n` and `Y` in MD), quality 0 and very high Q.**
6. **Reads with I = 0**: what score is printed, and are they kept at `--threshold 0` or below?
7. **`--ss`**: are G sites dropped, and how is y computed?
8. **`--CpG`, `--UDGhalf`, `--range`, `--refseq`**: their effect on which sites are scored.
9. **Filters**: mapping quality and `--requirebaseq` defaults, and whether they change scores or only which bases or reads are used.
10. **The `DS:Z:` tag** (issue #9): what PMDtools writes into it, and when it re-reads it.

## 3. What later work found about the PMD score

**Skoglund et al. 2014 (the original), main text:**
- PMDS > 5 kept ~15–20% of reads from Neandertals and Neolithic Scandinavians, but only ~0.01–0.02% from present-day individuals and 0.27% from 100-year-old remains.
- In mtDNA mixtures, PMDS ≥ 5 gave no detectable contamination from 1% to ~93% input contamination.
- PMDS ≥ 3 removed the visible effect of contamination on a PCA.
- Okladnikov 2 contamination fell from 10.2% to 1.3% at PMDS ≥ 0.
- They chose hard thresholds rather than posteriors, because the unknown starting contamination changes what any threshold means.

**Skoglund et al. 2014, SI:**
- A PMDS threshold of about **2** gave a 100-fold enrichment of ancient over present-day reads. A threshold of about **6** was needed for the same enrichment against 100-year-old "contamination".
- A simple rule ("C→T within n bp of the 5′ end") reached at most ~20-fold enrichment (n = 1). PMDS reached ~5,000-fold while keeping ~10% of authentic reads.
- So the right threshold depends on how damaged the contaminant is. This argues for reporting the score distribution rather than only "reads ≥ 3".
- (The SI's figure numbers in the text and captions are offset by one for these figures: the text cites Fig. S3 and S9, the captions are S9 and S10.)

**AuthentiCT** (Peyrégne & Peter 2020, [doi:10.1186/s13059-020-02123-y](https://doi.org/10.1186/s13059-020-02123-y)). Single-stranded only; an HMM with overhang, internal single-stranded and double-stranded states.
- Classification is similar to PMDtools at false-positive rates below 0.02.
- Where they differ, it is mostly on reads with one internal C→T. PMDtools calls more of those ancient; the authors argue many are polymorphisms or errors.
- Adjacent C→T substitutions are not independent, but the PMD model assumes they are.
- A simple rule (C→T in the first or last three positions) caught ~30% of endogenous reads at 0.4% false positives. Both models caught ~50% at the same rate.
- Not applicable to UDG-treated libraries.

**ngsBriggs** (Zhao et al. 2025, [doi:10.1111/1755-0998.14029](https://doi.org/10.1111/1755-0998.14029)). Per-read posterior from damage plus read length; needs an external contamination estimate.
- It beat PMDtools in simulated ROC curves when ancient and modern length distributions were well separated. PMDtools cannot classify undamaged ancient reads.
- Mismatch-matrix runtime: 14,000–36,000 s for PMDtools against 170–260 s for ngsBriggs. These are their numbers; PMDtools' options and version were not reported.
- Not for single-stranded or USER-treated data. The length advantage fades when lengths overlap.

**Cuesta-Aguirre et al. 2024** ([doi:10.1371/journal.pone.0311115](https://doi.org/10.1371/journal.pone.0311115)). Amplicon mtDNA pipeline.
- Chose **PMDS ≥ 1** because higher thresholds "removed almost all reads".
- Cites two PMDtools versions (0.23.2 and 0.60) without explanation.

**PyDamage** (Borry et al. 2021, [doi:10.7717/peerj.11845](https://doi.org/10.7717/peerj.11845)). Read-level PMD filtering before de novo assembly loses data and lowers assembly quality.

**mapDamage2** (Jónsson et al. 2013, [doi:10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)). Sample-level Bayesian damage estimates and quality rescaling. Relevant to the damage-summary design.

Seen in search results but not read, so not used as evidence:
- Martiniano et al. 2020 (Genome Biology; reference bias under simulated damage);
- a 2024/2025 Genome Biology paper comparing read-end masking with transversion-only analysis.

**No independent benchmark of the threshold.** I found no study that tests PMDtools' threshold of 3 across damage levels, read lengths, library types and quality profiles. "PMD ≥ 3" is a convention from the 2014 experiments, not a calibrated cut-off. Searched: PubMed, Consensus and web; not a systematic review.

## 4. PMDtools on GitHub (Doc/issue)

### 4.1 Repository facts (checked 2026-11-11)
- Licence: **GPL-3.0**. The LICENSE file was added in February 2018. A 2018 comment (eager #36) says the code was unlicensed before that.
- Last code push: April 2021. 23 stars, 9 forks, 11 open issues and PRs; none merged since 2018.
- **Versions are confusing:**
  - The only tag is `0.60` (commit `18b5901`, June 2018).
  - Master is 11 commits ahead. One of those commits ("Add files via upload", 2018-12-11) changed `pmdtools.0.60.py`; this is the change that makes `--first` work on master but crash in the tagged release (issue #8, eager #349).
  - `pmdtools --version` on the tagged release reports **0.50** (issue #7, eager #349).
- **Bioconda** packages the **tagged 0.60 tarball** unmodified (md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`, build 5, `noarch`). nf-core's `pmdtools/filter` uses `pmdtools:0.60--hdfd78af_5` with Python 3.9.7.

**Parity target:** what eager users actually run is the Bioconda tag 0.60 under Python 3. Parity should be defined against that exact archive, not GitHub master. This matters for `--first`.

### 4.2 Issues and PRs (all read)
| # | Type, state | Year | What it shows | Relevance to pmdshodh |
|---|---|---|---|---|
| 16 | PR, open | 2026 | `--platypus` puts 3′ CpG C→T (ss data) into the non-CpG group, so the CpG 3′ series is always zero; plot also drops a terminal point. Validated on the Bacho Kiro A11197 ss library. | Bug in damage tables, not scores. Our CpG-context tables must classify C by its downstream base. Bacho Kiro A11197 is a candidate public ss test dataset (check size and licence). |
| 15 | issue, open | 2025 | `--deamination` table empty on a Strobealign BAM. | Aligner-specific input (MD or flags) silently gives empty output. pmdshodh must count and report skipped reads by reason (CLAUDE.md §7). |
| 14 | issue, open | 2025 | `--platypus`/`--first`: variable `n` undefined (the tag-vs-master bug again). | Confirms the tagged release is broken for `--first`. |
| 13 | issue, open | 2024 | Python ≥ 3.12 `SyntaxWarning` from regex strings. | Python packaging debt. Not relevant to a Rust core. |
| 12 | issue, closed by discussion | 2023 | Merged-read damage plots showed only a 3′ peak. The author suspected the merging. | Damage plots depend on read processing. A summary should state the read types it counted (merged / single / paired). |
| 11, 10 | PRs, open/closed | 2021 | Attempts to fix #9. | — |
| 9 | issue, open | 2021 | PMDtools searches the **whole SAM line** for `DS:Z:` to re-use scores from earlier runs; one read in a 60 GB BAM whose quality string contained `DS:Z:` broke it. | PMDtools' score tag is `DS:Z:`. Tag parsing must be field-aware (htslib/noodles handle this). Tag name for pmdshodh is still to be agreed (CLAUDE.md §7). |
| 8 | issue, open | 2021 | `--first` errors in the tagged release, works on master. | See 4.1. |
| 7 | issue, open | 2019 | `-n N` exits early, breaking the samtools pipe (SIGPIPE, exit 141); eager had to ignore exit codes 141 and 143. The author recommends `--platypus` and `--first --CpG` for UDG data. | Reading BAM directly removes the pipe. "Stop after N reads" for quick damage estimates should exit cleanly with status 0. |
| 6 | issue, open | 2018–22 | Python 3 decoding errors on binary input; the author wanted minimal dependencies (no pysam). Two unofficial Python 3 ports exist (boulund, rsh249). | Users asked for a maintained version for four years. Also, nobody had a test dataset to check that a port's output was unchanged — exactly the parity test pmdshodh needs. |
| 5 | issue, open | 2018 | MD tags like `^CCTTTC` straight after a deletion, without the `0` separator (written by MALT), are misparsed: the deletion swallows the following mismatches and shifts the alignment. | The SAM spec requires the separator, so such MD is malformed. pmdshodh should check MD against CIGAR and report or skip mismatching reads by reason, never score them silently. |
| 1–4 | PRs/issue | 2018 | Python 3 conversion, Bioconda recipe, divide-by-zero fix on empty data. | Empty input must be handled. |

### 4.3 How nf-core/eager uses PMDtools
- **nf-core module `pmdtools/filter`** (MIT). The pipeline is `samtools calmd <bam> <ref> | pmdtools --threshold T --header [args] | samtools view -Sb`.
  - MD is always recomputed against the reference, so in eager every read has MD and scores reflect the supplied reference.
- **Options eager passes** (issues #324, #349, #844, #989):
  - `--threshold`, `--header`;
  - `--CpG` (always for full-UDG data) and `--UDGhalf` (for half-UDG);
  - `--deamination --range 10 -n 10000` for damage tables;
  - `--platypus` (added 2021);
  - `--refseq <masked reference>` for SNP-capture data. A reference with capture SNPs masked as `N` stops real transition SNPs near read ends being scored as damage, which reduces reference bias. eager #844 found this was silently not applied for a while.
- **eager #324** (open since 2019, seen again in 2023 with ~25–30 bp reads): PMDtools exits with status 1 on some reads. The only workaround offered was a minimum-length filter. Root cause unknown. Candidate for our test suite: very short reads, and IUPAC or `n` bases in MD.
- **eager #989** (DSL2 plan): damage *estimation* (range, max reads, platypus) is separated from damage *manipulation* (threshold, reference mask). Design hint: pmdshodh could expose these as separate subcommands, or one pass producing both outputs.
- **MultiQC** has modules for `damageprofiler` and `mapdamage`, but none for PMDtools, metaDMG or AuthentiCT. AuthentiCT issue #4 (2023) asks for JSON output "for reading into tools like MultiQC", so users want this kind of output.

## 5. Related tools (Doc/issue; licences checked 2026-11-11)

| Tool | Language, licence | Last push | What it does | Notes for pmdshodh |
|---|---|---|---|---|
| **metaDMG-cpp `pmd`** | C++, **no licence file** | 2026-11-09 | Per-read PMD score to stdout (`read<TAB>PMD: score`) | See 5.1. No licence means its code cannot be reused; only its documentation was read. |
| **ngsBriggs** | C++/htslib, no licence file | 2024-05 | Briggs-model fit; read posterior; writes BAM | v0.0.1, no release; issue #1 (2024) asks for a release so it can go on Bioconda, unanswered. |
| **AuthentiCT** | Python, GPL-3.0 | 2023-10 | ss-only contamination estimate; per-read likelihood | Requires SAM, not BAM. Install errors in issues #2/#5. JSON-for-MultiQC request #4 unanswered. |
| **condamage** | C/htslib, no licence file | 2022 | Damage patterns conditional on terminal deamination (Meyer 2016) | Idea for a summary metric. |
| **pmd-mask** | **Rust** (rust-htslib), GPL-3.0 | 2025-05 | Masks likely-damaged bases using mapDamage rates | A Rust aDNA tool. Its bug history is instructive (5.2). Masks bases, which pmdshodh must not do. |
| **libtaph** | C++, no licence | 2026-10 | Reference-free damage from FASTQ; **classifies libraries as ds / ss / unknown** | pmdshodh must never guess library type (CLAUDE.md §7). A future *warning* when observed patterns contradict the declared mode is a possible feature, not a default. |
| **DamageProfiler** | Java, GPL-3.0 | 2023-11 | Damage plots, length distributions | Has a MultiQC module. Useful for checking our C→T/G→A tables in phase 3. |
| **mapDamage2** | Python/R, MIT | 2025-07 (marked inactive) | Damage estimates, quality rescaling | Has a MultiQC module. |
| **MADERA** | Rust, MIT | 2026-03 | Reference-free metagenomic clustering using damage features | Not a read filter. |
| **CUNA** | Python | 2025-09 | Deep-learning deamination detection from **Nanopore** signal | A different data type. Out of scope, but shows aDNA moving beyond Illumina. |

No other PMDtools re-implementation was found on GitHub: searched "pmdtools", "postmortem damage score", "ancient DNA damage", "deamination", and Rust + "ancient DNA". The only forks are code.google exports and the two unofficial Python 3 ports.

### 5.1 metaDMG `pmd`: an active, overlapping implementation
From the metaDMG-cpp README and commit messages (code not read):
- On **2026-10-19** four commits revised `pmd`:
  - "Fix PMD likelihood and alignment handling";
  - "Match PMDtools informative mismatches";
  - "Run PMD model regression in test suite";
  - "Document PMD usage and fixes reported by @miwipe".
- The README says earlier versions had "incorrect handling of base-quality scores and strand orientation". This is evidence that strand handling and quality handling are where implementations go wrong (Open points 1 and 5).
- Still **single-threaded** (`-@` is "accepted for compatibility; PMD is currently single-threaded"). Output is scores only; no filtered BAM and no summary.
- Defaults: unmapped reads are skipped. `-a 1` (default) keeps only single-end or collapsed reads. `-q` sets minimum mapping quality (off by default).
- metaDMG issues #24, #31 and #44 show hard read-length limits (assertion on reads ≥ 200 bp; a 500 bp limit removed in 2026). pmdshodh must have no fixed read-length limit.

**Risk update:** metaDMG is moving towards PMDtools parity right now. If it adds BAM output and threads, pmdshodh's gap narrows. The background doc's suggestion to talk to its developers before building is more pressing. Contacting them is your decision.

### 5.2 Bugs other tools hit (to test for from day one)
- **pmd-mask:**
  - indels ignored in position counting (#10);
  - panic on reads extending past the reference length (#13);
  - segfault on a bad reference path (#8);
  - a rust-htslib FASTA reader memory leak (#15).
- **metaDMG:**
  - fixed read-length limits;
  - "could not recreate the sequence" from MD/CIGAR mismatches (#33);
  - silently empty output files (#48).
- **PMDtools:**
  - whole-line tag matching (#9);
  - malformed MD (#5);
  - empty outputs (#15);
  - pipe exit codes (#7);
  - crashes on short reads (eager #324).

## 6. Summary of evidence by question

| Question | Best evidence | Strength |
|---|---|---|
| Exact ds and ss likelihood | Skoglund 2014 Eqs. 1–5, SI Eq. S1 | Paper |
| Default parameters | π 0.001, P 0.3, C 0.01, ε from Q/3 | Paper |
| Behaviour on clips, indels, reverse strand, pairs, I = 0 | none | Open → black-box tests |
| Does threshold 3 work? | 2014 experiments; SI enrichment curves; AuthentiCT and ngsBriggs comparisons | Indirect; no independent benchmark |
| Which version users run | Bioconda tag 0.60, Python 3 | Doc/issue, verified recipe |
| Demand for a maintained replacement | PMDtools #6 (4 years of Python 3 requests), #7/#8/#14/#15, eager #324/#349/#567 | Doc/issue |
| Competitors moving | metaDMG `pmd` fixes, Oct 2026 | Doc/commit messages |

## 7. What this means for pmdshodh (inputs to phase 4 design)

1. **Parity target:** the Bioconda PMDtools 0.60 tag archive (md5 above) under Python 3. Document which `--first` behaviour we match (the tag crashes).
2. **Resolve the ten Open points in 2.7 by black-box tests** before writing scoring code. Prioritise reverse strand and base quality, where metaDMG went wrong.
3. **Make parameters visible.** P, C and π are fixed in the paper but sample-dependent in reality. Report them in the summary. Whether they are settable is a phase 4 scope decision.
4. **Report the score distribution**, not only reads kept. Thresholds in use range from 0 to 6 (2014 main text and SI; Cuesta-Aguirre 2024).
5. **Report library mode, read types counted, and terminal mismatch rates** (C→T by 5′ position, G→A by 3′ position; CpG/non-CpG split). Consider the conditional-damage rate (Meyer 2016) as a cheap contamination hint.
6. **Paired-end and clipped reads need an explicit, documented policy.** The paper did not score clipped or gapped alignments at all, and metaDMG defaults to single or collapsed reads only.
7. **Reference masking (`--refseq`) is used by eager** for SNP-capture data. It needs a decision: support it, or document how to get the same effect.
8. **Robustness requirements from other tools' bugs:**
   - field-aware tag parsing;
   - MD/CIGAR consistency checks with per-reason skip counts;
   - no read-length limit;
   - clean exit with "first N reads";
   - clear errors for bad paths;
   - test reads that extend past the contig end;
   - very short reads;
   - IUPAC/`n` in MD.
9. **Known model weaknesses are documentation, not silent fixes:** interior-mismatch sensitivity, the independence assumption, and blindness to undamaged ancient reads. Any alternative model is a separate, opt-in score, never "the PMD score".
10. **Speed claims need our own measurements.** The ngsBriggs runtime table is evidence that PMDtools is slow, not a number pmdshodh can quote.
11. **Candidate public test data:** Bacho Kiro A11197 (ss; used in PR #16). Size, licence and source still need checking before any download.

## 8. Gaps remaining
- metaDMG preprint Appendix 1 (its PMD definition) not read; bioRxiv returns 403. Its README and commits were used instead.
- Meyer et al. 2016 not read directly; described via the condamage README.
- No paper evaluates `--CpG` or `--UDGhalf`.
- PMDtools `--help` not yet captured. It needs PMDtools installed (phase 2), which also gives the full option list.
- Searches: PubMed, Consensus (3 results per query on this plan), web search, and GitHub repository/issue search. Not a systematic review.

## Sources

Papers:
- Skoglund et al. 2014, PNAS, and SI — [doi:10.1073/pnas.1318934111](https://doi.org/10.1073/pnas.1318934111) (PMC3926038)
- Briggs et al. 2007, PNAS — [doi:10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)
- Briggs et al. 2010, Nucleic Acids Res — [doi:10.1093/nar/gkp1163](https://doi.org/10.1093/nar/gkp1163)
- Gansauge & Meyer 2013, Nat Protoc — [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038)
- Rohland et al. 2015, Phil Trans R Soc B — [doi:10.1098/rstb.2013.0624](https://doi.org/10.1098/rstb.2013.0624)
- Meyer et al. 2016, Nature — [doi:10.1038/nature17405](https://doi.org/10.1038/nature17405) (not read; via condamage README)
- Jónsson et al. 2013, Bioinformatics — [doi:10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)
- Peyrégne & Peter 2020, Genome Biol — [doi:10.1186/s13059-020-02123-y](https://doi.org/10.1186/s13059-020-02123-y)
- Borry et al. 2021, PeerJ — [doi:10.7717/peerj.11845](https://doi.org/10.7717/peerj.11845)
- Zhao et al. 2025, Mol Ecol Resour — [doi:10.1111/1755-0998.14029](https://doi.org/10.1111/1755-0998.14029)
- Cuesta-Aguirre et al. 2024, PLoS One — [doi:10.1371/journal.pone.0311115](https://doi.org/10.1371/journal.pone.0311115)
- Michelsen et al. 2022, bioRxiv — [doi:10.1101/2022.12.06.519264](https://doi.org/10.1101/2022.12.06.519264) (abstract only)

Repositories:
- PMDtools — https://github.com/pontussk/PMDtools (README, issues, PR descriptions, commit messages; no source)
- Bioconda recipe — https://github.com/bioconda/bioconda-recipes/tree/master/recipes/pmdtools
- nf-core module — https://github.com/nf-core/modules/tree/master/modules/nf-core/pmdtools/filter
- nf-core/eager issues #36, #324, #349, #567, #597, #844, #989, #1029 — https://github.com/nf-core/eager/issues
- metaDMG-cpp — https://github.com/metaDMG-dev/metaDMG-cpp (README, issues, commit messages)
- ngsBriggs — https://github.com/RAHenriksen/ngsBriggs; AuthentiCT — https://github.com/StephanePeyregne/AuthentiCT; condamage — https://github.com/grahamgower/condamage; pmd-mask — https://github.com/MaelLefeuvre/pmd-mask; libtaph — https://github.com/genomewalker/libtaph; DamageProfiler — https://github.com/Integrative-Transcriptomics/DamageProfiler; mapDamage — https://github.com/ginolhac/mapDamage; MADERA — https://github.com/michi-sxc/MADERA; CUNA — https://github.com/iris1901/CUNA

Bibliographic details for PubMed-indexed papers were retrieved from PubMed.
