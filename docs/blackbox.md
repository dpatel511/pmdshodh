# Phase 3: black-box behaviour of PMDtools 0.60 — 2026-11-11

Purpose: settle how PMDtools actually scores reads, so pmdshodh can match it (parity) and document every deliberate difference. Answers the Open points in [literature.md](literature.md) §2.7 and [features.md](features.md) §8.

## Method

- **Tool under test.** `pmdtools.0.60.py` from the Bioconda source archive (`PMDtools/archive/0.60.tar.gz`, md5 `c059d7ab4c24e1ce00f1eee4c2a4b669`).
  - Run with Python 3.12.3 under WSL Ubuntu 24.04. It needs `/dev/null`, so it does not run on native Windows.
  - `--version` prints `v0.50`.
- **Clean room.** The script was only run; its source was not opened. Python tracebacks from crashes print the single failing source line. Those lines appear in the issue notes below as error evidence only.
- **Inputs.** Synthetic SAM records written by `docs/blackbox/harness.py`. MD is computed from (reference, read, CIGAR) per the SAM spec. Scores were read from `--printDS`.
- **Reproduce:** `py -3 -I docs/blackbox/probeN.py`. Set the `PMDTOOLS_PY` environment variable to the Windows path of `pmdtools.0.60.py` first.
- **Comparison.** Against two models in `harness.py`:
  - the paper's equations (Skoglund 2014 Eqs. 1–5, SI Eq. S1);
  - a *fitted* model describing what PMDtools does (section 1).
- **No real data** was used. Parity on real BAMs is still to be done (section 5).

## 1. The score PMDtools 0.60 actually computes

**Fit.** These rules reproduce `--printDS` on 300 random reads to within **6 × 10⁻¹⁴** (floating-point noise) in both default and `--ss` mode, and on 200 more in default mode. The reads:
- forward and reverse strand;
- lengths 25–90 bp;
- base qualities 2–40;
- random reference composition including CpGs;
- C→T and G→A damage on 25–30% of sites.

They also reproduce every single-site probe at Q = 0–93 (probe1, probe2, probe7, probe9).

### 1.1 Notation and output
- Q is the phred base quality (offset 33).
- ε = 10^(−Q/10) / 3, with **no cap** (checked up to Q = 93) and no floor (Q = 0 gives ε = 1/3).
- D(z) = (1 − P)^(z−1) · P + C. Defaults P = 0.3 (`--PMDpparam`) and C = 0.01 (`--PMDconstant`).
- π_a = `--polymorphism_ancient` (default 0.001); π_c = `--polymorphism_contamination` (default 0.001).
- B = **0.001**, a background damage rate in the null model. It does not change with `--PMDconstant`; no option to change it was found.
- `--printDS` prints four tab-separated columns per read, **without the read name**: L_PMD, L_null, L_PMD / L_null, PMDS = ln(ratio).

### 1.2 Which sites count, and how distance is measured
- Sites come from the reference rebuilt from SEQ + CIGAR + MD.
- **Default (double-stranded):**
  - reference C with read C or T, scored from the 5′ end;
  - reference G with read G or A, scored from the 3′ end;
  - **in SAM orientation regardless of the reverse flag.** The ds model is symmetric under reverse-complementing, so this gives the molecule-oriented result anyway.
- **Distance:** z = 1-based position **in SEQ**, counting soft-clipped bases and inserted bases, not counting deletions. The 3′ distance is len(SEQ) − position + 1.
- **Not scored:** other mismatches (e.g. C→A), `N` in the read, and non-ACGT reference letters from MD (`n`, `Y`).

### 1.3 Per-site terms (default, double-stranded)
P_match(D, π) = (1 − π)(1 − ε)(1 − D) + (1 − π) ε D + π ε (1 − D).

For each informative site, using its own 5′ or 3′ distance z:

- PMD model: p₁ = P_match(D(z), π_a)
- Null model: p₀ = P_match(**B**, π_c)
- Contribution: ln(p₁ / p₀) for a match; ln((1 − p₁) / (1 − p₀)) for a C→T / G→A mismatch.

PMDS is the sum over sites. **A read with no informative sites scores exactly 0.0** (all columns print 1.0, 1.0, 1.0, 0.0).

### 1.4 `--ss` (single-stranded)
- Only C sites in **molecule orientation** are scored, and G sites are ignored:
  - on forward reads: reference C → read C or T;
  - on reverse-flag reads (flag 16): reference G → read G or A in SAM terms.
- With z and y the distances from the molecule's 5′ and 3′ ends (SEQ positions, as in 1.2), **each site contributes two terms**:
  1. the double-stranded term of 1.3 with D(z);
  2. an Eq. S1 term: PMD model P_match_S1(D(z), D(y), π_a); null model P_match_S1(B, B, π_c), where P_match_S1(a, b, π) = (1−π)(1−ε)(1−a)(1−b) + (1−π)ε a (1−b) + (1−π)ε b (1−a) + π ε (1−a)(1−b).

### 1.5 Treatment modes (combine with ds or `--ss`)
| Mode | Sites scored |
|---|---|
| default, `--UDGminus` | all C/G sites as above (identical outputs) |
| `--CpG`, `--UDGplus` | only C followed by reference G, and G preceded by reference C. Sites at a read edge, whose neighbour is unknown, are dropped. |
| `--UDGhalf` | CpG sites, plus the first base (5′ C) and last base (3′ G) in any context |
| `--noCpG` | drops the **G** of a CpG but keeps the **C**. This looks like a bug (expected: drop both). Not v0.1 scope. |

### 1.6 Read and base selection defaults
| Behaviour | Observed |
|---|---|
| Mapping quality | **No default filter** (MAPQ 0 is scored). `-m N` keeps MAPQ ≥ N and drops other reads. |
| Base quality | **No default filter.** `-q N` drops bases below N, not reads; a read can then score 0. |
| SAM flags | **Ignored entirely.** Unmapped (4), secondary (256), supplementary (2048), QC-fail (512), duplicate (1024), and paired first/second mates are all scored. |
| Soft clips, indels | Scored, and they count towards z (1.2). `--noclips`, `--noindels`, `--noinsertions`, `--nodeletions` drop the whole read. |
| Hard clips (`H`) | Read dropped, with a stderr message ("only supports cigar operations M, I, S and D"), exit 0. |
| Read length | No upper limit seen (5,000 bp scored). A **1-bp read is dropped silently**; 2 bp and longer are scored. |
| Missing MD tag | **Read dropped silently**, exit 0. |
| MD length ≠ CIGAR length | Scored without warning. |

### 1.6b Context, CIGAR and limits (probes 12–14, added during build step 3)
- **CpG context is taken along SEQ positions:**
  - the "next reference base" of a C is the reference base at the next SEQ position, and the "previous" one for a G is the previous SEQ position;
  - deleted reference bases are skipped (C, deleted A, then G counts as CpG);
  - an inserted or soft-clipped neighbour has no reference base, so the site is not CpG.
- **`--UDGhalf` terminal exception** applies to SEQ position 1 (a 5′ C) and the last SEQ position (a 3′ G). After a soft clip, the first aligned base is not terminal.
  - With `--ss` the exception applies only to the molecule's 5′ first base, not its last base.
- **`=` and `X`** CIGAR operations are scored like `M`. **`N` and `H`** make PMDtools drop the read with a stderr message.
- **Lowercase read bases** are not recognised: a lowercase `t` is not counted as damage, and a lowercase `c` with an all-match MD drops the read.
- **QUAL `*` or SEQ `*`:** the read is dropped silently (exit 0).
- **Crashes (exit 1, whole run lost):**
  - **any informative site with z ≥ 1000** (IndexError). Single sites at z ≤ 999 work, for C and G alike;
  - **many mismatches in one read** (e.g. 400 C→T in 800 bp): ZeroDivisionError, because the product of probabilities underflows to 0. A sum of logs gives 687.96.

### 1.6c Reference reconstruction bug with combined CIGAR operations (probes 15–16, fixture)
Found when whole-read parity was tested on 1,500 random reads with mixed CIGARs (tests/data/reads_random.sam).

**The rule:**
1. If the CIGAR has **no `I` and no `S`**, PMDtools lays the MD reference (deleted bases removed) along the read. **This is correct**, for `M`, `=`, `X` and `D` alike.
2. Otherwise it builds a gap mask by walking the CIGAR with a position counter:
   - `M` and `X` advance it;
   - `I` and `S` mark gaps and advance it;
   - `D` **also advances it**, as if deleted bases were in the read;
   - `=` is **ignored**;
   - `H`/`N` reads were already dropped.
3. It then fills the non-gap positions in order with the MD reference, **deleted bases included**. Matches take the read base; anything past the end is dropped.

**Consequences, whenever a read has an insertion or soft clip together with a deletion or `=` ops:**
- deleted reference bases are scored against the read bases that follow the deletion (a deleted G over a read A looks like G→A damage);
- later mismatches shift by the deletion length;
- gaps move. A trailing soft clip after `=` ops is treated as if it were at the start, and clipped bases can be scored as matches.

Mismatch counts on the fixture against a correct reconstruction:

| CIGAR features | Wrong scores |
|---|---|
| deletion + insertion | 36/38 |
| insertion + `=`/`X` | 24/24 |
| clip + `=`/`X` | 35/58 |
| clip + deletion | 37/161 |
| any single feature alone | 0 |

The rule reproduces **all 1,500** reads' PMDtools scores to ≤ 2.3e-14 (ds mode).

**CpG context wraps around the read ends (probes 17–20).**
- PMDtools scores ds reads in molecule orientation (reverse-complemented for flag 16). The site terms are identical in SAM orientation, which is why section 1.2 looked strand-independent.
- In molecule orientation, a **G at the first base takes the last base as its "previous" base** (Python's index −1). So it counts as CpG in `--CpG`/`--UDGhalf` when the last reference base is C. A C at the last base never wraps.
- In SAM orientation that is the first G of a forward read, or the last C of a reverse read.
- `--ss` scores only C sites and is unaffected.
- With this rule and 1.6c, **all 8 modes** of the fixture (ds, ss, ×none/half/full, `-q 20`, changed parameters) match PMDtools for all 1,500 reads (≤ 8.3e-14).

### 1.6d Damage table (`--deamination`), found in build step 5
- Rows are z = 0…range−1 from the terminal base. Columns:
  - CT CA CG CC: the fraction of read bases A/C/G/T at reference C, from the **5′ end**;
  - GA GT GC GG: the same at reference G, from the **3′ end**.
- Counted in **molecule orientation** (reverse reads reverse-complemented), over the reference as rebuilt in 1.6c, **skipping `N` read bases**.
- `--ss` gives the **same table** as the default; the table always uses 5′ C and 3′ G.
- `--CpG`: a C site's context is the **next read base** (must be G), while a G site's context is the previous **reference** base (C, wrapping at the read start). Scoring uses reference context for both.
- With these rules pmdshodh `damage --pmdtools-compat` reproduces PMDtools' default and `--CpG` tables on the 1,500-read fixture (30 positions, 5-decimal output).
- **Open:** the `--UDGhalf` table differs from these rules at positions 1–5 (up to 0.009). Several candidate rules were tried; none matched. Scores in UDG-half mode do match (1.6c).

### 1.7 Filtering output (`--threshold`)
- **Keeps PMDS ≥ threshold** (a read exactly at the threshold passes). `--upperthreshold` keeps PMDS **<** its value.
- Reads with no informative sites (PMDS 0) pass `--threshold 0`.
- Output records are passed through **unchanged**, and in input order.
- **No `@PG` line.** The header is written only with `--header`.
- `--writesamfield` (with `--threshold`) appends `DS:Z:<score rounded to 3 decimals>`, e.g. `DS:Z:4.637` or `DS:Z:5.03`. It replaces an existing `DS:Z`.

## 2. Differences between PMDtools 0.60 and the paper

1. **Null model background damage.** The paper's null model has D = 0. PMDtools uses a constant 0.001.
   - Effect: every C→T/G→A mismatch scores about **0.68 lower** than the paper's equations (Q40, π = 0.001). Matches score about 0.001 higher.
   - A single terminal C→T scores **5.03**, against 5.71 from the paper (Q40).
   - My worked example in literature.md §2.3 (5.056 at Q30, paper equations) is therefore **not** what PMDtools gives. PMDtools gives 4.637 for the same read at Q40, against 5.310 from the paper.
2. **Two polymorphism rates** (π_a in the PMD model, π_c in the null model) where the paper has one. With defaults both are 0.001, so this only matters if a user changes one.
3. **`--ss` does not implement SI Eq. S1 alone.**
   - Each C site gets the ds 5′ term **and** the S1 term: in effect the site is counted twice.
   - A terminal C→T scores **9.68** in `--ss` mode, against 5.05 from Eq. S1 alone.
   - So `--ss` scores are not on the same scale as ds scores, and the same threshold (e.g. 3) is far looser for ss libraries. Whether this is intended is unknown. It is not described anywhere I have read.
4. **Reads with no informative sites** score 0 and are kept at threshold ≤ 0. The paper says such reads were discarded.
5. **Clipped and gapped alignments are scored** by default, with clipped bases counted in z. The paper's own analyses excluded clipped and gapped alignments (SI).

None of these is a "bug" in a proven sense. They are what the released program does. For parity, pmdshodh must reproduce 1–5 by default and document them. Any change (e.g. a paper-exact `--ss`) must be opt-in and never called "the PMDtools score" (CLAUDE.md §7).

## 3. Reported issues: reproduction status

| Issue | Status | Evidence | pmdshodh behaviour to adopt |
|---|---|---|---|
| #5 MD without `0` after a deletion (MALT) | **Reproduced** | Same read: spec MD `4^CCTT0T0C24` scores 3.408; MALT-style `4^CCTTT0C24` scores −0.085, exit 0, no warning | Parse MD with the CIGAR. Reject or skip malformed MD **by reason**, never score silently. |
| #8 / #14 `--first` crashes in tagged release | **Reproduced** | `NameError: name 'n' is not defined`, exit 1, with and without `--CpG` | Implement first-position rate + SE from its documented meaning; parity cannot be tested against the tag. |
| #9 `DS:Z:` matched anywhere in the line | **Reproduced** | `DS:Z:` in QUAL → `ValueError`, exit 1 | Field-aware tag handling (htslib/noodles). |
| New: input already carrying `DS:Z` | **Found** | Without `--writesamfield`, any read with an existing `DS:Z` tag crashes (`NameError`, exit 1). Re-running PMDtools on its own tagged output fails. | Existing tags are kept untouched. Overwriting our own tag needs a decision. |
| New: empty input | **Found** | SAM with no reads (empty, or header only) → `IndexError`, exit 1. With `--threshold --header` the header is printed first, then the crash. | Empty input gives empty output, exit 0, and a summary with zero counts. |
| New: long or heavily damaged reads | **Found** | Site distance ≥ 1000 → IndexError; ~400 mismatches → ZeroDivisionError; both exit 1 | Sum log terms; no distance limit. |
| New: silent drops | **Found** | Missing MD, 1-bp reads, hard-clipped reads (stderr message only) dropped with exit 0 | Count every skipped read by reason in the summary. Missing MD is an error (CLAUDE.md §7). |
| #7 `-n` breaks the samtools pipe | Not applicable here | `-n 5` works when reading a file; the problem is upstream SIGPIPE | Read BAM directly; "first N reads" ends cleanly. |
| #15 empty deamination table | Not reproduced | `--deamination` works on synthetic input; the cause is probably in the reporter's BAM (Strobealign) | Report skipped reads by reason so empty tables explain themselves. |
| eager #324 exit 1 on short reads | Not reproduced | Reads of 1–5 bp and `Y`/`n` in MD do not crash | Keep these cases in pmdshodh's tests anyway. |
| #16 `--platypus` 3′ CpG classification | Not tested | Damage-table option; later scope | — |

`--deamination --range N` works and prints, per position starting at 0, the fractions CT, CA, CG, CC, GA, GT, GC, GG. Its exact semantics are not yet checked.

## 4. Status of the earlier Open points

| Open point (literature.md §2.7 / features.md §8) | Answer |
|---|---|
| Reverse-strand reads | ds: SAM orientation, equivalent by symmetry. `--ss`: molecule orientation (strand-aware). |
| Insertions in z | Counted (z is a SEQ position) |
| Soft clips | Counted in z and scored by default |
| Unmerged pairs | Flags ignored; each mate scored as if it were a whole molecule |
| N, IUPAC, Q extremes | Not scored; no quality cap |
| Reads with I = 0 | PMDS 0.0; kept at threshold ≤ 0 |
| `--ss` model | Section 1.4 (differs from SI Eq. S1) |
| `--CpG`, `--UDGhalf` | Section 1.5 |
| Filters | No defaults; `-m` drops reads, `-q` drops bases |
| `DS:Z:` tag | Written by `--writesamfield`, 3 decimals; an existing tag crashes the run unless overwritten |
| `>` vs `≥` | `≥` for `--threshold`, `<` for `--upperthreshold` |
| One π or two | Two (π_a, π_c), both 0.001 by default |
| `--adjustss`, `--PMDSprim`, `--refseq` | **Not tested** (not needed for v0.1 parity) |

## 5. What this means for pmdshodh

1. **The parity specification is section 1.** It is complete for ds and `--ss` with `--CpG`, `--UDGhalf` and `--UDGminus`, on reads with M/I/D/S CIGARs. It should be turned into pmdshodh's reference tests, including the 500 random reads used here (seeds in probe3 and probe8).
2. **Parity is not yet shown on real data.** The next check, when the time comes, is a full comparison on a public aDNA BAM (number of reads, largest difference, tolerance; CLAUDE.md §7). That needs a dataset choice and your approval for its size.
3. **Decisions needed, because parity and CLAUDE.md conflict:**
   - **Missing MD:** PMDtools drops silently; CLAUDE.md requires an error. Proposed: error by default, documented as a departure.
   - **Flags:** PMDtools scores unmapped, secondary, supplementary, QC-fail and duplicate reads. Proposed: match by default for parity, but count them in the summary. An opt-in flag skips them.
   - **Hard-clipped and 1-bp reads:** PMDtools drops them. Proposed: score hard-clipped reads (hard clips are not in SEQ, so z is unaffected). Documented departure: it fixes data loss, and kept reads would differ from PMDtools output.
   - **`--ss` double counting:** match by default for parity. Offer a paper-exact S1 score only as an explicitly named alternative, if at all.
   - **No `@PG` line:** pmdshodh adds one (allowed by CLAUDE.md §7). Records stay byte-identical, so BAM content differs from PMDtools only in the header.
4. **Score tag:** PMDtools writes 3-decimal text (`DS:Z`). A float tag (`f`) keeps full precision. Still for you to decide (features.md §5).
5. **eager 2.5.3 `--CpG` for non-UDG libraries** (features.md §2.2). Section 1.5 now shows the consequence if the inference is right: only CpG C/G sites are scored, so most damage evidence on non-UDG reads is ignored. The Groovy behaviour itself was not executed. Telling eager is your call.

## 6. Limits of this phase
- Synthetic reads only. No real BAM, no CRAM, no comparison with metaDMG.
- Tested on PMDtools under Python 3.12 (WSL). Python 2 behaviour was not tested; nf-core runs 3.9.
- Untested options: `--adjustss`, `--PMDSprim`, `--refseq`, `--estimate`, `--platypus`, damage-table semantics.
- The fitted model explains all observed scores, but "fits every probe" is not proof of identical code paths for inputs we did not try. Real-data parity (5.2) is the stronger test.

## 7. Real-data parity (build step 7, 2026-11-11)

**Data:** nf-core/test-datasets, branch `eager`, MIT licence, downloaded 2026-11-11 into the session scratchpad (not committed).

| File | Library | SHA-256 | Reads scored |
|---|---|---|---|
| `testdata/Mammoth/bam/JK2782_TGGCCGATCAACGA_L008_R1_001.fastq.gz.tengrand.fq.combined.fq.mapped.bam` | mammoth mtDNA capture, ds, non-UDG, merged | `7bf408da…f0eb` | 1,197 mapped (8,695 unmapped skipped) |
| `testdata/Human/bam/JK2067_downsampled_s0.1.bam` | human 1240k capture, ds, UDG-treated (README says half, design sheet says full) | `65f68bde…5d4a` | 12,473 |

Both files carry MD tags, so no reference was needed.

**Method:**
- Each BAM was converted to SAM with `pmdshodh filter -t -1000000` (all mapped reads, records unchanged) and fed to PMDtools 0.60 under WSL.
- PMDtools gave per-read scores via `--printDS`, and read names via `--threshold=-500` (same order).
- pmdshodh gave per-read scores via `score`.
- Script: `docs/blackbox/realdata_compare.py`. pmdshodh built at commit `4acea79` (debug build, Windows).

**Results (every read compared; tolerance 1e-9):**

| Mode | `--pmdtools-compat` | default |
|---|---|---|
| ds | 13,670 / 13,670 equal, max diff 3.0e-14 | equal (max 3.0e-14) |
| `--CpG` | equal (max 2.0e-14) | 772 reads differ (26 mammoth, 746 human = 6.0%), max diff 1.62 |
| `--UDGhalf` | equal (max 2.0e-14) | same 772 reads differ |
| `--ss` | equal (max 2.9e-14) | equal (max 2.9e-14) |

- **Every** default-mode difference is a read whose rebuilt reference starts with G and ends with C: the CpG wrap-around in 1.6c (772 of 781 such reads; why the other 9 did not change was not investigated). No other difference exists.
- The CIGAR reconstruction bug (1.6c) did not affect these files: they contain M/I/D/S operations, but no read combines them in the triggering way.

**Side finding:** `--threshold=-1000000` makes PMDtools print no reads at all, while `--threshold=-5` and `-500` work. Very negative thresholds appear to be treated as "unset". Not investigated further.

**Not covered:** single-stranded real data (none small and public was found), CRAM input, and BAMs from aligners that produce the triggering CIGAR combinations (e.g. bwa mem with soft clips and indels).
