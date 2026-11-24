# Damage-curve study: results

The study plan is docs/damage-curve-study.md. This file collects its results as they come, with the commands that produced them. There are no conclusions yet.

## 1. The tools are checked

The tools live in tests/sim (Python 3, standard library only):
- **simulate.py:** writes damaged and undamaged reads with a truth file.
- **check_tables.py:** checks pmdshodh's damage tables against the simulator's rates.
- **score.py:** the oracle scorers.

The checks were run 2026-11-24 with pmdshodh 0.3.3 (PyPI) under WSL Ubuntu 24.04, Python 3.12.

### 1.1 pmdshodh's damage tables show the simulated rates

For each configuration below, 100,000 reads were simulated on a random 1 Mb contig. pmdshodh's `damage` tables of the file were then compared with the expected rate at every position within the range (30) from both ends:
- ds: 5′ C→T and 3′ G→A;
- ss: C→T at both ends.

The expectation mixes damaged and undamaged reads by their length distributions, counts sites as pmdshodh does under each treatment, and adds sequencing errors from the base qualities.

| Run | Library, treatment | Damage model | Undamaged share | Largest \|z\| over 60 positions |
|---|---|---|---|---|
| base | ds, none | geometric A 0.30, q 0.30, C 0.01 | 0.5 | 2.20 |
| ss | ss, none | geometric | 0.5 | 3.24 |
| half | ds, UDG-half | terminal excess (+0.10 at the terminal base) | 0.5 | 2.42 |
| full | ds, full UDG | geometric | 0.5 | 2.43 |
| obs | ds, none | observed: BER11_L_NOuser's profile | 0.5 | 2.30 |
| sshalf | ss, UDG-half | geometric; undamaged reads longer (mean 70 bp) | 0.8 | 2.45 |
| sshalf2 | ss, UDG-half | terminal excess | 0.2 | 2.12 |

Every position is within 4 standard errors. With 420 positions in all, a largest |z| of 3.24 is what chance gives.

```
py -3 -I tests/sim/simulate.py --out base --seed 1
pmdshodh damage --quiet base.sam -o base.damage.json
py -3 -I tests/sim/check_tables.py base base.damage.json
```
The other runs add `--library`, `--treatment`, `--model`, `--observed`, `--undamaged` and `--undamaged-length-mean` as in the table, with seeds 2 to 7. `damage` gets the same `--library` and `--treatment`.

An eighth run, **mt**, used the human mtDNA reference (rCRS, `--fasta`, with its N at position 3107), ds UDG-half, and `--divergence 0.001`. Its largest |z| was 3.54.

**Split by population.** Pooled rates look the same for "half the reads damaged" and "every read damaged at half the rate". So each file's damaged reads (`_D`) and undamaged reads (`_U`) were also run through `damage` separately, and checked with `--population`:
- Damaged reads match the damaged-only expectation: largest |z| 2.06–3.32 across the eight runs.
- Undamaged reads show sequencing errors (and divergence) only: largest |z| 1.77–3.91.

**CpG density.** The random contig (GC 0.42) has a CpG at 4.4% of positions; the rCRS has 2.6%, and the human nuclear genome about 1%. UDG-half and full runs on the random contig therefore score more CpG sites than real data would. The study's treated runs use the rCRS.

### 1.2 The oracle scorer reproduces pmdshodh's score

score.py with `--curve pmdtools` (P 0.3, C 0.01) was run against `pmdshodh score` with the same `--library` and `--treatment` on all eight files: 800,000 reads.
- Scores: largest difference 0.
- Site counts: no read differs.

The oracle therefore selects the same sites and computes the same terms as pmdshodh. On **base**, whose true curve is PMDtools' default (A = q = 0.3, C = 0.01), `--curve true` gives the same scores as `--curve pmdtools` (largest difference 0). That checks the true-curve path, including the swap of the 5′ and 3′ curves on reverse reads.

### 1.3 The ceiling scorer

`--curve true` alone keeps pmdshodh's null model (0.1% damage, π 0.001), which the simulated undamaged reads do not follow: they carry sequencing errors only (and divergence when it is set). A fitted curve could therefore beat it, so it is not a ceiling.

The ceiling is `--curve true --null exact --length`. This is the exact log-likelihood ratio of the simulator's own model on the sites pmdshodh scores, plus the read-length term.

An exact likelihood ratio s must satisfy E[exp(−s)] = 1 over damaged reads. Measured on each file (100,000 reads, about half of them damaged):

| Run | base | ss | half | full | obs | sshalf | sshalf2 | mt |
|---|---|---|---|---|---|---|---|---|
| Mean exp(−s), damaged reads | 1.007 | 1.009 | 0.997 | 1.004 | 1.000 | 0.974 | 1.001 | 1.000 |

The mirror check, mean exp(s) over undamaged reads, is heavy-tailed: a few undamaged reads with several sequencing errors at C sites carry very large ratios. Its values (0.86–2.14) are too noisy to check against.

Two limits of the ceiling:
- **It ignores divergence**, which the mt run has. Its mean of 1.000 shows the effect is small at 0.001.
- **It covers only the sites pmdshodh scores.** With ss UDG-half, the 3′ terminal C that is not at a CpG is damaged in the simulation but not scored, so it adds nothing.

The oracle in 1.3 had no divergence. In step 3 the exact model gained it: the sample carries the reference base, or with probability `--divergence` another one. Without it, an undamaged read that carries a diverged C/T looks like overwhelming evidence of damage. The ceiling then fell *below* the default on the first base-case run, which exposed the gap.

### 1.4 Open point for the observed-profile runs

Observed profiles come from pooled real libraries, so they already include their undamaged reads. The simulator applies such a profile to damaged reads only, and then adds undamaged reads, which dilutes it twice. Before the observed-profile runs (study step 4), the profile is either rescaled, or the runs are reported as testing its shape only. Decided 2026-11-24: shape only, scaled to a fixed terminal level (section 3).

## 2. Core runs (study step 3, 2026-11-24)

Each case has 400,000 simulated reads with divergence 0.001. The base case is:
- damage: A 0.30, q 0.30, C 0.01, geometric (PMDtools' default curve);
- undamaged share 0.5;
- reads: lognormal lengths, mean 45 bp and sd 15, the same for both populations;
- library: ds, no UDG, on a random contig.

Each other case changes one thing. UDG-half runs use the rCRS (1.1). The base case was run with three seeds.

The scorers (docs/damage-curve-study.md section 5) are:
- **default:** pmdshodh, P 0.3, C 0.01;
- **fitted:** pmdshodh with P = q, where q is fitted to the run's own damage tables by `tests/sim/fit.py`;
- **true curve:** the simulator's curve, with pmdshodh's null model;
- **oracle:** the exact likelihood ratio, the ceiling;
- **+length:** the oracle plus the read-length term.

The measure is the share of damaged reads kept when each scorer removes the same share of undamaged reads as the default does at threshold 3. The default removes 99.82–99.96% of undamaged reads at threshold 3, so this compares the scorers where filtering is actually used.

| Case | Default kept | Fitted | True curve | Oracle | +length | Fitted q |
|---|---|---|---|---|---|---|
| base, seed 11 | 31.1% | +0.0 | 0 | +0.7 | +0.7 | 0.331 |
| base, seed 12 | 30.9% | −0.3 | 0 | +1.2 | +1.2 | 0.328 |
| base, seed 13 | 30.6% | −0.4 | 0 | +0.9 | +0.9 | 0.328 |
| A 0.05 | 7.6% | −0.0 | +1.0 | +1.5 | +1.5 | 0.321 |
| A 0.15 | 17.4% | −0.1 | +0.6 | +1.2 | +1.2 | 0.325 |
| A 0.45 | 42.6% | −0.2 | −0.2 | +0.7 | +0.7 | 0.329 |
| q 0.20 | 36.9% | +0.7 | +0.9 | +2.4 | +2.4 | 0.236 |
| q 0.45 | 24.9% | −0.1 | +0.1 | +0.7 | +0.7 | 0.470 |
| q 0.60 | 20.5% | +0.4 | +0.5 | +0.8 | +0.8 | 0.613 |
| terminal excess +0.10 | 33.9% | −0.1 | +0.1 | +0.7 | +0.7 | 0.367 |
| undamaged 0.2 | 30.7% | +0.0 | 0 | +0.4 | +0.4 | 0.329 |
| undamaged 0.8 | 31.0% | −0.1 | 0 | +1.0 | +1.0 | 0.328 |
| long reads (mean 80 bp) | 31.0% | +0.0 | 0 | +1.8 | +1.8 | 0.302 |
| undamaged reads longer (mean 70 bp) | 30.8% | −0.2 | 0 | +0.4 | **+5.4** | 0.329 |
| ss | 30.2% | −0.2 | 0 | +1.6 | +1.6 | 0.341 |
| ds UDG-half, rCRS | 15.9% | −0.0 | 0 | +0.3 | +0.3 | 0.349 |
| ss UDG-half, rCRS | 9.9% | +0.0 | 0 | +0.7 | +0.7 | 0.348 |

The other columns are percentage points relative to the default. "0" means the true curve is the default curve.

How precise the measure is: across the three base seeds, the default ranges 30.6–31.1% and the oracle's gain 0.7–1.2 points. Differences under about 0.5 points are noise.

What the runs show:
1. **No change to the damage curve is worth making.** The ceiling for any change to the curve or the null model is the oracle. It gains 0.3–2.4 points over the default, at most 2.4 (shallow decay, q 0.2), against the agreed margin of 5 (S3). The true curve under pmdshodh's null gains at most 1.0 point.
2. **Fitting P gains nothing.** The fitted scorer is within ±0.7 points of the default in every case. It also recovers q only roughly (0.236 for 0.20, 0.613 for 0.60, about 0.33 for 0.30): the fit's background soaks up the floor C, as on real data (section 2.1 of the study).
3. **The read limits what any scorer can find, not the model.** At this removal rate even the oracle keeps only 9–43% of damaged reads. Most damaged reads carry too few damaged bases to be told from undamaged ones. The ROC areas of all scorers agree to within 0.01 in every case without a length difference.
4. **Read length is the only lever above the margin, and only when the populations differ in length.** With undamaged reads 25 bp longer on average, the length term adds 5.0 points over the oracle (+5.4 over the default). With equal length distributions it adds nothing. Its value in real data depends on how far endogenous and contaminant lengths differ, which this study does not know.

Commands (pmdshodh 0.3.3; WSL Ubuntu 24.04, Python 3.12):
```
PMDSHODH=pmdshodh tests/sim/run_case.sh base1 --seed 11 --reads 400000 --divergence 0.001
```
The other cases add `--amplitude`, `--decay`, `--model terminal`, `--undamaged`, `--length-mean 80 --length-sd 25`, `--undamaged-length-mean 70 --undamaged-length-sd 25`, `--library ss`, and `--treatment half --fasta rCRS.fa`, with seeds 12–82. `tests/sim/evaluate.py` gives each case's numbers.

## 3. Real damage shapes (study step 4, 2026-11-24)

The damaged reads take their damage from five real libraries' tables, by shape only. Each profile is scaled so the 5′ terminal base has 31% damage, the default curve's value (`--observed-terminal 0.31`). Everything else is as in the base case: 400,000 reads, divergence 0.001, undamaged share 0.5, mean length 45 bp.

The profile is the C→T rate (ds 3′: G→A) minus, at each position, the transition that damage does not cause there (`--observed-background transition`):
- ds: 5′ G→A and 3′ C→T;
- ss: G→A at both ends.

This also removes divergence from the reference. A first attempt subtracted only the transversion rate. It left BER11's interior G→A of about 1.7%, which is divergence from the horse reference, inside the "damage", and inflated BER11's ceiling to +9.2. That attempt is not reported further.

Single-stranded tables pool reads of every length, so a short read's 5′ table also holds damage from its 3′ end. `simulate.py` unmixes them using the library's own length histogram, solving for the per-end curves whose pooled expectation matches the tables. In the interior the two ends cannot be told apart: UKY001's interior damage moves wholly to the 3′ curve. Its total, about 7% C→T throughout the molecule against 0.12% G→A, is what matters for scoring.

| Profile (source) | Scaled shape, 5′ z = 1, 2, 3, 5, 10, 20, 30 | Default kept | Fitted | True curve | Oracle |
|---|---|---|---|---|---|
| BER11_L_NOuser (ds, no UDG) | 0.31, 0.21, 0.17, 0.14, 0.06, 0.04, 0.02 | 30.3% | +1.2 | +3.0 | +4.3 |
| JK2782 (ds, no UDG) | 0.31, 0.16, 0.14, 0.05, 0.05, 0.02, 0.00 | 29.6% | +0.1 | +1.8 | +2.3 |
| I8525 (ds, UDG-half, rCRS) | 0.31, 0.20, 0.08, 0.03, 0.01, 0.00, 0.00 | 14.5% | −0.0 | −0.0 | +0.1 |
| UKY001.B0102 MT (ss, no UDG) | unmixed 5′: 0.31, 0.15, 0.07, 0.06, 0.02, 0, 0; 3′: 0.26 at the end, 0.10 at z = 20–30 | 46.7% | +0.7 | +5.1 | **+6.2** |
| I19359 (ss, USER) | 0.31, then 0 (3′: 0.27, 0.17, then 0) | 14.7% | +0.1 | +0.3 | +0.4 |

The measure and the columns are as in section 2. The UKY001 shape is given in the scaled units (unmixed 0.223 at z = 1, times 0.31/0.223).

**Raising PMDtools' floor C** (`--PMDconstant`, an existing option) with P kept at 0.3:

| Profile | C 0.02 | C 0.035 | C 0.05 | C 0.07 | Oracle |
|---|---|---|---|---|---|
| UKY001 (ss) | +2.0 | +3.2 | +4.3 | +4.2 | +6.2 |
| BER11 (ds) | +1.6 | +1.7 | +1.7 | +0.9 | +4.3 |

What the real shapes show:
1. **Real shapes leave more room than the geometric curve,** because real profiles have long interior tails that PMDtools' C = 0.01 underweights.
2. **For ds libraries, the room stays under the margin.** The ceiling is +0.1 to +4.3 for the three ds profiles and +0.4 for the USER-treated ss one.
3. **One profile clears the margin: UKY001,** untreated single-stranded, at +6.2 (true curve +5.1). Its about 7% interior C→T is the reason. A higher floor alone, `--PMDconstant 0.05`, recovers +4.3 of it. Fitting P recovers +0.7.
4. **Fitting P stays within +1.2 everywhere.**

Caveats:
- Each shape is scaled to one terminal level.
- The ss unmixing is not unique in the interior.
- These are shapes in simulated reads, not real reads. Step 5, the labelled real-read mixture, can check the UKY001 result directly, since it is itself an untreated ss library.

Commands:
```
PMDSHODH=pmdshodh tests/sim/run_case.sh tUKY001 --seed 93 --reads 400000 --divergence 0.001 \
    --model observed --observed UKY001_MT.damage_ss.json --observed-terminal 0.31 \
    --observed-background transition --library ss
pmdshodh score --quiet --library ss --PMDconstant 0.05 uk.sam -o uk.c0.05.tsv
```
Seeds: BER11 91, JK2782 92, UKY001 93, I19359 94, I8525 95 (`--treatment half --fasta rCRS`).

## 4. Real reads with known labels (study step 5, 2026-11-24)

The mixture uses real ancient and real modern reads on the human mtDNA, so the label is each read's origin:
- **Ancient:** UKY001.B0102 MT, an untreated single-stranded library (ENA ERR3991018; 281,906 reads), used as it is.
- **Modern:** NA12878's mtDNA reads from 1000 Genomes phase 3. Only the `MT` region of `NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam` was fetched with `samtools view` (21.5 MB, 336,018 records, md5 `de47c7a8972efd40a74da6ef9ad1fdfd`). Its `MT` contig is the rCRS, the reference UKY001 was aligned to.

`tests/sim/mix_real.py` draws 100,000 reads of each:
- The modern reads with an all-match CIGAR are cut to one window each, with a length drawn from the ancient reads' own lengths.
- The two groups end up alike in length (mean 51.8 and 51.9 bp) and base quality (mean 38.9 and 38.4).
- MD tags are dropped, and pmdshodh takes reference bases from the rCRS (`-r`), in ss mode.

The fitted scorer's P is the decay fitted to UKY001's own damage table (q = 0.203). The measure is as before: ancient reads kept when each scorer removes the same share of modern reads as the default does at threshold 3. That share is 99.86% (seed 101) and 99.87% (seed 102).

| Scorer | Kept, seed 101 | Seed 102 |
|---|---|---|
| default (P 0.3, C 0.01) | 41.7% | 41.8% |
| fitted P = 0.203 | +1.0 | +1.2 |
| `--PMDconstant 0.02` | +1.4 | |
| `--PMDconstant 0.035` | +3.4 | +3.3 |
| `--PMDconstant 0.05` | +3.3 | +3.5 |
| `--PMDconstant 0.07` | +3.6 | +3.6 |
| `--PMDconstant 0.1` | +3.4 | |

The real reads agree with the simulated UKY001 shape (section 3: C 0.05 +4.3, fitted P +0.7). Raising the floor C to 0.035–0.07 keeps about 3.5 points more ancient reads at the same removal of modern reads. The gain is flat across that range. Fitting P keeps about 1 point more.

Caveats:
- The ancient library carries its own unknown share of modern contamination, labelled "ancient" here, so every kept rate understates. The comparisons between scorers stand.
- The modern windows come from anywhere in 101 bp reads, not from molecule ends.
- This is one library and one modern sample.

Commands:
```
samtools view -b -o NA12878.MT.bam https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam MT
samtools view UKY001.B0102.MT.bam > uky.sam; samtools view NA12878.MT.bam > na.sam
py -3 -I tests/sim/mix_real.py --ancient uky.sam --modern na.sam --contig MT --length 16569 --reads 100000 --seed 101 --out mix
pmdshodh score --quiet --library ss -r rCRS.fa --PMDconstant 0.05 mix.sam -o c0.05.tsv
py -3 -I tests/sim/evaluate.py default=default.tsv c0.05=c0.05.tsv
```
