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

Observed profiles come from pooled real libraries, so they already include their undamaged reads. The simulator applies such a profile to damaged reads only, and then adds undamaged reads, which dilutes it twice. Before the observed-profile runs (study step 4), the profile is either rescaled, or the runs are reported as testing its shape only.

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
