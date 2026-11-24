# Proposal: fitted damage parameters (roadmap milestone 6)

Status: **proposed 2026-11-24.** Decisions F1–F7 (section 6) are open. Nothing is implemented.

## 1. Goal and scope

Roadmap milestone 6 (docs/roadmap.md): estimate the damage parameters from the sample's own damage profile, then score reads with them. The fitted values and the fit quality are recorded in the summary. The default score stays PMDtools-equivalent; fitting is opt-in and labelled in every output.

Done when (roadmap): validated on simulated data with known parameters and on real ds and UDG-half samples, with the default unchanged.

In scope:
- the damage curve's parameters, fitted from the damage tables pmdshodh already counts;
- scoring with the fitted values (`score`, `filter`);
- a read simulator for validation.

Out of scope:
- contamination estimates (roadmap principle; see F2);
- the polymorphism rates π and the null model's background damage (0.001), which stay as they are (F5);
- recalibrating the threshold (F7).

Fitting methods come from papers only: Skoglund et al. 2014 for the score, Briggs et al. 2007 for the damage pattern. No code from PMDtools, metaDMG or ngsBriggs is read.

## 2. The model now

A read's PMD score sums one log-likelihood ratio per informative site (src/model.rs). The ancient model's damage at distance z from the relevant end is (paper Eq. 3):

  d(z) = P (1 − P)^(z−1) + C,  with P = 0.3 and C = 0.01 by default.

One parameter, P, sets both the damage at the terminal base and how fast it decays. C is a floor at every position. Sequencing error comes from each base's quality and is not part of d(z).

## 3. What the data say: a prototype fit

Before choosing a design, both candidate curves were fitted to the damage tables of eight real libraries. The fit is binomial maximum likelihood on the C→T counts by position (ds: 5′ C→T and 3′ G→A; ss: C→T at both ends), positions 1–30.

The background transition rate b (sequencing error and divergence) is estimated from the transition that damage does not cause on the same strand: 5′ G→A and 3′ C→T for ds, G→A at both ends for ss. Then:

- **M1, PMDtools' form:** f(z) = b + P (1 − P)^(z−1) + C (2 parameters).
- **M2, free amplitude:** f(z) = b + A (1 − q)^(z−1) + C (3 parameters). M1 is the case A = q = P.

| Library | Type | 5′ C→T z=1 | M1: P, C | M2: A, q, C | Deviance M1 / M2 (df 58 / 57) |
|---|---|---|---|---|---|
| BER11_L_NOuser | ds, no UDG | 18.4% | 0.093, 0.002 | 0.108, 0.225, 0.016 | 505,867 / 232,166 |
| JK2782 (eager test) | ds, no UDG | 13.9% | 0.102, 0 | 0.116, 0.314, 0.008 | 165 / 57 |
| JK2067 (eager test) | ds | 0.9% | **0.0001, 0** | 0.009, 0.615, 0 | 243 / 83 |
| Motala3 | ds | 1.4% | **0.0011, 0** | 0.012, 0.493, 0.0003 | 346,920 / 3,558 |
| UKY001.B0102 MT | ss, no UDG | 28.0% | 0.145, 0.059 | 0.145, 0.203, 0.067 | 6,861 / 5,165 |
| I8525 MT | ds, UDG-half | 8.8% | 0.053, 0 | 0.076, 0.434, 0.0001 | 26,256 / 364 |
| I0071 MT | ds, full UDG | 0.8% | **0.0001, 0** | 0.007, 0.560, 0 | 242 / 112 |

What this shows:
1. **PMDtools' form does not fit real profiles.** With one parameter for height and decay, a low-damage library (JK2067, Motala3, I0071) gets P ≈ 0: a flat curve with no damage at the ends. Scoring with that would give every read a score near 0. Where damage is high, M1 misses the terminal base (BER11: 11.0% fitted against 18.4% observed).
2. **A free amplitude fits much better** (deviance down 1.3- to 100-fold) and gives plausible decay rates (q 0.2–0.6).
3. **Neither fits exactly.** Deviances far above their degrees of freedom on the big libraries mean real profiles have more structure than a geometric curve. The terminal base in particular is often higher than M2 predicts (BER11 18.4% against 13.8%; UKY001 28.0% against 21.4%). Standard errors from a binomial fit will be too small; fit quality must be reported, not just the estimates.
4. **The treatment matters.** I19359 (ss, USER-treated) shows damage only at the terminal base (23.4% at z = 1, 0.4% at z = 2). Fitted as untreated, M1 and M2 both fail. The fit must use the same sites the score uses, i.e. follow `--treatment`.

## 4. Design

### 4.1 What is fitted

The ancient model's damage curve becomes d(z) = A (1 − q)^(z−1) + C (F1). With A = q = P this is exactly today's formula, so the default score is bit-identical and keeps the PMDtools tests passing.

The fit uses the summary's damage tables, which already follow `--library`, `--treatment` and `--min-baseq` (src/damage.rs). So the fit sees the same sites the score uses:
- **ds:** 5′ C→T and 3′ G→A, one shared curve.
- **ss:** C→T at both ends, one shared curve. The 5′ table mixes reads of all lengths, so in short reads it also carries some 3′ damage. Simulation will show whether this biases q. If it does, the fit gets its own table that counts only reads longer than twice the range.
- **UDG-half:** position 1 counts every C, later positions only CpG sites, as in scoring.

The background b is estimated as in section 3 and subtracted inside the fit. It is recorded but not used in the score, because the score already models sequencing error from base qualities.

Numerics: binomial maximum likelihood with Nelder–Mead from fixed starting points, A, q ∈ (0, 1), C ∈ [0, 0.2). The counts are integer sums that do not depend on thread count, so the same input always gives the same fit and the same summary.

### 4.2 What the fit describes

The damage tables pool every read, so contaminant reads dilute them. The fitted A is the damage of the sample's average read, not of its ancient molecules. With 30% contamination, A is about 70% of the endogenous value. Fitting the endogenous curve needs a two-component mixture, which yields a contamination-like fraction; the roadmap keeps contamination estimates out of scope (F2).

### 4.3 Where the fit runs

Input can be stdin, so the input cannot always be read twice. Recommended (F3):
- `pmdshodh damage` always fits and records the fit in its summary (or the reason it could not).
- `score` and `filter` take `--fit-from <summary.json>` and score with that summary's fitted values. They refuse a summary from another library type or treatment, or one without a usable fit.

A pipeline then runs `damage` and `filter` on the same file. The alternative, `--fit`, reads a file twice in one run and refuses stdin.

### 4.4 When the fit fails

pmdshodh stops with an error (F4) when:
- too few sites are counted at the terminal base (the threshold is set from the simulations);
- A is not clearly above 0, i.e. there is no damage signal;
- the optimizer does not converge.

Scoring a library with no damage signal would give scores near 0 for every read. Falling back silently to the defaults would score with parameters the user did not ask for.

### 4.5 Outputs

- **Summary:** `format_version` 5. `model` gains `"source": "default"` or `"fitted"`, and a `fit` object: A, q, C, b, their standard errors, the deviance and its degrees of freedom, the number of sites, and the summary the fit came from. A `damage` summary always carries `fit`.
- **BAM:** the `@PG` line already records the command, including `--fit-from`.
- **stderr:** one line with the fitted values when scoring with them.
- **`report`:** the fitted curve drawn over the observed damage.

### 4.6 What the scores mean

Fitted scores are not on PMDtools' scale. For a low-damage library the ancient model moves towards the null, so every score shrinks; for a high-damage library they grow. A threshold of 3 means something else with fitted values, and the summary says the scores are fitted (F7). This is the trade-off of the milestone: scores that suit the sample, against comparability with published thresholds.

## 5. Validation

**Simulator** (tests/sim, Python like tests/oracle): reads from a reference with a set length distribution, both strands, base qualities, sequencing errors from those qualities, and C→T (ds: also 3′ G→A) damage drawn from a chosen A, q and C. A fraction of the reads is undamaged, as contamination. Output is SAM with MD tags. UDG-half and full treatment remove damage except at the terminal base and CpG sites.

**Parameter recovery:** over a grid of A (0.02–0.4), q (0.2–0.7), C (0–0.02), library (ds, ss), treatment (none, half), and 10⁴–10⁶ reads, the fit recovers A (scaled by the undamaged fraction), q and C within a stated tolerance. The tolerance and the minimum site count (section 4.4) are set from these runs before any real-data claim.

**Separation:** on simulated mixtures of damaged and undamaged reads, compare default and fitted scores: the share of damaged reads kept and of undamaged reads removed at several thresholds. This is the only accuracy claim the milestone makes. It is measured on simulated data, and no claim is made about real data.

**Real libraries:** the fits for the libraries in section 3, rerun with the implemented fit, are reported in docs. The roadmap's real ds and UDG-half criterion is met by BER11_L_NOuser (ds, no UDG) and I8525 (ds, UDG-half; Narasimhan et al. 2019). For real data the check is that fits converge and match the observed profile, not that filtering improves.

**Default unchanged:** all existing PMDtools parity tests pass with bit-identical scores.

## 6. Decisions needed (my recommendation first)

| # | Decision | Options |
|---|---|---|
| F1 | Damage curve | **A (1 − q)^(z−1) + C: three parameters, PMDtools' form when A = q = P** / PMDtools' P and C only (section 3: P → 0 on low-damage libraries) / three parameters plus a separate terminal-base term |
| F2 | What the fit describes | **The sample's average read, documented as diluted by contamination** / a two-component mixture fit (gives a contamination-like fraction, out of scope unless agreed) |
| F3 | Where the fit runs | **`damage` records the fit; `score`/`filter --fit-from summary.json`** (works with pipes) / `--fit`, reading the file twice (no stdin) / both |
| F4 | Fit fails or finds no damage | **Error and stop** / warn and score with the defaults |
| F5 | Null model | **Unchanged (background 0.001, π 0.001); b is recorded only** / also fit the null model's background |
| F6 | Summary and MultiQC | **`format_version` 5; update MultiQC PR #3745 now so it reads newer versions (the fields it uses do not change) and shows A, q and C** / keep #3745 at versions 1–4 and send another PR after the release |
| F7 | Thresholds | **No recalibration; fitted scores are labelled in the summary and the docs say thresholds are not comparable** / also offer a threshold that keeps the same share of reads as the default |

## 7. Order of work

1. The simulator and its tests.
2. The fit as a function of damage counts (src/fit.rs), tested on simulated counts; the generalized curve in src/model.rs with a bit-identity test for the defaults.
3. Parameter-recovery runs; set the tolerance and the minimum site count.
4. CLI (`--fit-from`), summary version 5, docs, CHANGELOG.
5. Separation runs and the real-library fits, in docs/fit-validation.md.
6. The MultiQC module update (F6) and a 0.4.0 release, with your OK.
