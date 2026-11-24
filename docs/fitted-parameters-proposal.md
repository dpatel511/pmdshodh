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
- contamination estimates (roadmap principle);
- the polymorphism rates π and the null model's background damage (0.001), which stay as they are (F5);
- recalibrating the threshold (F7).

Fitting methods come from papers only: Skoglund et al. 2014 for the score, Briggs et al. 2007 for the damage pattern. No code from PMDtools, metaDMG or ngsBriggs is read.

## 2. The model now

A read's PMD score sums one log-likelihood ratio per informative site (src/model.rs). The ancient model's damage at distance z from the relevant end is (paper Eq. 3):

  d(z) = P (1 − P)^(z−1) + C,  with P = 0.3 and C = 0.01 by default.

One parameter, P, sets both the damage at the terminal base and how fast it decays. C is a floor at every position. Sequencing error comes from each base's quality and is not part of d(z).

## 3. What the data say: a prototype fit

Before choosing a design, two curves were fitted to the damage tables of eight real libraries: the seven in the table, and I19359 (finding 5). The fit is binomial maximum likelihood on the counts by position, positions 1–30:
- **ds:** 5′ C→T and 3′ G→A;
- **ss:** C→T at both ends.

A background rate b (sequencing error and divergence) is added to both curves. The prototype is a throwaway script and is not in the repository; the implemented fit replaces it (section 7).

- **M1, PMDtools' form:** f(z) = b + P (1 − P)^(z−1) + C.
- **M2, free amplitude:** f(z) = b + A (1 − q)^(z−1) + C. M1 is the case A = q = P.

M2 is also what PMDtools' form gives when only a fraction a of the reads is damaged: a·P (1 − P)^(z−1) + a·C is M2 with q = P and A = aP. So **A/q is the damaged fraction, if the damaged reads follow PMDtools' form.**

| Library | Type | Summary from | 5′ C→T z=1 | M1: P | M2: A, q | A/q | Deviance M1 / M2 (df 58 / 57) |
|---|---|---|---|---|---|---|---|
| BER11_L_NOuser | ds, no UDG | 0.1.0 | 18.4% | 0.093 | 0.108, 0.225 | 0.48 | 505,867 / 232,166 |
| JK2782 (eager test) | ds, no UDG | 0.2.4 | 13.9% | 0.102 | 0.116, 0.314 | 0.37 | 165 / 57 |
| JK2067 (eager test) | ds, presumed UDG-treated, fitted as untreated | 0.2.4 | 0.9% | **0.0001** | 0.009, 0.615 | 0.015 | 243 / 83 |
| Motala3 | ds, presumed UDG-treated, fitted as untreated | 0.1.0 | 1.4% | **0.0011** | 0.012, 0.493 | 0.024 | 346,920 / 3,558 |
| UKY001.B0102 MT | ss, no UDG | 0.2.4 | 28.0% | 0.145 | 0.145, 0.203 | 0.72 | 6,861 / 5,165 |
| I8525 MT | ds, UDG-half | 0.3.0 | 8.8% | 0.053 | 0.076, 0.434 | 0.17 | 26,256 / 364 |
| I0071 MT | ds, full UDG | 0.3.0 | 0.8% | **0.0001** | 0.007, 0.560 | 0.012 | 242 / 112 |

What this shows:
1. **PMDtools' two parameters alone do not fit the pooled profile.** Low-damage libraries get P ≈ 0: a flat curve with no damage at the ends, which would score every read near 0. Where damage is high, M1 misses the terminal base (BER11: 11.0% fitted against 18.4% observed).
2. **A free amplitude fits far better** (deviance 1.3 to 100 times lower), with decay rates q of 0.2–0.6.
   - This does not show that damaged reads break PMDtools' form: dilution by undamaged reads gives the same curve.
   - A/q is a plausible damaged fraction for BER11, JK2782, UKY001 and I8525 (0.17–0.72).
   - It is not plausible for JK2067, Motala3 and I0071 (0.01–0.02). There, UDG treatment rather than contamination lowers the damage.
   - The pooled profile cannot tell these two causes apart.
3. **Neither curve fits exactly.**
   - On the large libraries the deviance is far above its degrees of freedom, so real profiles have more structure than a geometric curve.
   - The terminal base in particular is often higher than M2 predicts (BER11 18.4% against 13.8%; UKY001 28.0% against 21.4%).
   - Binomial standard errors will therefore be too small. The fit quality must be reported, not only the estimates.
4. **Only C + b is identified, not C.** The table uses b from the transition that damage does not cause on the same strand (ds: 5′ G→A and 3′ C→T; ss: G→A).
   - Refitting with b = 0, or with b from transversions (5′ C→A and C→G), leaves A and q the same on five libraries. On the two lowest-damage libraries q moves by up to 0.08.
   - C, however, shifts by whatever b gives up. BER11's C is 0.016, 0.029 or 0.031 depending on the estimate.
   - The transition background is not clean either. In BER11, 5′ G→A rises from 0.26% at z = 1 to 1.7% at z = 20–30, which is the interior damage floor.
   - Without a separate background model, these tables cannot give C.
5. **The treatment matters.** I19359 (ss, USER-treated) shows damage only at the terminal base: 23.4% at z = 1 and 0.4% at z = 2. Fitted as untreated, both curves fail. The fit must use the sites the score uses, so it must follow `--treatment`.

## 4. Design

### 4.1 What is fitted and scored

Recommended (F1, F2):
- Fit A, q and a combined floor (C + b) to the pooled profile.
- Score with PMDtools' formula unchanged, with **P = q** and C at its default 0.01.

Why:
- **The formula does not change.** The scoring code stays as it is; only P changes, as if `--PMDpparam` had been set from the data.
- **q is what the data determine robustly** (finding 4).
- **q does not depend on contamination.** Dilution lowers A, not q.
- **C cannot be fitted from these tables** (finding 4).
- **A, A/q and the floor go in the summary** as diagnostics.

Caveat:
- P = q assumes that damaged reads have about q damage at the terminal base.
- In UDG-treated libraries, A/q is 0.01–0.02 because of the treatment. P = q may then overstate the damage of a typical read.
- The simulations of treated libraries test this before anything is released (section 5).

The alternative is to score with the fitted curve as it is: d(z) = A (1 − q)^(z−1) + C. That describes the sample's average read, so contamination dilutes it: with 30% contamination, A is 70% of the endogenous value. It also changes the formula, though A = q = P keeps the default bit-identical.

Neither option models contamination. A two-component fit per read would yield a contamination-like fraction, which the roadmap keeps out of scope.

The fit uses the summary's damage tables. These already follow `--library`, `--treatment` and `--min-baseq` (src/damage.rs), so the fit sees the same sites the score uses:
- **ds:** 5′ C→T and 3′ G→A, one shared curve.
- **ss:** C→T at both ends, one shared curve. The 5′ table mixes reads of all lengths, so in short reads it also carries some 3′ damage. Simulation will show whether this biases q. If it does, the fit gets its own table that counts only reads longer than twice the range.
- **UDG-half:** position 1 counts every C, later positions only CpG sites, as in scoring.

Numerics: binomial maximum likelihood with Nelder–Mead from fixed starting points, with A and q in (0, 1) and the floor in [0, 0.2). The counts are integer sums that do not depend on thread count. The same input therefore always gives the same fit and the same summary.

### 4.2 Where the fit runs

Input can be stdin, so the input cannot always be read twice. Recommended (F3):
- `pmdshodh damage` always fits and records the fit in its summary, or the reason it could not fit.
- `score` and `filter` take `--fit-from <summary.json>` and score with that summary's fitted P. They refuse a summary from another library type or treatment, or one without a usable fit.

A pipeline then runs `damage` and `filter` on the same file. The alternative, `--fit`, reads the file twice in one run and refuses stdin.

### 4.3 When the fit fails

pmdshodh stops with an error (F4) when:
- too few sites are counted at the terminal base (the threshold is set from the simulations);
- A is not clearly above 0, i.e. there is no damage signal;
- q's standard error is too large to use;
- the optimizer does not converge.

Falling back silently to the defaults would score with parameters the user did not ask for.

### 4.4 Outputs

- **Summary:** `format_version` 5.
  - `model` gains `"source": "default"` or `"fitted"`, and a `fit` object.
  - `fit` holds A, q, A/q, the floor, their standard errors, the deviance and its degrees of freedom, the number of sites, and the summary the fit came from.
  - A `damage` summary always carries `fit`.
- **BAM:** the `@PG` line already records the command, including `--fit-from`.
- **stderr:** one line with the fitted values when scoring with them.
- **`report`:** the fitted curve drawn over the observed damage.

### 4.5 What the scores mean

Fitted scores are not on PMDtools' default scale.
- With P = q, a library with steep decay (q around 0.5) gives more weight to terminal C→T and penalizes terminal matches more.
- A threshold of 3 therefore means something else with fitted values.
- The summary labels the scores as fitted (F7).

This is the trade-off of the milestone: scores that suit the sample, against comparability with published thresholds.

## 5. Validation

**Simulator** (tests/sim, Python like tests/oracle):
- reads from a reference, with a set length distribution, both strands and base qualities;
- sequencing errors drawn from those qualities;
- C→T damage (ds: also 3′ G→A) drawn from a chosen A, q and C;
- a fraction of undamaged reads, as contamination;
- UDG-half and full treatment, which remove damage except at the terminal base and CpG sites;
- output as SAM with MD tags.

The simulator must not only draw from the curve the fit assumes, since section 3 shows real profiles do not follow it exactly. It also has two misspecified cases:
- **terminal excess:** extra damage at the terminal base only;
- **observed profiles:** per-position damage rates taken from the real libraries in section 3.

**Parameter recovery** runs over a grid:
- A from 0.02 to 0.4, q from 0.2 to 0.7, C from 0 to 0.02;
- undamaged fractions from 0 to 0.7;
- library ds and ss, treatment none and half;
- 10⁴ to 10⁶ reads.

The fit must recover q within a stated tolerance. The tolerance and the minimum site count (section 4.3) are set from these runs before any real-data claim.

**Separation:** on simulated mixtures of damaged and undamaged reads, compare default and fitted scores. Measure the share of damaged reads kept and of undamaged reads removed, at several thresholds.
- This includes the misspecified cases and treated libraries (the caveat in 4.1).
- It is the only accuracy claim the milestone makes. It is measured on simulated data; no claim is made about real data.

**Real libraries:** the section 3 fits are rerun with the implemented fit and reported in docs.
- The roadmap's real ds and UDG-half criterion is met by BER11_L_NOuser (ds, no UDG) and I8525 (ds, UDG-half; Narasimhan et al. 2019).
- On real data the check is that fits converge and match the observed profile, not that filtering improves.

**Default unchanged:** all existing PMDtools parity tests pass with bit-identical scores.

## 6. Decisions needed (my recommendation first)

| # | Decision | Options |
|---|---|---|
| F1 | What the score uses | **PMDtools' formula with P = the fitted decay q** (formula unchanged; not affected by contamination) / the fitted curve A (1 − q)^(z−1) + C as it is (the average read, diluted by contamination; new formula) / PMDtools' P and C fitted directly (P → 0 on low-damage libraries, finding 1) |
| F2 | C | **Stays at 0.01; the fitted floor C + b is recorded only** (the tables cannot separate C from b, finding 4) / fit C with b from transversions (a lower bound for b, so C is an upper bound) |
| F3 | Where the fit runs | **`damage` records the fit; `score`/`filter --fit-from summary.json`** (works with pipes) / `--fit`, reading the file twice (no stdin) / both |
| F4 | Fit fails or finds no damage | **Error and stop** / warn and score with the defaults |
| F5 | Null model | **Unchanged (background 0.001, π 0.001)** / also fit the null model's background |
| F6 | Summary and MultiQC | **`format_version` 5; update MultiQC PR #3745 now so it reads newer versions (the fields it uses do not change) and shows the fitted values** / keep #3745 at versions 1–4 and send another PR after the release |
| F7 | Thresholds | **No recalibration; fitted scores are labelled in the summary and the docs say thresholds are not comparable** / also offer a threshold that keeps the same share of reads as the default |

## 7. Order of work

1. The simulator and its tests.
2. The fit as a function of damage counts (src/fit.rs), tested on simulated counts. If F1 changes the formula, the generalized curve goes in src/model.rs with a bit-identity test for the defaults.
3. Parameter recovery and separation on simulated data, including the misspecified cases, in docs/fit-validation.md.
   - **This is the go/no-go point.** If fitted scores do not separate damaged from undamaged reads better than the defaults, the milestone stops here and nothing user-facing changes.
4. CLI (`--fit-from`), summary version 5, docs, CHANGELOG.
5. The real-library fits, added to docs/fit-validation.md.
6. The MultiQC module update (F6) and a 0.4.0 release, with your OK.
