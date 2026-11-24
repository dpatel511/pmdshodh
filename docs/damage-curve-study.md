# Proposal: measure before fitting (roadmap milestone 6, re-scoped)

Status: **agreed 2026-11-24.** Milestone 6 is re-scoped from "fitted damage parameters" to a measurement study (section 1). S1–S5 (section 8) were accepted as recommended.

## 1. Why the milestone was re-scoped

Roadmap milestone 6 was to estimate the damage parameters from the sample's own profile and score reads with them. A prototype fit on eight real libraries (section 2) found:
- PMDtools' two parameters fit the pooled damage profile badly.
- Only the decay rate is robust.
- Scoring with P set to that decay changes scores by about 0.3–0.4 per read.

So the open question comes before any fitting: **how much better can any per-read damage score separate damaged from undamaged reads than PMDtools' defaults already do?**

This milestone answers that question with simulated and real reads, and the answer decides what follows (section 7):
- the small fitted-P option;
- a better per-read model;
- or closing the milestone.

The study changes no CLI option, score or summary field.

## 2. Evidence so far

### 2.1 Prototype fit

Two curves were fitted to the damage tables of eight real libraries: the seven in the table, and I19359 (point 5). The fit is binomial maximum likelihood by position, positions 1–30:
- **ds:** 5′ C→T and 3′ G→A;
- **ss:** C→T at both ends.

Both curves include a background rate b. The prototype was a throwaway script and is not in the repository.

- **M1, PMDtools' form:** f(z) = b + P (1 − P)^(z−1) + C.
- **M2, free amplitude:** f(z) = b + A (1 − q)^(z−1) + C. M1 is the case A = q = P.

M2 is also what PMDtools' form gives when only a fraction a of the reads is damaged (q = P, A = aP). So A/q is the damaged fraction, if damaged reads follow PMDtools' form.

| Library | Type | Summary from | 5′ C→T z=1 | M1: P | M2: A, q | A/q | Deviance M1 / M2 (df 58 / 57) |
|---|---|---|---|---|---|---|---|
| BER11_L_NOuser | horse, ds, no UDG | 0.1.0 | 18.4% | 0.093 | 0.108, 0.225 | 0.48 | 505,867 / 232,166 |
| JK2782 (eager test) | mammoth, ds, no UDG | 0.2.4 | 13.9% | 0.102 | 0.116, 0.314 | 0.37 | 165 / 57 |
| JK2067 (eager test) | ds, presumed UDG-treated, fitted as untreated | 0.2.4 | 0.9% | **0.0001** | 0.009, 0.615 | 0.015 | 243 / 83 |
| Motala3 | ds, presumed UDG-treated, fitted as untreated | 0.1.0 | 1.4% | **0.0011** | 0.012, 0.493 | 0.024 | 346,920 / 3,558 |
| UKY001.B0102 MT | ss, no UDG | 0.2.4 | 28.0% | 0.145 | 0.145, 0.203 | 0.72 | 6,861 / 5,165 |
| I8525 MT | ds, UDG-half | 0.3.0 | 8.8% | 0.053 | 0.076, 0.434 | 0.17 | 26,256 / 364 |
| I0071 MT | ds, full UDG | 0.3.0 | 0.8% | **0.0001** | 0.007, 0.560 | 0.012 | 242 / 112 |

1. **PMDtools' P and C fitted directly fail.** Low-damage libraries get P ≈ 0, a flat curve that would score every read near 0. Where damage is high, M1 misses the terminal base (BER11: 11.0% fitted against 18.4% observed).
2. **A free amplitude fits far better** (deviance 1.3 to 100 times lower; q 0.2–0.6). It cannot tell dilution by undamaged reads from a different damage shape. A/q is plausible as a damaged fraction for the untreated libraries (0.17–0.72), but not for the UDG-treated ones (0.01–0.02), where treatment lowers the damage.
3. **Neither fits exactly.** Real profiles have more structure than a geometric curve, above all at the terminal base (BER11 18.4% against 13.8% from M2; UKY001 28.0% against 21.4%).
4. **Only the decay q is robust.**
   - With b = 0, b from transitions or b from transversions, A and q stay the same on five libraries, and q moves by up to 0.08 on the two lowest-damage ones.
   - C shifts by whatever b gives up (BER11: 0.016, 0.029 or 0.031). The tables cannot separate C from the background.
5. **The treatment matters.** I19359 (ss, USER-treated) has damage only at the terminal base: 23.4% at z = 1, 0.4% at z = 2.

### 2.2 What fitted values would do to scores

The table shows per-site terms at base quality 30 and the scores of two made-up reads:
- **Read A:** C→T at positions 1 and 3, plus 6 matching C sites.
- **Read B:** 6 matching C sites and no C→T.

| Curve | Damage at z=1 | C→T at z=1 | C→T at z=5 | C match at z=1 | Read A | Read B |
|---|---|---|---|---|---|---|
| Default (P 0.3, C 0.01) | 0.310 | 4.89 | 3.58 | −0.37 | 8.64 | −0.68 |
| BER11, P = q = 0.225 | 0.235 | 4.62 | 3.68 | −0.27 | 8.30 | −0.60 |
| BER11, fitted curve as is | 0.124 | 3.98 | 3.18 | −0.13 | 7.27 | −0.34 |
| JK2067, P = q = 0.615 | 0.625 | 5.59 | 2.36 | −0.98 | 9.00 | −1.14 |
| JK2067, fitted curve as is | 0.009 | 1.49 | **−0.42** | −0.01 | **1.62** | 0.00 |

What this shows:
- **P = q moves scores little.** Most of a C→T's weight comes from comparing any real damage with the null model's 0.1%.
- **Scoring with the pooled curve as it is breaks on low-damage libraries.** It drops below the null's 0.1% after a few bases, so a C→T at position 5 counts *against* damage.

## 3. The question and the measures

For each read the study knows whether it came from a damaged (ancient) or an undamaged (modern) population. A read from the damaged population counts as damaged even if, by chance, it carries no C→T.

Each scorer (section 5) is measured on the same reads:
- **At fixed thresholds** (1 to 5): the share of damaged reads kept and the share of undamaged reads removed.
- **At equal removal:** the share of damaged reads kept when each scorer removes the same share of undamaged reads as the default does at threshold 3. This is the main comparison. It does not depend on each scorer's scale, so a scorer whose scores run higher or lower is not favoured.
- **ROC AUC** over all thresholds.

## 4. Simulator

tests/sim, in Python like tests/oracle. It writes SAM with MD tags, and a truth file with each read's population and its damaged positions. It is seeded, so a run can be repeated exactly.

Each read gets:
- **Origin:** a fragment of a reference (the rCRS mtDNA and a random-sequence contig), with a length from a chosen distribution, on either strand.
- **Base qualities:** from a chosen distribution, and sequencing errors drawn from them.
- **Damage**, for damaged reads only, from the chosen damage model:
  - **geometric:** A (1 − q)^(z−1) + C;
  - **terminal excess:** geometric plus extra damage at the terminal base;
  - **observed:** per-position rates taken from a real library's table (section 2.1).
- **Library type:** ds (5′ C→T, 3′ G→A) or ss (C→T at both ends).
- **Treatment:** none, UDG-half (damage kept at the terminal base and at CpG) or full (CpG only).

The undamaged reads get errors only. A test checks that pmdshodh's damage table of a simulated file recovers the rates the simulator used.

## 5. Scorers

All scorers use PMDtools' per-read log-likelihood ratio, with the same sites, base qualities and null model. Only the ancient model's damage curve differs:
- **Default:** P = 0.3, C = 0.01 (pmdshodh's default score).
- **Fitted:** P = q from the fit to the run's pooled damage table, C = 0.01 (the earlier plan).
- **True curve:** the curve the simulator used for damaged reads, including terminal excess, with pmdshodh's null model. This is the best any change to the damage curve alone can do.
- **Oracle:** the exact likelihood ratio of the simulator's own model on the same sites: the true curve, and a null with no damage and no polymorphism. No scorer of these sites can beat it, so it is the ceiling (docs/fit-validation.md 1.3).
- **Oracle with length:** the oracle plus the log-likelihood ratio of the read's length under the two populations' length distributions. It shows what read length adds.

The gap between the true curve and the oracle measures what the null model costs. Real data have divergence and polymorphism, which the null model is there for, so the core runs use `--divergence 0.001`.

Default and fitted scores come from pmdshodh itself (`--PMDpparam`). The oracle scorers are Python, built on tests/oracle/pmd_reference.py, which pmdshodh's own parity tests already check.

## 6. Conditions

The grid of simulated conditions:
- **Damage at the terminal base** (A): 0.05, 0.15, 0.30, 0.45.
- **Decay** (q): 0.2, 0.3, 0.45, 0.6.
- **Damage model:** geometric, terminal excess, observed (BER11, JK2782, UKY001, I8525).
- **Undamaged share:** 0, 0.2, 0.5, 0.8.
- **Read length:** short (mean 45 bp) and long (mean 80 bp); damaged and undamaged reads have the same length distribution, except in the length runs.
- **Library and treatment:** ds and ss, none and UDG-half.
- **Size:** 10⁵ reads per run.

The full grid is large, so the core runs vary one factor at a time from a base case: A 0.30, q 0.3, geometric, undamaged share 0.5, short reads, ds, untreated. A second set is run for each real library's observed profile.

**Real reads with known labels:**
- Ancient reads from an untreated ds human mtDNA library, mixed in set proportions with modern human mtDNA reads.
- The modern reads are trimmed to the ancient length distribution, so length does not separate them.
- The ancient library carries its own unknown contamination, so its labels are noisy. The real-read results are reported as such, next to the simulated ones.
- The datasets are found on ENA and 1000 Genomes and proposed to you before any download (S4).

## 7. What the results decide

- **The oracle keeps few more damaged reads than the default** at equal removal (below the margin in S3) in the conditions that resemble the real libraries. Then no change to the damage curve is worth making: close the milestone, record the result in docs, and keep the default.
- **The fitted scorer gets most of the oracle's gain.** Then implement fitted P (`damage` records the fit; `score`/`filter --fit-from`), as in the earlier plan.
- **The oracle gains, but the fitted scorer does not.** The gain then comes from the curve's shape. Propose a model change, such as a terminal-base term, as a new design note for your decision.
- **Length adds much.** Propose length as evidence, as a new design note for your decision.

Every result is reported, including negative ones, in docs/fit-validation.md, with the commands, seeds and versions.

## 8. Decisions needed (my recommendation first)

| # | Decision | Options |
|---|---|---|
| S1 | Re-scope milestone 6 as this study | **Yes; the roadmap row changes to "measure what a better damage curve could gain", and fitting follows only if section 7 says so** / keep the fitting milestone and run the study inside it |
| S2 | Simulator and oracle language | **Python in tests/sim** (readable, beside tests/oracle; only the oracle scores in Python, pmdshodh scores the rest) / Rust, as a test-only binary |
| S3 | Margin for a worthwhile gain | **5 percentage points more damaged reads kept at equal removal**, in the conditions resembling the real libraries / 2 points / no fixed margin, judged from the results |
| S4 | Real reads with known labels | **Yes, with datasets proposed to you first** / simulated data only |
| S5 | Where results go | **docs/fit-validation.md, linked from the README's validation section only if a user-facing change follows** / docs only, never linked |

## 9. Order of work

1. The simulator, with its truth file and the table-recovery test.
2. The oracle scorers, checked against pmdshodh on the default curve (they must agree to 1e-9).
3. The core runs, with the base case first, then one factor at a time.
4. The observed-profile runs.
5. The real-read mixture, after S4.
6. docs/fit-validation.md and a proposal for what follows (section 7), for your decision.
