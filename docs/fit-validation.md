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

### 1.2 The oracle scorer reproduces pmdshodh's score

score.py with `--curve pmdtools` (P 0.3, C 0.01) was run against `pmdshodh score` with the same `--library` and `--treatment` on all seven files: 700,000 reads.
- Scores: largest difference 0.
- Site counts: no read differs.

The oracle therefore selects the same sites and computes the same terms as pmdshodh. With `--curve true`, only the damage curve differs.
