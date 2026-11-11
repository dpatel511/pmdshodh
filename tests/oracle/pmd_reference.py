"""Reference model of the PMDtools 0.60 per-site score (readable oracle for pmdshodh).

Written from black-box tests of PMDtools 0.60 (docs/blackbox.md section 1), not from its
source. The paper is Skoglund et al. 2014, PNAS, doi:10.1073/pnas.1318934111.

Run directly to check this model against PMDtools' own outputs:
    py -3 -I tests/oracle/pmd_reference.py tests/data/pmdtools_site_terms.tsv
"""
import math
import sys

# Damage assumed by PMDtools' null model. The paper uses 0 here; PMDtools 0.60 uses 0.001,
# and no option changes it (docs/blackbox.md section 2).
NULL_BACKGROUND = 0.001


def error_rate(q):
    """Probability that a sequencing error gives this one mismatch type (phred Q, no cap)."""
    return 10 ** (-q / 10) / 3


def damage(z, p, c):
    """Damage probability at distance z (1 = terminal base) from the relevant end (paper Eq. 3)."""
    return (1 - p) ** (z - 1) * p + c


def p_match(d, pi, e):
    """Probability of a match at a C (or G) site with damage d (paper Eq. 1)."""
    return (1 - pi) * (1 - e) * (1 - d) + (1 - pi) * e * d + pi * e * (1 - d)


def p_match_two_ends(a, b, pi, e):
    """Match probability with damage from both ends, a and b (paper SI Eq. S1)."""
    return ((1 - pi) * (1 - e) * (1 - a) * (1 - b) + (1 - pi) * e * a * (1 - b)
            + (1 - pi) * e * b * (1 - a) + pi * e * (1 - a) * (1 - b))


def log_ratio(is_match, p_ancient, p_null):
    """Log-likelihood ratio for one observation, given the two models' match probabilities."""
    if is_match:
        return math.log(p_ancient / p_null)
    return math.log((1 - p_ancient) / (1 - p_null))


def site_term_ds(is_match, z, q, p=0.3, c=0.01, pi_ancient=0.001, pi_contamination=0.001):
    """Double-stranded score term of one informative site."""
    e = error_rate(q)
    return log_ratio(is_match,
                     p_match(damage(z, p, c), pi_ancient, e),
                     p_match(NULL_BACKGROUND, pi_contamination, e))


def site_term_ss(is_match, z, y, q, p=0.3, c=0.01, pi_ancient=0.001, pi_contamination=0.001):
    """Single-stranded (--ss) score term: PMDtools adds the ds term and the S1 term for the same site."""
    e = error_rate(q)
    two_ends = log_ratio(is_match,
                         p_match_two_ends(damage(z, p, c), damage(y, p, c), pi_ancient, e),
                         p_match_two_ends(NULL_BACKGROUND, NULL_BACKGROUND, pi_contamination, e))
    return site_term_ds(is_match, z, q, p, c, pi_ancient, pi_contamination) + two_ends


def main(path):
    worst = 0.0
    n = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("mode"):
                continue
            mode, z, y, q, is_match, p, c, pi_a, pi_c, expected = line.rstrip("\n").split("\t")
            args = dict(p=float(p), c=float(c), pi_ancient=float(pi_a), pi_contamination=float(pi_c))
            if mode == "ds":
                got = site_term_ds(is_match == "1", int(z), int(q), **args)
            else:
                got = site_term_ss(is_match == "1", int(z), int(y), int(q), **args)
            worst = max(worst, abs(got - float(expected)))
            n += 1
    print(f"rows {n}, largest |reference - PMDtools| = {worst:.3e}")


if __name__ == "__main__":
    main(sys.argv[1])
