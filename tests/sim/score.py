"""Score simulated reads with a chosen damage curve: the oracle scorers of the damage-curve study.

docs/damage-curve-study.md section 5. Uses pmdshodh's default score (no --pmdtools-compat): the
same sites, base qualities and null model, built on tests/oracle/pmd_reference.py. Only the
ancient model's damage curve is chosen here:

    --curve pmdtools   P (1 - P)^(z-1) + C, with --pmd-p and --pmd-c (default 0.3, 0.01)
    --curve true       the curve the simulator used for damaged reads (PREFIX.params.json)

--length adds the log-likelihood ratio of the read's length under the two populations'
length distributions. Writes read_name, pmds and n_sites, like `pmdshodh score`.

    py -3 -I tests/sim/score.py sim/base --curve true > sim/base.oracle.tsv
"""

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "oracle"))
from pmd_reference import NULL_BACKGROUND, error_rate, log_ratio, p_match, p_match_two_ends  # noqa: E402
from simulate import Damage, length_pmf  # noqa: E402

PI = 0.001  # pmdshodh's default polymorphism rate, in both models
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def read_fasta(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def ds_term(is_match, d, q):
    e = error_rate(q)
    return log_ratio(is_match, p_match(d, PI, e), p_match(NULL_BACKGROUND, PI, e))


def ss_term(is_match, a, b, q):
    e = error_rate(q)
    return log_ratio(
        is_match,
        p_match_two_ends(a, b, PI, e),
        p_match_two_ends(NULL_BACKGROUND, NULL_BACKGROUND, PI, e),
    )


DIVERGENCE = 0.0  # set from the simulation's parameters in main()


def exact_term(is_match, d, q):
    """The exact log-likelihood ratio of the simulator's own model at a reference C (or G).
    The sample carries the reference base, or with probability DIVERGENCE one of the other three;
    a sample C is damaged with probability d (undamaged reads: 0); then a sequencing error with
    probability 10^(-q/10) gives one of the other three bases. The read base is C or T (match or
    damage-type mismatch); CpG context changes at diverged sites are ignored."""
    eps, dv = 10 ** (-q / 10), DIVERGENCE

    def p_read(dmg):
        if is_match:
            return (1 - dv) * ((1 - dmg) * (1 - eps) + dmg * eps / 3) + dv * eps / 3
        return (1 - dv) * (dmg * (1 - eps) + (1 - dmg) * eps / 3) + dv / 3 * (1 - eps) + 2 * dv / 3 * eps / 3

    return math.log(p_read(d) / p_read(0.0))


def exact_ds_term(is_match, d, q):
    return exact_term(is_match, d, q)


def exact_ss_term(is_match, a, b, q):
    return exact_term(is_match, 1 - (1 - a) * (1 - b), q)


def score_read(seq, qual, ref, reverse, library, treatment, five, three, ds_term=ds_term, ss_term=ss_term):
    """src/sites.rs score_read, default mode, with damage curves `five(z)` and `three(z)` for
    the molecule's 5' and 3' ends. `seq`, `qual` and `ref` are in SAM orientation."""
    n = len(seq)
    total, sites = 0.0, 0
    if library == "ds":
        # In SAM orientation a forward read's C sites sit at the molecule's 5' end and its G
        # sites at the 3' end; a reverse read's are the other way round.
        c_curve, g_curve = (three, five) if reverse else (five, three)
        for i in range(n):
            r, b, q = ref[i], seq[i], qual[i]
            if r == "C" and b in "CT":
                g_next = i + 1 < n and ref[i + 1] == "G"
                if treatment == "half" and not (i == 0 or g_next):
                    continue
                if treatment == "full" and not g_next:
                    continue
                total += ds_term(b == "C", c_curve(i + 1), q)
                sites += 1
            elif r == "G" and b in "GA":
                c_prev = i > 0 and ref[i - 1] == "C"
                if treatment == "half" and not (i == n - 1 or c_prev):
                    continue
                if treatment == "full" and not c_prev:
                    continue
                total += ds_term(b == "G", g_curve(n - i), q)
                sites += 1
        return total, sites
    # Single-stranded: molecule orientation.
    if reverse:
        seq = "".join(COMPLEMENT[x] for x in reversed(seq))
        ref = "".join(COMPLEMENT[x] for x in reversed(ref))
        qual = qual[::-1]
    for i in range(n):
        r, b, q = ref[i], seq[i], qual[i]
        if r != "C" or b not in "CT":
            continue
        g_next = i + 1 < n and ref[i + 1] == "G"
        if treatment == "half" and not (i == 0 or g_next):
            continue
        if treatment == "full" and not g_next:
            continue
        total += ss_term(b == "C", five(i + 1), three(n - i), q)
        sites += 1
    return total, sites


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("prefix", help="the simulator's output prefix")
    p.add_argument("--curve", choices=["pmdtools", "true"], required=True)
    p.add_argument("--pmd-p", type=float, default=0.3)
    p.add_argument("--pmd-c", type=float, default=0.01)
    p.add_argument("--length", action="store_true", help="add the read-length log-likelihood ratio")
    p.add_argument(
        "--null",
        choices=["pmdshodh", "exact"],
        default="pmdshodh",
        help="pmdshodh: its null model and polymorphism rate (default); exact: the simulator's "
        "own model, the exact likelihood ratio and so the ceiling on these sites",
    )
    a = p.parse_args()
    params = json.load(open(f"{a.prefix}.params.json"))
    args = argparse.Namespace(**params)
    global DIVERGENCE
    DIVERGENCE = args.divergence
    if a.curve == "true":
        damage = Damage(args)
        five, three = damage.five_prime, damage.three_prime
    else:
        def five(z):
            return (1 - a.pmd_p) ** (z - 1) * a.pmd_p + a.pmd_c

        three = five
    if a.length:
        pmf_d = length_pmf(args.length_mean, args.length_sd, args.min_length, args.max_length)
        pmf_u = length_pmf(
            args.undamaged_length_mean, args.undamaged_length_sd, args.min_length, args.max_length
        )
    ref = read_fasta(f"{a.prefix}.fa")
    out = sys.stdout
    out.write("read_name\tpmds\tn_sites\n")
    with open(f"{a.prefix}.sam") as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            name, flag, pos, seq, qual = f[0], int(f[1]), int(f[3]), f[9], f[10]
            n = len(seq)
            score, sites = score_read(
                seq,
                [ord(c) - 33 for c in qual],
                ref[pos - 1 : pos - 1 + n],
                bool(flag & 16),
                args.library,
                args.treatment,
                five,
                three,
                *((exact_ds_term, exact_ss_term) if a.null == "exact" else ()),
            )
            if a.length:
                score += math.log(pmf_d[n] / pmf_u[n])
            out.write(f"{name}\t{score!r}\t{sites}\n")


if __name__ == "__main__":
    main()
