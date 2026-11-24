"""Check that pmdshodh's damage tables of a simulated file show the rates the simulator used.

docs/damage-curve-study.md section 4. For each position within the summary's range, the
observed C>T rate (and for ds the 3' G>A rate) is compared with its expectation, which mixes
the damaged and undamaged reads by their length distributions and adds sequencing errors
from the base qualities. Exits 1 if any position is more than 4 standard errors off.

    pmdshodh damage --quiet sim/base.sam -o sim/base.damage.json
    py -3 -I tests/sim/check_tables.py sim/base sim/base.damage.json
"""

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulate import Damage, length_pmf  # noqa: E402

BASES = "ACGT"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("prefix", help="the simulator's output prefix")
    p.add_argument("summary", help="pmdshodh damage summary of PREFIX.sam")
    p.add_argument("--max-z", type=float, default=4.0)
    a = p.parse_args()
    params = json.load(open(f"{a.prefix}.params.json"))
    summary = json.load(open(a.summary))
    args = argparse.Namespace(**params)
    damage = Damage(args)
    ss = args.library == "ss"
    quals = range(args.qual_min, args.qual_max + 1)
    eps = sum(10 ** (-q / 10) for q in quals) / len(quals)
    pmf_d = length_pmf(args.length_mean, args.length_sd, args.min_length, args.max_length)
    pmf_u = length_pmf(args.undamaged_length_mean, args.undamaged_length_sd, args.min_length, args.max_length)
    u = args.undamaged

    def site_damage(end, z, n):
        """True damage probability at distance z from `end` in a read of length n."""
        y = n - z + 1
        five, three = (z, y) if end == "five_prime" else (y, z)
        if ss:
            return 1 - (1 - damage.five_prime(five)) * (1 - damage.three_prime(three))
        return damage.five_prime(five) if end == "five_prime" else damage.three_prime(three)

    # Base composition of the reference, in molecule orientation (both strands): the share of
    # bases that are C, and of positions that start a CpG.
    ref = "".join(line.strip() for line in open(f"{a.prefix}.fa") if not line.startswith(">"))
    p_c = (ref.count("C") + ref.count("G")) / (2 * len(ref))
    p_cpg = sum(1 for i in range(len(ref) - 1) if ref[i : i + 2] == "CG") / (len(ref) - 1)

    def counted(end, z, n):
        """Expected number of sites pmdshodh counts at distance z from `end` in a read of
        length n (src/damage.rs): every C (ds: G at the 3' end) without treatment; with
        UDG-half the terminal base plus CpG sites; with full UDG CpG sites only."""
        if args.treatment == "none":
            return p_c
        i = z - 1 if end == "five_prime" else n - z  # molecule position, from the 5' end
        if ss or end == "five_prime":  # a C: terminal at i == 0, CpG needs a next base
            terminal, cpg_possible = i == 0, i + 1 < n
        else:  # a ds G read from the 3' end: terminal at i == n - 1, CpG needs a previous base
            terminal, cpg_possible = i == n - 1, i > 0
        if args.treatment == "half" and terminal:
            return p_c
        return p_cpg if cpg_possible else 0.0

    def expected(end, z):
        num = den = 0.0
        for n in range(max(z, args.min_length), args.max_length + 1):
            w = counted(end, z, n)
            wd, wu = w * (1 - u) * pmf_d[n], w * u * pmf_u[n]
            d = site_damage(end, z, n)
            num += wd * (d * (1 - eps) + (1 - d) * eps / 3) + wu * eps / 3
            den += wd + wu
        return num / den

    tables = {"five_prime": ("C", "T"), "three_prime": ("C", "T") if ss else ("G", "A")}
    worst = 0.0
    print("end\tz\tsites\tobserved\texpected\tz_score")
    for end, (ref, alt) in tables.items():
        t = summary["damage"][end][ref]
        for i in range(summary["damage"]["range"]):
            n = sum(t[b][i] for b in BASES)
            if n == 0:
                continue
            obs, exp = t[alt][i] / n, expected(end, i + 1)
            zs = (obs - exp) / math.sqrt(exp * (1 - exp) / n)
            worst = max(worst, abs(zs))
            print(f"{end}\t{i + 1}\t{n}\t{obs:.5f}\t{exp:.5f}\t{zs:+.2f}")
    print(f"largest |z| = {worst:.2f}")
    sys.exit(0 if worst <= a.max_z else 1)


if __name__ == "__main__":
    main()
