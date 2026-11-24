"""Simulate ancient (damaged) and modern (undamaged) reads, for the damage-curve study.

docs/damage-curve-study.md section 4. Writes, for an output prefix P:
    P.fa         the reference (a random contig, or the first contig of --fasta)
    P.sam        reads with MD tags, all-match CIGARs, sorted by nothing in particular
    P.truth.tsv  each read's population, length, strand and number of damaged bases
    P.params.json  every setting, so the scorers and checks know the true curve

Damage is drawn in molecule orientation, then reverse-strand reads are written
reverse-complemented, as an aligner would. Sequencing errors follow the base qualities.
The run is seeded: the same arguments give the same files.

    py -3 -I tests/sim/simulate.py --out sim/base --reads 100000 --seed 1
"""

import argparse
import json
import math
import random
import sys

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def observed_curves(path, library):
    """Per-position damage rates from a pmdshodh summary's damage tables, minus the
    transversion rate (C>A plus C>G at the 5' end) as a stand-in for errors. Returns the
    5' curve and the 3' curve as lists indexed from 0 (the terminal base)."""
    s = json.load(open(path))
    d = s["damage"]
    r = d["range"]

    def rates(end, ref, alt):
        t = d[end][ref]
        out = []
        for i in range(r):
            n = sum(t[b][i] for b in BASES)
            out.append(t[alt][i] / n if n else 0.0)
        return out

    c5 = d["five_prime"]["C"]
    n5 = sum(sum(c5[b][i] for b in BASES) for i in range(r))
    background = sum(c5["A"][i] + c5["G"][i] for i in range(r)) / n5 if n5 else 0.0
    five = rates("five_prime", "C", "T")
    three = rates("three_prime", "C", "T") if library == "ss" else rates("three_prime", "G", "A")
    return ([max(x - background, 0.0) for x in five], [max(x - background, 0.0) for x in three])


class Damage:
    """Damage probability at distance z (1 = terminal base) from the 5' or 3' end."""

    def __init__(self, args):
        self.model = args.model
        if self.model == "observed":
            self.five, self.three = observed_curves(args.observed, args.library)
        else:
            a, q, c = args.amplitude, args.decay, args.floor
            excess = args.terminal_excess if self.model == "terminal" else 0.0
            curve = [a * (1 - q) ** z + c for z in range(args.max_length)]
            curve[0] += excess
            self.five = self.three = [min(x, 1.0) for x in curve]

    @staticmethod
    def _at(curve, z):
        return curve[min(z, len(curve)) - 1]

    def five_prime(self, z):
        return self._at(self.five, z)

    def three_prime(self, z):
        return self._at(self.three, z)


def lognormal(rng, mean, sd):
    """A lognormal draw with the given mean and standard deviation (on the natural scale)."""
    sigma2 = math.log(1 + (sd / mean) ** 2)
    return rng.lognormvariate(math.log(mean) - sigma2 / 2, math.sqrt(sigma2))


def length_pmf(mean, sd, lo, hi):
    """The distribution `draw_length` draws from: a lognormal, rounded and truncated to lo..hi."""
    sigma2 = math.log(1 + (sd / mean) ** 2)
    mu, sigma = math.log(mean) - sigma2 / 2, math.sqrt(sigma2)

    def cdf(x):
        return 0.5 * (1 + math.erf((math.log(x) - mu) / (sigma * math.sqrt(2)))) if x > 0 else 0.0

    w = {n: cdf(n + 0.5) - cdf(n - 0.5) for n in range(lo, hi + 1)}
    total = sum(w.values())
    return {n: v / total for n, v in w.items()}


def draw_length(rng, mean, sd, lo, hi):
    while True:
        n = round(lognormal(rng, mean, sd))
        if lo <= n <= hi:
            return n


def damage_molecule(mol, damage, library, treatment, rng):
    """Apply damage to a molecule (5'->3'); return the damaged sequence and the event count."""
    n = len(mol)
    out = list(mol)
    events = 0
    for i, base in enumerate(mol):
        z, y = i + 1, n - i  # distance from the 5' and the 3' end
        if base == "C":
            cpg = i + 1 < n and mol[i + 1] == "G"
            if library == "ss":
                p = 1 - (1 - damage.five_prime(z)) * (1 - damage.three_prime(y))
                terminal = z == 1 or y == 1
            else:
                p = damage.five_prime(z)
                terminal = z == 1
            new = "T"
        elif base == "G" and library == "ds":
            cpg = i > 0 and mol[i - 1] == "C"
            p = damage.three_prime(y)
            terminal = y == 1
            new = "A"
        else:
            continue
        if treatment == "half" and not (terminal or cpg):
            continue
        if treatment == "full" and not cpg:
            continue
        if rng.random() < p:
            out[i] = new
            events += 1
    return "".join(out), events


def add_errors(seq, quals, rng):
    out = list(seq)
    for i, q in enumerate(quals):
        if rng.random() < 10 ** (-q / 10):
            out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def md_tag(read, ref):
    """MD for an all-match alignment of `read` against `ref` (same length, forward strand)."""
    parts, run = [], 0
    for r, b in zip(ref, read):
        if r == b:
            run += 1
        else:
            parts.append(f"{run}{r}")
            run = 0
    parts.append(str(run))
    return "".join(parts)


def reference(args, rng):
    if args.fasta:
        name, seq = None, []
        with open(args.fasta) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name is not None:
                        break
                    name = line[1:].split()[0]
                elif line:
                    seq.append(line.upper())
        seq = "".join(seq)
        if set(seq) - set(BASES + "N"):
            sys.exit("the reference must contain only A, C, G, T and N")
        return name, seq
    gc = args.gc / 2
    weights = [0.5 - gc, gc, gc, 0.5 - gc]
    return "sim", "".join(rng.choices(BASES, weights=weights, k=args.ref_length))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--reads", type=int, default=100_000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--fasta", help="use the first contig of this FASTA instead of a random one")
    p.add_argument("--ref-length", type=int, default=1_000_000)
    p.add_argument("--gc", type=float, default=0.42)
    p.add_argument("--library", choices=["ds", "ss"], default="ds")
    p.add_argument("--treatment", choices=["none", "half", "full"], default="none")
    p.add_argument("--model", choices=["geometric", "terminal", "observed"], default="geometric")
    p.add_argument("--amplitude", type=float, default=0.30, help="A: damage at the terminal base")
    p.add_argument("--decay", type=float, default=0.30, help="q: geometric decay")
    p.add_argument("--floor", type=float, default=0.01, help="C: damage at every position")
    p.add_argument("--terminal-excess", type=float, default=0.10)
    p.add_argument("--observed", help="pmdshodh summary whose damage tables give the curve")
    p.add_argument("--undamaged", type=float, default=0.5, help="share of undamaged reads")
    p.add_argument("--length-mean", type=float, default=45.0)
    p.add_argument("--length-sd", type=float, default=15.0)
    p.add_argument("--undamaged-length-mean", type=float, help="default: --length-mean")
    p.add_argument("--undamaged-length-sd", type=float, help="default: --length-sd")
    p.add_argument("--min-length", type=int, default=25)
    p.add_argument("--max-length", type=int, default=150)
    p.add_argument(
        "--divergence",
        type=float,
        default=0.0,
        help="share of reference positions where the sample differs (the same for every read)",
    )
    p.add_argument("--qual-min", type=int, default=20)
    p.add_argument("--qual-max", type=int, default=40)
    args = p.parse_args()
    if args.model == "observed" and not args.observed:
        p.error("--model observed needs --observed SUMMARY")
    if args.undamaged_length_mean is None:
        args.undamaged_length_mean = args.length_mean
    if args.undamaged_length_sd is None:
        args.undamaged_length_sd = args.length_sd

    rng = random.Random(args.seed)
    name, ref = reference(args, rng)
    if len(ref) < args.max_length:
        sys.exit("the reference is shorter than --max-length")
    damage = Damage(args)
    # The sample's genome: the reference with --divergence of its positions changed.
    sample = list(ref)
    if args.divergence > 0:
        for i, b in enumerate(sample):
            if b != "N" and rng.random() < args.divergence:
                sample[i] = rng.choice([x for x in BASES if x != b])
    sample = "".join(sample)

    with open(f"{args.out}.fa", "w") as fa:
        fa.write(f">{name}\n")
        for i in range(0, len(ref), 70):
            fa.write(ref[i : i + 70] + "\n")
    with open(f"{args.out}.params.json", "w") as pj:
        json.dump(vars(args) | {"contig": name}, pj, indent=1)

    with open(f"{args.out}.sam", "w") as sam, open(f"{args.out}.truth.tsv", "w") as truth:
        sam.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{name}\tLN:{len(ref)}\n")
        sam.write("@PG\tID:simulate\tPN:simulate.py\n")
        truth.write("read_name\tpopulation\tlength\tstrand\tdamage_events\n")
        for k in range(args.reads):
            damaged = rng.random() >= args.undamaged
            if damaged:
                n = draw_length(rng, args.length_mean, args.length_sd, args.min_length, args.max_length)
            else:
                n = draw_length(
                    rng, args.undamaged_length_mean, args.undamaged_length_sd, args.min_length, args.max_length
                )
            while True:  # fragments never span an N
                start = rng.randrange(len(ref) - n + 1)
                fragment = ref[start : start + n]
                if "N" not in fragment:
                    break
            source = sample[start : start + n]
            reverse = rng.random() < 0.5
            mol = revcomp(source) if reverse else source
            events = 0
            if damaged:
                mol, events = damage_molecule(mol, damage, args.library, args.treatment, rng)
            quals = [rng.randint(args.qual_min, args.qual_max) for _ in range(n)]
            mol = add_errors(mol, quals, rng)
            seq = revcomp(mol) if reverse else mol
            if reverse:
                quals.reverse()
            read = f"r{k}_{'D' if damaged else 'U'}"
            qual = "".join(chr(q + 33) for q in quals)
            nm = sum(1 for a, b in zip(seq, fragment) if a != b)
            sam.write(
                f"{read}\t{16 if reverse else 0}\t{name}\t{start + 1}\t37\t{n}M\t*\t0\t0\t{seq}\t{qual}"
                f"\tNM:i:{nm}\tMD:Z:{md_tag(seq, fragment)}\n"
            )
            truth.write(
                f"{read}\t{'damaged' if damaged else 'undamaged'}\t{n}\t{'-' if reverse else '+'}\t{events}\n"
            )


if __name__ == "__main__":
    main()
