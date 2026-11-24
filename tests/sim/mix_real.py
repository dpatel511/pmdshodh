"""Mix real ancient and real modern reads with known labels: step 5 of the damage-curve study.

docs/damage-curve-study.md section 6. Reads two headerless SAM files (`samtools view`):
- ancient reads are kept as they are, single-end, named NAME_D;
- modern reads with an all-match CIGAR are cut to one window each, whose length is drawn from
  the ancient reads' own lengths (so length does not separate the two), named NAME_U.
Flags keep only the strand; MD and NM tags are dropped, so score with `-r` on the reference.
Both sets are put on one contig. Writes PREFIX.sam; seeded.

    py -3 -I tests/sim/mix_real.py --ancient a.sam --modern m.sam --contig MT --length 16569 \
        --reads 100000 --out mix
"""

import argparse
import random
import re

KEEP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800  # unmapped, secondary, QC fail, duplicate, supplementary


def records(path):
    with open(path) as f:
        for line in f:
            r = line.rstrip("\n").split("\t")
            if int(r[1]) & KEEP_FLAGS or r[2] == "*" or r[5] == "*":
                continue
            yield r


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--ancient", required=True)
    p.add_argument("--modern", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--length", type=int, required=True, help="contig length")
    p.add_argument("--reads", type=int, default=100_000, help="reads per population")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    rng = random.Random(a.seed)

    ancient = list(records(a.ancient))
    modern = [r for r in records(a.modern) if re.fullmatch(r"\d+M", r[5])]
    rng.shuffle(ancient)
    rng.shuffle(modern)
    lengths = [len(r[9]) for r in ancient]
    out = []
    for i, r in enumerate(ancient[: a.reads]):
        out.append([f"{r[0]}_{i}_D", str(int(r[1]) & 16), a.contig, r[3], r[4], r[5], r[9], r[10]])
    made = 0
    for r in modern:
        if made == a.reads:
            break
        n = len(r[9])
        k = rng.choice(lengths)
        while k > n:
            k = rng.choice(lengths)
        off = rng.randrange(n - k + 1)
        out.append(
            [f"{r[0]}_{made}_U", str(int(r[1]) & 16), a.contig, str(int(r[3]) + off), r[4], f"{k}M",
             r[9][off : off + k], r[10][off : off + k]]
        )
        made += 1
    if made < a.reads or len(ancient) < a.reads:
        raise SystemExit(f"not enough reads: {min(len(ancient), a.reads)} ancient, {made} modern")
    rng.shuffle(out)
    with open(f"{a.out}.sam", "w") as f:
        f.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{a.contig}\tLN:{a.length}\n@PG\tID:mix_real\tPN:mix_real.py\n")
        for name, flag, rname, pos, mapq, cigar, seq, qual in out:
            f.write(f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")


if __name__ == "__main__":
    main()
