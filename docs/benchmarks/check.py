#!/usr/bin/env python3
"""Correctness checks and timing summary for a run of run.sh.

Usage: python3 check.py WORK_DIR RESULTS_DIR INPUT [INPUT ...]
Exits nonzero if any parity or thread-identity check fails.
"""
import csv
import itertools
import statistics
import subprocess
import sys
from collections import defaultdict

TOL = 1e-6  # absolute tolerance for PMD score differences


def records(bam, fields=None):
    """Yield the SAM record lines of a BAM (no header), optionally cut to the first `fields` columns."""
    p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        line = line.rstrip("\n")
        yield "\t".join(line.split("\t")[:fields]) if fields else line
    if p.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam}")


def compare(a, b):
    """Return (n_a, n_b, n_differing_positions) for two line iterators."""
    na = nb = diff = 0
    for x, y in itertools.zip_longest(a, b):
        na += x is not None
        nb += y is not None
        diff += x != y
    return na, nb, diff


def lines(path):
    with open(path) as f:
        yield from (l.rstrip("\n") for l in f)


def main():
    work, results, inputs = sys.argv[1], sys.argv[2], sys.argv[3:]
    ok = True
    out = []

    def report(name, passed, detail):
        nonlocal ok
        ok &= passed
        out.append(f"{'PASS' if passed else 'FAIL'}\t{name}\t{detail}")

    for i in inputs:
        w = f"{work}/{i}"
        # 1. Kept reads: PMDtools vs pmdshodh --pmdtools-compat, whole records (fields 1-11) and in order.
        na, nb, d = compare(records(f"{w}.pmdtools.bam", 11), records(f"{w}.compat.t1.bam", 11))
        report(f"{i} kept reads PMDtools vs compat", na == nb and d == 0,
               f"PMDtools {na}, compat {nb}, differing positions {d}")
        na, nb, d = compare(records(f"{w}.pmdtools.bam"), records(f"{w}.compat.t1.bam"))
        out.append(f"INFO\t{i} full-line (incl. tags) differences PMDtools vs compat\t{d} of {na}")
        # Default mode vs PMDtools, by read identity (informational: known bug fixes).
        key = lambda bam: set(records(bam, 4))
        p, s = key(f"{w}.pmdtools.bam"), key(f"{w}.filter.t1.bam")
        out.append(f"INFO\t{i} default vs PMDtools kept reads\tonly PMDtools {len(p - s)}, only default {len(s - p)}, both {len(p & s)}")

        # 2. Per-read scores: PMDtools --printDS (last column, input order) vs pmdshodh score --pmdtools-compat.
        ds = (l.split()[-1] for l in lines(f"{w}.pmdtools.ds") if l.strip())
        tsv = lines(f"{w}.score_compat.t8.tsv")
        next(tsv)  # header
        n = over = 0
        maxd = 0.0
        count_mismatch = False
        for x, y in itertools.zip_longest(ds, tsv):
            if x is None or y is None:
                count_mismatch = True
                break
            dlt = abs(float(x) - float(y.split("\t")[2]))
            n += 1
            maxd = max(maxd, dlt)
            over += dlt > TOL
        report(f"{i} scores PMDtools vs compat", not count_mismatch and over == 0,
               f"reads compared {n}, largest |diff| {maxd:.3g}, above {TOL:g}: {over}"
               + (", LINE COUNTS DIFFER" if count_mismatch else ""))

        # 3. Thread identity: records (headers differ in @PG CL) and score tables.
        for kind, threads in (("filter", (2, 4, 8, 12)), ("compat", (8,))):
            for t in threads:
                na, nb, d = compare(records(f"{w}.{kind}.t1.bam"), records(f"{w}.{kind}.t{t}.bam"))
                report(f"{i} {kind} -@ {t} records == -@ 1", na == nb and d == 0, f"{na} vs {nb}, differing {d}")
        na, nb, d = compare(lines(f"{w}.score.t1.tsv"), lines(f"{w}.score.t8.tsv"))
        report(f"{i} score -@ 8 == -@ 1", na == nb and d == 0, f"{na} vs {nb} lines, differing {d}")
        na, nb, d = compare(lines(f"{w}.metadmg.t1.txt"), lines(f"{w}.metadmg.t8.txt"))
        out.append(f"INFO\t{i} metaDMG -@ 8 vs -@ 1 output\t{na} vs {nb} lines, differing {d}")

    # 4. Timing summary: median and range per input and configuration.
    rows = defaultdict(list)
    with open(f"{results}/times.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[(r["input"], r["config"])].append(r)
    table = ["input\tconfig\truns\twall_median_s\twall_min_s\twall_max_s\tcpu_median_s\tmax_rss_median_mb"]
    for (inp, cfg), rs in rows.items():
        wall = [float(r["wall_s"]) for r in rs]
        cpu = [float(r["user_s"]) + float(r["sys_s"]) for r in rs]
        rss = [int(r["max_rss_kb"]) / 1024 for r in rs]
        table.append(f"{inp}\t{cfg}\t{len(rs)}\t{statistics.median(wall):.1f}\t{min(wall):.1f}\t{max(wall):.1f}"
                     f"\t{statistics.median(cpu):.1f}\t{statistics.median(rss):.0f}")
    with open(f"{results}/summary.tsv", "w") as f:
        f.write("\n".join(table) + "\n")
    with open(f"{results}/checks.tsv", "w") as f:
        f.write("\n".join(out) + "\n")
    print("\n".join(out))
    print()
    print("\n".join(table))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
