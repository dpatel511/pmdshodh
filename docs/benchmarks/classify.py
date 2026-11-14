#!/usr/bin/env python3
"""Classify reads whose default pmdshodh score differs from PMDtools 0.60's, by documented cause.

Usage: python3 classify.py WORK_DIR DATA_DIR INPUT [TOLERANCE]
Reads, all in input order: WORK_DIR/INPUT.score.t1.tsv (pmdshodh score, default mode, from run.sh),
WORK_DIR/INPUT.pmdtools.ds (PMDtools --printDS, score in the last column, from run.sh) and
DATA_DIR/INPUT.bam (needs samtools). TOLERANCE defaults to 1e-6, check.py's parity tolerance.
Causes are the departures in docs/design.md section 5; any other difference is UNEXPLAINED.
"""
import collections, itertools, re, subprocess, sys

work, data, inp = sys.argv[1:4]
TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
d = open(f"{work}/{inp}.score.t1.tsv")
assert next(d).startswith("read_name\t")
ds = (l for l in open(f"{work}/{inp}.pmdtools.ds") if l.strip())
sam = subprocess.Popen(["samtools", "view", f"{data}/{inp}.bam"], stdout=subprocess.PIPE, text=True)
causes, crossed = collections.Counter(), collections.Counter()
only_pmdtools = only_default = both = 0
unexplained = []
n = differ = 0
for a, p, s in itertools.zip_longest(d, ds, sam.stdout):
    assert a is not None and p is not None and s is not None, f"streams differ in length at read {n + 1}"
    n += 1
    a, f = a.split("\t"), s.split("\t", 12)
    assert a[0] == f[0], f"read {n}: {a[0]} vs {f[0]}"
    x, y = float(a[2]), float(p.split()[-1])
    kd, kp = x >= 3, y >= 3
    both += kd and kp; only_default += kd and not kp; only_pmdtools += kp and not kd
    if abs(x - y) <= TOL:
        continue
    differ += 1
    cig = f[5]
    md = next((t[5:].rstrip("\n") for t in s.rstrip("\n").split("\t")[11:] if t.startswith("MD:Z:")), "")
    rev = int(f[1]) & 16
    k = []
    if re.search(r"[IS]", cig) and re.search(r"[D=]", cig): k.append("CIGAR reconstruction fix (I/S with D/=)")
    if re.search(r"[a-z]", md): k.append("lowercase MD")
    if rev and re.search(r"[^0-9^ACGTNRYacgt]", md): k.append("reverse-read shift")
    key = " + ".join(k) or "UNEXPLAINED"
    causes[key] += 1
    crossed[key] += kd != kp
    if not k and len(unexplained) < 10: unexplained.append((f[0], f[1], cig, md, x, y))
assert sam.wait() == 0
print(f"{inp}: reads {n}, default differs from PMDtools (> {TOL}) {differ}")
print(f"  kept at PMD >= 3: both {both}, only PMDtools {only_pmdtools}, only default {only_default}")
for k, v in causes.most_common():
    print(f"  {v:8d}  {k}  (crossing PMD 3: {crossed[k]})")
for u in unexplained:
    print("  unexplained:", u)
