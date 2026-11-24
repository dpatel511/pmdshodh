"""Compare scorers on simulated reads: the measures of the damage-curve study.

docs/damage-curve-study.md section 3. Each scorer is a TSV with read_name and pmds columns
(`pmdshodh score` or score.py). A read is damaged if its name ends in _D. A read is kept when
its score is at least the threshold, as `pmdshodh filter` keeps it.

Measures, per scorer:
- kept_damaged and removed_undamaged at thresholds 1 to 5;
- kept_at_equal_removal: damaged reads kept when the scorer removes the same share of undamaged
  reads as the reference scorer does at threshold 3 (read off the ROC curve, ties as straight
  segments, so it does not depend on the scorer's scale);
- auc: the ROC area (ties count a half).

    py -3 -I tests/sim/evaluate.py default=a.tsv oracle=b.tsv --reference default
"""

import argparse
import json


def load(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        name, col = header.index("read_name"), header.index("pmds")
        return {r[name]: float(r[col]) for r in (line.rstrip("\n").split("\t") for line in f)}


def roc(scores, labels):
    """ROC points (FPR, TPR) from the highest threshold down, one point per distinct score."""
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    pos = sum(labels)
    neg = len(labels) - pos
    tp = fp = 0
    points = [(0.0, 0.0)]
    i = 0
    while i < len(pairs):
        s = pairs[i][0]
        while i < len(pairs) and pairs[i][0] == s:
            tp += pairs[i][1]
            fp += 1 - pairs[i][1]
            i += 1
        points.append((fp / neg, tp / pos))
    return points


def tpr_at(points, fpr):
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= fpr <= x1:
            return y0 if x1 == x0 else y0 + (y1 - y0) * (fpr - x0) / (x1 - x0)
    return 1.0


def auc(points):
    return sum((x1 - x0) * (y0 + y1) / 2 for (x0, y0), (x1, y1) in zip(points, points[1:]))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("scorers", nargs="+", help="NAME=TSV")
    p.add_argument("--reference", default="default", help="scorer whose threshold-3 removal is matched")
    p.add_argument("--threshold", type=float, default=3.0)
    a = p.parse_args()
    tables = {s.split("=", 1)[0]: load(s.split("=", 1)[1]) for s in a.scorers}
    names = sorted(next(iter(tables.values())))
    for t in tables.values():
        if sorted(t) != names:
            raise SystemExit("the scorers cover different reads")
    labels = [1 if n.endswith("_D") else 0 for n in names]
    n_dmg, n_und = sum(labels), len(labels) - sum(labels)
    ref = tables[a.reference]
    target_fpr = sum(1 for n, y in zip(names, labels) if not y and ref[n] >= a.threshold) / n_und
    out = {"reads": len(names), "damaged": n_dmg, "undamaged_removed_by_reference": 1 - target_fpr}
    for s, t in tables.items():
        scores = [t[n] for n in names]
        pts = roc(scores, labels)
        row = {"auc": auc(pts), "kept_at_equal_removal": tpr_at(pts, target_fpr)}
        for thr in (1, 2, 3, 4, 5):
            kd = sum(1 for x, y in zip(scores, labels) if y and x >= thr) / n_dmg
            ru = sum(1 for x, y in zip(scores, labels) if not y and x < thr) / n_und
            row[f"t{thr}"] = [kd, ru]
        out[s] = row
    print(json.dumps(out))


if __name__ == "__main__":
    main()
