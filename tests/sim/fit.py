"""Fit the damage curve to a pmdshodh summary's damage tables: the "fitted" scorer's P.

docs/damage-curve-study.md sections 2.1 and 5. Binomial maximum likelihood, positions 1 to the
summary's range, of f(z) = b + A (1 - q)^(z-1) + C on 5' C>T and 3' G>A (ds) or C>T at both
ends (ss). b is fixed beforehand from the transition damage does not cause on the same strand
(ds: 5' G>A and 3' C>T; ss: G>A at both ends). Only q is robust (section 2.1, point 4); the
fitted scorer uses P = q. Prints the fit as JSON.

    py -3 -I tests/sim/fit.py sample.damage.json
"""

import json
import math
import sys

BASES = "ACGT"


def counts(summary, end, ref, alt):
    t = summary["damage"][end][ref]
    r = summary["damage"]["range"]
    return [t[alt][i] for i in range(r)], [sum(t[b][i] for b in BASES) for i in range(r)]


def nelder_mead(f, x0, step=0.5, iters=4000, tol=1e-12):
    n = len(x0)
    pts = [list(x0)] + [[x0[j] + (step if j == i else 0) for j in range(n)] for i in range(n)]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        xr = [2 * cen[j] - pts[-1][j] for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [3 * cen[j] - 2 * pts[-1][j] for j in range(n)]
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [(cen[j] + pts[-1][j]) / 2 for j in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [(pts[0][j] + pts[i][j]) / 2 for j in range(n)]
                    vals[i] = f(pts[i])
    i = min(range(n + 1), key=lambda i: vals[i])
    return pts[i], vals[i]


def sig(x):
    return 1 / (1 + math.exp(-x)) if x > -700 else 0.0


def fit(summary):
    ss = summary["mode"]["library"] == "ss"
    if ss:
        curves = [counts(summary, "five_prime", "C", "T"), counts(summary, "three_prime", "C", "T")]
        bg = [counts(summary, "five_prime", "G", "A"), counts(summary, "three_prime", "G", "A")]
    else:
        curves = [counts(summary, "five_prime", "C", "T"), counts(summary, "three_prime", "G", "A")]
        bg = [counts(summary, "five_prime", "G", "A"), counts(summary, "three_prime", "C", "T")]
    b = sum(sum(k) for k, _ in bg) / sum(sum(n) for _, n in bg)

    def nll(theta):
        a, q, c = sig(theta[0]), sig(theta[1]), sig(theta[2]) * 0.2
        tot = 0.0
        for k, n in curves:
            for z in range(len(k)):
                if n[z]:
                    p = min(max(b + a * (1 - q) ** z + c, 1e-12), 1 - 1e-12)
                    tot -= k[z] * math.log(p) + (n[z] - k[z]) * math.log(1 - p)
        return tot

    starts = ([-1, -1, -3], [-3, 0, -5], [-2, -2, -4], [0, 1, -4])
    theta, val = min((nelder_mead(nll, x0) for x0 in starts), key=lambda r: r[1])
    a, q, c = sig(theta[0]), sig(theta[1]), sig(theta[2]) * 0.2
    sat = 0.0
    for k, n in curves:
        for z in range(len(k)):
            if n[z]:
                p = min(max(k[z] / n[z], 1e-12), 1 - 1e-12)
                sat -= k[z] * math.log(p) + (n[z] - k[z]) * math.log(1 - p)
    df = sum(1 for _, n in curves for x in n if x) - 3
    return {"A": a, "q": q, "floor": c, "background": b, "deviance": 2 * (val - sat), "df": df}


if __name__ == "__main__":
    print(json.dumps(fit(json.load(open(sys.argv[1])))))
