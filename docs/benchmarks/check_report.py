"""Check every plotted number in a pmdshodh report against the JSON summaries it was made from.

Usage: python3 -I check_report.py report.html
Rates are recomputed independently from the raw counts in each summary named in the report.
Prints one PASS/FAIL line per check and exits 1 on any failure.
"""
import html
import json
import re
import sys

TOL = 0.0005  # values are printed with 3 decimals of a percent

checks = {"PASS": 0, "FAIL": 0}


def check(ok, what):
    checks["PASS" if ok else "FAIL"] += 1
    print(("PASS " if ok else "FAIL ") + what)


def pct(num, den):
    return None if den == 0 else 100.0 * num / den


def parse_values(s):
    return [None if v == "n/a" else float(v) for v in s.split(" ")] if s else []


def same(got, want):
    if len(got) != len(want):
        return False
    for g, w in zip(got, want):
        if (g is None) != (w is None):
            return False
        if g is not None and abs(g - w) > TOL:
            return False
    return True


def sites(t, i):
    return sum(t[b][i] for b in "ACGT")


def expected(j):
    d = j["damage"]
    rng = d["range"]
    ss = j["mode"]["library"] == "ss"
    five_c, three = d["five_prime"]["C"], d["three_prime"]
    shown3 = ("C", "T") if ss else ("G", "A")

    def damage(end, ref, base):
        t = d[end][ref]
        return [pct(t[base][i], sites(t, i)) for i in range(rng)]

    def other(end, shown):
        out = []
        for i in range(rng):
            n = 0
            for ref, same_base in (("C", "C"), ("G", "G")):
                for b in "ACGT":
                    if b != same_base and (ref, b) != shown:
                        n += d[end][ref][b][i]
            out.append(pct(n, sites(d[end]["C"], i) + sites(d[end]["G"], i)))
        return out

    scored = j["reads"]["scored"]
    e = {
        "5′ C→T": [pct(five_c["T"][i], sites(five_c, i)) for i in range(rng)],
        "3′ " + ("C→T" if ss else "G→A"): damage("three_prime", *shown3),
        "5′ other substitutions": other("five_prime", ("C", "T")),
        "3′ other substitutions": other("three_prime", shown3),
        "PMD score": [pct(n, scored) for n in j["score_histogram"]["counts"]],
    }
    lengths = lambda m: " ".join(f"{k}:{v}" for k, v in sorted((int(k), v) for k, v in m.items()))
    e["scored reads"] = lengths(j["length_histogram"]["scored"])
    if j["length_histogram"]["kept"] is not None:
        e["kept reads"] = lengths(j["length_histogram"]["kept"])
    return e


def main(path):
    page = open(path, encoding="utf-8").read()
    series_re = re.compile(r'data-series="([^"]*)" data-values="([^"]*)"')
    sections = re.findall(r'<section id="sample-\d+" data-sample="([^"]*)">(.*?)</section>', page, re.S)
    check(len(sections) > 0, f"{path}: {len(sections)} sample sections")
    by_label = {}
    for label, body in sections:
        label = html.unescape(label)
        src = html.unescape(re.search(r"<th>Summary file</th><td>(.*?)</td>", body).group(1))
        j = json.load(open(src, encoding="utf-8"))
        by_label[label] = j
        want = expected(j)
        got = {html.unescape(n): v for n, v in series_re.findall(body)}
        check(set(got) == set(want), f"{label}: series {sorted(got)}")
        for name, w in want.items():
            g = got.get(name)
            if isinstance(w, str):
                check(g == w, f"{label}: {name} ({len(w.split())} lengths)")
            else:
                check(g is not None and same(parse_values(g), w), f"{label}: {name} ({len(w)} values)")
        fp = j["damage"]["first_position"]
        if fp["c_to_t"] is not None:
            text = f"First-base C→T: {100 * fp['c_to_t']:.2f}% ± {100 * fp['se']:.2f}% (standard error), {fp['n']:,} sites"
            check(text in body, f"{label}: {text}")
            check(abs(want["5′ C→T"][0] - 100 * fp["c_to_t"]) <= 1e-9, f"{label}: first plotted 5′ point equals first_position.c_to_t")
        r = j["reads"]
        for key, row in (("input", "Input"), ("scored", "Scored"), ("no_informative_sites", "No informative sites")):
            check(f"<th>{row}</th><td>{r[key]:,}</td>" in body, f"{label}: {row} {r[key]:,}")
        if r["kept"] is not None:
            kept = f"{r['kept']:,} of {r['scored']:,} scored reads ({100 * r['kept'] / r['scored']:.1f}%)"
            check(kept in body, f"{label}: kept {kept}")
        for t in ("threshold", "upper_threshold"):
            v = j["filters"][t]
            if v is not None:
                check(f'data-threshold="{v:g}"' in body, f"{label}: {t} {v:g} marked")

    ov = re.search(r'<section id="overlays">(.*?)</section>', page, re.S)
    if len(sections) > 1:
        found = [(html.unescape(n), v) for n, v in series_re.findall(ov.group(1))]
        n = len(sections)
        check(len(found) == 3 * n, f"overlays: {len(found)} series for {n} samples")
        labels = [html.unescape(l) for l, _ in sections]
        for k, label in enumerate(labels):
            e = expected(by_label[label])
            three = [x for x in e if x.startswith("3′ ") and "other" not in x][0]
            for (name, vals), (want_name, want) in zip(
                (found[k], found[n + k], found[2 * n + k]),
                ((label, e["5′ C→T"]), (f"{label} ({three[3:]})", e[three]), (label, e["PMD score"])),
            ):
                check(name == want_name and same(parse_values(vals), want), f"overlay: {want_name} ({len(want)} values)")
    else:
        check(ov is None, "no overlay section for one sample")
    print(f"{checks['PASS']} PASS, {checks['FAIL']} FAIL")
    sys.exit(1 if checks["FAIL"] else 0)


main(sys.argv[1])
