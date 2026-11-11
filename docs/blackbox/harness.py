"""Black-box harness for PMDtools 0.60 (clean room: runs the script, never reads it).

Builds synthetic SAM records from (reference, read, CIGAR) and computes MD per the SAM spec.
Also holds an independent implementation of the Skoglund 2014 model (paper Eqs. 1-5, SI S1).
"""
import math
import os
import subprocess
import sys

# Windows path to pmdtools.0.60.py (Bioconda 0.60 archive), e.g. set PMDTOOLS_PY=C:/tools/PMDtools-0.60/pmdtools.0.60.py
PMDTOOLS = os.environ.get("PMDTOOLS_PY", "C:/tools/PMDtools-0.60/pmdtools.0.60.py")
REF_NAME = "chrT"
REF_LEN = 10000

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def parse_cigar(cigar):
    ops, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""
    return ops


def make_md(ref, seq, cigar):
    """ref: reference bases covered by the alignment (M/D ops). seq: SAM SEQ (incl. clips/insertions)."""
    md, run, r, q = [], 0, 0, 0
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            for _ in range(n):
                if ref[r] == seq[q]:
                    run += 1
                else:
                    md.append(str(run)); md.append(ref[r]); run = 0
                r += 1; q += 1
        elif op == "D":
            md.append(str(run)); md.append("^" + ref[r:r + n]); run = 0
            r += n
        elif op in "IS":
            q += n
        elif op == "H":
            pass
    md.append(str(run))
    return "".join(md)


def sam_record(name, ref, seq, cigar="", flag=0, pos=100, mapq=37, qual=None, md=None, extra=""):
    if not cigar:
        cigar = f"{len(seq)}M"
    if qual is None:
        qual = "I" * len(seq)  # Q40
    if md is None:
        md = make_md(ref, seq, cigar)
    fields = [name, str(flag), REF_NAME, str(pos), str(mapq), cigar, "*", "0", "0", seq, qual, f"MD:Z:{md}"]
    if extra:
        fields.append(extra)
    return "\t".join(fields)


HEADER = f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{REF_NAME}\tLN:{REF_LEN}\n"


def run(records, args, header=True):
    text = (HEADER if header else "") + "\n".join(records) + "\n"
    # PMDtools needs /dev/null, so run it under WSL (Ubuntu 24.04, Python 3.12.3).
    wsl_path = "/mnt/c/" + PMDTOOLS.split("C:/", 1)[1]
    p = subprocess.run(["wsl", "-e", "python3", "-W", "ignore", wsl_path] + args,
                       input=text.encode(), capture_output=True)
    return p.returncode, p.stdout.decode(errors="replace"), p.stderr.decode(errors="replace")


def scores(records, extra_args=()):
    """Return {readname: score} from --printDS output (format determined empirically)."""
    rc, out, err = run(records, ["--printDS"] + list(extra_args))
    res = {}
    for line in out.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        res[f[0]] = f
    return rc, res, err


# ---------- independent model from the paper ----------
def D(z, P=0.3, C=0.01):
    return (1 - P) ** (z - 1) * P + C


def p_match_ds(z, Q, pmd, P=0.3, C=0.01, pi=0.001):
    e = 10 ** (-Q / 10) / 3
    d = D(z, P, C) if pmd else 0.0
    return (1 - pi) * (1 - e) * (1 - d) + (1 - pi) * e * d + pi * e * (1 - d)


def p_match_ss(z, y, Q, pmd, P=0.3, C=0.01, pi=0.001):
    e = 10 ** (-Q / 10) / 3
    a = D(z, P, C) if pmd else 0.0
    b = D(y, P, C) if pmd else 0.0
    return ((1 - pi) * (1 - e) * (1 - a) * (1 - b) + (1 - pi) * e * a * (1 - b)
            + (1 - pi) * e * b * (1 - a) + pi * e * (1 - a) * (1 - b))


def p_match_ds_null_pmdtools(Q, pi_c=0.001, bg=0.001):
    """Null model as observed from PMDtools 0.60 output: Eq. 1 with D replaced by a constant 0.001."""
    e = 10 ** (-Q / 10) / 3
    return (1 - pi_c) * (1 - e) * (1 - bg) + (1 - pi_c) * e * bg + pi_c * e * (1 - bg)


def llr_pmdtools_ds(match, z, Q, P=0.3, C=0.01, pi_a=0.001, pi_c=0.001):
    return site_llr(match, p_match_ds(z, Q, True, P, C, pi_a), p_match_ds_null_pmdtools(Q, pi_c))


def p_match_ss_null_pmdtools(Q, pi_c=0.001, bg=0.001):
    e = 10 ** (-Q / 10) / 3
    a = b = bg
    return ((1 - pi_c) * (1 - e) * (1 - a) * (1 - b) + (1 - pi_c) * e * a * (1 - b)
            + (1 - pi_c) * e * b * (1 - a) + pi_c * e * (1 - a) * (1 - b))


def llr_pmdtools_ss(match, z, y, Q, P=0.3, C=0.01, pi_a=0.001, pi_c=0.001):
    """Observed --ss behaviour: ds 5' term for the C site PLUS an Eq. S1 term for the same site."""
    t1 = site_llr(match, p_match_ds(z, Q, True, P, C, pi_a), p_match_ds_null_pmdtools(Q, pi_c))
    t2 = site_llr(match, p_match_ss(z, y, Q, True, P, C, pi_a), p_match_ss_null_pmdtools(Q, pi_c))
    return t1 + t2


def site_llr(match, pm, pn):
    return math.log(pm / pn) if match else math.log((1 - pm) / (1 - pn))
