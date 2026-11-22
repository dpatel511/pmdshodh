"""0.3 audit, section 7 check "default mode against compat": every read whose default score differs
from its --pmdtools-compat score, and every read only one mode scores, must be explained by a
departure listed in docs/design.md section 5. Run under WSL with samtools on PATH:
    python3 -I docs/benchmarks/audit_default_vs_compat.py <pmdshodh binary> <workdir> [names]
Input paths are those of the bench machine (/root/bench; docs/benchmarks.md and docs/audit-0.3.md
sections 7 and 8 list the sources). Prints one PASS/FAIL/INFO line per check. Compat itself is
compared with PMDtools by audit_matrix.py and audit_ena_libraries.py; this script does not run it.

For each read, the triggers of the score departures are found from its SAM fields:
  reconstruction   CIGAR combines I or S with D or = (blackbox.md 1.6c)
  lowercase MD     MD has a lowercase letter
  reverse shift    reverse read whose MD has a character other than ACGTNRY/acgt
  CpG wrap         --UDGhalf or --CpG, and the reference along SEQ starts with G and ends with C
  ss sites         --library ss and the read has a C site (one term per site, SI Eq. S1)
and the triggers of the read-set departures:
  placed unmapped  flag 4 with a CIGAR: compat scores it, default skips it (B4)
  hard clip        an H in the CIGAR: compat drops it (B6)
  1-base SEQ       compat drops it (B8)
  no ref bases     no uppercase A, C, G or T in the reference along SEQ: compat drops it (B5)
A differing read with no trigger is UNEXPLAINED, and so is a read scored only by compat that is not
placed unmapped, or only by default without one of the other three. Reads whose only triggers are
lowercase MD and the reverse shift are also checked by removing them (uppercasing the MD, and
writing N for the dropped characters) and scoring them again in compat: the score must then equal
default's."""
import collections, os, re, subprocess, sys

BIN, W = sys.argv[1], sys.argv[2]
B = "/root/bench"
os.makedirs(W, exist_ok=True)
counts = {"PASS": 0, "FAIL": 0}
TOL = 1e-9


def out(kind, msg):
    if kind in counts:
        counts[kind] += 1
    print(f"{kind} {msg}", flush=True)


def check(ok, msg):
    out("PASS" if ok else "FAIL", msg)


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, executable="/bin/bash").stdout


LIBS = [  # name, path, library
    ("BER11_L_NOuser", f"{B}/data/BER11_L_NOuser.bam", "ds"),
    ("Motala3", f"{B}/data/Motala3.bam", "ds"),
    ("JK2067", f"{B}/data/JK2067.bam", "ds"),
    ("JK2782", f"{B}/eager_report/JK2782.bam", "ds"),
    ("UKY001_1240K", f"{B}/data_ss/UKY001.B0102.1240K.bam", "ss"),
    ("UKY001_MT", f"{B}/data_ss/UKY001.B0102.MT.bam", "ss"),
    ("I19359_MT", f"{B}/audit/enadata/ERR7195621.I19359.MT.bam", "ss"),
    ("I8525_MT", f"{B}/audit/enadata/ERR4589292.I8525.MT.bam", "ds"),
    ("I0071_MT", f"{B}/audit/enadata/ERR1995180.I0071.MT.bam", "ds"),
    ("multiref", f"{B}/audit/newdata/multiref/multiref.bam", "ds"),
] + [
    (f"sweep_{f[:-4]}", f"{B}/sweep/aln/{f}", "ds")
    for f in sorted(os.listdir(f"{B}/sweep/aln"))
    if f.endswith(".bam") and f != "calmd_eq.bam"  # its `=` bases need -r in default mode
] + [  # The ss libraries in ds mode too, where the ss trigger does not cover every read.
    (f"{name}_as_ds", path, "ds")
    for name, path in (
        ("UKY001_1240K", f"{B}/data_ss/UKY001.B0102.1240K.bam"),
        ("UKY001_MT", f"{B}/data_ss/UKY001.B0102.MT.bam"),
        ("I19359_MT", f"{B}/audit/enadata/ERR7195621.I19359.MT.bam"),
    )
]
if len(sys.argv) > 3:  # optional: a comma-separated list of names to run
    LIBS = [x for x in LIBS if x[0] in sys.argv[3].split(",")]
TREATMENTS = (("none", []), ("half", ["--UDGhalf"]), ("full", ["--CpG"]))
OPS = re.compile(r"(\d+)([MIDNSHP=X])")
MD_TOK = re.compile(r"(\d+)|(\^[A-Za-z]+)|([A-Za-z])")


def reference_along_seq(cigar, seq, md):
    """Reference base at each SEQ position (None for clips and insertions), from CIGAR and MD.
    Returns None if they do not fit together (such reads are input errors in pmdshodh)."""
    md_bases = []  # reference base for each aligned (M/=/X) position: None = match
    for num, dele, mm in MD_TOK.findall(md):
        if num:
            md_bases += [None] * int(num)
        elif mm:
            md_bases.append(mm)
    refs, q, a = [], 0, 0
    for n, op in OPS.findall(cigar):
        n = int(n)
        if op in "M=X":
            for _ in range(n):
                if a >= len(md_bases) or q >= len(seq):
                    return None
                refs.append(md_bases[a] if md_bases[a] is not None else seq[q])
                a += 1
                q += 1
        elif op in "IS":
            refs += [None] * n
            q += n
    return refs


def triggers(f, md, lib, treatment):
    cigar, seq, flag = f[5], f[9], int(f[1])
    t = []
    if re.search(r"[IS]", cigar) and re.search(r"[D=]", cigar):
        t.append("reconstruction")
    if re.search(r"[a-z]", md):
        t.append("lowercase MD")
    if flag & 16 and re.search(r"[^0-9^ACGTNRYacgt]", md):
        t.append("reverse shift")
    if treatment != "none":
        refs = reference_along_seq(cigar, seq, md)
        if refs and refs[0] == "G" and refs[-1] == "C":
            t.append("CpG wrap")
    if lib == "ss":
        t.append("ss sites")
    return t


def set_triggers(f, md, only_compat):
    """Read-set triggers for the direction seen: compat scores placed unmapped reads; it drops the
    other three."""
    cigar, seq, flag = f[5], f[9], int(f[1])
    if only_compat:
        return ["placed unmapped"] if flag & 4 and cigar != "*" else []
    t = []
    if "H" in cigar:
        t.append("hard clip")
    if len(seq) == 1:
        t.append("1-base SEQ")
    refs = reference_along_seq(cigar, seq, md) if cigar != "*" else None
    if refs is not None and not any(r in ("A", "C", "G", "T") for r in refs):
        t.append("no ref bases")
    return t


def rows(path):
    with open(path) as fh:
        assert next(fh).startswith("read_name\t")
        for line in fh:
            name, flag, pmds, n = line.rstrip("\n").split("\t")
            yield name, flag, float(pmds), int(n)


def neutralise(f, md):
    # Uppercase the MD; then characters PMDtools drops on reverse reads become N, which both modes
    # read alike.
    new_md = re.sub(r"[^0-9^ACGTNRY]", "N", md.upper())
    g = list(f)
    g[11:] = [x if not x.startswith("MD:Z:") else "MD:Z:" + new_md for x in f[11:]]
    return "\t".join(g)


out("INFO", f"binary {sh(f'{BIN} --version').decode().strip()}; samtools {sh('samtools --version | head -1').decode().split()[1]}")
grand = collections.Counter()
for name, path, lib in LIBS:
    ss = ["--library", "ss"] if lib == "ss" else []
    header = sh(f"samtools view -H {path}").decode()
    for tname, topt in TREATMENTS:
        p = f"{W}/{name}.{tname}"
        opts = " ".join(ss + topt)
        sh(f"{BIN} score --quiet {opts} -@ 8 {path} -o {p}.default.tsv")
        sh(f"{BIN} score --quiet --pmdtools-compat {opts} -@ 8 {path} -o {p}.compat.tsv")
        d_rows, c_rows = rows(f"{p}.default.tsv"), rows(f"{p}.compat.tsv")
        d_next, c_next = next(d_rows, None), next(c_rows, None)
        causes, sets = collections.Counter(), collections.Counter()
        unexplained, retest = [], []
        n = both = differ = 0
        sam = subprocess.Popen(["samtools", "view", path], stdout=subprocess.PIPE, text=True)
        for line in sam.stdout:
            n += 1
            f = line.rstrip("\n").split("\t")
            key = (f[0], f[1])
            d = c = None
            if d_next is not None and (d_next[0], d_next[1]) == key:
                d, d_next = d_next, next(d_rows, None)
            if c_next is not None and (c_next[0], c_next[1]) == key:
                c, c_next = c_next, next(c_rows, None)
            if d is None and c is None:
                continue
            md = next((x[5:] for x in f[11:] if x.startswith("MD:Z:")), "")
            if d is None or c is None:
                t = set_triggers(f, md, d is None)
                k = ("only compat: " if d is None else "only default: ") + (" + ".join(t) or "UNEXPLAINED")
                sets[k] += 1
                if not t and len(unexplained) < 5:
                    unexplained.append((k, f[0], f[1], f[5], md))
                continue
            both += 1
            if abs(d[2] - c[2]) <= TOL:
                continue
            differ += 1
            t = triggers(f, md, lib, tname)
            if lib == "ss" and c[3] == 0:
                t.remove("ss sites")
            k = " + ".join(t) or "UNEXPLAINED"
            causes[k] += 1
            if not t and len(unexplained) < 5:
                unexplained.append((k, f[0], f[1], f[5], md, d[2], c[2]))
            if t and set(t) <= {"lowercase MD", "reverse shift"}:
                retest.append((neutralise(f, md), d[2]))
        assert sam.wait() == 0
        check(d_next is None and c_next is None, f"{name} {tname}: every scored read matched to its record")
        expl = sum(v for k, v in causes.items() if k != "UNEXPLAINED")
        check(causes["UNEXPLAINED"] == 0,
              f"{name} {tname}: {n:,} records, {both:,} scored by both, {differ:,} differ, {expl:,} explained"
              + "".join(f"; {k} {v:,}" for k, v in causes.most_common()))
        set_bad = sum(v for k, v in sets.items() if k.endswith("UNEXPLAINED"))
        check(set_bad == 0, f"{name} {tname}: reads scored by one mode only, {sum(sets.values()):,}"
              + "".join(f"; {k} {v:,}" for k, v in sets.most_common()))
        for u in unexplained:
            out("INFO", f"{name} {tname}: unexplained {u}")
        if retest:
            sam_path = f"{p}.retest.sam"
            with open(sam_path, "w") as fh:
                fh.write(header)
                fh.writelines(r + "\n" for r, _ in retest)
            sh(f"{BIN} score --quiet --pmdtools-compat {opts} {sam_path} -o {p}.retest.tsv")
            again = [x[2] for x in rows(f"{p}.retest.tsv")]
            bad = sum(abs(a - b) > TOL for a, (_, b) in zip(again, retest))
            check(len(again) == len(retest) and bad == 0,
                  f"{name} {tname}: {len(retest):,} reads whose only causes are lowercase MD or the reverse shift, re-scored in compat without them, equal default ({bad} differ)")
        for k, v in causes.items():
            grand[k] += v
        for k, v in sets.items():
            grand[k] += v
        os.remove(f"{p}.default.tsv")
        os.remove(f"{p}.compat.tsv")
out("INFO", "all runs: " + "; ".join(f"{k} {v:,}" for k, v in grand.most_common()))
print(f"{counts['PASS']} PASS, {counts['FAIL']} FAIL")
