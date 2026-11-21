"""0.3 audit, real-data matrix (docs/audit-0.3.md section 7). Run under WSL with samtools on PATH:
    python3 -I docs/benchmarks/audit_matrix.py <pmdshodh binary> <previous release binary> <workdir>
Input paths are those of the bench machine (/root/bench; docs/benchmarks.md lists the sources).
Prints one PASS/FAIL/INFO line per check. PMDtools is only run (clean room)."""
import hashlib, json, os, re, subprocess, sys, time

BIN, OLD, W = sys.argv[1], sys.argv[2], sys.argv[3]
B = "/root/bench"
PMDTOOLS = f"{B}/pmdtools/PMDtools-0.60/pmdtools.0.60.py"
EQ, HG = f"{B}/refs/equCab2.fa", f"{B}/refs/human_g1k_v37.fasta"
CHECKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_report.py")
os.makedirs(W, exist_ok=True)
counts = {"PASS": 0, "FAIL": 0}


def out(kind, msg):
    if kind in counts:
        counts[kind] += 1
    print(f"{kind} {msg}", flush=True)


def check(ok, msg):
    out("PASS" if ok else "FAIL", msg)


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, executable="/bin/bash", **kw).stdout


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def strip_command(path):
    """Output without its recorded command line (JSON `command`, SAM/BAM @PG)."""
    if path.endswith(".json"):
        j = json.load(open(path)); j.pop("command"); return json.dumps(j, sort_keys=True).encode()
    if path.endswith(".bam"):
        return sh(f"samtools view -h {path} | grep -v '^@PG'")
    return open(path, "rb").read()


LIBS = [
    # name, input, library, reference, PMDtools .ds, PMDtools kept BAM, md5
    ("BER11_L_NOuser", f"{B}/data/BER11_L_NOuser.bam", "ds", EQ, f"{B}/work/BER11_L_NOuser.pmdtools.ds", f"{B}/work/BER11_L_NOuser.pmdtools.bam", "7a73fd9c2895576d8741a6bbca80b10a"),
    ("Motala3", f"{B}/data/Motala3.bam", "ds", HG, f"{B}/work/Motala3.pmdtools.ds", f"{B}/work/Motala3.pmdtools.bam", "92987d313786940f801c9977fcf515f6"),
    ("JK2067", f"{B}/data/JK2067.bam", "ds", None, f"{B}/work/JK2067.pmdtools.ds", f"{B}/work/JK2067.pmdtools.bam", None),
    ("UKY001_1240K", f"{B}/data_ss/UKY001.B0102.1240K.bam", "ss", None, f"{B}/work_ss/UKY001.B0102.1240K.pmdtools_ss.ds", f"{B}/work_ss/UKY001.B0102.1240K.pmdtools_ss.bam", None),
    ("UKY001_MT", f"{B}/data_ss/UKY001.B0102.MT.bam", "ss", None, f"{B}/work_ss/UKY001.B0102.MT.pmdtools_ss.ds", f"{B}/work_ss/UKY001.B0102.MT.pmdtools_ss.bam", None),
    ("JK2782", f"{B}/eager_report/JK2782.bam", "ds", None, None, None, None),
]

out("INFO", f"binary {sh(f'{BIN} --version').decode().strip()}; samtools {sh('samtools --version | head -1').decode().split()[1]}; "
    f"PMDtools archive md5 {md5(f'{B}/pmdtools/0.60.tar.gz')}")

summaries = []
for name, path, lib, ref, ds, kept_bam, want_md5 in LIBS:
    t0 = time.monotonic()
    if want_md5 and len(want_md5) == 32:
        check(md5(path) == want_md5, f"{name}: input md5 {want_md5}")
    ss = ["--library", "ss"] if lib == "ss" else []
    p = f"{W}/{name}"
    # PMDtools for the small eager file is run here.
    if ds is None:
        opt = "--ss" if lib == "ss" else ""
        sh(f"samtools view -h {path} | python3 -W ignore {PMDTOOLS} --printDS {opt} > {p}.pmdtools.ds")
        sh(f"samtools view -h {path} | python3 -W ignore {PMDTOOLS} --threshold 3 --header {opt} | samtools view -b -o {p}.pmdtools.bam -")
        ds, kept_bam = f"{p}.pmdtools.ds", f"{p}.pmdtools.bam"

    # 1. Compat scores against PMDtools (same reads, same order).
    sh(f"{BIN} score --quiet --pmdtools-compat {' '.join(ss)} -@ 8 {path} -o {p}.compat.tsv")
    ours = [float(l.split("\t")[2]) for l in open(f"{p}.compat.tsv").readlines()[1:]]
    theirs = [float(l.split("\t")[-1]) for l in open(ds)]
    diff = max((abs(a - b) for a, b in zip(ours, theirs)), default=0.0)
    over = sum(abs(a - b) > 1e-6 for a, b in zip(ours, theirs))
    check(len(ours) == len(theirs) and over == 0, f"{name}: compat scores vs PMDtools, {len(ours):,} vs {len(theirs):,} reads, {over} above 1e-6, largest {diff:.2g}")

    # 2. Compat kept reads (threshold 3) against PMDtools.
    sh(f"{BIN} filter --quiet --pmdtools-compat {' '.join(ss)} -t 3 -@ 8 {path} -o {p}.compat3.bam")
    a = sh(f"samtools view {p}.compat3.bam | cut -f1,2").split(b"\n")
    b = sh(f"samtools view {kept_bam} | cut -f1,2").split(b"\n")
    check(a == b, f"{name}: compat kept reads vs PMDtools --threshold 3, {len(a) - 1:,} vs {len(b) - 1:,}")

    # 3. Read counts against samtools.
    sh(f"{BIN} damage --quiet {' '.join(ss)} -@ 8 {path} -o {p}.damage.json")
    j = json.load(open(f"{p}.damage.json"))["reads"]
    n = int(sh(f"samtools view -c {path}"))
    check(j["input"] == n, f"{name}: input {j['input']:,} = samtools view -c {n:,}")
    for flag, key in ((1024, "duplicate"), (256, "secondary"), (2048, "supplementary"), (512, "qc_fail"), (1, "paired"), (16, "reverse")):
        m = int(sh(f"samtools view -c -f {flag} {path}"))
        check(j["flags"][key] == m, f"{name}: {key} {j['flags'][key]:,} = samtools -f {flag} {m:,}")
    unmapped = int(sh(f"samtools view -c -f 4 {path}"))
    check(j["skipped"]["unmapped"] == unmapped, f"{name}: skipped unmapped {j['skipped']['unmapped']:,} = samtools -f 4 {unmapped:,}")
    check(j["input"] == j["scored"] + sum(j["skipped"].values()), f"{name}: input = scored + skipped")

    # 4. Threads: -@ 1 against -@ 8, default mode.
    for cmd, ext in (("score", "tsv"), ("damage", "json"), ("filter -t 3", "bam")):
        for t in (1, 8):
            sh(f"{BIN} {cmd} --quiet {' '.join(ss)} -@ {t} {path} -o {p}.t{t}.{ext}")
        check(strip_command(f"{p}.t1.{ext}") == strip_command(f"{p}.t8.{ext}"), f"{name}: {cmd.split()[0]} -@ 1 = -@ 8")

    # 5. Filter output: valid, and kept records byte-identical to their input records, in order.
    sh(f"{BIN} filter --quiet {' '.join(ss)} -t 3 -@ 8 {path} -o {p}.f3.bam --summary {p}.f3.json")
    check(subprocess.run(f"samtools quickcheck {p}.f3.bam", shell=True).returncode == 0, f"{name}: filter output passes samtools quickcheck")
    proc = subprocess.run(
        f"python3 -I -c \"import sys,subprocess\n"
        f"a=subprocess.Popen(['samtools','view','{path}'],stdout=subprocess.PIPE).stdout\n"
        f"k=0;bad=0\n"
        f"for line in subprocess.Popen(['samtools','view','{p}.f3.bam'],stdout=subprocess.PIPE).stdout:\n"
        f"    for x in a:\n"
        f"        if x==line: k+=1; break\n"
        f"    else: bad+=1; break\n"
        f"print(k,bad)\"", shell=True, capture_output=True, text=True)
    k, bad = map(int, proc.stdout.split())
    summaries.append(f"{p}.f3.json")
    check(bad == 0 and k == json.load(open(f"{p}.f3.json"))["reads"]["kept"], f"{name}: {k:,} kept records are input records, unchanged, in order")

    # 6. -r against MD (default): only reads whose MD disagrees with the reference may differ.
    if ref:
        sh(f"{BIN} score --quiet -@ 8 {path} -o {p}.md.tsv")
        sh(f"{BIN} score --quiet -@ 8 -r {ref} {path} -o {p}.r.tsv --summary {p}.r.json")
        dis = json.load(open(f"{p}.r.json"))["reference"]["md_tags"]
        nd = int(sh(f"diff {p}.md.tsv {p}.r.tsv | grep -c '^>' || true"))
        check(nd <= dis["disagree"], f"{name}: -r vs MD, {nd:,} reads differ, {dis['disagree']:,} MD tags disagree with the reference ({dis['agree']:,} agree)")

    # 7. Report against its JSON.
    sh(f"{BIN} report --quiet {p}.f3.json -o {p}.html")
    r = subprocess.run(f"python3 -I {CHECKER} {p}.html", shell=True, capture_output=True, text=True)
    check(r.returncode == 0, f"{name}: report numbers vs JSON: {r.stdout.strip().splitlines()[-1]}")
    out("INFO", f"{name}: done in {time.monotonic() - t0:.0f} s")

# Report on all libraries together.
sh(f"{BIN} report --quiet {' '.join(summaries)} -o {W}/all.html")
r = subprocess.run(f"python3 -I {CHECKER} {W}/all.html", shell=True, capture_output=True, text=True)
check(r.returncode == 0, f"all libraries: report numbers vs JSON: {r.stdout.strip().splitlines()[-1]}")

# 8. Damage tables, compat, against PMDtools --deamination on 200,000-read subsets.
for sub in ("BER11_L_NOuser", "Motala3"):
    path = f"{B}/udgreal/{sub}.200k.bam"
    for mode, ours_flag, theirs_flag in (("none", "", ""), ("half", "--UDGhalf", "--UDGhalf"), ("full", "--CpG", "--CpG")):
        table = sh(f"samtools view -h {path} | python3 -W ignore {PMDTOOLS} --deamination --range 30 {theirs_flag}").decode()
        rows = [l.split() for l in table.splitlines() if l[:1].isdigit()]
        d = json.loads(sh(f"{BIN} damage --quiet --pmdtools-compat --range 30 {ours_flag} {path}"))["damage"]
        worst = 0.0
        for row in rows:
            z = int(row[0])
            def frac(end, refb, readb):
                c = d[end][refb]; tot = sum(c[x][z] for x in "ACGT")
                return c[readb][z] / tot if tot else 0.0
            mine = [frac("five_prime", "C", x) for x in "TAGC"] + [frac("three_prime", "G", x) for x in "ATCG"]
            worst = max([worst] + [abs(m - float(t)) for m, t in zip(mine, row[1:])])
        check(len(rows) == 30 and worst <= 5e-6, f"{sub} 200k {mode}: damage table vs PMDtools --deamination, largest difference {worst:.1g}")

# 9. Aligner sweep: threads and counts on every file.
for f in sorted(os.listdir(f"{B}/sweep/aln")):
    if not f.endswith((".bam", ".cram")) or f.endswith(".bai"):
        continue
    path = f"{B}/sweep/aln/{f}"
    rr = ["-r", f"{B}/sweep/ref/mt.fa"] if f.endswith(".cram") else []
    try:
        for t in (1, 8):
            sh(f"{BIN} score --quiet {' '.join(rr)} -@ {t} {path} -o {W}/sweep.t{t}.tsv --summary {W}/sweep.t{t}.json")
        same = open(f"{W}/sweep.t1.tsv", "rb").read() == open(f"{W}/sweep.t8.tsv", "rb").read()
        j = json.load(open(f"{W}/sweep.t1.json"))["reads"]
        n = int(sh(f"samtools view -c {' '.join(['-T', rr[1]] if rr else [])} {path}"))
        check(same and j["input"] == n, f"sweep {f}: -@ 1 = -@ 8, input {j['input']:,} = samtools {n:,}")
    except subprocess.CalledProcessError as e:
        out("INFO", f"sweep {f}: exit {e.returncode}: {e.stderr.decode().strip()[:200]}")

# 10. BER11_L with the USER-treated libraries, whose MD tags do not match their CIGARs.
r = subprocess.run(f"{BIN} score --quiet -@ 8 {B}/data/BER11_L.bam -o {W}/berall.tsv", shell=True, capture_output=True, text=True)
out("INFO", f"BER11_L (all libraries): exit {r.returncode}: {r.stderr.strip()[:300]}")

# 11. Speed and memory against the previous release, same commands, alternating, median of 3.
for name, path in (("BER11_L_NOuser", f"{B}/data/BER11_L_NOuser.bam"), ("Motala3", f"{B}/data/Motala3.bam")):
    for t in (1, 8):
        res = {"old": [], "new": []}
        rss = {"old": 0, "new": 0}
        for _ in range(3):
            for k, binary in (("old", OLD), ("new", BIN)):
                t0 = time.monotonic()
                r = subprocess.run(f"/usr/bin/time -v {binary} filter --quiet -t 3 -@ {t} {path} -o {W}/speed.bam", shell=True, capture_output=True, text=True, check=True)
                res[k].append(time.monotonic() - t0)
                rss[k] = max(rss[k], int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", r.stderr).group(1)))
        med = {k: sorted(v)[1] for k, v in res.items()}
        out("INFO", f"speed {name} filter -t 3 -@ {t}: old {med['old']:.2f} s, new {med['new']:.2f} s "
            f"(old {[round(x, 2) for x in res['old']]}, new {[round(x, 2) for x in res['new']]}); "
            f"max RSS old {rss['old'] / 1024:.0f} MB, new {rss['new'] / 1024:.0f} MB")

print(f"{counts['PASS']} PASS, {counts['FAIL']} FAIL", flush=True)
