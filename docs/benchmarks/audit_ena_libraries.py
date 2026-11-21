"""Audit section 8: the three ENA libraries. PASS/FAIL lines. Run in WSL with samtools on PATH."""
import json, os, subprocess, sys

BIN, W = sys.argv[1], "/root/bench/audit/enadata"
PMD = "/root/bench/pmdtools/PMDtools-0.60/pmdtools.0.60.py"
REF = "/root/bench/audit/enadata/MT.fa"
CHECKER = "/mnt/e/Projects/pmdshodh/docs/benchmarks/check_report.py"
n = {"PASS": 0, "FAIL": 0}


def check(ok, msg):
    n["PASS" if ok else "FAIL"] += 1
    print(("PASS " if ok else "FAIL ") + msg, flush=True)


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, executable="/bin/bash").stdout


if not os.path.exists(REF):
    sh(f"samtools faidx /root/bench/refs/human_g1k_v37.fasta MT > {REF} && samtools faidx {REF}")

LIBS = [
    ("I19359_ss_USER", "ERR7195621.I19359.MT.bam", [("ss", "--library ss", "--ss"), ("ss half", "--library ss --UDGhalf", "--ss --UDGhalf")]),
    ("I8525_ds_half", "ERR4589292.I8525.MT.bam", [("half", "--UDGhalf", "--UDGhalf"), ("none", "", "")]),
    ("I0071_ds_full", "ERR1995180.I0071.MT.bam", [("full", "--CpG", "--CpG"), ("none", "", "")]),
]
summaries = []
for name, f, modes in LIBS:
    path = f"{W}/{f}"
    p = f"{W}/{name}"
    for label, ours, theirs in modes:
        # Compat scores against PMDtools, same reads in the same order.
        sh(f"{BIN} score --quiet --pmdtools-compat {ours} -@ 8 {path} -o {p}.compat.tsv")
        sh(f"samtools view -h {path} | python3 -W ignore {PMD} --printDS {theirs} > {p}.pmd.ds")
        # PMDtools' read names in --printDS order, from a run that keeps every read it scores.
        names = sh(f"samtools view -h {path} | python3 -W ignore {PMD} --threshold=-500 --header {theirs} 2>/dev/null | samtools view | cut -f1,2").decode().splitlines()
        b = [float(l.split("\t")[-1]) for l in open(f"{p}.pmd.ds")]
        pm = {}
        for k, v in zip(names, b):
            pm.setdefault(k, []).append(v)
        d, extra, extra_sites = [], 0, 0
        for l in open(f"{p}.compat.tsv").readlines()[1:]:
            c = l.rstrip("\n").split("\t")
            k = c[0] + "\t" + c[1]
            if pm.get(k):
                d.append(abs(float(c[2]) - pm[k].pop(0)))
            else:
                extra += 1
                extra_sites += int(c[3])
        check(len(names) == len(b) and not any(x > 1e-6 for x in d),
              f"{name} [{label}]: compat vs PMDtools by name, {len(d):,} reads compared, {sum(x > 1e-6 for x in d)} above 1e-6, largest {max(d):.2g}; scored only by pmdshodh: {extra} ({extra_sites} sites)")
        # Kept reads at threshold 3.
        sh(f"{BIN} filter --quiet --pmdtools-compat {ours} -t 3 {path} -o {p}.k.bam")
        sh(f"samtools view -h {path} | python3 -W ignore {PMD} --threshold 3 --header {theirs} 2>/dev/null | samtools view -b -o {p}.pk.bam -")
        ka, kb = sh(f"samtools view {p}.k.bam | cut -f1,2"), sh(f"samtools view {p}.pk.bam | cut -f1,2")
        check(ka == kb, f"{name} [{label}]: compat kept reads (t = 3) = PMDtools, {len(ka.splitlines()):,}")
        # Damage table against PMDtools --deamination (ds tables only: C at 5', G at 3').
        table = sh(f"samtools view -h {path} | python3 -W ignore {PMD} --deamination --range 30 {theirs}").decode()
        rows = [l.split() for l in table.splitlines() if l[:1].isdigit()]
        dmg = json.loads(sh(f"{BIN} damage --quiet --pmdtools-compat --range 30 {ours} {path}"))["damage"]
        worst = 0.0
        for row in rows:
            z = int(row[0])
            def frac(end, refb, readb):
                c = dmg[end][refb]; tot = sum(c[x][z] for x in "ACGT")
                return c[readb][z] / tot if tot else 0.0
            mine = [frac("five_prime", "C", x) for x in "TAGC"] + [frac("three_prime", "G", x) for x in "ATCG"]
            worst = max([worst] + [abs(m - float(t)) for m, t in zip(mine, row[1:])])
        check(len(rows) == 30 and worst <= 5e-6, f"{name} [{label}]: damage table vs PMDtools --deamination, largest {worst:.1g}")
    ours0 = modes[0][1]
    # Counts, threads, -r against MD, report.
    sh(f"{BIN} damage --quiet {ours0} {path} -o {p}.d.json")
    j = json.load(open(f"{p}.d.json"))["reads"]
    total = int(sh(f"samtools view -c {path}"))
    check(j["input"] == total and j["input"] == j["scored"] + sum(j["skipped"].values()),
          f"{name}: input {j['input']:,} = samtools {total:,}, scored {j['scored']:,}, skipped {sum(j['skipped'].values()):,}")
    for t in (1, 8):
        sh(f"{BIN} score --quiet {ours0} -@ {t} {path} -o {p}.t{t}.tsv")
    check(open(f"{p}.t1.tsv", "rb").read() == open(f"{p}.t8.tsv", "rb").read(), f"{name}: score -@ 1 = -@ 8")
    one_seq = int(sh(f"samtools view -H {path} | grep -c '^@SQ'")) == 1
    ref = REF if one_seq else "/root/bench/refs/human_g1k_v37.fasta"
    sh(f"{BIN} score --quiet {ours0} -r {ref} {path} -o {p}.r.tsv --summary {p}.r.json")
    md = json.load(open(f"{p}.r.json"))["reference"]["md_tags"]
    nd = int(sh(f"diff {p}.t1.tsv {p}.r.tsv | grep -c '^>' || true"))
    check(nd <= md["disagree"], f"{name}: -r vs MD, {nd:,} reads differ, {md['disagree']:,} MD tags disagree ({md['agree']:,} agree)")
    sh(f"{BIN} filter --quiet {ours0} -t 3 {path} -o {p}.f.bam --summary {p}.f.json")
    summaries.append(f"{p}.f.json")
    sh(f"{BIN} report --quiet {p}.f.json -o {p}.html")
    r = subprocess.run(f"python3 -I {CHECKER} {p}.html", shell=True, capture_output=True, text=True)
    check(r.returncode == 0, f"{name}: report vs JSON, {r.stdout.strip().splitlines()[-1]}")
    fp = json.load(open(f"{p}.f.json"))["damage"]["first_position"]
    print(f"INFO {name}: first-base C>T {100 * fp['c_to_t']:.2f}% ({fp['n']:,} sites); kept at t = 3 (default): {json.load(open(p + '.f.json'))['reads']['kept']:,}")
sh(f"{BIN} report --quiet {' '.join(summaries)} -o {W}/all.html")
print(f"{n['PASS']} PASS, {n['FAIL']} FAIL")
