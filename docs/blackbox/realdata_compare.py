"""Per-read parity of pmdshodh against PMDtools 0.60 on real BAMs (converted to SAM)."""
import subprocess
import sys

D = sys.argv[1]
PMDTOOLS = sys.argv[2]   # path inside WSL
BIN = sys.argv[3]        # pmdshodh.exe
WSL_D = "/mnt/c/" + D.split("C:/", 1)[1]

MODES = {"ds": [], "ds_full": ["--CpG"], "ds_half": ["--UDGhalf"], "ss": ["--ss"]}


def pmdtools(sam, opts):
    text = open(sam, "rb").read()
    run = lambda extra: subprocess.run(["wsl", "-e", "python3", "-W", "ignore", PMDTOOLS] + extra + opts,
                                       input=text, capture_output=True, check=True).stdout.decode()
    scores = [l.split("\t")[-1].strip() for l in run(["--printDS"]).splitlines() if l.strip()]
    names = [l.split("\t")[0] for l in run(["--threshold=-500"]).splitlines()
             if l and not l.startswith("@")]
    assert len(scores) == len(names), (len(scores), len(names))
    return list(zip(names, map(float, scores)))


def ours(sam, opts, compat):
    args = [BIN, "score", "--quiet", sam] + opts + (["--pmdtools-compat"] if compat else [])
    out = subprocess.run(args, capture_output=True, check=True).stdout.decode()
    return [(l.split("\t")[0], float(l.split("\t")[2])) for l in out.splitlines()[1:]]


for sample in ("mammoth", "human"):
    sam = f"{D}/{sample}.sam"
    for mode, opts in MODES.items():
        theirs = pmdtools(sam, opts)
        for compat in (True, False):
            mine = ours(sam, opts, compat)
            assert [n for n, _ in mine] == [n for n, _ in theirs], f"{sample} {mode}: read lists differ"
            diffs = [abs(a - b) for (_, a), (_, b) in zip(mine, theirs)]
            n_diff = sum(d > 1e-9 for d in diffs)
            label = "compat " if compat else "default"
            print(f"{sample:8} {mode:8} {label}: {len(diffs):6} reads, largest |diff| {max(diffs):.2e}, "
                  f"reads differing > 1e-9: {n_diff}")
