"""Generate the whole-read parity fixture by running PMDtools 0.60.

Writes:
- tests/data/reads_random.sam: synthetic reads with soft clips, insertions, deletions, =/X
  CIGARs, both strands, CpG-rich references, varied base qualities, damage and other mismatches.
- tests/data/reads_random_pmdtools.tsv: the score PMDtools 0.60 prints for every read in each mode.

Only reads PMDtools scores are generated (no hard clips, N/P CIGARs, missing MD, 1-bp reads or
sites 1000+ bases from an end), so every read has a PMDtools score to compare with.
Clean room: PMDtools is only run, never read.

Usage (PMDtools 0.60 path inside WSL; see gen_site_terms.py):
    py -3 -I tests/oracle/gen_read_fixture.py <pmdtools.0.60.py in WSL> tests/data
"""
import random
import subprocess
import sys

SEED = 20261111
N_READS = 1500

MODES = {
    # column name: PMDtools options
    "ds": [],
    "ds_half": ["--UDGhalf"],
    "ds_full": ["--CpG"],
    "ss": ["--ss"],
    "ss_half": ["--ss", "--UDGhalf"],
    "ss_full": ["--ss", "--CpG"],
    "ds_q20": ["--requirebaseq", "20"],
    "ds_params": ["--PMDpparam", "0.5", "--PMDconstant", "0.02",
                  "--polymorphism_ancient", "0.002", "--polymorphism_contamination", "0.003"],
}

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def random_reference(rng, n):
    # Enrich CpG so the CpG modes see many sites.
    out = []
    while len(out) < n:
        out.extend("CG" if rng.random() < 0.15 else rng.choice("ACGT"))
    return "".join(out[:n])


def damage_molecule(rng, ref_mol, ss):
    """Apply C→T (and G→A for ds) damage in molecule orientation, plus noise."""
    seq = list(ref_mol)
    n = len(seq)
    for i, r in enumerate(ref_mol):
        d5 = 0.4 * 0.7 ** i + 0.01
        d3 = 0.4 * 0.7 ** (n - 1 - i) + 0.01
        if r == "C" and rng.random() < (d5 + (d3 if ss else 0)):
            seq[i] = "T"
        elif r == "G" and not ss and rng.random() < d3:
            seq[i] = "A"
        elif rng.random() < 0.01:
            seq[i] = rng.choice([b for b in "ACGTN" if b != r])
    return "".join(seq)


def make_read(rng, i):
    ss_like = rng.random() < 0.5
    reverse = rng.random() < 0.5
    aligned_len = rng.randint(20, 150)
    # Alignment layout: [S] M (I|D M)* [S]
    ops = []
    left_clip = rng.choice([0, 0, 0, rng.randint(1, 6)])
    right_clip = rng.choice([0, 0, 0, rng.randint(1, 6)])
    remaining = aligned_len
    while remaining > 0:
        block = min(remaining, rng.randint(8, 60))
        ops.append(("M", block))
        remaining -= block
        if remaining > 8 and rng.random() < 0.3:
            ops.append((rng.choice("ID"), rng.randint(1, 3)))
    ref_len = sum(n for k, n in ops if k in "MD")
    ref = random_reference(rng, ref_len)

    # Build the read in SAM orientation from the reference, damaging in molecule orientation.
    aligned_ref = []  # reference bases at M positions, in order
    r = 0
    for k, n in ops:
        if k == "M":
            aligned_ref.append(ref[r:r + n])
        if k in "MD":
            r += n
    m_ref = "".join(aligned_ref)
    mol = revcomp(m_ref) if reverse else m_ref
    damaged = damage_molecule(rng, mol, ss_like)
    m_seq = revcomp(damaged) if reverse else damaged

    seq, cigar, md, run, r, m = [], [], [], 0, 0, 0
    if left_clip:
        seq.append("".join(rng.choice("ACGT") for _ in range(left_clip)))
        cigar.append(f"{left_clip}S")
    use_eqx = rng.random() < 0.15
    for k, n in ops:
        if k == "M":
            for j in range(n):
                rb, qb = ref[r], m_seq[m]
                if rb == qb:
                    run += 1
                else:
                    md.append(f"{run}{rb}")
                    run = 0
                seq.append(qb)
                r += 1
                m += 1
            if use_eqx:
                block = m_seq[m - n:m]
                rblock = ref[r - n:r]
                kinds = ["=" if a == b else "X" for a, b in zip(block, rblock)]
                j = 0
                while j < n:
                    k2 = j
                    while k2 < n and kinds[k2] == kinds[j]:
                        k2 += 1
                    cigar.append(f"{k2 - j}{kinds[j]}")
                    j = k2
            else:
                cigar.append(f"{n}M")
        elif k == "I":
            seq.append("".join(rng.choice("ACGT") for _ in range(n)))
            cigar.append(f"{n}I")
        else:
            md.append(f"{run}^{ref[r:r + n]}")
            run = 0
            r += n
            cigar.append(f"{n}D")
    md.append(str(run))
    if right_clip:
        seq.append("".join(rng.choice("ACGT") for _ in range(right_clip)))
        cigar.append(f"{right_clip}S")
    seq = "".join(seq)
    qual = "".join(chr(33 + rng.choice([2, 8, 15, 20, 25, 30, 35, 37, 40, 41])) for _ in seq)
    flag = (16 if reverse else 0) | rng.choice([0, 0, 1 | 64, 1 | 128, 1024])
    return "\t".join([f"read{i:05d}", str(flag), "chrT", str(1000 + i), str(rng.choice([0, 20, 37, 60])),
                      "".join(cigar), "*", "0", "0", seq, qual, f"MD:Z:{''.join(md)}"])


def main():
    pmdtools, outdir = sys.argv[1], sys.argv[2]
    rng = random.Random(SEED)
    header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:1000000\n"
    lines = [make_read(rng, i) for i in range(N_READS)]
    sam = header + "\n".join(lines) + "\n"
    with open(f"{outdir}/reads_random.sam", "w", newline="\n") as f:
        f.write(sam)

    columns = {}
    for name, opts in MODES.items():
        out = subprocess.run(["wsl", "-e", "python3", "-W", "ignore", pmdtools, "--printDS"] + opts,
                             input=sam.encode(), capture_output=True, check=True).stdout.decode()
        scores = [line.split("\t")[-1].strip() for line in out.splitlines() if line.strip()]
        if len(scores) != N_READS:
            sys.exit(f"{name}: expected {N_READS} scores, got {len(scores)}")
        columns[name] = scores

    with open(f"{outdir}/reads_random_pmdtools.tsv", "w", newline="\n") as f:
        f.write("# PMD scores printed by PMDtools 0.60 (--printDS) for tests/data/reads_random.sam.\n")
        f.write("# Column = mode; options: " + "; ".join(f"{k}: {' '.join(v) or '(none)'}"
                                                         for k, v in MODES.items()) + "\n")
        f.write("# Generated by tests/oracle/gen_read_fixture.py; do not edit by hand.\n")
        f.write("read_name\t" + "\t".join(MODES) + "\n")
        for i in range(N_READS):
            f.write(f"read{i:05d}\t" + "\t".join(columns[m][i] for m in MODES) + "\n")


if __name__ == "__main__":
    main()
