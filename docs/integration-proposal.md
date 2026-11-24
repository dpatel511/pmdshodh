# Proposal: pipeline integration (roadmap milestone 7)

Status: **agreed 2026-11-23.** I1–I6 were accepted as recommended (section 5). Both modules are built and tested locally (section 7); nothing has been sent upstream.

## 1. Goal and scope

Make pmdshodh usable where ancient DNA analyses actually run. Roadmap milestone 7 (docs/roadmap.md):
- an nf-core module, `pmdshodh/filter`;
- a MultiQC module that reads the JSON run summary;
- a proposal to nf-core/eager to offer pmdshodh next to PMDtools.

The milestone is done when the modules are accepted upstream.

In scope:
- the two modules, their tests and their documentation;
- the local tools needed to build and test them.

Out of scope:
- changes to pmdshodh's scores or its summary format (`format_version` 4 stays);
- writing eager's subworkflow ourselves before eager's maintainers agree;
- fitted damage parameters (milestone 6).

## 2. Where things stand (checked 2026-11-23)

| Item | State |
|---|---|
| Bioconda | `bioconda/pmdshodh` has 0.3.0 and 0.3.1. The bot's 0.3.2 bump (bioconda/bioconda-recipes#72119) has passed every check and waits for a member's review; 0.3.3 should follow the same way. |
| Containers | `quay.io/biocontainers/pmdshodh:0.3.1--hfa8f182_0` and the Galaxy depot Singularity image exist. These are built from Bioconda, so a module can use them now. |
| nf-core/modules | `pmdtools/filter` exists: `samtools calmd \| pmdtools --threshold N --header \| samtools view -Sb`, with inputs `(meta, bam, bai)`, `val threshold` and `path reference`, and a PMDtools 0.60 container. There is no pmdshodh module. |
| MultiQC | Has `damageprofiler` and `mapdamage` modules; none for PMDtools or pmdshodh. |
| nf-core/eager | Latest release 2.5.3. The `dev` branch (the DSL2 rewrite) uses the nf-core `pmdtools/filter` module. |
| Local tools (WSL) | Docker, Java and conda are installed; Nextflow, nf-core tools, nf-test and MultiQC are not. |

## 3. Design

### 3.1 nf-core module `pmdshodh/filter`

A drop-in for `pmdtools/filter`, without the `calmd` and `samtools view` steps:

```nextflow
input:
tuple val(meta), path(bam)                 // BAM, SAM or CRAM; no index needed
val(threshold)                             // as in pmdtools/filter
tuple val(meta2), path(fasta), path(fai)   // optional ([]); required for CRAM

output:
tuple val(meta), path("*.bam"),  emit: bam
tuple val(meta), path("*.json"), emit: json      // run summary, for MultiQC
// versions, as a topic channel (the current nf-core convention)
```

```bash
pmdshodh filter --threshold ${threshold} ${reference ? "-r $fasta" : ""} \
    --summary ${prefix}.pmdshodh.json -@ ${task.cpus} ${args} ${bam} -o ${prefix}.bam
```

- `ext.args` carries everything else: `--library ss`, `--treatment half`, `--tag` and so on. PMDtools spellings (`--UDGhalf`) also work, which eases moving from `pmdtools/filter`.
- **Tests:** nf-test with the nf-core test data the PMDtools module already uses (`test.paired_end.sorted.bam` and `genome.fasta`): with and without `-r`, and a stub test. The snapshot takes the BAM's read MD5 and only the stable JSON fields, because the summary records the command line and input path.

### 3.2 MultiQC module `pmdshodh`

It finds summaries by their content (`"tool": {"name": "pmdshodh"`) and reads `format_version` 1 to 4:
- **General statistics:** reads scored, percentage kept (for `filter`), and the 5′ first-base C→T rate.
- **Plots:** the score histogram, damage by position (5′ C→T; 3′ G→A for ds libraries or C→T for ss), and skipped reads by reason.

MultiQC computes rates from the summary's counts, as `pmdshodh report` does, and takes nothing from the HTML report. The module needs a sample file in MultiQC's test-data repository as well, a second PR.

### 3.3 nf-core/eager

Once the module is in nf-core/modules: an issue or discussion on nf-core/eager, proposing pmdshodh as an option next to PMDtools in the `dev` pipeline. It would cite the validation (README, docs/audit-0.3.md) and the speed (docs/benchmarks.md). Their maintainers decide whether, and how, it goes in.

## 4. Order

1. **Local tools:** a new conda env in WSL with Nextflow, nf-core tools, nf-test and MultiQC; Docker for containers.
2. **nf-core module:** build it in a local clone of a fork of nf-core/modules; run `nf-core modules lint` and nf-test locally, and compare its output with running pmdshodh directly on the same files.
3. **MultiQC module:** build it in a local clone of a fork of MultiQC, on summaries from the eager test BAMs and the audit libraries, with all four `format_version`s.
4. **Upstream:** the nf-core module PR, then the MultiQC PRs, then the eager proposal. Each is posted from your GitHub account, only after your OK on the final text.

Bioconda needs nothing from us unless the bot's bump stalls. The module pins whichever container is newest when its PR opens.

## 5. Decisions needed (my recommendation first)

| # | Decision | Options |
|---|---|---|
| I1 | Which module first | **nf-core module** (smaller, and eager consumes modules) / MultiQC first / both at once |
| I2 | Module set | **`pmdshodh/filter` only**, as the roadmap says / also `pmdshodh/damage` for QC without filtering |
| I3 | Threshold input | **`val threshold`, as `pmdtools/filter`**, for a drop-in swap / in `ext.args` only |
| I4 | Reference input | **Optional** (`[]` without it; reads then need MD, and CRAM needs it) / required, as `pmdtools/filter` requires it for `calmd` |
| I5 | Announcing the work | **Open the nf-core "new module" issue when the module passes lint and tests locally** / open it now to claim the name / no issue, PR only |
| I6 | Local tools | **A new conda env `nfcore` in WSL** (nothing outside WSL changes) / something else you prefer |

## 6. What needs your OK at each step

- Installing the local tools (I6): local only, but it downloads Nextflow and its plugins.
- Forking nf-core/modules and MultiQC on your GitHub account.
- Each upstream issue, PR and comment, with the final text shown to you first.
- Any message to the eager or MultiQC maintainers.

## 7. Progress

Both modules are built locally, in clones under WSL `/root/nfcore`, committed there and not pushed.

**nf-core `pmdshodh/filter`** (2026-11-23; nf-core/tools 4.2.2, nf-test 0.9.5, Nextflow 26.10.7, Docker):
- `nf-core modules lint` passes all 67 checks, with no warnings.
- nf-test passes 4 of 4 tests, on two runs with stable snapshots: BAM with MD tags, BAM with FASTA, CRAM with FASTA in single-stranded UDG-half mode, and a stub. The test data is nf-core's `test.paired_end.sorted.bam` and `genome.fasta`.
- The module's BAM records are identical to running `pmdshodh filter` directly with the same container, and the read counts match.
- It pins Bioconda 0.3.1, the newest Bioconda has. The bot's 0.3.3 bump (bioconda/bioconda-recipes#72119) passes its checks and waits for review. The module PR waits for it, so it can pin 0.3.3 (decided 2026-11-24).
- The Seqera containers that nf-core now expects (`nf-core modules containers create`) are built only just before the PR, with your OK.

**MultiQC `pmdshodh`** (2026-11-24; MultiQC 1.35 development version):
- It finds summaries by their second line, `"command": "…pmdshodh …"`, and then checks `tool.name`. The tool name comes after the damage tables and is past MultiQC's 1,000-line search limit when `--range` is above about 40.
- Samples are named after the summary's `input` file. A second summary of the same input is named after its own file, as `pmdshodh report` does.
- **General Statistics:** 5′ first-base C>T, % kept by the threshold, and reads scored (hidden by default).
- **Sections:** read counts (kept, removed, scored, skipped by reason), damage by position (5′ C>T, and 3′ G>A for ds or C>T for ss, computed from the counts as `report` does), the PMD score distribution, and read lengths (scored and kept).
- `--strict` runs without warnings on summaries written by pmdshodh 0.1.0, 0.2.4, 0.3.1 and 0.3.3 (`format_version` 1 to 4), from `filter`, `score` and `damage`, ds and ss.
- 6 unit tests, `ruff` and `mypy` pass.
- Upstream it also needs sample files in MultiQC's test-data repository, made with relative paths so no local path is published.
