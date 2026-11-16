# Changelog

All notable changes to pmdshodh are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Changes to the JSON run summary are listed with its `format_version`.

## [Unreleased]

### Added
- `-r, --reference FILE`: take reference bases from an indexed FASTA (`samtools faidx`; bgzip with `.gzi` too) instead of MD tags, so reads need no MD tag. MD tags that are present are compared with the reference; disagreeing (usually stale) tags are reported on stderr and counted in the summary, and those reads are scored from the reference. Every header sequence must be in the index with the same length; this is checked before any record is read. With `--pmdtools-compat`, the reference bases go through an MD tag written the way `samtools calmd` writes it (reference bases in uppercase, soft-masked ones included), so scores match PMDtools run on calmd output.
- Real-data checks of `-r`: samtools calmd and PMDtools 0.60 under WSL Ubuntu 24.04, every read compared.
  - BER11_L_NOuser (14,006,183 reads, soft-masked equCab2) and Motala3 (23,622,338 reads, human_g1k_v37): `-r` on the original BAM against the MD path on calmd output gives 0 reads differing in score or site count in default mode on both. In `--pmdtools-compat`, Motala3 also has 0 differing reads.
  - PMDtools on calmd output against `--pmdtools-compat -r`: Motala3 has 0 reads above 1e-6 (largest difference 4.97e-14). BER11 has 339 above 1e-6. All 339 have an MD tag in the input BAM with lowercase reference bases; calmd left those tags as they were, and PMDtools never counts a lowercase base as C or G. With MD removed from those reads (483 have lowercase MD) and recomputed by calmd, PMDtools and `--pmdtools-compat -r` agree on all 483 (largest difference 2.29e-14).
  - Motala3 has 107 reads whose MD tag disagrees with the reference.
- CRAM input for `filter`, `score` and `damage`, detected from the file contents and decoded with `-r` (required). Each reference sequence is checked against the header's `M5` checksum when it is loaded (and single-reference slices against their own MD5), so a different reference is an input error. Multi-reference slices carry no MD5 of their own; before the `M5` check, a wrong reference decoded them silently into wrong read bases. Containers are decoded on the thread pool; output is identical for any thread count. `filter` still writes SAM or BAM.
  - Kept records are written as stored in the CRAM: MD and NM tags that htslib left out (it regenerates them on decoding) are not added back.
  - CRAM decoding loads whole reference sequences: with `score -@ 8 -r`, peak RSS was 1,181 MB on BER11_L_NOuser.cram and 1,247 MB on Motala3.cram, against 110 MB and 105 MB for the same BAMs (median of 3 runs; MB = 1024 KiB, as in docs/benchmarks.md). The `M5` check added about 19% to CRAM wall time (median of 3: 35.6 → 42.3 s and 54.5 → 64.7 s at `-@ 8`).
- Real-data check of CRAM: samtools 1.25 CRAMs of BER11_L_NOuser and Motala3 give byte-identical `score -r` output to their BAMs (14,006,183 and 23,622,338 reads; default and `--pmdtools-compat` at `-@ 8`, default at `-@ 1`). `filter -r -t 3` keeps 1,029,380 and 301,706 reads, identical in every field and tag to `samtools view --input-fmt-option decode_md=0` of the same CRAMs. Unit tests decode CRAM 3.0 and 3.1 files written by samtools, including one with a multi-reference slice, a proper pair and an unmapped read, and an empty CRAM.
- JSON summary `format_version` 2 also records `input_format` (`sam`, `bam` or `cram`). For CRAM, the `reference.md_tags` counts cover only MD tags stored in the file. Version 2 has not been released, so this needs no further bump.
- JSON summary `format_version` 2: a `reference` section with the `-r` path, whether reference bases came from MD tags or the FASTA, and (with `-r`) how many scored reads had MD tags that agree, disagree or are absent.
- PyPI packaging: wheels that carry the `pmdshodh` binary (Linux x86_64 and aarch64 with glibc 2.17 or newer, macOS arm64 and x86_64, Windows x86_64), built and published by a workflow with PyPI trusted publishing. Version 0.1.0 was published to PyPI on 2026-11-15, built from the `v0.1.0` tag.

## [0.1.0] - 2026-11-14

First release.

### Added
- `filter` command: keep reads by PMD score and write SAM or BAM, with an `@PG` header line and an optional score tag.
- `score` command: per-read PMD scores as TSV.
- `damage` command: JSON run summary only.
- PMD score following the likelihood model of Skoglund et al. (2014), with PMDtools-compatible options for double-stranded and single-stranded libraries and CpG/UDG-treated data.
- `--pmdtools-compat` to reproduce PMDtools 0.60 exactly, including its known bugs, for comparisons.
- Multi-threaded scoring and BGZF compression; output is identical for any thread count.
- JSON run summary (`format_version` 1) with read counts, skipped reads by reason, a score histogram and damage tables.
- Parity tests against PMDtools on synthetic and real reads.
- Prebuilt binaries on GitHub Releases for Linux x86_64 and aarch64 (glibc 2.17 or newer), macOS arm64 and x86_64, and Windows x86_64, with SHA-256 checksums. Scores agree across these platforms to within 2e-15: Linux and macOS give byte-identical output on the test fixture, and Windows differs in the last digits.
- Minimum supported Rust version 1.91 (the minimum of noodles 0.119), declared in `Cargo.toml` and tested in CI. The README used to say 1.85, which was never true.
- README, design notes, benchmark plan and roadmap.

### Changed
- `--ss` scores each C site once, following the paper's single-stranded model (SI Eq. S1). PMDtools 0.60 adds the double-stranded term as well, which roughly doubles ss scores and makes them depend on which end the damage is at. Default ss scores are therefore lower than PMDtools'; `--pmdtools-compat` keeps PMDtools' scores.
- `filter --upper-threshold` no longer keeps reads with no informative sites; they carry no evidence of being undamaged. They are still counted as `no_informative_sites`. `--pmdtools-compat` keeps them, as PMDtools does.
- Faster scoring, with byte-identical output: double-stranded site terms are precomputed once per run, damage counting no longer copies each read, and BAM records are read directly instead of through the generic record interface. On the first 2,000,000 reads of ENA ERR1892229 (i5-13500H laptop, WSL, medians of 7 interleaved runs), `score -@ 1` went from 11.05 s to 8.19 s, `filter -@ 1` from 10.71 s to 7.29 s and `filter -@ 8` from 3.00 s to 2.18 s.
- With `-@` above 1, reading, scoring and output now overlap (a reader thread and a scoring thread run alongside output), so multi-threaded runs scale further; `-@ 1` is unchanged. Output is byte-identical; peak memory rises from 12–22 MB to 21–33 MB (at most three batches of 10,000 records in flight). On the first 2,000,000 reads of ENA ERR1892229 (i5-13500H laptop, WSL, medians of 5 interleaved runs, monotonic clock), `score -@ 8` went from 3.08 s to 1.75 s and `filter -t 3 -@ 8` from 2.03 s to 1.45 s.
- `-@ N` now means at most N threads doing work at once, so pipelines can pass their CPU allocation.
  - Before, scoring, BGZF decompression and BGZF compression each had their own pool of N threads, and `-@ 8` used about 9 cores.
  - They now share one pool. The new `src/bgzf.rs` compresses exactly as noodles-bgzf did, so output bytes are unchanged.
  - At the same `-@`, runs are slower than before because they use fewer cores. For the same number of cores they are as fast or faster.
  - Measured on the first 2,000,000 reads of ENA ERR1892229 (i5-13500H laptop, WSL, medians of 5 interleaved runs, monotonic clock):
    - `score -@ 8`: 1.71 s on 6.3 cores (before: 1.65 s on 7.9);
    - `filter -t 3 -@ 8`: 1.45 s on 7.1 cores (before: 1.38 s on 9.4);
    - `filter -t 3 -@ 12`: 1.25 s on 9.6 cores (before: 1.33 s on 10.2);
    - `filter -t 3 -@ 1`: 7.58 s → 8.43 s (it used to take 1.2 cores).
  - Peak memory at `-@ 8` is 33–38 MB, because up to four batches are in memory.
- Errors for MD tags that disagree with the alignment now suggest recomputing MD with `samtools calmd`. Such tags are usually stale after realignment.
- `score` and `filter` count damage tables only when `--summary` is given; the `damage` command always counts them. Output is byte-identical. On the first 2,000,000 reads of BER11_L_NOuser (i5-13500H laptop, WSL, medians of 5 interleaved runs, monotonic clock):
  - `score -@ 1`: 8.77 s → 6.18 s (metaDMG `pmd -@ 1`: 8.23 s; before the shared thread pool: 8.68 s);
  - `filter -t 3 -@ 1`: 8.77 s → 6.17 s (before the shared thread pool: 7.95 s);
  - `score -@ 8`: 1.60 s → 1.26 s; `filter -t 3 -@ 8`: 1.42 s → 1.04 s;
  - with `--summary`, unchanged (1.43 s for `filter -t 3 -@ 8`).
- Each command's `--help` groups its options under headings (Output, Filtering, Summary, Scoring model, Read filters, Performance) and ends with examples. A test runs every example. Usage lines say `pmdshodh` on every platform.
- `--range 0` is rejected as invalid input (exit 2) before any input is read; it used to give empty damage tables.

### Fixed
- `--pmdtools-compat` with several threads no longer slows down 2–5× from contention on glibc's allocator lock. The shared thread pool above made this worse, and the fix sizes compat's per-read buffers up front instead of growing them. On BER11_L_NOuser (14.0 M reads, 5 alternating runs), `filter --pmdtools-compat -@ 8` took 10.8–11.8 s instead of 23.9–66.6 s. Scores are unchanged.
- `score -o -` writes to stdout; it used to create a file named `-`.
- `--pmdtools-compat` now also reproduces PMDtools 0.60 for lowercase reference bases in MD (never counted as C or G) and for reverse-strand reads whose MD contains characters PMDtools drops when reverse-complementing (e.g. lowercase `n`). Default scoring is unchanged; both are listed as departures in docs/design.md.

[Unreleased]: https://github.com/dpatel511/pmdshodh/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/dpatel511/pmdshodh/releases/tag/v0.1.0
