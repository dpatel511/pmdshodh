# Changelog

All notable changes to pmdshodh are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Changes to the JSON run summary are listed with its `format_version`.

## [Unreleased]

## [0.2.3] - 2026-11-18

`filter` no longer adds a second CG tag to reads with more than 65,535 CIGAR operations read from BAM; kept records are written unchanged. Scores are unchanged. This release also records checks on the output of bwa, bowtie2 and samtools against PMDtools 0.60, and documents a `--pmdtools-compat` limitation on `=` read bases.

### Added
- Checks on the output of real aligners (WSL Ubuntu 24.04; bwa 0.7.19, bowtie2 2.5.5, samtools 1.25, PMDtools 0.60). 100,000 reads of ENA ERR3991018 (single-stranded, no UDG) were aligned to the rCRS (`MT` of human_g1k_v37), single-end and as fully overlapping pairs.
  - The inputs were:
    - `bwa aln` (eager's `-l 1024 -n 0.01`), `bwa mem` (supplementary alignments with hard clips) and `bowtie2` (end-to-end, `--local`, `-k 5`);
    - then `samtools` sort, merge, markdup (single-end and paired), name-sort, `calmd -e`, and CRAM 3.0 and 3.1;
    - and 14 crafted reads: a 66,000 bp read with 66,000 CIGAR operations, 1,500 and 5,000 bp reads, IUPAC codes in reads and reference, an N run, a soft-masked region, `X`/`=` and padding CIGARs, clips next to indels, a secondary read without SEQ, and a 1 bp read.
  - For each input, `score` and `filter` gave the same output at `-@ 1` and `-@ 4`, with no NaN scores, and every output passed `samtools quickcheck`.
  - `--pmdtools-compat` against PMDtools, on the reads PMDtools scores: 0 reads above 1e-6 and identical kept reads on 17 of the 18 aligner outputs, all but `calmd -e` (up to 199,865 reads each, largest difference 4.2e-14), and on the 7 crafted reads PMDtools can read. PMDtools crashed on the reads over 1,000 bp, a known PMDtools bug (docs/design.md §5).
  - The `calmd -e` output differs in compat mode; this is documented instead of reproduced (docs/design.md §5).
  - `-r` against MD: identical scores for bowtie2. For bwa output, about 1% of reads differ (968 of 99,935 single-end), all covering MT:3107. There the rCRS has `N`, which bwa replaces with a random base when indexing, so its MD tags give that base; `-r` does not score the site.

### Changed
- The stderr note on MD tags that disagree with `-r` names both usual causes: MD written before realignment, or against reference `N` bases that the aligner replaced, as bwa does. It used to say "usually stale after realignment".

### Fixed
- `filter`: a read with more than 65,535 CIGAR operations, which BAM stores in a CG tag, gained a second CG tag when read from BAM and written as BAM or SAM (264,008 bytes more for 66,000 operations), so samtools showed CG as an ordinary tag. Kept records are now written unchanged. SAM input and `--tag` were not affected.

## [0.2.2] - 2026-11-18

Fixes from an edge-case audit of the command line; scores are unchanged. Writing output no longer replaces FIFOs, devices or symlinks, and an output that is the input file is refused instead of overwriting it. Three behaviour changes: input with no data at all is an input error (exit 2); a closed stdout (`| head`) ends the run quietly with exit 141, as in samtools; and `--range` (at most 1000) and `-@` (at most 1024) have upper limits.

### Added
- Real-data check of single-stranded scoring, on library UKY001.B0102 (ENA PRJEB37007, single-stranded and not UDG-treated per its ENA record). Its damage pattern fits that: C→T at both ends (25.5% at the 5′ and 17.7% at the 3′ terminal base of the 1240K BAM), G→A at 1.3% or below. PMDtools 0.60 `--ss` was run under WSL Ubuntu 24.04, and every read was compared.
  - `score --library ss --pmdtools-compat` against `--printDS --ss`: 0 reads above 1e-6 on the mtDNA BAM (ERR3991018, 281,906 reads, largest difference 6.2e-14) and on the 1240K-capture BAM (ERR3977849, 14,461,790 reads, 8.2e-14).
  - `filter --library ss --pmdtools-compat --threshold 3` keeps the same reads as `--threshold 3 --ss`, identical including tags: 162,734 and 7,169,568.
  - In default mode, which scores each C site once, `--threshold 3` keeps 117,380 and 5,327,200 reads. All but 2 of those are also kept by PMDtools, and 1,842,370 reads kept by PMDtools are not. Every read with an informative site scores differently from PMDtools; the reads that score the same have none.
  - `filter --library ss` output is identical at `-@ 1` and `-@ 8`.

### Changed
- Input with no data at all (0 bytes, or empty standard input) is an input error (exit 2), as in samtools; it usually means an upstream step failed. Input with only a header is still valid and gives empty output.
- When the reader of standard output closes it early (`pmdshodh score in.bam | head`), pmdshodh stops without a message and exits 141, as samtools does when SIGPIPE ends it. Before, it printed "cannot write output: Broken pipe" and exited 1.
- Options are checked before any input is read:
  - `--pmd-p` + `--pmd-c`, the damage rate at the read end, must be at most 1. Larger values gave NaN scores, which `filter` silently dropped; `--pmd-p 1` with the default `--pmd-c 0.01` did this for 1,081 of 2,912 reads.
  - `--threshold` and `--upper-threshold` must be finite numbers: `-t nan` used to keep nothing.
  - `--range` is at most 1000 (`--range 100000000` was killed for lack of memory), and `-@` at most 1024 (`-@ 100000` tried to start threads for 4 minutes, then failed).
  - A negative value given with a space (`--pmd-p -0.1`) gets the range message instead of "unexpected argument '-0'".
- `filter` output format: `.BAM` and `.SAM` extensions are accepted in any case. Devices have no extension, so `--output-format` decides there as for stdout, BAM by default: `-o /dev/null --summary run.json` keeps only the summary. `--output-format` that contradicts the extension, and `-o x.cram`, are input errors instead of being ignored or reported as an unknown format.
- Error messages name the path given, not the temporary file (`out.bam.pmdshodh-tmp`). Input that is not SAM, BAM or CRAM is reported as such ("it looks like FASTQ", "it looks like FASTA", "it is binary data") instead of "invalid header: invalid record" or "record 1: unexpected EOL", and a BAM cut inside a BGZF block says that the file is probably truncated.

### Fixed
- Output to an existing FIFO or device was renamed over it: `filter -o x.sam` replaced a FIFO named `x.sam` (its reader got nothing), and as root `score -o /dev/full`, `--summary` or `damage -o` to a device replaced the device node with a regular file. A symlink output was replaced by a regular file instead of updating its target (for `score` since 0.2.1, which made `score -o` write by rename; `filter`, `damage` and `--summary` already did this), and `score -o >(...)` failed (since 0.2.1). `-o /dev/stdout` failed for users other than root, and `-o NUL` on Windows failed ("cannot move the file to a different disk drive"). Only regular files are now written under a temporary name and renamed; a symlink is followed, and anything else is written in place.
- An output that is the input file (also through `./`, and on Linux and macOS through a hard link) replaced the input with exit 0: `filter -t 3 same.bam -o same.bam` left 34 of 2,912 reads. This, the reference as output, and one file given for two outputs (`-o x.bam --summary x.bam`) are now input errors.
- `--summary` (and `damage -o`) are created before the first record is read. A summary in a missing directory used to fail only after the whole run, with the BAM already written and exit 1.
- `filter`: a kept read whose reference or mate's reference is not an `@SQ` line of the header is an input error that names the read (exit 2). It used to fail with "cannot write output: invalid reference sequence name" (exit 1). `score` and `damage` do not write records and still read SAM without a header (`samtools view in.bam | pmdshodh score -`).
- docs/benchmarks.md said 32 checks passed in the 0.2.1 round: `check.py` ran 24 checks, all passed, and printed 8 informational lines.
- `docs/benchmarks/check.py` fails on a NaN or infinite score; it used to count it as a small difference. No NaN or infinite score was found in any earlier benchmark or parity output.

## [0.2.1] - 2026-11-17

Bug fixes and speed-ups. Two inputs that used to be read silently are now input errors (exit 2): BGZF input that ends without its end-of-file marker (usually a truncated file), and reads with `=` bases when no reference is given with `-r`. Default scores change only for reads with `=` bases scored with `-r`.

### Added
- Benchmarks of `--library ss`, `-r` and CRAM input, and a new benchmark round on this release: [docs/benchmarks.md](docs/benchmarks.md).

### Changed
- Read bases written as `=` (same as the reference, as `samtools calmd -e` writes them) are scored in default mode. With `-r` they take their reference base; without `-r` a read with `=` bases is an input error naming `-r`, because MD does not give their reference bases. PMDtools never counts `=` as C or G: a read with nine `=` at reference C scored 5.03 from 1 site instead of 4.22 from 10. `--pmdtools-compat` keeps that. Kept records are written with their `=` bases.
- Damage tables (`--summary`, `damage`) count only positions within `--range` of the read ends instead of looking at every base. Output is identical. `score -@ 1 --summary` on 2 M BER11_L_NOuser reads: median 9.16 → 8.33 s (5 alternating runs each; WSL Ubuntu 24.04, i5-13500H).
- CRAM input: reference sequences are checked against the header `M5` on the thread pool instead of the reading thread. `score -@ 8 -r` on BER11_L_NOuser.cram: median 45.9 → 38.8 s (5 alternating runs each; WSL Ubuntu 24.04, i5-13500H). A wrong reference is still an input error before any record of the batch is used.
- Single-stranded scoring (`--library ss`) is faster. The damage, error and null terms are computed once per run instead of at every site (bit-identical), and reads are no longer copied into molecule orientation. `score -@ 1 --library ss` on 2 M BER11_L_NOuser reads: median 6.98 → 5.53 s (5 alternating runs each; WSL Ubuntu 24.04, i5-13500H).

### Fixed
- BGZF input (BAM, bgzipped SAM) that ends without the end-of-file marker is now an input error (exit 2), and `filter` leaves no output file. Before, a file cut at a block boundary was read silently: on a samtools-written BAM cut that way, 1,261 of 1,500 reads were scored with exit 0. samtools only warns in this case. CRAM input cut before its end-of-file container was already an error.
- `score -o FILE` writes under a temporary name and renames on success, as `filter` does, so a failed run no longer leaves a partial table.
- gzip-compressed input that is not BGZF (as `gzip` writes it) is reported as such, with a hint to recompress with bgzip, instead of "invalid BGZF header".
- CRAM input: the cache of whole reference sequences is also cleared past its limit (256 Mi bases) when a batch has only multi-reference containers, as samtools writes for many small contigs. Before, such input kept every sequence it had used, up to the whole reference. On a synthetic 1 Gb reference of 10,000 contigs (100 multi-reference containers, samtools 1.25), peak RSS of `score -@ 4 -r` fell from 1,010 MB to 392 MB, with identical scores.

## [0.2.0] - 2026-11-16

Reference-based scoring (`-r`) and CRAM input: `pmdshodh filter -r ref.fa in.bam` replaces `samtools calmd | pmdtools | samtools view`.

### Added
- `-r, --reference FILE`: take reference bases from an indexed FASTA (`samtools faidx`; bgzip with `.gzi` too) instead of MD tags, so reads need no MD tag. MD tags that are present are compared with the reference; disagreeing (usually stale) tags are reported on stderr and counted in the summary, and those reads are scored from the reference. Every header sequence must be in the index with the same length; this is checked before any record is read. With `--pmdtools-compat`, the reference bases go through an MD tag written the way `samtools calmd` writes it (reference bases in uppercase, soft-masked ones included), so scores match PMDtools run on calmd output.
- CRAM input for `filter`, `score` and `damage`, detected from the file contents and decoded with `-r` (required). Each reference sequence is checked against the header's `M5` checksum when it is loaded (and single-reference slices against their own MD5), so a different reference is an input error. Containers are decoded on the thread pool; output is identical for any thread count. `filter` still writes SAM or BAM.
  - Kept records are written as stored in the CRAM: MD and NM tags that htslib left out (it regenerates them on decoding) are not added back.
  - CRAM decoding loads whole reference sequences: with `score -@ 8 -r`, peak RSS was 1,181 MB on BER11_L_NOuser.cram and 1,247 MB on Motala3.cram, against 110 MB and 105 MB for the same BAMs (median of 3 runs; MB = 1024 KiB, as in docs/benchmarks.md). The `M5` check added about 19% to CRAM wall time (median of 3: 35.6 → 42.3 s and 54.5 → 64.7 s at `-@ 8`).
- PyPI packaging: wheels that carry the `pmdshodh` binary (Linux x86_64 and aarch64 with glibc 2.17 or newer, macOS arm64 and x86_64, Windows x86_64), built and published by a workflow with PyPI trusted publishing. 0.1.0 was published to PyPI on 2026-11-15 from the `v0.1.0` tag.
- Real-data checks of `-r`: samtools calmd and PMDtools 0.60 under WSL Ubuntu 24.04, every read compared.
  - BER11_L_NOuser (14,006,183 reads, soft-masked equCab2) and Motala3 (23,622,338 reads, human_g1k_v37): `-r` on the original BAM against the MD path on calmd output gives 0 reads differing in score or site count in default mode on both. In `--pmdtools-compat`, Motala3 also has 0 differing reads.
  - PMDtools on calmd output against `--pmdtools-compat -r`: Motala3 has 0 reads above 1e-6 (largest difference 4.97e-14). BER11 has 339 above 1e-6. All 339 have an MD tag in the input BAM with lowercase reference bases; calmd left those tags as they were, and PMDtools never counts a lowercase base as C or G. With MD removed from those reads (483 have lowercase MD) and recomputed by calmd, PMDtools and `--pmdtools-compat -r` agree on all 483 (largest difference 2.29e-14).
  - Motala3 has 107 reads whose MD tag disagrees with the reference.
- Real-data check of CRAM: samtools 1.25 CRAMs of BER11_L_NOuser and Motala3 give byte-identical `score -r` output to their BAMs (14,006,183 and 23,622,338 reads; default and `--pmdtools-compat` at `-@ 8`, default at `-@ 1`). `filter -r -t 3` keeps 1,029,380 and 301,706 reads, identical in every field and tag to `samtools view --input-fmt-option decode_md=0` of the same CRAMs. Unit tests decode CRAM 3.0 and 3.1 files written by samtools, including one with a multi-reference slice, a proper pair and an unmapped read, and an empty CRAM.

### Changed
- JSON run summary `format_version` 2 (was 1). New fields: `input_format` (`sam`, `bam` or `cram`), and a `reference` section with the `-r` path, whether reference bases came from MD tags or the FASTA, and (with `-r`) how many scored reads had MD tags that agree, disagree or are absent. For CRAM, the MD tag counts cover only MD tags stored in the file. Fields from version 1 are unchanged.
- The error for a mapped read without an MD tag now suggests `-r` as well as `samtools calmd`.

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

[Unreleased]: https://github.com/dpatel511/pmdshodh/compare/v0.2.3...HEAD
[0.2.3]: https://github.com/dpatel511/pmdshodh/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/dpatel511/pmdshodh/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/dpatel511/pmdshodh/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/dpatel511/pmdshodh/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/dpatel511/pmdshodh/releases/tag/v0.1.0
