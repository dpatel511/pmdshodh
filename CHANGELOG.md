# Changelog

All notable changes to pmdshodh are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Changes to the JSON run summary are listed with its `format_version`.

## [Unreleased]

Nothing has been released yet. The first release is planned as v0.1.0 (see [docs/roadmap.md](docs/roadmap.md)).

### Added
- `filter` command: keep reads by PMD score and write SAM or BAM, with an `@PG` header line and an optional score tag.
- `score` command: per-read PMD scores as TSV.
- `damage` command: JSON run summary only.
- PMD score following the likelihood model of Skoglund et al. (2014), with PMDtools-compatible options for double-stranded and single-stranded libraries and CpG/UDG-treated data.
- `--pmdtools-compat` to reproduce PMDtools 0.60 exactly, including its known bugs, for comparisons.
- Multi-threaded scoring and BGZF compression; output is identical for any thread count.
- JSON run summary (`format_version` 1) with read counts, skipped reads by reason, a score histogram and damage tables.
- Parity tests against PMDtools on synthetic and real reads.
- README, design notes, benchmark plan and roadmap.

### Changed
- `--ss` scores each C site once, following the paper's single-stranded model (SI Eq. S1). PMDtools 0.60 adds the double-stranded term as well, which roughly doubles ss scores and makes them depend on which end the damage is at. Default ss scores are therefore lower than PMDtools'; `--pmdtools-compat` keeps PMDtools' scores.
- `filter --upper-threshold` no longer keeps reads with no informative sites; they carry no evidence of being undamaged. They are still counted as `no_informative_sites`. `--pmdtools-compat` keeps them, as PMDtools does.
- Faster scoring, with byte-identical output: double-stranded site terms are precomputed once per run, damage counting no longer copies each read, and BAM records are read directly instead of through the generic record interface. On the first 2,000,000 reads of ENA ERR1892229 (i5-13500H laptop, WSL, medians of 7 interleaved runs), `score -@ 1` went from 11.05 s to 8.19 s, `filter -@ 1` from 10.71 s to 7.29 s and `filter -@ 8` from 3.00 s to 2.18 s.
- With `-@` above 1, reading, scoring and output now overlap (a reader thread and a scoring thread run alongside output), so multi-threaded runs scale further; `-@ 1` is unchanged. Output is byte-identical; peak memory rises from 12–22 MB to 21–33 MB (at most three batches of 10,000 records in flight). On the first 2,000,000 reads of ENA ERR1892229 (i5-13500H laptop, WSL, medians of 5 interleaved runs, monotonic clock), `score -@ 8` went from 3.08 s to 1.75 s and `filter -t 3 -@ 8` from 2.03 s to 1.45 s.

### Fixed
- `score -o -` writes to stdout; it used to create a file named `-`.
- `--pmdtools-compat` now also reproduces PMDtools 0.60 for lowercase reference bases in MD (never counted as C or G) and for reverse-strand reads whose MD contains characters PMDtools drops when reverse-complementing (e.g. lowercase `n`). Default scoring is unchanged; both are listed as departures in docs/design.md.
