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
