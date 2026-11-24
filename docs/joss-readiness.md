# JOSS readiness check (2026-11-24)

This file goes through the Journal of Open Source Software's submission requirements, pre-review screening gates, review checklist and paper requirements, and matches each against pmdshodh. The JOSS documentation was read from its source repository, openjournals/joss `docs/` (`submitting.md`, `review_criteria.md`, `review_checklist.md`, `paper.md`), last changed 2026-11-24. joss.readthedocs.io showed a bot check.

## 1. Verdict

**Not eligible yet. JOSS would desk-reject pmdshodh today, on timing and evidence of use rather than on quality.**

Three of the four pre-review gates fail:
- The repository has been public for 13 days, and JOSS requires more than six months.
- No research use is documented yet.
- Nearly all of the work is concentrated in two weeks.

JOSS calls such a rejection "not yet" rather than "never". The earliest a submission could pass the history gate is **after 2027-05-11**, and only with development continuing through that period and research use to show.

Most of the material for the paper and the review already exists.

## 2. Pre-review screening gates (all four required)

| Gate | JOSS requirement | pmdshodh | Status |
|---|---|---|---|
| 1. Public development history | Public for more than six months before submission, with active development spanning that period. "We run automated checks on commit distribution — a repo dump is not a history." | Public since creation on 2026-11-11 (13 days) | **Not met.** Eligible after 2027-05-11 at the earliest |
| 2. Demonstrated research impact | Evidence that the software is used for research, at least by the developers: papers or preprints citing it, adoption by other groups, or integration into research workflows. "Aspirational statements about future use are not sufficient." | No publication, analysis or outside user yet. The MultiQC module is under review, the nf-core module is not submitted, and there is no eager integration | **Not met** |
| 3. Good open-source practices | For a single author, multiple indicators: a meaningful commit history over time, tagged releases or a changelog, tests and CI, clear documentation, a CONTRIBUTING file, and stated support or governance expectations | 197 commits, 79 PRs, 10 tagged releases, a CHANGELOG, CI and 150 tests, a detailed README and docs. **No CONTRIBUTING file and no support statement** | **Mostly met.** Missing pieces are quick to add |
| 4. Iterative development over time | History showing refinement through use and feedback over time, "not a single burst of commits" | Commits on 12 days between 2026-11-11 and 2026-11-24 | **Not met.** This is the "concentrated window" JOSS names |

Counting in favour, but not required: engagement beyond the authors, such as issues, PRs or discussions by others. pmdshodh has none yet: no outside issues and no stars.

## 3. Submission requirements

| Requirement | pmdshodh | Status |
|---|---|---|
| OSI-approved open source licence | MIT (`LICENSE`) | Met |
| Browsable source, open issue tracker, PRs without approval or payment, clonable without registration | Public GitHub repository with issues enabled | Met |
| An obvious research application | Filtering ancient DNA reads by post-mortem damage, a standard step in aDNA pipelines (nf-core/eager uses PMDtools for it) | Met |
| A major contributor submits, with a GitHub account | Single author | Met |
| The paper does not focus on new research results | The paper would describe the software; the damage-curve study is supporting validation | Met, if written that way |
| `paper.md` and `paper.bib` in the repository | Not written | To do at submission |
| Research software, not a minor utility or thin wrapper | A full reimplementation of the PMD model, with reference-based scoring, CRAM, damage tables, reports and summaries | Met |
| Feature-complete and packaged to community standards | Released binaries, PyPI wheels and Bioconda; Rust crate layout | Met |

## 4. Review checklist

**General checks:**
- **Repository:** met.
- **Licence:** met.
- **Contribution and authorship:** met; single author.
- **Scope and significance:** not yet, because research impact is undemonstrated (gate 2).

**Development history and open-source practice:**
- **Development timeline:** not met. JOSS calls "all or most commits concentrated in the last few weeks before submission" not acceptable.
- **Open development:** met in spirit, since the repository was public from its first day. But six months of history are required.
- **Collaborative effort:** weak. JOSS looks for multiple contributors, or for community use and influence; pmdshodh has a single author and no outside engagement. JOSS's alternative is established research use, which pmdshodh does not have yet.
- **Good practices:** mostly met. Contribution guidelines and support channels are missing.

**Functionality:**
- **Installation:** met (binaries, PyPI, Bioconda).
- **Functionality:** met. The claims are backed by read-by-read parity (README, docs/audit-0.3.md).
- **Performance:** met. The claims carry measurements (docs/benchmarks.md: 32–39× on one thread).

**Documentation:**
- **Statement of need:** mostly met. The README states the problem and the audience; the paper would sharpen it.
- **Installation instructions:** met.
- **Example usage:** met. The README quick start uses real-world commands.
- **Functionality documentation:** met. Commands and options are in the README and `--help`, with design.md behind them.
- **Automated tests:** met. CI on Linux and macOS, unit, integration, property and fuzz tests.
- **Community guidelines** (how to contribute, report issues and get support): **not met.** The README has a Development section but no explicit "how to contribute / report / get help" guidance.

## 5. Paper requirements

The paper is 750–1,750 words in Markdown with YAML metadata, plus `paper.bib`. The required sections, and what already exists for each:

| Section | What JOSS asks | Material we have | Status |
|---|---|---|---|
| Summary | High-level function and purpose for a non-specialist | README overview | Draftable now |
| Statement of need | The problem, the audience, and the relation to other work | README; docs/background-investigation.md | Draftable now |
| State of the field | Comparison with common tools, and a **"build vs. contribute"** justification | docs/literature.md and background-investigation.md cover PMDtools, metaDMG `pmd`, ngsBriggs, DamageProfiler and mapDamage2 | Needs care (below) |
| Software design | Trade-offs, architecture, and why they matter | docs/design.md: Rust and noodles, streaming and threads, default vs `--pmdtools-compat`, reference-based scoring | Draftable now |
| Research impact statement | Realized impact, or credible near-term significance (benchmarks, reproducible materials, community readiness) | Benchmarks, parity validation, the damage-curve study, MultiQC and nf-core modules (pending) | Near-term significance only; realized impact is missing (gate 2) |
| AI usage disclosure | Tools and versions used; where they were used (code, docs, paper); the nature of the assistance; confirmation that the human authors reviewed and validated everything and made the core design decisions | Not written | **Required; needs your decision (below)** |
| References | Complete, in citation syntax | Skoglund et al. 2014, Briggs et al. 2007, the tools compared, datasets | Draftable now |
| Software archive | A link to an archive with a DOI (for example Zenodo) at acceptance | None yet | To do: set up Zenodo archiving |

### Points that need thought

**Build vs. contribute.** JOSS asks why pmdshodh is a new tool rather than a contribution to an existing one. The honest case:
- **PMDtools:** Python, unmaintained, with the bugs documented in the README's compatibility section. Fixing it would not give the speed, reference-based scoring or CRAM support.
- **metaDMG `pmd`:** overlaps the per-read score. It has no licence file, so its code cannot be reused, and it focuses on per-taxon damage, without BAM filtering output (docs/literature.md 5.1).
- **ngsBriggs:** a different model. It needs an external contamination estimate and has no release.
- Contact with metaDMG's developers is still undecided (roadmap). A reviewer may ask whether contributing BAM output to metaDMG was considered.

**AI usage disclosure.** pmdshodh's code, tests, docs and these analyses were written with an AI coding assistant (Claude Code), under your direction and review. JOSS requires the paper to disclose this:
- the tool and model versions;
- where it was used;
- the nature of the help;
- your confirmation that you reviewed and validated the outputs and made the core design decisions.

This is separate from your rule of no AI attribution in commits and PRs. JOSS makes the disclosure a condition of submission, and leaving it out can count as an ethics breach.

## 6. What would make a future submission pass

**Do now (cheap, and good practice anyway):**
1. **Add CONTRIBUTING.md and a support statement:** how to report bugs, ask questions and propose changes; issue templates; optionally a code of conduct. This closes the remaining parts of gate 3 and of community guidelines.
2. **Add CITATION.cff and Zenodo archiving of releases.** Each release then gets a DOI and pmdshodh is citable now, without a paper. Turning on the GitHub–Zenodo integration is an account action for you.

**Over the next six months or more (the gates JOSS cannot waive):**
3. **Keep developing openly over time:**
   - track work in public issues and PRs;
   - make releases as fixes and features land (0.3.4 with the CLI fixes is ready);
   - respond to MultiQC, nf-core and user feedback.

   JOSS checks the distribution of commits, so steady activity matters more than volume.
4. **Show research use:**
   - the MultiQC and nf-core modules accepted upstream;
   - pmdshodh offered in nf-core/eager;
   - and above all pmdshodh used in an analysis: yours, a collaborator's, or an outside group's, in a paper, preprint or public workflow.

   Outside issues or PRs would count strongly in its favour.
5. **Decide on contact with the metaDMG developers.** It strengthens the build-vs-contribute answer either way.

**At submission (not before 2027-05-12):**
6. Write `paper.md` and `paper.bib` from the material above, within 750–1,750 words, including the AI usage disclosure.
7. Have a colleague install and try pmdshodh, which JOSS asks for new software.
8. Submit through the JOSS website, with the archive DOI at acceptance.
