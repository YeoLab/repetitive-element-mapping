# Architecture Decisions Summary

**Stage:** 2-architect
**Prepared:** 2026-05-14
**Full architecture plan:** `architecture-plan.md`
**Full ADRs:** `adrs/`

---

## Top 5 Decisions

### 1. Pure Python + pybedtools (Approach A) — ADR-01

**Decision:** Four independent Python scripts using only the existing `workflow/envs/dropin.yaml` conda env. No new dependencies.

**Rejected:** HTSeq/pyranges (would require new packages, violates NFR-01); Bash/awk (fragile GTF parsing, locale-dependent sort).

**Risk:** pybedtools `BedTool.sequence()` strand/splice semantics must be validated against hg38 ground truth.

---

### 2. New scripts live in `bin/python/refdata_generation/` — ADR-02

**Decision:** Create a new Python package subdirectory for the four scripts + shared module + wrapper.

**Rejected:** Flat layout in `bin/python/` would mix new scripts with Perl-compat shims.

**Risk:** None significant.

---

### 3. Reproduction-test-first workflow — ADR-03

**Decision:** Tasks T-03, T-05, T-07, T-09 are blocking gates against hg38 references. Mouse generation (T-10, T-11) only proceeds after all four gates pass.

**Rejected:** Generate mouse first and validate later (would publish potentially-broken references and make debugging harder).

**Risk:** If the original hg38 generator was non-deterministic, AC-01 (100% identical) may be unreachable — escalation path documented (R-07).

---

### 4. Append `NR_046233.2.fasta` headers verbatim — ADR-04

**Decision:** Do not split or re-header the custom rRNA FASTA. The file is curated; preserve as-is.

**Rejected:** Synthetic 18S/28S/45S splitting without authoritative sub-region coordinates.

**Risk:** If headers don't match downstream conventions, surfaces in T-12 pipeline integration; user re-headers manually.

---

### 5. Bash wrapper for FR-10 — ADR-05

**Decision:** `bin/python/refdata_generation/run_assembly.sh` invokes the four Python scripts in order. Auto-detects optional inputs.

**Rejected:** Python wrapper (redundant argparse); Snakemake rule (over-engineered); Makefile (no precedent in repo).

**Risk:** None significant; matches existing `wf/eCLIP_repelement_SE` Bash launcher convention.

---

## Reading Order for Reviewers

1. `architecture-plan.md` §1-3 — mission and approaches
2. This document — high-level decision summary
3. `threat-model.md` — security boundaries (light)
4. `tasks/T-01.md` through `T-12.md` — task-level decomposition
5. `adrs/` — full ADR records for each decision
