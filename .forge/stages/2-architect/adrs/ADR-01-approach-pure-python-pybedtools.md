# ADR-01: Use Pure Python + pybedtools (Approach A)

Date: 2026-05-14
Status: Accepted

## Context

We must implement four reference-data generation scripts. Three architectural approaches were considered: (A) Pure Python with pybedtools, (B) HTSeq/pyranges-based single module, (C) Bash/awk pipelines.

## Decision

Adopt Approach A: four independent Python scripts that use only `argparse`, `gzip`, `subprocess`, `pybedtools`, and `pandas` — all already present in `workflow/envs/dropin.yaml`.

## Alternatives Considered

- **Approach B (HTSeq/pyranges):** rejected because `HTSeq` and `pyranges` are not in `dropin.yaml` — adding them violates NFR-01.
- **Approach C (Bash/awk pipelines):** rejected because GTF attribute-field parsing in awk is fragile, locale-dependent `sort` violates the deterministic-output guarantee (FR-09), and the existing codebase has no Bash precedents for non-trivial data processing.

## Consequences

- Positive: zero new dependencies; matches existing `bin/python/` style; each script independently testable; decoupled reproduction tests per output file.
- Negative: must hand-roll GTF attribute parsing (a ~10-line function, low complexity).
- Risks: pybedtools `BedTool.sequence()` semantics for strand and spliced-extraction must be carefully validated against hg38 ground truth.

## References

- Architecture plan section: §3 (Approaches Evaluated)
- Assumptions: A-09 (input files trusted), A-10 (samtools faidx available)
