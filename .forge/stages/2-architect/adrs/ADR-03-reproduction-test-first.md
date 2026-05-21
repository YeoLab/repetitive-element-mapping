# ADR-03: Reproduction-Test-First Workflow

Date: 2026-05-14
Status: Accepted

## Context

AC-01 requires the new `parsed_ucsc_tableformat` script to reproduce the hg38 reference file 100% line-for-line identically. ACs 04, 07, 10 require ≥99% similarity for the other three reference files. The original generator logic is not fully documented; sort order, CDS-derivation details, and disambiguation suffixes for duplicate names must be reverse-engineered.

## Decision

Each generation script has a sibling reproduction-test task (T-03, T-05, T-07, T-09) that runs the script against hg38 inputs and compares to the existing hg38 reference. **No mouse-reference generation (T-10, T-11) starts until all four reproduction gates pass.**

## Alternatives Considered

- **Generate mouse first, validate later:** rejected because debugging diffs is far harder when the only reference is the hg38 output that the script was meant to reproduce. Producing mouse outputs first risks publishing broken reference files.
- **Skip the 100% identity gate, target ≥99% for all four:** rejected because the parsed_ucsc_tableformat is a deterministic transformation; ≥99% leaves a 2,500-row tolerance window that is too wide to detect off-by-one errors.

## Consequences

- Positive: forces the script logic to be discovered against an authoritative ground truth before any new data is generated.
- Negative: adds four explicit reproduction tasks (T-03/05/07/09); ~2-4 iteration cycles per task.
- Risks: if the original hg38 generator was non-deterministic (e.g., relied on Python hash iteration), 100% identity may be unreachable. R-07 in the risk register; ADR-03 mandates escalation after 3 failed iterations.

## References

- Architecture plan section: §7 (Decision Register), §11 (Verification Strategy)
- Risk register: R-01, R-07
- Acceptance criteria: AC-01, AC-04, AC-07, AC-10
