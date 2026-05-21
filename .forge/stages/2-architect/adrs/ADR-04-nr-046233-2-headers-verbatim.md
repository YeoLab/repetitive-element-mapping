# ADR-04: Preserve NR_046233.2.fasta Headers Verbatim

Date: 2026-05-14
Status: Accepted

## Context

The mouse `NR_046233.2.fasta` is a custom RefSeq rRNA FASTA file with pre-curated headers. The Perl parse scripts use a `rRNA_extra_hash` that maps `RNA28S` / `RNA18S` / `RNA5-8S` to `RNA45S`. The architect must decide whether the new `generate_bowtie2_index.py` script should:

1. Append the FASTA records verbatim (preserve headers).
2. Split the file by sub-region and rename headers (e.g., `NR_046233.2-18S`, `NR_046233.2-28S`).

## Decision

Append the FASTA records verbatim. Do not split or rename. If the FASTA contains multiple records (one per sub-region) the curation has already been done; if it contains one record the implementer should not synthetically split it.

## Alternatives Considered

- **Split into 18S/28S/45S with renamed headers:** rejected because we don't have authoritative coordinates for the sub-regions and the file is short — its content is presumably the curated form expected by downstream code.
- **Conditional rename based on content inspection:** rejected as added complexity with unclear payoff. The Perl rRNA_extra_hash keys are matched by family label in MASTER_FILELIST col4, not by FASTA header.

## Consequences

- Positive: zero risk of corrupting pre-curated rRNA reference data; simplest implementation.
- Negative: if the FASTA's headers do not naturally line up with downstream expectations, the user may have to manually re-header it. The implementer can document this in `README.md` and the MASTER_FILELIST custom-fasta handling logic still groups by family label.
- Risks: if a single-record FASTA needs to be split, this decision delays surfacing the problem until pipeline integration (T-12).

## References

- Architecture plan section: §4.4 (Step 2 algorithm step 7), §6 (Open Technical Questions)
- Assumptions: A-07
- Edge case: Edge Case 5 (NR_046233.2 rRNA custom FASTA)
