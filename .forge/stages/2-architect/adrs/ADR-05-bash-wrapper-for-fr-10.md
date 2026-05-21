# ADR-05: Use Bash for `run_assembly.sh` Wrapper (FR-10)

Date: 2026-05-14
Status: Accepted

## Context

FR-10 calls for a wrapper that runs all four generation steps in order for a given assembly with configurable paths. The wrapper could be implemented as a Bash script, a Python script, a Snakemake rule, or a Makefile target.

## Decision

Use a Bash script (`bin/python/refdata_generation/run_assembly.sh`) accepting `--assembly` and `--date` flags. The script invokes the four Python scripts in dependency order with appropriate flags; auto-detects optional input files (trna, gff3, custom-fasta) and includes only the flags whose source files exist.

## Alternatives Considered

- **Python wrapper:** would add another Python entry point with no clear benefit; argparse plumbing for the four sub-scripts duplicates their CLIs.
- **Snakemake rule:** over-engineered for a single-assembly one-shot run; the user already has Snakemake available for the actual pipeline, not for reference generation.
- **Makefile target:** unusual for this codebase; would require dependency rules duplicating what we already encode in the task DAG.

## Consequences

- Positive: matches the style of existing `wf/eCLIP_repelement_SE` Bash launchers; no additional Python complexity; `set -euo pipefail` and trivial argument parsing.
- Negative: argument parsing in Bash is more verbose than argparse. The wrapper accepts only two flags (`--assembly`, `--date`), so this is acceptable.

## References

- Architecture plan section: §4.7 (Wrapper Script)
- Requirements: FR-10
- Existing Bash launchers: `wf/eCLIP_repelement_SE`, `wf/eCLIP_repelement_PE`
