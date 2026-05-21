# Security Requirements

## Scope

This is a bioinformatics reference data generation task running on a secure HPC cluster (TSCC). There is no network-facing component, no authentication layer, and no user data handling. Security concerns are minimal and focus on data integrity and execution environment.

## Data Integrity

- **No data modification outside write scope:** Scripts may only write to `examples/inputs/mm10/` and `examples/inputs/mm39/`. The hg38 reference files must not be overwritten (they serve as ground truth).
- **Reproducibility:** Scripts must be deterministic given the same inputs. Sort order in output files must be stable.
- **Verification gate:** hg38 reproduction must pass ≥99% similarity before mm10/mm39 outputs are trusted.

## Access Control

- Running on TSCC under account `bay001` / group `yeo-group`
- Source files are read-only group files; scripts run with standard user permissions
- No sudo or privileged operations required

## Dependency Trust

- Python packages used (pybedtools, pandas, numpy) are sourced from conda/pip with pinned versions in `workflow/envs/dropin.yaml`
- Perl scripts are in-repo and reviewed
- `module load ecliprepmap/1.0.0` loads a trusted internal module

## N/A Items

The following are not applicable to this task:
- Authentication / authorization
- Encryption at rest or in transit
- PII / sensitive data handling
- API keys or secrets management
- Input sanitization (all inputs are controlled bioinformatics files from trusted sources)
