# 07 - Security Requirements

## N/A — HPC Batch Pipeline

This pipeline runs as a batch job on a SLURM HPC cluster. It is not a web service, does not expose any network interfaces, and does not handle user credentials or sensitive personal data.

## File Permissions

- Input files are read-only (user owns them)
- Output files are written to user-specified directories
- No setuid, no privileged operations
- Snakemake lock files (`.snakemake/`) are in the working directory

## Secrets / Credentials

- No API keys, passwords, or tokens are used
- All file paths in config YAML — no hardcoded secrets in Snakefile or rules
- Commit: do not include paths to private/restricted data in example YAMLs committed to git; use placeholder paths or paths within the repository

## Data Sensitivity

- eCLIP sequencing data is research data, not PHI/PII
- No special data handling requirements
