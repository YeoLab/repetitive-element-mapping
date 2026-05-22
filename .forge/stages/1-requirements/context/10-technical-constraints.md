# 10 - Technical Constraints

## Environment

- **Cluster:** TSCC2 (San Diego Supercomputer Center)
- **Scheduler:** SLURM
- **SLURM partition:** gold
- **SLURM account:** csd792
- **Snakemake version:** 9.12.0 (pinned — use `conda activate snakemake9`)
- **Snakemake profile:** `profiles/tscc2_snakemake9/` (in repo root)
- **Python:** Snakemake uses Python 3 from the snakemake9 conda env
- **Perl:** MUST use `/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl` (NOT system perl)
- **bowtie2:** 2.2.6 (via `module load bowtie2/2.2.6`) or via conda env

## SLURM Profile

Located at `profiles/tscc2_snakemake9/` — currently empty (directory exists, config not yet written).

Default job requirements:
- Memory: 32 GB (32000 MB)
- CPUs: 1
- Wall time: 8 hours

Per-rule overrides:
- `map_repetitive_elements`: 8 CPUs, 16 GB (per CWL ResourceRequirement)
- `deduplicate`: 32 GB (per CWL ResourceRequirement)
- `getpair`: 1 GB (per CWL ResourceRequirement)

Retry policy: 2 retries; memory scales with attempt (e.g., `attempt * 32000` MB for deduplicate).

## Software Stack

- Conda env activation: `conda activate snakemake9`
- Module loads needed: `module load bowtie2/2.2.6` (if not in conda env)
- Python for fold change script: uses snakemake9 env python or ecliprepmap-0.1.0 env python
- Perl scripts: must be called with full path to perl 5.10.1

## Repository Constraints

- Do NOT modify Perl scripts (unless a script literally cannot be called from Snakemake, and even then ask the user first)
- Snakefile must live at repository root or in `workflow/`
- Rules modularized into: `workflow/rules/SE.smk`, `workflow/rules/PE.smk`, `workflow/rules/common.smk`
- Conda env file: `workflow/envs/dropin.yaml` (Python 3.11, bowtie2 ≥2.5, samtools ≥1.17, numpy, pandas)

## Key File Paths (on TSCC)

```
/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl
/tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/ecliprepmap-0.1.0/bin/python
/tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_PE_singleNode
/tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_SE_singleNode
```

## Git / Source Control

- Current branch: `refdata-generation`
- MUST commit and push before deleting any files (per cleanup step)
- Main branch: `master`

## Verification Test Data

- Reference CWL outputs: `test-provenance/tests/ecliprepmap-1.0.0/`
  - SE reference: `wf_ecliprepmap_se/wf_ecliprepmap_se/results/INV_B.IP.umi.r1.fqTrTr.sorted.fq.barcode1.{nopipes,withpipes}.tsv`
  - PE reference: `wf_ecliprepmap_pe/wf_ecliprepmap_pe/results/204_01_RBFOX2.{nopipes,withpipes}.tsv`
- Example data: `examples/example_data_for_repeat_mapping_hg38/EXAMPLE_{PE,SE}.*`
- Reference files: `examples/inputs/hg38/`
