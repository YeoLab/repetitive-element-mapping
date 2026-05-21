# Technical Constraints

## Execution Environment

- **Platform:** TSCC (Triton Shared Computing Cluster), SLURM scheduler
- **OS:** Linux (RHEL/CentOS-compatible)
- **Default shell:** bash

## Software Constraints

| Tool | Version / Source | Constraint |
|------|-----------------|------------|
| Python | 3.11 (dropin.yaml) | Use Python 3 only; no Python 2 compatibility needed |
| Perl | 5.10 (TSCC system perl via ecliprepmap module) | Perl ≤5.16 for deterministic hash iteration; do not require Perl 5.18+ |
| pybedtools | compatible with Python 3.11 | Must be importable in dropin conda env |
| bedtools | ≥2.30 | Underlying C tool for pybedtools |
| Bowtie2 | ≥2.5 | For bowtie2-build index generation |
| samtools | ≥1.17 | For FASTA indexing (samtools faidx) |
| pandas | in dropin.yaml | For TSV parsing if needed |
| numpy | in dropin.yaml | Statistical ops if needed |
| conda | Miniforge/Miniconda | Package manager |

## Module Load

```bash
module load ecliprepmap/1.0.0
```
This sets PATH for the correct Perl, Python, and bioinformatics tool versions on TSCC.

## Write Scope (Hard Constraint)

Scripts may ONLY write to:
- `examples/inputs/mm10/`
- `examples/inputs/mm39/`

Scripts must NOT modify:
- `examples/inputs/hg38/` (read-only ground truth)
- `bin/perl/*.pl` (except hardcoded values that block compatibility)
- Any CWL, Snakemake, or other workflow files

## Performance / Resource Constraints

- Reference generation for a single assembly is expected to complete within 2 hours on a standard TSCC compute node
- bowtie2-build may require up to 16 GB RAM for a large index; request appropriately if submitting via SLURM
- pybedtools getfasta may be slow on large FASTAs; indexing the FASTA first with `samtools faidx` is required
- No GPU resources required

## Conda Environment

All new Python scripts must be runnable within the existing dropin conda environment:
```
workflow/envs/dropin.yaml
```
No new conda packages may be added to dropin.yaml unless absolutely necessary. If a new package is required, it must be listed as a constraint for the architect to evaluate.

## File Format Constraints

- Output files must use tab (`\t`) as delimiter — no spaces
- Line endings must be Unix (`\n`) — no CRLF
- No trailing whitespace in BED or TSV output lines
- Coordinate system: all output BED files use 0-based half-open; all output parsed_ucsc_tableformat files use 0-based half-open (matching hg38 reference)

## CWL / Snakemake Compatibility Constraint

New reference files must work as drop-in replacements for the hg38 references in:
- `cwl/wf_ecliprepmap_se.cwl` / `cwl/wf_ecliprepmap_pe.cwl` (via YAML config)
- `workflow/rules/dropin_repelement.smk` (via `--config` CLI args)

No changes to CWL or Snakemake files are permitted.

## Naming Conventions

Output filenames must follow the established pattern:
- `gencode.{version}.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat`
- `bowtie2_index/` (directory)
- `UniqueGenomicElements.{assembly}.bed`
- `MASTER_FILELIST.{date}.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv`

The `{date}` token in MASTER_FILELIST must use the format `YYYYMMDD` (e.g., `20260514`).

## Similarity Threshold

- hg38 reproduction tests: ≥99% line similarity is the acceptance threshold
- 100% identical reproduction is the target for parsed_ucsc_tableformat (deterministic transformation)
- The 1% tolerance accounts for potential minor differences in repeat element extraction methodology
