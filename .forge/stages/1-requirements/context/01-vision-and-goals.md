# 01 - Vision and Goals

## Project Purpose

The eCLIP repetitive element mapping pipeline (ecliprepmap) identifies which repetitive genomic elements (transposons, repeats, rRNA, etc.) are bound by RNA-binding proteins, using enhanced CLIP (eCLIP) sequencing data.

The existing pipeline is implemented in CWL (Common Workflow Language) and runs via `cwltool` or the `eCLIP_repelement_PE/SE` launchers on TSCC. This feature translates the CWL workflows to Snakemake so the pipeline can run natively on SLURM via the `snakemake9` conda environment, removing the CWL runtime dependency and enabling better integration with existing Snakemake-based eCLIP infrastructure.

## Goals

1. **Exact output fidelity**: Snakemake outputs must match CWL outputs for `.nopipes.tsv` and `.withpipes.tsv` files (the primary scientific results). Minor differences due to random tie-breaking in Perl hash iteration are acceptable.

2. **Single Snakefile with SE/PE dispatch**: One Snakefile handles both single-end and paired-end datasets, dispatching to appropriate rules based on `se_or_pe` config parameter.

3. **Testable with downsampled data**: Small downsampled datasets enable fast iteration and CI-style validation without requiring full dataset runtimes.

4. **SLURM-native execution**: All rules run as SLURM jobs via the existing `profiles/tscc2_snakemake9/` profile (partition gold, account csd792).

5. **Drop-in replacement**: After validation, unused CWL files and Perl scripts are removed; README is updated with Snakemake usage instructions.

## Definition of Done

- Snakemake workflow runs end-to-end on both downsampled and full SE/PE datasets
- Outputs match CWL reference outputs in test-provenance/
- README updated with clear usage instructions for both dataset sizes
- All changes committed and pushed to git
