# Snakemake SE Foundation (Step 5a)

## Scope completed in this step
Converted and validated the following SE workflow stages in Snakemake:
1. Split SAM by UMI prefix (`split_bam_to_subfiles_SEorPE`) using Python conversion.
2. Merge parsed per-prefix outputs (`merge_multiple_parsed_files`) using Python conversion.

## Snakemake files
- `Snakefile`
- `config/config.yaml`
- `workflow/rules/se_foundation.smk`
- `workflow/scripts/verify_split_manifest.py`
- `workflow/scripts/verify_merge_vs_perl.py`

## Validation approach
- Split stage: compare Python split outputs against a Perl-derived manifest
  (`tests/fixtures/mini/expected/split.perl.manifest.tsv`).
- Merge stage: run Perl merge and Python merge on same mini parsed inputs and
  compare structured outputs with tolerant float checks.

## Run command
```bash
HOME=/Users/brianyee/Documents/github/repetitive-element-mapping \
/Users/brianyee/Documents/github/repetitive-element-mapping/.conda-env/bin/snakemake -j1 -p -F all
```

Note: setting `HOME` to workspace avoids sandbox cache permission issues for Snakemake.

## Pending for full SE migration
Still to implement in Snakemake for full SE parity:
- repetitive-element mapping (`parse_bowtie2...SE`)
- split of rmRep BAM + prefix pairing
- deduplication stage wiring
- concatenate/gzip final outputs
- fold-change stage wiring and expected output comparisons
