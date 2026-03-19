# Snakemake SE Foundation

## Implemented stages
1. Split SAM by UMI prefix (`split_bam_to_subfiles_SEorPE.py`).
2. Merge parsed statistics files (`merge_multiple_parsed_files.simplified_20191022.py`).

## Files
- `Snakefile`
- `config/config.yaml`
- `workflow/rules/se_foundation.smk`
- `workflow/scripts/verify_split_manifest.py`
- `workflow/scripts/verify_merge_against_expected.py`

## Validation strategy
- Split stage output files are validated against
  `tests/fixtures/mini/expected/split.expected.manifest.tsv`.
- Merge stage output is validated against
  `tests/fixtures/mini/expected/merged.expected.parsed`.

## Run
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all
```

## Next extensions
- Add mapping stage rules.
- Add deduplication stage rules.
- Add full end-to-end stage chaining for production data.
