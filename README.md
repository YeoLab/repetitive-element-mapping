# repetitive-element-mapping

Python + Snakemake pipeline for repetitive element mapping with deterministic mini fixtures for validation.

## What is in this branch
- Python implementation scripts in `bin/python/`.
- Snakemake workflow entrypoint in `Snakefile`.
- Rule modules in `workflow/rules/`.
- Validation scripts in `workflow/scripts/`.
- Mini fixtures and expected outputs in `tests/fixtures/mini/`.

## Repository layout
- `Snakefile`: top-level workflow entrypoint.
- `config/config.yaml`: runtime file paths for mini validation inputs and expected outputs.
- `bin/python/split_bam_to_subfiles_SEorPE.py`: split SAM/BAM by UMI prefix.
- `bin/python/merge_multiple_parsed_files.simplified_20191022.py`: merge parsed statistics files.
- `workflow/rules/se_foundation.smk`: currently implemented, testable workflow stages.
- `workflow/scripts/verify_split_manifest.py`: checksum validation for split outputs.
- `workflow/scripts/verify_merge_against_expected.py`: checksum validation for merged output.

## Requirements
- Python 3.11+
- Snakemake 9+
- samtools 1.23+
- bowtie2 2.5+
- pytest 9+

## Installation

### 1) Clone
```bash
git clone https://github.com/YeoLab/repetitive-element-mapping.git
cd repetitive-element-mapping
git checkout codex/python-conversion
```

### 2) Create environment
```bash
mamba create -y -p ./.conda-env -c conda-forge -c bioconda \
  python=3.11 snakemake bowtie2 samtools pandas numpy pyyaml pytest
```

### 3) Activate
```bash
mamba activate ./.conda-env
```

## Input specification

## CWL-style YAML Input Alignment
This branch supports two equivalent ways to provide run configuration:

1) Populate top-level keys in `config/config.yaml`.
2) Pass an existing CWL-style job YAML at runtime:
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p all \
  --config cwl_input_yaml=/absolute/path/to/job.yaml
```

Supported CWL-compatible keys:
- `barcode1r1FastqGz`, `barcode1rmRepBam`
- `barcode1Inputr1FastqGz`, `barcode1InputrmRepBam`
- `bowtie2_db`, `bowtie2_prefix`, `fileListFile1`, `gencodeGTF`, `gencodeTableBrowser`, `repMaskBEDFile`
- `dataset`, `prefixes`, `se_or_pe`

For `class: File` / `class: Directory` objects, the `path` value is extracted. Relative `path` values are resolved relative to the YAML file location.

The separate `mini_validation` block is only for deterministic fixture tests while conversion is in progress.
See `docs/config_mapping.md` for the one-to-one mapping details.
The current implemented Snakemake stages use mini fixture inputs configured in `config/config.yaml`:

- `mini_validation.source_sam_gz`: compressed SAM-like file used for split stage.
- `mini_validation.merge_input_1`: parsed stats input file #1 for merge stage.
- `mini_validation.merge_input_2`: parsed stats input file #2 for merge stage.
- `mini_validation.split_manifest`: expected checksums for 25 split output files.
- `mini_validation.merge_expected`: expected merged parsed output.

You can update these paths to point to your own test data as long as formats match.

## Output specification
Running the current workflow stages produces:

- `results/mini/split/python_tmp/*.tmp`: 25 prefix-split files (`AA`..`NN`).
- `results/mini/split/verified.ok`: split-stage validation success marker.
- `results/mini/merge/merged.python.parsed`: merged parsed output.
- `results/mini/merge/verified.ok`: merge-stage validation success marker.

## Run instructions (deployment)

### Local execution
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p all
```

`HOME=$(pwd)` keeps Snakemake cache files inside workspace and avoids host cache permission issues.

### Re-run all stages from scratch
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all
```

### Cluster deployment pattern
Use your site profile/launcher as usual, for example:
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake --profile <your-profile> all
```

## Testing

### Unit/integration tests
```bash
./.conda-env/bin/python -m pytest -q tests/test_mini_fixtures_manifest.py
```

### Workflow validation run
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all
```

## Regenerating expected artifacts
If you intentionally update core Python logic, regenerate mini expected artifacts and commit them:
```bash
python3 scripts/generate_python_expected_artifacts.py
```

This refreshes:
- `tests/fixtures/mini/expected/split.expected.manifest.tsv`
- `tests/fixtures/mini/expected/merged.expected.parsed`
