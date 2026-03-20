# repetitive-element-mapping

Python + Snakemake implementation of the repetitive element mapping workflow.

This branch contains no CWL execution path and no Perl runtime dependency.

## Overview

There are two execution profiles:

1. `mini` profile (default): deterministic fixture validation for core converted components.
2. `dropin` profile: full SE/PE Snakemake workflow driven by CWL-style YAML job inputs.

## Installation

```bash
git clone https://github.com/YeoLab/repetitive-element-mapping.git
cd repetitive-element-mapping
git checkout codex/python-conversion-python-only-cleanup

mamba create -y -p ./.conda-env -c conda-forge -c bioconda \
  python=3.11 snakemake bowtie2 samtools pandas numpy pyyaml pytest
```

## Inputs

The ETL parser accepts either:

1. Legacy CWL-style YAML (`class` + `path` objects).
2. New simplified YAML schema (`samples` + `references` blocks).

Both normalize to the same internal config keys.

Examples:
- Legacy: `examples/repeat_mapping_SE.yaml`, `examples/repeat_mapping_PE.yaml`
- Simplified: `examples/repeat_mapping_SE.simple.yaml`, `examples/repeat_mapping_PE.simple.yaml`

Use either with:

```bash
--config cwl_input_yaml=/absolute/path/to/job.yaml
```

Normalized keys include:
- `dataset`
- `barcode1r1FastqGz`, `barcode1r2FastqGz`, `barcode1rmRepBam`
- `barcode2r1FastqGz`, `barcode2r2FastqGz`, `barcode2rmRepBam`
- `barcode1Inputr1FastqGz`, `barcode1Inputr2FastqGz`, `barcode1InputrmRepBam`
- `bowtie2_db`, `bowtie2_prefix`, `fileListFile1`
- `gencodeGTF`, `gencodeTableBrowser`, `repMaskBEDFile`
- `prefixes`, `se_or_pe`

For both formats, relative paths are resolved relative to the YAML file location.

## Outputs

### Drop-in workflow (SE)
Generated in `results/dropin/` and exported by wrapper to `<job_name>/results`:
- `<ip>.preRmDup.sam.gz`
- `<input>.preRmDup.sam.gz`
- `<ip>.rmDup.sam.gz`
- `<input>.rmDup.sam.gz`
- `<ip>.parsed`
- `<input>.parsed`
- `<ip>.nopipes.tsv`
- `<ip>.withpipes.tsv`

Wrapper also writes:
- `REPELEMENTMAPPING_<job>_OUTPUT.json`
- `REPELEMENTMAPPING_<job>_VERSION-1.0.0`

### Drop-in workflow (PE)
Generated similarly, including combined IP outputs from barcode1+barcode2 and input outputs.

## Running

### 1) Mini validation profile

```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all
```

### 2) Drop-in profile directly with Snakemake

SE:
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all \
  --config pipeline_profile=dropin run_mode=SE cwl_input_yaml=/abs/job_se.yaml
```

PE:
```bash
HOME=$(pwd) ./.conda-env/bin/snakemake -j1 -p -F all \
  --config pipeline_profile=dropin run_mode=PE cwl_input_yaml=/abs/job_pe.yaml
```

### 3) Drop-in wrapper commands (legacy entrypoint style)

SE:
```bash
./wf/eCLIP_repelement_SE /abs/job_se.yaml
```

PE:
```bash
./wf/eCLIP_repelement_PE /abs/job_pe.yaml
```

## Testing

Unit tests:
```bash
./.conda-env/bin/python -m pytest -q \
  tests/test_mini_fixtures_manifest.py \
  tests/test_split_merge_python_expected.py \
  tests/test_parse_se_python.py \
  tests/test_parse_pe_python.py \
  tests/test_cwl_yaml_adapter.py
```

Drop-in synthetic fixtures:
- SE fixture YAML: `tests/fixtures/dropin/se_job.yaml`
- PE fixture YAML: `tests/fixtures/dropin/pe_job.yaml`

Example:
```bash
./wf/eCLIP_repelement_SE tests/fixtures/dropin/se_job.yaml
./wf/eCLIP_repelement_PE tests/fixtures/dropin/pe_job.yaml
```
