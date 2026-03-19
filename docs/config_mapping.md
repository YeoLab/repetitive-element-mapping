# Config Mapping: CWL YAML to Snakemake config.yaml

## Why they looked different
The previous `config/config.yaml` only covered mini validation fixtures for staged conversion tests.
It did not yet expose the full CWL-style runtime inputs.

## Current mapping approach
`config/config.yaml` includes CWL-compatible key names at top level.
You can also pass an existing CWL-style job YAML via:

```bash
snakemake --config cwl_input_yaml=/path/to/job.yaml ...
```

When `cwl_input_yaml` is provided, mapped keys from that file override defaults in `config/config.yaml`.

CWL YAML key -> Snakemake key:
- `dataset` -> `dataset`
- `barcode1r1FastqGz` -> `barcode1r1FastqGz`
- `barcode1rmRepBam` -> `barcode1rmRepBam`
- `barcode1Inputr1FastqGz` -> `barcode1Inputr1FastqGz`
- `barcode1InputrmRepBam` -> `barcode1InputrmRepBam`
- `bowtie2_db` -> `bowtie2_db`
- `bowtie2_prefix` -> `bowtie2_prefix`
- `fileListFile1` -> `fileListFile1`
- `gencodeGTF` -> `gencodeGTF`
- `gencodeTableBrowser` -> `gencodeTableBrowser`
- `repMaskBEDFile` -> `repMaskBEDFile`
- `prefixes` -> `prefixes`
- `se_or_pe` -> `se_or_pe`

## Validation-only section
The `mini_validation` subsection is separate and only used for deterministic fixture testing of currently implemented stages.

## `class/path` object handling
For any supported key with a value like:

```yaml
someKey:
  class: File
  path: relative/or/absolute/path
```

the adapter uses `path` as the runtime value.
Relative paths are resolved relative to the job YAML file directory.

## Practical note
For production runs, populate the top-level CWL-compatible keys with real file paths.
For conversion regression tests, keep using `mini_validation`.
