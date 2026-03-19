# Config Mapping: CWL YAML to Snakemake config.yaml

## Why they looked different
The previous `config/config.yaml` only covered mini validation fixtures for staged conversion tests.
It did not yet expose the full CWL-style runtime inputs.

## Current mapping approach
`config/config.yaml` now includes CWL-compatible key names at top level.

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

## Practical note
For production runs, populate the top-level CWL-compatible keys with real file paths.
For conversion regression tests, keep using `mini_validation`.
