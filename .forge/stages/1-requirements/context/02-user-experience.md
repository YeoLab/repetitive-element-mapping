# 02 - User Experience

## User Profile

Bioinformatics analysts at Yeo Lab (TSCC), comfortable with YAML config files and SLURM cluster job submission. Familiar with eCLIP data processing but not necessarily Snakemake internals.

## Primary Interface: Config YAML + Snakemake CLI

Users invoke the pipeline by:

1. Creating a config YAML file (modeled after existing `examples/repeat_mapping_PE.yaml` or `repeat_mapping_SE.yaml`)
2. Running Snakemake with the config and SLURM profile

Example invocation (SE):
```bash
module load singularitypro
conda activate snakemake9
snakemake -s Snakefile --configfile examples/repeat_mapping_SE_small.yaml --profile profiles/tscc2_snakemake9
```

Example invocation (PE):
```bash
snakemake -s Snakefile --configfile examples/repeat_mapping_PE_small.yaml --profile profiles/tscc2_snakemake9
```

## Config YAML Structure

The config must support both SE and PE modes. The existing CWL YAML uses flat keys:

```yaml
se_or_pe: SE   # or PE

dataset: seCLIP_example

barcode1r1FastqGz: /path/to/file.fq.gz
barcode1rmRepBam: /path/to/file.bam

barcode1Inputr1FastqGz: /path/to/input.fq.gz
barcode1InputrmRepBam: /path/to/input.bam

# PE-only additional fields:
barcode1r2FastqGz: /path/to/file.r2.fq.gz
barcode2r1FastqGz: /path/to/file.r1.fq.gz
barcode2r2FastqGz: /path/to/file.r2.fq.gz
barcode2rmRepBam: /path/to/file.bam
barcode1Inputr2FastqGz: /path/to/input.r2.fq.gz

# Reference files (same for SE and PE)
bowtie2_db: /path/to/bowtie2_index/
bowtie2_prefix: MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat
fileListFile1: /path/to/MASTER_FILELIST...tsv
gencodeGTF: /path/to/gencode.gtf
gencodeTableBrowser: /path/to/gencode.gtf.parsed_ucsc_tableformat
repMaskBEDFile: /path/to/UniqueGenomicElements.hg38.bed
```

## Config Validation

- If `se_or_pe: PE`, barcode2 fields (barcode2r1FastqGz, barcode2r2FastqGz, barcode2rmRepBam, barcode1Inputr2FastqGz) are required; validation error if missing.
- If `se_or_pe: SE`, barcode2 fields must not be present; validation error if found.
- Validation occurs at Snakemake startup (before any jobs are submitted).

## User-Visible Outputs

All outputs land in a configurable output directory (default: alongside the config YAML or in a `results/` subdirectory):

| File | Description |
|------|-------------|
| `<dataset>.nopipes.tsv` | Primary result: unambiguously mapped families with fold enrichment |
| `<dataset>.withpipes.tsv` | All families including multi-family reads |
| `<dataset>.barcode1.rmDup.sam.gz` | Deduplicated SAM (IP barcode1) |
| `<dataset>.barcode1.preRmDup.sam.gz` | Pre-dedup SAM (IP barcode1) |
| `<dataset>.barcode1.parsed` | Per-element read counts (IP barcode1) |
| `<dataset>.input.rmDup.sam.gz` | Deduplicated SAM (input) |
| `<dataset>.input.preRmDup.sam.gz` | Pre-dedup SAM (input) |
| `<dataset>.input.parsed` | Per-element read counts (input) |
| (PE only) `<dataset>.barcode2.*` | Same set for barcode2 |
| (PE only) `<dataset>.combined.parsed` | Merged barcode1+barcode2 parsed file |
