# 04 - Data Models

## Input Files

### FASTQ files
- SE: `<sample>.r1.fqTrTr.sorted.fq.gz` — Read1 after inline barcode trimming
- PE: `<sample>.r1.fqTrTr.sorted.fq.gz` and `<sample>.r2.fqTrTr.sorted.fq.gz`
- Read name format (SE): `K00180:223:...:10013:5587_CGCCTTGCCG` — UMI is after `_`, first 2 nt are the split prefix
- Read name format (PE): `AGAAA:SN1001:449:...:12732:17919` — UMI is before `:`, first 2 nt are the split prefix

### BAM files (rmRep)
- BAM file of reads that mapped to the unique genome (post-STAR 2nd pass, repetitive elements removed)
- Must contain the same read names as the corresponding FASTQ file

### Reference files

| File | Description |
|------|-------------|
| bowtie2 index (.bt2 files) | Repeat element database index |
| fileListFile1 (.tsv, 5 columns) | ENST/ENSG/name/chrom/genelist mapping |
| gencodeGTF (.gtf) | Gencode annotation |
| gencodeTableBrowser (.parsed_ucsc_tableformat) | GTF in UCSC table browser format |
| repMaskBEDFile (.bed) | BED file of unique genomic regions (UniqueGenomicElements.hg38.bed) |

## Intermediate Files

### Rep SAM-like file
- Output of bowtie2 parse step
- Named: `<dataset>.<barcode>.Rep.sam`
- SAM-like format with custom multi-family fields

### Split .tmp files
- Named by UMI prefix: `AA.tmp`, `AC.tmp`, ..., `NN.tmp` (25 files each)
- Two sets: rep .tmp and rmrep .tmp (from BAM)
- Must be in same directory for deduplicate script pairing

### Deduplicated SAM files (per prefix)
- `<prefix>.*combined_w_uniquemap.rmDup.sam` — post-dedup
- `<prefix>.*combined_w_uniquemap.prermDup.sam` — pre-dedup
- `<prefix>.*.parsed_v2.20201210.txt` — read count stats
- `<prefix>.*.done` — completion marker

### Parsed file
- Tab-separated, with `#READINFO` header lines
- Columns: element name | family | read count | RPM | ...
- `|` in family column = multi-family read (excluded from nopipes output)

## Output Files

### Primary outputs (scientific results)
- `<dataset>.nopipes.tsv` — fold enrichment for unambiguously mapped families
- `<dataset>.withpipes.tsv` — fold enrichment including multi-family reads

### SAM outputs
- `<dataset>.<barcode>.rmDup.sam.gz` — gzipped concatenated deduplicated SAM
- `<dataset>.<barcode>.preRmDup.sam.gz` — gzipped concatenated pre-dedup SAM

### Parsed statistics
- `<dataset>.<barcode>.parsed` — merged per-element counts from all 25 prefixes
- (PE only) `<dataset>.combined.parsed` — merged counts from barcode1 + barcode2

## File Naming Convention

Dataset name derives from the `dataset` config field.
Barcode suffix: `.barcode1`, `.barcode2`, `.input`
All outputs should be deterministic from the `dataset` name + config paths.
