# Implementation Report: mm10/mm39 Reference Data Generation
<!-- FORGE_STAGE: 3-implement -->
<!-- STATUS: IN_PROGRESS -->
<!-- STARTED_UTC: 2026-05-14T00:00:00Z -->
<!-- UPDATED_UTC: 2026-05-14T00:00:00Z -->

## Task Overview
| # | Task | Status | Verify Cycles | Last Updated |
|---|------|--------|---------------|-------------|
| T-01 | Shared utilities (_shared.py) | PENDING | 0 | - |
| T-02 | generate_parsed_ucsc_tableformat.py | PENDING | 0 | - |
| T-03 | hg38 parsed reproduction (AC-01) | PENDING | 0 | - |
| T-04 | generate_bowtie2_index.py | PENDING | 0 | - |
| T-05 | hg38 bowtie2 reproduction (AC-04) | PENDING | 0 | - |
| T-06 | generate_unique_genomic_elements.py | PENDING | 0 | - |
| T-07 | hg38 unique elements reproduction (AC-07) | PENDING | 0 | - |
| T-08 | generate_master_filelist.py | PENDING | 0 | - |
| T-09 | hg38 master filelist reproduction (AC-10) | PENDING | 0 | - |
| T-10 | Generate mm10 references | PENDING | 0 | - |
| T-11 | Generate mm39 references | PENDING | 0 | - |
| T-12 | Pipeline integration & wrapper | PENDING | 0 | - |

## Key Design Notes (Pre-Implementation Discoveries)

- The MASTER_FILELIST uses gene_name (not "Gencode") as family for Gencode transcripts
- Repeat element names are uppercased in MASTER_FILELIST (e.g., AluJb -> ALUJB)
- Repeat element family (e.g., "Alu") is derived from the uppercase gene_id prefix
- The col5 (genelist_label) for repeats is the RepBase class (SINE, LINE, DNA) NOT "genelists.{family}"
- The col5 for Gencode transcripts IS "genelists.{gene_name}"
- The hg38 tRNA entries in MASTER_FILELIST use a DIFFERENT naming format than what's in hg38.trna.tsv.gz
  - MASTER_FILELIST: tRNA-Ala-AGC-1-1 (from a separate tRNA FASTA file)
  - hg38.trna.tsv.gz: nm-tRNA-Tyr-GTA-chr1-142 format
  - mm10.trna.tsv.gz: chr1.tRNA1555-GluTTC format
- miRNA entries: col1=MI_ID, col2=col3=col4="miRNA", col5=name from gff3
- Simple repeat entries: col1=transcript_id (trf, trf_dup1, ...), col2=col3=col4="Simple_repeat", col5="Simple_repeat"
- rRNA/NR_ entries (hg38): come from a separate custom FASTA with specific naming like NR_046235.3-18S
- For mm10/mm39: NR_046233.2.fasta has a single header (not pre-split); appended verbatim per ADR-04

## Files Modified
| File | Action | Task | Notes |
|------|--------|------|-------|

## Task Reports
(Populated as tasks complete)
