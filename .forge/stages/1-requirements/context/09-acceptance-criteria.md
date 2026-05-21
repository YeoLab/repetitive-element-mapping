# Acceptance Criteria

## AC-01: hg38 parsed_ucsc_tableformat Reproduction

**Given** `gencode.v33.chr_patch_hapl_scaff.annotation.gtf` in `examples/inputs/hg38/downloaded/`  
**When** `generate_parsed_ucsc_tableformat.py` is run against that GTF  
**Then** the output is line-for-line identical to `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` (249,044 lines including header)

---

## AC-02: mm10 parsed_ucsc_tableformat Generated

**Given** `gencode.VM23.annotation.gtf.gz` in `examples/inputs/mm10/downloaded/`  
**When** `generate_parsed_ucsc_tableformat.py` is run against that GTF  
**Then** a valid parsed_ucsc_tableformat file is written to `examples/inputs/mm10/` with:
- Header line `#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds`
- One row per unique transcript_id from the VM23 GTF
- `exonCount` = number of exon features for each transcript
- `exonStarts` and `exonEnds` with trailing commas, 0-based coordinates

---

## AC-03: mm39 parsed_ucsc_tableformat Generated

**Given** `gencode.VM38.annotation.gtf.gz` in `examples/inputs/mm39/downloaded/`  
**When** `generate_parsed_ucsc_tableformat.py` is run  
**Then** same structural requirements as AC-02 but for VM38 annotation

---

## AC-04: hg38 bowtie2 FASTA Reproduced Within 99% Similarity

**Given** hg38 source files in `examples/inputs/hg38/downloaded/` and `GRCh38_no_alt_analysis_set.fasta`  
**When** `generate_bowtie2_index.py` is run for hg38  
**Then** the generated combined FASTA contains ≥99% of the sequences present in `examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa` (measured by sequence header count)

---

## AC-05: mm10 bowtie2_index Generated and Functional

**Given** mm10 source files (repeatmasker, simplerepeats, trna, mmu.gff3, NR_046233.2.fasta) and mm10 genome FASTA  
**When** `generate_bowtie2_index.py` is run for mm10  
**Then**:
- `examples/inputs/mm10/bowtie2_index/` contains 7 files (6 `.bt2` index files + 1 `.fa` source FASTA)
- `bowtie2-inspect --summary <prefix>` exits 0 confirming a valid index
- NR_046233.2 sequences appear in the FASTA with correct rRNA naming convention

---

## AC-06: mm39 bowtie2_index Generated and Functional

Same as AC-05 but for mm39, and without tRNA/GFF3 inputs (those are absent for mm39).

---

## AC-07: hg38 UniqueGenomicElements Reproduced Within 99% Similarity

**Given** hg38 source files  
**When** `generate_unique_genomic_elements.py` is run for hg38  
**Then** the generated `UniqueGenomicElements.hg38.bed` matches the existing reference within ±1% line count (5,618,483 expected) and column format is identical (6-column BED)

---

## AC-08: mm10 UniqueGenomicElements Generated

**Given** mm10 repeatmasker, simplerepeats, trna, mmu.gff3, and mm10 parsed_ucsc_tableformat  
**When** `generate_unique_genomic_elements.py` is run for mm10  
**Then**:
- `UniqueGenomicElements.mm10.bed` is written to `examples/inputs/mm10/`
- File is 6-column BED with no header
- Each element has two corresponding `-proximal` entries (500 bp upstream and downstream)
- No truncation at chromosome boundaries (or if truncated, `max(0, start-500)` is used)

---

## AC-09: mm39 UniqueGenomicElements Generated (Without tRNA/miRNA)

**Given** mm39 repeatmasker, simplerepeats, and mm39 parsed_ucsc_tableformat (no trna, no gff3)  
**When** `generate_unique_genomic_elements.py` is run without `--trna` and `--gff3`  
**Then**:
- Script logs two WARNINGs: one for missing trna, one for missing gff3
- `UniqueGenomicElements.mm39.bed` is generated successfully
- File contains repeatmasker + simplerepeats + parsed_ucsc_tableformat entries with proximal flanks
- File does NOT contain any tRNA or miRNA entries

---

## AC-10: hg38 MASTER_FILELIST Reproduced Within 99% Similarity

**Given** hg38 source files  
**When** `generate_master_filelist.py` is run for hg38  
**Then** generated file matches existing `MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` within ±1% line count (26,422 expected) and is 5 columns tab-separated with no header

---

## AC-11: mm10 MASTER_FILELIST Generated

**Given** mm10 parsed_ucsc_tableformat, repeatmasker, simplerepeats, trna, mmu.gff3, NR_046233.2.fasta  
**When** `generate_master_filelist.py` is run for mm10  
**Then**:
- TSV written to `examples/inputs/mm10/`
- 5 columns, tab-separated, no header
- NR_046233.2 entries present with correct family labels
- tRNA entries present (mm10 has trna file)
- miRNA entries present (mm10 has mmu.gff3)

---

## AC-12: mm39 MASTER_FILELIST Generated (Without tRNA/miRNA)

Same structural requirements as AC-11 but without tRNA and miRNA categories (those sources are absent for mm39). NR_046233.2 entries must still be present.

---

## AC-13: Perl Script Compatibility (No Modification Required)

**Given** mm10 or mm39 MASTER_FILELIST passed as ARGV[3] to `parse_bowtie2_output_realtime_includemultifamily_SE.pl`  
**When** a small test eCLIP BAM is mapped using Bowtie2 and piped through the Perl script  
**Then** the script runs without error and produces a valid `.parsed_v2` output file

---

## AC-14: Pipeline Dry-Run Success

**Given** all four mm10 reference files generated (parsed_ucsc_tableformat, bowtie2_index, UniqueGenomicElements.mm10.bed, MASTER_FILELIST.*.tsv)  
**When** `snakemake -s workflow/rules/dropin_repelement.smk --config ... -n` is run with mm10 reference paths  
**Then** the dry-run completes with exit code 0, listing all expected rules

---

## AC-15: Missing ID Warning Fires Correctly

**Given** a test scenario where >1% of expected FASTA IDs cannot be retrieved  
**When** `generate_bowtie2_index.py` is run  
**Then** a WARNING is printed to stderr and `missing_ids_report.txt` is written to the output directory, containing the list of missing IDs and their expected sources

---

## AC-16: Scripts Accept Both .gtf and .gtf.gz Inputs

**Given** either an uncompressed `.gtf` or a gzipped `.gtf.gz` file  
**When** any GTF-parsing script is run  
**Then** the script handles both formats without error (auto-detects from file extension)
