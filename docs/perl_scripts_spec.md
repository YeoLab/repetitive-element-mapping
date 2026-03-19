# Perl Script Specification (Current Repository)

## Conversion Scope Decision
Per user-approved scope:
- Convert scripts used by CWL workflow execution path.
- Exclude legacy/non-CWL wrappers if not in workflow path:
  - `bin/perl/RepElement_pipeline_1dataset.pl` (excluded)
  - `bin/perl/duplicate_removal.pl` (excluded)

## In-scope scripts

### `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl`
Purpose:
- Align SE reads with bowtie2 and emit SAM-like repetitive-element assignments with family/element annotations.

CLI (from script):
1. `fastq_file1`
2. `bowtie_db` prefix
3. `output` SAM-like file
4. `filelist_file` mapping ENST -> type/priority

Outputs:
- Main: `<output>`
- Side files:
  - `<output>.bowtieout`
  - `<output>.multimapping_deleted`
  - `<output>.done`

Behavior highlights:
- Streams bowtie2 output (`--sensitive -a --reorder`).
- Applies family prioritization and special rRNA handling (`RNA45S` vs `RNA18S/28S/5-8S` rules).
- Emits `ZZ:Z:` tag with multiple compatible ENSTs.

### `bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl`
Purpose:
- Same as SE parser but for paired-end reads.

CLI:
1. `fastq_file1`
2. `fastq_file2`
3. `bowtie_db` prefix
4. `output` SAM-like file
5. `filelist_file`

Outputs:
- Main: `<output>`
- Side files: `.bowtieout`, `.multimapping_deleted`, `.done`

Behavior highlights:
- Reads bowtie2 SAM stream in paired-record chunks.
- Enforces read-name pairing checks.
- Maintains multi-family assignment bookkeeping.

### `bin/perl/split_bam_to_subfiles_SEorPE.pl`
Purpose:
- Split a SAM/BAM into 25 temporary files by first two UMI nucleotides.

CLI:
1. `sam_fi` (`.sam` or `.bam`)
2. `se_or_pe_flag` (`SE` or `PE`)

Outputs:
- `AA.<sam_basename>.tmp` ... `NN.<sam_basename>.tmp` in current working dir.

Behavior highlights:
- For BAM input, uses `samtools view -h` pipe.
- Skips unmapped reads (`flag 4` for SE; `77/141` for PE).
- In PE mode, assumes paired records are adjacent and may reorder interpretation by flag.

### `bin/perl/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
Purpose:
- Merge repeat-mapped SAM-like records with genome-mapped records and remove PCR duplicates.
- Produce rmdup/pre-rmdup SAM-like outputs and parsed count summaries.

CLI:
1. `repfamily_sam`
2. `gabe_rmRep_sam`
3. `eCLIP_read_type_flag` (`SE`/`PE`)
4. `gencode_gtf_file`
5. `gencode_tablebrowser_file`
6. `repmask_bed_fi`
7. `filelist_file`

Outputs (derived from input basename):
- `<repfamily_sam_short>.combined_w_uniquemap.rmDup.sam`
- `<repfamily_sam_short>.combined_w_uniquemap.prermDup.sam`
- `<rmDup>.parsed_v2.20201210.txt`
- `<output>.done`

Behavior highlights:
- SE and PE logic diverges after shared annotation loading.
- Uses per-read keys including UMI, strand, coordinates for deduping.
- Produces #READINFO lines + TOTAL/ELEMENT summaries.

### `bin/perl/merge_multiple_parsed_files.simplified_20191022.pl`
Purpose:
- Merge multiple per-prefix parsed files into one dataset-level parsed file.

CLI:
1. `output_fi`
2. `parsed_file_1`
3. `parsed_file_2`
4. ...

Outputs:
- `output_fi`
- Appends to `failed_jobs_list.txt` in inferred working directory using a lock file.

Behavior highlights:
- Sums readinfo counters across parsed inputs.
- Aggregates `TOTAL` and `ELEMENT` rows, then sorts by descending counts.
- Emits four canonical `#READINFO` lines in merged output.

## Out-of-scope scripts in this phase

### `bin/perl/RepElement_pipeline_1dataset.pl`
- Legacy wrapper generating shell/PBS pipeline commands.
- Not called by CWL workflow files.

### `bin/perl/duplicate_removal.pl`
- Older variant of duplicate-removal logic.
- Not referenced by `cwl/deduplicate.cwl`.
