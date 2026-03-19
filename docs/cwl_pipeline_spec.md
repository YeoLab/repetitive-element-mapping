# CWL Pipeline Specification

## Scope
This document describes the current CWL graph with emphasis on SE (`wf_ecliprepmap_se.cwl`) because conversion priority is SE first, followed by PE.

## Top-level SE workflow: `cwl/wf_ecliprepmap_se.cwl`

### Inputs
- `dataset` (string)
- IP sample:
  - `barcode1r1FastqGz` (File)
  - `barcode1rmRepBam` (File)
- Input/control sample:
  - `barcode1Inputr1FastqGz` (File)
  - `barcode1InputrmRepBam` (File)
- Reference/config:
  - `bowtie2_db` (Directory)
  - `bowtie2_prefix` (string)
  - `fileListFile1` (File)
  - `gencodeGTF` (File)
  - `gencodeTableBrowser` (File)
  - `repMaskBEDFile` (File)
- Runtime defaults:
  - `prefixes` (string[], 25 two-base UMI prefixes from `AA` to `NN`)
  - `se_or_pe` (string, default `SE`)

### Steps
1. `step_ecliprepmap_barcode1`
- Runs subworkflow `wf_ecliprepmap_se_1barcode.cwl` for IP sample.
- Derives internal `dataset` name from input FASTQ nameroot with `.barcode1` suffix.

2. `step_ecliprepmap_input`
- Runs same subworkflow for Input sample.
- Derives internal `dataset` name from input FASTQ nameroot with `.input` suffix.

3. `step_calculate_fold_change_from_parsed_files`
- Runs `calculate_fold_change_from_parsed_files.cwl` over combined parsed outputs from IP vs Input.

### Outputs
- Pre-rmdup SAM-like gz:
  - `output_ip_concatenated_pre_rmDup_sam_file`
  - `output_input_concatenated_pre_rmDup_sam_file`
- Rmdup SAM-like gz:
  - `output_barcode1_concatenated_rmDup_sam_file`
  - `output_input_concatenated_rmDup_sam_file`
- Parsed stats and fold-change tables:
  - `output_ip_parsed`
  - `output_input_parsed`
  - `output_nopipes`
  - `output_withpipes`

## One-barcode SE subworkflow: `cwl/wf_ecliprepmap_se_1barcode.cwl`

### Step graph (execution order)
1. `step_map_repetitive_elements` -> `map_repetitive_elements_se.cwl`
2. `step_splitbam_repsam` -> `splitbam.cwl`
3. `step_splitbam_rmrepbam` -> `splitbam.cwl`
4. `step_getpair` (scatter over 25 `prefixes`) -> `getpair.cwl`
5. `step_deduplicate` (dotproduct scatter over prefix-matched rep/rmrep files) -> `deduplicate.cwl`
6. `step_concatenate_rmDup` -> `concatenate.cwl`
7. `step_concatenate_preRmDup` -> `concatenate.cwl`
8. `step_gzip_rmDup` -> `gzip.cwl`
9. `step_gzip_preRmDup` -> `gzip.cwl`
10. `step_combine_parsed` -> `combine.cwl`

### Subworkflow outputs
- `output_repeat_mapped_sam_file` (`.Rep.sam`)
- `output_rmDup_sam_files` (array of per-prefix rmdup SAM-like files)
- `output_pre_rmDup_sam_files` (array of per-prefix pre-rmdup SAM-like files)
- `output_concatenated_rmDup_sam_file` (gz)
- `output_concatenated_preRmDup_sam_file` (gz)
- `output_parsed_files` (array of parsed per-prefix files)
- `output_combined_parsed_file` (combined `.parsed`)

## CommandLineTool contracts used by SE

### `cwl/map_repetitive_elements_se.cwl`
- Base command: `parse_bowtie2_output_realtime_includemultifamily_SE.pl`
- Positional args:
  1. `read1`
  2. `bowtie2_db.path + "//" + bowtie2_prefix`
  3. `output_file` (default: `<read1.nameroot>.Rep.sam`)
  4. `file_list_file`
- Output: `rep_sam` (SAM-like file)

### `cwl/splitbam.cwl`
- Base command: `split_bam_to_subfiles_SEorPE.pl`
- Positional args: `sam_file`, `se_or_pe`
- Output glob: `*.tmp`
- Behavior depends on SE/PE flag and input extension (`.sam` vs `.bam`)

### `cwl/getpair.cwl` (ExpressionTool)
- Inputs: one `prefix`, `rep_s[]`, `rmrep_s[]`
- Returns the file in each array whose basename starts with prefix.
- Outputs: `prefixrep`, `prefixrmrep`

### `cwl/deduplicate.cwl`
- Base command: `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
- Positional args:
  1. `repFamilySam`
  2. `rmRepSam`
  3. `se_or_pe`
  4. `gencodeGTF`
  5. `gencodeTableBrowser`
  6. `repMaskBedFile`
  7. `fileList1`
- Outputs:
  - `*combined_w_uniquemap.rmDup.sam`
  - `*combined_w_uniquemap.prermDup.sam`
  - `*.parsed_v2.20201210.txt`
  - `*.done`

### `cwl/concatenate.cwl`
- Base command: `cat`
- Inputs: `files[]`
- `stdout` redirected to provided `concatenated_output` filename.

### `cwl/gzip.cwl`
- Base command: `gzip -c <input>`
- Output filename: `<input.basename>.gz`

### `cwl/combine.cwl`
- Base command: `merge_multiple_parsed_files.simplified_20191022.pl`
- Args: `outputFile`, then `files[]`
- Output glob: `outputFile`

### `cwl/calculate_fold_change_from_parsed_files.cwl`
- Base command: `calculate_fold_change_from_parsed_files.py`
- Required args: `--ip_parsed`, `--input_parsed`
- Defaults:
  - `--out_file_nopipes`: `<ip_parsed.nameroot>.nopipes.tsv`
  - `--out_file_withpipes`: `<ip_parsed.nameroot>.withpipes.tsv`

## PE workflow notes (for subsequent conversion phase)
- `cwl/wf_ecliprepmap_pe.cwl` runs `wf_ecliprepmap_pe_1barcode.cwl` for barcode1, barcode2, and input.
- It adds cross-barcode `combine`, `concatenate`, and gzip at top level.
- The one-barcode PE subworkflow mirrors SE structure but uses PE mapping parser (`map_repetitive_elements_pe.cwl`) and PE reads (R1/R2).
