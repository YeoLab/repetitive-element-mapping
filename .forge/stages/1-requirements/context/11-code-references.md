# 11 - Code References

## Key Repository Files

### CWL Workflows (source of truth for translation)

| File | Purpose |
|------|---------|
| `cwl/wf_ecliprepmap_pe.cwl` | Top-level PE workflow (2 barcodes + input) |
| `cwl/wf_ecliprepmap_se.cwl` | Top-level SE workflow (1 barcode + input) |
| `cwl/wf_ecliprepmap_pe_1barcode.cwl` | Per-barcode PE sub-workflow with scatter |
| `cwl/wf_ecliprepmap_se_1barcode.cwl` | Per-barcode SE sub-workflow with scatter |
| `cwl/map_repetitive_elements_pe.cwl` | Bowtie2 PE mapping step (8 CPU, 16 GB) |
| `cwl/map_repetitive_elements_se.cwl` | Bowtie2 SE mapping step (8 CPU, 16 GB) |
| `cwl/splitbam.cwl` | Split BAM/SAM by UMI prefix |
| `cwl/getpair.cwl` | ExpressionTool to match rep/rmrep .tmp pairs by prefix |
| `cwl/deduplicate.cwl` | Deduplication (32 GB) |
| `cwl/concatenate.cwl` | cat command wrapper |
| `cwl/gzip.cwl` | gzip -c wrapper |
| `cwl/combine.cwl` | merge parsed files (uses InitialWorkDirRequirement) |
| `cwl/calculate_fold_change_from_parsed_files.cwl` | Fold change calculation |

### Perl Scripts (bin/perl/)

| Script | Role |
|--------|------|
| `parse_bowtie2_output_realtime_includemultifamily_PE.pl` | PE bowtie2 + parse |
| `parse_bowtie2_output_realtime_includemultifamily_SE.pl` | SE bowtie2 + parse |
| `split_bam_to_subfiles_SEorPE.pl` | UMI prefix splitting |
| `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl` | Deduplication |
| `duplicate_removal.pl` | Softlink to above |
| `merge_multiple_parsed_files.simplified_20191022.pl` | Merge parsed stats |

### Python Scripts (bin/python/)

| Script | Role |
|--------|------|
| `calculate_fold_change_from_parsed_files.py` | Compute fold enrichment from parsed files |

### Example Config YAMLs

| File | Purpose |
|------|---------|
| `examples/repeat_mapping_PE.yaml` | Full PE config (3 barcodes) |
| `examples/repeat_mapping_SE.yaml` | Full SE config (2 samples) |

### Example Data

- `examples/example_data_for_repeat_mapping_hg38/EXAMPLE_PE.*` — PE FASTQs and BAMs (3 samples)
- `examples/example_data_for_repeat_mapping_hg38/EXAMPLE_SE.*` — SE FASTQs and BAMs (2 samples)

### Reference Data (hg38)

- `examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.*` — bowtie2 index
- `examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` — fileListFile1
- `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` — gencodeTableBrowser
- `examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf` — gencodeGTF
- `examples/inputs/hg38/UniqueGenomicElements.hg38.bed` — repMaskBEDFile

### Validation Reference Outputs

- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_se/wf_ecliprepmap_se/results/INV_B.IP.umi.r1.fqTrTr.sorted.fq.barcode1.nopipes.tsv`
- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_se/wf_ecliprepmap_se/results/INV_B.IP.umi.r1.fqTrTr.sorted.fq.barcode1.withpipes.tsv`
- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_pe/wf_ecliprepmap_pe/results/204_01_RBFOX2.nopipes.tsv`
- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_pe/wf_ecliprepmap_pe/results/204_01_RBFOX2.withpipes.tsv`

## Critical CWL Logic to Replicate in Snakemake

### getpair.cwl (ExpressionTool)
The CWL getpair tool matches rep and rmrep .tmp files by checking if the filename starts with the prefix:
```javascript
if (rep_s[i].basename.indexOf(prefix) == 0) { prefixrep = rep_s[i]; }
if (rmrep_s[i].basename.indexOf(prefix) == 0) { prefixrmrep = rmrep_s[i]; }
```
In Snakemake: use a `{prefix}` wildcard; the splitbam rule must name outputs as `{prefix}.rep.tmp` and `{prefix}.rmrep.tmp` (or similar), and the deduplicate rule uses `{prefix}` as a wildcard to match both.

### dataset name derivation in CWL
In CWL, the dataset name for a barcode is derived from the r1 FASTQ filename:
```javascript
return self.nameroot + ".barcode1";  // self = r1FastqGz
```
`nameroot` strips one extension (e.g., `.gz` → `EXAMPLE_PE.rep1_clip.A01.r1.fqTrTr.sorted.fq`).

In Snakemake: use the `dataset` config key directly as the base name (simpler and more explicit).

### PE final rmDup concatenation
In CWL `wf_ecliprepmap_pe.cwl`, the `step_concatenate_rmDup` concatenates rmDup files from BOTH barcode1 AND barcode2 (merge_flattened linkMerge). The final IP SAM output covers both barcodes.

### combine_parsed for PE
In CWL PE, `step_combine_parsed` receives parsed files from barcode1 AND barcode2 (merge_flattened), then the merged output is used for fold change against the input's combined parsed file.
