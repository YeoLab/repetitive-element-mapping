# 03 - User Flows

## Flow A: SE Pipeline (Single Barcode IP + Input)

```
Input:
  barcode1r1FastqGz  +  barcode1rmRepBam
  barcode1Inputr1FastqGz  +  barcode1InputrmRepBam
  reference files

Step 1 [map_repetitive_elements_SE]:
  parse_bowtie2_output_realtime_includemultifamily_SE.pl
    read1 → bowtie2 → streaming parse → <dataset>.barcode1.Rep.sam

Step 2 [splitbam - rep]:
  split_bam_to_subfiles_SEorPE.pl <dataset>.barcode1.Rep.sam SE
    → AA.tmp, AC.tmp, AG.tmp, ..., NN.tmp  (25 files)

Step 3 [splitbam - rmrep]:
  split_bam_to_subfiles_SEorPE.pl barcode1rmRepBam SE
    → AA.tmp, AC.tmp, AG.tmp, ..., NN.tmp  (25 files)

Step 4 [deduplicate x25 - scattered by prefix]:
  For each prefix in [AA, AC, AG, AT, AN, CA, CC, CG, CT, CN, GA, GC, GG, GT, GN,
                       TA, TC, TG, TT, TN, NA, NC, NG, NT, NN]:
    duplicate_removal_inline_paired...pl
      <prefix>.rep.tmp  <prefix>.rmrep.tmp  SE  gencodeGTF  gencodeTableBrowser
      repMaskBEDFile  fileListFile1
    →  <prefix>.*combined_w_uniquemap.rmDup.sam
       <prefix>.*combined_w_uniquemap.prermDup.sam
       <prefix>.*.parsed_v2.20201210.txt
       <prefix>.*.done

Step 5 [concatenate + gzip]:
  cat all 25 .rmDup.sam → <dataset>.barcode1.rmDup.sam → gzip → .rmDup.sam.gz
  cat all 25 .prermDup.sam → <dataset>.barcode1.preRmDup.sam → gzip → .preRmDup.sam.gz

Step 6 [combine_parsed]:
  merge_multiple_parsed_files.simplified_20191022.pl
    <dataset>.barcode1.parsed  [all 25 .parsed_v2 files]
  → <dataset>.barcode1.parsed

Repeat Steps 1-6 for Input (barcode1Inputr1FastqGz + barcode1InputrmRepBam)
  → <dataset>.input.parsed

Step 7 [calculate_fold_change]:
  calculate_fold_change_from_parsed_files.py
    --ip_parsed <dataset>.barcode1.parsed
    --input_parsed <dataset>.input.parsed
    --out_file_nopipes <dataset>.nopipes.tsv
    --out_file_withpipes <dataset>.withpipes.tsv
```

## Flow B: PE Pipeline (Two IP Barcodes + Input)

Identical to Flow A for each barcode (barcode1, barcode2, input) — all three run in parallel.

After barcode1 and barcode2 complete:

```
Step 6b [combine_parsed - merge barcodes]:
  merge_multiple_parsed_files.simplified_20191022.pl
    <dataset>.combined.parsed  [all parsed files from barcode1 AND barcode2]
  → <dataset>.combined.parsed

Step 7 [calculate_fold_change]:
  calculate_fold_change_from_parsed_files.py
    --ip_parsed <dataset>.combined.parsed
    --input_parsed <dataset>.input.parsed
    --out_file_nopipes <dataset>.nopipes.tsv
    --out_file_withpipes <dataset>.withpipes.tsv
```

Note: In the CWL PE workflow, `step_concatenate_rmDup` concatenates rmDup files from BOTH barcode1 and barcode2 before gzipping. The final `.rmDup.sam.gz` IP output is a merge of both barcodes.

## Flow C: Downsampled Dataset Generation (pre-pipeline)

```
1. Subsample PE FASTQ files to ~10,000 reads covering all chromosomes
2. Extract BAM reads whose read names exist in the subsampled FASTQ
3. Ensure ≥100 reads per UMI 2-nt prefix (all 25 prefixes covered)
4. Write to examples/inputs/downsampled/
5. Generate repeat_mapping_PE_small.yaml and repeat_mapping_SE_small.yaml
```
