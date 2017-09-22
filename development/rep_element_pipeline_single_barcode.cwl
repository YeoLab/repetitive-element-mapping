#!/usr/bin/env cwltool

cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement

inputs:
  # input R1.fastq.gz
  r1FastqGz:
    type: File
  # input R2.fastq.gz
  r2FastqGz:
    type: File
  # bowtie_reference_tar
  bowtieReferenceTar:
    type: File
  # perl script for parsing bowtie results inline
  parseBowtie2OutputRealtimeIncludefamilyPerlScript:
    type: File
  # fileListFile
  fileListFile1:
    type: File
  # fileListFile2
  fileListFile2:
    type: File
  # BAM file after removing repetitive elements (after 2nd STAR mapping)
  rmRepBAM:
    type: File
  gencodeGTF:
    type: File
  gencodeTableBrowser:
    type: File
  repMaskBEDFile:
    type: File
  outputSAMFile:
    type: string
  finalOutputSAMFile:
    type: string
  prefixes: string[]

outputs:
  # r1Fastq:
  #   type: File
  #   outputSource: gunzip_r1FastqGz/output
  # r2Fastq:
  #   type: File
  #   outputSource: gunzip_r2FastqGz/output
  map_repetitive_elements_output:
    type: File
    outputSource: map_repetitive_elements/output
  rep_split:
    type: File[]
    outputSource: split_rep_bam/outputs
  # rmrep_split:
  #   type: File[]
  #   outputSource: split_rmrep_bam/outputs
  rmDuped_sam:
    type: File[]
    outputSource: scattered_duplicate_removal_on_pair/combinedRmDupSam
  finalRmDupRepElementSam:
    type: File
    outputSource: combine_sam/output

steps:
  gunzip_r1_fastq_gz:
    run: gunzip.cwl
    in:
      input: r1FastqGz
    out:
      - output
  gunzip_r2_fastq_gz:
    run: gunzip.cwl
    in:
      input: r2FastqGz
    out:
      - output
  map_repetitive_elements:
    run: parse_bowtie2_output_realtime_includemultifamily_wrapper.cwl
    in:
      runScript: parseBowtie2OutputRealtimeIncludefamilyPerlScript
      read1: gunzip_r1_fastq_gz/output
      read2: gunzip_r2_fastq_gz/output
      indexTar: bowtieReferenceTar
      outputFile: outputSAMFile
      fileListFile: fileListFile1
    out:
      - output
  split_rep_bam:
    run: split_bam_to_subfiles.cwl
    in:
      inputFile: map_repetitive_elements/output
    out:
      - outputs
  split_rmrep_bam:
    run: split_bam_to_subfiles.cwl
    in:
      inputFile: rmRepBAM
    out:
      - outputs

  pair_matching_prefixes:
    run: get_pairs.cwl
    in:
      repInput: split_rep_bam/outputs
      rmRepInput: split_rmrep_bam/outputs
      prefix: prefixes
    scatter: prefix
    out:
      - rep
      - rmRep

  scattered_duplicate_removal_on_pair:
    run: duplicate_removal_inline_paired_count_region_other_reads.cwl
    in:
      repFamilySam: pair_matching_prefixes/rep
      rmRepSam: pair_matching_prefixes/rmRep
      gencodeGTF: gencodeGTF
      gencodeTableBrowser: gencodeTableBrowser
      repMaskBedFile: repMaskBEDFile
      fileList1: fileListFile1
      fileList2: fileListFile2
    scatter: [repFamilySam, rmRepSam]
    scatterMethod: dotproduct
    out:
      - combinedRmDupSam
      - combinedPreRmDupSam
      - parsedFile
      - doneFile

  combine_sam:
    run: combine_sam.cwl
    in:
      inputFiles: scattered_duplicate_removal_on_pair/combinedRmDupSam
      outputFile: finalOutputSAMFile
    out:
      - output

# input rmRep.bam (after second STAR mapping)
# input gencodeGTF
# input gencodeTableBrowser
# input repMaskBedFile
# input fileList1
# input fileList2

# unzip adaptertrim/polyAtrim read1.fastq.gz
# unzip adaptertrim/polyAtrim read2.fastq.gz

# split_bam_to_subfiles.pl rmRep.bam
# parse_bowtie2_output_realtime_includemultifamily.pl read1.fastq read2.fastq filelist.UpdatedSimpleRepeat working_dir/ samfile.sam
# split_bam_to_subfiles.pl samfile.sam

# for each AA/AC/AT/AT/CA/etc. perform duplicate removal and remove the tmp files.
# duplicate_removal_inline_paired.count_region_other_reads.pl samfile.sam.AA.tmp pre-merged-X1A-round2.rmRep.bam.AA.tmp

# cat XYZ.*.AA.tmp.combined_w_uniquemap.rmDup.sam >> XYZ.combined_w_uniquemap.rmDup.sam
# cat XYZ.*.AA.tmp.combined_w_uniquemap.prermDup.sam >> XYZ.combined_w_uniquemap.prermDup.sam
# gzip XYZ.combined_w_uniquemap.prermDup.sam

# merge_multiple_parsed_files.pl *.sam.parsed (both X1A and X1B)
