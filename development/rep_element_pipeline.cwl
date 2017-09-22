#!/usr/bin/env

cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:

  # input R1.fastq.gz
  barcode1R1FastqGz:
    type: File
  # input R2.fastq.gz
  barcode1R2FastqGz:
    type: File
  # BAM file after removing repetitive elements (after 2nd STAR mapping)
  barcode1RmRepBAM:
    type: File
  barcode1OutputSAMFile:
    type: string
  barcode1FinalOutputSAMFile:
    type: string

  # input R1.fastq.gz
  barcode2R1FastqGz:
    type: File
  # input R2.fastq.gz
  barcode2R2FastqGz:
    type: File
  barcode2RmRepBAM:
    type: File
  barcode2OutputSAMFile:
    type: string
  barcode2FinalOutputSAMFile:
    type: string

  bowtieReferenceTar:
    type: File
  # perl script for parsing bowtie results inline
  parseBowtie2OutputRealtimeIncludefamilyPerlScript:
    type: File
  gencodeGTF:
    type: File
  gencodeTableBrowser:
    type: File
  repMaskBEDFile:
    type: File

  # fileListFile
  fileListFile1:
    type: File

  # fileListFile2
  fileListFile2:
    type: File

  prefixes: string[]

  finalOutput:
    type: string

outputs:
  barcode1CombinedOutput:
    type: File
    outputSource: barcode_1/finalRmDupRepElementSam
  barcode2CombinedOutput:
    type: File
    outputSource: barcode_2/finalRmDupRepElementSam
  combinedOutput:
    type: File
    outputSource: combine_barcodes/output

steps:
  barcode_1:
    run: rep_element_pipeline_single_barcode.cwl
    in:
      r1FastqGz: barcode1R1FastqGz
      r2FastqGz: barcode1R2FastqGz
      bowtieReferenceTar: bowtieReferenceTar
      parseBowtie2OutputRealtimeIncludefamilyPerlScript: parseBowtie2OutputRealtimeIncludefamilyPerlScript
      fileListFile1: fileListFile1
      fileListFile2: fileListFile2
      rmRepBAM: barcode1RmRepBAM
      gencodeGTF: gencodeGTF
      gencodeTableBrowser: gencodeTableBrowser
      repMaskBEDFile: repMaskBEDFile
      outputSAMFile: barcode1OutputSAMFile
      finalOutputSAMFile: barcode1FinalOutputSAMFile
      prefixes: prefixes
    out:
      - map_repetitive_elements_output
      - rep_split
      - rmDuped_sam
      - finalRmDupRepElementSam

  barcode_2:
    run: rep_element_pipeline_single_barcode.cwl
    in:
      r1FastqGz: barcode2R1FastqGz
      r2FastqGz: barcode2R2FastqGz
      bowtieReferenceTar: bowtieReferenceTar
      parseBowtie2OutputRealtimeIncludefamilyPerlScript: parseBowtie2OutputRealtimeIncludefamilyPerlScript
      fileListFile1: fileListFile1
      fileListFile2: fileListFile2
      rmRepBAM: barcode2RmRepBAM
      gencodeGTF: gencodeGTF
      gencodeTableBrowser: gencodeTableBrowser
      repMaskBEDFile: repMaskBEDFile
      outputSAMFile: barcode2OutputSAMFile
      finalOutputSAMFile: barcode2FinalOutputSAMFile
      prefixes: prefixes
    out:
      - map_repetitive_elements_output
      - rep_split
      - rmDuped_sam
      - finalRmDupRepElementSam

  combine_barcodes:
    run: combine_sam.cwl
    in:
      inputFiles: [barcode_1/finalRmDupRepElementSam, barcode_2/finalRmDupRepElementSam]
      outputFile: finalOutput
    out:
      - output

