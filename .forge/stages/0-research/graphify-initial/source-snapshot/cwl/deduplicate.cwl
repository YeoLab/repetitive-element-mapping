#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: ResourceRequirement
    ramMin: 32000
        
hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:1.0.0

baseCommand: [duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl]

inputs:

  repFamilySam:
    type: File
    inputBinding:
      position: 1
    label: "input"
    doc: "input sam file"

  rmRepSam:
    type: File
    inputBinding:
      position: 2
    label: "rmRep.sam alignment file"
    doc: "rmRep sam alignment file"
  
  se_or_pe:
    type: string
    inputBinding:
      position: 3
    doc: "Either PE or SE (uppercase required)"
    
  gencodeGTF:
    type: File
    inputBinding:
      position: 4
    label: "gencode GTF File"
    doc: "gencode GTF File"

  gencodeTableBrowser:
    type: File
    inputBinding:
      position: 5
    label: "gencode GTF file in UCSC table browser format"
    doc: "gencode GTF file in UCSC table browser format"

  repMaskBedFile:
    type: File
    inputBinding:
      position: 6
    label: "repeatmasker bedfile"
    doc: "bedfile of repeat regions"

  fileList1:
    type: File
    inputBinding:
      position: 7
    doc: "tsv with 5 fields: ENST/ENSG/name/chrom/genelist.file"

outputs:

  deduplicatedRmDupSam:
    type: File
    outputBinding:
      glob: "*combined_w_uniquemap.rmDup.sam"
    label: "combined unique mapping + rep element mapping sam"

  deduplicatedPreRmDupSam:
    type: File
    outputBinding:
      glob: "*combined_w_uniquemap.prermDup.sam"
    label: "combined unique mapping + rep element mapping sam before rmdup"

  parsedFile:
    type: File
    outputBinding:
      glob: "*.parsed_v2.20201210.txt"
    label: "combined unique mapping + rep element mapping sam"

  doneFile:
    type: File
    outputBinding:
      glob: "*.done"
    label: "done file"