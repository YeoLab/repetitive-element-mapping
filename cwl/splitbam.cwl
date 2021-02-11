#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:1.0.0
    
baseCommand: [split_bam_to_subfiles_SEorPE.pl]

requirements:
  - class: InlineJavascriptRequirement

inputs:

  sam_file:
    type: File
    inputBinding:
      position: 1
  se_or_pe:
    type: string
    inputBinding:
      position: 2
    doc: "Either PE or SE (uppercase required)"
      
outputs:

  repsam_s:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tmp"
