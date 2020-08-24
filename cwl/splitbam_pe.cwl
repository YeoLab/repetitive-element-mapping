#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:0.1.0
    
baseCommand: [split_bam_to_subfiles_PE.pl]

requirements:
  - class: InlineJavascriptRequirement

inputs:

  sam_file:
    type: File
    inputBinding:
      position: 1

outputs:

  repsam_s:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tmp"
