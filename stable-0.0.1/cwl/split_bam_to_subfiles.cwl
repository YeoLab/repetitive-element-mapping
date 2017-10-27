#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [split_bam_to_subfiles.pl]

requirements:
  - class: InlineJavascriptRequirement

inputs:

  inputFile:
    type: File
    inputBinding:
      position: 1
    label: "input"
    doc: "input sam file"

outputs:
  outputs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tmp"
    label: "output"
    doc: "Files containing output"
