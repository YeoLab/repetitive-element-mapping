#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:0.1.0
    
baseCommand: [merge_multiple_parsed_files.pl]

inputs:

  outputFile: 
    type: string
    inputBinding: 
      position: 1

  files:
    type: File[]
    inputBinding:
      position: 2

outputs:

  output:
    type: File
    outputBinding:
      glob: $(inputs.outputFile)


