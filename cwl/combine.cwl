#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: $(inputs.files)
      
hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:1.0.0
    
baseCommand: [merge_multiple_parsed_files.simplified_20191022.pl]

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


