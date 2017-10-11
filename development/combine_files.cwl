#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [cat]

inputs:

  inputFiles:
    type: File[]
    inputBinding:
      position: 1
    label: "input"
    doc: "input file"
  outputFile:
    type: string

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.outputFile)

stdout: $(inputs.outputFile)
