#!/usr/bin/env python

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [/home/bay001/projects/codebase/repetitive-element-mapping/development/parse_bowtie2_output_realtime_includemultifamily_wrapper.sh]

inputs:

  runScript:
    type: File
    inputBinding:
      position: 1
    label: "runscript"
    doc: "runscript"

  read1:
    type: File
    inputBinding:
      position: 2
    label: "read1 trimmed fastq"
    doc: "read1 trimmed fastq"

  read2:
    type: File
    inputBinding:
      position: 3
    label: "read2 trimmed fastq"
    doc: "read2 trimmed fastq"

  indexDir:
    type: File
    inputBinding:
      position: 4
    label: "bowtie2 index"
    doc: "bowtie2 index"

  outputFile:
    type: string
    inputBinding:
      position: 5
    label: "output sam file"
    doc: "output sam file"

  fileListFile:
    type: File
    inputBinding:
      position: 6
    label: "file list file"
    doc: "dunno ask eric"

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.outputFile)
    label: "output"
    doc: "File containing output"


