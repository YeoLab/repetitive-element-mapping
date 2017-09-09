#!/usr/bin/env python

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [perl, /home/bay001/projects/codebase/repetitive-element-mapping/development/parse_bowtie2_output_realtime_includemultifamily.pl]

inputs:

  read1:
    type: File
    inputBinding:
      position: 1
    label: "read1 trimmed fastq"
    doc: "read1 trimmed fastq"

  read2:
    type: File
    inputBinding:
      position: 2
    label: "read2 trimmed fastq"
    doc: "read2 trimmed fastq"

  indexDir:
    type: Directory
    inputBinding:
      position: 3
    label: "bowtie2 index"
    doc: "bowtie2 index"

  tempDirectory:
    type: Directory
    inputBinding:
      position: 4
    label: "working directory"
    doc: "working directory"

  outputFile:
    type: string
    inputBinding:
      position: 5
    label: "output sam file"
    doc: "output sam file"

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.outputFile)
    label: "output"
    doc: "File containing output"


