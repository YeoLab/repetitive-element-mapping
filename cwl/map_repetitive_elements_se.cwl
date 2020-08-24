#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

hints: 
  - class: DockerRequirement
    dockerPull: brianyee/repetitive_element_mapping:0.1.0

# # wrapped perl script for parsing bowtie results inline
baseCommand: [parse_bowtie2_output_realtime_includemultifamily_SE.pl]


inputs:


  #read1 trimmed fastq
  read1:
    type: File
    inputBinding:
      position: 1

  bowtie2_db:
    type: Directory

  bowtie2_prefix:
    type: string

  output_file:
    default: ""
    type: string
    inputBinding:
      position: 3
      valueFrom: |
        ${
          if (inputs.output_file == "") {
            return inputs.read1.nameroot + ".Rep.sam";
          }
          else {
            return inputs.output_file;
          }
        }
  file_list_file:
    type: File
    inputBinding:
      position: 4

arguments:
  - valueFrom: $(inputs.bowtie2_db.path + "//" + inputs.bowtie2_prefix)
    position: 2
    shellQuote: false

outputs:

  rep_sam:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_file == "") {
            return inputs.read1.nameroot + ".Rep.sam";
          }
          else {
            return inputs.output_file;
          }
        }
