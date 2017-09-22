cwlVersion: v1.0

class: CommandLineTool

baseCommand: [perl, /home/bay001/projects/codebase/repetitive-element-mapping/development/split_bam_to_subfiles.pl]

requirements:
  - class: InlineJavascriptRequirement

inputs:

  inputFile:
    type: File
    inputBinding:
      position: 1
    label: "input"
    doc: "input sam file"

# outputs:
  # AA:
  #   type: File
  #   outputBinding:
  #     glob: "AA*.tmp"
  # AT:
  #   type: File
  #   outputBinding:
  #     glob: "AT*.tmp"
outputs:
  outputs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tmp"
    label: "output"
    doc: "Files containing output"
