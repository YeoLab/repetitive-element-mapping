#!/usr/bin/env cwltool

cwlVersion: v1.0

class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

inputs:
  input:
    type: File[]

outputs:
  outputs:
    type: File[]

expression: |
   ${
      var output = inputs.input;
      var sorted = output.sort(function(a,b) { return a.location > b.location } )
      return {'outputs': sorted }
    }