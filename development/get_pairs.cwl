#!/usr/bin/env cwltool

cwlVersion: v1.0

class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

inputs:
  prefix:
    type: string
  repInput:
    type: File[]
  rmRepInput:
    type: File[]

outputs:
  rep:
    type: File
  rmRep:
    type: File

expression: |
   ${
      var prefix = inputs.prefix;
      var rep = inputs.repInput;
      var rmRep = inputs.rmRepInput;
      var returnRep = '';
      var returnRmRep = '';
      for (var i = 0; i < rep.length; i++) {
        if (rep[i].basename.indexOf(prefix) == 0) {
          returnRep = rep[i];
        }
        if (rmRep[i].basename.indexOf(prefix) == 0) {
          returnRmRep = rmRep[i];
        }
      }
      return {'rep': returnRep, 'rmRep': returnRmRep}
    }