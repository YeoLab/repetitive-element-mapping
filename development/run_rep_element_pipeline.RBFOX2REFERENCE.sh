#!/bin/bash

cwltool \
--outdir /projects/ps-yeolab3/bay001/rep_element_reference/ \
--cachedir /home/bay001/scratch/codebase/rep_element_pipeline/RBFOX2_CACHE \
rep_element_pipeline.cwl \
rep_element_pipeline.RBFOX2REFERENCE.yaml
