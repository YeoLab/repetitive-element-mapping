#!/bin/bash

cwltool --cachedir splitbamcache \
split_bam_to_subfiles.cwl \
split_bam_to_subfiles_repsam.yaml
