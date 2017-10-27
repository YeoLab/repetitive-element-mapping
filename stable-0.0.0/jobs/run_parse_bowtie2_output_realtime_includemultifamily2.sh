#!/usr/bin/env bash

perl ../parse_bowtie2_output_realtime_includemultifamily.pl \
/projects/ps-yeolab3/bay001/rep_element_reference/RBFOX2-204-CLIP_S1_R1.A01_204_01_RBFOX2.adapterTrim.round2.fastq \
/projects/ps-yeolab3/bay001/rep_element_reference/RBFOX2-204-CLIP_S1_R2.A01_204_01_RBFOX2.adapterTrim.round2.fastq \
/home/bay001/projects/codebase/repetitive-element-mapping/data/bowtie_reference/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
/home/bay001/projects/codebase/repetitive-element-mapping/current/jobs \
original_outfile_run2.sam
