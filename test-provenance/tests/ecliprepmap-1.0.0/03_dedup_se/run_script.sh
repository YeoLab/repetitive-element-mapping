#!/bin/bash

duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl \
inputs/INV_B.IP.umi.r1.fqTrTr.sorted.fq.Rep.sam.AA.tmp \
inputs/INV_B.IP.umi.r1.fq.genome-mappedSoSo.bam.AA.tmp \
SE \
/projects/ps-yeolab4/genomes/hg38/gencode/v33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf \
/projects/ps-yeolab4/genomes/hg38/gencode/v33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat \
/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/hg38/UniqueGenomicElements.hg38.bed \
/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list
