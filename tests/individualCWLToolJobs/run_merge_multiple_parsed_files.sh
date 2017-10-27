#!/usr/bin/env bash

WD=/projects/ps-yeolab3/bay001/rep_element_reference

merge_multiple_parsed_files.pl \
RBFOX2.barcode1.final.parsed.new \
${WD}/AA.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/AC.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/AG.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/AN.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/AT.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/CA.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/CC.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/CG.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/CN.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/CT.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/GA.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/GC.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/GG.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/GN.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/GT.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/NA.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/NC.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/NG.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/NN.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/NT.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/TA.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/TC.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/TG.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/TN.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed \
${WD}/TT.tmp.barcode1.repmapped.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed

