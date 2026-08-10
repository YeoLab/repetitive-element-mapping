# Stage 5-QA Summary: Perl-to-Python Translation Verification

## Status: COMPLETE

## What Was Verified

**SE small pipeline**: Ran `examples/repeat_mapping_SE_small.yaml` end-to-end. Output `seCLIP_small.nopipes.tsv` has 62 rows with correct header. Top element: RNA28S (703 reads). SE deduplication was verified byte-identical to Perl 5.10.1 on 5 prefix bins (AA, AC, TT, NN, GG): identical READINFO counts, TOTAL counts, ELEMENT assignments, and rmDup.sam content.

**PE small pipeline**: Ran `examples/repeat_mapping_PE_small.yaml` end-to-end. Output `peCLIP_small.nopipes.tsv` has 76 rows with correct header. Top elements: unique_distintron (569), unique_proxintron (383), RNA28S (346). No `RG:Z:foo` contamination.

## Bugs Found and Fixed

1. **Stale PE results** (empty rep.sam): Old run used bowtie2-not-found path. Fixed by deleting stale `results/pe_small/` and re-running.

2. **RG:Z:foo in PE nopipes** (commit d44ca22): `split_bam_to_subfiles.py` was swapping the r1/r2 strings for flags 147/163/403/419 in the written output. The Perl reference only swaps internal arrays for UMI extraction, not the written lines. The swap caused pairing mismatches in dedup, leading to the "else" branch pulling `RG:Z:foo` (last SAM optional field) as an element name. Fix: removed the string swap.

## New Python Scripts

Five Perl scripts translated: `split_bam_to_subfiles.py`, `merge_parsed_files.py`, `map_repetitive_elements_se.py`, `map_repetitive_elements_pe.py`, `deduplicate.py`. All Snakemake rules updated to use Python.

## Remaining

Full dataset validation (Step 5 of translate_perl.md) requires SLURM job submission.
