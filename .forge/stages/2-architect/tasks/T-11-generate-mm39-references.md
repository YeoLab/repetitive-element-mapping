# T-11: Generate mm39 Reference Files (AC-03, AC-06, AC-09, AC-12)

<!-- DEPENDENCIES: T-03, T-05, T-07, T-09 -->
<!-- BLOCKS: T-12 -->

## Goal

Generate all four reference files for mm39 using the validated scripts, with `--trna` and `--gff3` omitted. WARNINGs must be logged for both absent inputs.

## Files Touched

- CREATE: `examples/inputs/mm39/gencode.VM38.annotation.gtf.parsed_ucsc_tableformat`
- CREATE: `examples/inputs/mm39/bowtie2_index/` (directory with 7 files)
- CREATE: `examples/inputs/mm39/UniqueGenomicElements.mm39.bed`
- CREATE: `examples/inputs/mm39/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv`
- CREATE: `examples/inputs/mm39/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list`

## Implementation Notes

Execution mirrors T-10 but without `--trna` and `--gff3` flags. The wrapper script `run_assembly.sh --assembly mm39` automatically detects the absence of those source files and omits the flags.

Expected WARNINGs:
- generate_bowtie2_index: `WARNING: --trna not provided; tRNA entries omitted` and `WARNING: --gff3 not provided; miRNA entries omitted`.
- generate_unique_genomic_elements: same two WARNINGs.
- generate_master_filelist: same two WARNINGs.

## Acceptance Criteria

- AC-T11-1: All four output files exist under `examples/inputs/mm39/`. (AC-03, AC-06, AC-09, AC-12)
- AC-T11-2: `bowtie2-inspect --summary` on the mm39 index returns exit 0.
- AC-T11-3: `UniqueGenomicElements.mm39.bed` contains NO tRNA-named rows and NO miRNA-named rows (the categories are entirely absent).
- AC-T11-4: `MASTER_FILELIST.*.tsv` contains NO `tRNA` family rows and NO `miRNA` family rows; DOES contain NR_046233.2 rRNA rows.
- AC-T11-5: At least two WARNINGs are logged to stderr by the wrapper (one for trna, one for gff3).
- AC-T11-6: `parsed_ucsc_tableformat` is a valid 11-column file with header.
- AC-T11-7: All scripts exit 0 with no traceback.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
bash bin/python/refdata_generation/run_assembly.sh --assembly mm39 --date 20260514 2>&1 | tee /tmp/mm39_run.log

grep -c WARNING /tmp/mm39_run.log
# Expect: >= 2

ls examples/inputs/mm39/
ls examples/inputs/mm39/bowtie2_index/
awk -F'\t' '$4 == "tRNA" || $4 == "miRNA"' examples/inputs/mm39/MASTER_FILELIST.20260514.*.tsv | wc -l
# Expect: 0
```
