# T-10: Generate mm10 Reference Files (AC-02, AC-05, AC-08, AC-11)

<!-- DEPENDENCIES: T-03, T-05, T-07, T-09 -->
<!-- BLOCKS: T-12 -->

## Goal

Generate all four reference files for mm10 using the validated scripts, placing outputs in `examples/inputs/mm10/`.

## Files Touched

- CREATE: `examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat`
- CREATE: `examples/inputs/mm10/bowtie2_index/` (directory with 7 files)
- CREATE: `examples/inputs/mm10/UniqueGenomicElements.mm10.bed`
- CREATE: `examples/inputs/mm10/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv`
- CREATE: `examples/inputs/mm10/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list`

## Implementation Notes

Date token: use `20260514` (today's date in YYYYMMDD form per tech-constraints §"Naming Conventions"). Document the date in NOTES so it can be reproduced exactly.

Execution order:
1. `generate_parsed_ucsc_tableformat.py` → writes parsed_ucsc to `examples/inputs/mm10/`.
2. `generate_bowtie2_index.py` → reads parsed_ucsc, writes index to `examples/inputs/mm10/bowtie2_index/`.
3. `generate_unique_genomic_elements.py` → reads parsed_ucsc, writes BED.
4. `generate_master_filelist.py` → reads parsed_ucsc, writes TSV+list.

## Acceptance Criteria

- AC-T10-1: All four output files exist under `examples/inputs/mm10/` with the expected names. (AC-02, AC-05, AC-08, AC-11)
- AC-T10-2: `parsed_ucsc_tableformat` has the correct 11-column header and ≥1 data row.
- AC-T10-3: `bowtie2_index/` directory contains the 6 `.bt2` files plus the source `.fa` (7 total). `bowtie2-inspect --summary` exits 0.
- AC-T10-4: `NR_046233.2` sequences are present in the FASTA file (`grep -c "^>NR_046233.2" <fa> >= 1`).
- AC-T10-5: `UniqueGenomicElements.mm10.bed` is 6-column tab-separated BED, ≥1 row.
- AC-T10-6: `MASTER_FILELIST.*.tsv` is 5-column tab-separated TSV, ≥1 row, AND `.list` sibling exists with identical content.
- AC-T10-7: MASTER_FILELIST contains at least one `NR_046233.2*` row, at least one tRNA row, and at least one miRNA row.
- AC-T10-8: Scripts refuse to overwrite `examples/inputs/hg38/` (assert_writable test from T-01 still passes).
- AC-T10-9: All scripts exit 0 with no traceback in stderr.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
bash bin/python/refdata_generation/run_assembly.sh --assembly mm10 --date 20260514
# OR invoke the four scripts manually (see Architecture Plan §4.7 for arguments)

ls examples/inputs/mm10/
ls examples/inputs/mm10/bowtie2_index/
bowtie2-inspect --summary examples/inputs/mm10/bowtie2_index/MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat | head
grep -c "^>NR_046233.2" examples/inputs/mm10/bowtie2_index/MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa
```
