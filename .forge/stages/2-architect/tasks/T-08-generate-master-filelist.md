# T-08: Implement `generate_master_filelist.py`

<!-- DEPENDENCIES: T-01, T-02 -->
<!-- BLOCKS: T-09, T-10, T-11 -->

## Goal

Implement the Python script that emits the 5-column TSV (and identical `.list`) used as `ARGV[3]` by the Perl parse scripts. Row order must follow the hg38 grouping: Gencode → RepeatMasker → simplerepeats → tRNA → miRNA → custom FASTA records.

## Files Touched

- CREATE: `bin/python/refdata_generation/generate_master_filelist.py`

## Implementation Notes

CLI signature:
```
--gtf PATH                  required, for transcript_id → gene_name lookup
--parsed-ucsc PATH          required
--repeatmasker PATH         required, .tsv.gz
--simplerepeats PATH        optional, .tsv.gz
--trna PATH                 optional, .tsv.gz
--gff3 PATH                 optional, miRBase gff3
--custom-fasta PATH         optional, repeatable (NR_046233.2.fasta)
--output PATH               required, destination .tsv (also writes .list with identical content)
```

Algorithm: see Architecture Plan §4.6.

Column schema (5-col tab-separated, NO header):
1. `sequence_id` — ENST..., repeat name (with `_dup{n}`), trna name, miRNA ID, NR_... header.
2. `gene_id` — ENSG for Gencode, same as sequence_id otherwise.
3. `short_name` — gene_name from GTF attr for Gencode; element name otherwise.
4. `family` — Gencode, repeat family from RepeatMasker col2 source, Simple_repeat, tRNA, miRNA, rRNA.
5. `genelist_label` — `genelists.{FAMILY}` (e.g. `genelists.Gencode`, `genelists.tRNA`, `genelists.miRNA`, `genelists.rRNA`, `genelists.Simple_repeat`, `genelists.{REPEAT_FAMILY}`).

Row order grouping:
1. Gencode transcripts (sorted within group by gene_id, transcript_id).
2. RepeatMasker repeats (deduped by name; `_dup{n}` suffix for collisions).
3. Simple repeats (using `transcript_id` field for col1).
4. tRNA (if provided).
5. miRNA (if provided).
6. Custom FASTA entries (in argument order).

After writing `.tsv`, write identical content to a sibling file with `.list` extension (or symlink). The Perl script reads the `.list`-extension file per OG-04.

## Acceptance Criteria

- AC-T08-1: Output is 5-column tab-separated TSV with no header line.
- AC-T08-2: Every row has exactly 5 tab-separated fields.
- AC-T08-3: Row order matches the grouping defined above (Gencode → RM → Simple → tRNA → miRNA → custom).
- AC-T08-4: Sibling `.list` file with identical content (or symlink) is created.
- AC-T08-5: When `--trna` is omitted, no tRNA rows appear (no error).
- AC-T08-6: When `--gff3` is omitted, no miRNA rows appear (no error).
- AC-T08-7: NR_046233.2 records (one per FASTA header) appear in the output with family=`rRNA`, label=`genelists.rRNA`.
- AC-T08-8: Each `sequence_id` (col1) appears exactly once (deduplicated).

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

python bin/python/refdata_generation/generate_master_filelist.py \
  --gtf examples/inputs/mm10/downloaded/gencode.VM23.annotation.gtf.gz \
  --parsed-ucsc examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat \
  --repeatmasker examples/inputs/mm10/downloaded/mm10.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/mm10/downloaded/mm10.simplerepeats.tsv.gz \
  --trna examples/inputs/mm10/downloaded/mm10.trna.tsv.gz \
  --gff3 examples/inputs/mm10/downloaded/mmu.gff3 \
  --custom-fasta examples/inputs/mm10/downloaded/NR_046233.2.fasta \
  --output /tmp/test_mm10_filelist.tsv

awk -F'\t' '{print NF}' /tmp/test_mm10_filelist.tsv | sort -u
# Expect: only "5"

ls /tmp/test_mm10_filelist.tsv /tmp/test_mm10_filelist.list
# Expect: both present
```
