# T-04: Implement `generate_bowtie2_index.py`

<!-- DEPENDENCIES: T-01, T-02 -->
<!-- BLOCKS: T-05, T-10, T-11 -->

## Goal

Implement the Python script that assembles a combined FASTA from Gencode transcripts, repeat elements, simple repeats, tRNA, miRNA, and custom FASTA inputs, then builds a Bowtie2 index. Must handle optional inputs gracefully and report missing IDs above 1% threshold.

## Files Touched

- CREATE: `bin/python/refdata_generation/generate_bowtie2_index.py`

## Implementation Notes

CLI signature:
```
--gtf PATH                  required, .gtf or .gtf.gz
--parsed-ucsc PATH          required, parsed_ucsc_tableformat from T-02
--repeatmasker PATH         required, .tsv.gz (9-col)
--simplerepeats PATH        required, .tsv.gz (9-col)
--fasta PATH                required, genome FASTA (must have .fai or will be auto-created)
--trna PATH                 optional, .tsv.gz
--gff3 PATH                 optional, miRBase gff3
--custom-fasta PATH         optional, repeatable
--output-dir PATH           required, directory for index files
--output-prefix STR         required, bowtie2 index prefix (e.g. MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat)
```

Algorithm: see Architecture Plan §4.4.

Sequence sources in this concatenation order:
1. Gencode transcript sequences (one entry per transcript, header `>{transcript_id}`).
2. RepeatMasker repeat elements (header from `gene_id`, with `_dup{n}` suffix for duplicates — pattern verified from hg38).
3. Simple repeats (header from `transcript_id` per Edge Case 9).
4. tRNA (header from `gene_id`).
5. miRNA (header from gff3 `Name=` attribute).
6. Custom FASTA records (headers preserved verbatim).

Missing IDs handling:
- Count: `len(expected_ids) - len(successfully_extracted_ids)`.
- If `missing / expected > 0.01`: write `{output_dir}/missing_ids_report.txt` with one ID per line; log WARNING.
- Else: log count to stderr only.

Bowtie2 build:
```python
subprocess.run(
    ["bowtie2-build", str(fa_path), str(output_dir / output_prefix)],
    check=True
)
subprocess.run(
    ["bowtie2-inspect", "--summary", str(output_dir / output_prefix)],
    check=True, stdout=subprocess.DEVNULL
)
```

Pre-check `samtools faidx` via `faidx_if_missing()` from `_shared.py`.

Write-scope guard: call `assert_writable(output_dir, allowed_prefixes=["examples/inputs/hg38/bowtie2_index", "examples/inputs/mm10", "examples/inputs/mm39", "/tmp"])` to allow hg38 reproduction validation under T-05 without modifying hg38 itself (write to `/tmp/test_hg38_bowtie2_index/`).

## Acceptance Criteria

- AC-T04-1: Script exits 0 when given valid mm10 inputs and writes 7 files (`<prefix>.1.bt2` through `.rev.2.bt2` plus `<prefix>.fa`) into `--output-dir`.
- AC-T04-2: `bowtie2-inspect --summary <output-dir>/<prefix>` exits 0 (index is valid).
- AC-T04-3: When `--trna` is omitted, script logs WARNING `"--trna not provided; tRNA entries omitted"` and continues.
- AC-T04-4: When `--gff3` is omitted, script logs WARNING `"--gff3 not provided; miRNA entries omitted"` and continues.
- AC-T04-5: When >1% of expected IDs are missing, `missing_ids_report.txt` is written to `--output-dir`.
- AC-T04-6: Script refuses to write to a path outside the allowed prefixes (exits non-zero with PermissionError message).
- AC-T04-7: Script handles `.gtf.gz` inputs (AC-16).
- AC-T04-8: NR_046233.2 records appear in the output `.fa` (when `--custom-fasta NR_046233.2.fasta` is provided) with headers preserved verbatim from the source FASTA.
- AC-T04-9: [SECURITY] Subprocess calls use list form (`["bowtie2-build", ...]`), never `shell=True`.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
mkdir -p /tmp/test_mm10_index

python bin/python/refdata_generation/generate_bowtie2_index.py \
  --gtf examples/inputs/mm10/downloaded/gencode.VM23.annotation.gtf.gz \
  --parsed-ucsc examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat \
  --repeatmasker examples/inputs/mm10/downloaded/mm10.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/mm10/downloaded/mm10.simplerepeats.tsv.gz \
  --fasta examples/inputs/mm10/downloaded/mm10.fa \
  --trna examples/inputs/mm10/downloaded/mm10.trna.tsv.gz \
  --gff3 examples/inputs/mm10/downloaded/mmu.gff3 \
  --custom-fasta examples/inputs/mm10/downloaded/NR_046233.2.fasta \
  --output-dir /tmp/test_mm10_index \
  --output-prefix MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat

ls /tmp/test_mm10_index/
bowtie2-inspect --summary /tmp/test_mm10_index/MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat | head
```

## Out of Scope

The ≥99% hg38 reproduction comparison lives in T-05.
