# T-06: Implement `generate_unique_genomic_elements.py`

<!-- DEPENDENCIES: T-01, T-02 -->
<!-- BLOCKS: T-07, T-10, T-11 -->

## Goal

Implement the Python script that emits a 6-column BED of repeat/tRNA/miRNA elements plus 500 bp proximal flanks. Must handle missing optional inputs gracefully.

## Files Touched

- CREATE: `bin/python/refdata_generation/generate_unique_genomic_elements.py`

## Implementation Notes

CLI signature:
```
--repeatmasker PATH         required, .tsv.gz
--simplerepeats PATH        optional, .tsv.gz
--trna PATH                 optional, .tsv.gz
--gff3 PATH                 optional, miRBase gff3
--parsed-ucsc PATH          optional, parsed_ucsc_tableformat
--assembly NAME             required, used in WARNING messages
--output PATH               required, destination BED
[--flank N]                 optional, default 500
```

Algorithm: see Architecture Plan §4.5.

Output format: tab-separated 6-column BED:
```
chrom    start    end    name    score    strand
```

No header. LF line endings. No trailing whitespace.

Proximal flank rule:
- Upstream: `(chrom, max(0, start - flank), start, name + "-proximal", 0, strand)`.
- Downstream: `(chrom, end, end + flank, name + "-proximal", 0, strand)`.

Sort order: `(chrom, start, end, name)` — natural string sort for chrom (e.g. chr1 < chr10 < chr2 lexicographic, matching `sort` default).

When an optional input is absent, log WARNING with `--assembly` value:
```
WARNING [generate_unique_genomic_elements]: --trna not provided; tRNA entries omitted from {assembly} UniqueGenomicElements
```

## Acceptance Criteria

- AC-T06-1: Output is a tab-separated 6-column BED file.
- AC-T06-2: For each input element, exactly two `-proximal` rows are written (upstream and downstream).
- AC-T06-3: Upstream proximal start is `max(0, start - 500)`.
- AC-T06-4: When `--trna` is omitted, script logs WARNING and produces a smaller BED file (no tRNA entries).
- AC-T06-5: When `--gff3` is omitted, script logs WARNING (no miRNA entries).
- AC-T06-6: Simple repeats use `transcript_id` (not `gene_id`) as the name field (Edge Case 9).
- AC-T06-7: Output sorted by (chrom, start, end, name).
- AC-T06-8: Script exits 0 on all required-input combinations from AC-T06-3 through AC-T06-7.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

# mm10 with all inputs
python bin/python/refdata_generation/generate_unique_genomic_elements.py \
  --repeatmasker examples/inputs/mm10/downloaded/mm10.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/mm10/downloaded/mm10.simplerepeats.tsv.gz \
  --trna examples/inputs/mm10/downloaded/mm10.trna.tsv.gz \
  --gff3 examples/inputs/mm10/downloaded/mmu.gff3 \
  --parsed-ucsc examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat \
  --assembly mm10 \
  --output /tmp/test_mm10.bed
awk -F'\t' '{print NF}' /tmp/test_mm10.bed | sort -u
# Expect: only "6"

# mm39 without tRNA/gff3 — should log two WARNINGs
python bin/python/refdata_generation/generate_unique_genomic_elements.py \
  --repeatmasker examples/inputs/mm39/downloaded/mm39.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/mm39/downloaded/mm39.simplerepeats.tsv.gz \
  --parsed-ucsc examples/inputs/mm39/gencode.VM38.annotation.gtf.parsed_ucsc_tableformat \
  --assembly mm39 \
  --output /tmp/test_mm39.bed 2>&1 | grep WARNING
# Expect: 2 WARNING lines
```
