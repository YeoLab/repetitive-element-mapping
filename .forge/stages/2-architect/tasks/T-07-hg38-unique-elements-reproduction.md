# T-07: hg38 UniqueGenomicElements Reproduction (AC-07)

<!-- DEPENDENCIES: T-06 -->
<!-- BLOCKS: T-10, T-11 -->

## Goal

Verify that `generate_unique_genomic_elements.py` reproduces `examples/inputs/hg38/UniqueGenomicElements.hg38.bed` within ±1% line count (target: 5,618,483 lines).

## Files Touched

- CREATE: `.forge/stages/2-architect/notes/T-07-hg38-unique-elements-diff.log`

## Acceptance Criteria

- AC-T07-1: `wc -l /tmp/test_hg38_unique.bed` is within ±1% of 5,618,483 (i.e. between 5,562,298 and 5,674,668).
- AC-T07-2: Generated file is 6-column tab-separated BED (`awk -F'\t' '{print NF}' | sort -u` yields only `6`).
- AC-T07-3: All proximal entries have name suffix `-proximal`.
- AC-T07-4: No row has start < 0.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

python bin/python/refdata_generation/generate_unique_genomic_elements.py \
  --repeatmasker examples/inputs/hg38/downloaded/hg38.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/hg38/downloaded/hg38.simplerepeats.tsv.gz \
  --trna examples/inputs/hg38/downloaded/hg38.trna.tsv.gz \
  --gff3 examples/inputs/hg38/downloaded/hsa.gff3 \
  --parsed-ucsc examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat \
  --assembly hg38 \
  --output /tmp/test_hg38_unique.bed

NEW=$(wc -l < /tmp/test_hg38_unique.bed)
REF=5618483
RATIO=$(python -c "print(min($NEW,$REF)/max($NEW,$REF))")
echo "new=$NEW ref=$REF ratio=$RATIO"
# Expect: ratio >= 0.99
```
