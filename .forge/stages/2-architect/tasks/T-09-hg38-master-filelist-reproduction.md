# T-09: hg38 MASTER_FILELIST Reproduction (AC-10)

<!-- DEPENDENCIES: T-08 -->
<!-- BLOCKS: T-10, T-11 -->

## Goal

Verify that `generate_master_filelist.py` reproduces `examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` within ±1% line count (target: 26,422 lines).

## Files Touched

- CREATE: `.forge/stages/2-architect/notes/T-09-hg38-master-filelist-diff.log`

## Acceptance Criteria

- AC-T09-1: `wc -l /tmp/test_hg38_filelist.tsv` is within ±1% of 26,422 (between 26,158 and 26,686).
- AC-T09-2: Output is 5-column tab-separated, no header.
- AC-T09-3: First 100 lines of the new output match the family/labeling convention of the hg38 reference (sample inspection).
- AC-T09-4: `.list` sibling exists with identical content (cmp returns 0).

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

python bin/python/refdata_generation/generate_master_filelist.py \
  --gtf examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf \
  --parsed-ucsc examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat \
  --repeatmasker examples/inputs/hg38/downloaded/hg38.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/hg38/downloaded/hg38.simplerepeats.tsv.gz \
  --trna examples/inputs/hg38/downloaded/hg38.trna.tsv.gz \
  --gff3 examples/inputs/hg38/downloaded/hsa.gff3 \
  --output /tmp/test_hg38_filelist.tsv

NEW=$(wc -l < /tmp/test_hg38_filelist.tsv)
REF=$(wc -l < examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv)
RATIO=$(python -c "print(min($NEW,$REF)/max($NEW,$REF))")
echo "new=$NEW ref=$REF ratio=$RATIO"
# Expect: ratio >= 0.99

cmp /tmp/test_hg38_filelist.tsv /tmp/test_hg38_filelist.list && echo "TSV/list match"
```
