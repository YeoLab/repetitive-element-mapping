# T-05: hg38 bowtie2 FASTA Reproduction (AC-04)

<!-- DEPENDENCIES: T-04 -->
<!-- BLOCKS: T-10, T-11 -->

## Goal

Verify that `generate_bowtie2_index.py` reproduces the hg38 bowtie2 source FASTA at ≥99% similarity, measured by FASTA header count.

## Files Touched

- CREATE: `.forge/stages/2-architect/notes/T-05-hg38-bowtie2-headers.log`

## Implementation Notes

The hg38 reference FASTA is at:
`examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa`

Reproduction test:
1. Run `generate_bowtie2_index.py` with hg38 inputs, writing to a /tmp directory.
2. Compare FASTA header counts:
   - `grep -c "^>" <new>.fa`
   - `grep -c "^>" <reference>.fa`
3. Compute ratio: `min(new, ref) / max(new, ref) >= 0.99`.

Iteration points (if ratio < 0.99):
- Verify header naming convention for repeats (`_dup{n}` suffix pattern).
- Verify miRNA naming from `hsa.gff3` `Name=` attribute.
- Verify simple repeat naming uses `transcript_id` not `gene_id`.

## Acceptance Criteria

- AC-T05-1: `grep -c "^>" /tmp/test_hg38_index/<prefix>.fa` ≥ 0.99 * `grep -c "^>" examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa`. (AC-04)
- AC-T05-2: Header set intersection: at least 99% of reference headers appear in new output (`comm -12` of sorted header sets).
- AC-T05-3: `bowtie2-inspect --summary` on the new index returns exit 0.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

# Identify hg38 inputs
HG38_GTF=examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf
HG38_FA=examples/inputs/hg38/downloaded/hg38.fasta

mkdir -p /tmp/test_hg38_index

python bin/python/refdata_generation/generate_bowtie2_index.py \
  --gtf "$HG38_GTF" \
  --parsed-ucsc examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat \
  --repeatmasker examples/inputs/hg38/downloaded/hg38.repeatmasker.tsv.gz \
  --simplerepeats examples/inputs/hg38/downloaded/hg38.simplerepeats.tsv.gz \
  --fasta "$HG38_FA" \
  --trna examples/inputs/hg38/downloaded/hg38.trna.tsv.gz \
  --gff3 examples/inputs/hg38/downloaded/hsa.gff3 \
  --output-dir /tmp/test_hg38_index \
  --output-prefix MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat

NEW=$(grep -c "^>" /tmp/test_hg38_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa)
REF=$(grep -c "^>" examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa)
echo "new=$NEW ref=$REF ratio=$(python -c "print(min($NEW,$REF)/max($NEW,$REF))")"
# Expect: ratio >= 0.99
```
