# T-03: hg38 parsed_ucsc_tableformat Reproduction (AC-01)

<!-- DEPENDENCIES: T-02 -->
<!-- BLOCKS: T-10, T-11 -->

## Goal

Verify that `generate_parsed_ucsc_tableformat.py` reproduces `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` exactly (line-for-line identical, 249,044 lines including header). Iterate on the script until `diff` is empty.

## Files Touched

- MODIFY: `bin/python/refdata_generation/generate_parsed_ucsc_tableformat.py` (iteratively, as needed)
- CREATE: `.forge/stages/2-architect/notes/T-03-hg38-parsed-diff.log` (capture each iteration's diff sample)

## Implementation Notes

Standard reproduction-test-first loop:
1. Run script against hg38 GTF.
2. `diff /tmp/test_hg38_parsed.tsv examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat`.
3. Read first 20 differing rows; identify pattern (sort order? cds calc? trailing whitespace? exon count?).
4. Patch script. Re-run.
5. Repeat until diff is empty or implementer has iterated 3+ times — at which point escalate to user with diff sample.

Likely iteration points based on Risk Register:
- **R-01:** sort order — check whether hg38 sorts by `(gene_id, transcript_id)` or some other key (e.g., chromosomal position).
- **R-02:** `cdsStart`/`cdsEnd` — confirm CDS-feature aggregation matches.
- Trailing/leading whitespace.
- Newline at end of file (LF terminator).

## Acceptance Criteria

- AC-T03-1: **`diff <new_output> <hg38_reference>` produces zero output (empty diff).** (AC-01)
- AC-T03-2: Output line count is exactly 249,044.
- AC-T03-3: MD5 of new output matches MD5 of reference: `md5sum` on both yields identical hashes.
- AC-T03-4: Diff iteration log captured in `.forge/stages/2-architect/notes/T-03-hg38-parsed-diff.log` showing the path from first attempt to convergence.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
python bin/python/refdata_generation/generate_parsed_ucsc_tableformat.py \
  --gtf examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf \
  --output /tmp/test_hg38_parsed.tsv

diff -q /tmp/test_hg38_parsed.tsv \
        examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
# Expect: no output (files identical)

md5sum /tmp/test_hg38_parsed.tsv examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
# Expect: identical hashes
```

## Escalation Path

If after 3 implementation iterations the diff is still non-empty:
1. Stop iterating.
2. Capture a sample of differing rows in the log.
3. Surface to user: "AC-01 requires 100% identical reproduction; after N iterations residual diff is {pattern}. Proceed at ≥99% similarity threshold (AC-04/07/10 standard), or investigate the original hg38 generator?"
