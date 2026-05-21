# T-12: Pipeline Integration, Wrapper, and Dry-Run Validation (AC-13, AC-14, AC-16)

<!-- DEPENDENCIES: T-10, T-11 -->
<!-- BLOCKS: none -->

## Goal

Verify the new mm10 and mm39 references work as drop-in replacements in the Snakemake dropin workflow, and confirm Perl `parse_bowtie2_output_realtime_includemultifamily_SE.pl` accepts the mm10 MASTER_FILELIST without modification. Provide the `run_assembly.sh` wrapper (FR-10) if not already created.

## Files Touched

- CREATE: `bin/python/refdata_generation/run_assembly.sh` (if not already produced in T-10)
- CREATE: `bin/python/refdata_generation/README.md` (usage notes; documents the RepElement_pipeline_1dataset.pl hardcoded-species caveat per Edge Case 7)
- READ ONLY: workflow/rules/dropin_repelement.smk, bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl

## Implementation Notes

`run_assembly.sh` content:
```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: run_assembly.sh --assembly {hg38|mm10|mm39} --date YYYYMMDD
# Generates all four reference files for the given assembly.
```

The script resolves source paths from `examples/inputs/<assembly>/downloaded/` and invokes the four Python scripts in dependency order. It auto-detects optional inputs (trna, gff3, custom fasta) and includes flags only when files exist.

Pipeline dry-run command (from architect-prompt §10):
```bash
snakemake -s workflow/rules/dropin_repelement.smk \
  --config se_or_pe=SE \
           bowtie2_db=examples/inputs/mm10/bowtie2_index \
           bowtie2_prefix=MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
           fileListFile1=examples/inputs/mm10/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv \
           gencodeGTF=examples/inputs/mm10/downloaded/gencode.VM23.annotation.gtf.gz \
           gencodeTableBrowser=examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat \
           repMaskBEDFile=examples/inputs/mm10/UniqueGenomicElements.mm10.bed \
  -n
```

Perl compat test (AC-13):
- Construct a tiny test SAM (5-10 reads) mapped against the mm10 bowtie2 index.
- Pipe through `parse_bowtie2_output_realtime_includemultifamily_SE.pl <sam> <bowtie2_db> /tmp/test_parsed.out <mm10_filelist>`.
- Assert exit 0 and that `/tmp/test_parsed.out` is non-empty.

## Acceptance Criteria

- AC-T12-1: `bin/python/refdata_generation/run_assembly.sh` exists, is executable, and accepts `--assembly` and `--date` flags.
- AC-T12-2: `run_assembly.sh --assembly hg38 --date 20260514` regenerates all four hg38 outputs into a `/tmp` scratch dir (does not overwrite `examples/inputs/hg38/`).
- AC-T12-3: `snakemake -n` against mm10 references exits 0 with a non-empty rule listing. (AC-14)
- AC-T12-4: `snakemake -n` against mm39 references exits 0.
- AC-T12-5: `parse_bowtie2_output_realtime_includemultifamily_SE.pl` runs to exit 0 against a small test SAM and the mm10 MASTER_FILELIST. (AC-13)
- AC-T12-6: All four scripts also accept `.gtf` (uncompressed) inputs in addition to `.gtf.gz` (AC-16) — verified by re-running mm10 wrapper after `gunzip --keep` on the GTF.
- AC-T12-7: README.md documents the RepElement_pipeline_1dataset.pl hardcoded-species caveat (Edge Case 7) and the date-token convention.
- AC-T12-8: [SECURITY] `run_assembly.sh` uses `set -euo pipefail`; refuses unknown `--assembly` values.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

# AC-T12-1, AC-T12-2
ls -l bin/python/refdata_generation/run_assembly.sh
bash bin/python/refdata_generation/run_assembly.sh --assembly hg38 --date 20260514 --output-root /tmp/hg38_scratch

# AC-T12-3
snakemake -s workflow/rules/dropin_repelement.smk \
  --config se_or_pe=SE \
           bowtie2_db=examples/inputs/mm10/bowtie2_index \
           bowtie2_prefix=MASTER_FILELIST.20260514.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
           fileListFile1=examples/inputs/mm10/MASTER_FILELIST.20260514.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv \
           gencodeGTF=examples/inputs/mm10/downloaded/gencode.VM23.annotation.gtf.gz \
           gencodeTableBrowser=examples/inputs/mm10/gencode.VM23.annotation.gtf.parsed_ucsc_tableformat \
           repMaskBEDFile=examples/inputs/mm10/UniqueGenomicElements.mm10.bed \
  -n 2>&1 | tail -20
# Expect exit 0

# AC-T12-5 (Perl compat) — requires a small test SAM; can be generated from any prior eCLIP test fixture
# perl bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl <test_sam> <bowtie2_db> /tmp/test.out <mm10_filelist>
```
