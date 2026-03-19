# Python + Snakemake Conversion Roadmap

Base branch: `codex/python-conversion`

## Goals
1. Understand and document CWL pipeline and Perl scripts in this repository.
2. Understand and document the provided inputs in `data/` and the PPA1 example outputs in `data/PPA1_rep1*`.
3. Generate small deterministic test files derived from the PPA1 example.
4. Convert each Perl script in `bin/perl/` to Python.
5. Convert CWL workflow steps to Snakemake incrementally and test each step output vs expected.

## Current Observations (from repository)
- Main SE workflow path for PPA1 example:
  - `cwl/wf_ecliprepmap_se.cwl`
  - nested `cwl/wf_ecliprepmap_se_1barcode.cwl`
  - tools: `map_repetitive_elements_se.cwl`, `splitbam.cwl`, `getpair.cwl`, `deduplicate.cwl`, `concatenate.cwl`, `gzip.cwl`, `combine.cwl`, `calculate_fold_change_from_parsed_files.cwl`
- Main PE workflow path:
  - `cwl/wf_ecliprepmap_pe.cwl`, `cwl/wf_ecliprepmap_pe_1barcode.cwl`
- Perl scripts to convert (`bin/perl/`):
  - `RepElement_pipeline_1dataset.pl`
  - `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
  - `duplicate_removal.pl`
  - `merge_multiple_parsed_files.simplified_20191022.pl`
  - `parse_bowtie2_output_realtime_includemultifamily_SE.pl`
  - `parse_bowtie2_output_realtime_includemultifamily_PE.pl`
  - `split_bam_to_subfiles_SEorPE.pl`
- Data location in this workspace:
  - `data` is a symlink to `/Volumes/X9Pro/Yeo/repetitive-element-pipeline`
  - Example SE run exists at `data/PPA1_rep1/` with `results/*.gz`, input YAML, and workflow log/output JSON.

## Step-by-Step Execution Plan

### Step 1: Pipeline + Script Documentation
Deliverables:
- `docs/cwl_pipeline_spec.md`: tool-by-tool and step dependency documentation for SE and PE workflows.
- `docs/perl_scripts_spec.md`: each Perl script CLI contract, inputs, outputs, side effects, and edge-case behavior notes.
- `docs/output_format_spec.md`: `.parsed`, `.reparsed.tsv`, `.nopipes.tsv`, `.withpipes.tsv`, SAM-like and mpileup-short field-level definitions inferred from files + code.

Validation:
- Cross-check docs against CWL inputs/outputs and real files in `data/PPA1_rep1/results`.

### Step 2: Input + Example Characterization
Deliverables:
- `docs/data_inventory.md`: every input artifact used by PPA1 example (paths, file types, sample line structure, compression).
- `docs/ppa1_expected_outputs_manifest.tsv`: canonical list of expected files and lightweight fingerprints (size, sha256).

Validation:
- deterministic fingerprint generation scripts and reproducible inventory command.

### Step 3: Create Small Test Fixtures
Deliverables:
- `tests/fixtures/mini/` with reduced-size SE inputs derived from PPA1 (small FASTQ/BAM/reference slices).
- `scripts/make_mini_fixtures.py` to regenerate fixtures deterministically.
- `tests/fixtures/mini/expected/` expected outputs for mini run.

Validation:
- fixture-generation test and manifest checksums.

### Step 4: Perl -> Python Conversion
Deliverables:
- New Python equivalents under `bin/python/` (one file per Perl script) with CLI parity.
- shared utilities module for parsing SAM/BAM and table records where appropriate.
- regression tests for each converted script on synthetic + mini fixture cases.

Migration strategy:
- convert one script at a time in this order:
  1) `split_bam_to_subfiles_SEorPE.pl`
  2) `merge_multiple_parsed_files.simplified_20191022.pl`
  3) `parse_bowtie2_output_realtime_includemultifamily_SE.pl`
  4) `parse_bowtie2_output_realtime_includemultifamily_PE.pl`
  5) `duplicate_removal_inline_..._simple.pl`
  6) `duplicate_removal.pl`
  7) `RepElement_pipeline_1dataset.pl` (legacy wrapper)

Validation:
- script-level golden tests compare normalized outputs with Perl baseline.

### Step 5: CWL -> Snakemake Conversion (Incremental)
Deliverables:
- `Snakefile`, `workflow/rules/*.smk`, `config/config.yaml`, `workflow/envs/*.yaml`.
- first target scope: SE pipeline matching `wf_ecliprepmap_se.cwl` + `_1barcode` behavior.
- second target scope: PE parity.

Incremental conversion/testing sequence (SE first):
1. map repetitive elements
2. split rep sam
3. split rmRep bam
4. prefix pairing logic (`getpair` equivalent)
5. deduplicate
6. concatenate + gzip outputs
7. combine parsed files
8. fold-change tables

Validation at each step:
- compare produced files to expected for mini fixture and (where feasible) PPA1 reference outputs.
- compare row counts, key columns, and checksum (exact when deterministic; tolerant numeric comparison where needed).

## Branching / PR Execution Strategy
- Keep `codex/python-conversion` as integration base branch.
- For each completed roadmap step, create branch `codex/python-conversion-stepN-<shortname>` from latest `codex/python-conversion`.
- Commit/push branch and open PR targeting `codex/python-conversion`.
- After PR creation, fast-forward/merge locally into `codex/python-conversion` before starting next step branch.

## Tooling Gaps Identified Locally
Currently found in PATH:
- present: `perl`, `samtools`, `gzip`
- missing: `snakemake`, `bowtie2`, `cwltool`

These are required to fully execute conversion/verification end-to-end.
