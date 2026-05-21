# Architect Prompt — Reference Data Generation for mm10 and mm39

**STATUS:** READY_FOR_ARCHITECT  
**Stage:** 1-requirements  
**Feature:** snakemake workflow for repetitive element mapping  
**Prepared:** 2026-05-14  

---

## 1. Architect Role Definition

You are a senior bioinformatics software architect. Your task is to design a complete, production-ready implementation plan for generating reference data for two new mouse genome assemblies (mm10 and mm39) in the eCLIP repetitive element mapping pipeline. The pipeline already ships reference data for hg38; this work extends it to support mouse.

You will produce:
1. An architecture decision record explaining your approach choices
2. A decomposed task list with clear acceptance criteria per task
3. An implementer-ready prompt with all context embedded

---

## 2. User Request

> Generate reference files for mice species mm10 and mm39 — specifically: `parsed_ucsc_tableformat`, `bowtie2_index`, `UniqueGenomicElements.{assembly}.bed`, and `MASTER_FILELIST.*.tsv` — using the same pipeline logic as hg38 but adapted for mouse annotations.
> 
> New Python scripts should be reusable for any future assembly. Perl scripts may be modified only if hardcoded values prevent mm10/mm39 compatibility.

Source prompt: `prompts/generate_refdata.md`

---

## 3. Mission Brief

### In Scope
- Four Python scripts (one per reference file type) that accept assembly-generic CLI arguments
- Reproduction of all four hg38 reference files using the new scripts (validation step)
- Generation of all four reference files for mm10 and mm39
- Minimal Perl script modifications (only if hardcoded values block compatibility)
- Documentation of any missing source files or IDs

### Out of Scope
- Changes to CWL workflow files
- Changes to Snakemake dropin rules
- Adding new conda packages to `workflow/envs/dropin.yaml` (unless unavoidable)
- Reference generation for assemblies other than mm10 and mm39
- New pipeline features or output formats

### Success Definition
1. All 4 reference files generated for mm10 — scripts exit 0
2. All 4 reference files generated for mm39 — scripts exit 0
3. hg38 reproduction test: parsed_ucsc_tableformat = 100% identical; MASTER_FILELIST and UniqueGenomicElements ≥99% line similarity; bowtie2 FASTA header count ≥99%
4. `snakemake -n` (dry-run) passes with mm10 and mm39 references
5. Any Perl modifications are documented and minimal

---

## 4. Technical Context

### Repository Facts
- **Repo root:** `/tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping`
- **Language:** Python 3.11 + Perl (≤5.16 for determinism)
- **Pipeline implementations:** CWL (production) and Snakemake dropin (both call same Perl scripts)
- **Scheduler:** SLURM on TSCC, profile at `profiles/tscc2_snakemake9/`
- **Module:** `module load ecliprepmap/1.0.0` provides all tools in PATH
- **Conda env:** `workflow/envs/dropin.yaml` (Python 3.11, bowtie2 ≥2.5, samtools ≥1.17, pybedtools, pandas, numpy)
- **Docker image:** `brianyee/repetitive_element_mapping:1.0.0` (for CWL external runs)

### Key File Paths
See `context/11-code-references.md` for the complete file inventory. Critical paths:
- Ground truth hg38 parsed: `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` (249,044 lines)
- Ground truth MASTER_FILELIST: `examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` (26,422 lines)
- Ground truth UniqueGenomicElements: `examples/inputs/hg38/UniqueGenomicElements.hg38.bed` (5,618,483 lines)
- Ground truth bowtie2 FASTA: `examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa`

### Source Input Files Available

| File type | hg38 | mm10 | mm39 |
|-----------|------|------|------|
| GTF | v33 (uncompressed symlink) | VM23.annotation.gtf.gz | VM38.annotation.gtf.gz |
| Genome FASTA | symlinked | MISSING | MISSING |
| repeatmasker.tsv.gz | present | present | present |
| simplerepeats.tsv.gz | present | present | present |
| trna.tsv.gz | present | present | ABSENT |
| .gff3 (miRNA) | hsa.gff3 | mmu.gff3 | ABSENT |
| NR_046233.2.fasta | absent | present | present |

### Perl Scripts — Hardcoded Value Assessment
- `parse_bowtie2_output_realtime_includemultifamily_SE/PE.pl`: hg38-specific paths are **commented out**; active code takes MASTER_FILELIST path via ARGV[3]. **No modification needed.**
- `split_bam_to_subfiles_SEorPE.pl`: hg38 mentioned only in a comment on line 10. **No modification needed.**
- `RepElement_pipeline_1dataset.pl`: `my $species = "hg38"` hardcoded on line 4; hardcoded bowtie_db path for hg38. **This script is NOT called by CWL or Snakemake dropin.** Modification not required for pipeline compatibility.
- `duplicate_removal_inline_paired...pl`: No hg38 hardcoding detected. **No modification needed.**

### Coordinate System
- All input TSV files use 1-based inclusive coordinates (GTF-style)
- All output files use 0-based half-open coordinates (BED/UCSC-style)
- Conversion: `start_0 = gtf_start - 1`, `end_0 = gtf_end`

---

## 5. Requirements

### Functional Requirements

**FR-01 (CRITICAL):** `generate_parsed_ucsc_tableformat.py` must accept `--gtf` (`.gtf` or `.gtf.gz`) and `--output` arguments; produce the 11-column tab-separated UCSC table format with header line `#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds`; one row per unique transcript_id where GTF feature == "transcript"; `exonCount` = count of "exon" features; `exonStarts`/`exonEnds` comma-delimited with trailing comma; 0-based coordinates.

**FR-02 (CRITICAL):** `generate_bowtie2_index.py` must accept `--gtf`, `--parsed-ucsc`, `--repeatmasker`, `--simplerepeats`, `--fasta` (required), `--trna` (optional), `--gff3` (optional), `--custom-fasta` (optional, repeatable), `--output-dir`, `--output-prefix`; produce combined FASTA and bowtie2 index; handle missing IDs per Edge Case 2 (see `context/08-edge-cases.md`).

**FR-03 (CRITICAL):** `generate_unique_genomic_elements.py` must accept `--repeatmasker` (required), `--simplerepeats` (optional), `--trna` (optional), `--gff3` (optional), `--parsed-ucsc` (optional), `--assembly`, `--output`; produce 6-column BED with 500 bp proximal flanks per element; log WARNINGs for absent optional inputs; never silently drop categories.

**FR-04 (CRITICAL):** `generate_master_filelist.py` must accept `--parsed-ucsc` (required), `--repeatmasker` (required), `--simplerepeats` (optional), `--trna` (optional), `--gff3` (optional), `--custom-fasta` (optional), `--output`; produce 5-column tab-separated TSV with no header; NR_ custom sequences mapped to correct family labels.

**FR-05 (HIGH):** All scripts must handle `.gtf.gz` inputs by auto-detecting from file extension and decompressing via `gzip.open()` or `subprocess`.

**FR-06 (HIGH):** `generate_bowtie2_index.py` must call `bowtie2-build` as a subprocess after assembling the FASTA; verify the index with `bowtie2-inspect --summary` and log success/failure.

**FR-07 (HIGH):** All scripts must print progress to stderr (source file read, records processed, output written); all errors must produce a descriptive message and exit non-zero.

**FR-08 (MEDIUM):** Missing IDs exceeding 1% must produce `missing_ids_report.txt` in the output directory; missing IDs under 1% must be logged to stderr only.

**FR-09 (MEDIUM):** All four scripts must produce deterministic output (same input → same output, regardless of Python dict iteration order or OS).

**FR-10 (LOW):** A wrapper script or Makefile target that runs all four steps in order for a given assembly, with configurable paths.

### Non-Functional Requirements

**NFR-01:** Scripts run within `workflow/envs/dropin.yaml` conda environment without new package additions.

**NFR-02:** Total reference generation time for one assembly ≤ 2 hours on a standard TSCC compute node (20 GB RAM, 4 CPUs).

**NFR-03:** Output files use tab (`\t`) delimiters, Unix line endings (`\n`), no trailing whitespace.

**NFR-04:** Scripts must be runnable non-interactively (no prompts); all parameters via CLI flags with argparse.

**NFR-05:** All output files include the generating script name and timestamp in a header comment where the format permits (e.g., `# Generated by generate_parsed_ucsc_tableformat.py on YYYY-MM-DD`). Exception: formats with strict no-header requirements (e.g., plain 5-column TSV for MASTER_FILELIST) should log this to stderr instead.

### Security Requirements

See `context/07-security-requirements.md`. Summary: no network access, no sensitive data, standard HPC file permissions. Scripts must write only to `examples/inputs/mm10/` or `examples/inputs/mm39/`.

---

## 6. Acceptance Criteria

See `context/09-acceptance-criteria.md` for the full 16 acceptance criteria. Summary:

| ID | Description | Priority |
|----|-------------|----------|
| AC-01 | hg38 parsed_ucsc_tableformat reproduced 100% identically | CRITICAL |
| AC-02 | mm10 parsed_ucsc_tableformat generated with correct structure | CRITICAL |
| AC-03 | mm39 parsed_ucsc_tableformat generated | CRITICAL |
| AC-04 | hg38 bowtie2 FASTA reproduced ≥99% by header count | CRITICAL |
| AC-05 | mm10 bowtie2_index generated, functional | CRITICAL |
| AC-06 | mm39 bowtie2_index generated, functional | CRITICAL |
| AC-07 | hg38 UniqueGenomicElements reproduced within ±1% line count | HIGH |
| AC-08 | mm10 UniqueGenomicElements generated with correct 6-col BED format | CRITICAL |
| AC-09 | mm39 UniqueGenomicElements generated without tRNA/miRNA, with WARNINGs | HIGH |
| AC-10 | hg38 MASTER_FILELIST reproduced within ±1% line count | HIGH |
| AC-11 | mm10 MASTER_FILELIST generated with NR_046233.2, tRNA, miRNA entries | CRITICAL |
| AC-12 | mm39 MASTER_FILELIST generated with NR_046233.2, without tRNA/miRNA | CRITICAL |
| AC-13 | Perl SE/PE parse scripts run successfully with mm10 MASTER_FILELIST | HIGH |
| AC-14 | `snakemake -n` dry-run passes with mm10 references | HIGH |
| AC-15 | Missing ID warning fires and report written when >1% missing | MEDIUM |
| AC-16 | Scripts accept both .gtf and .gtf.gz inputs | HIGH |

---

## 7. Assumptions and Defaults

| Assumption | Rationale | If Wrong |
|-----------|-----------|----------|
| mm10/mm39 genome FASTAs must be provided by the user as CLI arguments | Not present in `downloaded/`; too large to store in repo | Scripts fail with clear error message pointing to genome FASTA requirement |
| mm39 will not have tRNA or miRNA entries in its reference files | `mm39.trna.tsv.gz` and `mm39.gff3` are absent from `downloaded/` | If files appear later, scripts already accept them as optional args |
| `RepElement_pipeline_1dataset.pl` need not be modified | It is not invoked by CWL or Snakemake dropin | If user runs it directly for mm10/mm39, they'll need to update the hardcoded species variable |
| Output MASTER_FILELIST column order follows existing hg38 convention | Perl parse scripts load the full 5-column file into a hash indexed by col1 | If Perl script reads specific column positions differently, column order must be re-verified |
| `parsed_ucsc_tableformat` sort order is by gene_id then transcript_id | Matches visual inspection of hg38 reference | Run diff to confirm before writing mouse outputs |
| NR_046233.2 represents the same rRNA sequences for mouse as NR_046235.3 does for human | Both are ribosomal RNA precursor sequences | Validate by checking NCBI records for NR_046233.2 |
| bowtie2 index prefix naming: `MASTER_FILELIST.{YYYYMMDD}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat` | Matches hg38 convention | Adjust if CWL/Snakemake YAML fields require a fixed prefix |

---

## 8. Open Gaps Ledger

**Critical (blocks implementation): 1**

| ID | Gap | Impact | Resolution path |
|----|-----|--------|-----------------|
| OG-01 | Mouse genome FASTA (mm10, mm39) files are not present in `examples/inputs/{assembly}/downloaded/`. The genome FASTA is required by `generate_bowtie2_index.py` for pybedtools getfasta to extract transcript sequences. | Blocks Step 2 (bowtie2 index generation) for both assemblies | User must provide paths to mm10 and mm39 reference FASTAs. Suggested sources: UCSC GRCm38/GRCm39 no-alt analysis set. The scripts should accept the path as a required `--fasta` CLI argument and fail clearly if not found. |

**High (meaningful risk, workaround exists): 2**

| ID | Gap | Impact | Resolution path |
|----|-----|--------|-----------------|
| OG-02 | Exact sort order and cdsStart/cdsEnd derivation logic for `parsed_ucsc_tableformat` must be reverse-engineered from the hg38 reference file. | Script may produce correct column content but incorrect row order or CDS boundaries. | Implementer must run reproduction test against hg38 and iterate until diff = 0 before proceeding to mouse. |
| OG-03 | MASTER_FILELIST row construction logic (which elements get included, in what order, with what family labels) is inferred from the 26,422-row hg38 reference. If the construction logic has assembly-specific elements (e.g., a hardcoded list of RNA families), those need identification. | mm10/mm39 MASTER_FILELIST may be structurally valid but fail Perl script compatibility. | Implementer must trace from bowtie2 FASTA headers → MASTER_FILELIST rows → Perl parse script's `read_in_filelists()` function to confirm the mapping is general. |

**Low (cosmetic/documentation): 1**

| ID | Gap | Impact | Resolution path |
|----|-----|--------|-----------------|
| OG-04 | The `.list` and `.tsv` extensions for MASTER_FILELIST appear to have identical content in hg38. Need to confirm which extension the Perl script's `read_in_filelists()` requires. | If Perl expects `.list`, the generated `.tsv` must also be copied or symlinked. | Check `ARGV[3]` usage in `parse_bowtie2_output_realtime_includemultifamily_SE.pl`'s `read_in_filelists` sub; generate both or symlink. |

---

## 9. Architect Decision Checklist

- [ ] Choose script placement: `bin/python/` (alongside shims) or `bin/python/refdata_generation/` (separate subdir)?
- [ ] Choose GTF parsing approach: pure Python vs. pybedtools vs. HTSeq/pyranges?
- [ ] Choose FASTA extraction approach for transcript sequences: pybedtools getfasta (BED→FASTA per region then concatenate by transcript) vs. HTSeq vs. bedtools getfasta subprocess?
- [ ] Determine NR_046233.2 naming convention: extract sub-sequences for 18S/28S/45S regions, or include the full sequence once?
- [ ] Determine MASTER_FILELIST family label format for mouse-specific tRNA/repeat families (may differ from hg38 labeling)
- [ ] Determine whether a wrapper script (FR-10) is a shell script, Python, or Snakemake rule
- [ ] Confirm `cdsStart`/`cdsEnd` derivation: use `CDS` feature rows from GTF, or use `start_codon`/`stop_codon`, or fall back to txStart/txEnd?
- [ ] Confirm sort order for parsed_ucsc_tableformat by diffing with hg38 reference
- [ ] Assess whether `pybedtools` getfasta is strand-aware by default (verify `-s` flag behavior for sense/antisense transcript extraction)
- [ ] Determine proximal flank behavior at chromosome boundaries (cap at 0 for upstream, cap at chrom_size for downstream?)

---

## 10. Verification Environment

```bash
# Environment setup
module load ecliprepmap/1.0.0
# OR
conda activate <dropin-env-name>

# Tool availability check
python --version        # expect 3.11.x
bowtie2 --version       # expect ≥2.5
bedtools --version      # expect ≥2.30
samtools --version      # expect ≥1.17
python -c "import pybedtools; print(pybedtools.__version__)"

# Reference file line count verification
wc -l examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
# expected: 249044

wc -l examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv
# expected: 26422

wc -l examples/inputs/hg38/UniqueGenomicElements.hg38.bed
# expected: 5618483

grep "^>" examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa | wc -l
# record this count as the baseline for AC-04

# Snakemake dry-run (once references are generated)
snakemake -s workflow/rules/dropin_repelement.smk \
  --config se_or_pe=SE \
  --config bowtie2_db=examples/inputs/mm10/bowtie2_index \
  -n 2>&1 | tail -20
```

---

## 11. Context Files

| File | Contents |
|------|---------|
| `context/01-vision-and-goals.md` | Project vision, goals, success criteria, out-of-scope |
| `context/02-user-experience.md` | Primary user, invocation pattern, UX expectations |
| `context/03-user-flows.md` | 6 user flows covering all scenarios |
| `context/04-data-models.md` | All input/output file formats with column schemas |
| `context/05-business-logic.md` | Step-by-step logic for all 5 steps from generate_refdata.md |
| `context/06-api-integrations.md` | External tools (bowtie2, pybedtools, samtools) |
| `context/07-security-requirements.md` | Minimal HPC security constraints |
| `context/08-edge-cases.md` | 9 edge cases with expected behavior |
| `context/09-acceptance-criteria.md` | 16 acceptance criteria with GIVEN/WHEN/THEN |
| `context/10-technical-constraints.md` | Environment, performance, naming, format constraints |
| `context/11-code-references.md` | All key file paths, script inventory, verification commands |
