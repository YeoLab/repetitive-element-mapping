# Architect Prompt: Snakemake Workflow for eCLIP Repetitive Element Mapping

STATUS: INITIALIZING

---

## 1. Architect Role Definition

You are a senior bioinformatics workflow engineer. Your task is to design a complete Snakemake implementation of the eCLIP repetitive element mapping pipeline, translating from existing CWL workflows while preserving exact output fidelity. You will produce an implementer-ready task decomposition with precise acceptance criteria.

---

## 2. User Request

Translate the CWL-based eCLIP repetitive element mapping pipeline to a Snakemake workflow. The pipeline maps trimmed eCLIP FASTQ files to a curated repeat element database, deduplicates by UMI prefix, resolves conflicts between unique genomic and repeat-family mappings, and produces per-element read counts and fold enrichment (IP vs. Input). Both paired-end (PE) and single-end (SE) workflows must be supported from a single Snakefile with mode selection via config.

---

## 3. Mission Brief

**Scope:**
- Translate all CWL tools and workflows to Snakemake rules/modules
- Generate downsampled test datasets from example data
- Produce config.yaml with SE/PE mode selection and enforce appropriate input structure
- Run full pipeline on both downsampled and full datasets
- Validate that Snakemake outputs match CWL reference outputs

**In scope:**
- Snakefile with SE/PE dispatch logic
- workflow/rules/SE.smk, PE.smk, common.smk modularization
- Config schema with SE vs PE input validation
- Downsampled test datasets in examples/inputs/downsampled/
- SLURM profiles (profiles/tscc2_snakemake9/)
- Conda environment (workflow/envs/)
- README and changelog updates

**Out of scope:**
- Modifying Perl scripts (unless strictly required and approved)
- Docker/Singularity container changes
- New reference data generation (covered by existing prompts)

**Success criteria:**
- Snakemake workflow completes successfully on downsampled and full datasets
- Output .nopipes.tsv and .withpipes.tsv match CWL reference outputs exactly (or within random tie-breaking tolerance)
- SE and PE configs both validated and tested

---

## 4. Technical Context

See context files 10-technical-constraints.md and 11-code-references.md for full details.

**Key facts:**
- Language: Python/Perl; Framework: Snakemake 9.12.0 + CWL (existing)
- SLURM cluster: TSCC2, partition gold, account csd792
- Perl version sensitivity: use system perl 5.10.1 only (at /tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl)
- 25 UMI prefix bins (AA..NN) drive scatter; deduplication may need 32GB memory
- CWL reference outputs in test-provenance/tests/ecliprepmap-1.0.0/

---

## 5. Requirements

### Functional Requirements

**FR-01 [P0]:** Implement `step_map_repetitive_elements` — run `parse_bowtie2_output_realtime_includemultifamily_PE.pl` or `_SE.pl` via bowtie2, producing a Rep.sam file. Inputs: r1 (and r2 for PE), bowtie2_db directory, bowtie2_prefix string, fileListFile1.

**FR-02 [P0]:** Implement `step_splitbam` — run `split_bam_to_subfiles_SEorPE.pl` on both Rep.sam and rmRep.bam, producing 25 `.tmp` files per input split by first 2 nt of UMI.

**FR-03 [P0]:** Implement `step_deduplicate` — run `duplicate_removal_inline_paired...pl` once per UMI prefix (25 scatter jobs), producing `.rmDup.sam`, `.prermDup.sam`, and `.parsed_v2.20201210.txt` per prefix. Inputs: matched rep.tmp + rmrep.tmp pair, se_or_pe, gencodeGTF, gencodeTableBrowser, repMaskBEDFile, fileListFile1.

**FR-04 [P0]:** Implement `step_concatenate` + `step_gzip` — cat all 25 `.rmDup.sam` files into one, gzip; repeat for `.prermDup.sam`.

**FR-05 [P0]:** Implement `step_combine_parsed` — run `merge_multiple_parsed_files.simplified_20191022.pl` to merge 25 `.parsed_v2.20201210.txt` into one `.parsed` file.

**FR-06 [P0]:** Implement `step_calculate_fold_change` — run `calculate_fold_change_from_parsed_files.py`, producing `.nopipes.tsv` and `.withpipes.tsv` from IP and Input `.parsed` files.

**FR-07 [P0]:** PE workflow — run the per-barcode sub-workflow for barcode1, barcode2, and input in parallel, then merge barcode1+barcode2 parsed files with `step_combine_parsed` before fold change calculation.

**FR-08 [P0]:** SE workflow — run per-barcode sub-workflow for barcode1 and input, then fold change directly from their `.parsed` files.

**FR-09 [P0]:** Config-driven SE/PE dispatch — single Snakefile with logic to select SE or PE sub-rules from config. If `se_or_pe: PE`, barcode2 fields are required; if `se_or_pe: SE`, barcode2 fields are forbidden.

**FR-10 [P1]:** Generate downsampled test datasets — PE and SE, all chromosomes, ≥100 reads per UMI 2-nt prefix, BAM reads must exist in corresponding FASTQ.

**FR-11 [P1]:** Generate small config YAMLs — `repeat_mapping_PE_small.yaml` and `repeat_mapping_SE_small.yaml` pointing to downsampled data.

**FR-12 [P1]:** Profile memory/CPU for `step_deduplicate` on full dataset; if ≤32GB, remove scatter and re-test; if >32GB, retain scatter.

**FR-13 [P2]:** Clean up — remove unused CWL and Perl scripts after successful validation; commit and push before deleting.

**FR-14 [P2]:** Update README and changelog with Snakemake usage and downsampled data instructions.

### Non-Functional Requirements

**NFR-01:** Snakemake version == 9.12.0 exactly.
**NFR-02:** Use system perl 5.10.1 (`/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl`) for all Perl scripts.
**NFR-03:** Default resources: 1 CPU, 32 GB memory, 8h walltime per rule.
**NFR-04:** Retry limit: 2; memory scales with attempt (`attempt * 32000` MB or similar).
**NFR-05:** Downsampled test data must be lightweight (complete within default resources).
**NFR-06:** SLURM profile: partition gold, account csd792.

### Security Requirements

**SEC-01:** No credentials in config files or Snakefiles.
**SEC-02:** All file paths in config (not hardcoded in rules).

---

## 6. Acceptance Criteria

**AC-01:** Snakemake dry-run (`snakemake -n`) completes without error for both SE and PE configs.

**AC-02:** Downsampled PE pipeline completes end-to-end; output `.nopipes.tsv` and `.withpipes.tsv` match CWL reference outputs from test-provenance/ (or are within random tie-breaking tolerance).

**AC-03:** Downsampled SE pipeline completes end-to-end; output `.nopipes.tsv` and `.withpipes.tsv` match CWL SE reference outputs.

**AC-04:** Full (non-downsampled) PE pipeline produces output that matches CWL PE reference outputs.

**AC-05:** Full (non-downsampled) SE pipeline produces output that matches CWL SE reference outputs.

**AC-06:** All 25 UMI prefix scatter jobs complete; `.parsed` files merge correctly.

**AC-07:** PE config with missing barcode2 fields raises a config validation error before execution.

**AC-08:** SE config with barcode2 fields present raises a config validation error.

**AC-09:** Downsampled datasets contain ≥100 reads per UMI 2-nt prefix (all 25 prefixes covered) in both SE and PE data.

**AC-10:** All BAM reads in downsampled files exist in their corresponding downsampled FASTQ file.

---

## 7. Assumptions and Defaults

- Perl scripts are used as-is; no modifications unless the script literally cannot work with Snakemake.
- `getpair` CWL ExpressionTool is replaced by Snakemake wildcard logic matching `.tmp` filenames by prefix.
- Working directory for each rule is the Snakemake working directory; scripts are referenced by absolute path or in PATH via the conda env.
- If deduplication memory profiling shows ≤32GB for the full dataset, scatter is removed from the workflow (25 serial rules replaced by 1 rule processing all prefixes).
- Outputs land in a configurable `output_dir` (or alongside inputs); exact structure TBD during architecture.

---

## 8. Open Gaps Ledger

Critical gaps requiring architect decision:

| ID | Gap | Category |
|----|-----|----------|
| G-01 | Output directory structure: flat per-barcode dirs or nested? | Architecture |
| G-02 | How to handle `InitialWorkDirRequirement` from combine.cwl (files must be in cwd) | Implementation |
| G-03 | Wildcard naming convention for the 25 tmp files — how to match rep and rmrep pairs by prefix | Implementation |
| G-04 | Whether to profile memory for deduplication before or after implementing scatter | Sequencing |
| G-05 | Config schema: nested keys (barcode1: {r1: ...}) vs flat keys (barcode1r1FastqGz:) | UX |

---

## 9. Architect Decision Checklist

- [ ] Directory/naming convention for all intermediate and final outputs
- [ ] Config schema design (flat vs. nested; SE vs. PE required fields)
- [ ] How `split_bam_to_subfiles_SEorPE.pl` names its output .tmp files (verify by inspection)
- [ ] Whether scatter jobs run in SLURM array or as sequential Snakemake jobs
- [ ] Memory profiling strategy for step_deduplicate
- [ ] Conda environment contents for workflow/envs/dropin.yaml (or equivalent)

---

## 10. Verification Environment

```bash
# Activate environment
module load singularitypro; conda activate snakemake9

# Dry-run SE
snakemake -s Snakefile --config se_or_pe=SE -n

# Dry-run PE
snakemake -s Snakefile --config se_or_pe=PE -n

# Run downsampled SE
snakemake -s Snakefile --configfile examples/inputs/downsampled/repeat_mapping_SE_small.yaml --profile profiles/tscc2_snakemake9

# Compare outputs
diff <snakemake_output>.nopipes.tsv test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_se/wf_ecliprepmap_se/results/INV_B.IP.umi.r1.fqTrTr.sorted.fq.barcode1.nopipes.tsv

# System perl check
/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl --version
```

---

## 11. Context Files Reference

- `context/01-vision-and-goals.md` — Project purpose and success definition
- `context/02-user-experience.md` — User-facing interface (config YAML, CLI invocation)
- `context/03-user-flows.md` — Step-by-step execution flows for SE and PE
- `context/04-data-models.md` — Input/output file formats and naming conventions
- `context/05-business-logic.md` — Pipeline step details, resource requirements, scatter logic
- `context/06-api-integrations.md` — External tool interfaces (bowtie2, Perl scripts, Python scripts)
- `context/07-security-requirements.md` — N/A (HPC batch pipeline, no network services)
- `context/08-edge-cases.md` — Empty prefix bins, missing files, Perl version sensitivity
- `context/09-acceptance-criteria.md` — Detailed AC with verification commands
- `context/10-technical-constraints.md` — Environment, SLURM, software versions
- `context/11-code-references.md` — Key file paths and code snippets
