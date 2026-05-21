# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

This is the **eCLIP repetitive element mapping** pipeline (ecliprepmap). Given trimmed eCLIP FASTQ files and a BAM of reads already mapped to the unique genome (rmRep BAM), it:
1. Maps reads to a curated repeat element database using Bowtie2
2. Splits reads by UMI prefix (AA, AC, ..., NN — 25 bins) for memory-efficient deduplication
3. Deduplicates using randomer UMIs, resolving conflicts between unique genomic and repeat-family mappings
4. Outputs per-element read counts and fold enrichment (IP vs. Input)

## Running the Pipeline

### CWL Workflow (Yeo Lab / TSCC)

```bash
module load ecliprepmap
# Paired-end (2 IP barcodes + 1 Input)
eCLIP_repelement_PE examples/repeat_mapping_PE.yaml

# Single-end (1 IP barcode + 1 Input)
eCLIP_repelement_SE examples/repeat_mapping_SE.yaml
```

The launcher creates `<job_name>/results/` next to the YAML file and removes `.tmp/` intermediates on success. Logs go to `<job_name>/<PIPELINE>_<job_name>_LOG.txt`.

### CWL Workflow (External)

```bash
cwltool --debug \
  --outdir <outdir> \
  --cachedir <cachedir> \
  cwl/wf_ecliprepmap_pe.cwl \
  examples/repeat_mapping_PE.yaml
```

Docker image: `brianyee/repetitive_element_mapping:1.0.0`

### Snakemake Workflow (TSCC)

```
TODO
```

On TSCC with SLURM, use the profile:

```bash
snakemake \
  --snakefile /path/to/Snakefile \
  --profile profiles/tscc2_snakemake9 \
  --config ...
```

The SLURM profile (`profiles/tscc2_snakemake9/`) defaults to partition `gold`, account `csd792`, 20 GB memory, 30-minute wall time. The dropin workflow uses `--use-conda` with `workflow/envs/dropin.yaml` (Python 3.11, bowtie2 ≥2.5, samtools ≥1.17, numpy, pandas).


## Workflow Architecture

### CWL Entry Points

| Launcher (`wf/`) | Top-level CWL | Use case |
|---|---|---|
| `eCLIP_repelement_PE` | `cwl/wf_ecliprepmap_pe.cwl` | 2 IP barcodes + 1 Input (PE) |
| `eCLIP_repelement_SE` | `cwl/wf_ecliprepmap_se.cwl` | 1 IP barcode + 1 Input (SE) |
| `eCLIP_repelement_PE_singleNode` | `cwl/wf_ecliprepmap_pe_1barcode.cwl` | 1 barcode only |
| `eCLIP_repelement_SE_singleNode` | `cwl/wf_ecliprepmap_se_1barcode.cwl` | 1 barcode only |

### Per-barcode sub-workflow (`wf_ecliprepmap_*_1barcode.cwl`)

```
map_repetitive_elements  →  splitbam (rep SAM)   ┐
                                                   ├→ getpair (25×) → deduplicate (25×, scattered) → concatenate → gzip
rmRepBam input           →  splitbam (rmrep BAM)  ┘                                                   combine_parsed
```

1. **`map_repetitive_elements_pe/se.cwl`** — Runs Bowtie2 via `parse_bowtie2_output_realtime_includemultifamily_PE/SE.pl` in streaming mode (`-q --sensitive -a -p 3 --no-mixed --reorder`), producing a SAM-like file with the best-scoring assignment per read pair.
2. **`splitbam.cwl`** — `split_bam_to_subfiles_SEorPE.pl` partitions reads into 25 `.tmp` files by the first 2 nt of the UMI randomer (AA, AC, ..., NN).
3. **`getpair.cwl`** — JavaScript expression tool that returns matched rep/rmrep `.tmp` file pairs for each UMI prefix.
4. **`deduplicate.cwl`** — `duplicate_removal_inline_paired...pl` per prefix: deduplicates using UMIs, resolves genome vs. repeat-family conflicts, outputs `.rmDup.sam`, `.prermDup.sam`, and `.parsed_v2.20201210.txt`.
5. **`combine.cwl`** — `merge_multiple_parsed_files.simplified_20191022.pl` merges the 25 per-prefix `.parsed` files into one.

`wf_ecliprepmap_pe.cwl` runs two IP barcodes in parallel (via the 1-barcode sub-workflow), then merges them with `combine.cwl` and computes fold change. `wf_ecliprepmap_se.cwl` is structurally identical but with one IP barcode.

## Required Reference Files

Each YAML job file must provide paths to:

| Field | Description |
|---|---|
| `bowtie2_db` | Directory containing Bowtie2 index |
| `bowtie2_prefix` | Index prefix (e.g. `MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat`) |
| `fileListFile1` | TSV mapping ENST/ENSG IDs → repeat families (5 columns) |
| `gencodeGTF` | Gencode annotation GTF |
| `gencodeTableBrowser` | GTF in UCSC table browser format (`.parsed_ucsc_tableformat`) |
| `repMaskBEDFile` | BED file of unique genomic regions (e.g. `UniqueGenomicElements.hg38.bed`) |

Reference data for hg38: [repeat-family-mapping-grch38.tar.gz](https://external-collaborator-data.s3-us-west-1.amazonaws.com/reference-data/repeat-family-mapping-grch38.tar.gz). Reference generation scripts and logic for new assemblies (mm10, mm39) are described in `prompts/generate_refdata.md`.

## Key Implementation Details

### Perl version sensitivity
Perl 5.18+ introduced non-deterministic hash iteration, affecting tie-breaking in deduplication. Validated with perl ≤5.16. TSCC system perl (5.10) is deterministic. Perl 5.18+ is partially mitigated by iterating over sorted hash keys in `duplicate_removal_inline_paired...pl`.

### UMI-based splitting
Deduplication is split into 25 serial jobs (one per UMI 2-nt prefix AA..NN) to cap memory usage. The `prefixes` input in CWL (and `PREFIXES` list in the Snakemake rules) controls this list and drives scatter.

### Conflict resolution (repeat vs. unique genome)
In `duplicate_removal_inline_paired...pl`: a read mapping to both the unique genome and a repeat element is assigned to the unique genome only if the genome alignment score is more than `2 * 2 * 6 = 24` alignment score units better. Otherwise the repeat-element assignment is kept.

### rRNA special handling
In `parse_bowtie2_output_realtime_includemultifamily_PE/SE.pl`: RNA28S, RNA18S, RNA5-8S are treated as members of the RNA45S precursor family (`rRNA_extra_hash`). Reads mapping to any of these count toward RNA45S.

### Multi-family reads
Elements in the second column of `.parsed` files that contain `|` characters are ambiguously mapped to multiple repeat families. `calculate_fold_change_from_parsed_files.py` separates these into `.withpipes.tsv` vs. `.nopipes.tsv`. Multifamily reads should generally be excluded from downstream analysis.

## Output Files

| Extension | Description |
|---|---|
| `.rmDup.sam.gz` | Deduplicated SAM-like file (repeat + unique genome assignments) |
| `.preRmDup.sam.gz` | Pre-deduplication SAM-like file |
| `.parsed` | Per-element read counts and reads-per-million (RPR) with `#READINFO` header lines |
| `.nopipes.tsv` | Unambiguously mapped families with fold enrichment and information content |
| `.withpipes.tsv` | All families including multi-family reads |

## Verification / Testing

Reference output files for validating a run:
- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_pe/` — PE `.nopipes.tsv` and `.withpipes.tsv`
- `test-provenance/tests/ecliprepmap-1.0.0/wf_ecliprepmap_se/` — SE `.nopipes.tsv` and `.withpipes.tsv`

Compare your output `.nopipes.tsv` and `.withpipes.tsv` against these reference files to verify a run.

## Subagent

A specialized bioinformatics pipeline subagent definition lives at `.claude/agents/bioinformatics-pipeline-developer.md`. It has expertise in CWL, Snakemake, Nextflow, and scripting in Perl/Python/R for genomics workflows.


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
