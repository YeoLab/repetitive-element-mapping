# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

This is the **eCLIP repetitive element mapping** pipeline (ecliprepmap). Given trimmed eCLIP FASTQ files and a BAM of reads already mapped to the unique genome (rmRep BAM), it:
1. Maps reads to a curated repeat element database using Bowtie2
2. Splits reads by UMI prefix (AA, AC, ..., NN — 25 bins) for memory-efficient deduplication
3. Deduplicates using randomer UMIs, resolving conflicts between unique genomic and repeat-family mappings
4. Outputs per-element read counts and fold enrichment (IP vs. Input)

There are two parallel implementations: the original **CWL workflow** (production-ready, orchestrated by bash scripts in `wf/`) and a newer **Snakemake "dropin" workflow** (in `workflow/rules/`) that calls the same Perl scripts via Python shim wrappers.

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

### Snakemake Dropin Workflow (TSCC)

The dropin workflow accepts the same YAML inputs as the CWL workflow. Run from the repo root:

```bash
snakemake \
  --snakefile workflow/rules/dropin_repelement.smk \
  --config barcode1r1FastqGz=/path/to/r1.fq.gz \
           barcode1r2FastqGz=/path/to/r2.fq.gz \
           barcode1rmRepBam=/path/to/rmrep.bam \
           barcode1Inputr1FastqGz=/path/to/input_r1.fq.gz \
           barcode1Inputr2FastqGz=/path/to/input_r2.fq.gz \
           barcode1InputrmRepBam=/path/to/input_rmrep.bam \
           bowtie2_db=/path/to/bowtie2_index \
           bowtie2_prefix=MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
           fileListFile1=/path/to/MASTER_FILELIST.tsv \
           gencodeGTF=/path/to/gencode.gtf \
           gencodeTableBrowser=/path/to/gencode.gtf.parsed_ucsc_tableformat.tsv \
           repMaskBEDFile=/path/to/UniqueGenomicElements.hg38.bed \
           se_or_pe=PE

# For SE runs, set se_or_pe=SE and only provide barcode1r1FastqGz/barcode1rmRepBam
```

On TSCC with SLURM, use the profile:

```bash
snakemake \
  --snakefile workflow/rules/dropin_repelement.smk \
  --profile profiles/tscc2_snakemake9 \
  --config ...
```

The SLURM profile (`profiles/tscc2_snakemake9/`) defaults to partition `gold`, account `csd792`, 20 GB memory, 30-minute wall time. The dropin workflow uses `--use-conda` with `workflow/envs/dropin.yaml` (Python 3.11, bowtie2 ≥2.5, samtools ≥1.17, numpy, pandas).

### SE Foundation Tests (mini-dataset)

```bash
snakemake \
  --snakefile workflow/rules/se_foundation.smk \
  --config mini.source_sam_gz=tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz \
           mini.split_manifest=tests/fixtures/mini/expected/split.perl.manifest.tsv \
           mini.merge_input_1=tests/fixtures/mini/source/merge_input_1.parsed_v2.txt \
           mini.merge_input_2=tests/fixtures/mini/source/merge_input_2.parsed_v2.txt
```

This verifies that the Python port of `split_bam_to_subfiles_SEorPE` produces the same split files as the Perl original, and that `merge_multiple_parsed_files.py` matches Perl output.

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

**Final step (both):** `calculate_fold_change_from_parsed_files.cwl` → `bin/calculate_fold_change_from_parsed_files.py` → `.nopipes.tsv` and `.withpipes.tsv`.

### Snakemake Dropin vs. CWL: Script Layer

The Snakemake dropin runs the **same Perl scripts** as the CWL workflow but invokes them through thin Python shims in `bin/python/`:

| `bin/python/` script | Delegates to |
|---|---|
| `parse_bowtie2_output_realtime_includemultifamily_PE.py` | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl` |
| `parse_bowtie2_output_realtime_includemultifamily_SE.py` | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl` |
| `duplicate_removal_inline_paired...py` | `bin/perl/duplicate_removal_inline_paired...pl` |
| `split_bam_to_subfiles_SEorPE.py` | Pure Python re-implementation (does NOT delegate to Perl) |
| `merge_multiple_parsed_files.simplified_20191022.py` | Pure Python re-implementation (does NOT delegate to Perl) |

The shims use `_perl_compat.py` which resolves the Perl script relative to the repo root and calls it with `subprocess`. `split_bam_to_subfiles_SEorPE.py` and `merge_multiple_parsed_files.py` are full Python ports that the `se_foundation.smk` tests verify against Perl output.

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
