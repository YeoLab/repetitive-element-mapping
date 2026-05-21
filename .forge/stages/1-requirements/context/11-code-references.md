# Code References

## Key Files for the Implementer

### Ground Truth Reference Files (hg38, read-only)

| File | Path | Purpose |
|------|------|---------|
| parsed_ucsc_tableformat | `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` | Column schema and row format ground truth (249,044 lines) |
| MASTER_FILELIST TSV | `examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` | 5-column TSV ground truth (26,422 lines) |
| MASTER_FILELIST list | `examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list` | Identical content to TSV; .list extension is what Perl reads |
| UniqueGenomicElements | `examples/inputs/hg38/UniqueGenomicElements.hg38.bed` | 6-column BED ground truth (5,618,483 lines) |
| bowtie2_index FASTA | `examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa` | Combined FASTA showing sequence header naming conventions |
| bowtie2_index dir | `examples/inputs/hg38/bowtie2_index/` | Contains 7 files: 6 .bt2 index + 1 .fa |

### Source Input Files (per assembly)

**hg38:**
- `examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf` (symlink)
- `examples/inputs/hg38/downloaded/hg38.fasta` (symlink to GRCh38_no_alt... fasta)
- `examples/inputs/hg38/downloaded/hg38.repeatmasker.tsv.gz`
- `examples/inputs/hg38/downloaded/hg38.simplerepeats.tsv.gz`
- `examples/inputs/hg38/downloaded/hg38.trna.tsv.gz`
- `examples/inputs/hg38/downloaded/hsa.gff3`

**mm10:**
- `examples/inputs/mm10/downloaded/gencode.VM23.annotation.gtf.gz`
- `examples/inputs/mm10/downloaded/mm10.repeatmasker.tsv.gz`
- `examples/inputs/mm10/downloaded/mm10.simplerepeats.tsv.gz`
- `examples/inputs/mm10/downloaded/mm10.trna.tsv.gz`
- `examples/inputs/mm10/downloaded/mmu.gff3`
- `examples/inputs/mm10/downloaded/NR_046233.2.fasta`

**mm39:**
- `examples/inputs/mm39/downloaded/gencode.VM38.annotation.gtf.gz`
- `examples/inputs/mm39/downloaded/mm39.repeatmasker.tsv.gz`
- `examples/inputs/mm39/downloaded/mm39.simplerepeats.tsv.gz`
- `examples/inputs/mm39/downloaded/NR_046233.2.fasta`
- **Missing:** mm39.trna.tsv.gz, mm39.gff3 (not available; scripts must handle absence gracefully)
- **Missing:** mm39 genome FASTA (must be sourced externally)

### Perl Scripts (read; modify only if hardcoded values block compatibility)

| Script | Path | Notes |
|--------|------|-------|
| parse_bowtie2_SE | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl` | Reads MASTER_FILELIST via ARGV[3]; hardcoded paths are commented out |
| parse_bowtie2_PE | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl` | Same pattern as SE |
| deduplicate | `bin/perl/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl` | Check for assembly assumptions |
| split_bam | `bin/perl/split_bam_to_subfiles_SEorPE.pl` | Comment mentions hg38 but logic is assembly-agnostic |
| RepElement_pipeline | `bin/perl/RepElement_pipeline_1dataset.pl` | Has `my $species = "hg38"` hardcoded on line 4; NOT called by CWL or Snakemake dropin |

### Python Shims (read-only context)

| Script | Path | Notes |
|--------|------|-------|
| perl_compat | `bin/python/_perl_compat.py` | Resolves Perl script paths relative to repo root |
| split_bam py | `bin/python/split_bam_to_subfiles_SEorPE.py` | Pure Python port verified against Perl output |
| merge_parsed py | `bin/python/merge_multiple_parsed_files.simplified_20191022.py` | Pure Python port |

### CWL Workflow Files (read-only for compatibility verification)

| File | Path |
|------|------|
| SE top-level | `cwl/wf_ecliprepmap_se.cwl` |
| PE top-level | `cwl/wf_ecliprepmap_pe.cwl` |
| Map step | `cwl/map_repetitive_elements_se.cwl`, `cwl/map_repetitive_elements_pe.cwl` |
| Deduplicate step | `cwl/deduplicate.cwl` |

### Snakemake Files (read-only for compatibility verification)

| File | Path |
|------|------|
| Dropin workflow | `workflow/rules/dropin_repelement.smk` |
| Conda env | `workflow/envs/dropin.yaml` |
| SLURM profile | `profiles/tscc2_snakemake9/` |

### Task Specification

| File | Path |
|------|------|
| Primary task prompt | `prompts/generate_refdata.md` |
| Project context | `CLAUDE.md` |
| Knowledge graph | `.forge/stages/0-research/graphify-initial/GRAPH_REPORT.md` |

### New Scripts to Create (by implementer)

All new scripts must be placed in a logical location (e.g., `bin/python/refdata_generation/` or `bin/python/`) and accept assembly-generic CLI arguments:

1. `generate_parsed_ucsc_tableformat.py` — Step 1
2. `generate_bowtie2_index.py` — Step 2
3. `generate_unique_genomic_elements.py` — Step 3
4. `generate_master_filelist.py` — Step 4

### Verification Commands

```bash
# Dry-run Snakemake with SE mm10 references
snakemake -s workflow/rules/dropin_repelement.smk \
  --config se_or_pe=SE \
    barcode1r1FastqGz=<path> \
    barcode1rmRepBam=<path> \
    barcode1Inputr1FastqGz=<path> \
    barcode1InputrmRepBam=<path> \
    bowtie2_db=examples/inputs/mm10/bowtie2_index \
    bowtie2_prefix=MASTER_FILELIST.<date>.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
    fileListFile1=examples/inputs/mm10/MASTER_FILELIST.<date>.*.tsv \
    gencodeGTF=<mm10_gtf> \
    gencodeTableBrowser=examples/inputs/mm10/gencode.VM23.*.parsed_ucsc_tableformat \
    repMaskBEDFile=examples/inputs/mm10/UniqueGenomicElements.mm10.bed \
  -n

# Line count comparison (example)
wc -l examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
# Expected: 249044

wc -l examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv
# Expected: 26422

wc -l examples/inputs/hg38/UniqueGenomicElements.hg38.bed
# Expected: 5618483
```
