# Repetitive Element Mapping Pipeline - Structural Codemap

**Project**: eCLIP Repetitive Element Mapping  
**Language**: Perl (primary scripts), Bash (launchers), CWL (workflows), Python (analysis), Snakemake (new dropin)  
**Architecture**: Parallel workflow orchestration with dual CWL + Snakemake implementations  
**Files Indexed**: 40+ source files | **Est. Tokens**: ~45K vs ~500K reading all | **Compression**: ~10x

---

## 1. Project Overview

This pipeline maps eCLIP sequencing reads to repetitive element databases and performs deduplication with UMI-based conflict resolution. Two parallel implementations exist:
- **CWL Workflow** (production, tested): Bash launcher scripts + CWL job descriptions
- **Snakemake Dropin** (newer): Snakemake rules calling same Perl scripts

### Core Data Flow
```
FASTQ files + rmRep BAM → Bowtie2 mapping (25 UMI bins) 
→ Deduplication + QC → Per-element read counts → Fold enrichment (IP/Input)
```

---

## 2. Directory Structure & Key Files

### Root Configuration
```
/
├── CLAUDE.md                          # Project guidance + running instructions
├── README.md                          # Methods, requirements, output format
├── config/                            # Runtime configuration stubs
├── examples/                          # YAML example inputs for both SE/PE
│   ├── repeat_mapping_SE.yaml         # Single-end example config
│   └── repeat_mapping_PE.yaml         # Paired-end example config
└── profiles/                          # Snakemake profiles (TSCC SLURM)
    └── tscc2_snakemake9/              # SLURM profile config
```

### CWL Workflow (Production)
```
wf/                                    # Bash launcher entry points
├── eCLIP_repelement_SE                # SE launcher → cwl/wf_ecliprepmap_se.cwl
├── eCLIP_repelement_PE                # PE launcher → cwl/wf_ecliprepmap_pe.cwl
├── eCLIP_repelement_SE_singleNode     # Single barcode variant
└── eCLIP_repelement_PE_singleNode     # Single barcode variant

cwl/                                   # CWL job definitions (Docker containers)
├── wf_ecliprepmap_se.cwl              # Top-level SE workflow
├── wf_ecliprepmap_pe.cwl              # Top-level PE workflow
├── wf_ecliprepmap_se_1barcode.cwl     # Single-barcode SE variant
├── wf_ecliprepmap_pe_1barcode.cwl     # Single-barcode PE variant
├── map_repetitive_elements_se.cwl     # Bowtie2 mapping (SE)
├── map_repetitive_elements_pe.cwl     # Bowtie2 mapping (PE)
├── deduplicate.cwl                    # UMI deduplication & conflict resolution
├── splitbam.cwl                       # Split BAM by UMI prefix (AA..NN)
├── combine.cwl                        # Combine deduplicated results
├── concatenate.cwl                    # Concatenate final outputs
├── calculate_fold_change_from_parsed_files.cwl  # IP/Input fold enrichment
├── getpair.cwl                        # Utility: extract paired reads
├── gzip.cwl                           # Utility: compress outputs
└── [others]                           # Helper/deprecated CWLs
```

### Snakemake Workflow (New Dropin)
```
workflow/                              # Snakemake-based implementation
├── rules/
│   ├── dropin_repelement.smk          # Main dropin workflow
│   └── se_foundation.smk              # Validation workflow for SE pipeline
└── envs/
    └── dropin.yaml                    # Conda environment (Python 3.11, Bowtie2, samtools)
```

### Processing Scripts
```
bin/
├── calculate_fold_change_from_parsed_files.py    # Compute fold enrichment IP/Input
└── perl/                              # Core Perl scripts (invoked by CWL)
    ├── parse_bowtie2_output_realtime_includemultifamily_SE.pl  (481 lines)
    │   └─ Parses Bowtie2 SAM; handles multi-family/multi-element mapping
    ├── parse_bowtie2_output_realtime_includemultifamily_PE.pl  (520 lines)
    │   └─ PE variant with pair validation logic
    ├── duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl  (1795 lines)
    │   └─ Main deduplication; UMI conflict resolution, genomic vs. repeat priority
    ├── duplicate_removal.pl            # Symlink to above
    ├── split_bam_to_subfiles_SEorPE.pl   (138 lines)
    │   └─ Split SAM/BAM by first 2 nucleotides of UMI (25 bins: AA..NN)
    ├── merge_multiple_parsed_files.simplified_20191022.pl  (116 lines)
    │   └─ Merge parsed results from 25 bins; compute per-element RPM
    ├── RepElement_pipeline_1dataset.pl  (173 lines)
    │   └─ Single-dataset orchestration (legacy/reference)
    └── ._reparse_samfile_updatedchrM_fixmultenstsort_SE.pl  (archived)

test-provenance/                       # Test data & expected outputs (git-tracked)
```

### Documentation & Examples
```
documentation/                         # Protocol & reference docs
examples/                              # YAML configurations
    ├── repeat_mapping_SE.yaml         # SE: 1 IP barcode + 1 Input
    └── repeat_mapping_PE.yaml         # PE: 2 IP barcodes + 1 Input
archived_workflow_definitions/         # Pre-CWL shell scripts + older variants
```

---

## 3. Data Processing Pipeline (Logical Flow)

### Input
- `barcode1r1FastqGz`, `barcode1r2FastqGz`: Trimmed eCLIP FASTQ (gzipped)
- `barcode1rmRepBam`: Unique genomic mapping BAM (reads mapped to hg38 minus repeats)
- `bowtie2_db` + `bowtie2_prefix`: Repeat element Bowtie2 index
- `fileListFile1`: Master repeat annotation (priority ordering)

### Processing Stages

#### 1. **Bowtie2 Mapping** → `parse_bowtie2_output_realtime_includemultifamily_*.pl`
- Command: `bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder -x $db -1 $r1 -2 $r2`
- Outputs: SAM-like file with multi-family annotated reads
- Logic:
  - Multi-element mappings per family → select best (mismatch score + quality)
  - Multi-family mappings → flagged but retained (excluded downstream)

#### 2. **UMI-Based Splitting** → `split_bam_to_subfiles_SEorPE.pl`
- Splits parsed/BAM by first 2 nucleotides of randomer UMI
- 25 bins: AA, AC, AG, AT, AN, CA, ..., NN
- Purpose: Process serially to avoid memory explosion

#### 3. **Deduplication** → `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
- For each UMI bin (AA..NN):
  - Parse repeat & genomic mappings
  - Conflict resolution:
    - If genomic alignment > 2 mismatches better (12 alignment score): use genomic
    - Else: use repeat element
  - rRNA handling: 18S/28S/5-8S merged into 45S family
  - Multi-element within family: keep shorter transcript
  - Output: deduplicated .parsed file

#### 4. **Result Merging** → `merge_multiple_parsed_files.simplified_20191022.pl`
- Concatenate 25 .parsed outputs
- Compute per-element RPM (reads per million usable)
- Output: combined .parsed + .reparsed + .nopipes.tsv + .withpipes.tsv

#### 5. **Fold Enrichment** → `calculate_fold_change_from_parsed_files.py`
- Input: IP .parsed + Input .parsed
- Compute: log2(fold change), information content
- Output: .withpipes.tsv (all elements) + .nopipes.tsv (unambiguous only)

---

## 4. CWL Workflow Graph

### SE Workflow: `cwl/wf_ecliprepmap_se.cwl`
```
Inputs: dataset, barcode1r1FastqGz, barcode1Inputr1FastqGz, bowtie2_db, fileListFile1, etc.

map_ip_barcode1 (map_repetitive_elements_se.cwl)
  ├─ Bowtie2 + parse_bowtie2_output_*_SE.pl
  └─ Outputs: .parsed (SAM-like)

map_input (map_repetitive_elements_se.cwl)
  └─ Same as above for Input control

splitbam_ip (splitbam.cwl) ─┐
splitbam_input (splitbam.cwl) ├─ Splits .parsed by UMI (AA..NN)
                             └─ Outputs: 25 files each

deduplicate (deduplicate.cwl)
  ├─ Loop over 25 UMI bins
  ├─ Calls duplicate_removal_*.pl for each
  └─ Outputs: 25 .parsed files per sample

combine (combine.cwl)
  ├─ Merges 25 deduplicated .parsed files (IP & Input separately)
  └─ Outputs: .parsed + .reparsed + .nopipes.tsv + .withpipes.tsv

calculate_fold_change (calculate_fold_change_from_parsed_files.cwl)
  ├─ Inputs: IP .nopipes.tsv + Input .nopipes.tsv
  └─ Outputs: .nopipes.tsv with fold change, information content

Final outputs in: <job_name>/results/
```

### PE Workflow: `cwl/wf_ecliprepmap_pe.cwl`
- Same as SE but with 2 IP barcodes
- Processes each barcode in parallel then merges all 3 samples (2 IP + 1 Input)

---

## 5. Snakemake Dropin Workflow

### Entry Point: `workflow/rules/dropin_repelement.smk`
- Accepts same YAML config as CWL
- Calls Perl scripts via Python wrapper rules
- Uses `workflow/envs/dropin.yaml`: Python 3.11, bowtie2 ≥2.5, samtools ≥1.17, numpy, pandas
- Profile: `profiles/tscc2_snakemake9/` (SLURM, gold partition, 20 GB, 30 min)

### Validation Workflow: `workflow/rules/se_foundation.smk`
- Tests that Python ports of `split_bam_to_subfiles_SEorPE` and `merge_multiple_parsed_files` match Perl
- Mini test data: `tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz`
- Validates output equivalence with Perl reference

---

## 6. Key Script Signatures

### Perl: `parse_bowtie2_output_realtime_includemultifamily_SE.pl`
```perl
Usage: parse_bowtie2_output_realtime_includemultifamily_SE.pl <SAM_file> <fileListFile> <output_file>

Parses Bowtie2 SAM output:
  - Reads: QNAME, bitflag, RNAME (repeat), MAPQ, CIGAR, TLEN, AS, NM
  - Annotations: RepBase family, transcript, gene
  - Multi-mapping: resolves by (mismatch score, quality, input ordering)
  - rRNA special: 18S/28S/5-8S → RNA45S family

Output: SAM-like with family|transcript|gene annotations
```

### Perl: `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
```perl
Usage: duplicate_removal_*.pl \
  <repeat_mapped_SAM> \
  <unique_genomic_BAM> \
  <fileListFile> \
  <output_parsed> \
  <output_reparsed> \
  <se_or_pe_flag> \
  <prefixes_string> \
  ...

Main deduplication logic:
  1. Load UMI (randomer) from read name
  2. Compare repeat vs. genomic alignment quality
  3. If genomic is > 12 score better (2 mismatches × 2 reads × 6): use genomic
  4. Else: use repeat element
  5. Multi-family within same element: keep shorter transcript
  6. Compute RPM per family and element

Output: .parsed (4 READINFO lines + TOTAL/ELEMENT rows)
```

### Perl: `split_bam_to_subfiles_SEorPE.pl`
```perl
Usage: split_bam_to_subfiles_SEorPE.pl <input_SAM> <output_dir> <se_or_pe>

Splits by first 2 nucleotides of UMI:
  - UMI extracted from QNAME (format: sample:barcode:umi:other)
  - 25 output files: AA.sam, AC.sam, ..., NN.sam
  - Reduces memory footprint for deduplication

Output: 25 split SAM files per input
```

### Perl: `merge_multiple_parsed_files.simplified_20191022.pl`
```perl
Usage: merge_multiple_parsed_files.simplified_20191022.pl <dir_with_25_parsed_files> <output_file>

Merges 25 .parsed files from deduplication:
  1. Aggregate READINFO totals
  2. Sum TOTAL and ELEMENT read counts
  3. Compute per-element RPM
  4. Output: .parsed (Perl format) or .reparsed (without AllReads info)

Output: Combined .parsed, .nopipes.tsv, .withpipes.tsv
```

### Python: `calculate_fold_change_from_parsed_files.py`
```python
Usage: python3 calculate_fold_change_from_parsed_files.py \
  --ip_parsed <IP.parsed> \
  --input_parsed <Input.parsed> \
  --output_nopipes <output.nopipes.tsv> \
  --output_withpipes <output.withpipes.tsv>

Computes enrichment:
  - Reads .parsed using pandas
  - For each element:
    - IP_rpr = IP_reads / IP_total_usable_reads
    - Input_rpr = Input_reads / Input_total_usable_reads (+ 1 pseudocount)
    - Fold_enrichment = IP_rpr / Input_rpr
    - Info_content = IP_rpr * log2(fold_enrichment)
  
Output: .nopipes.tsv (unambiguous) + .withpipes.tsv (all)
Columns: element, IP_read_num, IP_clip_rpr, Input_read_num, Input_clip_rpr, Fold_enrichment, Information_content
```

---

## 7. Configuration & Inputs

### Example YAML: `examples/repeat_mapping_SE.yaml`
```yaml
dataset: "eCLIP_SE_example"

barcode1r1FastqGz: /path/to/barcode1.r1.fq.gz
barcode1rmRepBam: /path/to/barcode1.rmRep.bam

barcode1Inputr1FastqGz: /path/to/barcode1_input.r1.fq.gz
barcode1InputrmRepBam: /path/to/barcode1_input.rmRep.bam

bowtie2_db: /path/to/bowtie2_index/
bowtie2_prefix: MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat

fileListFile1: /path/to/MASTER_FILELIST.tsv
gencodeGTF: /path/to/gencode.gtf
gencodeTableBrowser: /path/to/gencode.gtf.parsed_ucsc_tableformat.tsv
repMaskBEDFile: /path/to/UniqueGenomicElements.hg38.bed

se_or_pe: SE
```

### Example YAML: `examples/repeat_mapping_PE.yaml`
- Same as above but with `barcode2r1FastqGz`, `barcode2r2FastqGz`, `barcode2rmRepBam` for 2nd IP barcode
- `se_or_pe: PE`

### Bowtie2 Command Template
```bash
bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder \
  -x <bowtie2_db>/<bowtie2_prefix> \
  -1 <read1.fq.gz> -2 <read2.fq.gz> \
  2> <bowtie_stderr>
```

---

## 8. Output Format

### .parsed File (4-6 columns)
```
#READINFO	AllReads	<total_reads>
#READINFO	UsableReads	<dedup_reads>
#READINFO	GenomicReads	<unique_genome_reads>
#READINFO	RepFamilyReads	<repeat_element_reads>
TOTAL	ELEMENT_NAME	read_count	RPM	annotation	gene_names
ELEMENT	ELEMENT_NAME	read_count	RPM	family|transcript|...	gene_names
ELEMENT	ELEMENT_NAME	read_count	RPM	family1|family2|...	gene_names  # multi-family (pipes)
```

### .nopipes.tsv (No multi-family)
```
element	IP_read_num	IP_clip_rpr	Input_read_num	Input_clip_rpr	Fold_enrichment	Information_content
```

### .withpipes.tsv (All, including multi-family)
```
Same columns as .nopipes.tsv
```

---

## 9. Running the Pipeline

### CWL (Yeo Lab)
```bash
module load ecliprepmap
eCLIP_repelement_SE examples/repeat_mapping_SE.yaml
# Results: repeat_mapping_SE/results/
# Log: repeat_mapping_SE/REPELEMENTMAPPING_repeat_mapping_SE_LOG.txt
```

### CWL (External)
```bash
cwltool --outdir <outdir> --cachedir <cachedir> \
  cwl/wf_ecliprepmap_se.cwl examples/repeat_mapping_SE.yaml
```

### Snakemake (TSCC)
```bash
snakemake \
  --snakefile workflow/rules/dropin_repelement.smk \
  --profile profiles/tscc2_snakemake9 \
  --config barcode1r1FastqGz=/path/to/r1.fq.gz \
           barcode1rmRepBam=/path/to/rmrep.bam \
           ... [other params from YAML]
```

### Snakemake Validation
```bash
snakemake \
  --snakefile workflow/rules/se_foundation.smk \
  --config mini.source_sam_gz=tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz \
           mini.split_manifest=tests/fixtures/mini/expected/split.perl.manifest.tsv \
           mini.merge_input_1=tests/fixtures/mini/source/merge_input_1.parsed_v2.txt \
           mini.merge_input_2=tests/fixtures/mini/source/merge_input_2.parsed_v2.txt
```

---

## 10. External Tools & Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| Bowtie2 | ≥2.2.6 | Sequence alignment to repeat DB |
| Perl | 5.10.1+ (5.18+ has hash randomness warnings) | Main processing scripts |
| Python | 3.6+ (tested 3.11) | Fold enrichment, Snakemake rules |
| Samtools | ≥1.17 | BAM processing (Snakemake) |
| pandas/numpy | 1.1.3+/1.18+ | Data manipulation (Python) |
| CWL Runner | cwltool 1.0.20180306140409+ | CWL orchestration |
| Snakemake | (profile-specified) | Workflow orchestration (dropin) |
| Docker | (multi-stage) | brianyee/repetitive_element_mapping:1.0.0 |

---

## 11. Key Design Decisions

1. **UMI-based splitting** (25 bins): Mitigates memory pressure during deduplication by processing serially
2. **Conflict resolution priority**: Genomic mapping preferred if significantly better (>2 mismatches), else repeat family
3. **rRNA consolidation**: 18S/28S/5-8S → RNA45S to handle polycistronic precursor
4. **Multi-family flagging**: Reads mapping to multiple families retained but flagged as "pipes" (excluded from no-pipes output)
5. **Perl version variance**: Hash iteration randomness in Perl 5.18+ mitigated by sorting keys (v0.0.2+)
6. **Dual implementations**: CWL (stable, tested) + Snakemake dropin (modern, same Perl scripts)

---

## 12. File Statistics

| Category | Count | Total Lines | Language |
|----------|-------|-------------|----------|
| CWL Workflows | 14 | ~1000 | YAML |
| Bash Launchers | 4 | ~500 | Bash |
| Core Perl Scripts | 6 | ~4000 | Perl |
| Python Scripts | 1 | ~200 | Python |
| Snakemake Rules | 2+ | TBD | Snakemake |
| YAML Configs | 2 | ~40 | YAML |

**Total Active Source**: ~40 files | **Estimated Tokens**: 45-50K | **Compression vs. reading all**: 10x

---

## 13. Quick Navigation

- **To understand the pipeline**: Read `CLAUDE.md` → `README.md` → trace SE/PE workflow in `cwl/wf_ecliprepmap_se.cwl`
- **To run locally**: Use `cwltool` with `examples/repeat_mapping_SE.yaml` (requires Bowtie2, docker/singularity)
- **To run on TSCC**: `module load ecliprepmap` + `eCLIP_repelement_SE examples/repeat_mapping_SE.yaml`
- **To modify Perl logic**: Edit `bin/perl/parse_bowtie2_output_*.pl` or `duplicate_removal_*.pl`, rebuild CWL image
- **To test Snakemake dropin**: Run `snakemake ... workflow/rules/se_foundation.smk` with test fixtures
- **To understand deduplication**: Read `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl` (1795 lines, core algorithm)

---

Generated: 2026-05-14
Codemap Version: 1.0
