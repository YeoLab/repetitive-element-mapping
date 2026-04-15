# Test Documentation for Repetitive Element Mapping Pipeline

**Generated:** 2026-04-15
**Purpose:** Complete documentation of all tests in `tests/` and `examples/` directories to enable reproduction after deletion

---

## Overview

This document provides comprehensive documentation of all tests performed in the `tests/` and `examples/` directories. The total size of these directories is approximately **147GB**, consisting of:
- `tests/`: 130GB
- `examples/`: 17GB

---

## Table of Contents

1. [Python Unit Tests](#python-unit-tests)
2. [Test Fixtures](#test-fixtures)
3. [CWL/ecliprepmap-1.0.0 Integration Tests](#cwlecliprepmap-100-integration-tests)
4. [Example Data and Workflows](#example-data-and-workflows)
5. [Reproduction Instructions](#reproduction-instructions)
6. [File Manifests](#file-manifests)

---

## Python Unit Tests

### Location
`tests/*.py`

### Test Files

#### 1. `test_split_merge_python_expected.py`
**Purpose:** Validate Python implementations of BAM splitting and parsed file merging
**Tests:**
- `test_split_matches_expected_manifest()`: Verifies that `split_bam_to_subfiles_SEorPE.py` produces exactly the expected output files with correct sizes and SHA256 hashes
- `test_merge_matches_expected_output()`: Verifies that `merge_multiple_parsed_files.simplified_20191022.py` produces deterministic output

**Dependencies:**
- `bin/python/split_bam_to_subfiles_SEorPE.py`
- `bin/python/merge_multiple_parsed_files.simplified_20191022.py`
- `tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz`
- `tests/fixtures/mini/source/merge_input_1.parsed_v2.txt`
- `tests/fixtures/mini/source/merge_input_2.parsed_v2.txt`
- `tests/fixtures/mini/expected/split.expected.manifest.tsv`
- `tests/fixtures/mini/expected/merged.expected.parsed`

**Execution:**
```bash
pytest tests/test_split_merge_python_expected.py -v
```

---

#### 2. `test_mini_fixtures_manifest.py`
**Purpose:** Ensure integrity of mini test fixtures
**Tests:**
- `test_mini_fixture_manifest_matches_files()`: Validates that all files listed in `tests/fixtures/mini/manifest.tsv` exist with correct sizes and SHA256 checksums

**Dependencies:**
- `tests/fixtures/mini/manifest.tsv`

**Execution:**
```bash
pytest tests/test_mini_fixtures_manifest.py -v
```

---

#### 3. `test_cwl_yaml_adapter.py`
**Purpose:** Test configuration adapters that convert CWL-style and simplified YAML inputs to internal config format
**Tests:**
- `test_map_cwl_job_to_config_normalizes_path_objects()`: Validates path normalization for CWL File/Directory objects
- `test_merge_cwl_job_yaml_into_config_uses_fixture()`: Tests merging CWL YAML with base config
- `test_map_simple_job_to_config_se()`: Tests simplified SE mode YAML parsing
- `test_merge_job_yaml_into_config_with_simple_fixture()`: Tests simplified YAML merging

**Dependencies:**
- `workflow/config_adapter.py`
- `tests/fixtures/mini/cwl_job_example.yaml`
- `tests/fixtures/dropin/simple_se_job.yaml`

**Execution:**
```bash
pytest tests/test_cwl_yaml_adapter.py -v
```

---

#### 4. `test_parse_pe_python.py`
**Purpose:** Test paired-end bowtie2 parser
**Tests:**
- `test_parse_pe_runs_on_sam_input()`: Validates that the PE parser correctly processes SAM input and produces expected output files

**Dependencies:**
- `bin/python/parse_bowtie2_output_realtime_includemultifamily_PE.py`

**Execution:**
```bash
pytest tests/test_parse_pe_python.py -v
```

---

#### 5. `test_parse_se_python.py`
**Purpose:** Test single-end bowtie2 parser (inferred from naming pattern)

**Dependencies:**
- `bin/python/parse_bowtie2_output_realtime_includemultifamily_SE.py` (inferred)

**Execution:**
```bash
pytest tests/test_parse_se_python.py -v
```

---

## Test Fixtures

### Minimal Test Fixtures (`tests/fixtures/mini/`)

**Purpose:** Lightweight fixtures for fast unit testing

**Structure:**
```
tests/fixtures/mini/
├── manifest.tsv                                    # File integrity manifest
├── cwl_job_example.yaml                           # CWL-style job configuration
├── source/
│   ├── ip.preRmDup.sam.mini.gz                   # IP sample SAM (gzipped)
│   ├── input.preRmDup.sam.mini.gz                # Input control SAM (gzipped)
│   ├── merge_input_1.parsed_v2.txt               # Parsed file for merge testing
│   └── merge_input_2.parsed_v2.txt               # Parsed file for merge testing
├── expected/
│   ├── split.expected.manifest.tsv               # Expected split output manifest
│   ├── merged.expected.parsed                    # Expected merged output
│   ├── ip.parsed.mini.gz                         # Expected IP parsed output
│   ├── input.parsed.mini.gz                      # Expected input parsed output
│   ├── ip.rmDup.sam.mini.gz                      # Expected IP deduplicated SAM
│   ├── input.rmDup.sam.mini.gz                   # Expected input deduplicated SAM
│   └── ip.reparsed.nopipes.mini.tsv.gz          # Expected reparsed output
└── refs/
    ├── bowtie2_index/                            # Mini bowtie2 index
    ├── filelist.tsv                              # Repeat element annotation
    ├── gencode.gtf                               # Gene annotations
    ├── gencode.table.tsv                         # Gene table
    └── repmask.bed                               # RepeatMasker annotations
```

**Configuration:**
- Dataset: `fixture_dataset`
- Mode: SE
- Prefixes tested: `[AA, AC, TT]`

---

### Drop-in Test Fixtures (`tests/fixtures/dropin/`)

**Purpose:** Simplified YAML format testing

**Files:**
- `se_job.yaml`: Mock SE job with simplified format
- `simple_se_job.yaml`: Simple SE job configuration
- `pe_job.yaml`: Mock PE job configuration
- `simple_pe_job.yaml`: Simple PE job configuration

**Example Configuration (`se_job.yaml`):**
```yaml
dataset: mock_se
barcode1r1FastqGz:
  class: File
  path: mock.fastq.gz
barcode1rmRepBam:
  class: File
  path: mock_rmrep.sam
# ... (abbreviated)
prefixes: [AA]
se_or_pe: SE
```

---

## CWL/ecliprepmap-1.0.0 Integration Tests

### Location
`tests/ecliprepmap-1.0.0/`

**Purpose:** Integration tests for individual CWL tools and complete workflows from the ecliprepmap 1.0.0 pipeline

### Individual Tool Tests

#### SE (Single-End) Pipeline Steps

##### 1. `00_map_repetitive_elements_se/`
**Tool:** Map repetitive elements (SE)
**Config:** `map_repetitive_elements_se.yaml`
**Validates:** Initial mapping of reads to repetitive elements database

##### 2. `01_split_rep_bam_se/`
**Tool:** Split repeat-mapped BAM (SE)
**Config:** `splitbam.yaml`
**Validates:** Splitting repeat-mapped reads by family

##### 3. `02_split_rmrep_bam_se/`
**Tool:** Split genome-mapped (non-repeat) BAM (SE)
**Config:** `splitbam.yaml`
**Validates:** Splitting genome-mapped reads by prefix

##### 4. `03_dedup_se/`
**Tool:** Deduplication (SE)
**Config:** `deduplicate_se.yaml`
**Script:** `run_script.sh`
**Validates:** PCR duplicate removal for SE data

##### 5. `04_merge_parsed/`
**Tool:** Merge parsed files
**Config:** `combine.yaml`
**Validates:** Combining parsed outputs from multiple prefixes

---

#### PE (Paired-End) Pipeline Steps

##### 6. `05_map_repetitive_elements_pe/`
**Tool:** Map repetitive elements (PE)
**Config:** `map_repetitive_elements_pe.yaml`
**Validates:** Initial mapping of PE reads to repetitive elements database

##### 7. `06_split_rep_sam_pe/`
**Tool:** Split repeat-mapped SAM (PE)
**Config:** `splitbam.yaml`
**Validates:** Splitting repeat-mapped PE reads by family

##### 8. `07_split_rmrep_bam_pe/`
**Tool:** Split genome-mapped (non-repeat) BAM (PE)
**Config:** `splitbam.yaml`
**Validates:** Splitting genome-mapped PE reads by prefix

##### 9. `08_dedup_pe/`
**Tool:** Deduplication (PE)
**Config:** `deduplicate_pe.yaml`
**Validates:** PCR duplicate removal for PE data

---

### Complete Workflow Tests

#### SE Workflow: `wf_ecliprepmap_se/`
**Config:** `wf_ecliprepmap_se.yaml`
**Job Input:** `wf_ecliprepmap_se/REPELEMENTMAPPING_wf_ecliprepmap_se_INPUT.yaml`
**Validation Script:** `run_diff.sh`
**Purpose:** End-to-end test of SE repetitive element mapping workflow

#### PE Workflow: `wf_ecliprepmap_pe/`
**Config:** `wf_ecliprepmap_pe.yaml`
**Job Input:** `wf_ecliprepmap_pe/REPELEMENTMAPPING_wf_ecliprepmap_pe_INPUT.yaml`
**Purpose:** End-to-end test of PE repetitive element mapping workflow

---

### Legacy Outputs

#### `_eric_outputs/`
**Purpose:** Reference outputs from original developer (Eric Van Nostrand) for regression testing
**Location:** `_eric_outputs/20201213_PE_FINAL/`
**Size:** ~162KB directory

#### `_eric_scripts/`
**Purpose:** Original CWL/Bash scripts for comparison and provenance
**Size:** ~3.5KB directory

---

## Example Data and Workflows

### Location
`examples/`

### Reference Data Archives

#### 1. `data_for_repeat_mapping_hg19.tar.gz` (4.8GB)
**Genome:** hg19
**Contents:** Test data for hg19-based repeat mapping

#### 2. `repeat-mapping-hg19-refdata.tar.gz` (102MB)
**Genome:** hg19
**Contents:** Reference databases for hg19
- Bowtie2 index for repetitive elements
- MASTER_FILELIST annotations
- Gencode GTF/table
- RepeatMasker BED

#### 3. `example_data_for_repeat_mapping_hg38.tar.gz` (5.2GB)
**Genome:** hg38/GRCh38
**Contents:** Test data for hg38-based repeat mapping

#### 4. `example_data_for_repeat_mapping_hg38/` (extracted, 6.5KB metadata)
**Files:**
- `EXAMPLE_PE.rep1_clip.A01.r1.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_clip.A01.r2.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_clip.A01.r1.fq.genome-mapped.bam`
- `EXAMPLE_PE.rep1_clip.B06.r1.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_clip.B06.r2.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_clip.B06.r1.fq.genome-mapped.bam`
- `EXAMPLE_PE.rep1_input.NIL.r1.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_input.NIL.r2.fqTrTr.sorted.fq.gz`
- `EXAMPLE_PE.rep1_input.NIL.r1.fq.genome-mapped.bam`

#### 5. `repeat-mapping-grch38-refdata/` (3KB metadata)
**Contents:**
- `bowtie2_index/` - Bowtie2 indices for GRCh38
- MASTER_FILELIST annotations
- Gencode v33 GTF and table format
- UniqueGenomicElements.hg38.bed

---

### Example Workflows

#### SE Workflow: `repeat_mapping_SE.yaml`
**Dataset:** `seCLIP_example`
**Mode:** SE
**Data Source:** `/tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/`
**Samples:**
- IP: `INV_B.IP.umi.r1.fqTrTr.sorted.fq.gz` + `INV_B.IP.umi.r1.fq.genome-mappedSoSo.bam`
- Input: `INV_B.IN.umi.r1.fqTrTr.sorted.fq.gz` + `INV_B.IN.umi.r1.fq.genome-mappedSoSo.bam`

**References:**
- Bowtie2 DB: `/tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/bowtie2_index`
- Prefix: `MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat`
- FileList: `MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list`
- Gencode: `/tscc/projects/ps-yeolab4/genomes/hg38/gencode/v33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf`
- RepMask: `UniqueGenomicElements.hg38.bed`

#### SE Simplified: `repeat_mapping_SE.simple.yaml`
**Dataset:** `seCLIP_example`
**Mode:** SE
Simplified format version of the above

---

#### PE Workflow: `repeat_mapping_PE.yaml`
**Dataset:** `peCLIP_example`
**Mode:** PE
**Data Source:** `/home/centos/example_data_for_repeat_mapping_hg38/`
**Samples:**
- Barcode 1 (A01):
  - R1: `EXAMPLE_PE.rep1_clip.A01.r1.fqTrTr.sorted.fq.gz`
  - R2: `EXAMPLE_PE.rep1_clip.A01.r2.fqTrTr.sorted.fq.gz`
  - rmRep BAM: `EXAMPLE_PE.rep1_clip.A01.r1.fq.genome-mapped.bam`
- Barcode 2 (B06):
  - R1: `EXAMPLE_PE.rep1_clip.B06.r1.fqTrTr.sorted.fq.gz`
  - R2: `EXAMPLE_PE.rep1_clip.B06.r2.fqTrTr.sorted.fq.gz`
  - rmRep BAM: `EXAMPLE_PE.rep1_clip.B06.r1.fq.genome-mapped.bam`
- Input:
  - R1: `EXAMPLE_PE.rep1_input.NIL.r1.fqTrTr.sorted.fq.gz`
  - R2: `EXAMPLE_PE.rep1_input.NIL.r2.fqTrTr.sorted.fq.gz`
  - rmRep BAM: `EXAMPLE_PE.rep1_input.NIL.r1.fq.genome-mapped.bam`

**References:** Same as SE workflow (GRCh38)

#### PE Simplified: `repeat_mapping_PE.simple.yaml`
**Dataset:** `peCLIP_example`
**Mode:** PE
Simplified format version of the above

---

## Reproduction Instructions

### Prerequisites

1. **Python Environment:**
   ```bash
   conda env create -f conda-env/repmap.yaml
   conda activate repmap
   ```

2. **Reference Data:**
   - Download or obtain access to:
     - Gencode GTF (v33 for hg38 or appropriate version for hg19)
     - RepeatMasker BED files
     - MASTER_FILELIST annotations
     - Bowtie2 indices for repetitive elements

### Reproducing Python Unit Tests

```bash
# 1. Create mini test fixtures (if needed)
python scripts/make_mini_fixtures.py

# 2. Run all unit tests
pytest tests/test_*.py -v

# 3. Run specific test suites
pytest tests/test_split_merge_python_expected.py -v
pytest tests/test_mini_fixtures_manifest.py -v
pytest tests/test_cwl_yaml_adapter.py -v
pytest tests/test_parse_pe_python.py -v
pytest tests/test_parse_se_python.py -v
```

### Reproducing CWL Integration Tests

**Note:** These tests require the full CWL workflow engine (cwltool or Toil)

```bash
# For individual tool tests (example: SE mapping)
cd tests/ecliprepmap-1.0.0/00_map_repetitive_elements_se/
cwltool <path_to_cwl_tool> map_repetitive_elements_se.yaml

# For complete workflow tests
cd tests/ecliprepmap-1.0.0/wf_ecliprepmap_se/
bash run_diff.sh  # Runs workflow and compares output to expected
```

### Reproducing Snakemake Workflow Tests

**Current Implementation:** The repository now uses Snakemake instead of CWL

```bash
# SE workflow
snakemake --snakefile workflow/Snakefile \
  --configfile examples/repeat_mapping_SE.simple.yaml \
  --cores 4

# PE workflow
snakemake --snakefile workflow/Snakefile \
  --configfile examples/repeat_mapping_PE.simple.yaml \
  --cores 4
```

### Obtaining Test Data

**Option 1: Download from TSCC**
```bash
# If you have access to TSCC storage
rsync -avP /tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/ ./test_data/
rsync -avP /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/ ./refdata/
```

**Option 2: Use archived test data**
```bash
# Extract example archives
tar -xzf examples/data_for_repeat_mapping_hg19.tar.gz
tar -xzf examples/repeat-mapping-hg19-refdata.tar.gz
```

**Option 3: Generate mini fixtures from full data**
```bash
# If full test data is available
python scripts/make_mini_fixtures.py --source <full_data_path> --output tests/fixtures/mini/
```

---

## File Manifests

### Mini Fixtures Manifest

The file `tests/fixtures/mini/manifest.tsv` contains SHA256 checksums and sizes for all fixtures:

**Sample entries:**
```tsv
relative_path                                          size_bytes    sha256
tests/fixtures/mini/expected/input.parsed.mini.gz      6492          63a47234eff4e54612065eb4386bb24c2a817165884c15b3aebb5245e199148d
tests/fixtures/mini/expected/input.rmDup.sam.mini.gz   174824        479a3774b8867e41292aff93a23c539cd220549522a2244a49e1b6d509f7f4a1
tests/fixtures/mini/expected/ip.parsed.mini.gz         6292          fdc33e2b6f6a6b7df23f21c5bbffc87bfcd0e682d1309e541d703edf9e5c74dd
tests/fixtures/mini/expected/ip.reparsed.nopipes.mini.tsv.gz  7188  ab76d446267a23dad0f5d815ef0db2bc8d152155cbd7ac5d6ecf0be089ec5c72
```

### Expected Split Output Manifest

The file `tests/fixtures/mini/expected/split.expected.manifest.tsv` defines expected outputs from BAM splitting:

**Format:**
```tsv
filename              size_bytes    sha256
<output_file_1>       <size>        <checksum>
<output_file_2>       <size>        <checksum>
...
```

---

## Test Execution Summary

### Running All Tests

```bash
# Python unit tests (fast, ~seconds)
pytest tests/ -v

# Integration tests (slow, depends on data size)
# Must be run individually as documented above
```

### Expected Test Outcomes

1. **Python Unit Tests:** All should pass with matching SHA256 checksums
2. **Integration Tests:** Outputs should match expected files byte-for-byte
3. **Workflow Tests:** Final outputs should be scientifically valid (proper BAM format, expected read counts, etc.)

---

## Key Scientific Validations

1. **Deterministic Output:** All tools produce identical output given identical input
2. **SHA256 Integrity:** Critical files are checksummed to detect corruption
3. **Read Count Conservation:** Reads are not lost during splitting/merging
4. **Format Compliance:** Output BAM/SAM files pass validation
5. **Annotation Consistency:** Repetitive element families are correctly assigned

---

## Dependencies Graph

```
Python Tools:
├── split_bam_to_subfiles_SEorPE.py
│   ├── Input: SAM/BAM file + SE/PE mode
│   └── Output: Multiple split SAM files by prefix
├── merge_multiple_parsed_files.simplified_20191022.py
│   ├── Input: Multiple parsed TSV files
│   └── Output: Single merged parsed file
├── parse_bowtie2_output_realtime_includemultifamily_PE.py
│   ├── Input: Bowtie2 SAM + filelist
│   └── Output: Annotated SAM with repeat families
└── parse_bowtie2_output_realtime_includemultifamily_SE.py
    ├── Input: Bowtie2 SAM + filelist
    └── Output: Annotated SAM with repeat families

Reference Files:
├── bowtie2_index/ (repetitive elements database)
├── MASTER_FILELIST (repeat family annotations)
├── gencode.gtf (gene annotations)
├── gencode.table.tsv (gene table format)
└── repmask.bed (RepeatMasker annotations)
```

---

## Notes for Reproduction

1. **Paths:** All absolute paths in YAML configs must be updated to match your environment
2. **Python Version:** Requires Python 3.7+ (type hints and pathlib)
3. **External Tools:** Bowtie2, samtools may be required for some tests
4. **Temporary Files:** Tests use `tmp_path` pytest fixture for clean ephemeral storage
5. **Gzipped Inputs:** Many fixtures are gzipped to save space; tests handle decompression

---

## Storage Optimization Notes

- **Mini fixtures are essential:** Keep `tests/fixtures/mini/` (~few MB)
- **Large test data is reproducible:** The 130GB in `tests/ecliprepmap-1.0.0/` can be regenerated by re-running workflows
- **Example archives contain redundant data:** Extracted and compressed versions both exist
- **Legacy outputs are for reference only:** `_eric_outputs/` can be regenerated if needed

---

## Metadata

**Total test coverage:**
- 5 Python unit test files
- 13 CWL integration test directories
- 2 complete workflow tests (SE + PE)
- 4 example workflow configurations
- 2 genome builds tested (hg19, hg38)

**Test modes:**
- SE (single-end)
- PE (paired-end)
- Mock/mini fixtures (fast)
- Full-scale integration (slow)

---

## Contact & Provenance

**Original Pipeline:** ecliprepmap 1.0.0
**Developer:** Eric Van Nostrand (Yeo Lab)
**Conversion to Snakemake:** 2026 (current branch: codex/python-conversion-python-only-cleanup)
**Repository:** /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
**Documentation Generated:** 2026-04-15
