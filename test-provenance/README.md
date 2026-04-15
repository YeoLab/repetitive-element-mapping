# Test Provenance Archive

**Created:** 2026-04-15
**Purpose:** Preserve essential metadata and configurations to enable complete reproduction of tests and examples after deletion of large data directories

---

## Contents

This archive contains **all configuration files, scripts, manifests, and metadata** from the `tests/` and `examples/` directories, but **excludes large binary data files** (147GB total).

### Directory Structure

```
test-provenance/
├── README.md                          # This file
├── RECONSTRUCTION.md                  # Instructions to rebuild test data
├── configs/                           # (reserved for additional configs)
├── manifests/                         # (reserved for file manifests)
├── reference-metadata/                # (reserved for reference file metadata)
├── scripts/                           # (reserved for helper scripts)
├── tests/                             # All test configurations and code
│   ├── *.py                          # Python unit test files
│   ├── ecliprepmap-1.0.0/            # CWL integration test configs
│   └── fixtures/                     # Test fixture configs
└── examples/                          # Example workflow configurations
    ├── repeat_mapping_SE.yaml
    ├── repeat_mapping_SE.simple.yaml
    ├── repeat_mapping_PE.yaml
    └── repeat_mapping_PE.simple.yaml
```

### Archive Size

**Total:** ~692KB (compared to 147GB of original data)
**Compression ratio:** ~212,000:1

---

## What Was Excluded

- Binary data files (FASTQ, BAM, SAM)
- Large reference archives (*.tar.gz)
- Generated test outputs
- Bowtie2 indices
- Jupyter notebook checkpoints

---

## What Is Preserved

✅ All YAML/YML configuration files
✅ All Python test scripts
✅ All shell scripts
✅ All TSV manifests with SHA256 checksums
✅ Directory structure and naming conventions
✅ Complete documentation (see TEST_DOCUMENTATION.md in repo root)

---

## Reconstruction Workflow

To fully reproduce the test environment:

1. **Install dependencies:**
   ```bash
   conda env create -f ../conda-env/repmap.yaml
   conda activate repmap
   ```

2. **Generate mini test fixtures:**
   ```bash
   # Option A: From full data (if available)
   python ../scripts/make_mini_fixtures.py --source <full_data> --output tests/fixtures/mini/

   # Option B: Download pre-generated fixtures
   # See RECONSTRUCTION.md for download instructions
   ```

3. **Run Python unit tests:**
   ```bash
   pytest tests/test_*.py -v
   ```

4. **Run integration tests:**
   ```bash
   # See TEST_DOCUMENTATION.md for detailed instructions
   snakemake --configfile examples/repeat_mapping_SE.simple.yaml
   ```

---

## File Inventory

### Python Unit Tests (tests/)
- `test_split_merge_python_expected.py` - Validates split/merge tools
- `test_mini_fixtures_manifest.py` - Validates fixture integrity
- `test_cwl_yaml_adapter.py` - Tests YAML config parsing
- `test_parse_pe_python.py` - Tests PE parser
- `test_parse_se_python.py` - Tests SE parser

### Test Fixtures (tests/fixtures/)
- `mini/cwl_job_example.yaml` - CWL-style job config
- `mini/manifest.tsv` - File checksums
- `dropin/se_job.yaml` - Simplified SE config
- `dropin/pe_job.yaml` - Simplified PE config
- `dropin/simple_se_job.yaml` - Mock SE config
- `dropin/simple_pe_job.yaml` - Mock PE config

### CWL Integration Tests (tests/ecliprepmap-1.0.0/)
Each directory contains one or more `.yaml` config files:
- `00_map_repetitive_elements_se/` - SE mapping
- `01_split_rep_bam_se/` - SE repeat split
- `02_split_rmrep_bam_se/` - SE genome split
- `03_dedup_se/` - SE deduplication
- `04_merge_parsed/` - Merge parsed files
- `05_map_repetitive_elements_pe/` - PE mapping
- `06_split_rep_sam_pe/` - PE repeat split
- `07_split_rmrep_bam_pe/` - PE genome split
- `08_dedup_pe/` - PE deduplication
- `wf_ecliprepmap_se/` - Complete SE workflow
- `wf_ecliprepmap_pe/` - Complete PE workflow

### Example Workflows (examples/)
- `repeat_mapping_SE.yaml` - Full SE workflow config
- `repeat_mapping_SE.simple.yaml` - Simplified SE workflow
- `repeat_mapping_PE.yaml` - Full PE workflow config
- `repeat_mapping_PE.simple.yaml` - Simplified PE workflow

---

## Data Sources

### Original Test Data Locations (TSCC)

**SE Test Data:**
```
/tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/
├── INV_B.IP.umi.r1.fqTrTr.sorted.fq.gz
├── INV_B.IP.umi.r1.fq.genome-mappedSoSo.bam
├── INV_B.IN.umi.r1.fqTrTr.sorted.fq.gz
└── INV_B.IN.umi.r1.fq.genome-mappedSoSo.bam
```

**Reference Data:**
```
/tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/
├── bowtie2_index/
├── MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list
└── UniqueGenomicElements.hg38.bed

/tscc/projects/ps-yeolab4/genomes/hg38/gencode/v33/
├── gencode.v33.chr_patch_hapl_scaff.annotation.gtf
└── gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
```

---

## Checksums and Validation

Original fixtures include SHA256 checksums in:
- `tests/fixtures/mini/manifest.tsv`
- `tests/fixtures/mini/expected/split.expected.manifest.tsv`

These can be used to validate regenerated fixtures.

---

## Version Control

**Git branch:** codex/python-conversion-python-only-cleanup
**Archived from commit:** (see git log at time of archival)
**Original repository:** /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

---

## Contact

For questions about test reproduction or data access:
- Repository: https://github.com/[your-org]/repetitive-element-mapping
- Documentation: See TEST_DOCUMENTATION.md in repository root
