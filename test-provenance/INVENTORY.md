# Test Provenance Archive Inventory

**Generated:** 2026-04-15
**Archive Size:** ~692KB
**Original Size:** ~147GB
**Compression Ratio:** ~212,000:1

---

## Archive Contents

### Documentation Files

| File | Size | Purpose |
|------|------|---------|
| README.md | - | Archive overview and quick reference |
| RECONSTRUCTION.md | - | Step-by-step reconstruction instructions |
| BEADS.md | - | Beads/DVC usage guide |
| INVENTORY.md | - | This file - complete inventory |
| beads.input | - | Beads manifest with checksums |

### Preserved Test Code

| File | Lines | Purpose |
|------|-------|---------|
| tests/test_split_merge_python_expected.py | 61 | Tests split/merge Python tools |
| tests/test_mini_fixtures_manifest.py | 34 | Validates fixture integrity |
| tests/test_cwl_yaml_adapter.py | 83 | Tests YAML config adapters |
| tests/test_parse_pe_python.py | 53 | Tests PE bowtie2 parser |
| tests/test_parse_se_python.py | ? | Tests SE bowtie2 parser |

### Preserved Test Configurations

#### Mini Fixtures (`tests/fixtures/mini/`)
- `cwl_job_example.yaml` - CWL-style test configuration
- `manifest.tsv` - File integrity manifest with checksums

#### Drop-in Fixtures (`tests/fixtures/dropin/`)
- `se_job.yaml` - Mock SE job (CWL style)
- `simple_se_job.yaml` - Simple SE job (new style)
- `pe_job.yaml` - Mock PE job (CWL style)
- `simple_pe_job.yaml` - Simple PE job (new style)

#### CWL Integration Tests (`tests/ecliprepmap-1.0.0/`)

**SE Pipeline Steps:**
| Directory | Config File | Tool |
|-----------|-------------|------|
| 00_map_repetitive_elements_se/ | map_repetitive_elements_se.yaml | Map repeats (SE) |
| 01_split_rep_bam_se/ | splitbam.yaml | Split repeat BAM (SE) |
| 02_split_rmrep_bam_se/ | splitbam.yaml | Split genome BAM (SE) |
| 03_dedup_se/ | deduplicate_se.yaml, run_script.sh | Dedup (SE) |

**PE Pipeline Steps:**
| Directory | Config File | Tool |
|-----------|-------------|------|
| 05_map_repetitive_elements_pe/ | map_repetitive_elements_pe.yaml | Map repeats (PE) |
| 06_split_rep_sam_pe/ | splitbam.yaml | Split repeat SAM (PE) |
| 07_split_rmrep_bam_pe/ | splitbam.yaml | Split genome BAM (PE) |
| 08_dedup_pe/ | deduplicate_pe.yaml | Dedup (PE) |

**Shared:**
| Directory | Config File | Tool |
|-----------|-------------|------|
| 04_merge_parsed/ | combine.yaml | Merge parsed files |

**Complete Workflows:**
| Directory | Config File | Job Input |
|-----------|-------------|-----------|
| wf_ecliprepmap_se/ | wf_ecliprepmap_se.yaml | REPELEMENTMAPPING_wf_ecliprepmap_se_INPUT.yaml |
| wf_ecliprepmap_pe/ | wf_ecliprepmap_pe.yaml | REPELEMENTMAPPING_wf_ecliprepmap_pe_INPUT.yaml |

### Preserved Example Configurations

| File | Dataset | Mode | Purpose |
|------|---------|------|---------|
| repeat_mapping_SE.yaml | seCLIP_example | SE | Full SE workflow (CWL style) |
| repeat_mapping_SE.simple.yaml | seCLIP_example | SE | Simplified SE workflow |
| repeat_mapping_PE.yaml | peCLIP_example | PE | Full PE workflow (CWL style) |
| repeat_mapping_PE.simple.yaml | peCLIP_example | PE | Simplified PE workflow |

---

## What Was Excluded (147GB)

### Test Data Files

**Mini Fixtures - Source:** (~50MB)
- `tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz` - 7.6MB
- `tests/fixtures/mini/source/input.preRmDup.sam.mini.gz` - 7.3MB
- `tests/fixtures/mini/source/merge_input_*.parsed_v2.txt` - Small

**Mini Fixtures - Expected:** (~15MB)
- `tests/fixtures/mini/expected/*.gz` - Various parsed/dedup outputs
- `tests/fixtures/mini/expected/split.expected.manifest.tsv`

**Mini Fixtures - References:** (~1GB estimated)
- `tests/fixtures/mini/refs/bowtie2_index/` - Bowtie2 indices
- `tests/fixtures/mini/refs/*.gtf` - Gene annotations
- `tests/fixtures/mini/refs/*.bed` - RepeatMasker annotations

**CWL Integration Test Data:** (~129GB estimated)
- `tests/ecliprepmap-1.0.0/*/inputs/` - Test input files
- `tests/ecliprepmap-1.0.0/*/outputs/` - Expected outputs
- `tests/ecliprepmap-1.0.0/_eric_outputs/` - Reference outputs

**Example Data Archives:**
- `examples/data_for_repeat_mapping_hg19.tar.gz` - 4.8GB
- `examples/repeat-mapping-hg19-refdata.tar.gz` - 102MB
- `examples/example_data_for_repeat_mapping_hg38.tar.gz` - 5.2GB
- `examples/example_data_for_repeat_mapping_hg38/` - Extracted PE data
- `examples/repeat-mapping-grch38-refdata/` - Extracted references

---

## Reconstruction Requirements

### Minimal (Unit Tests Only)
**Size:** ~100MB
**Files needed:**
- Mini fixture source files (4 files, ~15MB)
- Mini fixture expected outputs (11 files, ~35MB)
- Mini fixture references (~50MB)

**What works:**
- All Python unit tests
- Fixture integrity validation

### Medium (Basic Workflows)
**Size:** ~10GB
**Files needed:**
- Minimal requirements above
- SE test data (4 files, ~5GB)
- Reference bowtie2 index (~4GB)
- Annotation files (~1GB)

**What works:**
- Unit tests
- SE workflow end-to-end
- Basic integration tests

### Complete (All Tests)
**Size:** ~150GB
**Files needed:**
- Medium requirements above
- PE test data
- All CWL integration test data
- All example archives
- Legacy reference outputs

**What works:**
- All unit tests
- All integration tests
- All workflow tests
- Full regression testing

---

## File Checksums (Critical Files)

### Mini Fixture Sources
```
f6749227444bc327255cd98fe7e54be5c7c19daaacc90dd15ec2fbb996f9f9e6  tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz
24f3307ebfce353325bdea0813a54b2201ef83da12059af8aa3292a4fe705f2b  tests/fixtures/mini/source/input.preRmDup.sam.mini.gz
3c0a90323c04fd3fb853f1acb45206390e68908e1042e7073bf4b6b1ae587b42  tests/fixtures/mini/source/merge_input_1.parsed_v2.txt
7e1faae81ff307f2522d2576ad4c4bde583e18be7e1eaa3b6f9d57eaf9338ef5  tests/fixtures/mini/source/merge_input_2.parsed_v2.txt
```

### Mini Fixture Expected Outputs
```
fdc33e2b6f6a6b7df23f21c5bbffc87bfcd0e682d1309e541d703edf9e5c74dd  tests/fixtures/mini/expected/ip.parsed.mini.gz
63a47234eff4e54612065eb4386bb24c2a817165884c15b3aebb5245e199148d  tests/fixtures/mini/expected/input.parsed.mini.gz
920ce5732417fcb46d9fbcf3c348e2dbf99f4d8fd74fe8c3248231f9855d77e8  tests/fixtures/mini/expected/ip.rmDup.sam.mini.gz
479a3774b8867e41292aff93a23c539cd220549522a2244a49e1b6d509f7f4a1  tests/fixtures/mini/expected/input.rmDup.sam.mini.gz
ab76d446267a23dad0f5d815ef0db2bc8d152155cbd7ac5d6ecf0be089ec5c72  tests/fixtures/mini/expected/ip.reparsed.nopipes.mini.tsv.gz
db976f71def48a10632aeaaa71226a6c19260c254da85c7f9e0afb0de3339255  tests/fixtures/mini/expected/ip.reparsed.withpipes.mini.tsv.gz
bce73da7c0297c605f646b17ed4b9ec0a2c2f2ecfc693d03601ff0799cbcbdbc  tests/fixtures/mini/expected/merged.expected.parsed
```

### Configuration Files
```
0897b8769657d6c5b5fbc7ca34af6759c4903c30550eb2eb6375ee79efad26ce  tests/fixtures/mini/cwl_job_example.yaml
873c9c964d70097c2ccd3299500ff6c9a2fd753c4727137a4885aaadd7cc2450  tests/fixtures/mini/manifest.tsv
e272596ba514c5b61355d405b9457fca595f516fc7646d6f3072e82967e5c642  tests/fixtures/mini/expected/split.expected.manifest.tsv
```

---

## Dependencies

### Python Packages
See `../conda-env/repmap.yaml` for complete list.

**Key dependencies:**
- pytest >= 5.0
- pyyaml
- biopython (for some tests)

### External Tools
- bowtie2 >= 2.3.0
- samtools >= 1.9 (optional)
- snakemake >= 6.0 (for workflow tests)

### Reference Data Sources

**Gencode:**
- v33 for hg38: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/

**Bowtie2 Index:**
- Custom repetitive element database
- Contact: Yeo Lab, UCSD
- TSCC path: `/tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/bowtie2_index/`

**Annotation Files:**
- MASTER_FILELIST: Custom repeat family annotations
- RepeatMasker: Derived from UCSC Genome Browser
- TSCC paths documented in beads.input

---

## Version Information

### Git Metadata
- **Branch:** codex/python-conversion-python-only-cleanup
- **Original Repository:** /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
- **Archival Date:** 2026-04-15

### Pipeline Versions
- **ecliprepmap:** 1.0.0 (CWL version, legacy)
- **Current:** Snakemake conversion (in progress)

### Test Data Versions
- **hg19 data:** 2020-08-25
- **hg38 data:** 2021-02-13
- **Mini fixtures:** Derived from full data, last updated 2021

---

## Usage Statistics

### Directory Sizes (Original)
```
130G    tests/
  ~100G   tests/ecliprepmap-1.0.0/
  ~25G    tests/fixtures/
  ~5G     tests/ecliprepmap-1.0.0/_eric_outputs/

17G     examples/
  4.8G    examples/data_for_repeat_mapping_hg19.tar.gz
  5.2G    examples/example_data_for_repeat_mapping_hg38.tar.gz
  102M    examples/repeat-mapping-hg19-refdata.tar.gz
  ~6G     examples/example_data_for_repeat_mapping_hg38/
  ~1G     examples/repeat-mapping-grch38-refdata/
```

### File Counts (Preserved)
- **Python test files:** 5
- **YAML config files:** ~30
- **Shell scripts:** ~5
- **TSV manifests:** ~3
- **Documentation:** 5 (including this file)

**Total files in archive:** ~50 lightweight files

---

## Verification

### Quick Verification

```bash
# Check archive structure
tree test-provenance/ -L 2

# Verify file count
find test-provenance/ -type f | wc -l

# Check documentation
cat test-provenance/README.md
```

### Full Verification

```bash
# Verify all checksums
grep -v "^#" test-provenance/beads.input | while IFS="|" read -r path checksum source; do
    path=$(echo "$path" | xargs)
    checksum=$(echo "$checksum" | xargs)
    if [ -f "$path" ] && [ "$checksum" != "CHECKSUM_NEEDED" ]; then
        current=$(sha256sum "$path" | awk '{print $1}')
        if [ "$current" == "$checksum" ]; then
            echo "OK: $path"
        else
            echo "FAIL: $path"
        fi
    fi
done
```

---

## Contact

For questions about this archive or test reconstruction:
- See TEST_DOCUMENTATION.md in repository root
- Contact repository maintainers
- File issue on GitHub

---

**Archive Status:** Complete and ready for distribution
**Last Verified:** 2026-04-15
