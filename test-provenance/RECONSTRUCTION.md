# Test Data Reconstruction Guide

**Purpose:** Step-by-step instructions to rebuild the complete test environment from this provenance archive

---

## Prerequisites

### Software Requirements

1. **Conda/Mamba:**
   ```bash
   # Install miniconda if needed
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. **Git:**
   ```bash
   git --version  # Should be 2.0+
   ```

3. **Python 3.7+** (provided by conda env)

4. **Optional: Snakemake** (for workflow tests)
   ```bash
   conda install -c conda-forge -c bioconda snakemake
   ```

---

## Reconstruction Steps

### Step 1: Set Up Repository

```bash
# Clone the repository
git clone <repository_url>
cd repetitive-element-mapping

# Check out the correct branch
git checkout codex/python-conversion-python-only-cleanup

# Verify test-provenance archive is present
ls -lh test-provenance/
```

---

### Step 2: Install Python Environment

```bash
# Create conda environment
conda env create -f conda-env/repmap.yaml

# Activate environment
conda activate repmap

# Verify pytest is available
pytest --version
```

---

### Step 3: Reconstruct Test Directories

```bash
# Create test directories from provenance
cp -r test-provenance/tests .
cp -r test-provenance/examples .
```

---

### Step 4: Obtain Reference Data

#### Option A: Access TSCC Storage (If Available)

```bash
# Create reference directory
mkdir -p references/hg38

# Copy bowtie2 index
rsync -avP /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/bowtie2_index/ \
  references/hg38/bowtie2_index/

# Copy annotation files
rsync -avP /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/MASTER_FILELIST.* \
  references/hg38/

rsync -avP /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/UniqueGenomicElements.hg38.bed \
  references/hg38/

# Copy Gencode annotations
rsync -avP /tscc/projects/ps-yeolab4/genomes/hg38/gencode/v33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf* \
  references/hg38/
```

#### Option B: Download from Public Sources

```bash
# Download Gencode v33 (hg38/GRCh38)
cd references/hg38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz

# Generate UCSC table format
# (script needed, or use existing from TSCC)

# Download RepeatMasker (if publicly available)
# Or build from UCSC Genome Browser data
```

#### Option C: Build Bowtie2 Index from Scratch

```bash
# This requires the MASTER_FILELIST.*.fa file
# Contact repository maintainers for access or generation script

bowtie2-build MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat \
  bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat
```

---

### Step 5: Obtain Test Data

#### Option A: Use Existing TSCC Data

```bash
# SE test data
mkdir -p test-data/se
rsync -avP /tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/*.fq.gz \
  test-data/se/
rsync -avP /tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/*.bam \
  test-data/se/
```

#### Option B: Generate Mini Fixtures

```bash
# If you have full test data, generate mini versions
python scripts/make_mini_fixtures.py \
  --source test-data/se/ \
  --output tests/fixtures/mini/ \
  --subsample 10000  # Create small test files
```

#### Option C: Use Pre-Generated Mini Fixtures

```bash
# Download from repository releases or shared storage
# (Coordinates to be provided by repository maintainers)

# Example:
wget https://example.com/test-fixtures/mini-fixtures.tar.gz
tar -xzf mini-fixtures.tar.gz -C tests/fixtures/
```

---

### Step 6: Update Configuration Paths

Many YAML configs contain hardcoded absolute paths. Update them:

```bash
# Option A: Manual editing
vim examples/repeat_mapping_SE.simple.yaml
# Update paths to match your environment

# Option B: Automated replacement (careful!)
find examples/ tests/ -name "*.yaml" -type f -exec sed -i \
  's|/tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs|references/hg38|g' {} \;

find examples/ tests/ -name "*.yaml" -type f -exec sed -i \
  's|/tscc/projects/ps-yeolab4/genomes/hg38/gencode/v33|references/hg38|g' {} \;
```

---

### Step 7: Validate Mini Fixtures

```bash
# Run manifest validation
pytest tests/test_mini_fixtures_manifest.py -v

# Expected output:
# tests/test_mini_fixtures_manifest.py::test_mini_fixture_manifest_matches_files PASSED
```

---

### Step 8: Run Python Unit Tests

```bash
# Run all unit tests
pytest tests/test_*.py -v

# Run individually for debugging
pytest tests/test_split_merge_python_expected.py -v
pytest tests/test_cwl_yaml_adapter.py -v
pytest tests/test_parse_pe_python.py -v
pytest tests/test_parse_se_python.py -v
```

**Expected Results:**
- All tests should PASS
- SHA256 checksums should match expected values
- Output files should be byte-for-byte identical

---

### Step 9: Run Snakemake Workflows

```bash
# Test SE workflow (simplified config)
snakemake --snakefile workflow/Snakefile \
  --configfile examples/repeat_mapping_SE.simple.yaml \
  --cores 4 \
  --dry-run  # Check workflow plan first

# Run for real
snakemake --snakefile workflow/Snakefile \
  --configfile examples/repeat_mapping_SE.simple.yaml \
  --cores 4

# Test PE workflow (simplified config)
snakemake --snakefile workflow/Snakefile \
  --configfile examples/repeat_mapping_PE.simple.yaml \
  --cores 4
```

---

## Troubleshooting

### Issue: Missing Mini Fixtures

**Solution:**
```bash
# Check if fixture creation script exists
ls scripts/make_mini_fixtures.py

# If not, mini fixtures must be obtained from maintainers or TSCC
```

---

### Issue: Bowtie2 Index Missing

**Error:** `bowtie2-inspect: not found` or index files missing

**Solution:**
```bash
# Verify bowtie2 is installed
conda install -c bioconda bowtie2

# Check index files exist
ls references/hg38/bowtie2_index/*.bt2
```

---

### Issue: Hardcoded Paths in Configs

**Error:** `FileNotFoundError: /tscc/projects/...`

**Solution:**
```bash
# Use provided path update script (if available)
python scripts/update_config_paths.py \
  --old-prefix /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs \
  --new-prefix references/hg38 \
  examples/*.yaml tests/ecliprepmap-1.0.0/*/*.yaml
```

---

### Issue: SHA256 Mismatch

**Error:** `AssertionError: sha256 mismatch: tests/fixtures/mini/expected/ip.parsed.mini.gz`

**Cause:** Fixture files were corrupted or incorrectly generated

**Solution:**
```bash
# Regenerate fixtures from known-good source
python scripts/make_mini_fixtures.py --source <original_data> --output tests/fixtures/mini/ --force
```

---

## Data Size Requirements

### Minimal (Unit Tests Only)
- **Size:** ~100MB
- **Contents:** Mini fixtures only
- **Tests:** Python unit tests

### Medium (Workflow Tests)
- **Size:** ~10GB
- **Contents:** Mini fixtures + subsampled test data
- **Tests:** Unit tests + basic workflow tests

### Complete (Full Integration)
- **Size:** ~150GB
- **Contents:** All original test data + references
- **Tests:** All tests including CWL integration tests

---

## Validation Checklist

After reconstruction, verify:

- [ ] Conda environment activates successfully
- [ ] All Python unit tests pass
- [ ] Mini fixture manifest validation passes
- [ ] Bowtie2 index is accessible
- [ ] Reference annotations are present
- [ ] At least one workflow runs to completion
- [ ] Output files match expected formats

---

## Alternative: Using Beads (Reproducibility Tool)

If the repository includes a `beads.input` file:

```bash
# Install beads
pip install beads-bio

# Reconstruct all inputs
beads make

# This will download/symlink all required data files
# based on checksums and metadata
```

See `BEADS.md` for more details.

---

## Contact and Support

For issues with reconstruction:
1. Check TEST_DOCUMENTATION.md for detailed test descriptions
2. Review git commit history for recent changes
3. Contact repository maintainers
4. Open an issue on GitHub with reconstruction logs

---

## Success Criteria

You have successfully reconstructed the test environment when:

1. `pytest tests/` runs without errors
2. At least one Snakemake workflow completes
3. Output files pass format validation
4. Key checksum files match expected values

---

## Appendix: File Checksums

### Critical Files (for validation)

```bash
# Generate checksums after reconstruction
sha256sum tests/fixtures/mini/expected/*.gz > checksums.txt

# Compare with original manifest
diff checksums.txt tests/fixtures/mini/manifest.tsv
```

### Expected Checksums (from manifest)

```
63a47234eff4e54612065eb4386bb24c2a817165884c15b3aebb5245e199148d  tests/fixtures/mini/expected/input.parsed.mini.gz
479a3774b8867e41292aff93a23c539cd220549522a2244a49e1b6d509f7f4a1  tests/fixtures/mini/expected/input.rmDup.sam.mini.gz
fdc33e2b6f6a6b7df23f21c5bbffc87bfcd0e682d1309e541d703edf9e5c74dd  tests/fixtures/mini/expected/ip.parsed.mini.gz
ab76d446267a23dad0f5d815ef0db2bc8d152155cbd7ac5d6ecf0be089ec5c72  tests/fixtures/mini/expected/ip.reparsed.nopipes.mini.tsv.gz
```

---

**Last Updated:** 2026-04-15
**Reconstruction Tested:** No (awaiting first reconstruction attempt)
