# Using Beads for Reproducible Test Data

**Beads** is a tool for reproducible computational research that manages input data files using checksums and manifests.

---

## Installation

```bash
pip install beads-bio
```

**Note:** As of 2026-04-15, beads may not be available as a standard package. Alternative approaches:
- Use manual reconstruction (see RECONSTRUCTION.md)
- Implement a custom data management solution
- Use DVC (Data Version Control) or similar tools

---

## Concept

The `beads.input` file lists all required input data with:
1. **Output path** - Where the file should be placed
2. **Checksum** - SHA256 hash for validation
3. **Source** - URL or path to download/copy from

Beads ensures:
- Files are fetched from source only if missing or corrupted
- All files match expected checksums
- Reproducibility across different environments

---

## Usage

### Basic Commands

```bash
# Fetch all required data files
beads make

# Check status of data files
beads status

# Clean generated files
beads clean
```

---

## Customization for TSCC

The `beads.input` file includes custom `tscc:` protocol for TSCC-specific paths:

```
tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz | CHECKSUM | tscc:/tscc/projects/...
```

To enable this, you may need to configure a custom fetcher:

### Option 1: Symlink (if running on TSCC)

```bash
# Create custom beads config
cat > .beads_config <<EOF
[fetchers]
tscc = symlink
EOF

# Now "beads make" will create symlinks instead of copying
```

### Option 2: Copy (if running on TSCC with writable destination)

```bash
# Create custom fetcher script
cat > .beads/fetchers/tscc.sh <<'EOF'
#!/bin/bash
# Custom fetcher for TSCC paths
SOURCE_PATH="${1#tscc:}"  # Remove "tscc:" prefix
DEST_PATH="$2"
cp "$SOURCE_PATH" "$DEST_PATH"
EOF
chmod +x .beads/fetchers/tscc.sh
```

### Option 3: Rsync (if running remotely)

```bash
# Create custom fetcher for remote TSCC access
cat > .beads/fetchers/tscc.sh <<'EOF'
#!/bin/bash
SOURCE_PATH="${1#tscc:}"
DEST_PATH="$2"
rsync -avP username@tscc-login.sdsc.edu:"$SOURCE_PATH" "$DEST_PATH"
EOF
chmod +x .beads/fetchers/tscc.sh
```

---

## Alternative: DVC (Data Version Control)

If beads is not available, consider using DVC:

```bash
# Install DVC
pip install dvc

# Initialize DVC
dvc init

# Add test data to DVC tracking
dvc add tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz
dvc add tests/fixtures/mini/expected/

# Commit DVC metadata
git add tests/fixtures/mini/source/*.dvc tests/fixtures/mini/expected/*.dvc .dvc/
git commit -m "Add test data to DVC"

# Configure remote storage
dvc remote add -d tscc /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping/dvc-cache

# Push data to remote
dvc push

# On another machine, pull data
dvc pull
```

---

## Manual Reconstruction (Without Beads)

If neither beads nor DVC is available:

```bash
# Create a simple manifest-based fetcher
cat > fetch_test_data.sh <<'EOF'
#!/bin/bash
# Read beads.input and fetch files

grep -v "^#" test-provenance/beads.input | while IFS="|" read -r dest checksum source; do
    dest=$(echo "$dest" | xargs)  # trim whitespace
    checksum=$(echo "$checksum" | xargs)
    source=$(echo "$source" | xargs)

    # Skip if checksum not yet computed
    if [[ "$checksum" == "CHECKSUM_NEEDED" ]]; then
        echo "SKIP (checksum needed): $dest"
        continue
    fi

    # Check if file exists and matches checksum
    if [ -f "$dest" ]; then
        current_checksum=$(sha256sum "$dest" | awk '{print $1}')
        if [ "$current_checksum" == "$checksum" ]; then
            echo "OK: $dest"
            continue
        fi
    fi

    # Fetch file
    mkdir -p "$(dirname "$dest")"

    if [[ "$source" == tscc:* ]]; then
        source_path="${source#tscc:}"
        echo "COPY: $source_path -> $dest"
        cp "$source_path" "$dest"
    elif [[ "$source" == http* ]]; then
        echo "DOWNLOAD: $source -> $dest"
        wget -O "$dest" "$source"
        # If source is .gz and dest is not, decompress
        if [[ "$source" == *.gz ]] && [[ "$dest" != *.gz ]]; then
            gunzip "$dest"
        fi
    else
        echo "ERROR: Unknown source type: $source"
    fi

    # Validate checksum
    current_checksum=$(sha256sum "$dest" | awk '{print $1}')
    if [ "$current_checksum" != "$checksum" ]; then
        echo "ERROR: Checksum mismatch for $dest"
        echo "  Expected: $checksum"
        echo "  Got:      $current_checksum"
    fi
done
EOF

chmod +x fetch_test_data.sh
./fetch_test_data.sh
```

---

## Generating Missing Checksums

Some entries in `beads.input` have `CHECKSUM_NEEDED` placeholders. Generate them:

```bash
# On TSCC (where original files exist)
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping

# Generate checksums for example data
sha256sum /tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/*.fq.gz
sha256sum /tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/*.bam

# Generate checksums for reference data
sha256sum /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/bowtie2_index/*.bt2
sha256sum /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/*.list
sha256sum /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/*.bed
sha256sum /tscc/projects/ps-yeolab4/genomes/hg38/gencode/v33/*.gtf*

# Update beads.input with computed checksums
```

---

## Verifying Data Integrity

After fetching data:

```bash
# Verify all files match checksums
grep -v "^#" test-provenance/beads.input | grep -v "CHECKSUM_NEEDED" | while IFS="|" read -r dest checksum source; do
    dest=$(echo "$dest" | xargs)
    checksum=$(echo "$checksum" | xargs)

    if [ ! -f "$dest" ]; then
        echo "MISSING: $dest"
        continue
    fi

    current=$(sha256sum "$dest" | awk '{print $1}')
    if [ "$current" == "$checksum" ]; then
        echo "OK: $dest"
    else
        echo "FAIL: $dest (expected $checksum, got $current)"
    fi
done
```

---

## Best Practices

1. **Always verify checksums** after fetching data
2. **Store large files externally** (not in git)
3. **Use beads.input as source of truth** for data provenance
4. **Document custom fetchers** if using non-standard protocols
5. **Keep beads.input updated** when adding new test data

---

## Troubleshooting

### Error: "Checksum mismatch"

**Cause:** File was corrupted during transfer or source file changed

**Solution:**
```bash
# Re-fetch the file
rm <problematic_file>
beads make
```

---

### Error: "Cannot fetch tscc: protocol"

**Cause:** Custom fetcher not configured

**Solution:**
- Configure custom fetcher (see above)
- Or manually copy files from TSCC
- Or update beads.input with accessible URLs

---

## See Also

- RECONSTRUCTION.md - Full manual reconstruction guide
- TEST_DOCUMENTATION.md - Complete test documentation
- README.md - Provenance archive overview

---

**Last Updated:** 2026-04-15
