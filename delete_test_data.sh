#!/bin/bash
# Script to safely delete 147GB of test data after documentation
# Generated: 2026-04-15
# Run with: bash delete_test_data.sh --dry-run
# Then: bash delete_test_data.sh --confirm

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

DRY_RUN=true

# Parse arguments
if [ "$1" == "--confirm" ]; then
    DRY_RUN=false
elif [ "$1" == "--dry-run" ]; then
    DRY_RUN=true
else
    echo -e "${YELLOW}Usage: $0 [--dry-run|--confirm]${NC}"
    echo "  --dry-run  : Show what would be deleted (default)"
    echo "  --confirm  : Actually delete files"
    exit 1
fi

# Check we're in the right directory
if [ ! -f "TEST_DOCUMENTATION.md" ] || [ ! -d "test-provenance" ]; then
    echo -e "${RED}ERROR: Must run from repository root${NC}"
    echo "Expected files: TEST_DOCUMENTATION.md, test-provenance/"
    exit 1
fi

# Check git status
if ! git diff-index --quiet HEAD -- 2>/dev/null; then
    echo -e "${YELLOW}WARNING: You have uncommitted changes${NC}"
    git status --short
    echo ""
    read -p "Continue anyway? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Check if commit exists
if ! git log --oneline | grep -q "Add comprehensive test documentation"; then
    echo -e "${RED}ERROR: Required commit not found${NC}"
    echo "Expected commit: 'Add comprehensive test documentation and provenance archive'"
    echo "Current commits:"
    git log --oneline -5
    exit 1
fi

echo -e "${GREEN}Pre-flight checks passed${NC}"
echo ""

# Define directories/files to delete
TARGETS=(
    # Test data directories
    "tests/ecliprepmap-1.0.0/00_map_repetitive_elements_se/inputs"
    "tests/ecliprepmap-1.0.0/00_map_repetitive_elements_se/outputs"
    "tests/ecliprepmap-1.0.0/01_split_rep_bam_se/inputs"
    "tests/ecliprepmap-1.0.0/01_split_rep_bam_se/outputs"
    "tests/ecliprepmap-1.0.0/02_split_rmrep_bam_se/inputs"
    "tests/ecliprepmap-1.0.0/02_split_rmrep_bam_se/outputs"
    "tests/ecliprepmap-1.0.0/03_dedup_se/inputs"
    "tests/ecliprepmap-1.0.0/03_dedup_se/outputs"
    "tests/ecliprepmap-1.0.0/04_merge_parsed/inputs"
    "tests/ecliprepmap-1.0.0/04_merge_parsed/outputs"
    "tests/ecliprepmap-1.0.0/05_map_repetitive_elements_pe/inputs"
    "tests/ecliprepmap-1.0.0/05_map_repetitive_elements_pe/outputs"
    "tests/ecliprepmap-1.0.0/06_split_rep_sam_pe/inputs"
    "tests/ecliprepmap-1.0.0/06_split_rep_sam_pe/outputs"
    "tests/ecliprepmap-1.0.0/07_split_rmrep_bam_pe/inputs"
    "tests/ecliprepmap-1.0.0/07_split_rmrep_bam_pe/outputs"
    "tests/ecliprepmap-1.0.0/08_dedup_pe/inputs"
    "tests/ecliprepmap-1.0.0/08_dedup_pe/outputs"
    "tests/ecliprepmap-1.0.0/_eric_outputs"
    "tests/ecliprepmap-1.0.0/_eric_scripts"

    # Mini fixture data
    "tests/fixtures/mini/source"
    "tests/fixtures/mini/expected"
    "tests/fixtures/mini/refs"

    # Example data archives and directories
    "examples/data_for_repeat_mapping_hg19.tar.gz"
    "examples/repeat-mapping-hg19-refdata.tar.gz"
    "examples/example_data_for_repeat_mapping_hg38.tar.gz"
    "examples/example_data_for_repeat_mapping_hg38"
    "examples/repeat-mapping-grch38-refdata"
)

# Calculate sizes
echo "Calculating sizes..."
total_size=0
for target in "${TARGETS[@]}"; do
    if [ -e "$target" ]; then
        size=$(du -sb "$target" 2>/dev/null | cut -f1)
        total_size=$((total_size + size))
    fi
done

echo -e "${GREEN}Total size to be deleted: $(numfmt --to=iec-i --suffix=B $total_size)${NC}"
echo ""

# Show what will be deleted
echo "Files/directories to be deleted:"
for target in "${TARGETS[@]}"; do
    if [ -e "$target" ]; then
        size=$(du -sh "$target" 2>/dev/null | cut -f1)
        echo "  - $target ($size)"
    else
        echo "  - $target (already deleted)"
    fi
done
echo ""

if [ "$DRY_RUN" = true ]; then
    echo -e "${YELLOW}DRY RUN MODE - No files will be deleted${NC}"
    echo "Run with --confirm to actually delete files"
    exit 0
fi

# Final confirmation
echo -e "${RED}WARNING: This will permanently delete $(numfmt --to=iec-i --suffix=B $total_size) of data${NC}"
echo ""
read -p "Are you sure you want to continue? Type 'DELETE' to confirm: " -r
echo
if [ "$REPLY" != "DELETE" ]; then
    echo "Aborted"
    exit 1
fi

# Perform deletion
echo "Deleting files..."
deleted_count=0
for target in "${TARGETS[@]}"; do
    if [ -e "$target" ]; then
        echo "  Deleting: $target"
        rm -rf "$target"
        deleted_count=$((deleted_count + 1))
    fi
done

echo ""
echo -e "${GREEN}Deletion complete!${NC}"
echo "  Files/directories deleted: $deleted_count"
echo "  Space recovered: ~$(numfmt --to=iec-i --suffix=B $total_size)"
echo ""

# Verify provenance is intact
echo "Verifying test-provenance archive..."
if [ -d "test-provenance" ] && [ -f "TEST_DOCUMENTATION.md" ]; then
    echo -e "${GREEN}✓ Provenance archive intact${NC}"
    echo "  - test-provenance/ directory: $(du -sh test-provenance | cut -f1)"
    echo "  - TEST_DOCUMENTATION.md: $(du -sh TEST_DOCUMENTATION.md | cut -f1)"
else
    echo -e "${RED}✗ WARNING: Provenance files missing!${NC}"
fi

# Show remaining test files
echo ""
echo "Remaining test configuration files:"
find tests/ -name "*.yaml" -o -name "*.py" | head -10
echo "  ... (see full list with: find tests/ -type f)"

echo ""
echo -e "${GREEN}Done!${NC}"
echo ""
echo "Next steps:"
echo "  1. Verify test-provenance/ is intact: ls -lh test-provenance/"
echo "  2. Push to GitHub if not already done: git push origin codex/python-conversion-python-only-cleanup"
echo "  3. Test reconstruction: See test-provenance/RECONSTRUCTION.md"
echo ""
echo "To restore data, see: test-provenance/BEADS.md"
