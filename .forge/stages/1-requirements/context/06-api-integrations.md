# API Integrations

## External Tools Required

### Bowtie2
- **Version:** ≥2.5 (available via `workflow/envs/dropin.yaml` conda env; also available via `module load ecliprepmap/1.0.0`)
- **Usage:** `bowtie2-build <fasta> <prefix>` to generate the index
- **Called by:** `generate_bowtie2_index.py` (new script)
- **No network access required**

### pybedtools
- **Version:** compatible with Python 3.x (already in dropin.yaml env)
- **Usage:** `pybedtools.BedTool.sequence(fi=fasta, s=True)` for strand-aware FASTA extraction
- **Called by:** `generate_bowtie2_index.py` (for getfasta) and `generate_unique_genomic_elements.py` (for coordinate operations)
- **Depends on:** bedtools binary in PATH

### samtools
- **Version:** ≥1.17 (already in dropin.yaml)
- **Usage:** Only if BAM manipulation is needed (not expected for reference generation)

## Module Environment

The `ecliprepmap/1.0.0` module provides the correct PATH for Perl, Python, and Bowtie2 on TSCC:
```bash
module load ecliprepmap/1.0.0
```

The dropin conda environment (`workflow/envs/dropin.yaml`) can also be used:
```bash
conda activate <dropin-env>
```

## No External API Calls

This task does not require network access. All source files are pre-downloaded:
- Gencode GTFs: already in `examples/inputs/{assembly}/downloaded/`
- RepeatMasker/SimpleRepeats/tRNA TSVs: already downloaded
- miRNA GFF3s: already downloaded
- NR_046233.2 custom FASTA: already provided in mm10/mm39 downloaded directories

**No calls to UCSC, Ensembl, NCBI, or miRBase APIs are needed.**

## File Format Dependencies

### pybedtools getfasta
- Input: BED3+ file with coordinates, reference FASTA (must be indexed with `.fai`)
- Output: FASTA sequences keyed by `chrom:start-end` or custom name
- Requires: `samtools faidx` or `bedtools faidx` pre-run on the reference FASTA

### bowtie2-build
- Input: FASTA file
- Output: 6-file index (`.1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`)
- Runtime estimate: ~10-30 minutes for a reference of this size
