# 06 - API Integrations

## External Tools

### bowtie2
- Called internally by `parse_bowtie2_output_realtime_includemultifamily_PE/SE.pl`
- Not invoked directly in Snakemake rules (the Perl script handles the subprocess)
- Must be available in PATH: `module load bowtie2/2.2.6` or via conda env
- Version used in CWL context: bowtie2 ≥2.5 (per CLAUDE.md dropin env), but historical CWL runs used 2.2.6

### samtools
- Used by `split_bam_to_subfiles_SEorPE.pl` to read BAM files
- Must be in PATH: samtools ≥1.17 (per CLAUDE.md dropin env)

### Perl scripts (in bin/perl/)

| Script | Invocation |
|--------|-----------|
| `parse_bowtie2_output_realtime_includemultifamily_PE.pl` | Positional args: r1, r2, bowtie2_db_path, output, fileListFile1 |
| `parse_bowtie2_output_realtime_includemultifamily_SE.pl` | Positional args: r1, bowtie2_db_path, output, fileListFile1 |
| `split_bam_to_subfiles_SEorPE.pl` | Positional args: sam/bam, SE_or_PE |
| `duplicate_removal_inline_paired...pl` (softlinked as `duplicate_removal.pl`) | Positional args: repFamilySam, rmRepSam, SE_or_PE, gencodeGTF, gencodeTableBrowser, repMaskBedFile, fileList1 |
| `merge_multiple_parsed_files.simplified_20191022.pl` | Positional args: output_file, input_files... |

All Perl scripts must be invoked with:
```bash
/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl <script_path> <args>
```

### Python scripts (in bin/python/)

| Script | Invocation |
|--------|-----------|
| `calculate_fold_change_from_parsed_files.py` | `--ip_parsed`, `--input_parsed`, `--out_file_nopipes`, `--out_file_withpipes` |

Python script must be invoked with the ecliprepmap conda environment python:
```bash
/tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/ecliprepmap-0.1.0/bin/python
```
OR via the Snakemake conda env that includes compatible numpy/pandas.

## Conda Environment

The existing `workflow/envs/dropin.yaml` (referenced in CLAUDE.md) specifies:
- Python 3.11
- bowtie2 ≥2.5
- samtools ≥1.17
- numpy
- pandas

Note: The original CWL used bowtie2/2.2.6 via `module load`. The Snakemake workflow should use the conda env for portability, but may need to confirm version compatibility with the Perl scripts.

## No Network Integrations

This is a batch HPC pipeline. No web services, APIs, or databases are accessed at runtime.
