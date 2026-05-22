## Pre-Execution Plan: 3-implement-perl2py

### Three Most Likely Failure Modes

1. **Dedup hash-ordering divergence** — The Perl 5.10 dedup script iterates sorted
   hash keys (Perl 5.18+ mitigation already in place). If Python sorts UMI keys or
   read-name keys in a different collation order, tie-breaking in deduplication will
   differ. Signal: >5% deviation in per-element counts in `.parsed` files.
   Watch for: any dict iteration over UMI barcodes or read IDs in the dedup script.

2. **Bowtie2 subprocess pipe truncation** — The Perl map scripts open bowtie2 via
   a pipe (`open(BOWTIE2, "bowtie2 ... |")`). Python's subprocess.Popen must drain
   stdout fully before close; otherwise EOF arrives early, truncating the SAM output.
   Signal: rep.sam shorter than Perl-produced reference or missing terminal reads.
   Watch for: any SIGPIPE errors or incomplete output on large input files.

3. **Absolute-path CWD coupling** — Both split_bam and dedup scripts cd into a
   working directory and write outputs relative to CWD. Python equivalents must
   replicate this exactly (output filenames, not paths, written to CWD).
   Signal: FileNotFoundError in downstream Snakemake rules.

### Verification Steps (after each script)

For each translated script:
  1. Run Snakemake SE small (`examples/repeat_mapping_SE_small.yaml`) with
     `--forcerun <rule>` to regenerate only that rule's outputs.
  2. Diff Python-generated output vs Perl-generated reference (from previous run
     or from test-provenance/).
  3. Threshold: ≤5% deviation in numeric counts; zero deviation in headers/structure.

### Context Dependencies

Files to read before starting:
  - bin/perl/split_bam_to_subfiles_SEorPE.pl
  - bin/perl/merge_multiple_parsed_files.simplified_20191022.pl
  - bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl
  - bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl
  - bin/perl/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl
  - workflow/rules/common.smk  (rules to update)
  - workflow/rules/SE.smk      (map_rep_se rule)
  - workflow/rules/PE.smk      (map_rep_pe, merge_ip_parsed rules)
  - Snakefile                  (PYTHON_ECLIP, PERL constants, shell.prefix)
  - examples/repeat_mapping_SE_small.yaml  (test config)
  - workflow/envs/dropin.yaml  (available packages)

### Implementation Order

Translate simplest-first to fail fast on easy scripts:
  1. split_bam_to_subfiles_SEorPE.pl  → workflow/scripts/split_bam_to_subfiles.py
  2. merge_multiple_parsed_files.simplified_20191022.pl → workflow/scripts/merge_parsed_files.py
  3. parse_bowtie2_output_realtime_includemultifamily_SE.pl → workflow/scripts/map_repetitive_elements_se.py
  4. parse_bowtie2_output_realtime_includemultifamily_PE.pl → workflow/scripts/map_repetitive_elements_pe.py
  5. duplicate_removal_inline_paired...pl → workflow/scripts/deduplicate.py

After each: update the corresponding Snakemake rule to call Python.

### Acceptance Criteria (from translate_perl.md)

- Python+Snakemake outputs match Perl+Snakemake outputs for SE small dataset
- Python+Snakemake outputs match Perl+Snakemake outputs for PE small dataset
- .nopipes.tsv and .withpipes.tsv match test-provenance/ reference files
- Tolerance: ≥95% row-level match on numeric columns; 100% structure match
