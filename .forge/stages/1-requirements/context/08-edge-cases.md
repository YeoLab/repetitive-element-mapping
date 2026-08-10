# 08 - Edge Cases

## EC-01: Empty UMI prefix bin

**Scenario:** A `.tmp` file for a given prefix (e.g., `NN.tmp`) contains no reads.
**CWL behavior:** The deduplicate script receives an empty file and produces empty output files (including an empty `.parsed_v2` file).
**Snakemake requirement:** The rule must still create all output files (including empty ones) to satisfy Snakemake's output expectations. The concatenate and combine steps must handle empty inputs gracefully.

## EC-02: Perl version sensitivity

**Scenario:** Using Perl ≥5.18 introduces non-deterministic hash iteration, causing different tie-breaking in deduplication.
**Requirement:** Always invoke Perl scripts with `/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl`. Do not use the system default `perl` command (which may point to a newer version).
**Verification:** Check `perl --version` shows 5.10.1 in the environment used by the rule.

## EC-03: BAM reads not in FASTQ

**Scenario:** The rmRep BAM contains read names that do not appear in the FASTQ. This is an upstream data preparation error.
**Requirement:** Document that this is a user error; add a validation step or warning. The pipeline should fail with a clear error message if the BAM contains reads not in the FASTQ.

## EC-04: Missing barcode2 in PE mode

**Scenario:** User sets `se_or_pe: PE` but omits barcode2 fields in config.
**Requirement:** Snakefile validates config at startup and raises a descriptive error before any rules run.

## EC-05: PE fields present in SE mode

**Scenario:** User sets `se_or_pe: SE` but includes barcode2 fields.
**Requirement:** Snakefile raises a validation error.

## EC-06: split_bam_to_subfiles_SEorPE.pl writes .tmp files to cwd

**Scenario:** The Perl script writes output .tmp files to the current working directory (not the directory of the input file).
**Requirement:** Each Snakemake rule invocation for splitbam must use `params.cwd` or `shadow:` to ensure the .tmp files land in the correct location and do not collide between concurrent barcode runs.

## EC-07: merge_multiple_parsed_files.pl requires files in cwd

**Scenario:** The Perl script (combine_parsed) may expect input files in the current working directory due to CWL `InitialWorkDirRequirement`.
**Requirement:** Either run the script from the output directory (using `shell: "cd {params.dir} && ..."`) or pass absolute paths. Verify which mode the script uses.

## EC-08: Bowtie2 version differences

**Scenario:** CWL originally used bowtie2/2.2.6 via module load; conda env may install a newer version.
**Requirement:** Test that outputs match between versions. If they differ, pin bowtie2 to 2.2.6 in the conda env.

## EC-09: Very large datasets

**Scenario:** Full (non-downsampled) PE datasets may produce deduplication jobs exceeding 32GB.
**Requirement:** Profile memory usage on full dataset before deciding whether to keep or remove scatter. If scatter is kept, each SLURM job is independently submitted and can be monitored.

## EC-10: .tmp file naming collisions

**Scenario:** Both splitbam_repsam and splitbam_rmrepbam produce files named `AA.tmp`, `AC.tmp`, etc. If they run in the same directory, they will overwrite each other.
**Requirement:** Architect must design separate working directories (e.g., `<barcode>/rep/` and `<barcode>/rmrep/`) or use distinct naming (e.g., `AA.rep.tmp` and `AA.rmrep.tmp`).
