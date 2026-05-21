# User Flows

## Flow 1: Generate All References for a New Assembly (mm10 or mm39)

**Preconditions:** Source files exist in `examples/inputs/{assembly}/downloaded/`

```
1. Decompress GTF if gzipped (.gtf.gz)
2. Run generate_parsed_ucsc_tableformat.py → parsed_ucsc_tableformat
3. Run generate_bowtie2_index.py → bowtie2_index/ directory + .fa FASTA
4. Run generate_unique_genomic_elements.py → UniqueGenomicElements.{assembly}.bed
5. Run generate_master_filelist.py → MASTER_FILELIST.{date}.*.tsv
6. Verify output counts (line counts, non-empty files)
7. Run pipeline dry-run: snakemake -n with new references
```

**Output location:** `examples/inputs/{assembly}/`

## Flow 2: Reproduce hg38 References for Validation

**Preconditions:** hg38 source files exist in `examples/inputs/hg38/downloaded/`

```
1. Run each script against hg38 sources
2. Diff generated files against existing hg38 reference files
3. Compute similarity metric: lines matching / total lines ≥ 99%
4. Report any discrepancies
```

**This flow is mandatory before producing mm10/mm39 references.**

## Flow 3: Handle Optional Inputs Missing

**Example:** mm39 has no `mm39.trna.tsv.gz` and no `mm39.gff3`

```
1. User runs generate_unique_genomic_elements.py without --trna and --gff3
2. Script logs: "WARNING: --trna not provided; tRNA elements will be omitted"
3. Script logs: "WARNING: --gff3 not provided; miRNA proximal regions will be omitted"
4. Script produces UniqueGenomicElements.mm39.bed without tRNA/miRNA entries
5. Script logs final element count per source category
```

**Important:** mm10 has `mm10.trna.tsv.gz` and `mmu.gff3`. mm39 has neither.

## Flow 4: Handle Missing IDs During Bowtie2 Index Generation

```
1. Script identifies all FASTA header IDs needed from GTF + repeat sources
2. Script attempts to retrieve each ID from the reference FASTA
3. If missing IDs > 1% of total:
   a. Print WARNING to stderr
   b. Write missing_ids_report.txt listing all missing IDs + suggested sources
   c. Proceed with remaining IDs (do NOT abort)
4. User reviews report, decides whether to provide additional custom FASTA sequences
```

## Flow 5: Perl Script Compatibility Check

```
1. Run CWL or Snakemake pipeline (or dry-run) with new mm10/mm39 references
2. Check parse_bowtie2_output Perl scripts accept new MASTER_FILELIST format
3. Check duplicate_removal Perl script accepts new references
4. If hardcoded hg38 values found: modify only those specific lines
5. Document any Perl modifications in NOTES or changelog
```

## Flow 6: Custom FASTA Inclusion

**Example:** NR_046233.2 exists in RefSeq but not in the genome FASTA

```
1. User provides --custom-fasta NR_046233.2.fasta to generate_bowtie2_index.py
2. Script appends custom sequences to the combined FASTA before indexing
3. Script logs: "Added 1 custom FASTA sequence(s): NR_046233.2"
4. Custom sequence is tracked in MASTER_FILELIST with appropriate family label
```
