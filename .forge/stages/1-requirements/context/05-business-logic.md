# 05 - Business Logic

## Step 1: map_repetitive_elements

**Script:** `parse_bowtie2_output_realtime_includemultifamily_PE.pl` (PE) or `_SE.pl` (SE)
**Located at:** `bin/perl/`
**What it does:** Internally spawns bowtie2 in streaming mode, then parses alignments in real-time to assign reads to repeat families. Best-scoring assignment per read pair is kept.

**bowtie2 flags used:** `-q --sensitive -a -p 3 --no-mixed --reorder`

**Positional arguments:**
- arg1: read1 fastq.gz
- arg2 (PE only): read2 fastq.gz
- arg3 (SE: arg2): bowtie2 db path (directory + prefix joined)
- arg4 (SE: arg3): output file name (`<r1_nameroot>.Rep.sam`)
- arg5 (SE: arg4): fileListFile1

**rRNA special handling:** RNA28S, RNA18S, RNA5-8S are assigned to RNA45S precursor family (`rRNA_extra_hash`).

**Multi-family reads:** If a read maps to multiple families equally, family column contains pipe-separated names (e.g., `AluJb|AluSx`).

**Resource requirements (CWL):** 8 cores, 16 GB RAM

## Step 2: splitbam

**Script:** `split_bam_to_subfiles_SEorPE.pl`
**Positional arguments:**
- arg1: SAM/BAM file path
- arg2: `PE` or `SE` (uppercase)

**Output:** 25 `.tmp` files in the current working directory, each named by UMI 2-nt prefix (e.g., `AA.tmp`). Files contain reads whose UMI starts with that prefix.

**UMI location:**
- SE: after `_` in read name (e.g., `...5587_CGCCTTGCCG` → UMI starts with `CG`)
- PE: before `:` in read name (e.g., `AGAAA:SN1001:...` → UMI starts with `AG`)

**Resource requirements:** Not specified in CWL (uses defaults).

## Step 3: deduplicate (scattered x25)

**Script:** `duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl`
**Located at:** `bin/perl/` (softlinked from full name to `duplicate_removal.pl`)
**Resource requirements (CWL):** 32 GB RAM

**Positional arguments:**
1. repFamilySam: `<prefix>.rep.tmp`
2. rmRepSam: `<prefix>.rmrep.tmp`
3. se_or_pe: `PE` or `SE`
4. gencodeGTF
5. gencodeTableBrowser
6. repMaskBedFile
7. fileList1

**Conflict resolution logic:** A read mapping to both unique genome and repeat element is assigned to the unique genome only if the genome alignment score exceeds the repeat score by more than `2 * 2 * 6 = 24` units. Otherwise repeat assignment is kept.

**Perl version sensitivity:** Non-deterministic hash iteration in Perl ≥5.18 can cause different tie-breaking results. Must use Perl 5.10.1 (system perl on TSCC). The script partially mitigates this by sorting hash keys in critical loops.

**Output file naming:** The script writes output files to the current working directory, using input filenames as prefixes. Exact glob patterns:
- `*combined_w_uniquemap.rmDup.sam`
- `*combined_w_uniquemap.prermDup.sam`
- `*.parsed_v2.20201210.txt`
- `*.done`

## Step 4: concatenate

**Command:** `cat <file1> <file2> ... > <output>`

For rmDup: concatenate all 25 `*combined_w_uniquemap.rmDup.sam` files into `<dataset>.<barcode>.rmDup.sam`
For preRmDup: concatenate all 25 `*combined_w_uniquemap.prermDup.sam` files into `<dataset>.<barcode>.preRmDup.sam`

In PE pipeline: for final IP rmDup output, concatenate files from BOTH barcode1 AND barcode2.

## Step 5: gzip

**Command:** `gzip -c <file> > <file>.gz`

Applied to concatenated rmDup.sam and preRmDup.sam.

## Step 6: combine_parsed

**Script:** `merge_multiple_parsed_files.simplified_20191022.pl`
**Positional arguments:**
- arg1: output filename
- arg2+: all input .parsed files

**Important:** In CWL this uses `InitialWorkDirRequirement` to stage input files into the working directory. In Snakemake, ensure the script is run from a directory where it can read all input files (or use absolute paths).

## Step 7: calculate_fold_change

**Script:** `calculate_fold_change_from_parsed_files.py`
**Located at:** `bin/python/calculate_fold_change_from_parsed_files.py`
**Arguments:**
- `--ip_parsed <file>`
- `--input_parsed <file>`
- `--out_file_nopipes <file>`
- `--out_file_withpipes <file>`

**Output logic:** Rows with `|` in family column go to `.withpipes.tsv`; rows without go to `.nopipes.tsv`. Both files also exist as `.withpipes.tsv` (all rows).

## Scatter Logic

The 25 UMI prefix scatter is required because deduplication for the full dataset may exceed 32 GB if run on all reads at once. CWL handles this with `scatter`; Snakemake uses wildcards over the `PREFIXES` list.

**Per the PRD:** Run the full (non-downsampled) dataset through deduplication to profile memory. If ≤32GB, remove scatter (1 rule processes all prefixes sequentially or via a loop inside the rule). If >32GB, keep scatter.

## UMI Prefix List

```python
PREFIXES = [
    "AA","AC","AG","AT","AN",
    "CA","CC","CG","CT","CN",
    "GA","GC","GG","GT","GN",
    "TA","TC","TG","TT","TN",
    "NA","NC","NG","NT","NN"
]
```
