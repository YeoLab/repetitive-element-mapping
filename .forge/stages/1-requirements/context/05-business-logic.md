# Business Logic

## Step 1: generate_parsed_ucsc_tableformat.py

### Rules
- Read input GTF; skip lines starting with `#`
- Process only rows where column 3 (feature) == `"transcript"` to build the table row
- Process only rows where feature == `"exon"` to count and collect exon coordinates
- Group exons by transcript_id; sort by start position
- Coordinates in GTF are 1-based inclusive → convert to 0-based half-open (subtract 1 from start only)
- `exonStarts` and `exonEnds` are comma-delimited lists with a trailing comma (matches hg38 reference format)
- `cdsStart` / `cdsEnd`: derive from UTR/CDS features if available; fall back to txStart/txEnd for non-coding transcripts
- Output sorted by `#ENSG` (gene_id) then by `name` (transcript_id) — verify exact sort from hg38 reference

### Validation Requirement
Run against `gencode.v33.chr_patch_hapl_scaff.annotation.gtf` → output must match existing `gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` with 100% line-for-line identity (this is the ground truth, not the 99% threshold).

---

## Step 2: generate_bowtie2_index.py

### Sequence Assembly Logic
1. Extract Gencode transcript sequences:
   - Parse parsed_ucsc_tableformat for all transcript IDs, chrom, strand, exon coordinates
   - Use pybedtools getfasta with the assembly FASTA to extract per-transcript sequence
   - Name FASTA headers by transcript_id (e.g. `>ENST00000383925.1`)
2. Extract repeat element sequences:
   - From repeatmasker TSV: extract unique element positions per gene_id
   - From simplerepeats TSV: extract positions per transcript_id
   - Use pybedtools getfasta to extract sequences; name headers by element name
3. Extract tRNA sequences (if trna TSV provided):
   - From trna TSV: extract positions per transcript_id
   - Use pybedtools getfasta
4. Include rRNA sequences (RefSeq NR_ entries):
   - The NR_ sequences in hg38 are named with suffixes: NR_046235.3-18S, NR_046235.3-28S, NR_046235.3-45S
   - These correspond to the rRNA_extra_hash in the Perl parse scripts (RNA18S, RNA28S, RNA45S)
5. Append custom FASTA sequences (e.g. NR_046233.2.fasta) verbatim

### Bowtie2 Index Build
- Concatenate all FASTA sequences into a single `.fa` file
- Run `bowtie2-build` on the concatenated FASTA
- Output index prefix must match the `.fa` filename base (without `.fa`)

### Missing ID Handling
- Count IDs that fail pybedtools getfasta retrieval
- If missing > 1% of total: write `missing_ids_report.txt`, log WARNING, continue
- If missing ≤ 1%: log count, continue silently

---

## Step 3: generate_unique_genomic_elements.py

### Element Inclusion Logic
1. **repeatmasker entries**: chrom, start-1, end, gene_id, score, strand → then add ±500 bp proximal pairs
2. **trna entries** (optional): same extraction from trna TSV
3. **simplerepeats entries**: chrom, start-1, end, transcript_id, score, strand → then add proximal pairs
4. **parsed_ucsc_tableformat entries**: chrom, txStart, txEnd, name, 0, "-" (strand always "-") → proximal pairs
5. **gff3 entries** (optional): extract Name attribute, chrom, start-1, end, ., strand → proximal pairs

### Proximal Flanking Rule (500 bp)
For element at `[s, e)` on chromosome:
- Upstream proximal: `[max(0, s-500), s)` named `{element_id}-proximal`
- Downstream proximal: `[e, e+500)` named `{element_id}-proximal`

### Output Sort
The BED output must be sorted — verify sort order from hg38 reference before applying to mouse.

---

## Step 4: generate_master_filelist.py

### Row Construction Logic
The 5-column format:
1. `sequence_id`: ENST transcript IDs (from parsed_ucsc_tableformat), repeat names (from repeatmasker gene_id), tRNA names, miRNA IDs (from gff3 Name attribute), NR_ custom entries
2. `gene_id`: ENSG for Gencode entries; same as sequence_id for non-Gencode entries
3. `short_name`: gene name for Gencode (from GTF gene_name attribute); element name for repeats
4. `family`: repeat class/family (from repeatmasker); "tRNA" for tRNA entries; "miRNA" for miRNA; "rRNA" for rRNA
5. `genelist_label`: `genelists.{FAMILY}` — must match the labeling convention used in existing hg38 file

### Ordering
Must produce rows in same order as hg38 reference (verify grouping order: first Gencode transcripts, then repeat elements, then tRNAs, etc.). Exact order determines compatibility with Perl scripts.

---

## Step 5: Perl Script Compatibility

### Scripts to Check
1. `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl` — reads MASTER_FILELIST via `ARGV[3]`; the hardcoded hg38 path is commented out; active code accepts any path via argument → **no modification needed**
2. `bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl` — same pattern → **no modification needed**
3. `bin/perl/duplicate_removal_inline_paired...pl` — check for any hardcoded assembly-specific paths or chromosome name assumptions
4. `bin/perl/split_bam_to_subfiles_SEorPE.pl` — comment mentions hg38 but actual logic appears assembly-agnostic
5. `bin/perl/RepElement_pipeline_1dataset.pl` — has `my $species = "hg38"` hardcoded on line 4; this is an **older orchestration script**, not called by CWL or Snakemake dropin — assess whether modification is needed based on usage trace

### Modification Threshold
Modify only if a hardcoded value causes incorrect behavior when mm10/mm39 references are provided as CLI arguments. Commented-out code is not a blocker.

### Chromosome Naming
Mouse assemblies use `chr1`-`chrX` format identical to hg38 conventions in Gencode annotations → no chromosome naming conflicts expected.
