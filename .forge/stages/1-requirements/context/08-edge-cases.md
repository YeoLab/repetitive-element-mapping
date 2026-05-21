# Edge Cases

## Edge Case 1: mm39 Missing tRNA and miRNA GFF3

**Situation:** `mm39/downloaded/` has no `mm39.trna.tsv.gz` and no `mm39.gff3`.

**Expected behavior:**
- `generate_unique_genomic_elements.py`: when `--trna` is absent, log WARNING and skip tRNA entries entirely. When `--gff3` is absent, log WARNING and skip miRNA proximal entries. Do not raise an exception; produce a valid (though smaller) BED file.
- `generate_master_filelist.py`: when tRNA and miRNA inputs are absent, produce TSV without those entry categories.
- The pipeline must still run correctly on a MASTER_FILELIST and UniqueGenomicElements that lack tRNA/miRNA entries — these element categories will simply have zero reads mapping to them.

## Edge Case 2: Missing IDs During FASTA Extraction

**Situation:** Some transcript IDs in the GTF annotation cannot be extracted from the genome FASTA (e.g., due to chromosome naming mismatch, scaffold-only entries, or alternate loci).

**Expected behavior:**
- Count missing IDs after pybedtools getfasta attempt
- If missing ≤ 1% of total: log count to stderr, continue without reporting file
- If missing > 1% of total: write `missing_ids_report.txt` in output directory, log WARNING, continue without those IDs
- Never silently drop entries without any logging

## Edge Case 3: GTF in .gz Format

**Situation:** mm10 and mm39 GTFs are gzipped (`gencode.VM23.annotation.gtf.gz`, `gencode.VM38.annotation.gtf.gz`); hg38 GTF is a symlink to an uncompressed file.

**Expected behavior:** `generate_parsed_ucsc_tableformat.py` must accept both `.gtf` and `.gtf.gz` inputs, detecting compression from the file extension.

## Edge Case 4: Duplicate Element Names in RepeatMasker/SimpleRepeats

**Situation:** Multiple genomic intervals may share the same `gene_id` in the RepeatMasker TSV (e.g., many AluY loci all named "AluY").

**Expected behavior for bowtie2 index:** Each unique locus gets its own FASTA entry. Deduplicate by (chrom, start, end, name) to avoid redundant sequences. Name disambiguation strategy: `AluY`, `AluY_dup1`, `AluY_dup2`, etc. — or whatever strategy reproduces the hg38 reference header format exactly.

**Expected behavior for MASTER_FILELIST:** Each unique sequence_id appears only once. If repeat elements appear once per family (not per locus), confirm this against hg38 reference.

## Edge Case 5: NR_046233.2 rRNA Custom FASTA

**Situation:** mm10 and mm39 include `NR_046233.2.fasta` as a custom sequence (mouse rRNA). This corresponds to the rRNA_extra_hash entries in the Perl parse scripts.

**Expected behavior:** `generate_bowtie2_index.py` must accept `--custom-fasta NR_046233.2.fasta` and append these sequences with headers matching the naming convention used for rRNA entries in the MASTER_FILELIST (e.g., `NR_046233.2-18S`, `NR_046233.2-28S`, `NR_046233.2-45S`).

**Validation:** The MASTER_FILELIST for mm10/mm39 must include NR_046233.2 entries with the correct family labels (`rRNA` or equivalent).

## Edge Case 6: Coordinate System Differences Between GTF and BED

**Situation:** GTF uses 1-based inclusive coordinates; BED uses 0-based half-open.

**Expected behavior:** All scripts must convert GTF coordinates to 0-based internally:
- `start_0based = gtf_start - 1`
- `end_0based = gtf_end` (no change)

Off-by-one errors here would break the 99% similarity threshold against hg38.

## Edge Case 7: Perl Script RepElement_pipeline_1dataset.pl

**Situation:** Line 4 has `my $species = "hg38"` hardcoded. Lines 95-96 branch on `$species eq "hg38"` to set bowtie_db path.

**Assessment:** This script is a standalone job orchestrator that hard-codes paths. It is NOT called by the CWL workflow or Snakemake dropin; it exists in `bin/perl/` as a legacy convenience script. The CWL and Snakemake paths pass all reference file paths via YAML/config arguments. Modifying this script is NOT required unless the user explicitly runs it directly.

**Expected behavior:** Document the hardcoded species in NOTES. Do not modify unless the user requests standalone RepElement_pipeline_1dataset.pl support for mm10/mm39.

## Edge Case 8: Mouse Genome FASTA Availability

**Situation:** `mm10/downloaded/` and `mm39/downloaded/` do not currently contain genome FASTA files (only hg38 has `hg38.fasta` as a symlink). The genome FASTA is needed for pybedtools getfasta to extract transcript sequences.

**Expected behavior:** The `generate_bowtie2_index.py` script requires `--fasta` as a mandatory argument. The script must verify the FASTA exists and is readable before proceeding. If the mouse genome FASTA is not yet present, the script must fail with a clear error: "Genome FASTA not found at {path}. Please provide a reference genome FASTA for {assembly}."

**Resolution:** The architect must determine whether sourcing the mouse genome FASTA is in-scope, or whether the user will provide it. The architect-prompt captures this as a critical open gap.

## Edge Case 9: Simple Repeat Naming Disambiguation

**Situation:** The hg38 simplerepeats TSV has `gene_id "trf"` for all entries, with `transcript_id "trf"`, `"trf_dup1"`, `"trf_dup2"`, etc. for disambiguation.

**Expected behavior:** The script must use `transcript_id` (not `gene_id`) from simplerepeats when naming simple repeat entries, to match the hg38 convention.
