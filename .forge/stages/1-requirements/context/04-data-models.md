# Data Models

## Input Source Files (per assembly)

### Downloaded Source Files

| File | hg38 | mm10 | mm39 | Notes |
|------|------|------|------|-------|
| `gencode.{version}.annotation.gtf[.gz]` | v33, uncompressed symlink | VM23, gzipped | VM38, gzipped | GTF format; skip `#` comment lines |
| `{assembly}.fasta` | symlink to GRCh38_no_alt... | absent (needs sourcing) | absent (needs sourcing) | Reference genome for pybedtools getfasta |
| `{assembly}.repeatmasker.tsv.gz` | present | present | present | GTF-like 9-col format (see below) |
| `{assembly}.simplerepeats.tsv.gz` | present | present | present | GTF-like 9-col format |
| `{assembly}.trna.tsv.gz` | present | present | absent for mm39 | GTF-like 9-col format; optional |
| `{species}.gff3` | hsa.gff3 (miRBase v22) | mmu.gff3 (miRBase v22) | absent for mm39 | miRNA coordinates; optional |
| `NR_046233.2.fasta` | absent | present | present | Custom RefSeq rRNA FASTA |

### TSV Source Column Format (repeatmasker / trna / simplerepeats)

All three source files share the same 9-column GTF-like format:
```
col1: chrom
col2: source (e.g. hg38_rmsk, hg38_tRNAs, hg38_simpleRepeat)
col3: feature (always "exon")
col4: start (1-based)
col5: end (1-based, inclusive)
col6: score (float)
col7: strand (+, -, .)
col8: frame (always .)
col9: attributes  → gene_id "NAME"; transcript_id "NAME";
```
Column 6 in the generate_refdata.md context refers to the gene_id in the attribute field (extracted from col9), used as the element name.

---

## Output Reference Files

### 1. parsed_ucsc_tableformat

**Filename pattern:** `gencode.{version}.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat`

**Format:** Tab-separated, 11 columns, first line is header:
```
#ENSG   name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds
```

**Derivation rules:**
- One row per unique `transcript_id` where GTF feature == "transcript"
- `#ENSG`: gene_id from GTF attributes
- `name`: transcript_id from GTF attributes
- `chrom`: GTF col1
- `strand`: GTF col7
- `txStart`: GTF col4 minus 1 (convert to 0-based)
- `txEnd`: GTF col5 (0-based, exclusive)
- `cdsStart`: set to start of first coding exon (or txStart if non-coding)
- `cdsEnd`: set to end of last coding exon (or txEnd if non-coding)
- `exonCount`: count of "exon" feature rows for this transcript
- `exonStarts`: comma-delimited 0-based start positions of exons, trailing comma
- `exonEnds`: comma-delimited 0-based end positions of exons, trailing comma

**Validation:** Reproducing from `gencode.v33.chr_patch_hapl_scaff.annotation.gtf` must produce file identical to existing `gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` (249,044 rows incl. header).

### 2. bowtie2_index/

**Directory contents:**
```
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.1.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.2.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.3.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.4.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.rev.1.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.rev.2.bt2
MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa   ← combined source FASTA
```

**FASTA sequence sources (hg38 reference composition):**
- Gencode transcript sequences (ENST... IDs) — extracted via pybedtools getfasta from genome FASTA
- RepBase / repeatmasker-derived repeat element sequences (AluY, L1P5, MIRB, etc.)
- tRNA sequences (nm-tRNA..., tRNA-... IDs)
- Simple repeat sequences (trf, trf_dup1, etc.)
- RefSeq rRNA sequences (NR_046235.3-18S, NR_046235.3-28S, NR_046235.3-45S, etc.)
- Custom sequences from NR_046233.2.fasta (for mm10/mm39)

**Missing ID threshold:** If >1% of expected IDs cannot be retrieved from source FASTA, produce a `missing_ids_report.txt` and continue. Do not abort.

### 3. UniqueGenomicElements.{assembly}.bed

**Format:** 6-column BED (standard):
```
chrom   start   end     name    score   strand
```

**Derivation sources (all optional except repeatmasker):**
- repeatmasker (gene_id from col9 attributes) → positions + proximal flanks
- trna (gene_id from col9 attributes) → positions + proximal flanks (optional)
- simplerepeats (transcript_id from col9 attributes) → positions + proximal flanks
- parsed_ucsc_tableformat (strand always "-") → positions + proximal flanks
- gff3 "Name" attribute (miRNA IDs like MI0000112) → positions + proximal flanks

**500bp proximal flanks:** For each element at `[start, end)`, add:
- Upstream proximal: `[start - 500, start)`
- Downstream proximal: `[end, end + 500)`
Proximal entry name appended with `-proximal`.

**Reference size:** hg38 UniqueGenomicElements.hg38.bed has 5,618,483 lines.

### 4. MASTER_FILELIST.{date}.*.tsv (and .list)

**Format:** Tab-separated, 5 columns, no header:
```
col1: sequence_id (ENST..., repeat_name, tRNA_name, miRNA_ID, NR_...)
col2: gene_id / ENSG (or same as col1 for non-Gencode entries)
col3: short name / symbol
col4: family / repeat class
col5: genelist label (e.g. genelists.RNU1, genelists.SINE, genelists.tRNA)
```

**Reference size:** hg38 MASTER_FILELIST has 26,422 lines.

**Note:** The `.list` and `.tsv` files in hg38 have identical content (same 5 columns, same rows). The `.list` filename is used by the Perl parse scripts via ARGV[3].

---

## Assembly-Specific Notes

| Assembly | GTF version | Gencode release | trna | gff3 (miRNA) | Custom FASTA |
|----------|-------------|-----------------|------|--------------|--------------|
| hg38 | v33 | GRCh38.p13 | hg38.trna.tsv.gz | hsa.gff3 | none |
| mm10 | VM23 | GRCm38.p6 | mm10.trna.tsv.gz | mmu.gff3 | NR_046233.2.fasta |
| mm39 | VM38 | GRCm39 | absent | absent | NR_046233.2.fasta |
