# Output Format Specification (Observed + Code-Aligned)

## Source basis
- Observed files under `data/PPA1_rep1/results/*.gz`.
- CWL tool wiring and Perl/Python scripts.

## 1) `.parsed` and `.reparsed.tsv`-style summary files

### Header lines
Common readinfo-like rows begin with `#READINFO`.
Observed keys:
- `AllReads` (in `.parsed`, not always present in reparsed files)
- `UsableReads`
- `GenomicReads`
- `RepFamilyReads`

### Body rows
Two row classes:
- `TOTAL\t<family>\t<read_count>\t<fraction_or_rpm>`
- `ELEMENT\t<primary_label>\t<read_count>\t<fraction_or_rpm>\t<element_bundle>\t<gene_bundle>`

Notes:
- Some files labeled `.reparsed.tsv.gz` retain `TOTAL` rows but may omit initial `AllReads` line.
- Floating precision differs by script/formatting (`sprintf` in Perl; pandas/float formatting in Python).

## 2) `.nopipes.tsv` and `.withpipes.tsv`
Tab-separated with header:
- `element`
- `IP_read_num`
- `IP_clip_rpr`
- `Input_read_num`
- `Input_clip_rpr`
- `Fold_enrichment`
- `Information_content`

Meaning:
- `clip_rpr` terms are read proportion metrics.
- `Fold_enrichment = IP_clip_rpr / Input_clip_rpr` (with pseudocount handling in upstream script).
- `Information_content = IP_clip_rpr * log2(IP_clip_rpr / Input_clip_rpr)`.
- `.nopipes.tsv` excludes ambiguous multi-family mappings (contains no `|` in element family labels).
- `.withpipes.tsv` includes both ambiguous and unambiguous families.

## 3) SAM-like pre-rmdup and rmdup outputs
Examples:
- `*.preRmDup.sam.gz`
- `*.rmDup.sam.gz`

Characteristics:
- Tab-delimited SAM fields preserved.
- Additional right-side annotations are appended beyond canonical SAM tags.
- `ZZ:Z:` tags may include pipe-delimited ENST sets.
- rmdup files append classification fields such as:
  - `RepFamily` / `UniqueGenomic`
  - normalized family/element label string

## 4) mpileup-short outputs
Examples:
- `*.rmDup.sam.gz.tmp.RNA18S.bam.sorted.bam.mpileup.short.gz`
- `*.rmDup.sam.gz.tmp.RNA28S.bam.sorted.bam.mpileup.short.gz`

Observed columns:
1. transcript/accession id
2. 1-based position
3. reference base
4. depth/count

## 5) Determinism/tolerance for validation
Validation should be tolerant for:
- floating point formatting differences
- row ordering differences where hash/dict iteration affects output order

Validation should be strict for:
- schema/column set
- primary key identity (element/family labels)
- integer read counts
