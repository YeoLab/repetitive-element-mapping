# Vision and Goals

## Project Vision

Extend the eCLIP repetitive element mapping pipeline (ecliprepmap) to support two additional mouse genome assemblies: mm10 (GRCm38) and mm39 (GRCm39). The pipeline currently ships with reference data only for hg38 (GRCh38). This work produces scripts that reproducibly generate all four required reference files for any supported assembly.

## Primary Goal

Generate the following four reference artifacts for mm10 and mm39 such that they are structurally and format-identical to the existing hg38 references:

1. `gencode.{version}.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` — tab-separated transcript table derived from a Gencode GTF
2. `bowtie2_index/` — Bowtie2 index built from a curated multi-source FASTA
3. `UniqueGenomicElements.{assembly}.bed` — 6-column BED of repeat/tRNA/miRNA regions plus 500 bp proximal flanks
4. `MASTER_FILELIST.{date}.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv` — 5-column TSV mapping sequence IDs to repeat families

## Secondary Goal

Produce reusable Python scripts (one per reference type) that accept assembly-generic inputs so that future assemblies can be added with minimal effort.

## Out of Scope

- Changes to the CWL or Snakemake pipeline execution logic
- Changes to Perl scripts except where hardcoded values prevent acceptance of new assemblies
- Adding new pipeline features (new output formats, new metrics, etc.)
- Reference generation for assemblies other than mm10 and mm39

## Success Criteria

1. All four reference files generated for mm10 and mm39.
2. All four hg38 reference files faithfully reproduced by the new scripts with ≥99% similarity to existing hg38 references.
3. Format of generated files is identical to hg38 references (same column count, same header lines, same sort order).
4. Perl scripts accept mm10/mm39 references without modification, or any required modifications are minimal and documented.
5. Pipeline dry-run (`snakemake -n`) completes without error when pointed at mm10/mm39 references.
