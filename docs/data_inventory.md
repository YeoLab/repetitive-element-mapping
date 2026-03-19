# Data Inventory (PPA1 Example + Required References)

## Location and path semantics
- In this workspace, `data/` is a symlink to `/Volumes/X9Pro/Yeo/repetitive-element-pipeline`.
- All paths below are documented using repository-relative paths for reproducibility.

## PPA1 run metadata
- Input YAML:
  - `data/PPA1_rep1.yaml`
  - `data/PPA1_rep1/REPELEMENTMAPPING_PPA1_rep1_INPUT.yaml`
- Workflow run metadata:
  - `data/PPA1_rep1/REPELEMENTMAPPING_PPA1_rep1_OUTPUT.json`
  - `data/PPA1_rep1/REPELEMENTMAPPING_PPA1_rep1_LOG.txt`
  - `data/PPA1_rep1/REPELEMENTMAPPING_PPA1_rep1_WORKFLOW-wf_ecliprepmap_se`
  - `data/PPA1_rep1/REPELEMENTMAPPING_PPA1_rep1_VERSION-0.1.0`

## Required SE reference inputs

| Path | Purpose | Size (bytes) |
|---|---|---:|
| `data/bowtie_reference/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.{1,2,3,4,rev.1,rev.2}.bt2` | Bowtie2 index shards | varied |
| `data/MASTER_filelist.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat` | ENST-to-family/type mapping list | 616168 |
| `data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf` | Gencode annotation | 1199970581 |
| `data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat` | Parsed Gencode table format | 41595212 |
| `data/RepeatMask.bed` | RepeatMasker regions | 180955350 |

Additional provided references (not required for minimal SE CWL execution path):
- `data/ALLRepBase_elements.id_table.FULL`
- `data/genelists.chrM.wchr.txt`
- `data/mirbase.v20.hg19.gff3`

## PPA1 example result artifacts
`data/PPA1_rep1/results/` contains expected compressed outputs from the SE workflow, including:
- Parsed summaries: `*.parsed.gz`, `*.reparsed.tsv.gz`
- Dedup outputs: `*.preRmDup.sam.gz`, `*.rmDup.sam.gz`
- Fold-change tables: `*.nopipes.tsv.gz`, `*.withpipes.tsv.gz`
- mpileup-short files: `*.mpileup.short.gz`

A full expected file fingerprint table is recorded in:
- `docs/ppa1_expected_outputs_manifest.tsv`

## Observed line-level structures (sampled)

### SAM-like rmDup output (`*.rmDup.sam.gz`)
- Contains standard SAM columns plus appended annotations such as `RepFamily`/`UniqueGenomic` and normalized family labels.

### Parsed summary (`*.parsed.gz`)
- Starts with `#READINFO` header lines, e.g. `AllReads`, `UsableReads`, `GenomicReads`, `RepFamilyReads`.

### Reparsed fold-change input (`*.reparsed.nopipes.tsv.gz`)
- Header columns:
  - `element`, `IP_read_num`, `IP_clip_rpr`, `Input_read_num`, `Input_clip_rpr`, `Fold_enrichment`, `Information_content`

## Reproducible manifest generation
Command:
```bash
python3 scripts/generate_ppa1_manifest.py \
  --input-dir data/PPA1_rep1/results \
  --output docs/ppa1_expected_outputs_manifest.tsv
```

This records `relative_path`, `size_bytes`, and `sha256` for each expected output file.
