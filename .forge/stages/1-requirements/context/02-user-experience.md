# User Experience

## Primary User

A bioinformatician (Yeo Lab member) who:
- Runs the eCLIP repetitive element mapping pipeline on TSCC (SLURM)
- Wants to analyze mouse eCLIP data (mm10 or mm39) using the same pipeline as human data
- Is comfortable running Python scripts and Snakemake on the command line
- Expects reference generation to be a one-time setup task, not repeated per experiment

## User Interaction Model

The user runs Python scripts to generate reference files. They do not interact with a GUI or web interface.

### Invocation Pattern (per script)

```bash
python generate_parsed_ucsc_tableformat.py \
    --gtf gencode.VM38.annotation.gtf.gz \
    --output gencode.VM38.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat

python generate_bowtie2_index.py \
    --gtf gencode.VM38.annotation.gtf.gz \
    --repeatmasker mm39.repeatmasker.tsv.gz \
    --simplerepeats mm39.simplerepeats.tsv.gz \
    --trna mm39.trna.tsv.gz \
    --fasta mm39.fasta \
    --custom-fasta NR_046233.2.fasta \
    --output-dir bowtie2_index/ \
    --output-prefix MASTER_FILELIST.{date}.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat

python generate_unique_genomic_elements.py \
    --repeatmasker mm39.repeatmasker.tsv.gz \
    --trna mm39.trna.tsv.gz \            # optional
    --simplerepeats mm39.simplerepeats.tsv.gz \
    --parsed-ucsc gencode.VM38...parsed_ucsc_tableformat \
    --gff3 mmu.gff3 \                    # optional
    --assembly mm39 \
    --output UniqueGenomicElements.mm39.bed

python generate_master_filelist.py \
    --parsed-ucsc gencode.VM38...parsed_ucsc_tableformat \
    --repeatmasker mm39.repeatmasker.tsv.gz \
    --simplerepeats mm39.simplerepeats.tsv.gz \
    --trna mm39.trna.tsv.gz \            # optional
    --gff3 mmu.gff3 \                    # optional
    --custom-entries NR_046233.2.fasta \
    --output MASTER_FILELIST.{date}.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv
```

## UX Expectations

- Scripts must print progress to stderr (not suppress output)
- Scripts must exit non-zero on fatal errors with a descriptive message
- Optional inputs (trna, gff3) must be skipped gracefully when absent — no silent data corruption
- Missing IDs exceeding 1% of total must produce a visible warning and a report file, not a silent omission
- Runtime on TSCC for mm10/mm39 is expected to be minutes to hours (acceptable for one-time setup)
- No interactive prompts; all parameters via CLI flags

## Verification UX

The user verifies generated files by:
1. Running the hg38 reproduction test (scripts regenerate hg38 references → diff against originals)
2. Running `snakemake -n` (dry-run) with mm10/mm39 references to verify pipeline acceptance
3. Optionally running a full pipeline test with a small mm10/mm39 eCLIP dataset
