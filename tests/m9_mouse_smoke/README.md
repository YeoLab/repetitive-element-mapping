# M-9: the pipeline end-to-end on real mouse eCLIP data (`475.14`)

Phase 6 of `docs/ROADMAP-refdata.md`. Every mouse check before this one compared **static
files** — M-8's acceptance suite proves the three artifacts are internally consistent, not that
the pipeline can read them. This is the first run of the actual workflow on a mouse reference
set.

## The sample

`EV245` from `eric_ifit_vsv_clips` (Eric Van Nostrand, IFIT2/IFIT3 seCLIP, mm10+VSV) — the only
mouse eCLIP in the lab tree carrying **both** inputs this pipeline needs, the trimmed fastq and
the genome-mapped rmRep BAM. Config: `examples/repeat_mapping_SE_mm10.yaml`.

| | IP (`EV245_CLIP`) | INPUT (`EV245_INPUT`) |
|---|---|---|
| reads | 12,310,292 | 14,804,186 |
| bowtie2 to the mm10 repeat index | 36.64% | 56.44% |

Two properties of the sample, neither of which invalidates the run: the BAM was mapped to an
mm10+VSV hybrid, so its VSV contig is not covered by `UniqueGenomicElements.mm10.bed` and those
reads resolve to no unique element; and its STAR annotation was gencode vM15 against the
refdata's vM23 — both mm10, so coordinates agree.

## Result: PASS, 70/70 steps, exit 0

**The repeat arm is alive at mouse scale.** This is the check `475.16` exists for — full-dataset
hg38 runs once assigned `RepFamilyReads 0` while downsampled ones were fine, so a mouse run at
12 M reads is the only way to know the defect does not recur on a different assembly:

```
#READINFO  AllReads        6,451,354          9,997,487        (IP / INPUT)
#READINFO  UsableReads     3,950,961  0.612   7,378,033  0.738
#READINFO  GenomicReads    1,278,195  0.324   1,398,056  0.189
#READINFO  RepFamilyReads  2,672,766  0.676   5,979,977  0.811
```

**Multifamily fraction is comparable to human**, which is the M-9 criterion that needed an hg38
comparator:

| run | withpipes elements | nopipes | IP reads | multifamily | fraction |
|---|---|---|---|---|---|
| hg38 SE reference (`test-provenance/`) | 1,915 | 182 | 17,691,900 | 186,438 | **1.05%** |
| mm10 EV245 | 1,748 | 313 | 3,950,961 | 35,938 | **0.91%** |

Both carry 13 `unique_*` rows. Mouse resolves more distinct families (313 vs 182) because mm10's
rmsk names more of them, not because anything is being split.

**No family-resolution warnings** — every job log is clean. Two element names that look wrong for
mouse and are not: `antisense_Alu` is legitimate, because UCSC's mm10 `rmsk` files B1 SINEs
(`B1F`, `B1F1`, `B1F2` …) under `repFamily = Alu`; `antisense_L1`, `antisense_B4` and `MTB_MM_LTR`
are mouse-specific as expected.

**`-9j2` confirmed on real data.** The 10 U5 transcripts recovered by the per-assembly ambiguity
fix are not inert annotation: `RNU5G` carries **702 IP / 4,766 INPUT reads** here. Before the fix
mouse U5 rested on a single transcript.

Top families (IP reads):

```
RNA28S             1,615,534      antisense_L1        90,859
RNA18S               475,436      unique_3utr         66,079
unique_distintron    467,287      Simple_repeat       46,402
RNA45S               249,495      antisense_Alu       40,051
unique_proxintron    118,354      tRNA                33,667
unique_CDS           105,682      antisense_B4        26,621
RNA5S                 94,123      RNU5G                  702
```

## Files

| file | what it is |
|---|---|
| `SE.mm10.nopipes.tsv` | unambiguous families, fold enrichment and information content |
| `SE.mm10.withpipes.tsv` | all families including the `|`-joined multifamily rows |
| `SE.mm10.{IP,INPUT}.readinfo` | the `#READINFO` header, which the TSVs do not carry |

There is **no reference to diff against** — mouse has no gold standard, which is the whole
premise of M-8. These are a provenance record and a regression baseline for the next mouse run,
not a correctness proof.

## Reproducing

```bash
module load singularitypro
conda activate snakemake9

# From inside an interactive job, this is required first -- see below.
unset SLURM_JOB_ID
snakemake --configfile examples/repeat_mapping_SE_mm10.yaml \
  --profile profiles/tscc2_snakemake9
```

**`unset SLURM_JOB_ID` before using the SLURM profile from a compute node.** With it set,
Snakemake's SLURM executor dies immediately on `resolve_ctls_from_dns_srv: res_nsearch error:
Unknown host` / `DNS SRV lookup failed`. That is not `sbatch` being unusable there — bare
`sbatch --wrap` and `sbatch --test-only` both submit fine with the variable set. It is specific to
the executor, which takes a different path when it believes it is running inside a job.

The run recorded here predates that diagnosis and used the local fallback, which is also fine:

```bash
snakemake --configfile examples/repeat_mapping_SE_mm10.yaml \
  --cores 8 --use-conda --conda-prefix conda-env --resources mem_mb=32000
```

~28 min wall. The memory cap is required — `dedup` asks 32,000 `mem_mb`, so without it Snakemake
runs eight of them at once on a 32 GB node.
