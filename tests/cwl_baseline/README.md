# CWL baseline outputs

Reference outputs from the **original CWL pipeline**, generated fresh so the
Snakemake translation can be checked against a run of known provenance. Unlike
`test-provenance/`, everything about these runs is recorded.

Compare with:

```bash
python3 bin/python/compare_pipeline_outputs.py \
    <snakemake>.nopipes.tsv \
    tests/cwl_baseline/ecliprepmap-1.0.0-SE/INV_B.IP.umi.r1.fqTrTr.sorted.fq.barcode1.nopipes.tsv
```

Read counts must match **exactly**; derived float columns are compared with a
relative tolerance. Row order is ignored — elements with equal read counts come
out in Perl hash-iteration order, which is not reproducible across
implementations or, post-5.18, across runs.

## `ecliprepmap-1.0.0-SE/` — INV_B, hg38

Generated 2026-08-10.

| | |
|---|---|
| Module | `ecliprepmap/1.0.0` (bowtie2 2.2.6, perl 5.10.1, toil 5.12.0) |
| Workflow | `eCLIP_repelement_SE` → `wf_ecliprepmap_se`, SLURM/toil |
| Job file | `/tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/repeat_mapping_SE.yaml` |
| Dataset | `INV_B` |
| Reference | `MASTER_FILELIST.20201203...`, Gencode v33, `UniqueGenomicElements.hg38.bed` |
| Result | 22,576,144 AllReads / 17,691,900 usable / 9,477,858 RepFamilyReads |

This is the same dataset behind `test-provenance/tests/ecliprepmap-1.0.0/`, so
it also documents the provenance that directory is missing.

### Reproducing

The shipped job file does **not** run as-is on TSCC2 (see issue `-475.28`).
Three fixes are needed:

1. Every `path:` uses the pre-migration `/projects/...` prefix — prepend `/tscc`.
2. The four `INV_B.*` files under `examples/inputs/` are broken symlinks into
   `/projects/ps-yeolab4/software/eclip/0.7.0/...`. The real files are at the
   `/tscc/...` equivalent.
3. **Run from a shared filesystem.** `eCLIP_repelement_SE` distributes steps as
   separate SLURM jobs, but `/scratch` is node-local NVMe — the toil file store
   lands on one node and the steps run on others. Use
   `/tscc/lustre/ddn/scratch/$USER`. The failure otherwise surfaces as a
   misleading `no <file>` from a Perl script, not a filesystem error.

```bash
module load ecliprepmap/1.0.0        # do NOT pipe this; a pipeline runs it in
                                      # a subshell and the PATH changes are lost
mkdir -p /tscc/lustre/ddn/scratch/$USER/cwl_baseline && cd $_
sed -e "s#/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/inputs/INV_B#/tscc/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/INV_B#g" \
    -e "s#: /projects/#: /tscc/projects/#g" \
    /tscc/projects/ps-yeolab4/software/ecliprepmap/1.0.0/examples/repeat_mapping_SE.yaml \
    > repeat_mapping_SE.yaml
chmod +x repeat_mapping_SE.yaml && ./repeat_mapping_SE.yaml
```

### Files

| File | What |
|---|---|
| `*.nopipes.tsv` | 182 unambiguous elements |
| `*.withpipes.tsv` | 1,915 elements including multi-family |
| `*.parsed.readinfo` | `#READINFO` header only — the full `.parsed` is 32 MB |

The `.parsed.readinfo` header is kept because it caught a defect the TSV
comparison could not: `UsableReads` was computed as `usable/usable` (always
`1.0`) instead of `usable/AllReads` (`-475.29`). Fold enrichment derives from
the `TOTAL` rows, so that field was wrong in every run while the TSVs compared
clean.
