# P-5: hg38 run on fully regenerated reference data (`475.22`)

Phase 5 of `docs/ROADMAP-refdata.md` — **the last point at which a refdata defect meets
ground truth.** Mouse has no reference, so any error surviving this comparison is invisible
from phase 6 onward.

All four hg38 artifacts were regenerated from source by the committed generators, then SE and
PE were run on them and compared against the P-1 baseline (`results/{se,pe}_full`, which
themselves reproduce the CWL 1.0.0 references exactly — SE 182/182 read counts, PE 169/169).

| file | what it is |
|---|---|
| `SE.baseline.nopipes.tsv`, `PE.baseline.nopipes.tsv` | P-1 baseline, stock 2020 reference data |
| `SE.regenerated.nopipes.tsv`, `PE.regenerated.nopipes.tsv` | same pipeline, fully regenerated reference data |
| `SE.{baseline,regenerated}.readinfo` | the `#READINFO` header, which the TSVs do not carry |

Compare with:

```bash
python3 bin/python/compare_pipeline_outputs.py \
  tests/p5_regenerated_hg38/SE.regenerated.nopipes.tsv \
  tests/p5_regenerated_hg38/SE.baseline.nopipes.tsv
```

It reports FAIL, and that is the expected, documented result — `compare_pipeline_outputs.py`
requires read counts to match **exactly**, which is the right rule for pipeline equivalence and
the wrong rule for a refdata-version comparison. The measured divergence is in
`docs/CHANGELOG-refdata-validation-2026-07-25.md` §10.

The regenerated MASTER_FILELIST is not committed (1.7 MB, and it needs RepBase 18.05, which is
TSCC-local). It is reproducible with the command in changelog §10 and has md5
`2b6d6807fe4ebd7d822562d9cab1155c`.
