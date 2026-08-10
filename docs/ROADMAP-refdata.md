# Reference Data Roadmap — hg38 provenance → mm10/mm39 generation

**Status as of 2026-08-10.** Tracked in beads under `repetitive-element-mapping-475`.
Run `bd show repetitive-element-mapping-475` for live status; this document is the rationale,
the evidence, and the phase structure.

## Why this document exists

The goal is mouse reference data (mm10/mm39) for the eCLIP repetitive element pipeline. Mouse
has **no reference outputs to validate against**, so mouse correctness cannot be established
directly. It can only be inherited — from a chain of evidence built on human, where reference
outputs *do* exist.

That chain is six phases. Every phase exists to make the next one interpretable:

| # | Phase | Why the next phase needs it |
|---|---|---|
| 1 | Run the pipeline on stock hg38 refdata, reproduce the reference outputs | Establishes that the *pipeline* is correct, so later discrepancies can be blamed on refdata |
| 2 | Document the structure of the four hg38 artifacts | You cannot regenerate a file whose schema and invariants are unstated |
| 3 | Identify the upstream sources that can regenerate each artifact | Names what regeneration is even made of |
| 4 | Gap analysis: auto-derivable vs manually curated | Isolates the parts that *cannot* transfer to mouse |
| 5 | Run the pipeline on regenerated hg38 refdata, compare to phase 1 | **The last checkpoint against ground truth.** Any refdata defect surviving here is invisible forever after |
| 6 | Generate mouse refdata by the same method | The deliverable |

The failure mode this structure guards against: regenerate mouse refdata, get plausible-looking
output, and have no way to know it is wrong. Phase 5 is the load-bearing step — it is the last
moment a generator bug meets a known-correct answer.

## Current status

| Phase | Status | Evidence |
|---|---|---|
| 1 | **BLOCKED — not established** | Full-dataset runs assign zero repeat-family reads (`475.16`) |
| 2 | Partial | Knowledge exists as defect narratives in `.forge/stages/2-architect/notes/`, not as a schema (`475.18`) |
| 3 | Essentially done, unrecorded | RepBase 18.05 source confirmed at 100% coverage; needs a source manifest (`475.19`) |
| 4 | Partial | Curated constants known but not enumerated as a class (`475.20`) |
| 5 | Not started | Blocked on phase 1 + both generators (`475.22`) |
| 6 | 0 of 4 mouse artifacts valid | All mm10/mm39 outputs predate a known fix (`475.6`–`475.15`) |

### Phase 1 is broken, and this was not previously known

Measured 2026-08-10 across the four committed runs in `results/`:

| Run | AllReads | RepFamilyReads |
|---|---|---|
| `se_small` | 2,077 / 2,206 | 1,162 (0.559) / 1,784 (0.809) |
| `pe_small` | 1,941 | 1,444 (0.744) |
| `se_full` | 10,272,097 | **0 (0)** |
| `pe_full` | 3,954,802 | **0 (0)** |

Small datasets assign 56–81% of reads to repeat families. Full datasets assign **exactly zero**,
in both SE and PE. The repeat-mapping arm — the arm every reference file in this project feeds —
is dead at scale.

The reference data is not the cause. `examples/inputs/hg38/bowtie2_index/` is the authentic
downloaded index (7,606 headers, zero `::coords` suffixes) and `repeat_mapping_SE_full.yaml`
points at it plus the 2020 reference MASTER_FILELIST. The pipeline fails with known-good inputs.

Comparing `results/se_full/seCLIP_example.nopipes.tsv` to the SE reference:

- Total assigned IP reads: **8,403,956 vs 17,505,462**
- `RNA28S`: **1,457 vs 5,399,580** (~3,700× low)
- 30 element classes absent entirely: `RNA45S`, `RNA5-8S`, `MTTRNA`, `MTRNR1/2`, `RNU5A`–`RNU5F`,
  `snRNA`, `Other`, `HSFAU`, `Eutr`, and their antisense partners
- Zero elements present in ours but not the reference
- `unique_*` genomic rows agree to ~0.005% (`unique_distintron` 2,444,022 vs 2,443,909)

That last row is doubly informative: it confirms `EXAMPLE_SE` is the same underlying dataset as
the reference `INV_B`, **and** that the genomic arm is fine while only the repeat arm fails.

`.forge/debug/causal-chain.json` diagnosed an identical 0-byte-`Rep.sam` symptom on CWL in
May 2026 as a *transient* bowtie2 failure that "is NOT reproducible under identical conditions
today." **That conclusion no longer holds.** It reproduces deterministically on full datasets
under Snakemake. It is scale-dependent, not transient. The prime suspect from that chain still
stands: bowtie2 dies without stdout, and the pipe reader treats EOF as success with no
exit-status check — so the pipeline reports `RepFamilyReads 0` instead of failing.

### The reference outputs have no provenance

`test-provenance/` holds the SE and PE reference `.nopipes.tsv`/`.withpipes.tsv`, but
`configs/`, `manifests/`, `scripts/` and `reference-metadata/` are **all empty**. Nothing records
which inputs, which reference data, which tool versions, or which command produced them. That
`EXAMPLE_SE` corresponds to `INV_B` is inferred from matching read counts above — not documented.
Until that is pinned down, "reproduce the reference" is not a well-defined target (`475.17`).

## Phase detail

### Phase 1 — human baseline (`475.21`)

Reproduce the reference outputs using the original downloaded hg38 refdata.

- `475.16` **P0** — fix the zero-repeat-reads failure at scale. Capture bowtie2's stderr this time
  rather than inferring; add an exit-status check so silent death fails loudly.
- `475.17` — reconstruct reference-output provenance; confirm the `EXAMPLE_SE`↔`INV_B` and
  `EXAMPLE_PE`↔`204_01_RBFOX2` correspondence by checksum.

Decide the acceptance tolerance explicitly. Byte equality may be unachievable: Perl 5.18+
non-deterministic hash iteration affects dedup tie-breaking (partially mitigated by sorting hash
keys). State a per-column band and justify it, rather than assuming exactness and discovering
otherwise.

### Phase 2 — artifact structure (`475.18`)

Consolidate into one schema document. Two invariants already learned the hard way:

- **MASTER_FILELIST row order is behaviorally significant.** `read_in_filelists`
  (`parse_bowtie2_output_realtime_includemultifamily_SE.pl:463`) assigns `priority_n` in file
  order, and that priority breaks equal-score mapping ties. Never sort it.
- **`parsed_ucsc_tableformat` is order-independent** — its consumer keys by transcript id.

### Phase 3 — sources (`475.19`)

Largely settled; needs recording, not rediscovery.

| Artifact | Source |
|---|---|
| Repeat consensus | RepBase 18.05 `species_specific/homo_sapiens_repbase_fixed_v2.fasta` |
| Gene models | Gencode v33 |
| Repeat instances | UCSC RepeatMasker + simpleRepeat tracks |
| tRNA | gtRNAdb |
| miRNA | miRBase gff3 |
| RNA45S | RefSeq `NR_046235.1` |

The RepBase file has 1,356 records covering **1,224/1,224** index families (100%), verified by
header prefix *and* canonical sequence digest: 700 byte-exact, 524 matching after IUPAC→N (which
is precisely the `.fixed.fa` step).

**Rejected, do not revisit:** the RepBase 24.01 two-library method and the RepeatMasker-Edition
EMBL parse. The 97.96% coverage figure that motivated them was an artifact of missing
`pseudo.ref`/`vrtrep.ref`/`invrep.ref`, not a real library gap.

### Phase 4 — gap analysis (`475.20`)

The pivotal phase for mouse feasibility. Known curated residue that falls out of no source file:

1. The 12-entry `DROP_EXACT` list (`refdata/repbase18.05-drop-exact.txt`) — derived by diffing
   against the hg38 reference index
2. The `gene_name`→family map with its `RNU`/`YRNA`/`RN7SL`/`SNORD` special cases
3. The chromosome allowlist, pinned to an assembly release
4. `GENCODE_EXCLUDED_FAMILIES` = {`RNA5S`, `RNA5-8S`, `MTTRNA`, `MTRNR1`, `MTRNR2`}
5. MASTER_FILELIST row order

Each is a place where hg38 could be reproduced **by fitting to the answer** — and where mouse has
no answer to fit to. Classify every one as derivable-from-source, derivable-from-rule, or
irreducibly-curated; for the last category, state how mouse gets a value and how it is checked.

### Phase 5 — regenerated-hg38 validation (`475.22`)

Regenerate all four hg38 artifacts, rerun SE and PE, compare to the phase-1 baseline. Requires
T-05 (`-7ee`) and T-09 (`-gz2`) to land first. The residual divergence measured here becomes the
documented expectation ceiling for mouse — mouse cannot be held to a standard human did not meet.

### Phase 6 — mouse generation (`475.6`–`475.15`)

All four mm10/mm39 artifacts are currently invalid or unverified; the files on disk date to
2026-05-14 and each predates a known fix.

| Artifact | Status | Why |
|---|---|---|
| bowtie2 index | Invalid | T-05 bug on disk: 5,874/5,875 (mm10) and 22,356/22,357 (mm39) headers carry `::chr:start-end(strand)` — genomic instances, not family consensus |
| MASTER_FILELIST | Invalid | T-09 bug: col4 holds `misc_RNA`/`Mt_tRNA` (biotype) where hg38 holds `RNU1` (family). Also mm10=11,036 vs mm39=22,368 rows — inconsistent runs |
| UniqueGenomicElements | Stale | 275/277 MB vs hg38's 208 MB; predates the T-07 five-rule fix, so still carries ~1M over-generated simple-repeat rows |
| parsed_ucsc_tableformat | Unverified | Generator is validated (T-03 PASS), but outputs predate its commit by a week |

Sequence:

```
M-1  475.6   mm39 tRNA + GRCm39 miRNA inputs        [ready]
M-2  475.7   mm10/mm39 chrom allowlists             [ready]
M-7  475.8   re-run parsed_ucsc + checksum          [ready]
T-05 -7ee    index generator rewrite                [ready]
T-09 -gz2    filelist generator fix          <- T-05
M-3  475.9   mouse index generation          <- T-05, M-1, P-5
M-4  475.10  mouse MASTER_FILELIST           <- T-09, M-3, P-5
M-5  475.11  mouse DROP_EXACT re-derive      <- M-4
M-6  475.12  mouse UniqueGenomicElements     <- M-1, M-2, M-4
M-8  475.13  mouse acceptance suite          <- M-5, M-6, M-7
M-9  475.14  end-to-end pipeline smoke       <- M-8, P-1a
M-10 475.15  package + document              <- M-9
```

Three gaps found 2026-08-10 that were in no document or issue:

1. **mm39 is missing two inputs entirely** — no `mm39.trna.tsv.gz`, and no GRCm39-coordinate
   miRBase gff3 (mm10's `mmu.gff3` is GRCm38). T-07 rules 4 and 5 need both.
2. **No mouse chromosome allowlists** — `refdata/hg38.chrom-allowlist.txt` exists, the generator
   takes `--chrom-allowlist` as required, mm10/mm39 have none.
3. **No mouse acceptance criteria existed.** Every hg38 criterion is "diff against the reference."
   Mouse has none. `M-8` replaces reference-diff with cross-artifact consistency checks — and its
   acceptance criterion requires those checks to *fail* on the known-bad hg38 outputs before they
   are trusted on mouse.

## Standing risks

- **Line-count acceptance criteria mask content bugs.** T-09 passed a ±1% row-count check
  (26,252 vs 26,422 = 0.9936) while col4 was entirely the wrong field. Any new criterion must
  test content, not shape.
- **`M-5` is a correctness dependency, not bookkeeping.** The mouse `DROP_EXACT` prediction of
  1,079 families assumes the human list transfers. The actual rule is "drop whatever Gencode
  supplies," which cannot be checked until a mouse MASTER_FILELIST exists. Getting it wrong
  double-counts reads across the repeat and Gencode portions.
- **Genomic-instance fallback for unresolved families is prohibited** (FIX-PLAN §4a). It
  reintroduces the T-05 error, and mouse has no reference to catch it.
- **Fitting to hg38 is the systemic hazard.** Rules derived by diffing against the hg38 reference
  carry no independent evidence for mouse. Phase 4 exists to make that explicit, artifact by
  artifact.

## Pointers

- `.forge/stages/2-architect/notes/FIX-PLAN-bowtie2-index.md` — T-05 rewrite plan (§2–§4b, §7)
- `.forge/stages/2-architect/notes/T-07-per-source-selection.log` — the five UniqueGenomicElements rules
- `.forge/debug/causal-chain.json` — prior (now superseded) 0-repeat-reads diagnosis
- `docs/CHANGELOG-refdata-validation-2026-07-25.md` — hg38 reproduction validation results
- `bd memories` — `repbase-build-method-confirmed`, `repbase-consensus-source`,
  `refdata-validation-policy`, `refdata-review-corrections`, `refdata-hg38-reproduction-results`
