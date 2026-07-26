# Fix Plan — `generate_bowtie2_index.py` (T-05), with knock-on to T-09

Status: DRAFT for review · 2026-07-25 · issue repetitive-element-mapping-475

## 1. Confirmed diagnosis

The hg38 bowtie2 index FASTA (7,606 seqs) = 5,002 Gencode transcripts + 864 tRNA +
15 rRNA (NR_) + **1,725 repeat entries** (1,224 RepBase family consensus + 501
SimpleRepeat kmers). The repeat consensus sequences come from **RepBase 24.01**, not the
genome. Current script does `pybedtools getfasta` on every genomic repeat *instance* →
26,353 seqs with coordinate-suffixed headers. Wrong unit (instance vs family), wrong
header format, wrong source.

## 2. Source (located, owned by bay001)

`/tscc/projects/ps-yeolab4/genomes/RepBase24.01.fasta/RepBase24.01.fasta/*.ref`
Header format `>NAME<TAB>CLASS<TAB>species`. Also `RepBaseRepeatMaskerEdition-20181026.tar.gz`
(more families) in the parent dir.

## 3. Proposed rewrite

Keep the working parts, replace the repeat part:

| Portion | Reference count | Current | Plan |
|---|---|---|---|
| Gencode transcripts | 5,002 | getfasta ✅ but header has `::coords` | getfasta, **strip coord suffix** → clean `>ENST…` |
| tRNA | 864 | ? | from tRNA source, clean names |
| rRNA (NR_) | 15 | NR_ fasta | keep |
| Repeat family consensus | 1,224 | getfasta per instance ❌ | **read from RepBase `.ref`; one seq/family; UPPERCASE name header** |
| SimpleRepeat kmers | 501 | ? | keep existing UpdatedSimpleRepeat logic (not RepBase simple.ref — 0 overlap) |

New CLI: `--repbase-refs <dir>` + `--species {human,mouse}` selecting the ref-file set:
- human: humrep, humsub, prirep, prisub, mamrep, mamsub
- mouse: rodrep, rodsub, mousub, ratsub, mamrep, mamsub

## 4. Source method — CONFIRMED by verbatim sequence matching

Compared reference-index sequences against candidate sources (byte-level):

| Family | in RepBase24.01 .ref | in RM Edition | index matches |
|---|---|---|---|
| LTR18B, MSTC | yes | yes | both verbatim |
| MER5A | yes (2 IUPAC diffs→N) | 5 real diffs | **.ref** (primary) |
| EULOR1 | no | yes (verbatim) | **RM Edition** (fallback) |

Confirmed build method (mirrors the original hg38 pipeline):
1. **Primary source: RepBase 24.01 `.ref` consensus** (taxon-split, species-selected).
2. **Fallback: RepeatMasker Edition `Libraries/RMRBSeqs.embl`** (EMBL, 49,011 seqs) for
   families absent from RepBase 24.01 (EULOR/DEUSINE-class).
3. **IUPAC→N normalization** (the `.fixed.fa` step): R/Y/M/… → N.

Coverage of the 1,224 hg38 index families:
- RepBase 24.01 alone: 1116; RM Edition alone: 864; **union: 1199/1224 = 97.96%**.
- **This does NOT yet meet the ≥99% gate** (needs ≥1212/1224). At least **13 of the 25
  residual** families (CR1L, DEUSINE, E1/2/3, UCON*, L7/L23/L28…) must be resolved first.
- **Required before implementation: an alias-resolution pass.** The residuals are candidate
  name variants; resolve via (a) suffix/case normalization (`DEUSINE`↔`DEUSINE1`), (b) the
  RepBase↔RepeatMasker name-mapping tables in the RM Edition `Libraries/`. Re-measure
  coverage after aliasing; only then is ≥99% claimable.

### 4a. Consensus-source policy (no genomic fallback)

A repeat-family sequence in this index is a **family consensus**. A genomic instance is not
a consensus — substituting one reintroduces exactly the source error that caused T-05, and
mouse has no gold reference that would catch it. Therefore:

- **Permitted sources, in priority order, and nothing else:**
  1. RepBase 24.01 `.ref` consensus (species-selected taxon files).
  2. RepeatMasker Edition `Libraries/RMRBSeqs.embl` consensus.
  3. An explicit, recorded alias into (1) or (2).
- **`getfasta` of a representative RepeatMasker instance is prohibited** as a source for any
  repeat-family sequence. (It remains the correct method for the *Gencode transcript*
  portion of the FASTA, which is genomic by definition.)
- **Every emitted repeat family must carry provenance.** The generator writes a sidecar TSV
  `<output>.provenance.tsv` with `family, source ∈ {repbase24.01, rm_edition}, source_file,
  resolved_identifier, alias_of_or_'-'`. A family with no provenance row is a build error.
- **Unresolved families cannot pass silently.** A family selected by §4b but absent from
  both libraries after aliasing is written to `<output>.unresolved.txt` and **omitted** from
  the FASTA. The generator exits non-zero unless the caller passes
  `--allow-unresolved <n>` acknowledging a count, and the coverage gate is computed over the
  *selected* set — so omissions lower coverage rather than being hidden by a smaller
  denominator.
- **Mouse uses the identical policy.** No hg38-only leniency: mm10/mm39 emit the same
  provenance sidecar and the same non-zero exit on unacknowledged unresolved families.

## 4b. Family-SELECTION rule (which families to include)

Sequence source is settled; the remaining knob is *which* families make the list. The hg38
index = species-primary refs ~wholesale + ancestral/RM-Edition families present in the
genome. Implement as:

    select = (species_primary_refs)                         # human: humrep,humsub
                                                             # mouse: rodrep,rodsub,mousub,ratsub
           ∪ {fam in genome repeatmasker : consensus exists in (RepBase24.01 ∪ RM Edition)}
    sequence(fam) = RepBase24.01[fam]  if present  else RM_Edition[fam]   # then IUPAC→N

Validate on hg38: iterate this rule to ≥99% family-header overlap with the reference index,
then apply with mouse ref-set to mm10/mm39.

## 5. Validation (hg38, before touching mouse)

Header-level checks are necessary but **not sufficient** — the MER5A example proves two
FASTAs can share a header yet differ in sequence. Validate at the **sequence** level:

Checks 1–3 are implemented by `bin/python/refdata_generation/validate_fasta_equivalence.py`
(tests: `tests/refdata_generation/test_validate_fasta_equivalence.py`). Do **not** use an
ad-hoc AWK one-liner: the obvious sketch substitutes `[RYMKSWBDHV]` *before* uppercasing, so
lowercase ambiguity codes (`rymk`) survive as `RYMK` instead of `NNNN` and soft-masked
records mis-compare. The validator uppercases first, then maps every non-ACGT byte to `N`.

1. **No duplicate headers** — the validator raises and exits 2 on any repeated header.
2. **Header set overlap** — reports `shared`, `missing` (in reference, not in new) and
   `extra`; `--max-missing` gates the missing count (default 0).
3. **Per-header sequence match** — SHA-256 of each canonicalized record is compared, and
   `--min-match` gates `matching/shared` (default 0.99). Non-matching headers are listed,
   not just counted. Records are streamed, so memory is O(headers), not O(sequence).
4. `bowtie2-inspect --summary` exit 0.
5. Clean header format (no `::coords` suffix) across ALL headers, not a spot-check.

```bash
python bin/python/refdata_generation/validate_fasta_equivalence.py \
  new.fa examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa \
  --min-match 0.99 --max-missing 0
# new=7606 reference=7606 shared=7606
# missing=0 extra=0
# matching=7606 mismatching=0 match_fraction=1.00000 (threshold 0.99)
# PASS
```

Confirmed on real data: the reference index has 7,606 unique headers (no duplicates), and
the MER5A `.ref`-vs-index IUPAC-only difference is correctly absorbed as a match, while a
genuine single-base difference is reported as `mismatching`.

Then regenerate mm10/mm39 with `--species mouse` and re-run checks 1–5.

## 6. Coupling / scope

- **T-09 (MASTER_FILELIST)** — reusing the T-05 repeat-family set fixes only the *repeat*
  rows. It does NOT fix the **Gencode** rows, which are the larger problem. Two derivations
  are needed, plus order preservation:
  1. **Gencode cols 3–5** (currently `generate_master_filelist.py:195` wrongly uses the
     biotype for col4). Reference convention, e.g. `RNU1-1 → RNU1`:
     `col3 = gene_name`; `col4 = curated family` (strip the copy-number suffix from
     gene_name: `RNU6-2→RNU6`, group `Y_RNA→YRNA`); `col5 = genelists.{family}`.
     A gene_name→family map (with the RNU/YRNA/RN7SL/SNORD/… special cases) must be defined.
  2. **Repeat cols 4–5** need the RepBase family→class mapping (`ALUY → Alu → SINE`), not the
     current `(name,name,name,name,name)`; col5 = RepBase class (SINE/LINE/DNA).
  3. **Row order** must match the reference (or be validated behaviorally) — it sets mapping
     priority via `priority_n` (see §1). Do NOT sort.
  Sequence: do T-05 first for the repeat family set, then T-09 adds the two mapping tables.
- **T-07 (UniqueGenomicElements)** — schema is correct; the fix is **selection only**:
  restrict the Gencode contribution from all 249,043 transcripts to the curated 4,670
  subset. Separate, coordinate-based fix; do not touch the column layout.

## 7. Effort / risk

- Rewrite repeat portion + RepBase `.ref` parser + EMBL parser + species ref-set selection:
  moderate.
- Risk: the residual 25 families (§4) — must run the alias-resolution pass to reach the ≥99%
  gate before implementation is accepted.
- Env to run: pybedtools+bedtools (snakemake738), samtools (ecliprepmap-0.1.0),
  bowtie2-build (bowtie2-2.5.4).
