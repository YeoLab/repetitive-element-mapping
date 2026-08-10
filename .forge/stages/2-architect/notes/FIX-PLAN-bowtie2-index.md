# Fix Plan — `generate_bowtie2_index.py` (T-05), with knock-on to T-09

Status: DRAFT for review · 2026-07-25 · issue repetitive-element-mapping-475

## 1. Confirmed diagnosis

The hg38 bowtie2 index FASTA (7,606 seqs) = 5,002 Gencode transcripts + 864 tRNA +
15 rRNA (NR_) + **1,725 repeat entries** (1,224 RepBase family consensus + 501
SimpleRepeat kmers). The repeat consensus sequences come from a **RepBase consensus
library** (identified in §2 as the RepBase 18.05 species-specific FASTA), not the
genome. Current script does `pybedtools getfasta` on every genomic repeat *instance* →
26,353 seqs with coordinate-suffixed headers. Wrong unit (instance vs family), wrong
header format, wrong source.

## 2. Source — RepBase 18.05 species-specific FASTA (SUPERSEDES the 24.01 plan)

```
/tscc/projects/ps-yeolab4/genomes/RepBase18.05.fasta/species_specific/
├── homo_sapiens_repbase_fixed_v2.fasta        1,356 records   (human)
├── mus_musculus_repbase_u1_fixed_v2.fastq     1,194 records   (mouse; FASTA content
│                                                               despite the .fastq name)
├── rat_rattus_repbase_fixed_v2.fasta                          (rat, unused here)
└── README                                     "*V2 indices contain a longer
                                                 NR_046235.1 (RNA45S). Using this one."
```

Header format: `>NAME_CLASS_TAXON` (e.g. `ALUY_SINE1/7SL_Primates`, `MER5A_hAT_Eutheria`).
Sequences are lowercase; the `_fixed_v2` naming matches the `.fixed.fa` in the reference
index filename, and the Dec-2017 mtime fits a 2020-12-03 index build.

**This one file covers 1,224/1,224 hg38 index repeat families (100%)**, verified two
independent ways: by `<FAMILY>_` header prefix, and by SHA-256 of the canonicalized
sequence. 700 are byte-exact; the remaining 524 match after IUPAC→N — i.e. exactly the
`.fixed.fa` normalization step, with nothing left over.

> **Superseded approach.** An earlier draft reconstructed the repeat set from
> `RepBase24.01.fasta/*.ref` (9 taxon files, species ref-set selection) with a
> `RepBaseRepeatMaskerEdition-20181026` EMBL fallback, reaching 1,199/1,224 = 97.96% and
> requiring an alias-resolution pass to clear a ≥99% gate. That is no longer needed: the
> 25 "missing" families were not absent from RepBase, only from the *human ref-set file
> selection* (18 live in `pseudo.ref`, 6 in `vrtrep.ref`, 1 in `invrep.ref`). The 18.05
> species file makes the second library, the EMBL parser, the alias table, and the
> coverage gate all unnecessary — this is an exact reproduction, not a thresholded one.

## 3. Proposed rewrite

Keep the working parts, replace the repeat part:

| Portion | Reference count | Current | Plan |
|---|---|---|---|
| Gencode transcripts | 5,002 | getfasta ✅ but header has `::coords` | getfasta, **strip coord suffix** → clean `>ENST…` |
| tRNA | 864 | ? | from tRNA source, clean names |
| rRNA (NR_) | 15 | NR_ fasta | keep |
| Repeat family consensus | 1,224 | getfasta per instance ❌ | **read from the 18.05 species FASTA; one seq/family; UPPERCASE name header** |
| SimpleRepeat kmers | 501 | ? | keep existing UpdatedSimpleRepeat logic |

New CLI: `--repbase-species-fasta <path>` (one file). No `--repbase-refs` directory, no
species ref-set list, no EMBL fallback path.

Emitting a repeat family is then: parse `NAME_CLASS_TAXON` → apply §4b filter →
uppercase → IUPAC→N → write `>NAME`.

## 4. Source method — CONFIRMED by 100% canonical sequence match

| Check (human, 1,224 index repeat families) | Result |
|---|---|
| resolvable by `<FAMILY>_` header prefix | 1,224 / 1,224 |
| sequence match after IUPAC→N canonicalization | **1,224 / 1,224 (100%)** |
| byte-exact (uppercase only, no IUPAC→N) | 700 / 1,224 |

Build method:
1. **Single source: RepBase 18.05 species-specific FASTA** (`_fixed_v2`, per README).
2. **IUPAC→N normalization** (the `.fixed.fa` step) — accounts for all 524 non-byte-exact
   families and nothing else.

There is no fallback source and no residual: coverage is exact, so the ≥99% gate that the
previous draft imposed no longer applies. Any future shortfall is therefore a real defect,
not an expected tolerance.

### 4a. Consensus-source policy (no genomic fallback)

A repeat-family sequence in this index is a **family consensus**. A genomic instance is not
a consensus — substituting one reintroduces exactly the source error that caused T-05, and
mouse has no gold reference that would catch it. Therefore:

- **Permitted source, and nothing else:** the RepBase 18.05 species-specific FASTA (§2).
  One file per species; no second library, no alias table.
- **`getfasta` of a representative RepeatMasker instance is prohibited** as a source for any
  repeat-family sequence. (It remains the correct method for the *Gencode transcript*
  portion of the FASTA, which is genomic by definition.)
- **Every emitted repeat family must carry provenance.** The generator writes a sidecar TSV
  `<output>.provenance.tsv` with `family, source_file, source_header, dropped_reason_or_'-'`.
  A family with no provenance row is a build error. `source_header` retains the full
  `NAME_CLASS_TAXON` string, which is also what T-09 reads for the family→class map.
- **Unresolved families cannot pass silently.** A family selected by §4b but absent from the
  source is written to `<output>.unresolved.txt` and **omitted** from the FASTA. The
  generator exits non-zero unless the caller passes `--allow-unresolved <n>` acknowledging a
  count. For hg38 the expected count is **0** — coverage is exact (§4), so any unresolved
  family is a real defect rather than an accepted tolerance.
- **Mouse uses the identical policy.** No hg38-only leniency: mm10/mm39 emit the same
  provenance sidecar and the same non-zero exit on unacknowledged unresolved families.

## 4b. Family-SELECTION rule (which families to include) — DERIVED, exact on hg38

The source file is a *superset* of the index (1,356 vs 1,224 for human). The 132 excluded
records are exactly those supplied by another portion of the FASTA, so the rule is
"drop what a dedicated source already provides":

    keep(record) unless:
      class ∈ {Simple_Repeat, tRNA, rRNA}      # -> SimpleRepeat kmers / tRNA file / RefSeq NR_
        EXCEPT headers starting 'MamSINE1_tRNA_'   # tRNA-derived SINE: a real repeat family
      OR header ∈ DROP_EXACT                   # -> supplied by Gencode/MASTER_FILELIST

    DROP_EXACT (12):
      U1_snRNA_Homo_sapiens     U2_snRNA_Homo_sapiens    U4B_snRNA_Homo_sapiens
      U5B1_snRNA_Homo_sapiens   U6_snRNA_Homo_sapiens    U7_snRNA_Vertebrata
      7SK_Pseudogene_Homo_sapiens   7SL_Pseudogene_Homo_sapiens
      HY1_Pseudogene_Homo_sapiens   HY3_Pseudogene_Homo_sapiens
      GGAAT_SAT_Homo_sapiens        NR_046235.1

**Verified: human 1,356 → 1,224, exact match to the reference index, zero symmetric
difference.**

`DROP_EXACT` is not arbitrary — those ten families are precisely what MASTER_FILELIST
supplies from Gencode (`RNU1, RNU2, RNU4, RNU5, RNU6, RNU7, RN7SK, RN7SL, YRNA`). It is the
same Gencode-vs-repeat partition found in T-07: whatever Gencode provides is removed from
the repeat portion so reads are not double-counted. Two traps this encodes:

- A blanket `_snRNA_` class drop is **wrong** — the index *keeps* `U3, U8, U14, UHG` and
  `E1/E2/E3`. Only the six Gencode-supplied U-snRNAs are dropped.
- Family names are not `header.split('_')[0]`: `HY1_SINE_Pseudogene_*` (kept) would be
  swallowed by a prefix match on `HY1` (dropped). Match `DROP_EXACT` on the **full header**.

Because it is a curated 12-entry constant rather than a general rule, commit it as data
(`refdata/repbase18.05-drop-exact.txt`) with this derivation recorded alongside.

### Mouse

The same rule transfers with **no mouse-specific additions**: `mus_musculus_repbase_u1_fixed_v2.fastq`
1,194 → **1,079** (dropped 115 = 68 Simple_Repeat + 39 tRNA + 6 rRNA + 2 exact). The mouse
file contains no `U1/U2/U4B/U5B1/U6/7SK/7SL/HY1/HY3` records at all, so only
`U7_snRNA_Vertebrata` and `NR_046235.1` apply. Its lone Pseudogene record,
`L32_Pseudogene_Mus_musculus`, is correctly **kept** — the analogue of human `L7/L23A/L28`,
which are index families.

**Caveat:** 1,079 is a *prediction*. There is no mouse reference index to validate against,
which is exactly why §4a's provenance sidecar and non-zero-exit rules are load-bearing for
mouse. Re-verify the mouse DROP_EXACT set once a mouse MASTER_FILELIST exists (T-09): any
family the mouse MASTER_FILELIST supplies from Gencode must be added to it.

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
   `--min-match` gates `matching/shared`. Non-matching headers are listed, not just counted.
   Records are streamed, so memory is O(headers), not O(sequence).
   **Use `--min-match 1.0 --max-missing 0` for hg38**: §4 establishes the repeat portion
   reproduces exactly, so the old 0.99 tolerance would now mask a regression.
4. `bowtie2-inspect --summary` exit 0.
5. Clean header format (no `::coords` suffix) across ALL headers, not a spot-check.

```bash
python bin/python/refdata_generation/validate_fasta_equivalence.py \
  new.fa examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa \
  --min-match 1.0 --max-missing 0
# new=7606 reference=7606 shared=7606
# missing=0 extra=0
# matching=7606 mismatching=0 match_fraction=1.00000 (threshold 1.0)
# PASS
```

Confirmed on real data: the reference index has 7,606 unique headers (no duplicates), and
the MER5A source-vs-index IUPAC-only difference is correctly absorbed as a match, while a
genuine single-base difference is reported as `mismatching`.

Then regenerate mm10/mm39 with `--repbase-species-fasta mus_musculus_repbase_u1_fixed_v2.fastq`
and re-run checks 1, 4, 5 (2 and 3 need a reference, which mouse lacks — see §4b caveat).

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
     **No separate mapping table is needed** — the class is the middle field of the 18.05
     header (`ALUY_SINE1/7SL_Primates`), so T-05's provenance sidecar (§4a,
     `source_header`) already carries it. Parse it there rather than building a second map.
  3. **Row order** must match the reference (or be validated behaviorally) — it sets mapping
     priority via `priority_n` (see §1). Do NOT sort.
  Sequence: do T-05 first for the repeat family set, then T-09 adds the Gencode derivation.
- **T-07 (UniqueGenomicElements)** — ~~pending~~ **DONE 2026-07-26.** Schema was correct;
  the fix was selection across five sources. See
  `.forge/stages/2-architect/notes/T-07-per-source-selection.log`. Note the Gencode rule
  there (MASTER_FILELIST minus rRNA/mito families) is the *same partition* as §4b's
  `DROP_EXACT` here — the two files must stay consistent.

## 7. Effort / risk

- Rewrite repeat portion + one `NAME_CLASS_TAXON` FASTA parser: **small**. The 18.05 source
  (§2) removes the `.ref` taxon-set selection, the EMBL parser, the fallback rule, and the
  alias table that the previous draft required.
- Risk is now concentrated in **mouse**, not human: hg38 is an exact reproduction
  (1,224/1,224), but the mouse 1,079 is unvalidated — no mouse reference index exists. §4a's
  provenance sidecar and non-zero exit are the only guards there.
- Residual risk: `DROP_EXACT` is a curated 12-entry constant derived from hg38; its mouse
  equivalent cannot be confirmed until a mouse MASTER_FILELIST exists (T-09).
- Env to run: pybedtools+bedtools (snakemake738), samtools (ecliprepmap-0.1.0),
  bowtie2-build (bowtie2-2.5.4).
