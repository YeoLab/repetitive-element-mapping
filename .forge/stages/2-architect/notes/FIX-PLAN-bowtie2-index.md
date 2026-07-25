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
  RepBase↔RepeatMasker name-mapping tables in the RM Edition `Libraries/`, (c) genomic
  `getfasta` of a representative instance for any family present in `repeatmasker.tsv.gz`
  but in neither library. Re-measure coverage after aliasing; only then is ≥99% claimable.

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

1. Header count within 99% of 7,606, and **no duplicate headers** (`grep '^>' | sort | uniq -d`
   must be empty).
2. Family-header set overlap ≥99% vs reference (`comm -12` on sorted headers).
3. **Per-header sequence match**: for each shared header, compare a normalized sequence
   hash (uppercase, IUPAC→N applied to both) — require ≥99% of shared headers to hash-match
   the reference sequence. Report the non-matching headers, not just a count.
4. `bowtie2-inspect --summary` exit 0.
5. Clean header format (no `::coords` suffix) across ALL headers, not a spot-check.

```bash
# sequence-level overlap sketch
norm() { awk '/^>/{if(h)print h"\t"s; h=$0; s=""; next}{gsub(/[RYMKSWBDHV]/,"N",$0); s=s toupper($0)} END{print h"\t"s}' "$1"; }
join -j1 <(norm new.fa|sort) <(norm ref.fa|sort) | awk '$2==$3{ok++} END{print ok" seq-identical"}'
```

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
