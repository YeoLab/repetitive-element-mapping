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
- RepBase 24.01 alone: 1116; RM Edition alone: 864; **union: 1199/1224 (98%)**.
- 25 residual (CR1L, DEUSINE, E1/2/3, UCON*…) — likely name variants; within the
  content-equivalence bar (≥99% overlap target per 2026-07-25 policy).

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

1. Header count within 99% of 7,606.
2. Family-header set overlap ≥99% vs reference (`comm -12` on sorted headers).
3. `bowtie2-inspect --summary` exit 0.
4. Spot-check clean header format (no `::coords`).

Then regenerate mm10/mm39 with `--species mouse` and sanity-check header counts +
`bowtie2-build`.

## 6. Coupling / scope

- **T-09 (MASTER_FILELIST)** shares this repeat-family list and currently mislabels col4
  with Gencode biotypes. Fixing the family set here feeds the T-09 fix (repeat rows should
  carry repeat-family names + RepBase class in col4/col5). Do T-05 first, then T-09 reuses
  the family set.
- **T-07 (UniqueGenomicElements)** is genomic-coordinate based (independent); separate fix.

## 7. Effort / risk

- Rewrite repeat portion + RepBase `.ref` parser + species ref-set selection: moderate.
- Risk: the last ~9% of families (RM Edition) — mitigated by choosing option A/B/C.
- Env to run: pybedtools+bedtools (snakemake738), samtools (ecliprepmap-0.1.0),
  bowtie2-build (bowtie2-2.5.4).
