# P-4: what is derivable, what was curated, and what was fit to the answer

> Phase 4 of `docs/ROADMAP-refdata.md` (`475.20`). The central question for whether mouse can be
> generated at all. Companion to `docs/MASTER_FILELIST-decisions.md`, which records *why* each rule
> is what it is; this file records *where each rule's authority comes from* and whether that
> authority extends to mouse.

Three classifications are used throughout:

| | meaning |
|---|---|
| **derivable-from-source** | falls out of a public input file. Regenerating it for a new assembly is a download. |
| **derivable-from-rule** | a stated rule applied to a public input. Transfers if the rule holds; the rule is the thing to check. |
| **irreducibly-curated** | a human choice with no source to regenerate it from. Mouse needs its own, and its own validation. |

A fourth label, **fit-to-hg38**, is orthogonal and is the real hazard: a constant obtained by
diffing against the hg38 reference reproduces hg38 *by construction* and therefore carries **no
independent evidence** that it is right — for hg38 or for mouse. Every such item is flagged.

---

## Summary

| # | Item | Class | Fit to hg38? | Mouse status |
|---|---|---|---|---|
| 1 | `repbase18.05-drop-exact.txt` (12) | irreducibly-curated | **yes** | own list derived by rule (`475.11`), 9 entries |
| 2 | `GENE_NAME_FAMILY_RULES` (14) | irreducibly-curated | partly | transfers after `-4ee`; validated no-op for human |
| 3 | chromosome allowlists | derivable-from-rule | hg38 yes | mouse pinned to UCSC `chrom.sizes` (`475.7`) |
| 4 | `GENCODE_EXCLUDED_FAMILIES` (5) | irreducibly-curated | **yes** | 5/5 on mm39, 4/5 on mm10 |
| 5 | MASTER_FILELIST row order | derivable-from-rule | n/a | behaviourally pinned by test, not by the reference |
| 6 | `repbase18.05-classes.txt` (76) | derivable-from-rule | **yes** | validated: 0 unparsed on mouse |
| 7 | `RMSK_REPNAME_TO_FAMILY` (21) | irreducibly-curated | **yes** | 16/21 fire on mouse; no mouse gap found |
| 8 | **`RMSK_AMBIGUOUS_REPNAMES` (2)** | irreducibly-curated | **yes** | **costs mouse 10 transcripts — see §7** |
| 9 | `hg38.repbase-class-family.tsv` (1,224) | irreducibly-curated | **yes** | mouse has its own 37-row table (`475.43`) |
| 10 | `*.rfam-family.tsv` | derivable-from-source | no | regenerated for mouse from RNAcentral |
| 11 | `*.rmsk-{class-family,smallrna}` | derivable-from-source | no | regenerated per assembly |
| 12 | `DROP_CLASSES` + `KEEP_PREFIX` | derivable-from-rule | **yes** | validated: the one exception exists in mouse too |
| 13 | index structural constants | derivable-from-rule | **yes** | identical output shape on mouse (501 k-mers) |
| 14 | `REFERENCE_FAMILY_OVERLAP` (2) | irreducibly-curated | **yes, by design** | it is the mouse gate; see §8 |
| 15 | `MIN_RMSK_FRACTION` = 0.99 | **invented** | no | see §8 |

**One item actively harms mouse (§7). One is a gate defined by the thing it gates (§8). The rest
either transfer with evidence or have a mouse counterpart derived independently.**

---

## 1. `refdata/repbase18.05-drop-exact.txt` — irreducibly-curated, fit-to-hg38

12 RepBase headers excluded from the repeat portion. Obtained by diffing the generated index
against the hg38 reference until the count reached 1,224/1,224. **Pure fitting**: no source states
that `U1_snRNA_Homo_sapiens` should be dropped.

The underlying *rule* — "drop from the repeat portion whatever Gencode/MASTER_FILELIST supplies" —
is derivable, and `475.11` applied it to mouse from scratch rather than transferring the list. That
was the right call: transferring would have been wrong in both directions. Mouse needed `U1 U2 U4B
U6 U7 U7_snRNA_Vertebrata U8 U14 NR_046235.1` (9 entries, of which only `U7_snRNA_Vertebrata` and
`NR_046235.1` appear in the human list), and `U5B1` — which human drops — had to be decided twice.

**Mouse validation:** the M-8 check that no family appears in both the RepBase and Gencode blocks
beyond the reference's own overlap. Currently `{SNORD}` against hg38's `{SNORD, YRNA}`.

## 2. `GENE_NAME_FAMILY_RULES` — irreducibly-curated, partly fit-to-hg38

14 ordered regex→family rules. The *patterns* are human gene-symbol conventions (HGNC `RNU5A`,
`SNORD13`, `MT-TF`), which is domain knowledge rather than reference-fitting; the *ordering* and
the special cases (`RNVU1-`→`RNU1`, ATAC before the generic `RNU\d+`) were tuned against hg38.

Until 2026-08-12 this was the worst hidden case of hg38-fitting in the codebase: the rules matched
literally, so mouse's title-case symbols (`Rnu5g`) never matched and **13 of 14 rules were dead on
mouse**, resolving 22 rows where human got 755. Fixed by `-4ee` (uppercase before matching);
verified a strict no-op for human — 0 of 62,629 hg38 gene names change family.

**Mouse validation:** 18 parametrised mouse cases in `test_master_filelist.py`, plus the tier
totals (mm10 173, mm39 177) being the right order of magnitude against human's 755.

## 3. Chromosome allowlists — derivable-from-rule

hg38's 674-scaffold list was **derived from the scaffold set of the reference BED** — fit-to-hg38.
The rule behind it ("the assembly release the annotation tracks were built against") is derivable,
and mouse uses it directly: `mm10.chrom-allowlist.txt` and `mm39.chrom-allowlist.txt` are UCSC
`chrom.sizes` for the base release (66 and 61 sequences, `475.7`).

**Mouse validation:** applied to RepeatMasker the rule drops 186,003 rows on 173 mm10 patch
scaffolds and 0 on mm39 — and mm39's 61 names are *exactly* its rmsk scaffold set, an independent
confirmation the list is neither over- nor under-inclusive.

## 4. `GENCODE_EXCLUDED_FAMILIES` — irreducibly-curated, fit-to-hg38

`{RNA5S, RNA5-8S, MTTRNA, MTRNR1, MTRNR2}`, obtained by observing that 582 of the 591
MASTER_FILELIST transcripts the reference BED omits belong to these 5 families. Fit to the answer,
though with a coherent biological reading: multi-copy rRNA and mitochondrial transcripts belong to
their repeat family, not to a "unique" genomic locus (cf. the `RNA45S` `rRNA_extra_hash` handling
in the mapper perl).

**Mouse:** 5/5 present on mm39; **4/5 on mm10, which has no `RNA5-8S` family**. The exclusion is a
no-op there rather than a silent error, but it means mm10's BED cannot be checked against this rule
the way mm39's can.

## 5. MASTER_FILELIST row order — derivable-from-rule

`read_in_filelists` assigns `priority_n` in file order, so order is behaviourally significant. It
would be easy to treat "reproduce the reference's order" as the requirement — that would be
fitting. Instead the requirement was made behavioural: four tests reimplement the priority
assignment and pin that priority is file order across blocks, that pipe-joined ids each consume
one, and that **permuting rows inside a family cannot move a count between families** — it only
changes which element represents the family.

`refdata/gencode-family-order.txt` (32 families) *is* taken from the reference, and is
fit-to-hg38 — but the priority comparison in the mapper is keyed by family, so between-family order
changes no count. It exists to keep generated files diffable. **Mouse inherits it harmlessly**;
families it does not name are appended in first-seen order.

## 6. `repbase18.05-classes.txt` — derivable-from-rule, fit-to-hg38

76 class tokens, derived by taking each of the 1,224 hg38 reference family names and subtracting
NAME and TAXON from the header. Fit to hg38 by construction.

**Mouse validation, and it is strong:** the vocabulary is a property of the RepBase 18.05 header
grammar, not of human. Applied to the mouse library it yields **1,071 families with 0 unparsed**
(human: 1,224 with 0 unparsed). `select_families` treats a non-empty unparsed list as fatal, so a
class token missing for mouse could not pass silently.

## 7. `RMSK_AMBIGUOUS_REPNAMES` — irreducibly-curated, fit-to-hg38, **and it costs mouse**

This is the finding P-4 exists to produce.

`RMSK_REPNAME_TO_FAMILY` maps 21 rmsk small-RNA repNames to families. Two names, `U5` and `U17`,
map to more than one family in human, so they are listed in `RMSK_AMBIGUOUS_REPNAMES` and
deliberately fall through to tier 3 (`gene_name`). For human that costs nothing:

| | transcripts overlapping rmsk `U5` | resolved by tier 3 | dropped |
|---|---|---|---|
| hg38 | 31 | **31** (`RNU5A-4P`→`RNU5A`, `RNU5D-1`→`RNU5D`, …) | 0 |
| mm39 | 11 | **1** (`Rnu5g`→`RNU5G`) | **10** |

Human gene symbols encode the subfamily, so the fallthrough always lands. Mouse names the same
genes `Gm24043`, `Gm23102`, `Gm22365` … — tier 3 returns nothing, Rfam has no U5 family, and the
transcripts are dropped. **Ten mouse U5 snRNA transcripts are lost to a rule whose justification
does not exist in mouse**: mouse has exactly one U5 family (`RNU5G`), so `U5` is not ambiguous
there at all.

This is the measured cause of the thinness recorded in `-9j2` — mouse U5 resting on a single
transcript against human's 34 — and it is 10 of mm39's 569 unresolved small-RNA transcripts.

**Strategy for mouse:** make the ambiguity assembly-aware rather than global — if an assembly's
Gencode supplies exactly one family for an ambiguous repName, that repName is not ambiguous for
that assembly and tier 2 should use it. mm39 would map all 11 U5-overlapping transcripts to
`RNU5G`; hg38 is unchanged, since it genuinely has five U5 families.
**Validation:** hg38 tier-2/3/4 totals must not move (4,104 / 755 / 298), and mm39's unresolved
residue must fall by exactly the 10 U5 and 4 U17 transcripts. Tracked in `-9j2`.

The other 19 map entries carry no mouse penalty, measured: of mm39's 569 unresolved small-RNA
transcripts, **554 overlap no rmsk small-RNA locus at all** and **zero** overlap a mouse-specific
repName the map lacks (`BC1_Mm`, `4.5SRNA`). The map's human origin is not, by itself, a mouse gap.
Six entries (`ACA64`, `U11`, `U12`, `U4atac`, `U6atac`, and `U13_` on mm39) match no locus in *any*
assembly, hg38 included — dead weight, not a species issue.

## 8. The two constants that gate mouse and were defined from hg38

`REFERENCE_FAMILY_OVERLAP = {SNORD, YRNA}` in `validate_refdata_set.py` is the M-8 check that mouse
does not invent a repeat/Gencode family overlap human does not have. It is fit-to-hg38 **by
design** — it is a statement about the reference — but it is worth being explicit that a mouse
build is being judged against a human artifact's blemish. If the hg38 reference were rebuilt, this
constant would have to be re-measured.

`MIN_RMSK_FRACTION = 0.99` is not curated or derived: **it is invented**, from three observations
(hg38 .9981, mm10 .99867, mm39 .99869). It caught the T-07 BEDs at 0.744 with enormous margin, so
nothing rests on its exact value — but it has no authority beyond "well below everything seen and
well above the known-bad case", and it should not be tightened without a reason.

## 9. What is genuinely derivable and needs no mouse decision

- **`*.rfam-family.tsv`** (`475.41`) — built by `build_rfam_family_table.py` from RNAcentral +
  Rfam. Explicitly not reference-derived. Validated on hg38 at 255/255 correct, and it carries far
  more weight for mouse (973/1,208 rows) than for human (298), precisely because mouse gene symbols
  are uninformative.
- **`*.rmsk-class-family.tsv.gz`, `*.rmsk-smallrna.bed.gz`** — cuts of UCSC `rmsk.txt.gz`,
  regenerated per assembly. Their *drift* is the dominant P-5 residual (83 family reassignments),
  which is a property of the source, not of a curated choice.
- **`DROP_CLASSES` + `KEEP_PREFIX`** — "drop RepBase class Simple_Repeat/tRNA/rRNA, except
  `MamSINE1_tRNA_`". The exception is a single hg38-derived special case; validated for mouse, whose
  library contains `MamSINE1_tRNA_Mammalia` and 113 records dropped by the class rule.
- **Index structural constants** — `SIMPLE_REPEAT_MAX_K=6`, `SIMPLE_REPEAT_LENGTH=60`,
  `TRNA_PAD/FLANK/CCA`, `RRNA_ANCHOR_K=40`. Read off the hg38 index, but they describe a
  construction, and the construction reproduces on mouse: 501 synthesized k-mers on all three
  assemblies, and the tRNA plain+flanked pairing holds.

---

## Conclusion

Of the 15 items, **9 are fit-to-hg38 in whole or in part.** Eight of those have since been either
validated against mouse by an independent check (§3, §6, §12, §13), given a mouse-specific
counterpart derived by rule rather than transferred (§1, §9), or shown to be harmless (§5, §7's
remaining 19 entries).

**One is not: `RMSK_AMBIGUOUS_REPNAMES`.** It is the only curated constant measured to actively
degrade the mouse artifacts, and it does so silently — the transcripts are omitted, not
mis-assigned, so nothing downstream reports a problem. It is exactly the failure mode this phase
was created to find: a rule that reproduces hg38 perfectly, is justified entirely by a human
property, and costs mouse real annotation.
