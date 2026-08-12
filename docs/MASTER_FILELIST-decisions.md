# MASTER_FILELIST: every decision, and why

> Companion to `docs/ROADMAP-refdata.md`. The function docstrings in
> `bin/python/refdata_generation/generate_master_filelist.py` carry the same rationale in short
> form and cross-reference the section numbers here. Keep the two in step.

Decision record for `bin/python/refdata_generation/generate_master_filelist.py` (T-09, `-gz2`).
Each rule below was measured against the hg38 reference before being adopted; each rejected
alternative was measured too, and the number that killed it is recorded. Function docstrings
carry the same rationale in short form — this file is the long form.

Reference file:
`examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv`
(26,422 rows, built 2020-12-03).

---

## 1. What the file is

**It is a concatenation of source lists, not a table.** Column 5 records which list each row came
from — `genelists.RNU1`, `hg38-tRNAs.fa.list.wgenome_flank.flankNs`, `hsa-mir-6859-1`. The
original was assembled by concatenating per-family gene lists from a directory that no longer
exists; the surviving evidence is the commented-out call at
`bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl:43`:

```perl
#&read_in_filelists(".../RNA_type_analysis/genelists.RNA5S","RNA5S");
```

Two consequences follow, and both are enforced by tests:

- **Never deduplicate by column 1.** 61 names appear twice with *different* annotations —
  `DNA1_MAM` is `TcMar` in the RepBase block (line 5880) and `TcMar-Tc1` in the rmsk block
  (line 20674). The previous generator deduplicated globally and destroyed this.
  → `test_reference_repeats_column_one_so_dedup_would_be_wrong`
- **Never sort.** Row order is assigned as `priority_n` by `read_in_filelists`. See §6.

### Which column matters

`read_in_filelists` parses each row as
`($allenst, $allensg, $gid, $type_label, $typefile)`, and **column 4 is `$type_label`** — the
family a read gets counted under, since `print_output` accumulates `$count{$ensttype_join}++`.
Columns 2, 3 and 5 are descriptive. Getting column 4 wrong changes results; getting column 5
wrong does not.

---

## 2. Build order

The index selects its Gencode portion by MASTER_FILELIST membership, so:

```
generate_master_filelist.py   →   generate_bowtie2_index.py
```

`generate_master_filelist.py` takes `--repbase-species-fasta` directly rather than the index's
`.repbase_provenance.tsv` sidecar. Both call `repbase.select_families()`, so the RepBase block
matches the index exactly, and the filelist never has to wait for an index to exist.
**Decision:** avoid the circular dependency rather than document a two-stage build around it.

---

## 3. Block layout and `main()` call order

`main()` emits six blocks in this order. hg38 counts shown.

| # | Block | n | Builder | Columns |
|---|---|---|---|---|
| 1 | Gencode + rRNA | 5,276 | `build_gencode_rows`, `build_rrna_rows` | ENST │ ENSG │ gene_name │ **FAMILY** │ `genelists.FAMILY` |
| 2 | RepBase families | 1,224 | `build_repbase_rows` | NAME │ FAM │ FAM │ FAM │ CLASS |
| 3 | tRNA | 864 | `build_trna_rows` | name │ anticodon │ anticodon │ `tRNA` │ source list |
| 4 | SimpleRepeat k-mers | 501 | `build_simple_repeat_rows` | `AT_SimpleRepeat` │ `Simple_repeat` ×4 |
| 5 | miRNA | 3,765 | `build_mirna_rows` | `MI0022705` │ `miRNA` ×3 │ mirbase name |
| 6 | rmsk leftovers | 14,220+ | `build_rmsk_rows` | NAME │ FAM │ FAM │ FAM │ CLASS, and `(XXX)N` |

Blocks 1–4 are exactly the index's contents plus the 259 filelist-only scaffold transcripts the
index drops. Blocks 5–6 are filelist-only.

Within block 1, Gencode and rRNA rows are merged and re-grouped by family into the order given by
`refdata/gencode-family-order.txt` (`read_family_order`). The reference's Gencode region is 32
contiguous family runs — one per `genelists.<FAMILY>` source file — and the rRNA runs
(`RNA18S`, `RNA28S`, `RNA45S`) sit inside that sequence between `MTTRNA` and `RNA5S`.

---

## 4. Gencode family assignment (`build_gencode_rows`)

The hardest decision in the file, because the original values were curated by hand. Four tiers,
tried in order:

| Tier | Mechanism | Function | hg38 | mm10 | mm39 |
|---|---|---|---|---|---|
| 1 | `--family-override` TSV | `read_family_overrides` | (operator-supplied) | | |
| 2 | rmsk small-RNA overlap ≥ 50% | `overlapping_repname` | 4,104 | 1,993 | 1,972 |
| 3 | `gene_name` pattern rules | `family_from_gene_name` | 755 | 173 | 177 |
| 4 | Rfam family via RNAcentral | `read_rfam_families` | 298 | 973 | 1,208 |
| — | none of them | — | **559** | 785 | 569 |

The mouse tier-3 numbers were 22 and 22 until 2026-08-12 (`-4ee`). The rules are written against
human symbols (`RNU5A`, `SNORD13`, `MT-TF`) and were matched literally, so every one of them was
dead against mouse's title-case symbols (`Rnu5g`, `Snord13`, `mt-Tf`) except the single `^Mt-t`
rule that carried `re.I`. `family_from_gene_name` now uppercases before matching, which is a
verified no-op for human — 0 of hg38's 62,629 gene names change family — and resolves 177 mouse
names that previously fell through to Rfam or to nothing. It is what gives mouse the RNU5G,
RNU4ATAC, RNU11, RNU12, MTRNR1 and MTRNR2 families at all, and RNU5G is what settled `U5B1` in
`-475.11`.

**Tier 2 detail.** `refdata/<asm>.rmsk-smallrna.bed.gz` holds rmsk loci whose `repClass` is one of
`srpRNA, scRNA, snRNA, tRNA, rRNA, RNA`. A transcript takes the `repName` of the feature covering
≥ 50% of it, mapped through `RMSK_REPNAME_TO_FAMILY`. The mapping is effectively a function — 19
of 21 observed repNames are unambiguous (`U6`→`RNU6`, `7SLRNA`→`RN7SL`, `5S`→`RNA5S`,
`7SK`→`RN7SK`, `HY1`/`HY3`→`YRNA`). The two that are not, `U5` (splits RNU5A/B/D/E/F) and `U17`
(SNORA vs RNU105), are listed in `RMSK_AMBIGUOUS_REPNAMES` and deliberately fall through to tier 3.

**Why a shipped BED digest rather than `rmsk.txt.gz`:** the full table is 155 MB; the small-RNA
subset is 12,753 rows / 125 KB and is all this rule needs.

**Tier 3 detail.** `GENE_NAME_FAMILY_RULES` is an ordered list and **the order is load-bearing**:
`^(RNU\d+)` would swallow `RNU6ATAC`, so the ATAC rules precede it.
→ `test_rnu6atac_is_not_swallowed_by_the_rnu6_rule`

**Tier 4 detail (`475.41`).** Rfam names its families after the box class — `SNORA70`,
`SNORD56`, `SCARNA20` — which is exactly the distinction the clone-named snoRNA rows are missing.
`refdata/<asm>.rfam-family.tsv` is built by `build_rfam_family_table.py`, which streams
RNAcentral's `ensembl.tsv` (transcript → URS) and `rfam.tsv` (URS → Rfam accession) plus Rfam's
`family.txt.gz` (accession → Rfam ID), then maps the ID through `RFAM_PREFIX_TO_FAMILY`.

Validated against the 365 residue transcripts the 2020 reference does label: **255 predicted, 255
correct, zero wrong.** Unlike the curated RepBase table this is **not** reference-derived, so it is
not circular — and it transfers to mouse, where it matters far more: mouse gene symbols barely
match the human `RNU`/`SNORD` patterns, so tier 3 fires on only 22 rows and tier 4 carries
1,075 (mm10) and 1,337 (mm39).

Three Rfam IDs are listed as `AMBIGUOUS` and deliberately predict nothing: `snoU2_19` and
`snoU2-30`, which the reference itself splits between SNORD and SCARNA, and `Vault`, which cannot
choose between VTRNA1/2/3. `SNORA73` (U17) also splits in the reference, 5 SNORA / 2 RNU105, but
both RNU105 rows carry gene_name `RNU105C` and are settled at tier 3, so Rfam is never consulted
for them — which is why tier order is what makes the validation above the right comparison.

**Decision on what tiers 2–4 still miss: omit and count, never guess.** Column 4 is the counting label,
so a plausible-but-wrong family is worse than an absent row. Filling them with the Gencode biotype
was rejected — that is precisely the defect T-09 exists to fix. They are logged at build time and
can be supplied through `--family-override`. What remains after tier 4 has no Rfam family at all,
so no source currently to hand can recover it.

**Selection is a side effect of resolution.** A transcript is in the filelist iff it resolves to a
family. There is no separate type filter — the discredited `DEFAULT_TRANSCRIPT_TYPES` approach
(killed for the index in `475.35`) is not used here either. `SMALL_RNA_GENE_TYPES` exists **only**
to decide whether an unresolved transcript is worth warning about; it never selects.

---

## 5. Repeat annotation (`resolve_repeat`)

Shared by blocks 2 and 6. Resolution order, most authoritative first:

```
1. RMSK_REPNAME_TO_FAMILY_UPPER   curated small-RNA map   →  'small_rna'
2. assembly rmsk table            read_rmsk_class_family  →  'rmsk'
3. curated RepBase table          read_curated_class_family → 'curated'
4. name aliases                   repeat_name_aliases     →  '..._int' / '..._stem'
5. fall back to the name itself                           →  'unresolved'
```

`repeat_row` wraps it into `NAME │ FAM │ FAM │ FAM │ CLASS`. The `how` string is tallied per build
so the source mix is visible in the log.

**Step 1 — why a curated small-RNA map at all.** The reference has `5S` as `RNA5S` in every
annotation column, never rmsk's `rRNA`/`rRNA`. These repNames *are* the curated families.

**Step 2 — class and family are UCSC `repClass`/`repFamily`, not RepBase.** `L2B_CR1_Eutheria`
carries RepBase class `CR1` but the reference says family `L2`, class `LINE`. This killed the
original T-09 premise that the class could be parsed from the RepBase header.
`read_rmsk_class_family` also
  - strips a trailing `?` (`DNA?` is a provisional RepeatMasker call) — **but not for the reason
    originally recorded here.** The old justification was "the reference has the settled label";
    P-5 (`475.22`) measured it and it is false. The reference filelist has **53** rows whose family
    ends in `?`, and the `?` does reach the output: the dedup perl strips it at `read_in_filelists`,
    but the *mapper* perl strips only a trailing `_`, so the repeat arm reports `ERVL?` as its own
    element, exactly as the reference outputs do.

    Stripping survives because it is measurably the closer approximation, not because it is
    faithful. The real problem is upstream — **71 repNames in the modern UCSC table carry
    conflicting families, 40 of them a `?`/non-`?` pair for the same name** (`EULOR5A` is both
    `Crypton-A` and `DNA?`), so `first occurrence wins` picks by sort order. Measured against the
    reference:

    | policy | reassignments | reference-has-`?` | regen-has-`?` | genuine drift |
    |---|---|---|---|---|
    | strip `?` | **93** | 10 | 0 | 83 |
    | keep `?` | 112 | 8 | 21 | 83 |

    The 83 genuine reclassifications (`DNA`→`Crypton-A`, `hAT`→`hAT-Ac`) are invariant under either
    policy: they are era drift between the 2020 reference and today's rmsk, and only an era-matched
    snapshot would remove them. Stripping collapses the conflicting pairs and removes the arbitrary
    choice; keeping the `?` reproduces 2 more of the reference's own `?` rows at the cost of 21 new
    wrong-direction ones. Revisit only with a 2020-era `rmsk.txt.gz`.

    Also,
  - indexes UCSC's slash-qualified names by their bare prefix, because UCSC writes alpha satellite
    as `ALR/Alpha` while the index and filelist use `ALR`.

The `examples/inputs/<asm>.repeatmasker.tsv.gz` files are stripped to `gene_id` and **cannot**
supply class/family, which is why `refdata/<asm>.rmsk-class-family.tsv.gz` exists.

**Step 3 — the curated table (`475.42`).** `refdata/hg38.repbase-class-family.tsv` holds all 1,224
curated pairs lifted from the reference's RepBase block. It exists because 354 of those families
have no genomic instances in the modern rmsk track, so column 4 would otherwise fall back to the
family's own name — splitting alpha-satellite variants out of `centr` instead of pooling them.

Two deliberate properties:
- **It is consulted after rmsk**, so current RepeatMasker calls win and post-2020 reclassification
  (99 hg38 families) stays visible instead of being frozen at 2020.
  → `test_rmsk_wins_over_the_curated_table`
- **It is circular** — derived from the reference — so it must be excluded when measuring how well
  an autogenerated hg38 filelist reproduces that reference (`475.22`). Same derivation style, and
  same caveat, as `refdata/repbase18.05-classes.txt` and `repbase18.05-drop-exact.txt`.

It doubles as the **cross-species** source: pass the *human* file to a mouse build and mouse
families inherit by exact family name. Safe because 167 of 168 shared names carry the identical
RepBase class token in both species (only `UHG` differs, and its mouse token is empty).

`--repbase-class-family` is repeatable and **later files win**, so mouse passes the human table
first and then `refdata/mm.repbase-class-family.tsv` (`475.43`), 36 mouse-specific rows led by
`B1 → Alu/SINE` and `B2 → B2/SINE`. Those two are the major mouse SINEs and had resolved to
nothing, because UCSC has no bare `B1`/`B2` repName — the instances are `B1F`, `B1_Mus1`,
`B2_Mm1a`. B1's subfamilies alone carry 108 Mb.

That table was derived **parent → children**, the opposite direction from the longest-prefix rule
rejected above: for an unresolved family X, collect every rmsk repName that *extends* X and adopt
their annotation only if all of them agree; disagreement rejects the row. On hg38 ground truth the
same rule scores 44/47 with one subfamily and 28/29 with two — good, but not good enough to run
automatically, so it is used as **evidence for a curated table, never as a fallback**. Every row
was additionally cross-checked against the family's own RepBase class token (`ERV2`→LTR,
`SINE1/7SL`→SINE, …): 36 rows, zero conflicts. Columns 4 and 5 of the file record the subfamily
count and their genomic bases so any row can be audited.

**Step 4 — two systematic name aliases** (`repeat_name_aliases`), each validated before adoption:

| Alias | Rationale | Recovered | Validation |
|---|---|---|---|
| `NAME_I` → `NAME_I-int` | RepeatMasker's LTR internal-segment convention | 8 mm10 | 8/8 agree with the family's RepBase class token; 1 human ground-truth case, exact |
| `_MM`, `_HS` stripped | RepBase tags the species in the family name | 27 mouse | 27/27 agree with the RepBase class token, 0 contradictions |

Only those two species tags. `_LTR`, `_DNA`, `_II`, `_MAM` are class tokens *inside* the family
name (`ERVB3_1-LTR_MM`) and stripping them would corrupt it.
→ `test_class_tokens_in_names_are_not_treated_as_species_tags`

The exact name is always tried before any alias.
→ `test_exact_name_beats_every_alias`

---

## 6. Row order

**Between-family order is reproduced; within-family order is not, and provably need not be.**

`read_in_filelists` assigns `priority_n` in file order, but the comparison that uses it,
`$enstpriority < $read_hash{$r1name}{flags}{$ensttype}`, is keyed **by family**. So priority only
chooses which element is named as a family's representative (`master_enst`) in the SAM output.
Counts are accumulated on the family join, `$count{$ensttype_join}++`, so **no count and no fold
change depends on within-family order.**

That matters because within-family order is not recoverable. Measured against the reference, the
number of the 29 Gencode families whose rows follow a candidate order:

| Candidate order | Families matching |
|---|---|
| `parsed_ucsc` file order | 7 / 29 |
| `gene_name` lexicographic | 6 / 29 |
| `gene_name` natural sort | 6 / 29 |
| ENST id | 11 / 29 |
| ENSG id | 9 / 29 |

It came from the internal order of the `genelists.*` files. Between-family order *is* reproduced,
via `refdata/gencode-family-order.txt`, which keeps generated files diffable against the reference.
→ `test_family_order_file_matches_the_reference_block_order`

The RepBase block order is likewise not reproducible: it is close to rmsk abundance-descending but
not equal to it (position 6 breaks monotonicity), and matches neither rmsk first-occurrence order
nor the RepBase FASTA order.

---

## 7. Per-block decisions

**rRNA (`build_rrna_rows`, `rrna_copy_suffix`).** Column 3 is `RNA18SN5` — the rDNA copy
identifier — read from the GenBank `/gene="RNA45SN5"` qualifier, *not* a positional index. Rows
are grouped 18S, then 28S, then 45S, matching the reference.

**tRNA (`build_trna_rows`).** Each locus is followed immediately by its `_withgenomeflank` twin,
in gtRNAdb file order — verified to reproduce the reference block exactly. Columns 2 and 3 are the
anticodon-level name (`tRNA-Ala-AGC`), column 4 is the constant `tRNA`, column 5 is
`<gtrnadb file>.list.wgenome_flank.flankNs`. Source is the gtRNAdb FASTA, **not** the assembly
tRNA track (`475.36`).

**SimpleRepeat (`build_simple_repeat_rows`).** The same 501 synthesized canonical k-mers the index
carries, taken from `simple_repeat_fasta()` so the two artifacts cannot drift apart. Deriving them
from the RepeatMasker track was the defect that gave mm39 17,181 rows against the correct 501.

**miRNA (`build_mirna_rows`, `read_mirna_gff3`).** Only `miRNA_primary_transcript` features. Each
emits two rows, the entry and a `-proximal` twin.

**RepBase (`build_repbase_rows`).** Family set comes from `repbase.select_families()`, identical to
the index. `--repbase-drop-exact` overrides the drop list; mouse needs
`refdata/repbase18.05-drop-exact.mm.txt` because `U7`, `U8` and `U14` each appear twice in the
mouse file, once well-formed and once as an unparseable bare header. Unparseable families are
**fatal**, never silently dropped.

**rmsk leftovers (`build_rmsk_rows`, `read_repeatmasker_names`).** Everything in the assembly track
the RepBase block did not already cover, in first-occurrence order, uppercased. `(XXX)N` simple
repeats are emitted verbatim with constant `Simple_repeat` annotation.

---

## 8. Rejected alternatives

| Alternative | Why rejected |
|---|---|
| Parse repeat class from the RepBase header | `L2B_CR1_Eutheria` is family `L2`, class `LINE`, not `CR1` |
| Deduplicate rows by column 1 | 61 reference names carry two different annotations |
| Longest-prefix inheritance for unresolved repeats | 18 of 161 matches wrong on hg38: `MER21`/`MER22`/`MER25` all collapse onto `MER2`, giving `TcMar-Tigger/DNA` where the truth is `ERVL/LTR`, `centr/Satellite`, `L1/LINE` |
| Fill unresolved Gencode families with the biotype | Reintroduces the exact defect T-09 fixes, in the counting column |
| Dfam API for unresolved families | `dfam.org` does not resolve from TSCC (NCBI, UCSC, gtRNAdb, EBI all do) |
| Select Gencode rows by transcript type | Discredited in `475.35`; the reference spans types no plausible set admits |
| Take the RepBase family set from the index provenance sidecar | Creates a circular filelist↔index dependency for no gain |

---

## 9. Known residue

| Issue | Residue |
|---|---|
| `475.41` | **Closed.** Tier 4 (Rfam) cut the hg38 residue 857 → 559 and the mouse residue 1,907 → 833 (mm10) and 1,928 → 592 (mm39). What is left has no Rfam family at all |
| `475.43` | **Closed.** `B1` and `B2` resolve through `refdata/mm.repbase-class-family.tsv`. `B1-DID` does not: no rmsk name extends it, so there is no evidence to curate from |
| — | 112 mm10 / 114 mm39 RepBase families unresolved, none with a populated UCSC subfamily group left to derive from |
| — | 99 hg38 repeat rows differ from the reference through post-2020 RepeatMasker reclassification. Accepted on purpose (rmsk outranks the curated table) |

## 10. Fidelity against the reference

Per block, on shared column 1, after every rule above:

| Block | shared | identical on all 5 columns |
|---|---|---|
| Gencode | 4,890 | 4,890 |
| rRNA | 15 | 15 |
| tRNA | 864 | 864 |
| SimpleRepeat | 501 | 501 |
| `(XXX)N` | 14,218 | 14,218 |
| miRNA | 3,702 | 3,700 (2 = miRBase drift) |
| repeat | 1,718 | 1,619 (99 = RepeatMasker drift) |

Gencode family block order matches the reference exactly.

Building the hg38 index **from the generated filelist** rather than the reference one gives
`shared=7494, missing=112, extra=267, mismatching=0` — every shared record sequence-identical, so
the whole divergence is Gencode membership. Tier 4 cut `missing` from 369 to 112. The 267 extra are
real Gencode small RNAs the 2020 curated lists omitted; `475.22` must account for them separately.
