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

The hardest decision in the file, because the original values were curated by hand. Three tiers,
tried in order, measured on hg38's 5,261 reference rows:

| Tier | Mechanism | Function | Resolved |
|---|---|---|---|
| 1 | `--family-override` TSV | `read_family_overrides` | (operator-supplied) |
| 2 | rmsk small-RNA overlap ≥ 50% | `overlapping_repname` | 3,932 |
| 3 | `gene_name` pattern rules | `family_from_gene_name` | 831 |
| — | neither | — | **473 unresolved** |

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

**Decision on the 473 (`475.41`): omit and count, never guess.** Column 4 is the counting label,
so a plausible-but-wrong family is worse than an absent row. Filling them with the Gencode biotype
was rejected — that is precisely the defect T-09 exists to fix. They are logged at build time and
can be supplied through `--family-override`. Their only Gencode annotation is `gene_type "snoRNA"`
with a clone-style `gene_name` like `AC020634.1`, so no rule over the GTF can recover them; Rfam
or snoDB is the likely source.

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
  - strips a trailing `?` (`DNA?` is a provisional RepeatMasker call; the reference has the settled
    label), and
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

It doubles as the **cross-species** source: pass the *human* file to a mouse build and 215 mouse
families inherit by exact family name. Safe because 167 of 168 shared names carry the identical
RepBase class token in both species (only `UHG` differs, and its mouse token is empty).

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
| `475.41` | 473 hg38 Gencode rows resolved by neither rmsk nor gene_name (229 SNORA, 92 SNORD, 54 YRNA, …) |
| `475.43` | Mouse `B1`, `B1-DID`, `B2` — the major mouse SINEs. UCSC has no bare `B1`/`B2` repName (instances are `B1_Mus1`, `B2_Mm1a`, …), and prefix matching is the rejected rule above |
| — | 147 mm10 / 149 mm39 RepBase families unresolved. Every mouse family with a nonzero genomic footprint now resolves; the remainder have zero footprint, `B1`/`B2` excepted |
| — | 99 hg38 repeat rows differ from the reference through post-2020 RepeatMasker reclassification. Accepted on purpose (rmsk outranks the curated table) |

## 10. Fidelity against the reference

Per block, on shared column 1, after every rule above:

| Block | shared | identical on all 5 columns |
|---|---|---|
| Gencode | 4,633 | 4,633 |
| rRNA | 15 | 15 |
| tRNA | 864 | 864 |
| SimpleRepeat | 501 | 501 |
| `(XXX)N` | 14,218 | 14,218 |
| miRNA | 3,702 | 3,700 (2 = miRBase drift) |
| repeat | 1,718 | 1,619 (99 = RepeatMasker drift) |

Gencode family block order matches the reference exactly.
