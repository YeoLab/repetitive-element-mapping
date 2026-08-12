# Reference Data Roadmap — hg38 provenance → mm10/mm39 generation

**Status as of 2026-08-10.** Tracked in beads under `repetitive-element-mapping-475`.
Run `bd show repetitive-element-mapping-475` for live status; this document is the rationale,
the evidence, and the phase structure.

## Why this document exists

The goal is mouse reference data (mm10/mm39) for the eCLIP repetitive element pipeline. Mouse
has **no reference outputs to validate against**, so mouse correctness cannot be established
directly. It can only be inherited — from a chain of evidence built on human, where reference
outputs *do* exist.

That chain is seven phases. Every phase exists to make the next one interpretable:

| # | Phase | Why the next phase needs it |
|---|---|---|
| 0 | Verify the Snakemake pipeline is equivalent to the CWL/Perl one | Everything downstream assumes the engine computes the same thing; it was never tested |
| 1 | Run the pipeline on stock hg38 refdata, reproduce the reference outputs | Establishes that the *pipeline* is correct, so later discrepancies can be blamed on refdata |
| 2 | Document the structure of the four hg38 artifacts | You cannot regenerate a file whose schema and invariants are unstated |
| 3 | Identify the upstream sources that can regenerate each artifact | Names what regeneration is even made of |
| 4 | Gap analysis: auto-derivable vs manually curated | Isolates the parts that *cannot* transfer to mouse |
| 5 | Run the pipeline on regenerated hg38 refdata, compare to phase 1 | **The last checkpoint against ground truth.** Any refdata defect surviving here is invisible forever after |
| 6 | Generate mouse refdata by the same method | The deliverable |

The failure mode this structure guards against: regenerate mouse refdata, get plausible-looking
output, and have no way to know it is wrong. Phase 5 is the load-bearing step — it is the last
moment a generator bug meets a known-correct answer.

Phase 0 was added on 2026-08-10 after the first end-to-end comparison of the two engines found
a defect that silently corrupted 21% of the output rows. The lesson generalizes: **a translated
pipeline needs equivalence testing against the original, not smoke testing.** Both defects found
so far ran to completion and produced plausible numbers.

## Current status

| Phase | Status | Evidence |
|---|---|---|
| 0 | **PASS on SE and PE.** Three defects found and closed | Read counts reproduce the CWL references exactly (SE 182/182, PE 169/169); residual ≤1.6e-11 (`475.16`, `475.23`, `475.29` all closed) |
| 1 | Blocker fixed, full SE rerun reproduces the reference | Zero repeat-family reads traced to bowtie2 missing from PATH + a swallowed exit code; fixed (`475.16`) |
| 2 | Partial | Knowledge exists as defect narratives in `.forge/stages/2-architect/notes/`, not as a schema (`475.18`) |
| 3 | Essentially done, unrecorded | RepBase 18.05 source confirmed at 100% coverage; needs a source manifest (`475.19`) |
| 4 | Partial | Curated constants known but not enumerated as a class (`475.20`) |
| 5 | Not started | Blocked on phase 1 + both generators (`475.22`) |
| 6 | 0 of 4 mouse artifacts valid | All mm10/mm39 outputs predate a known fix (`475.6`–`475.15`) |

### Phase 1 blocker: zero repeat-family reads on full datasets

Measured 2026-08-10 across the four committed runs in `results/`:

| Run | AllReads | RepFamilyReads |
|---|---|---|
| `se_small` | 2,077 / 2,206 | 1,162 (0.559) / 1,784 (0.809) |
| `pe_small` | 1,941 | 1,444 (0.744) |
| `se_full` | 10,272,097 | **0 (0)** |
| `pe_full` | 3,954,802 | **0 (0)** |

Small datasets assign 56–81% of reads to repeat families. Full datasets assign **exactly zero**,
in both SE and PE.

**Root cause, found 2026-08-10 (`475.16`): bowtie2 was not on `PATH` for the full runs.**
`results/se_full/barcode1/mapped/rep.sam.bowtieout` — 67 bytes, present on disk, not deleted as
the earlier investigation assumed — reads:

```
stdbuf: failed to run command 'bowtie2': No such file or directory
```

Identical in all five full-run samples. The small runs' `.bowtieout` files contain real bowtie2
alignment summaries. The full runs were simply launched without `--use-conda`. **This is not
scale-dependent** — dataset size correlated only because the full and small runs were launched
differently.

The real defect is the *silent* failure. `map_repetitive_elements_{se,pe}.py` called
`proc.wait()` and discarded the return code, then wrote the `.done` file and exited 0. `stdbuf`
exits 127; the pipe reader saw an immediate EOF and treated it as a successful zero-alignment
run, so the empty `rep.sam` flowed through splitbam, dedup and combine untouched. An obvious
environment problem was thereby converted into a months-long misdiagnosis.

The reference data was never the cause. `examples/inputs/hg38/bowtie2_index/` is the authentic
downloaded index (7,606 headers, zero `::coords` suffixes) and `repeat_mapping_SE_full.yaml`
points at it plus the 2020 reference MASTER_FILELIST.

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
May 2026 as a *transient* bowtie2 failure whose stderr "was deleted by cwltool after the run."
Both halves were wrong: the stderr file was on disk the whole time, and the cause was a missing
binary, not a transient fault. Its one correct call was the mechanism — "Perl script exits with
code 0 (no error checking on pipe return value)" — which is exactly what was fixed.

**Fixed** on `fix/bowtie2-silent-failure`: `check_bowtie_exit()` in both mappers surfaces
bowtie2's own stderr and exits non-zero *without* writing `.done`, so downstream rules cannot
consume an empty `rep.sam` as a valid result. Four regression tests in `tests/workflow_scripts/`
fail against the pre-fix code and pass after.

### The reference outputs have no provenance

`test-provenance/` holds the SE and PE reference `.nopipes.tsv`/`.withpipes.tsv`, but
`configs/`, `manifests/`, `scripts/` and `reference-metadata/` are **all empty**. Nothing records
which inputs, which reference data, which tool versions, or which command produced them. That
`EXAMPLE_SE` corresponds to `INV_B` is inferred from matching read counts above — not documented.
Until that is pinned down, "reproduce the reference" is not a well-defined target (`475.17`).

## Phase 0 — is the Snakemake pipeline equivalent to the CWL one?

Every phase below assumes the Snakemake implementation computes what the CWL/Perl one
computed. That assumption was never tested end to end, and testing it found a second
output-corrupting defect. It is now its own phase, ahead of phase 1.

### Result: the conversion is sound, after three fixes

**Validated on both SE and PE.** All three defects are closed; the full suite is 45 tests.

Validated against a **fresh CWL 1.0.0 run** of the same dataset (`INV_B`), generated
2026-08-10 on the same machine — not against the 2020 artifact of unrecorded provenance in
`test-provenance/`. Baseline committed under `tests/cwl_baseline/ecliprepmap-1.0.0-SE/`.

| | `.nopipes.tsv` | `.withpipes.tsv` |
|---|---|---|
| Elements | 182 / 182 shared | 1,915 / 1,915 shared |
| Read counts | **exactly identical** | **exactly identical** |
| Max relative deviation (derived floats) | 4.6e-12 | 1.6e-11 |
| Verdict | **PASS** | **PASS** |

The `.parsed` files agree too: all four `#READINFO` totals and every `TOTAL` read count
identical (22,576,144 all / 17,691,900 usable / 8,214,042 genomic / 9,477,858 rep-family).

PE, against the `test-provenance` PE reference: **169 / 169 elements, read counts exactly
identical, max relative deviation 2.0e-13.**

**The SE and PE evidence are not equivalent.** The PE reference is `.nopipes.tsv`/`.withpipes.tsv`
only — those files carry no `#READINFO` lines, and no CWL PE `.parsed` exists anywhere. So PE's
four header fields have never been compared against the Perl implementation; only SE's have.
That matters because `475.29` was a defect in one of those fields that never reached the TSVs.
The PE workflow also has code the SE path lacks (two-barcode merge, PE pairing logic — the
r1/r2 flag-swap bug in `d44ca22` was PE-only), so a PE-specific defect of that class is not
hypothetical. Tracked as `475.31`.

### Full SE surface comparison

Every SE output was compared against the fresh CWL 1.0.0 run, not just the TSVs:

| Surface | Result |
|---|---|
| `.nopipes.tsv` | 182 / 182 elements, read counts **exact** |
| `.withpipes.tsv` | 1,915 / 1,915 elements, read counts **exact** |
| `.parsed` `#READINFO` | all 4 totals match |
| `.parsed` `TOTAL` | all 1,915 rows match |
| `.parsed` `ELEMENT` | all **181,342** rows match on id / readnum / enst / ensg |
| `rmDup.sam.gz` | 17,691,900 lines, identical **except column 5** |
| `preRmDup.sam.gz` | identical **except column 5** |

One difference remains, and it is a tool-version artifact rather than a translation defect:
bowtie2 2.2.6 (CWL) emits `MAPQ 1` on secondary alignments where 2.5.5 emits `MAPQ 255`,
confirmed by running both aligners on the same fastq and index. Blanking column 5 makes both
SAM files byte-identical after sorting, and no pipeline logic reads MAPQ — the mapper keys on
`AS:i`. Accepted deliberately; dependency versions are now pinned exactly in
`workflow/envs/dropin.yaml` (`475.32`).

Two things are *more* reproducible in the Snakemake version than in the original: row order
among equal-count elements (Perl uses hash iteration, which post-5.18 varies run to run), and
float precision (full float64 versus Perl's 15-digit truncation). Neither is a gap to close.

### The three defects, and what each one teaches

| Issue | Defect | Symptom it presented as | Reached the TSVs? |
|---|---|---|---|
| `475.16` | bowtie2 exit code discarded | `RepFamilyReads 0` — looked like bad reference data | Yes — everything zero |
| `475.23` | RPR truncated to 5 decimals | `inf` fold enrichment — looked like division by zero | Yes — 21% of rows |
| `475.29` | usable fraction as `usable/usable` | a plausible `1.0` | **No** |

The third is the instructive one. It produced a value nothing would flag, never reached the
TSVs at all (fold enrichment derives from the `TOTAL` rows), and survived every output
comparison. It was caught only by diffing `.parsed` **headers** against a CWL run — which is
why `tests/cwl_baseline/` keeps the `.parsed` `#READINFO` header alongside the TSVs, and why
`475.27` should compare both.

`475.30` records the same unchecked-pipe defect as `475.16` still present in the Perl mappers.
It is deliberately unpatched: the CWL resolves the script from `PATH` to the installed module
copy, not this repo's, so patching here would change nothing while perturbing the reference
implementation that equivalence testing depends on.

Read counts — the quantity the pipeline actually measures — reproduce the CWL reference
exactly. The residual ~1e-12 is floating-point noise from the two paths reaching the derived
columns by different arithmetic.

Row **order** is deliberately not compared. Elements with equal read counts are emitted in
Perl hash-iteration order, which is not reproducible across implementations (and, post-5.18,
not reproducible across *runs*). Ordering among ties carries no information.

### The defect this uncovered (`475.23`, P0)

`merge_parsed_files.py` wrote the RPR column and the `#READINFO` fractions with `:.5f`; the
Perl original writes a bare double. `calculate_fold_change_from_parsed_files.py` reads
`clip_rpr` straight out of that column, so the truncation landed in the pipeline's primary
output.

On `se_full` this was **data corruption, not a rounding nicety**:

- **39 of 182 elements (21%)** had `Input_clip_rpr` truncated to exactly `0.0`
- a zero denominator makes `Fold_enrichment` and `Information_content` `inf`
- where `IP_clip_rpr` also truncated to zero, `0/0` left the cells **empty**
- 86 destroyed cells; 636 float columns wrong by up to **37%** relative
- `5S-Deu-L2` fold enrichment: `1.4999` vs the correct `1.9228`
- `antisense_Crypton`: `12  0.0  5.0  0.0  <empty>  <empty>` vs
  `12  6.78e-07  5.0  2.82e-07  2.404  8.58e-07`

Fixed by removing the format specs. Three regression tests in
`tests/workflow_scripts/test_merge_parsed_files_precision.py` fail against the pre-fix code.

### Why this keeps happening

Both defects found so far — the swallowed bowtie2 exit code (`475.16`) and this truncation —
share a shape: **the pipeline kept running and produced plausible output**. Neither raised an
error; both were found only by comparing numbers against a reference. A translated pipeline
needs equivalence testing, not smoke testing, and `475.26` audits the remaining scripts for
the same two classes.

### Remaining phase-0 work

| Issue | Work |
|---|---|
| `475.24` | Fresh CWL v0.1.0 baseline run (EV136, hg19) — same-machine, same-day |
| `475.25` | Snakemake on identical v0.1.0 inputs, compared to that baseline |
| `475.26` | Audit remaining translated scripts for silent-failure/precision defects |
| `475.27` | Wire the equivalence comparison into the test suite |

The v0.1.0 example is **not** the repo's hg38 data: it uses EV136, Gencode v19,
`MASTER_filelist.wrepbaseandtRNA`, mirbase v20 hg19 and `RepeatMask.bed`. Comparing against it
requires pointing Snakemake at those same hg19 inputs.

Operational note: do not pipe `module load` (`module load x | head`) — the pipeline runs it in
a subshell and the `PATH` changes are lost.

## Phase detail

### Phase 1 — human baseline (`475.21`)

Reproduce the reference outputs using the original downloaded hg38 refdata.

- `475.16` **P0** — **fixed.** bowtie2 missing from `PATH`, plus a swallowed exit code that
  reported the failure as a successful zero-alignment run. Full-dataset rerun pending.
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
| tRNA | gtRNAdb `hg38-tRNAs.fa` (`gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/`) |
| miRNA | miRBase gff3 |
| rRNA | RefSeq GenBank `NR_046235.3`, `NR_145819.1`, `NR_146117.1`, `NR_146144.1`, `NR_146151.1` |

The tRNA and rRNA rows above are exact, verified 2026-08-11 (`475.36`, `475.37`). The tRNA portion
is the gtRNAdb FASTA rendered twice per locus: `N*10 + genomic + CCA + N*10`, and a
`_withgenomeflank` twin `N*10 + locus±50bp + N*10` with no CCA. The rRNA portion is the whole
GenBank record as `-45S` plus the `18S`/`28S` `misc_feature` spans; the annotated `5.8S` span is
not emitted. **The assembly tRNA track (`hg38.trna.tsv.gz`) is not a substitute for the gtRNAdb
FASTA** — it carries different names and different sequences, and using it is what produced the
long-standing "the reference used a different gtRNAdb release" dead end.

Mouse rRNA does **not** follow the human rule (`475.40`). Its precursor `NR_046233.2` annotates
only `source`/`gene`/`rRNA` over the full 13,400 bp — no `misc_feature` to read — so the subunit
spans are located by anchoring the standalone RefSeq records `NR_003278.3` (18S) and
`NR_003279.1` (28S) on their terminal 40-mers, supplied via `--rrna-subunit LABEL=PATH`. The
emitted bases are always the precursor's: 18S is an exact copy of Rn18s, while the precursor's
28S is 4,727 bp against Rn28s1's 4,730 — a real 3 bp indel between rDNA copies, which is why a
whole-substring match is not required.

The RepBase file has 1,356 records covering **1,224/1,224** index families (100%), verified by
header prefix *and* canonical sequence digest: 700 byte-exact, 524 matching after IUPAC→N (which
is precisely the `.fixed.fa` step).

**Rejected, do not revisit:** the RepBase 24.01 two-library method and the RepeatMasker-Edition
EMBL parse. The 97.96% coverage figure that motivated them was an artifact of missing
`pseudo.ref`/`vrtrep.ref`/`invrep.ref`, not a real library gap.

### MASTER_FILELIST structure (T-09, settled 2026-08-11)

It is a **concatenation of source lists**, not a table. Column 5 names the list each row came
from. Blocks, in file order:

| Block | hg38 n | Columns |
|---|---|---|
| Gencode + rRNA | 5,276 | ENST │ ENSG │ gene_name │ **FAMILY** │ `genelists.FAMILY` |
| RepBase families | 1,224 | NAME │ FAM │ FAM │ FAM │ **CLASS** |
| tRNA | 864 | name │ anticodon │ anticodon │ `tRNA` │ source list |
| SimpleRepeat k-mers | 501 | `AT_SimpleRepeat` │ `Simple_repeat` ×4 |
| miRNA | 3,765 | `MI0022705` │ `miRNA` ×3 │ mirbase name (+ `-proximal` twin) |
| rmsk leftovers | 14,220+ | NAME │ FAM │ FAM │ FAM │ CLASS, and `(XXX)N` |

Three rules that are easy to get wrong:

- **Never deduplicate by column 1.** 61 names appear twice with *different* annotations
  (`DNA1_MAM` is `TcMar` in the RepBase block and `TcMar-Tc1` in the rmsk block).
- **Repeat class/family are UCSC `repClass`/`repFamily`,** not the RepBase header's class field —
  `L2B_CR1_Eutheria` is family `L2`, class `LINE`. The `examples/inputs` RepeatMasker GTFs are
  stripped to `gene_id`, so this needs `refdata/<asm>.rmsk-class-family.tsv.gz`.

Repeat annotation resolves in this order (`475.42`), most authoritative first:

1. the curated small-RNA family map — `5S` is `RNA5S` everywhere, never `rRNA/rRNA`
2. the assembly rmsk table, including UCSC's slash-qualified aliases (`ALR/Alpha` → `ALR`)
3. `refdata/hg38.repbase-class-family.tsv` — the 1,224 curated pairs lifted from the reference's
   RepBase block. It covers the 354 hg38 families with no genomic instances in the modern rmsk
   track, and doubles as the **cross-species** source: pass the human file to a mouse build and
   215 mouse families inherit by exact family name. It is consulted *after* rmsk so current
   RepeatMasker calls win and post-2020 reclassification stays visible. Being reference-derived it
   is circular, so exclude it when measuring hg38 reproduction fidelity (`475.22`).
4. two systematic name aliases: `NAME_I` → `NAME_I-int` (RepeatMasker's LTR internal-segment
   convention) and the species tags `_MM`/`_HS`. Only those two tags — `_LTR`, `_DNA`, `_II`,
   `_MAM` are class tokens inside the family name and stripping them would corrupt it.
- **Column 4 is behaviorally significant.** `read_in_filelists` reads it as `$type_label`, and
  `print_output` counts with `$count{$ensttype_join}++`, so col4 *is* the family reads are
  counted under. Between-family block order and within-family row order are not: the priority
  comparison is keyed by family (`flags{$ensttype}`), so `priority_n` only chooses which element
  is named as a family's representative in the SAM output.

**Gencode family assignment**, measured on hg38: rmsk small-RNA overlap ≥50% resolves 3,932 of
5,261 rows and `repName → family` is effectively a function; gene_name resolves another 831; 473
resolve to neither and need `--family-override` (`475.41`). The hg38 originals came from a
per-family `genelists.*` directory that no longer exists — see the commented `read_in_filelists`
call at `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl:43`.

Build order is **filelist → index**: the index selects its Gencode portion by filelist
membership, so `generate_master_filelist.py` takes `--repbase-species-fasta` directly rather than
the index's provenance sidecar.

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
| UniqueGenomicElements | **Regenerated 2026-08-11 (`475.12`)** | was 275/277 MB against hg38's 208 MB; now 189 MB (mm10, 5,154,546 rows) and 195 MB (mm39, 5,327,711). See below |
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

1. ~~**mm39 is missing two inputs entirely**~~ — **closed 2026-08-11 (`475.6`)**. Neither the
   index nor the filelist reads the tRNA track any more; both take the gtRNAdb FASTA
   (`mm39-tRNAs.fa`, 407 loci). The GRCm39 miRBase gff3 is `https://mirbase.org/download/mmu.gff3`
   — miRBase **v23**, GRCm39, 1,190 primary transcripts — now at
   `examples/inputs/mm39/downloaded/mmu.gff3`, giving the mm39 filelist 2,380 miRNA rows.
   mm10 stays on its v22/GRCm38 copy, because v23 ships mouse coordinates on GRCm39 only.
   **The v23 file mixes seqid conventions** — 1,164 Ensembl-style rows (`1`, `X`) and 26 already
   `chr`-prefixed — so `parse_gff3_mirna` normalizes through `ucsc_chrom`; everything else in the
   repo is UCSC-named and un-normalized seqids match no chromosome. ~~T-07 still reads the tRNA
   track.~~ **T-07 now takes `--gtrnadb-fasta` too** (`-dvy`) — see below.
2. ~~**No mouse chromosome allowlists**~~ — **closed 2026-08-11 (`475.7`)**. Both are now pinned
   to the UCSC base release (`{mm10,mm39}.chrom.sizes`), 66 and 61 sequences. The earlier mm39
   placeholder was derived from `GRCm39.primary_assembly.genome.fa.fai`, which names unplaced
   scaffolds Ensembl-style (`GL456210.1`) where every track the generators read is UCSC-named
   (`chr1_GL456210v1_random`); it silently dropped 8,292 RepeatMasker rows on 39 scaffolds.
   Applied to RepeatMasker the rule drops **186,003 rows on 173 patch scaffolds for mm10** and
   **0 for mm39** (UCSC's mm39 rmsk snapshot carries nothing that postdates the release — its
   61 scaffolds are exactly `mm39.chrom.sizes`).
3. ~~**No mouse acceptance criteria existed.**~~ **Closed 2026-08-12 (`475.13`).**
   `bin/python/refdata_generation/validate_refdata_set.py` runs six cross-artifact consistency
   checks over an assembly's filelist + index + BED and exits non-zero on failure. hg38
   reference PASSES (5 checks, 1 skipped — no provenance sidecar for a downloaded 2020 index),
   mm10 and mm39 PASS 6/6. Its own acceptance requirement — that the checks *fail* on known-bad
   input — is met against all three real defect classes: the T-05 `20260514` indices fail A+C,
   the T-07 `.stale-20260514` BEDs fail B+E, the M-5 `.pre-475.11` filelists fail D. Full
   results and the two checks whose wording had to change: changelog §8.

### M-6 — mouse UniqueGenomicElements, regenerated 2026-08-11 (`475.12`)

Two generator defects had to be fixed first; both were found by measuring, not by inspection,
and both are written up in `docs/CHANGELOG-refdata-validation-2026-07-25.md` §T-07:

- `-dvy` — the tRNA source. The UCSC `{assembly}_tRNAs` track gives 222 spurious / 23 missing
  rows on hg38; the gtRNAdb FASTA gives 432/432 exactly. mm39 has no UCSC track at all.
- `-72c` — `read_gencode_allowed_transcripts` matched `ENST` only, returning **0** allowed
  transcripts for mouse and dropping the whole Gencode contribution.

With both fixed, hg38 re-measures at recall 0.999937 / precision 0.999957 (was .99993/.99991),
so the mouse run is not resting on a rule that human never passed.

| | hg38 (reference) | mm10 | mm39 |
|---|---|---|---|
| rows | 5,618,483 | 5,154,546 | 5,327,711 |
| size | 208 MB | 189 MB | 195 MB |
| RepeatMasker | 5,607,738 | 5,147,736 | 5,320,771 |
| Gencode | 4,670 | 2,721 | 2,963 |
| tRNA | 432 | 408 | 407 |
| miRNA (+proximal) | 1,881 (+3,762) | 1,227 (+2,454) | 1,190 (+2,380) |
| simple-repeat (trf) rows | 0 | 0 | 0 |
| distinct col4 names | 24,436 | 7,137 | 24,234 |
| **names unresolved in MASTER_FILELIST** | **0** | **0** | **0** |

The last row is the acceptance criterion that needs no gold standard: `read_peakfi` uppercases
col4 and looks it up in `convert_enst2type`, which `read_in_filelists` fills from column 1 of the
MASTER_FILELIST — an unresolved name is an untyped peak. The superseded 2026-05-14 mouse BEDs
score **1,826,959 unresolved rows (mm10)** and **1,916,058 (mm39)** on the same test, almost all
of them the trf simple-repeat track they should never have carried (1,687,263 / 1,641,063 rows)
plus the unrestricted Gencode dump (142,351 / 278,326). They are kept as
`UniqueGenomicElements.{mm10,mm39}.bed.stale-20260514` until `M-8` passes.

mm39's stale BED also had **zero** tRNA and zero miRNA rows — neither input existed before
`475.6`.

Reproduce (`$A` = mm10 with `gencode.vM23`, mm39 with `gencode.vM38`; ~40 s each):

```bash
PYTHONPATH=bin/python/refdata_generation \
/tscc/nfs/home/bay001/miniconda3/envs/marine_environment/bin/python \
  bin/python/refdata_generation/generate_unique_genomic_elements.py \
  --repeatmasker    examples/inputs/$A/downloaded/$A.repeatmasker.tsv.gz \
  --gtrnadb-fasta   examples/inputs/$A/downloaded/$A-tRNAs.fa \
  --gff3            examples/inputs/$A/downloaded/mmu.gff3 \
  --parsed-ucsc     examples/inputs/$A/gencode.vM##.annotation.gtf.parsed_ucsc_tableformat \
  --master-filelist examples/inputs/$A/MASTER_FILELIST.20260811.*.list \
  --chrom-allowlist refdata/$A.chrom-allowlist.txt \
  --assembly $A --output examples/inputs/$A/UniqueGenomicElements.$A.bed
```

miRBase pairing is not interchangeable: mm10 takes the v22/GRCm38 `mmu.gff3`, mm39 the
v23/GRCm39 one.

### M-5 — mouse DROP_EXACT, re-derived 2026-08-11 (`475.11`)

The `_u1_` in `mus_musculus_repbase_u1_fixed_v2.fastq` is a block of 11 bare-header records
appended to the mouse library — `U1 U2 U3 U4B U5B1 U6 UHG U13 U7 U14 U8` — standing in for the
human library's `U<n>_snRNA_Homo_sapiens` records. Until now all but three were kept, and the
MASTER_FILELIST showed the double-count directly: its RepBase block carried

```
U1   RNU1 RNU1 RNU1 RNU1        <- family RNU1, which genelists.RNU1 already supplies (204 rows)
U2   RNU2 ...                   <- RNU2 (48)
U6   RNU6 ...                   <- RNU6 (936) + RNU6ATAC (22)
U7   RNU7 ...                   <- RNU7 (15)
U4B  U4B  ...                   <- unresolved, so its own family; Gencode supplies RNU4 (38)
NR   NR   ...                   <- NR_046235.1, the HUMAN 45S RefSeq, resolving to family "NR"
```

Those six are now dropped, plus `U7_snRNA_Vertebrata` (the well-formed twin — the human list
drops that form, so both go) and, after `-4ee`, `U5B1`.

`U5B1` was **kept** on 2026-08-11 and **dropped on 2026-08-12**, once `-4ee` was fixed. It was
kept because the mouse MASTER_FILELIST then had no RNU5 family at all — which turned out to be a
symptom of `-4ee`, not a fact about mouse: `family_from_gene_name` matched human symbols
literally, so vM38's `Rnu5g` never resolved. With the fix Gencode supplies **RNU5G** and the
ordinary rule applies, exactly as the human list drops U5B1 against RNU5A–RNU5F. Mouse's U5
annotation is thin — one transcript against human's 34, and the Rfam tier supplies no U5 family
either — but that is a statement about the annotation, not a reason to keep a second competing
representation of U5 in the repeat portion. Mouse RepBase selection: **1,078 → 1,072 → 1,071**.

Measured overlap between the RepBase block's families and the Gencode block's:

| | before | after | hg38 reference |
|---|---|---|---|
| mm10 | RNU1 RNU2 RNU6 RNU7 SNORD | **SNORD** | — |
| mm39 | RNU1 RNU2 RNU6 RNU7 SNORD | **SNORD** | — |
| hg38 | | | SNORD, YRNA |

(unchanged by the `-4ee` rebuild: still `{SNORD}` on both, with RNU5G now in the Gencode block
and U5B1 gone from the repeat block, so U5 is represented once.)

The residual `SNORD` is inherited from the human rule, not a mouse exception: the kept
U3/U8/U13/U14 snoRNA records carry family SNORD, and the hg38 reference does the same (and also
carries YRNA). Mouse now overlaps strictly less than the validated human artifact.

Note that `U1`, `U2`, `U6` and `U7` still appear *later* in the file, in the rmsk-leftovers
block with families RNU1/RNU2/RNU6/RNU7 — because RepeatMasker's `repName` for those loci is
literally `U1`, `U2`, … The hg38 reference has exactly the same rows (lines 11786–12345), so
this is the intended shape; the rule is about the RepBase block only.

All four mouse artifacts were rebuilt, twice — once for M-5 and again after `-4ee`, which
changed the Gencode block and so propagated into the index and the BED as well:

| | filelist rows | index records | BED rows | Gencode in BED |
|---|---|---|---|---|
| mm10 before | 11,526 | 5,488 | 5,154,546 | 2,721 |
| mm10 after M-5 | 11,520 | 5,482 | — | — |
| **mm10 after `-4ee`** | **11,568** | **5,530** | **5,154,593** | **2,768** |
| mm39 before | 25,871 | 5,727 | 5,327,711 | 2,963 |
| mm39 after M-5 | 25,869 | 5,721 | — | — |
| **mm39 after `-4ee`** | **25,894** | **5,746** | **5,327,735** | **2,987** |

Contracts re-checked after the final rebuild: every index header resolves in its filelist
(5,530/5,530 and 5,746/5,746), the M-6 BED check returns **0 unresolved rows** on both, and the
BEDs still carry zero trf simple-repeat rows.

The filelist invocation was first validated by reproducing the committed 2026-08-11 mm39 file
**byte-identically** before changing anything. It needs both curated tables, human first:
`--repbase-class-family refdata/hg38.repbase-class-family.tsv --repbase-class-family
refdata/mm.repbase-class-family.tsv`.

## Standing risks

- **Line-count acceptance criteria mask content bugs.** T-09 passed a ±1% row-count check
  (26,252 vs 26,422 = 0.9936) while col4 was entirely the wrong field. Any new criterion must
  test content, not shape.
- ~~**`M-5` is a correctness dependency, not bookkeeping.**~~ **Confirmed and closed 2026-08-11.**
  The prediction that the human list transfers was wrong in both directions: four families
  (RNU1, RNU2, RNU6, RNU7, plus U4B and the human `NR_046235.1` record) were being double-counted
  across the repeat and Gencode portions, and one (`U5B1`) must be kept where human drops it.
  Realized 1,072 families, not the predicted 1,079. See the M-5 section above.
- **Genomic-instance fallback for unresolved families is prohibited** (FIX-PLAN §4a). It
  reintroduces the T-05 error, and mouse has no reference to catch it.
- **Fitting to hg38 is the systemic hazard.** Rules derived by diffing against the hg38 reference
  carry no independent evidence for mouse. Phase 4 exists to make that explicit, artifact by
  artifact.

## Pointers

- `.forge/stages/2-architect/notes/FIX-PLAN-bowtie2-index.md` — T-05 rewrite plan (§2–§4b, §7)
- `docs/MASTER_FILELIST-decisions.md` — **the full decision record for the filelist generator**:
  every rule, the function that implements it, the measurement that justified it, and the
  alternatives that were rejected with the numbers that killed them
- `.forge/stages/2-architect/notes/T-07-per-source-selection.log` — the five UniqueGenomicElements rules
- `.forge/debug/causal-chain.json` — prior (now superseded) 0-repeat-reads diagnosis
- `docs/CHANGELOG-refdata-validation-2026-07-25.md` — hg38 reproduction validation results
- `bd memories` — `repbase-build-method-confirmed`, `repbase-consensus-source`,
  `refdata-validation-policy`, `refdata-review-corrections`, `refdata-hg38-reproduction-results`
