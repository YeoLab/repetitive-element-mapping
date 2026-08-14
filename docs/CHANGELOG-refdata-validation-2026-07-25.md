# Changelog — mm10/mm39 reference-data validation (hg38 reproduction)

- **Date:** 2026-07-25
- **Issue:** `repetitive-element-mapping-475` (Reference data generation, mm10/mm39)
- **Scope:** Ran the hg38 reproduction-test-first validation gate (tasks T-03/T-05/T-07/T-09)
  for the four `bin/python/refdata_generation/` generators, diagnosed failures, and located
  + confirmed the original repeat-consensus source.
- **Outcome:** 1 of 4 generators passes; 3 fail on content. Root cause of the bowtie2-index
  failure identified and the correct source method proven by byte-level sequence matching.
  No generator code changed yet — this is the validation + diagnosis record and fix plan.

---

## 1. Validation policy

> **Content-equivalence = PASS — for `parsed_ucsc_tableformat` ONLY.** A sorted-content MD5
> match satisfies the "faithful reproduction" AC for that file, because its consumer reads
> it order-independently. **This policy does NOT generalize to the other three files.**

Verified order-independence of the parsed-table consumer (rows keyed by transcript id):

```perl
# bin/perl/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl
sub read_gencode {
    open(F,$fi) || die;
    while (<F>) {                      # line order irrelevant:
        my @tmp = split(/\t/,$_);
        my $enst = $tmp[1];            # every row keyed by transcript id
        ...
    }
}
```

**MASTER_FILELIST order IS behaviorally significant** and must be reproduced (or validated
behaviorally), NOT treated as order-independent. `read_in_filelists` assigns a per-line
priority in file order, and that priority breaks ties between equal-score mappings:

```perl
# bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl:463
sub read_in_filelists {
    my $priority_n = 0;
    for my $line (<F>) {                                   # file order == priority order
        my ($allenst,$allensg,$gid,$type_label,$typefile) = split(/\t/,$line);
        ...
        $convert_enst2type{$enst} = $type_label.":".$priority_n;   # priority encoded per enst
        $priority_n++;
    }
}
```

Order-dependence of the remaining two files (`bowtie2_index` FASTA, `UniqueGenomicElements`
BED) has not been assumed either way; validate content, not sorted content, for those.

## 2. Results summary

| Task | Generator | Count/line check | Content check | Verdict |
|---|---|---|---|---|
| T-03 | `generate_parsed_ucsc_tableformat.py` | 249,044 = 249,044 | sorted-MD5 identical | **PASS** |
| T-05 | `generate_bowtie2_index.py` | 26,353 vs 7,606 (0.29) | wrong headers + wrong source | **FAIL** |
| T-07 | `generate_unique_genomic_elements.py` | 6,988,833 vs 5,618,483 (0.80) | multi-source over-selection; schema OK | **FAIL → FIXED 2026-07-26** (recall .99993 / precision .99991) |
| T-09 | `generate_master_filelist.py` | 26,252 vs 26,422 (0.99) | biotype vs repeat-family labels | **FAIL (content)** |

Per-task logs: `.forge/stages/2-architect/notes/T-0{3,5,7,9}-*.log`.

### T-03 — parsed_ucsc_tableformat: PASS

```bash
python bin/python/refdata_generation/generate_parsed_ucsc_tableformat.py \
  --gtf examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf \
  --output /tmp/test_hg38_parsed.tsv

# Raw diff is huge (row order differs) but content is identical:
sort /tmp/test_hg38_parsed.tsv | md5sum      # 01244eef2f97da25128feff32e79fc06
sort examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat | md5sum
                                             # 01244eef2f97da25128feff32e79fc06  -> IDENTICAL
```

The only difference is that the script emits `sorted((gene_id, transcript_id))` while the
reference preserves the UCSC tool's internal order. Content-equivalent → PASS.

### T-05 — bowtie2 index: FAIL (two bugs)

```bash
# Generated headers carry a pybedtools getfasta coordinate suffix:
>ENST00000516207.1::chrX:88148791-88148918(+)     # generated (WRONG)
>ENST00000384010.1                                # reference   (clean id)

grep -c '^>' <generated>.fa   # 26,353
grep -c '^>' <reference>.fa   #  7,606   -> 3.5x over-generation
```

Cause: the script does `pybedtools.getfasta` on **every genomic repeat instance** instead
of using **one consensus sequence per repeat family**. Wrong unit (instance vs family) and
wrong header format. See §3–§4.

### T-07 — UniqueGenomicElements: FAIL (over-selection, NOT schema)

```bash
wc -l <generated>.bed   # 6,988,833
wc -l <reference>.bed   # 5,618,483   ratio 0.80
```

**Correction (2026-07-25 review):** an earlier draft claimed the BED *schema* was wrong
(that Gencode rows should carry `-` as the name). That was an analysis error — it compared a
generated Gencode row against an unrelated RepeatMasker reference row. The reference BED
**does** contain Gencode rows in the generator's exact format:

```bash
awk -F'\t' '$4 ~ /^ENST/' <reference>.bed | wc -l    # 4,670 rows
# e.g.  chr6  159785593 159785733  ENST00000384183.1  -  -   (col4=ENST, col5='-', col6=strand)
```

matching `parse_parsed_ucsc` in `generate_unique_genomic_elements.py`. The **schema is
correct.** The real
defect is **selection** — and it is **multi-source**, not Gencode-only.

**Update (2026-07-26, issue -475.2 — RESOLVED).** A Gencode-only fix cannot close the gap:
restricting 249,043 transcripts to 4,670 removes 244,373 of the 1,370,350 excess rows.
A per-source diff attributes the delta completely:

| source | over-generated | rule |
|---|---|---|
| simple repeats | 1,049,715 | **exclude entirely** — reference has 0 such rows |
| Gencode | 244,608 | restrict to MASTER_FILELIST minus rRNA/mito families |
| RepeatMasker | 75,952 | chromosome allowlist pinned to the assembly release |
| tRNA | 631 | col5 = family, not GTF score (a derivation bug) |
| miRNA | 231 | miRBase version drift |

With all five rules implemented the generator reproduces hg38 at **recall 0.99993 /
precision 0.99991** (5,618,567 vs 5,618,483 rows); the 5.6M RepeatMasker rows — 99.8% of
the file — match exactly. Residuals are source-version drift (miRBase, gtRNAdb, scaffold
naming), not selection error. Full derivation, rejected alternatives, and the mm10/mm39
application: `.forge/stages/2-architect/notes/T-07-per-source-selection.log`.
Tests: `tests/refdata_generation/test_unique_genomic_elements.py` (19, passing).

**Update (2026-08-11, issues `-dvy` and `-72c`).** Two of those residuals were not drift.

*tRNA was the wrong source, not a version skew.* Rule 4 fixed col5 but kept reading the
UCSC `{assembly}_tRNAs` GTF track, which on hg38 yields 631 rows against the reference's
432: **409 shared, 222 spurious, 23 missing**. The 222 are named in a convention no
MASTER_FILELIST carries (`nm-tRNA-Tyr-GTA-chr1-142`, `tRNA-Und-NNN-chr1-1`). Reading the
same rows from the **gtRNAdb genomic FASTA** — coordinates taken from the header's
trailing `chr6:28795964-28796035 (-)` field, the source the index and filelist generators
already use — reproduces the reference tRNA rows **exactly: 432/432, zero spurious, zero
missing**. `--gtrnadb-fasta` is now the preferred source and takes precedence over
`--trna`. This was also a hard mouse blocker: mm39 has no UCSC tRNA track at all, and
mm10's names its rows `chr1.tRNA1555-GluTTC`.

*The Gencode allowlist was human-only.* `read_gencode_allowed_transcripts` filtered on
`col1.startswith('ENST')`. Mouse ids are `ENSMUST`, so it returned **0** allowed
transcripts for both mouse MASTER_FILELISTs and silently dropped the entire Gencode
contribution. Now matched with `^ENS[A-Z]*T[0-9]`; human selection is unchanged (4,679).

Re-measured hg38 with both fixes: **recall 0.999937 / precision 0.999957** (5,618,368 vs
5,618,483 rows) — spurious rows 507 → 240, and the tRNA residual is now zero. What is
left is genuine version drift: miRNA 231 spurious / 120 missing (miRBase) and Gencode
9 / 235 (v33 coordinates).

A cross-artifact check that needs no reference, and the one mouse is held to: **every
col4 name in the BED, uppercased, must resolve in column 1 of the MASTER_FILELIST** —
that is what `read_peakfi` / `read_in_filelists` in
`duplicate_removal_inline_paired...pl` require to type a peak at all. The hg38 reference
scores 24,436/24,436 distinct names, 0 unresolved.

### T-09 — MASTER_FILELIST: FAIL (content)

```bash
wc -l <generated>.tsv   # 26,252   (ratio 0.9936 — line-count AC PASSES)

# but col4 (family) is wrong:
cut -f4 <generated>.tsv | sort | uniq -c | sort -rn | head   # misc_RNA, snRNA, snoRNA ...  (Gencode BIOTYPES)
cut -f4 <reference>.tsv | sort | uniq -c | sort -rn | head   # RNU6, YRNA, RN7SL, SNORD ... (REPEAT FAMILIES)
```

Transcript-id overlap is only ~42% even after stripping version suffixes. The ±1% line-count
AC passed while the content is substantially wrong — a weak proxy.

## 3. Root cause of the bowtie2-index failure

The reference index FASTA (7,606 seqs) decomposes as:

| Component | Count | Correct source |
|---|---|---|
| Gencode transcripts `>ENST…` | 5,002 | genome `getfasta` (OK; must strip `::coords`) |
| tRNA | 864 | tRNA fasta |
| rRNA `>NR_…-18S/28S/45S` | 15 | RefSeq |
| Repeat-family **consensus** `>ALUY`,`>L2C` | 1,224 | **RepBase consensus (NOT genomic)** |
| SimpleRepeat kmers `>AAAAAC_SimpleRepeat` | 501 | curated kmer set |

Only **868 of the 15,600** families in `repeatmasker.tsv.gz` appear in the index → it is a
curated, consensus-driven set, not a mechanical genomic extraction.

## 4. Repeat-consensus source — LOCATED and CONFIRMED

Original author trail (from perl comments): `elvannostrand` /
`.../RNA_type_analysis/`. The consensus library is present and readable:

```
/tscc/projects/ps-yeolab4/genomes/RepBase24.01.fasta/
├── RepBase24.01.fasta/*.ref                       # RepBase 24.01, taxon-split FASTA
│      humrep humsub prirep prisub mamrep mamsub    #   human set
│      rodrep rodsub mousub ratsub                  #   mouse set
│      simple.ref
└── RepBaseRepeatMaskerEdition-20181026.tar.gz
       └── Libraries/RMRBSeqs.embl                  # RepeatMasker Edition, 49,011 EMBL seqs
```

RepBase `.ref` header format: `>NAME<TAB>CLASS<TAB>species` (e.g. `>ALU  SINE1/7SL  Primates`).

### Confirmed build method (byte-level sequence match)

Sequences were extracted by family name from each source and compared to the reference index:

```python
# getseq.py — extract uppercased concatenated sequence by header first-token (FASTA)
import sys
fa, name = sys.argv[1], sys.argv[2].upper()
seq, grab = [], False
for line in open(fa):
    if line.startswith('>'):
        grab = (line[1:].split('\t')[0].split()[0].strip().upper() == name)
    elif grab:
        seq.append(line.strip().upper())
print(''.join(seq))
```

```python
# getembl.py — extract sequence by ID from EMBL (RMRBSeqs.embl)
import sys
embl, name = sys.argv[1], sys.argv[2].upper()
grab = inseq = False; seq = []
for line in open(embl):
    if line.startswith('ID'):
        grab = (line.split()[1].rstrip(';').upper() == name); inseq = False
        if grab: seq = []
    elif grab and line.startswith('SQ'):
        inseq = True
    elif grab and line.startswith('//'):
        if seq:
            break
    elif grab and inseq:
        seq.append(''.join(c for c in line if c.isalpha()))
print(''.join(seq).upper())
```

| Family | RepBase 24.01 `.ref` | RM Edition EMBL | Index matches |
|---|---|---|---|
| `LTR18B`, `MSTC` | verbatim | verbatim | both |
| `MER5A` | 2 diffs (IUPAC R/Y → N) | 5 real diffs | **`.ref`** (primary) |
| `EULOR1` | absent | verbatim | **RM Edition** (fallback) |

**Method that mirrors the original pipeline:**
1. **Primary:** RepBase 24.01 `.ref` consensus (taxon-split → species selection).
2. **Fallback:** RepeatMasker Edition `RMRBSeqs.embl` for families absent from 24.01.
3. **Normalize** IUPAC ambiguity codes → `N` (the `.fixed.fa` step).

Coverage of the 1,224 index families: RepBase 24.01 = 1,116; RM Edition = 864;
**union = 1,199/1,224 = 97.96%** — below the ≥99% gate.

**Superseded (2026-07-27).** The two-library reconstruction above is no longer the plan.
The 25 "missing" families were never absent from RepBase — only from the *human ref-set
file selection* (18 are in `pseudo.ref`, 6 in `vrtrep.ref`, 1 in `invrep.ref`; adding those
three files takes RepBase 24.01 alone to 1,221/1,224 and the union to 1,224/1,224).

The actual source is simpler and exact:

```
/tscc/projects/ps-yeolab4/genomes/RepBase18.05.fasta/species_specific/
    homo_sapiens_repbase_fixed_v2.fasta       1,356 records
    mus_musculus_repbase_u1_fixed_v2.fastq    1,194 records  (FASTA content)
```

One file covers **1,224/1,224 = 100%** of the hg38 index repeat families, verified by both
`<FAMILY>_` header prefix and canonical sequence digest (700 byte-exact; the other 524 match
after IUPAC→N — precisely the `.fixed.fa` step). The `_fixed_v2` naming and Dec-2017 date fit
the 2020-12-03 index build, and a mouse counterpart exists.

This removes the second library, the EMBL parser, the alias-resolution pass, and the ≥99%
coverage gate — the repeat portion is an exact reproduction, not a thresholded one. Full
derivation, including the family-selection rule (1,356 → 1,224, exact): FIX-PLAN §2–§4b.
Substituting a genomic instance for an unresolved family remains prohibited (FIX-PLAN §4a).

### Reproduction of the analysis

```bash
DIR=/tscc/projects/ps-yeolab4/genomes/RepBase24.01.fasta/RepBase24.01.fasta
# index repeat-family headers (uppercased, no ENST/tRNA/NR/SimpleRepeat)
grep '^>' <reference>.fa | grep -v '^>ENST' | grep -vi 'tRNA' | grep -v '^>NR_' \
  | grep -v '_SimpleRepeat' | sed 's/^>//' | tr a-z A-Z | sort -u > idx_fams.txt   # 1,224
# RepBase 24.01 names (human set)
cat $DIR/{humrep,humsub,prirep,prisub,mamrep,mamsub,simple}.ref \
  | grep '^>' | cut -f1 | sed 's/^>//' | tr a-z A-Z | sort -u > repbase_hg38.txt
# RM Edition names
tar -xzf $DIR/../RepBaseRepeatMaskerEdition-20181026.tar.gz -C /tmp Libraries/RMRBSeqs.embl
grep '^ID' /tmp/Libraries/RMRBSeqs.embl | awk '{print $2}' | sed 's/;//' | tr a-z A-Z | sort -u > rmedition_names.txt
# union coverage
comm -12 idx_fams.txt <(sort -u repbase_hg38.txt rmedition_names.txt) | wc -l   # 1,199 / 1,224
```

## 5. Fix plan (next phase)

Full plan: `.forge/stages/2-architect/notes/FIX-PLAN-bowtie2-index.md`. Summary:

- Rewrite the repeat portion of `generate_bowtie2_index.py` to read consensus from RepBase
  24.01 `.ref` (primary) → RM Edition `RMRBSeqs.embl` (fallback), IUPAC→N normalized;
  add `--species {human,mouse}` to pick the ref set. Strip the `::coords` suffix from
  Gencode headers.
- Family selection: species-primary refs ~wholesale ∪ genome-present families with a
  consensus. Current union coverage is 97.96% — an **alias-resolution pass is required** to
  reach the ≥99% gate before implementation is accepted.
- **Validate at the sequence level**, not just headers (per-header normalized sequence hash
  + no-duplicate-header check); headers alone can match while sequences differ.
- **T-09 (MASTER_FILELIST)** needs more than the repeat family set: explicit Gencode
  column derivation (`col3=gene_name`, `col4=family` by stripping the copy-number suffix
  e.g. `RNU6-2→RNU6`, `col5=genelists.{family}`), a repeat family→class map
  (`ALUY→Alu→SINE`), **and row-order preservation** (order sets mapping priority — see §1).
- **T-07 (UniqueGenomicElements)**: ~~pending~~ **DONE 2026-07-26.** Schema was correct; the
  fix was selection, across five sources (not Gencode alone) — see §T-07 above.
  ~~mm10/mm39 still blocked on T-09~~ **mm10/mm39 regenerated 2026-08-11** (`475.12`), after
  two further fixes to the tRNA source (`-dvy`) and the human-only Gencode id match (`-72c`);
  hg38 re-measures at recall 0.999937 / precision 0.999957.

## 6. Environment notes (tools needed to run the generators)

The `ecliprepmap/1.0.0` runtime env lacks the refdata-generation dependencies. Use:

| Tool | Source env |
|---|---|
| `pybedtools` 0.9.0, `bedtools` | `~/miniconda3/envs/snakemake738` |
| `samtools` | `.../miniconda_tscc2/envs/ecliprepmap-0.1.0` |
| `bowtie2-build` | `.../miniconda_tscc2/envs/bowtie2-2.5.4` |

```bash
export PATH="$HOME/miniconda3/envs/snakemake738/bin:\
/tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/bowtie2-2.5.4/bin:\
/tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/ecliprepmap-0.1.0/bin:$PATH"
```

## 7. Status of generated mouse files

Because mm10/mm39 were produced by the three failing generators, their bowtie2 index,
UniqueGenomicElements, and MASTER_FILELIST are wrong in the same ways; only the
`parsed_ucsc_tableformat` files are trustworthy. The mm10 "32% missing IDs" report is a
downstream symptom of the T-05 genomic-getfasta approach. Reference generation is **not
~95% complete** — the earlier "~95%" figure counted scripts written, not outputs validated.

**Status as of 2026-07-26:**

| Generator | hg38 reproduction | mouse outputs |
|---|---|---|
| `generate_parsed_ucsc_tableformat.py` | PASS | valid |
| `generate_unique_genomic_elements.py` | **PASS** (recall .999937 / precision .999957 after `-dvy`, `-72c`, 2026-08-11) | **regenerated 2026-08-11** (`475.12`); 0 names unresolved against the mouse MASTER_FILELIST |
| `generate_bowtie2_index.py` | ~~FAIL~~ **rewritten** (T-05, `-7ee`) | **regenerated 2026-08-12** |
| `generate_master_filelist.py` | ~~FAIL~~ **rewritten** (T-09, `-gz2`) | **regenerated 2026-08-12** |

~~Remaining order of work: T-05 → T-09 → regenerate mm10/mm39 UniqueGenomicElements →
integration test.~~ All four mouse artifacts are now generated and pass the §8 acceptance
suite. Remaining: `M-7` (`475.8`, re-run parsed_ucsc + checksum), `M-9` (`475.14`, end-to-end
smoke run on mouse eCLIP data), `M-10` (`475.15`, package + document), and `P-5` (`475.22`),
which is the last checkpoint against human ground truth and still gates M-4.

---

## 8. The mouse acceptance suite (2026-08-12, `475.13` / M-8)

Every criterion in sections 1–7 is "reproduce the reference file". Mouse has no reference, so
that method does not transfer, and until now the mouse artifacts had no acceptance gate at all —
the single largest risk in the effort. `bin/python/refdata_generation/validate_refdata_set.py`
replaces reference-diff with six internal-consistency checks across the three artifacts of one
assembly, and exits non-zero on any failure.

| | check | catches |
|---|---|---|
| A | every index FASTA header resolves in MASTER_FILELIST column 1 | an indexed sequence `read_in_filelists` cannot type |
| B | every UniqueGenomicElements column-4 name resolves there too | an untyped peak in `read_peakfi` |
| C | no index header carries a `::` coordinate suffix | T-05, the per-genomic-instance index |
| D | the RepBase-block/Gencode-block family overlap is within the hg38 reference's own | M-5, a family counted in both portions |
| E | BED source profile: 0 trf rows, rmsk ≥ 99%, miRNA rows == 3 × entries, no source empty | T-07, and any optional input silently omitted |
| F | the provenance sidecar covers every indexed RepBase family | T-05 from the other side |

**Results.**

| set | outcome |
|---|---|
| hg38 reference | **PASS** (5 passed, 1 skipped — no provenance sidecar exists for a downloaded 2020 index) |
| mm10 generated | **PASS** (6/6) |
| mm39 generated | **PASS** (6/6) |

**The checks have teeth**, demonstrated against the three known-bad artifact classes this
project actually produced rather than against synthetic damage — all six runs exit 1:

| known-bad input | fails |
|---|---|
| mm10/mm39 index `20260514` (T-05) | A (0/5,875 and 0/22,357 headers resolve), C (5,874 and 22,356 headers carry `::`) |
| mm10/mm39 BED `.stale-20260514` (T-07) | B (1.83M and 1.92M unresolved rows), E (rmsk 0.744 and 0.735; 1,687,263 and 1,641,063 trf rows) |
| mm10/mm39 filelist `.pre-475.11` (M-5) | D (overlap `RNU1 RNU2 RNU6 RNU7 SNORD`), plus A and B as collateral |

Unit tests in `tests/refdata_generation/test_validate_refdata_set.py` (13) pin each check's
failure case independently, so a check cannot silently stop being able to fail.

**Two of the six could not be stated as the issue worded them, and the wording was wrong rather
than the artifacts:**

- *"every MASTER_FILELIST id resolves to an index header **and vice versa**."* Only
  index→filelist is an invariant. 18,748 of hg38's 26,354 filelist ids (71%) have no index
  sequence **by design** — the rmsk-leftover and SimpleRepeat blocks contribute names the dedup
  perl needs in order to type a peak, not sequences to align against. The reverse direction
  fails on the reference set.
- *"no family appears in both the repeat and Gencode portions."* The hg38 reference has two,
  SNORD and YRNA, from the RepBase U3/U8/U13/U14 snoRNA records it deliberately keeps. The check
  is therefore that an assembly's overlap is a **subset** of the reference's: mouse may not
  invent an overlap human does not have. Mouse sits at `{SNORD}`.

One incidental finding, upstream data rather than a defect of ours: mm10's `MamSINE1` row
carries family `tRNA` where mm39's carries `tRNA-RTE`, each taken from its own assembly's rmsk
table. Delimiting the RepBase block on column 4 therefore cut mm10's block at 387 rows instead
of 1,071; the suite delimits on the tRNA block's column-5 source-list name instead.

---

## 9. parsed_ucsc_tableformat provenance (2026-08-12, `475.8` / M-7)

The mouse `.parsed_ucsc_tableformat` files were dated 2026-05-14, a week before
`generate_parsed_ucsc_tableformat.py` was committed, so nothing recorded that they came from the
committed generator. Re-ran it against the same GTFs:

| assembly | source GTF | rows | byte MD5 | sorted MD5 | vs on-disk |
|---|---|---|---|---|---|
| mm10 | `gencode.vM23.annotation.gtf.gz` | 142,351 | `f7b9c484a64a305098f892c4f39430f2` | same | **byte-identical** |
| mm39 | `gencode.vM38.annotation.gtf.gz` | 278,326 | `8baf63a7b437579bd1e30f68a7fdfb77` | same | **byte-identical** |
| hg38 | `gencode.v33.chr_patch_hapl_scaff.annotation.gtf` | 249,043 | `cc390c9a8de9fc00e0c49334c05d798c` | `01244eef2f97da25128feff32e79fc06` | sorted-identical |

Both mouse files reproduce byte for byte, so the 2026-05-14 artifacts *were* produced by the
committed script and no downstream artifact needed rebuilding. Byte MD5 equals sorted MD5 for
them because the generator emits `sorted(transcripts.keys())`; hg38's on-disk file is the 2020
UCSC-tool output, which preserves that tool's own row order — hence the byte difference and the
sorted-content criterion of §T-03. Its sorted MD5 is unchanged from the value recorded there on
2026-07-25, which independently confirms the generator has not drifted.

**The mm10 "duplicate" was not a duplicate.** `gencode.vM23.*` and `gencode.VM23.*` differ only
in one letter but are different kinds of file:

- `gencode.vM23.annotation.gtf.gz` (28 MB) is the real GENCODE release — `##description:
  evidence-based annotation of the mouse genome (GRCm38), version M23`, with `gene`/`transcript`/
  `exon`/`CDS` features and full attributes.
- `gencode.VM23.annotation.gtf.gz` (11 MB) is a **UCSC table-browser export** of `mm10_knownGene`
  — exon rows only, `gene_id` identical to `transcript_id`, no `gene_name` or `gene_type`. The
  same shape as `mm10.repeatmasker.tsv.gz` and `mm10.trna.tsv.gz`.

Parsed through the generator, that difference is visible in column 1: the `VM23`-derived file
carries `ENSMUST00000000001.4` where every other assembly carries a gene id
(`ENSMUSG00000000001.4`, `ENSG00000164054.15`). It is the wrong file to keep, and nothing in the
repo referenced it. Canonical name is the lowercase `gencode.vM23.annotation.gtf.parsed_ucsc_tableformat`;
the other is set aside as `.superseded-not-gencode`. mm39 never had the equivalent stray, though
the `gencode.VM38.annotation.gtf.gz` UCSC export is on disk in the same way.

---

## 10. P-5 — the pipeline on fully regenerated hg38 reference data (2026-08-12, `475.22`)

Phase 5 of `docs/ROADMAP-refdata.md`: regenerate all four hg38 artifacts from source with the
committed generators, run SE and PE on them, compare against the P-1 baseline. **This is the last
point at which a refdata defect meets ground truth** — mouse has no reference, so anything
surviving here is invisible from phase 6 onward.

The baseline is trustworthy: `results/se_full` and `results/pe_full` reproduce the CWL 1.0.0
references *exactly* (SE 182/182 read counts, max float deviation 4.6e-12; PE 169/169, 2.0e-13).
Inputs were identical between the two runs — bowtie2 reports the same 25,958,442 reads on both
sides — so every difference below is refdata and nothing else.

### Artifact-level agreement

| artifact | agreement with the 2020 reference |
|---|---|
| `parsed_ucsc_tableformat` | sorted-identical (row order only) |
| `MASTER_FILELIST` | 25,908 / 26,354 ids, recall 0.983 / precision 0.984 |
| index FASTA | 7,494 shared headers, **every shared sequence byte-identical** (match_fraction 1.00000) |
| `UniqueGenomicElements` | recall 0.999917 / precision 0.999913 |

No sequence differs anywhere. No Gencode transcript is mislabelled: all 4,905 shared Gencode
transcripts get the identical family in both files.

### Pipeline-level divergence

| | SE | PE |
|---|---|---|
| total assigned IP reads | 17,505,462 → 17,507,036 (**+0.009%**) | 5,590,588 → 5,591,351 (**+0.014%**) |
| elements | 182 → 183 | 169 → 170 |
| shared elements identical | 74 / 174 | 71 / 161 |
| reads moved between shared elements | 23,629 (**0.135%**) | 8,645 (**0.155%**) |
| reads on dropped / new elements | 79 / 2,792 | 17 / 1,927 |

The families that carry the data are stable. Of the ten largest, nine move by **< 0.05%**:

```
RNA28S             5,399,580 -> 5,398,240   -0.025%
unique_distintron  2,443,909 -> 2,443,792   -0.005%
RNA18S             1,976,076 -> 1,975,587   -0.025%
unique_proxintron  1,017,390 -> 1,017,170   -0.022%
unique_CDS           638,276 ->   638,076   -0.031%
RNA45S               547,301 ->   547,303   +0.000%
antisense_Alu        368,380 ->   365,578   -0.761%
Alu                  271,398 ->   269,691   -0.629%
```

Large *relative* swings exist (`hAT` 807 → 1,609) but only on elements of a few hundred reads,
where a handful of reassigned repeat names moves a visible fraction.

### Where the divergence comes from, and why it is accepted not fixed

**1. rmsk-era drift — 83 family reassignments (38.5% of the moved reads).** The 2020 reference and
today's UCSC `rmsk` disagree about the family of 83 repeat names: `DNA`→`Crypton-A`,
`hAT`→`hAT-Ac`, `ERV3`→`ERVL`. These are genuine reclassifications in the upstream database, not
generator error. **Accepted:** removing them needs a 2020-era `rmsk.txt.gz`, which is not
available; the modern classification is the better one to carry forward.

**2. Gencode membership — 371 reference-only and 267 regenerated-only transcripts (61.5%).**
Dominated by the 559-row unresolved small-RNA residue tracked as `475.41`: transcripts the rmsk
overlap, gene_name and Rfam tiers all fail to name, which the generator omits rather than guesses.
**Accepted:** the omission rule is deliberate (§4 of `docs/MASTER_FILELIST-decisions.md`) — column
4 is a counting label, and a wrong one silently misattributes reads.

**3. A false premise, tested and corrected.** `read_rmsk_class_family` strips a trailing `?` on the
recorded grounds that "the reference carries the settled label." That is wrong — the reference has
53 rows ending in `?`, and the `?` reaches the output because the *mapper* perl (unlike the dedup
perl) does not strip it. Keeping the `?` was implemented and measured: it makes agreement **worse**
(112 reassignments against 93), because 71 repNames in the modern table carry conflicting families,
40 of them a `?`/non-`?` pair, so `first occurrence wins` picks by sort order. The behaviour is
unchanged and the rationale is now the measurement.

### The expectation ceiling for mouse

Mouse cannot be held to a standard human did not meet. Regenerating reference data from source and
rerunning costs, on hg38:

- total assigned reads within **0.02%**
- **≈0.15%** of reads landing on a different family label
- the ten largest families stable to **< 1%**, most to **< 0.05%**
- ~1 element appearing or disappearing per 180, always with a handful of reads

A mouse build inside those bounds is behaving as human does. Outside them, something is wrong that
hg38 would have caught. Comparison artifacts: `tests/p5_regenerated_hg38/`.

---

## 11. repName ambiguity is a per-assembly property (2026-08-14, `-9j2`)

The one finding of P-4 (§`docs/P4-curated-vs-derivable.md` §7) that had teeth: `RMSK_AMBIGUOUS_REPNAMES
= {U5, U17}` in `generate_master_filelist.py` was a human observation applied globally, and it was
the only curated constant measured to actively degrade the mouse artifacts.

hg38's `U5` loci really do span `RNU5A/B/D/E/F` and its `U17` loci span `SNORA` and `RNU105`, so
tier 2 cannot name a family from the repName and must fall through to tier 3 — which always lands,
because human symbols encode the subfamily (`RNU5A-4P` → `RNU5A`). Mouse names the same genes
`Gm24043`, `Gm23102`, `Gm22365`, so tier 3 returns nothing, tier 4 has no U5 family, and the
transcripts are dropped for a conflict mouse does not have: it carries one U5 family (`RNU5G`) and
one U17 family (`SNORA`).

`resolve_ambiguous_repnames` now decides ambiguity from each assembly's own data. A listed repName
whose transcripts resolve to exactly one family downstream is not ambiguous there and tier 2 uses
that family; zero or several leaves the fallthrough unchanged. Overrides count as evidence, so an
operator-supplied second family restores ambiguity. Tier 2 runs over every candidate before any row
is emitted, because the question is a property of the whole assembly.

**Measured**, transcripts overlapping each ambiguous repName and the families tiers 3/4 give them:

| | `U5` | `U17` | verdict |
|---|---|---|---|
| hg38 | 31 → `RNU5A` 8, `RNU5B` 5, `RNU5D` 2, `RNU5E` 9, `RNU5F` 7 | 8 → `SNORA` 5, `RNU105` 2 | ambiguous, unchanged |
| mm10 | 11 → `RNU5G` 1, unresolved 10 | 13 → `SNORA` 9, unresolved 4 | one family each, tier 2 used |
| mm39 | identical to mm10 | identical to mm10 | one family each, tier 2 used |

**hg38 is untouched, verified rather than argued.** The regenerated filelist is *byte-identical* to
P-5's — md5 `2b6d6807fe4ebd7d822562d9cab1155c` — and the tier totals are the §10 ones exactly:
tier 2 4,104 / tier 3 755 / tier 4 298, residue 559.

**Mouse gains 14 rows per assembly and every one is an addition** — no row changed family, none was
lost:

| | MASTER_FILELIST | index FASTA | UniqueGenomicElements | small-RNA residue |
|---|---|---|---|---|
| mm10 | 11,568 → **11,582** | 5,530 → **5,544** | 5,154,593 → **5,154,607** | 785 → **771** |
| mm39 | 25,894 → **25,908** | 5,746 → **5,760** | 5,327,735 → **5,327,749** | 569 → **555** |

The 14 are 10 U5 (`Gm22265`, `Gm24871`, `Gm24043`, `Gm23102`, `Gm22365`, `Gm25099`, `Gm25313`,
`Gm23793`, `Gm23287`, `Gm23143` → `RNU5G`) and 4 U17 (→ `SNORA`), the exact counts `-9j2` predicted.
Mouse U5 now rests on **11** transcripts, not 1.

All four mouse artifacts were rebuilt on the new filelists — the Gencode block feeds the index and
the BED too — and the §8 acceptance suite is **6/6 PASS** on both assemblies, with the
RepBase/Gencode family overlap still `{SNORD}`. The RepBase `U5B1` record stays dropped;
reinstating it would reintroduce the double-count `475.11` removed.

What this does *not* close: the residue is still 555 (mm39) and 771 (mm10) against hg38's 559. The
remaining clone-named `Gm*` small RNAs are `475.41`'s scope.
