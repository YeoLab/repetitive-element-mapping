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
| `generate_bowtie2_index.py` | FAIL — rewrite planned, not implemented | invalid |
| `generate_master_filelist.py` | FAIL — depends on T-05 | invalid |

Remaining order of work: T-05 (alias resolution → ≥99% coverage gate → rewrite) → T-09
(mouse MASTER_FILELIST) → regenerate mm10/mm39 UniqueGenomicElements → integration test.
