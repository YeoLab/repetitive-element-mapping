# Architecture Plan — mm10/mm39 Reference Data Generation

**Stage:** 2-architect
**Status:** READY_FOR_BUILD
**Prepared:** 2026-05-14
**Feature:** Generate `parsed_ucsc_tableformat`, `bowtie2_index/`, `UniqueGenomicElements.{assembly}.bed`, and `MASTER_FILELIST.*.tsv` for mm10 and mm39.

---

## 1. Mission Restatement

Produce four reusable, assembly-generic Python scripts that:
1. Reproduce the existing hg38 reference files at ≥99% similarity (100% for parsed_ucsc_tableformat).
2. Generate the same four reference files for mm10 and mm39 from their respective source inputs.
3. Run within the existing `workflow/envs/dropin.yaml` conda environment with no new dependencies.
4. Drop into the CWL and Snakemake pipelines without any workflow file changes.

All four reference files for mm10 and mm39 land in `examples/inputs/mm10/` and `examples/inputs/mm39/`. The four Python scripts land in `bin/python/refdata_generation/`. No Perl modifications are required.

## 2. Constraints Summary

- Python 3.11 only, conda env `workflow/envs/dropin.yaml`, no new packages.
- Write scope: only `examples/inputs/mm10/`, `examples/inputs/mm39/`, and `bin/python/refdata_generation/`.
- Deterministic output: explicit sort, tab delimiter, LF line endings.
- Coordinate conversion: GTF 1-based inclusive → BED/UCSC 0-based half-open.
- Performance: ≤2 hr per assembly on 20 GB / 4 CPU TSCC node.
- mm39 has no tRNA and no miRNA gff3; scripts log WARNINGs and continue.

## 3. Approaches Evaluated

### Approach A: Pure Python + pybedtools (per-step scripts)

Four independent Python scripts, each owning a single output file type. GTF parsing via streaming text I/O (no external GTF library). FASTA extraction via `pybedtools.BedTool.sequence()`. Subprocess calls limited to `bowtie2-build` and `samtools faidx`.

**Pros:**
- Minimal dependency surface — uses only what is already in `dropin.yaml`.
- Each script is independently testable against hg38 ground truth (decoupled reproduction tests).
- Streaming GTF parser avoids loading 11 M GTF rows into memory.
- Clear 1:1 mapping from script to reference file → easier to debug a single output regression.

**Cons:**
- Hand-rolled GTF parsing requires explicit handling of attribute field quoting and comment lines.
- pybedtools requires a `.fai` index; we must ensure `samtools faidx` is invoked if absent.

**Effort:** ~600-900 LOC across four scripts plus ~150 LOC shared utility.

### Approach B: HTSeq / pyranges-based pipeline (single integrated module)

Use `HTSeq` or `pyranges` for GTF parsing and interval operations; single integrated CLI with `--step` flag selecting which output to generate.

**Pros:**
- Mature GTF parsing handles edge cases (multi-line records, malformed attributes).
- Integrated module shares parsed intermediate state across steps (one GTF read serves all four outputs).

**Cons:**
- **`pyranges` and `HTSeq` are NOT in `workflow/envs/dropin.yaml`** — violates NFR-01.
- Single-CLI design couples four otherwise-independent outputs; failure in one step blocks regeneration of others.
- Larger code surface; harder to validate per-output.

**Effort:** ~400-600 LOC but requires dependency addition (rejected) or substituting back to pybedtools/pandas (collapses into Approach A).

### Approach C: Shell + bedtools/awk pipelines

Implement each step as a Bash script that pipes `zcat → awk → sort → bedtools → bowtie2-build`. Python only as a glue layer for the bowtie2 wrapper.

**Pros:**
- Minimal new code; leverages stable UNIX tools.
- Trivially parallelizable per step.
- Lowest runtime overhead.

**Cons:**
- Complex multi-column GTF attribute parsing in awk is fragile and hard to read.
- Difficult to maintain assembly-agnostic CLI semantics (option parsing in pure shell).
- Cannot meet FR-04/FR-09 deterministic-output guarantees easily (sort behavior locale-dependent).
- Misses Python ecosystem for unit testing.
- The existing pipeline conventions are Python (`bin/python/`) and Perl (`bin/perl/`) — Bash glue would be an outlier.

**Effort:** ~400 LOC shell but with high maintenance cost and risk of locale-specific bugs.

### Selected: Approach A

**Rationale:** Approach A is the only option that:
1. Honors NFR-01 (no new conda packages) without compromise.
2. Matches the existing `bin/python/` style established by `_perl_compat.py`, `split_bam_to_subfiles_SEorPE.py`, and `merge_multiple_parsed_files.simplified_20191022.py`.
3. Decouples the four outputs so each can be regenerated and validated independently against its hg38 ground truth.
4. Is the simplest sufficient design (Approach B over-engineers; Approach C trades clarity for marginal performance).

## 4. High-Level Design

### 4.1 File Layout

```
bin/python/refdata_generation/
├── __init__.py                          # Marks as a package
├── _shared.py                           # Shared utilities (see §4.2)
├── generate_parsed_ucsc_tableformat.py  # Step 1
├── generate_bowtie2_index.py            # Step 2
├── generate_unique_genomic_elements.py  # Step 3
├── generate_master_filelist.py          # Step 4
└── run_assembly.sh                      # FR-10 wrapper (Bash; runs all 4 in order)
```

### 4.2 Shared Utilities (`_shared.py`)

Single small module providing helpers used by all four scripts:

```
def open_maybe_gz(path: str) -> TextIO: ...        # auto-detect .gz and open accordingly
def parse_gtf_attributes(attrs_field: str) -> dict[str, str]: ...
def setup_logger(name: str) -> Logger: ...         # writes to stderr, timestamped
def assert_writable(out_path: str, allowed_prefixes: list[str]) -> None: ...
def faidx_if_missing(fasta: str) -> None: ...      # invoke `samtools faidx` if .fai absent
```

Each of these is ≤30 LOC. They eliminate duplication across the four scripts and centralize the write-scope security check (Security req).

### 4.3 Step 1 — `generate_parsed_ucsc_tableformat.py`

**CLI:**
```
generate_parsed_ucsc_tableformat.py \
  --gtf <path>             # .gtf or .gtf.gz
  --output <path>          # destination file (in examples/inputs/{assembly}/)
  [--include-cds]          # if set, derive cdsStart/cdsEnd from CDS features (default: ON)
```

**Algorithm:**
1. Stream-read GTF with `open_maybe_gz()`; skip `#`-prefixed lines.
2. For each row:
   - feature == `transcript`: record `(gene_id, transcript_id, chrom, strand, txStart-1, txEnd)` keyed by `transcript_id`.
   - feature == `exon`: append `(start-1, end)` to `exons[transcript_id]`.
   - feature == `CDS`: track min start and max end per `transcript_id` for `cdsStart`/`cdsEnd`.
3. For each transcript: assemble row. If no CDS feature observed, fall back to `cdsStart=txStart`, `cdsEnd=txEnd`.
4. Sort rows by `(gene_id, transcript_id)`.
5. Within each transcript, sort exons by start position.
6. Write header `#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds` then data rows.
7. Exon `exonStarts`/`exonEnds` use comma-delimited list with trailing comma.

**Validation loop:** Run against hg38 GTF, `diff` against `examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat`, iterate until empty diff before producing mouse outputs.

### 4.4 Step 2 — `generate_bowtie2_index.py`

**CLI:**
```
generate_bowtie2_index.py \
  --gtf <path>                 # for transcript IDs (cross-check)
  --parsed-ucsc <path>         # transcript coordinate source
  --repeatmasker <path.gz>     # required
  --simplerepeats <path.gz>    # required
  --fasta <path>               # genome FASTA, required
  --trna <path.gz>             # optional
  --gff3 <path>                # optional, miRNA
  --custom-fasta <path>        # optional, repeatable (NR_046233.2.fasta etc.)
  --output-dir <dir>           # bowtie2_index/ destination
  --output-prefix <prefix>     # base name for .fa and .bt2 files
```

**Algorithm:**
1. Ensure `<fasta>.fai` exists via `faidx_if_missing()`.
2. **Gencode transcript sequences:** Read `parsed_ucsc_tableformat`. For each transcript build a BED (chrom, exonStart, exonEnd, transcript_id, 0, strand) per exon. Group exons by transcript and use `pybedtools.BedTool().sequence(fi=fasta, s=True, name=True, split=True)` to extract spliced transcript sequence. Header: `>{transcript_id}`.
3. **Repeat element sequences:** Read repeatmasker TSV (gzipped 9-col). Convert to BED. Extract sequences with pybedtools getfasta (strand-aware). Header naming: use `gene_id` from attributes; disambiguate duplicates by appending `_dup1`, `_dup2`, etc. (matching hg38 convention — verified by inspecting hg38 `MASTER_FILELIST.*.fa` headers).
4. **Simple repeats:** Read simplerepeats TSV. Use `transcript_id` (not `gene_id`) for naming, per Edge Case 9.
5. **tRNA (if `--trna` given):** same pattern as repeatmasker; header from `gene_id`.
6. **miRNA (if `--gff3` given):** parse gff3 `Name` attribute; extract sequences; header is the miRNA ID (e.g., `mmu-let-7a-5p`).
7. **Custom FASTAs (`--custom-fasta`):** for `NR_046233.2.fasta`, append verbatim. If multiple records inside the FASTA represent 18S/28S/45S separately, preserve their headers as-is. If the FASTA is a single record but rRNA sub-region naming is required, the implementer must split as documented in NR_046233.2 file (see §6 below).
8. Concatenate all sequences (in this order: Gencode → repeatmasker → simplerepeats → tRNA → miRNA → custom) to `<output-dir>/<output-prefix>.fa`.
9. Index missing IDs: any transcript/element that pybedtools could not extract → count and compare to total.
   - If > 1% missing: write `<output-dir>/missing_ids_report.txt` and log WARNING.
   - If ≤ 1%: log count to stderr only.
10. Run `subprocess.run(["bowtie2-build", "<output.fa>", "<output-prefix>"], cwd=<output-dir>, check=True)`.
11. Verify with `subprocess.run(["bowtie2-inspect", "--summary", "<output-prefix>"], ...)`. Exit non-zero on failure.

### 4.5 Step 3 — `generate_unique_genomic_elements.py`

**CLI:**
```
generate_unique_genomic_elements.py \
  --repeatmasker <path.gz>          # required
  --simplerepeats <path.gz>         # optional
  --trna <path.gz>                  # optional
  --gff3 <path>                     # optional
  --parsed-ucsc <path>              # optional (Gencode transcripts as "elements")
  --assembly <name>                 # e.g. mm10 / mm39 / hg38
  --output <path>                   # destination BED file
  [--flank 500]                     # proximal flank size (default 500)
```

**Algorithm:**
1. Initialize empty list of BED rows.
2. For each provided input source, parse and produce 6-column BED rows:
   - **repeatmasker / trna**: `(chrom, start-1, end, gene_id, score, strand)`.
   - **simplerepeats**: `(chrom, start-1, end, transcript_id, score, strand)` (per Edge Case 9).
   - **parsed_ucsc_tableformat**: `(chrom, txStart, txEnd, name, 0, "-")` — strand hardcoded to "-" per hg38 convention.
   - **gff3 miRNA**: parse `Name=` attribute; `(chrom, start-1, end, name, 0, strand)`.
3. For each element, emit two additional proximal rows:
   - Upstream: `(chrom, max(0, start-flank), start, name + "-proximal", 0, strand)`.
   - Downstream: `(chrom, end, end+flank, name + "-proximal", 0, strand)`.
4. Sort all rows by `(chrom, start, end, name)` — verify against hg38 sort by sampling.
5. Write tab-separated 6-column BED to `--output`. No header.
6. For each optional input that was not provided, log WARNING: `"<source>.tsv not provided; <category> entries omitted from {assembly} UniqueGenomicElements"`.

### 4.6 Step 4 — `generate_master_filelist.py`

**CLI:**
```
generate_master_filelist.py \
  --parsed-ucsc <path>              # required
  --repeatmasker <path.gz>          # required
  --simplerepeats <path.gz>         # optional
  --trna <path.gz>                  # optional
  --gff3 <path>                     # optional
  --custom-fasta <path>             # optional, repeatable
  --gtf <path>                      # for gene_name attribute lookup
  --output <path>                   # destination TSV (and matching .list copy)
```

**Algorithm:**
1. Parse GTF for `transcript_id → gene_name` map (single pass; from `gene_name "..."` attribute).
2. Construct rows in this order:
   1. **Gencode transcripts** (from parsed_ucsc_tableformat):
      `(transcript_id, gene_id, gene_name, "Gencode", "genelists.Gencode")`.
   2. **RepeatMasker repeats** (deduplicated by `gene_id` with `_dup_N` suffix when needed):
      `(disambiguated_name, gene_id, gene_id, repeat_family_from_col2, "genelists.{FAMILY}")`. Map `col2` (source field) → family (e.g., `hg38_rmsk` → derive family from attributes).
   3. **Simple repeats** (using `transcript_id` per Edge Case 9):
      `(transcript_id, "trf", "trf", "Simple_repeat", "genelists.Simple_repeat")`.
   4. **tRNA**:
      `(gene_id, gene_id, gene_id, "tRNA", "genelists.tRNA")`.
   5. **miRNA**: from gff3 `Name=` attribute:
      `(name, name, name, "miRNA", "genelists.miRNA")`.
   6. **Custom FASTA entries** (e.g. NR_046233.2 records):
      `(fasta_header, fasta_header, fasta_header, "rRNA", "genelists.rRNA")` — header value preserved as-is.
3. Deduplicate by col1 (sequence_id) keeping first occurrence.
4. Write tab-separated, no header, to `--output`.
5. Also write identical content to `--output` with `.list` suffix (per OG-04: the `.list` filename is what Perl's `read_in_filelists()` reads).

### 4.7 Step 5 — `run_assembly.sh` (FR-10 wrapper)

A Bash script accepting `--assembly {mm10|mm39|hg38}` and `--date YYYYMMDD`. Resolves paths from `examples/inputs/{assembly}/downloaded/`, invokes the four Python scripts in order with appropriate flags, and writes outputs into `examples/inputs/{assembly}/`. Chosen as Bash (not Snakemake or Python) because:
- It is a thin invocation harness, not a workflow.
- Existing pipeline launchers in `wf/` are Bash; this matches established style.
- No new dependencies.

## 5. Data Models

**Inputs:** Documented in `.forge/stages/1-requirements/context/04-data-models.md` §"Input Source Files".
**Outputs:** Same file, §"Output Reference Files".

This architecture does NOT introduce new file formats; it produces files matching the existing hg38 reference schema verbatim.

## 6. Open Technical Questions Resolved By Architecture

| Question | Decision | Rationale |
|----------|---------|----------|
| Script placement | `bin/python/refdata_generation/` (new subdir) | Groups related new scripts; doesn't clutter `bin/python/` root which holds Perl-compat shims |
| GTF parsing approach | Pure Python streaming | No new dependency; satisfies NFR-01 |
| FASTA extraction approach | pybedtools `BedTool.sequence(s=True)` | Already in env; strand-aware; batched (one call per category) |
| NR_046233.2 sub-region naming | Preserve FASTA headers as-is (don't split) — let the FASTA file itself dictate naming | Files are pre-curated; splitting risks header mismatch against Perl `rRNA_extra_hash` keys |
| `cdsStart`/`cdsEnd` derivation | Use CDS feature rows (min start, max end per transcript); fall back to `txStart`/`txEnd` for non-coding transcripts | Standard UCSC table conversion; matches hg38 pattern |
| Sort order for parsed | `(gene_id, transcript_id)` — verify against hg38 reference and iterate | OG-02 explicitly calls for iterative verification |
| Wrapper form | Bash script | Lightest weight; matches `wf/` launchers |
| Strand-aware getfasta | `s=True` flag | Required for transcripts; repeats can use same flag (strand info is in TSV col7) |
| Chromosome boundary truncation | `max(0, start - flank)` for upstream; no upper-bound cap on downstream (chromosome size lookup is unnecessary because BED tools downstream tolerate over-extent coordinates) | Simplification; aligned with hg38 reference behavior |

## 7. Decision Register (Top Decisions)

See `DECISIONS_SUMMARY.md` and `adrs/` for full ADRs.

| # | Decision | Alternative Rejected | ADR |
|---|----------|---------------------|-----|
| 1 | Pure Python + pybedtools (Approach A) | HTSeq/pyranges (B), Bash/awk (C) | ADR-01 |
| 2 | Place scripts in `bin/python/refdata_generation/` subdir | Place in `bin/python/` root | ADR-02 |
| 3 | hg38 reproduction-test-first iteration | Generate mouse first, validate later | ADR-03 |
| 4 | Preserve NR_046233.2 FASTA headers verbatim | Split into 18S/28S/45S sub-sequences | ADR-04 |
| 5 | Bash wrapper for FR-10 | Python or Snakemake wrapper | ADR-05 |

## 8. Risk Register

| ID | Risk | Likelihood | Impact | Mitigation | Owner |
|----|------|-----------|--------|-----------|-------|
| R-01 | Sort order mismatch breaks AC-01 (100% identical hg38) | HIGH | CRITICAL | Reproduction-test-first; iterate; capture sort key by inspecting hg38 file head | Implementer |
| R-02 | cdsStart/cdsEnd derivation produces wrong values | MEDIUM | HIGH | Test against known-coding transcript (e.g. ENST00000456328); compare to hg38 row | Implementer |
| R-03 | FASTA header naming for repeats does not match hg38 → AC-04 fails | MEDIUM | HIGH | Inspect 100 random hg38 headers, document pattern, reproduce exactly | Implementer |
| R-04 | pybedtools getfasta runtime exceeds 2 hr (NFR-02) | LOW | MEDIUM | Use single batched call per category; ensure `.fai` index present; profile with `time` | Implementer |
| R-05 | Edge Case 4 disambiguation suffix format differs from hg38 | MEDIUM | MEDIUM | Inspect duplicate names in hg38 FASTA; reproduce `_dup{n}` or whatever pattern is observed | Implementer |
| R-06 | mm39 MASTER_FILELIST without tRNA/miRNA breaks downstream Perl | LOW | LOW | AC-13 requires test run; Perl uses hash lookup, missing categories are silent no-ops | Implementer |
| R-07 | The architecture-plan assumes hg38 reference is genuinely deterministic; if the original generator had random behavior, 100% reproduction is impossible | LOW | HIGH | If diff is non-empty after 3 iterations, escalate to user with diff sample | Implementer |
| R-08 | RepElement_pipeline_1dataset.pl users break on mm10/mm39 | LOW | LOW | Document hardcoded `$species` in NOTES; do not modify (out of scope per Edge Case 7) | Documented |

## 9. Trust Boundaries

For a local HPC reference-generation script, the trust boundary surface is minimal. See `threat-model.md` for the STRIDE analysis. Summary:
- **Filesystem-input boundary:** GTF/TSV/FASTA files trusted (read-only group files from internal sources).
- **Subprocess boundary:** `bowtie2-build`, `bowtie2-inspect`, `samtools faidx` — all trusted module-loaded binaries.
- **Filesystem-output boundary:** Scripts must enforce write scope (`examples/inputs/{mm10,mm39}/` only) via `assert_writable()` helper.

## 10. Task Decomposition

See `tasks/` directory. Twelve tasks total, organized as:

```
T-01: shared/_shared.py utility module
T-02: generate_parsed_ucsc_tableformat.py
T-03: generate_parsed_ucsc_tableformat — hg38 reproduction (AC-01)
T-04: generate_bowtie2_index.py
T-05: generate_bowtie2_index — hg38 reproduction (AC-04)
T-06: generate_unique_genomic_elements.py
T-07: generate_unique_genomic_elements — hg38 reproduction (AC-07)
T-08: generate_master_filelist.py
T-09: generate_master_filelist — hg38 reproduction (AC-10)
T-10: Generate mm10 references (AC-02, AC-05, AC-08, AC-11)
T-11: Generate mm39 references (AC-03, AC-06, AC-09, AC-12)
T-12: Pipeline integration & dry-run validation (AC-13, AC-14, AC-16) + Bash wrapper (FR-10)
```

Each task has explicit acceptance criteria, dependency declarations, file touchpoints, and verification commands.

## 11. Verification Strategy

| Layer | Method | Run when |
|-------|--------|---------|
| Unit | Per-script: assert exit 0 on hg38 inputs; compare line counts and column counts | After each script lands |
| Integration | `diff` against hg38 ground truth file (parsed: line-for-line; others: ≥99% line count) | After each reproduction task |
| Pipeline | `snakemake -s workflow/rules/dropin_repelement.smk --config ... -n` for mm10 | T-12 |
| Perl compat | Pipe small test FASTQ through bowtie2 + parse_bowtie2 SE with mm10 references | T-12 |

## 12. Out of Scope

- Reference generation for assemblies other than mm10, mm39, hg38.
- Modification of CWL / Snakemake workflow files.
- Modification of Perl scripts (no hardcoded values block our work).
- Updating `RepElement_pipeline_1dataset.pl`'s hardcoded `$species` (legacy script, not in active workflow).
- New output formats or new metrics.
- Adding conda packages to `dropin.yaml`.
