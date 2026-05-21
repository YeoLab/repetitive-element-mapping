# Codebase Exploration Report — mm10/mm39 Reference Data Generation

**Stage:** 2-architect (sub-phase 2a: explore)
**Prepared:** 2026-05-14
**Scope:** Identify files to modify, integration points, patterns, and constraints for generating mm10 and mm39 reference data.

---

## 1. Files Analyzed

| Category | Path | Role |
|----------|------|------|
| Source prompt | `prompts/generate_refdata.md` | Original user request (5 steps) |
| Requirements pkg | `.forge/stages/1-requirements/architect-prompt.md` + `context/` | FRs, NFRs, ACs, edge cases |
| Codemap | `.forge/codemap.md` | Repo structure, data flow |
| Perl SE parser | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl` | Consumes MASTER_FILELIST via ARGV[3] |
| Perl PE parser | `bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl` | Same ARGV[3] convention |
| Snakemake dropin | `workflow/rules/dropin_repelement.smk` | Consumes references via `--config` |
| Conda env | `workflow/envs/dropin.yaml` | Runtime environment |
| Existing Python ports | `bin/python/_perl_compat.py`, `bin/python/split_bam_to_subfiles_SEorPE.py`, `bin/python/merge_multiple_parsed_files.simplified_20191022.py` | Style reference for new scripts |
| hg38 ground truth | `examples/inputs/hg38/{parsed_ucsc_tableformat, MASTER_FILELIST.*.tsv, MASTER_FILELIST.*.list, UniqueGenomicElements.hg38.bed, bowtie2_index/*.fa}` | Reproduction targets |
| mm10 inputs | `examples/inputs/mm10/downloaded/{*.gtf.gz, *.tsv.gz, mmu.gff3, NR_046233.2.fasta, mm10.fa, mm10.fa.fai}` | Source data |
| mm39 inputs | `examples/inputs/mm39/downloaded/{*.gtf.gz, *.tsv.gz, NR_046233.2.fasta, GRCm39.primary_assembly.genome.fa, *.fa.fai}` | Source data; no tRNA/miRNA |

**Total**: ~25 source files inspected. No new conda packages needed; all required tools (pybedtools, pandas, bowtie2, bedtools, samtools) are present in `workflow/envs/dropin.yaml`.

## 2. Integration Points

1. **Perl `read_in_filelists()` (parse_bowtie2 SE/PE)** — reads `ARGV[3]` MASTER_FILELIST as a 5-column TSV. Indexed by col1 (sequence_id). Hardcoded paths in the script are commented out. **No script modification required.**
2. **Snakemake `dropin_repelement.smk`** — accepts reference paths via `--config` flags: `bowtie2_db`, `bowtie2_prefix`, `fileListFile1`, `gencodeGTF`, `gencodeTableBrowser`, `repMaskBEDFile`. New mm10/mm39 references are drop-in replacements at the config layer.
3. **CWL `wf_ecliprepmap_se.cwl` / `wf_ecliprepmap_pe.cwl`** — accept the same reference paths via YAML input. No CWL file changes.
4. **pybedtools** — `BedTool.sequence(fi=fasta, s=True)` for strand-aware FASTA extraction. Requires `.fa.fai` index (already present for mm10/mm39).
5. **bowtie2-build** — invoked as subprocess from Step-2 script; writes 6 `.bt2` files plus the source `.fa`.

## 3. Constraints Identified

| Constraint | Source | Effect on design |
|-----------|--------|-----------------|
| Write scope limited to `examples/inputs/mm10/` and `examples/inputs/mm39/` | Security req | Scripts must validate `--output` paths |
| No new conda packages | NFR-01 | Use only pybedtools, pandas, numpy, bowtie2 already installed |
| Deterministic output | FR-09, NFR-03 | Must sort all collections explicitly; tab delimiters; Unix LF |
| GTF 1-based inclusive → BED 0-based half-open | Edge Case 6, Data models | `start_0 = gtf_start - 1`, `end_0 = gtf_end` |
| `.gtf.gz` support | FR-05, AC-16 | Auto-detect by extension; use `gzip.open()` |
| hg38 reproduction at ≥99% (parsed = 100%) | AC-01,04,07,10 | Reproduction-test loop is required before mouse outputs are trusted |
| Missing IDs threshold = 1% | FR-08, AC-15 | Write `missing_ids_report.txt` only when >1% missing |
| Performance ≤ 2 hr / 20 GB / 4 CPU | NFR-02 | pybedtools getfasta must use indexed FASTA; avoid per-record subprocess |
| mm39 has no tRNA, no miRNA gff3 | Edge Case 1 | Optional CLI flags; WARNING on absence; never abort |
| Sort order for parsed_ucsc_tableformat must match hg38 | Business logic, OG-02 | Implementer iterates against hg38 ground truth until diff = 0 |
| Perl scripts ≤5.16 | Tech constraints | Not relevant for new Python scripts |

## 4. Patterns To Follow

Observed from existing Python scripts in `bin/python/`:
- **CLI**: `argparse` with `--required` flags; no positional args. All scripts non-interactive.
- **Progress logging**: `print(..., file=sys.stderr)` for status; exit non-zero on error with descriptive message.
- **Path resolution**: `_perl_compat.py` style — derive paths relative to repo root via `Path(__file__).resolve().parents[N]`.
- **Style**: ~Python 3.11, no type hints in existing ports, but new scripts can use them (additive).
- **Testing**: Snakemake-based validation in `workflow/rules/se_foundation.smk`; we extend the same pattern by reproducing hg38 outputs.

Observed from existing Perl scripts:
- All file paths via `ARGV`; hardcoded paths are commented out.
- Output naming follows `<input_basename>.<transformation_suffix>` convention.

## 5. Risk Areas

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|-----------|
| Sort order mismatch breaks 100% hg38 parsed reproduction | HIGH | CRITICAL | Iterate against hg38 ground truth with `diff` until line-for-line identical; capture sort key from inspection of hg38 file |
| `cdsStart`/`cdsEnd` derivation logic ambiguous | MEDIUM | HIGH | Inspect hg38 reference rows for known non-coding transcripts; choose CDS-feature-based derivation; fall back to txStart/txEnd for non-coding |
| FASTA header naming for repeat/tRNA/rRNA differs subtly from hg38 | MEDIUM | HIGH (could break AC-04 ≥99%) | Inspect hg38 `MASTER_FILELIST.*.fa` headers; reproduce naming exactly; especially `NR_046235.3-18S` / `-28S` / `-45S` pattern |
| MASTER_FILELIST row order differs from hg38 | MEDIUM | MEDIUM (AC-10 ±1%) | Preserve grouping: Gencode transcripts → repeats → simple repeats → tRNA → miRNA → NR_ |
| pybedtools getfasta memory/time on large repeatmasker BED | MEDIUM | MEDIUM (NFR-02) | Pre-index FASTA with `samtools faidx`; batch all intervals into single getfasta call |
| Edge Case 4: duplicate element names (AluY) | HIGH | MEDIUM | Inspect hg38 naming convention (likely positional disambiguation); reproduce exactly |
| Missing genome FASTA blocking work | RESOLVED | — | OG-01 marked critical in requirements, but mm10.fa and GRCm39 fa symlinks already present in `downloaded/` |

## 6. Key Findings

1. **OG-01 is resolved**: Both `mm10.fa` and `GRCm39.primary_assembly.genome.fa` symlinks (with `.fai` indices) are already in `examples/inputs/{mm10,mm39}/downloaded/`. No external sourcing needed; scripts pass-through to existing files.
2. **No Perl modifications required**: All three core Perl scripts (`parse_bowtie2_*_SE.pl`, `parse_bowtie2_*_PE.pl`, `duplicate_removal_*.pl`) read MASTER_FILELIST via `$ARGV[3]` with no assembly-specific hardcoding. `RepElement_pipeline_1dataset.pl` is a legacy orchestrator not invoked by CWL/Snakemake.
3. **Sort order is the dominant risk** for AC-01 (100% identical hg38 reproduction). Implementer must run reproduction-tests early and iterate.

## 7. Recommended Script Placement

`bin/python/refdata_generation/` (a new subdirectory) — keeps the four new scripts grouped together and out of the way of the existing Perl-compatible shims in `bin/python/`. Existing pattern: `bin/python/_perl_compat.py` is a sibling to `bin/python/split_bam_to_subfiles_SEorPE.py`. A subdirectory makes the new artifacts independently testable.

## 8. Conclusion

The exploration confirms:
- All input files are present (including the previously-flagged missing genome FASTAs).
- No new conda packages needed; no Perl modifications needed; no CWL/Snakemake changes needed.
- All four scripts can be cleanly implemented as a self-contained Python subpackage.
- The single highest-risk technical decision is `parsed_ucsc_tableformat` sort order and CDS derivation — this drives a reproduction-test-first workflow.

Proceed to design sub-phase.
