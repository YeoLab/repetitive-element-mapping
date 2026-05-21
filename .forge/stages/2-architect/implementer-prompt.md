# Implementer Prompt — mm10/mm39 Reference Data Generation

**Stage:** 3-implement
**Status:** READY_FOR_IMPLEMENTATION
**Architect handoff date:** 2026-05-14

---

## 1. Mission

Implement four Python scripts (and a Bash wrapper) that generate reference data files for the mm10 and mm39 mouse genome assemblies. The scripts must reproduce the existing hg38 reference files (≥99% similarity; 100% for `parsed_ucsc_tableformat`) and then produce mm10 and mm39 outputs. The pipeline must accept the new references via existing CWL YAML and Snakemake `--config` interfaces without any workflow file changes.

## 2. Constraints

1. Python 3.11, conda env `workflow/envs/dropin.yaml`. **No new conda packages.**
2. Write scope: `examples/inputs/mm10/`, `examples/inputs/mm39/`, `bin/python/refdata_generation/`, and `/tmp/` (for tests). DO NOT modify `examples/inputs/hg38/`, CWL files, Snakemake workflow files, or Perl scripts.
3. Deterministic output: explicit sort orders, tab delimiters, Unix LF line endings.
4. Performance: each assembly's reference generation must complete within 2 hours on a 20 GB / 4 CPU node.
5. mm39 has no tRNA TSV and no miRNA gff3 — handle as optional with WARNINGs.
6. All four scripts must accept both `.gtf` and `.gtf.gz` inputs (auto-detect by extension).

## 3. Reading List (in this order)

1. `.forge/stages/2-architect/architecture-plan.md` — design overview.
2. `.forge/stages/2-architect/tasks/T-01.md` through `T-12.md` — twelve tasks in dependency order.
3. `.forge/stages/1-requirements/context/04-data-models.md` — file formats.
4. `.forge/stages/1-requirements/context/05-business-logic.md` — algorithm specs.
5. `.forge/stages/1-requirements/context/08-edge-cases.md` — gotchas.
6. `.forge/stages/2-architect/assumptions.json` — assumptions you may rely on (or invalidate).
7. `.forge/stages/2-architect/threat-model.md` — security boundaries.

## 4. Execution Order

Tasks form this dependency DAG (run in topological order):

```
T-01 (shared utils)
 ├─→ T-02 (parsed_ucsc script) ──→ T-03 (hg38 parsed reproduction, AC-01)
 ├─→ T-04 (bowtie2 script)    ──→ T-05 (hg38 bowtie2 reproduction, AC-04)
 ├─→ T-06 (unique elements)   ──→ T-07 (hg38 elements reproduction, AC-07)
 └─→ T-08 (master filelist)   ──→ T-09 (hg38 filelist reproduction, AC-10)
                                                    │
                                                    ▼
                              T-10 (generate mm10) ←┤ ←─ T-11 (generate mm39)
                                                    │
                                                    ▼
                                                  T-12 (integration, dry-run, Perl compat)
```

T-03, T-05, T-07, T-09 are reproduction-test gates. **Do not start T-10/T-11 until all four gates pass.**

## 5. Critical Decisions Already Made by Architect

- **Script placement:** `bin/python/refdata_generation/` (new subdir).
- **GTF parsing:** pure-Python streaming (no new deps).
- **FASTA extraction:** `pybedtools.BedTool.sequence(s=True)`.
- **Wrapper:** Bash (`run_assembly.sh`), not Snakemake or Python.
- **`cdsStart`/`cdsEnd`:** derived from CDS feature rows; fall back to `(txStart, txEnd)` for non-coding.
- **`NR_046233.2.fasta`:** append verbatim — do not split or rename headers.
- **`.list` file:** identical content to `.tsv` (Perl reads the `.list` extension).
- **Reproduction-test-first:** AC-01 (100% identical hg38 parsed) is the first gate; iterate on T-02 until diff is empty.

If you find an ADR-recorded decision wrong, **stop and surface to user**; do not silently change architecture.

## 6. Success Definition

All 16 ACs from `.forge/stages/1-requirements/context/09-acceptance-criteria.md` pass, plus all 12 task-level ACs from this stage's `tasks/`.

The fast path to "done":
1. All four scripts in place under `bin/python/refdata_generation/`.
2. All four hg38 reproduction gates green.
3. All four mm10 reference files in `examples/inputs/mm10/`.
4. All four mm39 reference files in `examples/inputs/mm39/`.
5. `snakemake -n` passes for mm10 and mm39.
6. Perl SE parse script accepts mm10 MASTER_FILELIST without error.
7. `run_assembly.sh` available as the user-facing entry point.

## 7. When to Escalate

Escalate to user (do not silently proceed) if:
- Hg38 parsed_ucsc_tableformat diff is non-empty after 3 implementation iterations on T-03.
- Any reproduction-gate similarity drops below 95% (we can't get within 4% of the ≥99% threshold by ordinary debugging).
- A required source file is missing or unreadable (e.g., the mm10.fa symlink target disappears).
- A subprocess (`bowtie2-build`, `samtools faidx`) exits non-zero with an unrecoverable error.
- You discover an assumption from `assumptions.json` is invalid.

---

**End of implementer prompt.** Begin with T-01.
