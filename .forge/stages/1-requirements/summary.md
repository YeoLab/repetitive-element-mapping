# Stage Summary: 1-requirements

**Completed:** 2026-05-14

## What Was Done
Translated `prompts/generate_refdata.md` into a structured 11-file context package for the architect. 
Inventoried existing hg38 reference files and mm10/mm39 source files.

## Key Findings
- **4 Python scripts** to be written (one per reference file type): generate_parsed_ucsc_tableformat.py, generate_bowtie2_index.py, generate_unique_genomic_elements.py, generate_master_filelist.py
- **No Perl modifications needed** — parse scripts already accept MASTER_FILELIST via ARGV
- **Critical gap (OG-01):** mm10 and mm39 genome FASTAs missing — scripts must accept --fasta as required CLI arg
- **mm39 has no tRNA or miRNA source files** — scripts must handle these as optional
- **All scripts must be assembly-generic** (no hardcoded assembly names)

## Acceptance Criteria
16 ACs defined. CRITICAL: AC-01 (hg38 reproduction 100% identical), AC-05/06 (bowtie2 indexes functional), AC-08 (mm10 UniqueGenomicElements), AC-11/12 (mm10/mm39 MASTER_FILELIST with NR_046233.2).

## Artifacts
- `architect-prompt.md` (269 lines, includes all requirements, ACs, gaps)
- `context/01-11` (11 context files covering vision, data models, business logic, edge cases, ACs)
