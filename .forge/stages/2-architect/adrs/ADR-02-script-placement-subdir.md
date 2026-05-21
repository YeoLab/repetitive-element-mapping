# ADR-02: Place New Scripts in `bin/python/refdata_generation/` Subdirectory

Date: 2026-05-14
Status: Accepted

## Context

Existing Python scripts live directly in `bin/python/` (e.g., `_perl_compat.py`, `split_bam_to_subfiles_SEorPE.py`). Adding four new scripts plus a shared utility module to the same flat directory would obscure their relationship.

## Decision

Create `bin/python/refdata_generation/` as a Python package containing all four scripts, the `_shared.py` utility module, an `__init__.py`, the `run_assembly.sh` wrapper, and a `README.md`.

## Alternatives Considered

- **Flat layout in `bin/python/`:** rejected because the four new scripts logically group together and are independent of the Perl-compat shims already in that directory.
- **Top-level `refdata_generation/`:** rejected because the existing repo convention places Python under `bin/python/`.

## Consequences

- Positive: clear ownership boundary; easier to navigate; isolates shared utility (`_shared.py`) from the broader `bin/python/` namespace.
- Negative: callers must use a slightly longer path (e.g., `python bin/python/refdata_generation/generate_master_filelist.py`). Negligible.

## References

- Architecture plan section: §4.1 (File Layout)
- Codemap: `.forge/codemap.md` §2 ("Directory Structure & Key Files")
