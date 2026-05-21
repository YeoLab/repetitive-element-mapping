# T-01: Shared Utilities Module

<!-- DEPENDENCIES: none -->
<!-- BLOCKS: T-02, T-04, T-06, T-08 -->

## Goal

Create the `bin/python/refdata_generation/_shared.py` module providing helpers used by all four generation scripts. Eliminates duplication and centralizes write-scope security checks.

## Files Touched

- CREATE: `bin/python/refdata_generation/__init__.py` (empty)
- CREATE: `bin/python/refdata_generation/_shared.py`

## Implementation Notes

Module must expose at minimum:

```python
def open_maybe_gz(path: str | Path) -> TextIO:
    """Open .gtf or .gtf.gz transparently. Uses gzip.open(..., 'rt')."""

def parse_gtf_attributes(attrs_field: str) -> dict[str, str]:
    """Parse a GTF column-9 attribute string into a dict. Handles 'key "value";' format."""

def setup_logger(name: str) -> logging.Logger:
    """Return a logger that writes ISO-timestamped messages to stderr."""

def assert_writable(out_path: str | Path, allowed_prefixes: list[str]) -> None:
    """Raise PermissionError if out_path is not under any of allowed_prefixes.
    Default allowed_prefixes should include examples/inputs/mm10/, examples/inputs/mm39/."""

def faidx_if_missing(fasta: str | Path) -> None:
    """If <fasta>.fai is missing, invoke samtools faidx as subprocess."""
```

## Acceptance Criteria

- AC-T01-1: Module imports cleanly under Python 3.11 inside the dropin conda env.
- AC-T01-2: `open_maybe_gz("test.gtf")` returns a text file; `open_maybe_gz("test.gtf.gz")` returns a decoded text stream (verified via a 5-line round-trip test).
- AC-T01-3: `parse_gtf_attributes('gene_id "ENSG1"; transcript_id "ENST1";')` returns `{"gene_id": "ENSG1", "transcript_id": "ENST1"}`.
- AC-T01-4: `assert_writable("/etc/passwd", ["examples/inputs/"])` raises `PermissionError`.
- AC-T01-5: `faidx_if_missing(path)` is a no-op when `path + ".fai"` exists; runs `samtools faidx` when missing.
- AC-T01-6: All public functions have a one-line docstring.

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
python -c "from bin.python.refdata_generation._shared import open_maybe_gz, parse_gtf_attributes, setup_logger, assert_writable, faidx_if_missing; print('OK')"
python -c "from bin.python.refdata_generation._shared import parse_gtf_attributes as p; assert p('gene_id \"E1\"; transcript_id \"T1\";') == {'gene_id':'E1','transcript_id':'T1'}; print('OK')"
```

## Out of Scope

Type hints are optional. Test framework wiring (pytest) is optional for this task — verification commands above are sufficient.
