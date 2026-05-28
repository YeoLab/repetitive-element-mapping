# Mini Fixture Set

## Purpose
These mini fixtures provide deterministic, small test data derived from real `PPA1_rep1` outputs so we can validate Perl-to-Python and CWL-to-Snakemake conversion steps quickly.

## Layout
- Source mini inputs:
  - `tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz`
  - `tests/fixtures/mini/source/input.preRmDup.sam.mini.gz`
- Expected mini outputs:
  - `tests/fixtures/mini/expected/ip.rmDup.sam.mini.gz`
  - `tests/fixtures/mini/expected/input.rmDup.sam.mini.gz`
  - `tests/fixtures/mini/expected/ip.parsed.mini.gz`
  - `tests/fixtures/mini/expected/input.parsed.mini.gz`
  - `tests/fixtures/mini/expected/ip.reparsed.nopipes.mini.tsv.gz`
  - `tests/fixtures/mini/expected/ip.reparsed.withpipes.mini.tsv.gz`
- Checksum manifest:
  - `tests/fixtures/mini/manifest.tsv`

## Generation
```bash
python3 scripts/make_mini_fixtures.py
```

Defaults:
- SAM-like files: first 3000 lines
- Parsed files: all `#READINFO` + first 60 `TOTAL` + first 200 `ELEMENT`
- TSV files: header + first 300 data rows

## Validation
```bash
python3 -m pytest -q tests/test_mini_fixtures_manifest.py
```
