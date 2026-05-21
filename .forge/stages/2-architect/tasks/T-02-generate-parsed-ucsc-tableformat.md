# T-02: Implement `generate_parsed_ucsc_tableformat.py`

<!-- DEPENDENCIES: T-01 -->
<!-- BLOCKS: T-03, T-04, T-06, T-08 -->

## Goal

Implement the Python script that converts a Gencode GTF to the 11-column `parsed_ucsc_tableformat` file. The script must be assembly-generic and accept both `.gtf` and `.gtf.gz` inputs.

## Files Touched

- CREATE: `bin/python/refdata_generation/generate_parsed_ucsc_tableformat.py`

## Implementation Notes

CLI signature (argparse):
```
--gtf PATH       (required)  .gtf or .gtf.gz
--output PATH    (required)  destination file
```

Algorithm: see Architecture Plan §4.3.

Header line: `#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds`

Sort key: `(gene_id, transcript_id)`. Within each transcript, exons sorted by start.

CDS derivation: track `min(start-1)` and `max(end)` over rows where feature == "CDS"; fall back to `(txStart, txEnd)` if no CDS rows seen.

Coordinate conversion: GTF is 1-based inclusive → BED is 0-based half-open: `start_0 = gtf_start - 1`, `end_0 = gtf_end`.

Exon list formatting: comma-delimited with TRAILING comma, e.g. `100,200,300,`.

Output format: tab-separated, Unix LF, no trailing whitespace, no comment header beyond the column header line.

## Acceptance Criteria

- AC-T02-1: Script runs to exit 0 on `examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf` and produces an output file.
- AC-T02-2: Output file's first line equals exactly `#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds`.
- AC-T02-3: Output has exactly 11 tab-separated columns on every data row.
- AC-T02-4: `exonStarts` and `exonEnds` strings always end with `,` (trailing comma).
- AC-T02-5: Script writes progress (records read, transcripts emitted) to stderr.
- AC-T02-6: Script exits non-zero with a descriptive error if `--gtf` does not exist.
- AC-T02-7: Script handles `--gtf path.gtf.gz` (auto-detected gzip) without error.
- AC-T02-8: Output line count equals `transcripts_in_gtf + 1` (header + one row per transcript).

## Verification

```bash
cd /tscc/projects/ps-yeolab3/bay001/codebase/repetitive-element-mapping
python bin/python/refdata_generation/generate_parsed_ucsc_tableformat.py \
  --gtf examples/inputs/hg38/downloaded/gencode.v33.chr_patch_hapl_scaff.annotation.gtf \
  --output /tmp/test_hg38_parsed.tsv

head -2 /tmp/test_hg38_parsed.tsv | cut -f1-4
wc -l /tmp/test_hg38_parsed.tsv
# Expect: 249044 (matches AC-T02-8 against hg38)
```

## Out of Scope

100%-identical reproduction validation lives in T-03, not here. T-02 only requires the script structure and AC checks above.
