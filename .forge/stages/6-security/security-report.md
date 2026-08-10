# Security Report: Perl-to-Python Translation

**Date**: 2026-05-22  
**Scope**: `workflow/scripts/{split_bam_to_subfiles,merge_parsed_files,map_repetitive_elements_se,map_repetitive_elements_pe,deduplicate}.py`  
**Tool**: bandit 1.9.4 + manual review  
**Gate result**: PASS_WITH_NOTES

---

## Summary

No critical or high security issues found. Bandit reports 0 medium/high issues across 1,832 lines of code. All subprocess calls use list arguments (no shell injection possible). Three low-severity resource management issues documented below — acceptable for a trusted-environment HPC pipeline.

---

## OWASP Top 10 Coverage

| Category | Status | Notes |
|---|---|---|
| A01 Broken Access Control | N/A | No auth, no access control |
| A02 Cryptographic Failures | N/A | No crypto |
| A03 Injection | **PASS** | All Popen calls use list args, no shell=True |
| A04 Insecure Design | N/A | CLI tool, trusted input |
| A05 Security Misconfiguration | N/A | No server config |
| A06 Vulnerable Components | LOW | conda env pins versions; routine updates apply |
| A07-A10 | N/A | No web, auth, or logging surface |

---

## Findings

### LOW-001: File descriptor not explicitly closed in Popen stderr
**Files**: `map_repetitive_elements_pe.py:263`, `map_repetitive_elements_se.py:272`  
**Pattern**: `stderr=open(bowtie_out, "w")` passed directly to Popen without a context manager.  
**Risk**: File handle leaks until GC runs. Not a security issue; a resource hygiene issue.  
**Recommendation**: Use `with open(bowtie_out, "w") as err_fh: proc = Popen(..., stderr=err_fh)`. Not blocking.

### LOW-002: samtools subprocess not waited in split_bam_to_subfiles.py
**File**: `split_bam_to_subfiles.py:56-61`  
**Pattern**: `proc = Popen(["samtools", ...])` with `infile = proc.stdout` — no `proc.wait()` after reading.  
**Risk**: Zombie process on Linux until parent exits. Minimal impact for short-lived Snakemake rules.  
**Recommendation**: Add `proc.wait()` after the read loop. Not blocking.

### LOW-003: Open file handle without context manager in deduplicate.py
**File**: `deduplicate.py:631`, `deduplicate.py:710`  
**Pattern**: `fh = open(sam_file)` in the `.sam`/`.tmp` branch, no `with` statement.  
**Risk**: Descriptor leak on exception. File closes on GC; no data loss risk.  
**Recommendation**: Wrap in `with open(sam_file) as fh:`. Not blocking.

### LOW-004: Partial executable paths (bandit B607)
**Files**: All scripts using `"samtools"`, `"bowtie2"`, `"stdbuf"` without full paths.  
**Pattern**: Standard for conda-env tools. PATH is controlled by Snakemake's `--use-conda` activation.  
**Risk**: PATH hijacking if the conda env is compromised. Accepted in HPC bioinformatics.  
**Recommendation**: None — full paths would require hardcoding conda env location.

---

## Non-Findings (explicitly verified)

- **No shell=True** anywhere in any script
- **No f-string interpolation into shell commands** (all Popen args are list literals with variable elements)
- **No eval/exec/os.system**
- **Argument count validated** in all scripts (sys.argv length checked, usage printed on error)
- **se_or_pe validated** in split_bam_to_subfiles.py — invalid value prints fatal error and exits
- **No network calls**, no credentials, no secrets handling
- **No SQL**, no template rendering, no deserialization of untrusted data

---

## Conclusion

**Gate**: PASS_WITH_NOTES  
Three low-severity resource management issues. None exploitable in the trusted Snakemake/HPC environment. Recommend addressing LOW-001 and LOW-003 in a future cleanup pass; not blocking for this translation milestone.
