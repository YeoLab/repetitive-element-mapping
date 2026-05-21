# Threat Model — mm10/mm39 Reference Data Generation

**Stage:** 2-architect
**Prepared:** 2026-05-14
**Scope:** Four Python scripts and one Bash wrapper that read trusted bioinformatics files and write reference data into a controlled directory tree on a private HPC cluster (TSCC).

---

## 1. Context

This is a single-user, single-cluster, no-network workload. There is no authentication, no incoming traffic, and no PII. The traditional attack surface (web request handlers, secrets, user uploads) is absent. The threat model is therefore narrow and focuses on data integrity, write-scope containment, and subprocess safety.

## 2. Trust Boundaries

```
┌─────────────────────────────────────────────────────────────────────┐
│  TSCC user session (user: bay001, group: yeo-group)                 │
│                                                                     │
│   ┌────────────┐    ┌──────────────────────┐    ┌──────────────┐    │
│   │  source    │ →  │  Python scripts in   │ →  │   output     │    │
│   │  files     │    │  refdata_generation/ │    │   directory  │    │
│   │ (read-only)│    │                      │    │ (write-scope)│    │
│   └────────────┘    └──────────┬───────────┘    └──────────────┘    │
│                                │                                    │
│                                ▼                                    │
│                  ┌─────────────────────────┐                        │
│                  │  subprocesses:          │                        │
│                  │   bowtie2-build         │                        │
│                  │   bowtie2-inspect       │                        │
│                  │   samtools faidx        │                        │
│                  └─────────────────────────┘                        │
└─────────────────────────────────────────────────────────────────────┘
```

### Identified Boundaries

| # | Name | From | To | Crossing Data |
|---|------|------|----|---------------|
| TB-1 | Filesystem-input | `examples/inputs/*/downloaded/` | Python scripts | GTF, TSV, FASTA, GFF3 files |
| TB-2 | Filesystem-output | Python scripts | `examples/inputs/{mm10,mm39}/` | Generated reference files |
| TB-3 | Subprocess | Python scripts | `bowtie2-build`, `bowtie2-inspect`, `samtools faidx` binaries | argv-encoded file paths |
| TB-4 | Module/env | `module load ecliprepmap/1.0.0` or conda env | Python script process | PATH, library locations |

## 3. STRIDE per Trust Boundary

### TB-1: Filesystem-Input Boundary

| Threat | Category | Applicability | Notes |
|--------|----------|---------------|-------|
| Spoofing | S | LOW | Files in `examples/inputs/*/downloaded/` are group-readable bioinformatics files. No identity assertion required by the scripts. |
| Tampering | T | LOW | The scripts trust group-readable files. The TSCC group is small and trusted (yeo-group). Out of scope for code-level mitigation. |
| Repudiation | R | N/A | No user actions tracked. |
| Information Disclosure | I | LOW | Source files contain only public bioinformatics data (Gencode, RepeatMasker, miRBase). |
| Denial of Service | D | LOW | Adversary would need group write access; out of scope. A malformed GTF could in theory loop the parser — mitigated by streaming reads and progress logging. |
| Elevation of Privilege | E | N/A | Scripts run as user, no privilege elevation. |

### TB-2: Filesystem-Output Boundary

| Threat | Category | Applicability | Notes |
|--------|----------|---------------|-------|
| Spoofing | S | N/A | |
| Tampering | T | **MEDIUM** | The most realistic risk: a buggy script could overwrite hg38 reference files (ground truth), invalidating future reproduction tests. **DREAD = D5/R10/E1/A5/D10 = 31.** |
| Repudiation | R | LOW | Filesystem journaling preserves history. |
| Information Disclosure | I | N/A | |
| Denial of Service | D | LOW | Disk-full could block writes; routine. |
| Elevation of Privilege | E | N/A | |

**Mitigation:** `_shared.assert_writable()` enforces an allow-list of output prefixes (`examples/inputs/mm10/`, `examples/inputs/mm39/`, `/tmp/`). Every script calls this guard before opening any output stream. AC-T10-8 verifies the guard rejects hg38 paths.

### TB-3: Subprocess Boundary

| Threat | Category | Applicability | Notes |
|--------|----------|---------------|-------|
| Spoofing | S | LOW | PATH binaries are from a trusted module / conda env. |
| Tampering | T | LOW | argv content originates from controlled CLI args; no user-supplied strings concatenated into shells. |
| Repudiation | R | N/A | |
| Information Disclosure | I | N/A | |
| Denial of Service | D | LOW | bowtie2-build can take significant time; out-of-scope behavior. |
| Elevation of Privilege | E | LOW | **DREAD = D2/R3/E5/A1/D3 = 14.** If a malicious file path containing shell metacharacters were passed via CLI and `shell=True` were used, command injection could occur. Mitigated by mandating `subprocess.run([...], shell=False)` list-form invocation. |

**Mitigation:** AC-T04-9 requires list-form subprocess invocation. Code review must reject any `shell=True`.

### TB-4: Module/Env Boundary

| Threat | Category | Applicability | Notes |
|--------|----------|---------------|-------|
| Spoofing | S | LOW | Module is administrator-curated. |
| Tampering | T | LOW | Module loads from controlled paths. |
| Others | — | N/A | Single-user environment. |

## 4. DREAD-Scored Threats (Selected for AC Promotion)

| ID | Threat | D | R | E | A | D | Total | Promoted to AC? |
|----|--------|---|---|---|---|---|-------|-----------------|
| T-OUT-1 | Bug overwrites `examples/inputs/hg38/` reference files | 5 | 10 | 1 | 5 | 10 | **31** | YES → AC-T10-8 |
| T-SUB-1 | Command injection via shell=True subprocess | 2 | 3 | 5 | 1 | 3 | 14 | YES → AC-T04-9 |
| T-IN-1 | Malformed GTF causes infinite loop or OOM | 1 | 3 | 3 | 1 | 3 | 11 | NO (below threshold 30) |

Two threats with DREAD ≥ 30 have been promoted to acceptance criteria (prefix `[SECURITY]`):
- **AC-T10-8:** Scripts refuse to overwrite `examples/inputs/hg38/` (via `assert_writable`).
- **AC-T04-9:** Subprocess calls use list form (`["bowtie2-build", ...]`), never `shell=True`.

## 5. MAESTRO Analysis

**Not applicable.** This project does not contain AI/ML model components. The MAESTRO framework (Model Input, Model Output, Training Pipeline, Inference Pipeline, Human-AI Interface, Agent Autonomy) addresses LLM/ML system risks that do not exist here.

## 6. Residual Risk

The remaining unmitigated risk is operator error: a user invoking `run_assembly.sh --assembly hg38` could regenerate hg38 outputs into the canonical hg38 directory. We mitigate by routing hg38 reproduction tests through `/tmp/` in T-03, T-05, T-07, T-09. Production hg38 outputs remain immutable.

## 7. Out of Scope

- Network-layer threats (no network access).
- Authentication / authorization (HPC-shell-level user trust).
- Cryptographic considerations (no secrets, no PII).
- Supply-chain attacks on pybedtools, pandas, numpy, bowtie2 (delegated to conda env pinning).
