import os

# ---------------------------------------------------------------------------
# eCLIP Repetitive Element Mapping — Snakemake 9.12.0
#
# Usage:
#   snakemake --configfile examples/repeat_mapping_SE_small.yaml
#   snakemake --configfile examples/repeat_mapping_PE_small.yaml --profile profiles/tscc2_snakemake9
# ---------------------------------------------------------------------------

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PERL = "/tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl"

# Perl 5.10.1 is a Singularity container wrapper; load singularitypro for all rules
shell.prefix("module load singularitypro; ")
PYTHON_ECLIP = "/tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/ecliprepmap-0.1.0/bin/python"
PREFIXES = [
    "AA", "AC", "AG", "AT", "AN",
    "CA", "CC", "CG", "CT", "CN",
    "GA", "GC", "GG", "GT", "GN",
    "TA", "TC", "TG", "TT", "TN",
    "NA", "NC", "NG", "NT", "NN",
]

# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

SE_OR_PE = config.get("se_or_pe", "SE")
if SE_OR_PE not in ("SE", "PE"):
    raise ValueError(f"se_or_pe must be 'SE' or 'PE', got: {SE_OR_PE!r}")

OUTPUT_DIR = config.get("output_dir", "results")
DATASET    = config.get("dataset", "dataset")


def _require(key, section):
    if section not in config:
        raise ValueError(f"Config missing required section: {section!r}")
    if key not in config[section] or not config[section][key]:
        raise ValueError(
            f"se_or_pe={SE_OR_PE!r} requires config[{section!r}][{key!r}] to be set"
        )


# barcode1 and input always required
for _s in ("barcode1", "input"):
    for _k in ("r1", "bam"):
        _require(_k, _s)

if SE_OR_PE == "PE":
    if "barcode2" not in config or not config.get("barcode2"):
        raise ValueError("se_or_pe=PE requires a 'barcode2' section in config")
    for _s in ("barcode1", "barcode2", "input"):
        for _k in ("r1", "r2", "bam"):
            _require(_k, _s)
else:  # SE
    if config.get("barcode2"):
        raise ValueError("se_or_pe=SE config must not contain 'barcode2'")
    for _s in ("barcode1", "input"):
        if config.get(_s, {}).get("r2"):
            raise ValueError(f"se_or_pe=SE config[{_s!r}] must not have 'r2'")

# ---------------------------------------------------------------------------
# Derived globals (used by common.smk, SE.smk, PE.smk)
# ---------------------------------------------------------------------------

# For SE: IP is barcode1's parsed file.
# For PE: IP is the merged barcode1+barcode2 combined.parsed (from merge_ip_parsed).
if SE_OR_PE == "PE":
    IP_PARSED = f"{OUTPUT_DIR}/{DATASET}.combined.parsed"
else:
    IP_PARSED = f"{OUTPUT_DIR}/barcode1/{DATASET}.barcode1.parsed"

# ---------------------------------------------------------------------------
# rule all — must be in the root Snakefile (Snakemake 9 bug: rule all in
# included files is not recognized as the default target rule).
# ---------------------------------------------------------------------------

if SE_OR_PE == "SE":
    _SAMPLES = ["barcode1", "input"]
    rule all:
        input:
            expand("{outdir}/{dataset}.nopipes.tsv",   outdir=OUTPUT_DIR, dataset=DATASET),
            expand("{outdir}/{dataset}.withpipes.tsv", outdir=OUTPUT_DIR, dataset=DATASET),
            expand(
                "{outdir}/{sample}/{dataset}.{sample}.rmDup.sam.gz",
                outdir=OUTPUT_DIR, dataset=DATASET, sample=_SAMPLES,
            ),
            expand(
                "{outdir}/{sample}/{dataset}.{sample}.preRmDup.sam.gz",
                outdir=OUTPUT_DIR, dataset=DATASET, sample=_SAMPLES,
            ),

elif SE_OR_PE == "PE":
    _SAMPLES = ["barcode1", "barcode2", "input"]
    rule all:
        input:
            expand("{outdir}/{dataset}.nopipes.tsv",   outdir=OUTPUT_DIR, dataset=DATASET),
            expand("{outdir}/{dataset}.withpipes.tsv", outdir=OUTPUT_DIR, dataset=DATASET),
            expand(
                "{outdir}/{sample}/{dataset}.{sample}.rmDup.sam.gz",
                outdir=OUTPUT_DIR, dataset=DATASET, sample=_SAMPLES,
            ),
            expand(
                "{outdir}/{sample}/{dataset}.{sample}.preRmDup.sam.gz",
                outdir=OUTPUT_DIR, dataset=DATASET, sample=_SAMPLES,
            ),

# ---------------------------------------------------------------------------
# Include shared and mode-specific rules
# ---------------------------------------------------------------------------

include: "workflow/rules/common.smk"

if SE_OR_PE == "SE":
    include: "workflow/rules/SE.smk"
elif SE_OR_PE == "PE":
    include: "workflow/rules/PE.smk"
