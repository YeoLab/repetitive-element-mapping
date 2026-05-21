import os

# ---------------------------------------------------------------------------
# Wildcard constraints — shared by all rules in this file and included files
# ---------------------------------------------------------------------------
wildcard_constraints:
    sample  = "barcode1|barcode2|input",
    prefix  = "|".join(PREFIXES),
    dataset = "[^/]+",


# ---------------------------------------------------------------------------
# prepare_rmrep_sam
# Convert the user-supplied rmRep BAM to SAM with a fixed basename ("rmrep.sam")
# so that split_bam produces predictably-named .tmp files.
# ---------------------------------------------------------------------------
rule prepare_rmrep_sam:
    input:
        bam = lambda wc: config[wc.sample]["bam"],
    output:
        sam = "{outdir}/{sample}/mapped/rmrep.sam",
    conda:
        "../envs/dropin.yaml"
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        "{outdir}/{sample}/logs/prepare_rmrep_sam.log",
    shell:
        """
        mkdir -p $(dirname {output.sam})
        samtools view -h {input.bam} > {output.sam} 2>{log}
        """


# ---------------------------------------------------------------------------
# splitbam_rep
# Split the repeat-mapped SAM (rep.sam) by UMI 2-nt prefix into 25 .tmp files.
# The Perl script writes {prefix}.rep.sam.tmp to CWD; we rename to {prefix}.rep.tmp.
# ---------------------------------------------------------------------------
rule splitbam_rep:
    input:
        rep_sam = "{outdir}/{sample}/mapped/rep.sam",
    output:
        expand("{{outdir}}/{{sample}}/split_rep/{prefix}.rep.tmp", prefix=PREFIXES),
    params:
        outdir      = lambda wc: os.path.abspath(f"{wc.outdir}/{wc.sample}/split_rep"),
        rep_sam_abs = lambda wc, input: os.path.abspath(input.rep_sam),
        perl        = PERL,
        script      = os.path.join(workflow.basedir, "bin/perl/split_bam_to_subfiles_SEorPE.pl"),
        se_or_pe    = SE_OR_PE,
        rename      = " && ".join(f"mv {p}.rep.sam.tmp {p}.rep.tmp" for p in PREFIXES),
    resources:
        mem_mb  = lambda wc, attempt: attempt * 8000,
        runtime = 120,
    log:
        "{outdir}/{sample}/logs/splitbam_rep.log",
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        ( cd {params.outdir}
          {params.perl} {params.script} {params.rep_sam_abs} {params.se_or_pe}
          {params.rename}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# splitbam_rmrep
# Split the rmRep SAM (rmrep.sam) by UMI 2-nt prefix into 25 .tmp files.
# The Perl script writes {prefix}.rmrep.sam.tmp; we rename to {prefix}.rmrep.tmp.
# ---------------------------------------------------------------------------
rule splitbam_rmrep:
    input:
        rmrep_sam = "{outdir}/{sample}/mapped/rmrep.sam",
    output:
        expand("{{outdir}}/{{sample}}/split_rmrep/{prefix}.rmrep.tmp", prefix=PREFIXES),
    params:
        outdir        = lambda wc: os.path.abspath(f"{wc.outdir}/{wc.sample}/split_rmrep"),
        rmrep_sam_abs = lambda wc, input: os.path.abspath(input.rmrep_sam),
        perl          = PERL,
        script        = os.path.join(workflow.basedir, "bin/perl/split_bam_to_subfiles_SEorPE.pl"),
        se_or_pe      = SE_OR_PE,
        rename        = " && ".join(f"mv {p}.rmrep.sam.tmp {p}.rmrep.tmp" for p in PREFIXES),
    resources:
        mem_mb  = lambda wc, attempt: attempt * 8000,
        runtime = 120,
    log:
        "{outdir}/{sample}/logs/splitbam_rmrep.log",
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        ( cd {params.outdir}
          {params.perl} {params.script} {params.rmrep_sam_abs} {params.se_or_pe}
          {params.rename}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# dedup
# Deduplicate one UMI prefix bin using duplicate_removal_inline_paired...pl.
# The script writes ALL outputs to CWD named by the basename of arg1:
#   {prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam
#   {prefix}.rep.tmp.combined_w_uniquemap.prermDup.sam
#   {prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt
#   {prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt.done
# We run from a per-prefix subdir to isolate each job's CWD output.
# ALL paths passed to the script must be absolute (because we cd first).
# ---------------------------------------------------------------------------
rule dedup:
    input:
        rep_tmp      = "{outdir}/{sample}/split_rep/{prefix}.rep.tmp",
        rmrep_tmp    = "{outdir}/{sample}/split_rmrep/{prefix}.rmrep.tmp",
        gtf          = config["reference"]["gencode_gtf"],
        table_browser= config["reference"]["gencode_table_browser"],
        rep_bed      = config["reference"]["rep_mask_bed"],
        file_list    = config["reference"]["file_list"],
    output:
        rmdup_sam    = "{outdir}/{sample}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam",
        prermdup_sam = "{outdir}/{sample}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.prermDup.sam",
        parsed_txt   = "{outdir}/{sample}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt",
        done_file    = "{outdir}/{sample}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt.done",
    params:
        dedup_dir      = lambda wc: os.path.abspath(f"{wc.outdir}/{wc.sample}/dedup/{wc.prefix}"),
        rep_tmp_abs    = lambda wc, input: os.path.abspath(input.rep_tmp),
        rmrep_tmp_abs  = lambda wc, input: os.path.abspath(input.rmrep_tmp),
        gtf_abs        = lambda wc, input: os.path.abspath(input.gtf),
        table_abs      = lambda wc, input: os.path.abspath(input.table_browser),
        repbed_abs     = lambda wc, input: os.path.abspath(input.rep_bed),
        filelist_abs   = lambda wc, input: os.path.abspath(input.file_list),
        perl           = PERL,
        script         = os.path.join(workflow.basedir,
                             "bin/perl/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl"),
        se_or_pe       = SE_OR_PE,
    resources:
        mem_mb  = lambda wc, attempt: attempt * 32000,
        runtime = 480,
    threads: 1
    log:
        "{outdir}/{sample}/logs/dedup_{prefix}.log",
    shell:
        """
        mkdir -p {params.dedup_dir}
        ( cd {params.dedup_dir}
          {params.perl} {params.script} \
            {params.rep_tmp_abs} \
            {params.rmrep_tmp_abs} \
            {params.se_or_pe} \
            {params.gtf_abs} \
            {params.table_abs} \
            {params.repbed_abs} \
            {params.filelist_abs}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# concatenate_rmdup
# Cat all 25 per-prefix rmDup SAMs into one per-sample file.
# ---------------------------------------------------------------------------
rule concatenate_rmdup:
    input:
        expand(
            "{{outdir}}/{{sample}}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam",
            prefix=PREFIXES,
        ),
    output:
        "{outdir}/{sample}/{dataset}.{sample}.rmDup.sam",
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        "{outdir}/{sample}/logs/{dataset}.concatenate_rmdup.log",
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} > {output} 2>{log}
        """


# ---------------------------------------------------------------------------
# concatenate_prermdup
# ---------------------------------------------------------------------------
rule concatenate_prermdup:
    input:
        expand(
            "{{outdir}}/{{sample}}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.prermDup.sam",
            prefix=PREFIXES,
        ),
    output:
        "{outdir}/{sample}/{dataset}.{sample}.preRmDup.sam",
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        "{outdir}/{sample}/logs/{dataset}.concatenate_prermdup.log",
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} > {output} 2>{log}
        """


# ---------------------------------------------------------------------------
# gzip_rmdup
# ---------------------------------------------------------------------------
rule gzip_rmdup:
    input:
        "{outdir}/{sample}/{dataset}.{sample}.rmDup.sam",
    output:
        "{outdir}/{sample}/{dataset}.{sample}.rmDup.sam.gz",
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        "{outdir}/{sample}/logs/{dataset}.gzip_rmdup.log",
    shell:
        "gzip -c {input} > {output} 2>{log}"


# ---------------------------------------------------------------------------
# gzip_prermdup
# ---------------------------------------------------------------------------
rule gzip_prermdup:
    input:
        "{outdir}/{sample}/{dataset}.{sample}.preRmDup.sam",
    output:
        "{outdir}/{sample}/{dataset}.{sample}.preRmDup.sam.gz",
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        "{outdir}/{sample}/logs/{dataset}.gzip_prermdup.log",
    shell:
        "gzip -c {input} > {output} 2>{log}"


# ---------------------------------------------------------------------------
# combine_parsed_per_sample
# Merge 25 per-prefix parsed_v2 files into one .parsed per sample.
# merge_multiple_parsed_files.pl accepts absolute paths.
# ---------------------------------------------------------------------------
rule combine_parsed_per_sample:
    input:
        expand(
            "{{outdir}}/{{sample}}/dedup/{prefix}/{prefix}.rep.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt",
            prefix=PREFIXES,
        ),
    output:
        "{outdir}/{sample}/{dataset}.{sample}.parsed",
    params:
        perl        = PERL,
        script      = os.path.join(workflow.basedir,
                          "bin/perl/merge_multiple_parsed_files.simplified_20191022.pl"),
        out_abs     = lambda wc, output: os.path.abspath(output[0]),
        inputs_abs  = lambda wc, input: " ".join(os.path.abspath(f) for f in input),
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 30,
    log:
        "{outdir}/{sample}/logs/{dataset}.combine_parsed.log",
    shell:
        """
        mkdir -p $(dirname {output})
        {params.perl} {params.script} \
            {params.out_abs} \
            {params.inputs_abs} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# fold_change
# Calculate fold enrichment (IP vs Input) from merged parsed files.
# IP_PARSED is set in Snakefile: barcode1.parsed for SE, combined.parsed for PE.
# ---------------------------------------------------------------------------
rule fold_change:
    input:
        ip_parsed    = IP_PARSED,
        input_parsed = expand(
            "{outdir}/input/{dataset}.input.parsed",
            outdir=OUTPUT_DIR, dataset=DATASET,
        ),
    output:
        nopipes   = expand("{outdir}/{dataset}.nopipes.tsv",   outdir=OUTPUT_DIR, dataset=DATASET),
        withpipes = expand("{outdir}/{dataset}.withpipes.tsv", outdir=OUTPUT_DIR, dataset=DATASET),
    params:
        python       = PYTHON_ECLIP,
        script       = os.path.join(workflow.basedir, "bin/calculate_fold_change_from_parsed_files.py"),
        nopipes_abs  = lambda wc, output: os.path.abspath(output.nopipes[0]),
        withpipes_abs= lambda wc, output: os.path.abspath(output.withpipes[0]),
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 30,
    log:
        expand("{outdir}/logs/fold_change.log", outdir=OUTPUT_DIR)[0],
    shell:
        """
        mkdir -p $(dirname {params.nopipes_abs})
        {params.python} {params.script} \
            --ip_parsed {input.ip_parsed} \
            --input_parsed {input.input_parsed} \
            --out_file_nopipes {params.nopipes_abs} \
            --out_file_withpipes {params.withpipes_abs} \
            > {log} 2>&1
        """
