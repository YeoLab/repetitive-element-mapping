rule map_rep_pe:
    """
    Map trimmed PE FASTQ pair to the repeat element database.
    parse_bowtie2_output_realtime_includemultifamily_PE.pl args:
      1: r1 fastq.gz
      2: r2 fastq.gz
      3: bowtie2_db/prefix
      4: output rep.sam path
      5: fileListFile1
    """
    input:
        r1       = lambda wc: config[wc.sample]["r1"],
        r2       = lambda wc: config[wc.sample]["r2"],
        filelist = config["reference"]["file_list"],
    output:
        rep_sam  = "{outdir}/{sample}/mapped/rep.sam",
    params:
        db_path  = config["reference"]["bowtie2_db"] + "/" + config["reference"]["bowtie2_prefix"],
        perl     = PERL,
        script   = os.path.join(workflow.basedir,
                       "bin/perl/parse_bowtie2_output_realtime_includemultifamily_PE.pl"),
        out_abs  = lambda wc, output: os.path.abspath(output.rep_sam),
    conda:
        "../envs/dropin.yaml"
    resources:
        mem_mb  = lambda wc, attempt: attempt * 16000,
        runtime = 480,
    threads: 4
    log:
        "{outdir}/{sample}/logs/map_rep_pe.log",
    shell:
        """
        mkdir -p $(dirname {output.rep_sam})
        {params.perl} {params.script} \
            {input.r1} \
            {input.r2} \
            {params.db_path} \
            {params.out_abs} \
            {input.filelist} \
            > {log} 2>&1
        """


rule merge_ip_parsed:
    """
    PE only: merge barcode1 + barcode2 parsed files into a single combined.parsed
    used as the IP input to fold_change.
    """
    input:
        bc1 = expand(
            "{outdir}/barcode1/{dataset}.barcode1.parsed",
            outdir=OUTPUT_DIR, dataset=DATASET,
        ),
        bc2 = expand(
            "{outdir}/barcode2/{dataset}.barcode2.parsed",
            outdir=OUTPUT_DIR, dataset=DATASET,
        ),
    output:
        combined = expand(
            "{outdir}/{dataset}.combined.parsed",
            outdir=OUTPUT_DIR, dataset=DATASET,
        )[0],
    params:
        perl    = PERL,
        script  = os.path.join(workflow.basedir,
                      "bin/perl/merge_multiple_parsed_files.simplified_20191022.pl"),
        out_abs = lambda wc, output: os.path.abspath(output.combined),
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4000,
        runtime = 60,
    log:
        expand("{outdir}/logs/merge_ip_parsed.log", outdir=OUTPUT_DIR)[0],
    shell:
        """
        mkdir -p $(dirname {output.combined})
        {params.perl} {params.script} \
            {params.out_abs} \
            {input.bc1} \
            {input.bc2} \
            > {log} 2>&1
        """
