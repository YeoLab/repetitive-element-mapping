rule map_rep_se:
    """
    Map trimmed SE FASTQ to the repeat element database.
    parse_bowtie2_output_realtime_includemultifamily_SE.pl args:
      1: r1 fastq.gz
      2: bowtie2_db/prefix  (directory + "/" + prefix, no extension)
      3: output rep.sam path
      4: fileListFile1 (5-column TSV mapping ENST/ENSG to repeat families)
    Bowtie2 is invoked internally by the Perl script with -p 3.
    """
    input:
        r1       = lambda wc: config[wc.sample]["r1"],
        filelist = config["reference"]["file_list"],
    output:
        rep_sam  = "{outdir}/{sample}/mapped/rep.sam",
    params:
        db_path  = config["reference"]["bowtie2_db"] + "/" + config["reference"]["bowtie2_prefix"],
        script   = os.path.join(workflow.basedir,
                       "workflow/scripts/map_repetitive_elements_se.py"),
        out_abs  = lambda wc, output: os.path.abspath(output.rep_sam),
    conda:
        "../envs/dropin.yaml"
    resources:
        mem_mb  = lambda wc, attempt: attempt * 16000,
        runtime = 480,
    threads: 4
    log:
        "{outdir}/{sample}/logs/map_rep_se.log",
    shell:
        """
        mkdir -p $(dirname {output.rep_sam})
        python {params.script} \
            {input.r1} \
            {params.db_path} \
            {params.out_abs} \
            {input.filelist} \
            > {log} 2>&1
        """
