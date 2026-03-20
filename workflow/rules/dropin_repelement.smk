RUN_MODE = (config.get('run_mode') or config.get('se_or_pe') or 'SE').upper()
SE_OR_PE = (config.get('se_or_pe') or RUN_MODE).upper()
PREFIXES = config.get('prefixes') or [
    'AA', 'AC', 'AG', 'AT', 'AN',
    'CA', 'CC', 'CG', 'CT', 'CN',
    'GA', 'GC', 'GG', 'GT', 'GN',
    'TA', 'TC', 'TG', 'TT', 'TN',
    'NA', 'NC', 'NG', 'NT', 'NN',
]


def cfg_required(key: str) -> str:
    v = config.get(key, '')
    if not v:
        raise ValueError(f'Missing required config key for dropin profile: {key}')
    return v


def sample_r1(sample: str) -> str:
    if sample == 'barcode1':
        return cfg_required('barcode1r1FastqGz')
    if sample == 'barcode2':
        return cfg_required('barcode2r1FastqGz')
    if sample == 'input':
        return cfg_required('barcode1Inputr1FastqGz')
    raise ValueError(sample)


def sample_r2(sample: str) -> str:
    if sample == 'barcode1':
        return cfg_required('barcode1r2FastqGz')
    if sample == 'barcode2':
        return cfg_required('barcode2r2FastqGz')
    if sample == 'input':
        return cfg_required('barcode1Inputr2FastqGz')
    raise ValueError(sample)


def sample_rmrep(sample: str) -> str:
    if sample == 'barcode1':
        return cfg_required('barcode1rmRepBam')
    if sample == 'barcode2':
        return cfg_required('barcode2rmRepBam')
    if sample == 'input':
        return cfg_required('barcode1InputrmRepBam')
    raise ValueError(sample)


def split_rep_tmp(prefix: str) -> str:
    return f'{prefix}.rep.sam.tmp'


def split_rmrep_tmp(sample: str, prefix: str) -> str:
    return f'{prefix}.rmrep.sam.tmp'


def dedup_parsed(sample: str, prefix: str) -> str:
    return f'results/dropin/{sample}/dedup/{prefix}/{split_rep_tmp(prefix)}.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt'


def dedup_rmdup(sample: str, prefix: str) -> str:
    return f'results/dropin/{sample}/dedup/{prefix}/{split_rep_tmp(prefix)}.combined_w_uniquemap.rmDup.sam'


def dedup_prermdup(sample: str, prefix: str) -> str:
    return f'results/dropin/{sample}/dedup/{prefix}/{split_rep_tmp(prefix)}.combined_w_uniquemap.prermDup.sam'


if RUN_MODE == 'SE':
    rule all:
        input:
            'results/dropin/barcode1/final/preRmDup.sam.gz',
            'results/dropin/input/final/preRmDup.sam.gz',
            'results/dropin/barcode1/final/rmDup.sam.gz',
            'results/dropin/input/final/rmDup.sam.gz',
            'results/dropin/barcode1/final/combined.parsed',
            'results/dropin/input/final/combined.parsed',
            'results/dropin/barcode1/final/nopipes.tsv',
            'results/dropin/barcode1/final/withpipes.tsv',
else:
    rule all:
        input:
            'results/dropin/pe_ip/final/preRmDup.sam.gz',
            'results/dropin/input/final/preRmDup.sam.gz',
            'results/dropin/pe_ip/final/rmDup.sam.gz',
            'results/dropin/input/final/rmDup.sam.gz',
            'results/dropin/pe_ip/final/combined.parsed',
            'results/dropin/input/final/combined.parsed',
            'results/dropin/pe_ip/final/nopipes.tsv',
            'results/dropin/pe_ip/final/withpipes.tsv',


if RUN_MODE == 'SE':
    rule map_repetitive_elements:
        input:
            r1=lambda wc: sample_r1(wc.sample),
            filelist=lambda wc: cfg_required('fileListFile1'),
        output:
            'results/dropin/{sample}/mapped/rep.sam'
        wildcard_constraints:
            sample='barcode1|input'
        params:
            bowtie=lambda wc: f"{cfg_required('bowtie2_db')}/{cfg_required('bowtie2_prefix')}",
        shell:
            (
                'mkdir -p $(dirname {output}) && '
                '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/parse_bowtie2_output_realtime_includemultifamily_SE.py '
                '{input.r1} {params.bowtie} {output} {input.filelist}'
            )
else:
    rule map_repetitive_elements:
        input:
            r1=lambda wc: sample_r1(wc.sample),
            r2=lambda wc: sample_r2(wc.sample),
            filelist=lambda wc: cfg_required('fileListFile1'),
        output:
            'results/dropin/{sample}/mapped/rep.sam'
        wildcard_constraints:
            sample='barcode1|barcode2|input'
        params:
            bowtie=lambda wc: f"{cfg_required('bowtie2_db')}/{cfg_required('bowtie2_prefix')}",
        shell:
            (
                'mkdir -p $(dirname {output}) && '
                '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/parse_bowtie2_output_realtime_includemultifamily_PE.py '
                '{input.r1} {input.r2} {params.bowtie} {output} {input.filelist}'
            )


rule split_repsam:
    input:
        'results/dropin/{sample}/mapped/rep.sam'
    output:
        expand('results/dropin/{{sample}}/split_rep/{prefix}.rep.sam.tmp', prefix=PREFIXES)
    shell:
        r'''
        rm -rf results/dropin/{wildcards.sample}/split_rep
        mkdir -p results/dropin/{wildcards.sample}/split_rep
        cd results/dropin/{wildcards.sample}/split_rep
        {workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/split_bam_to_subfiles_SEorPE.py ../mapped/rep.sam {SE_OR_PE}
        '''


rule prepare_rmrep_sam:
    input:
        lambda wc: sample_rmrep(wc.sample)
    output:
        'results/dropin/{sample}/mapped/rmrep.sam'
    shell:
        r'''
        mkdir -p $(dirname {output})
        if [[ "{input}" == *.bam ]]; then
          samtools view -h {input} > {output}
        else
          cat {input} > {output}
        fi
        '''


rule split_rmrepbam:
    input:
        'results/dropin/{sample}/mapped/rmrep.sam'
    output:
        expand('results/dropin/{{sample}}/split_rmrep/{prefix}.rmrep.sam.tmp', prefix=PREFIXES)
    shell:
        r'''
        rm -rf results/dropin/{wildcards.sample}/split_rmrep
        mkdir -p results/dropin/{wildcards.sample}/split_rmrep
        cd results/dropin/{wildcards.sample}/split_rmrep
        {workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/split_bam_to_subfiles_SEorPE.py ../mapped/rmrep.sam {SE_OR_PE}
        '''


rule deduplicate_prefix:
    input:
        rep='results/dropin/{sample}/split_rep/{prefix}.rep.sam.tmp',
        rmrep='results/dropin/{sample}/split_rmrep/{prefix}.rmrep.sam.tmp',
        gtf=lambda wc: cfg_required('gencodeGTF'),
        table=lambda wc: cfg_required('gencodeTableBrowser'),
        bed=lambda wc: cfg_required('repMaskBEDFile'),
        filelist=lambda wc: cfg_required('fileListFile1'),
    output:
        rmdup='results/dropin/{sample}/dedup/{prefix}/{prefix}.rep.sam.tmp.combined_w_uniquemap.rmDup.sam',
        prermdup='results/dropin/{sample}/dedup/{prefix}/{prefix}.rep.sam.tmp.combined_w_uniquemap.prermDup.sam',
        parsed='results/dropin/{sample}/dedup/{prefix}/{prefix}.rep.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt',
        done='results/dropin/{sample}/dedup/{prefix}/_done.ok'
    shell:
        r'''
        rep_abs=$(realpath {input.rep})
        rmrep_abs=$(realpath {input.rmrep})
        gtf_abs=$(realpath {input.gtf})
        table_abs=$(realpath {input.table})
        bed_abs=$(realpath {input.bed})
        filelist_abs=$(realpath {input.filelist})
        outdir=$(dirname {output.done})
        mkdir -p "$outdir"
        cd "$outdir"
        {workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.py \
          "$rep_abs" "$rmrep_abs" {SE_OR_PE} "$gtf_abs" "$table_abs" "$bed_abs" "$filelist_abs"
        touch _done.ok
        '''


rule combine_parsed_per_sample:
    input:
        dedup_done=expand('results/dropin/{{sample}}/dedup/{prefix}/_done.ok', prefix=PREFIXES),
        parsed=lambda wc: [dedup_parsed(wc.sample, p) for p in PREFIXES],
    output:
        'results/dropin/{sample}/final/combined.parsed'
    shell:
        (
            'mkdir -p $(dirname {output}) && '
            '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/merge_multiple_parsed_files.simplified_20191022.py '
            '{output} {input.parsed}'
        )


rule concatenate_rmdup_per_sample:
    input:
        dedup_done=expand('results/dropin/{{sample}}/dedup/{prefix}/_done.ok', prefix=PREFIXES),
        sam=lambda wc: [dedup_rmdup(wc.sample, p) for p in PREFIXES],
    output:
        'results/dropin/{sample}/final/rmDup.sam'
    shell:
        'mkdir -p $(dirname {output}) && cat {input.sam} > {output}'


rule concatenate_prermdup_per_sample:
    input:
        dedup_done=expand('results/dropin/{{sample}}/dedup/{prefix}/_done.ok', prefix=PREFIXES),
        sam=lambda wc: [dedup_prermdup(wc.sample, p) for p in PREFIXES],
    output:
        'results/dropin/{sample}/final/preRmDup.sam'
    shell:
        'mkdir -p $(dirname {output}) && cat {input.sam} > {output}'


rule gzip_rmdup_per_sample:
    input:
        'results/dropin/{sample}/final/rmDup.sam'
    output:
        'results/dropin/{sample}/final/rmDup.sam.gz'
    shell:
        'gzip -c {input} > {output}'


rule gzip_prermdup_per_sample:
    input:
        'results/dropin/{sample}/final/preRmDup.sam'
    output:
        'results/dropin/{sample}/final/preRmDup.sam.gz'
    shell:
        'gzip -c {input} > {output}'


if RUN_MODE == 'SE':
    rule foldchange:
        input:
            ip='results/dropin/barcode1/final/combined.parsed',
            ctrl='results/dropin/input/final/combined.parsed',
        output:
            nopipes='results/dropin/barcode1/final/nopipes.tsv',
            withpipes='results/dropin/barcode1/final/withpipes.tsv',
        shell:
            (
                '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/calculate_fold_change_from_parsed_files.py '
                '--ip_parsed {input.ip} --input_parsed {input.ctrl} '
                '--out_file_nopipes {output.nopipes} --out_file_withpipes {output.withpipes}'
            )
else:
    rule combine_pe_ip_parsed:
        input:
            dedup_done=expand('results/dropin/barcode1/dedup/{prefix}/_done.ok', prefix=PREFIXES)
            + expand('results/dropin/barcode2/dedup/{prefix}/_done.ok', prefix=PREFIXES),
            parsed=[dedup_parsed('barcode1', p) for p in PREFIXES] + [dedup_parsed('barcode2', p) for p in PREFIXES],
        output:
            'results/dropin/pe_ip/final/combined.parsed'
        shell:
            (
                'mkdir -p $(dirname {output}) && '
                '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/python/merge_multiple_parsed_files.simplified_20191022.py '
                '{output} {input.parsed}'
            )

    rule concatenate_pe_ip_rmdup:
        input:
            dedup_done=expand('results/dropin/barcode1/dedup/{prefix}/_done.ok', prefix=PREFIXES)
            + expand('results/dropin/barcode2/dedup/{prefix}/_done.ok', prefix=PREFIXES),
            sam=[dedup_rmdup('barcode1', p) for p in PREFIXES] + [dedup_rmdup('barcode2', p) for p in PREFIXES],
        output:
            'results/dropin/pe_ip/final/rmDup.sam'
        shell:
            'mkdir -p $(dirname {output}) && cat {input.sam} > {output}'

    rule concatenate_pe_ip_prermdup:
        input:
            dedup_done=expand('results/dropin/barcode1/dedup/{prefix}/_done.ok', prefix=PREFIXES)
            + expand('results/dropin/barcode2/dedup/{prefix}/_done.ok', prefix=PREFIXES),
            sam=[dedup_prermdup('barcode1', p) for p in PREFIXES] + [dedup_prermdup('barcode2', p) for p in PREFIXES],
        output:
            'results/dropin/pe_ip/final/preRmDup.sam'
        shell:
            'mkdir -p $(dirname {output}) && cat {input.sam} > {output}'

    rule gzip_pe_ip_rmdup:
        input:
            'results/dropin/pe_ip/final/rmDup.sam'
        output:
            'results/dropin/pe_ip/final/rmDup.sam.gz'
        shell:
            'gzip -c {input} > {output}'

    rule gzip_pe_ip_prermdup:
        input:
            'results/dropin/pe_ip/final/preRmDup.sam'
        output:
            'results/dropin/pe_ip/final/preRmDup.sam.gz'
        shell:
            'gzip -c {input} > {output}'

    rule foldchange:
        input:
            ip='results/dropin/pe_ip/final/combined.parsed',
            ctrl='results/dropin/input/final/combined.parsed',
        output:
            nopipes='results/dropin/pe_ip/final/nopipes.tsv',
            withpipes='results/dropin/pe_ip/final/withpipes.tsv',
        shell:
            (
                '{workflow.basedir}/.conda-env/bin/python {workflow.basedir}/bin/calculate_fold_change_from_parsed_files.py '
                '--ip_parsed {input.ip} --input_parsed {input.ctrl} '
                '--out_file_nopipes {output.nopipes} --out_file_withpipes {output.withpipes}'
            )
