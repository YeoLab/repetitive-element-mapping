shell.executable('/bin/bash')

CONDA_MINI = f'{workflow.basedir}/workflow/envs/mini.yaml'
MODULE_INIT = config.get('module_init', 'source /etc/profile.d/modules.sh')
MODULES = config.get('modules', {})


def module_shell(*module_keys: str) -> str:
    names = [MODULES.get(k, '') for k in module_keys if MODULES.get(k, '')]
    if not names:
        return ':'
    return f'{MODULE_INIT} && module purge && module load {" ".join(names)}'


# Mini-profile validation targets for split and merge parity checks.
rule all:
    input:
        'results/mini/split/verified.ok',
        'results/mini/merge/verified.ok',


# Inflate the mini fixture SAM so split parity tests can run on plain SAM input.
rule unzip_mini_sam:
    message:
        'Unzipping mini fixture SAM'
    input:
        config['mini_validation']['source_sam_gz']
    output:
        'results/mini/split/ip.preRmDup.sam.mini.sam'
    conda:
        CONDA_MINI
    shell:
        'gzip -cd {input} > {output}'


# Split the mini SAM with the Python splitter for SE behavior validation.
rule split_mini_sam_python:
    message:
        'Splitting mini fixture SAM with Python implementation'
    input:
        'results/mini/split/ip.preRmDup.sam.mini.sam'
    output:
        directory('results/mini/split/python_tmp')
    params:
        modules=lambda wc: module_shell('samtools'),
    conda:
        CONDA_MINI
    shell:
        r'''
        {params.modules}
        rm -rf {output}
        mkdir -p {output}
        cd {output}
        python {workflow.basedir}/bin/python/split_bam_to_subfiles_SEorPE.py ../ip.preRmDup.sam.mini.sam SE
        '''


# Confirm split output manifest parity against expected fixture listing.
rule verify_split_against_manifest:
    message:
        'Verifying split output manifest against expected fixture'
    input:
        split_dir='results/mini/split/python_tmp',
        manifest=config['mini_validation']['split_manifest'],
    output:
        'results/mini/split/verified.ok'
    conda:
        CONDA_MINI
    script:
        '../scripts/verify_split_manifest.py'


# Merge parsed fixture shards with the Python merger implementation.
rule merge_parsed_python:
    message:
        'Merging parsed fixture shards with Python implementation'
    input:
        config['mini_validation']['merge_input_1'],
        config['mini_validation']['merge_input_2'],
    output:
        'results/mini/merge/merged.python.parsed'
    conda:
        CONDA_MINI
    shell:
        'python {workflow.basedir}/bin/python/merge_multiple_parsed_files.simplified_20191022.py {output} {input[0]} {input[1]}'


# Validate merged parsed output against expected fixture values.
rule verify_merge_against_expected:
    message:
        'Verifying merged parsed output against expected fixture'
    input:
        produced='results/mini/merge/merged.python.parsed',
        expected=config['mini_validation']['merge_expected'],
    output:
        'results/mini/merge/verified.ok'
    conda:
        CONDA_MINI
    script:
        '../scripts/verify_merge_against_expected.py'
