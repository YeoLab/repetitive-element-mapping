rule all:
    input:
        'results/mini/split/verified.ok',
        'results/mini/merge/verified.ok',


rule unzip_mini_sam:
    input:
        config['mini']['source_sam_gz']
    output:
        'results/mini/split/ip.preRmDup.sam.mini.sam'
    shell:
        'gzip -cd {input} > {output}'


rule split_mini_sam_python:
    input:
        'results/mini/split/ip.preRmDup.sam.mini.sam'
    output:
        directory('results/mini/split/python_tmp')
    shell:
        r'''
        rm -rf {output}
        mkdir -p {output}
        cd {output}
        python3 {workflow.basedir}/bin/python/split_bam_to_subfiles_SEorPE.py ../ip.preRmDup.sam.mini.sam SE
        '''


rule verify_split_against_manifest:
    input:
        split_dir='results/mini/split/python_tmp',
        manifest=config['mini']['split_manifest'],
    output:
        'results/mini/split/verified.ok'
    script:
        '../scripts/verify_split_manifest.py'


rule merge_parsed_python:
    input:
        config['mini']['merge_input_1'],
        config['mini']['merge_input_2'],
    output:
        'results/mini/merge/merged.python.parsed'
    shell:
        'python3 {workflow.basedir}/bin/python/merge_multiple_parsed_files.simplified_20191022.py {output} {input[0]} {input[1]}'


rule verify_merge_against_expected:
    input:
        produced='results/mini/merge/merged.python.parsed',
        expected=config['mini']['merge_expected'],
    output:
        'results/mini/merge/verified.ok'
    script:
        '../scripts/verify_merge_against_expected.py'
