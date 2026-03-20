from __future__ import annotations

from pathlib import Path

from workflow.config_adapter import (
    map_cwl_job_to_config,
    map_job_to_config,
    merge_cwl_job_yaml_into_config,
    merge_job_yaml_into_config,
)


def test_map_cwl_job_to_config_normalizes_path_objects() -> None:
    base = Path('tests/fixtures/mini')
    mapped = map_cwl_job_to_config(
        {
            'barcode1r1FastqGz': {'class': 'File', 'path': 'source/ip.preRmDup.sam.mini.gz'},
            'bowtie2_db': {'class': 'Directory', 'path': 'refs/bowtie2_index'},
            'se_or_pe': 'SE',
            'unrelated': 123,
        },
        base_dir=base,
    )

    assert mapped['barcode1r1FastqGz'].endswith('tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz')
    assert mapped['bowtie2_db'].endswith('tests/fixtures/mini/refs/bowtie2_index')
    assert mapped['se_or_pe'] == 'SE'
    assert 'unrelated' not in mapped


def test_merge_cwl_job_yaml_into_config_uses_fixture() -> None:
    merged = merge_cwl_job_yaml_into_config(
        {'dataset': 'base_dataset', 'mini_validation': {'source_sam_gz': 'x'}},
        Path('tests/fixtures/mini/cwl_job_example.yaml'),
    )

    assert merged['dataset'] == 'fixture_dataset'
    assert merged['barcode1r1FastqGz'].endswith('tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz')
    assert merged['barcode1rmRepBam'] == '/abs/path/ip.bam'
    assert merged['mini_validation']['source_sam_gz'] == 'x'


def test_map_simple_job_to_config_se() -> None:
    base = Path('tests/fixtures/mini')
    mapped = map_job_to_config(
        {
            'mode': 'SE',
            'dataset': 'simple_ds',
            'samples': {
                'ip': {'r1': 'source/ip.preRmDup.sam.mini.gz', 'rmrep_bam': 'source/ip.preRmDup.sam.mini.gz'},
                'input': {
                    'r1': 'source/input.preRmDup.sam.mini.gz',
                    'rmrep_bam': 'source/input.preRmDup.sam.mini.gz',
                },
            },
            'references': {
                'bowtie2_db': 'refs/bowtie2_index',
                'bowtie2_prefix': 'mock',
                'fileListFile1': 'refs/filelist.tsv',
                'gencodeGTF': 'refs/gencode.gtf',
                'gencodeTableBrowser': 'refs/gencode.table.tsv',
                'repMaskBEDFile': 'refs/repmask.bed',
            },
        },
        base_dir=base,
    )

    assert mapped['run_mode'] == 'SE'
    assert mapped['se_or_pe'] == 'SE'
    assert mapped['barcode1r1FastqGz'].endswith('tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz')
    assert mapped['barcode1Inputr1FastqGz'].endswith('tests/fixtures/mini/source/input.preRmDup.sam.mini.gz')


def test_merge_job_yaml_into_config_with_simple_fixture() -> None:
    merged = merge_job_yaml_into_config(
        {'dataset': 'base_dataset', 'mini_validation': {'source_sam_gz': 'x'}},
        Path('tests/fixtures/dropin/simple_se_job.yaml'),
    )
    assert merged['dataset'] == 'simple_se'
    assert merged['run_mode'] == 'SE'
    assert merged['barcode1r1FastqGz'].endswith('tests/fixtures/dropin/mock.fastq.gz')
    assert merged['mini_validation']['source_sam_gz'] == 'x'
