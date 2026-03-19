from __future__ import annotations

from pathlib import Path

from workflow.config_adapter import map_cwl_job_to_config, merge_cwl_job_yaml_into_config


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
