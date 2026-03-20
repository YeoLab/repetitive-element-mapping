from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


# Keys accepted from CWL-style job YAML and merged into Snakemake config.
CWL_COMPAT_KEYS = {
    'dataset',
    'pipeline_profile',
    'run_mode',
    'barcode1r1FastqGz',
    'barcode1r2FastqGz',
    'barcode1rmRepBam',
    'barcode2r1FastqGz',
    'barcode2r2FastqGz',
    'barcode2rmRepBam',
    'barcode1Inputr1FastqGz',
    'barcode1Inputr2FastqGz',
    'barcode1InputrmRepBam',
    'bowtie2_db',
    'bowtie2_prefix',
    'fileListFile1',
    'fileListFile2',
    'gencodeGTF',
    'gencodeTableBrowser',
    'repMaskBEDFile',
    'chrM_genelist_file',
    'mirbase_gff3_file',
    'prefixes',
    'se_or_pe',
}


def _normalize_value(value: Any, *, base_dir: Path) -> Any:
    if isinstance(value, dict):
        if 'path' in value and isinstance(value['path'], str):
            candidate = Path(value['path'])
            if candidate.is_absolute():
                return str(candidate)
            return str((base_dir / candidate).resolve())
        return {k: _normalize_value(v, base_dir=base_dir) for k, v in value.items()}
    if isinstance(value, list):
        return [_normalize_value(v, base_dir=base_dir) for v in value]
    return value


def load_yaml(path: str | Path) -> dict[str, Any]:
    yaml_path = Path(path)
    with yaml_path.open('r', encoding='utf-8') as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f'Expected mapping at YAML root: {yaml_path}')
    return data


def map_cwl_job_to_config(job_data: dict[str, Any], *, base_dir: str | Path) -> dict[str, Any]:
    based = Path(base_dir)
    mapped: dict[str, Any] = {}
    for key, value in job_data.items():
        if key in CWL_COMPAT_KEYS:
            mapped[key] = _normalize_value(value, base_dir=based)
    return mapped


def merge_cwl_job_yaml_into_config(config: dict[str, Any], cwl_job_yaml: str | Path) -> dict[str, Any]:
    yaml_path = Path(cwl_job_yaml)
    job_data = load_yaml(yaml_path)
    mapped = map_cwl_job_to_config(job_data, base_dir=yaml_path.parent)
    merged = dict(config)
    merged.update(mapped)
    return merged
