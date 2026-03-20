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


def _norm_pathish(value: Any, *, base_dir: Path) -> str:
    norm = _normalize_value(value, base_dir=base_dir)
    if not isinstance(norm, str):
        raise ValueError(f'Expected path-like value, got {type(norm).__name__}')
    p = Path(norm)
    if p.is_absolute():
        return str(p)
    return str((base_dir / p).resolve())


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


def _map_simple_job_to_config(job_data: dict[str, Any], *, base_dir: Path) -> dict[str, Any]:
    mapped: dict[str, Any] = {}

    if 'dataset' in job_data:
        mapped['dataset'] = _normalize_value(job_data['dataset'], base_dir=base_dir)
    if 'prefixes' in job_data:
        mapped['prefixes'] = _normalize_value(job_data['prefixes'], base_dir=base_dir)

    mode = job_data.get('mode') or job_data.get('run_mode') or job_data.get('se_or_pe') or 'SE'
    mode = str(mode).upper()
    mapped['run_mode'] = mode
    mapped['se_or_pe'] = mode

    refs = job_data.get('references', {})
    if refs:
        mapped['bowtie2_db'] = _norm_pathish(refs['bowtie2_db'], base_dir=base_dir)
        mapped['bowtie2_prefix'] = _normalize_value(refs['bowtie2_prefix'], base_dir=base_dir)
        mapped['fileListFile1'] = _norm_pathish(refs['fileListFile1'], base_dir=base_dir)
        mapped['gencodeGTF'] = _norm_pathish(refs['gencodeGTF'], base_dir=base_dir)
        mapped['gencodeTableBrowser'] = _norm_pathish(refs['gencodeTableBrowser'], base_dir=base_dir)
        mapped['repMaskBEDFile'] = _norm_pathish(refs['repMaskBEDFile'], base_dir=base_dir)

    samples = job_data.get('samples', {})
    if not samples:
        return mapped

    ip = samples.get('ip') or samples.get('barcode1') or {}
    inp = samples.get('input') or {}
    ip2 = samples.get('ip2') or samples.get('barcode2') or {}

    if ip:
        mapped['barcode1r1FastqGz'] = _norm_pathish(ip['r1'], base_dir=base_dir)
        if 'r2' in ip:
            mapped['barcode1r2FastqGz'] = _norm_pathish(ip['r2'], base_dir=base_dir)
        mapped['barcode1rmRepBam'] = _norm_pathish(ip['rmrep_bam'], base_dir=base_dir)
    if inp:
        mapped['barcode1Inputr1FastqGz'] = _norm_pathish(inp['r1'], base_dir=base_dir)
        if 'r2' in inp:
            mapped['barcode1Inputr2FastqGz'] = _norm_pathish(inp['r2'], base_dir=base_dir)
        mapped['barcode1InputrmRepBam'] = _norm_pathish(inp['rmrep_bam'], base_dir=base_dir)
    if ip2:
        mapped['barcode2r1FastqGz'] = _norm_pathish(ip2['r1'], base_dir=base_dir)
        mapped['barcode2r2FastqGz'] = _norm_pathish(ip2['r2'], base_dir=base_dir)
        mapped['barcode2rmRepBam'] = _norm_pathish(ip2['rmrep_bam'], base_dir=base_dir)

    return mapped


def map_job_to_config(job_data: dict[str, Any], *, base_dir: str | Path) -> dict[str, Any]:
    based = Path(base_dir)
    if 'samples' in job_data or 'references' in job_data:
        return _map_simple_job_to_config(job_data, base_dir=based)
    return map_cwl_job_to_config(job_data, base_dir=based)


def merge_job_yaml_into_config(config: dict[str, Any], job_yaml: str | Path) -> dict[str, Any]:
    yaml_path = Path(job_yaml)
    job_data = load_yaml(yaml_path)
    mapped = map_job_to_config(job_data, base_dir=yaml_path.parent)
    merged = dict(config)
    merged.update(mapped)
    return merged


def merge_cwl_job_yaml_into_config(config: dict[str, Any], cwl_job_yaml: str | Path) -> dict[str, Any]:
    # Backward-compatible name; now accepts either CWL-style or simplified YAML.
    yaml_path = Path(cwl_job_yaml)
    job_data = load_yaml(yaml_path)
    mapped = map_job_to_config(job_data, base_dir=yaml_path.parent)
    merged = dict(config)
    merged.update(mapped)
    return merged
