from __future__ import annotations

import importlib.util
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SHIM = REPO / 'bin/python/_perl_compat.py'


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location('perl_compat_test', path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_resolve_perl_script_path() -> None:
    mod = load_module(SHIM)
    resolved = mod.resolve_perl_script('split_bam_to_subfiles_SEorPE.pl')
    assert resolved.exists()
    assert resolved.name == 'split_bam_to_subfiles_SEorPE.pl'


def test_missing_perl_script_returns_nonzero() -> None:
    mod = load_module(SHIM)
    rc = mod.main('__definitely_missing__.pl')
    assert rc == 2
