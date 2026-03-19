#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_perl_script(perl_script: Path, argv: list[str]) -> int:
    cmd = ['perl', str(perl_script), *argv]
    completed = subprocess.run(cmd)
    return completed.returncode


def resolve_perl_script(perl_basename: str) -> Path:
    here = Path(__file__).resolve().parent
    repo = here.parents[1]
    return repo / 'bin' / 'perl' / perl_basename


def main(perl_basename: str) -> int:
    perl_script = resolve_perl_script(perl_basename)
    if not perl_script.exists():
        print(f'Could not find perl script: {perl_script}', file=sys.stderr)
        return 2
    return run_perl_script(perl_script, sys.argv[1:])
