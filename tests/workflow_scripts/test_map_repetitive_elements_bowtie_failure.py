"""bowtie2 failure must abort the mapper, not produce an empty rep.sam (issue -475.16).

Regression for the defect that made se_full/pe_full report
'#READINFO RepFamilyReads 0 0': bowtie2 was absent from PATH, stdbuf exited 127,
and the pipe reader saw an immediate EOF. `proc.wait()` discarded the exit code,
so the mapper wrote its '.done' file and exited 0. The failure only became visible
as zero repeat-family reads in the final .nopipes.tsv, and was misdiagnosed for
months as a transient bowtie2 fault.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"

# The real stderr stdbuf writes when bowtie2 is not on PATH. Quoting around the
# command name is locale-dependent (straight vs typographic), so match either.
STDBUF_MISSING_CMD = "failed to run command"


@pytest.fixture
def filelist(tmp_path):
    """Minimal 5-column filelist; read_filelists runs before bowtie2 launches."""
    p = tmp_path / "filelist.tsv"
    p.write_text("ENST00000384010.1\tENSG00000199568.1\tRNU1-1\tRNU1\tgenelists.RNU1\n")
    return p


@pytest.fixture
def fastq(tmp_path):
    p = tmp_path / "reads.fq"
    p.write_text("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
    return p


def run_mapper(script, args, tmp_path, path_env):
    """Run a mapper with a PATH that does or does not contain bowtie2."""
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *args],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env={"PATH": path_env},
    )


@pytest.fixture
def bowtie_free_path(tmp_path):
    """A PATH with coreutils (for stdbuf) but no bowtie2."""
    stub = tmp_path / "bin"
    stub.mkdir()
    for tool in ("stdbuf",):
        real = subprocess.run(["which", tool], capture_output=True, text=True).stdout.strip()
        if not real:
            pytest.skip(f"{tool} not available")
        (stub / tool).symlink_to(real)
    return str(stub)


def test_se_mapper_exits_nonzero_when_bowtie2_missing(
    tmp_path, filelist, fastq, bowtie_free_path
):
    out = tmp_path / "rep.sam"
    result = run_mapper(
        "map_repetitive_elements_se.py",
        ["index_prefix", str(fastq), str(out), str(filelist)],
        tmp_path,
        bowtie_free_path,
    )
    assert result.returncode != 0, (
        "mapper exited 0 despite bowtie2 being absent; this is the -475.16 defect"
    )
    assert STDBUF_MISSING_CMD in result.stderr, (
        f"bowtie2's own stderr must be surfaced, got: {result.stderr!r}"
    )


def test_pe_mapper_exits_nonzero_when_bowtie2_missing(
    tmp_path, filelist, fastq, bowtie_free_path
):
    out = tmp_path / "rep.sam"
    result = run_mapper(
        "map_repetitive_elements_pe.py",
        ["index_prefix", str(fastq), str(fastq), str(out), str(filelist)],
        tmp_path,
        bowtie_free_path,
    )
    assert result.returncode != 0
    assert STDBUF_MISSING_CMD in result.stderr and "bowtie2" in result.stderr


@pytest.mark.parametrize(
    "script,args",
    [
        ("map_repetitive_elements_se.py", ["index_prefix", "reads.fq", "rep.sam", "filelist.tsv"]),
        ("map_repetitive_elements_pe.py",
         ["index_prefix", "reads.fq", "reads.fq", "rep.sam", "filelist.tsv"]),
    ],
)
def test_done_file_not_written_on_bowtie_failure(
    tmp_path, filelist, fastq, bowtie_free_path, script, args
):
    """Downstream rules key off '.done'; writing it on failure is what let the
    empty rep.sam flow through splitbam, dedup and combine as if it were real."""
    run_mapper(script, args, tmp_path, bowtie_free_path)
    assert not (tmp_path / "rep.sam.done").exists(), (
        "'.done' written despite bowtie2 failing — downstream steps will treat "
        "an empty rep.sam as a valid zero-alignment result"
    )
