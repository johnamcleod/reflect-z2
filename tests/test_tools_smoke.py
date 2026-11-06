"""Smoke tests: the independently built AlVin and CoxIter binaries run and agree on a
known compact case (d = 1 + sqrt2, n = 5: Bugaenko's polytope, compact)."""
import pathlib
import subprocess
import tempfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
ALVIN = ROOT / "tools/src/AlVin/build/alvin"
COXITER = ROOT / "tools/src/CoxIter/build/coxiter"

needs_tools = pytest.mark.skipif(not (ALVIN.exists() and COXITER.exists()),
                                 reason="run tools/build.sh first")


@needs_tools
def test_alvin_and_coxiter_smoke():
    with tempfile.TemporaryDirectory() as td:
        (pathlib.Path(td) / "output").mkdir()
        out = subprocess.run([str(ALVIN), "-k=Q[sqrt 2]", "-qf", "-1-T,1,1,1,1,1"],
                             cwd=td, capture_output=True, text=True, timeout=300).stdout
        assert "Algorithm ended" in out
        files = list((pathlib.Path(td) / "output").glob("*.coxiter"))
        assert len(files) == 1
        graph = files[0].read_text()
        assert graph.splitlines()[0].split() == ["10", "5"]
        res = subprocess.run([str(COXITER), "-c", "-fv"], input=graph,
                             capture_output=True, text=True, timeout=300).stdout
        assert "Cocompact: yes" in res
        assert "Finite covolume: yes" in res
