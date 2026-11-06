"""Calibration (d): stored isometries <-d,1,1> ~ Mark's entries re-checked exactly, and the
grid results are consistent with Mark's classification."""
import json
import pathlib

import pytest

from reflect.linalg import inverse, matmul, transpose
from reflect.mark import load
from reflect.params import enumerate_params
from reflect.ring import QI, ZS2

ROOT = pathlib.Path(__file__).resolve().parent.parent
ISO = sorted((ROOT / "data/isometries").glob("*.json"))
ENTRIES = load()
CANDS = {c["d"]: c for c in json.loads((ROOT / "data/candidates.json").read_text())}


@pytest.mark.parametrize("path", ISO, ids=lambda p: p.stem)
def test_isometries(path):
    rec = json.loads(path.read_text())
    d = QI.parse(rec["d"])
    F = [[(-d if i == j == 0 else QI(1) if i == j else QI(0)) for j in range(3)] for i in range(3)]
    assert rec["matches"], "reflective <-d,1,1> must be isometric to some entry of Mark's list"
    for m in rec["matches"]:
        e = ENTRIES[m["entry"]]
        phi = [[QI.parse(x) for x in r] for r in m["phi"]]
        assert all(ZS2.contains(x) for r in phi for x in r)
        pinv = inverse(phi)
        assert all(ZS2.contains(x) for r in pinv for x in r)
        assert matmul(transpose(phi), matmul(F, phi)) == e["G"]


def test_grid_consistent_with_mark():
    """Grid d reflective at n = 2 must be candidates; p47b (certified non-reflective at
    n = 2) need not be; nothing certified non-reflective may match an entry."""
    for p in enumerate_params():
        run = ROOT / "data/runs"
        cert = (ROOT / f"data/certs/{p.label}_n2.json").exists()
        if not cert:
            assert str(p.d) in CANDS, p.label
        iso = ROOT / f"data/isometries/{p.label}.json"
        if cert:
            assert not iso.exists() or not json.loads(iso.read_text())["matches"]
