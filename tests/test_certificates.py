"""Every stored certificate passes the standalone verifier; corrupted ones fail."""
import copy
import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import verify_certificate as V  # noqa: E402
from reflect.ring import QI  # noqa: E402


def _bump(s):
    return str(QI.parse(s) + 1)

CERTS = sorted((ROOT / "data/certs").glob("*.json"))


@pytest.mark.parametrize("path", CERTS, ids=lambda p: p.stem)
def test_certificate_verifies(path):
    assert V.verify(json.loads(path.read_text()), check_facets=False)


@pytest.mark.parametrize("path", CERTS[:: max(1, len(CERTS) // 8)], ids=lambda p: p.stem)
def test_certificate_facets(path):
    """Slower: re-run Vinberg to confirm the target vertex roots are facets of P."""
    assert V.verify(json.loads(path.read_text()), check_facets=True)


def _first():
    return json.loads(CERTS[0].read_text())


def _expect_fail(cert):
    with pytest.raises(AssertionError):
        V.verify(cert, check_facets=False)


@pytest.mark.skipif(not CERTS, reason="no certificates yet")
def test_corrupted_entry_fails():
    c = _first()
    c["factors"][0]["g"][0][0] = _bump(c["factors"][0]["g"][0][0])
    _expect_fail(c)


@pytest.mark.skipif(not CERTS, reason="no certificates yet")
def test_corrupted_bracket_fails():
    c = _first()
    c["infinite_order"]["bracket"] = ["1/2", "1"]      # a <= 1 is not allowed
    _expect_fail(c)
    c = _first()
    c["infinite_order"]["charpoly"][0] = _bump(c["infinite_order"]["charpoly"][0])
    _expect_fail(c)


@pytest.mark.skipif(not CERTS, reason="no certificates yet")
def test_wrong_vertex_fails():
    c = _first()
    c["factors"][0]["vertex_to"] = list(reversed(c["factors"][0]["vertex_to"]))
    _expect_fail(c)


@pytest.mark.skipif(not CERTS, reason="no certificates yet")
def test_wrong_d_fails():
    c = _first()
    c["d"] = "2+r2" if c["d"] != "2+r2" else "3+r2"
    _expect_fail(c)


def _stab(n):
    rows = [["0"] * (n + 1) for _ in range(n)]
    for j in range(1, n):
        rows[j - 1][j], rows[j - 1][j + 1] = "-1", "1"
    rows[n - 1][n] = "-1"
    return rows


def _is_general(path):
    c = json.loads(path.read_text())
    return any(f["vertex_from"] != _stab(c["n"]) for f in c["factors"])


GENERAL = [p for p in CERTS if _is_general(p)]


@pytest.mark.parametrize("path", GENERAL, ids=lambda p: p.stem)
def test_general_base_certificate_checks_its_base_facets(path):
    """Certificates whose base vertex is not the v0-stabiliser vertex: the verifier must
    re-establish S_from as facets of P, and a fake base root must be rejected."""
    c = json.loads(path.read_text())
    assert V.verify(c, check_facets=False)       # facets check runs automatically
    bad = copy.deepcopy(c)
    f = bad["factors"][0]
    f["vertex_from"] = [r for r in f["vertex_from"]]
    # replace a base root by a non-facet root of L: twice a facet direction is not primitive
    f["vertex_from"][-1] = [str(QI.parse(x) * 2) for x in f["vertex_from"][-1]]
    with pytest.raises(AssertionError):
        V.verify(bad, check_facets=False)
