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
