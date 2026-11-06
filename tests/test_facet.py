"""Calibration (e): for every terminated pair (d, n), (d, n+1) the facet -v_{n+1}^perp of
the (d, n+1) polytope reproduces the (d, n) polytope, with no odd angle at that facet
(docs/monotonicity.md)."""
import json

import pytest

from reflect.facet import check_pair
from reflect.params import enumerate_params
from reflect.run import run_path


def _reflective(d, n):
    p = run_path("Zr2", d, n)
    return p.exists() and json.loads(p.read_text())["status"] == "reflective"


PAIRS = [(p.label, p.d, n) for p in enumerate_params() for n in range(2, 10)
         if _reflective(p.d, n) and _reflective(p.d, n + 1)]


def test_pair_count():
    # unit 2-6, ram2 2-5, p7a 2-5, p7b 2-4, inert3 2-3, p23a 2-3 -> 4+3+3+2+1+1
    assert len(PAIRS) == 14


@pytest.mark.parametrize("label,d,n", PAIRS, ids=lambda x: str(x))
def test_facet_reproduces_lower_dimension(label, d, n):
    assert check_pair(d, n) is None
