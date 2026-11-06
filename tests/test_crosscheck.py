"""(c) Independent cross-checks: root lists against AlVin, polytope data against CoxIter."""
import pytest

from reflect.bridges import run_alvin, run_coxiter, tools_available
from reflect.ring import QI, ZS2, ZZ
from reflect.vinberg import vinberg

pytestmark = pytest.mark.skipif(not tools_available(), reason="run tools/build.sh first")

CASES = [
    (ZS2, QI(1, 1), 3), (ZS2, QI(1, 1), 5), (ZS2, QI(1, 1), 6),
    (ZZ, QI(5), 6), (ZZ, QI(5), 8), (ZZ, QI(3), 9), (ZZ, QI(3), 12), (ZZ, QI(13), 2),
    (ZZ, QI(11), 4),
]


@pytest.mark.parametrize("ring,d,n", CASES, ids=lambda x: str(x))
def test_root_lists_match_alvin(ring, d, n):
    ours = vinberg(ring, d, n, delta_max=20000)
    assert ours.status == "reflective"
    theirs = run_alvin(ring, d, n, timeout=600)
    assert theirs["ended"]
    assert theirs["signature"] == (n, 1)
    # same set of roots (both normalise to k0 > 0, canonical norm classes)
    assert set(r.k for r in ours.roots) == set(theirs["roots"])


@pytest.mark.parametrize("ring,d,n", CASES, ids=lambda x: str(x))
def test_polytope_matches_coxiter(ring, d, n):
    ours = vinberg(ring, d, n, delta_max=20000)
    cx = run_coxiter(ours.coxiter)
    assert cx["finite_covolume"] == "yes"
    assert cx["cocompact"] == ("yes" if ring is ZS2 or ours.termination["ideal_vertices"] == 0
                               else "no")
    t = ours.termination
    assert cx["f_vector"][-2] == len(ours.roots)                       # facets
    assert cx["f_vector"][0] == t["finite_vertices"] + t["ideal_vertices"]
    assert cx["f_vector"][1] == t["edges"]
    assert int(cx["vertices_at_infinity"]) == t["ideal_vertices"]
