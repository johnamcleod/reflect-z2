"""Calibration against the literature (must pass before trusting the Z[sqrt2] grid).

(a) K = Q: Mark (p = 5, reflective exactly for 2 <= n <= 8) and Mcleod (p = 3,
    reflective exactly for 2 <= n <= 13), including face counts and, for p = 3, the
    numbers of ideal vertices of Mcleod (arXiv:1007.2299) Corollary 3.
(b) Bugaenko (via Belolipetsky's survey): d = 1+sqrt2 reflective for n <= 6, and the
    algorithm does not terminate for n = 7.
Non-termination is checked "within caps" only; the proof of non-reflectivity is the
business of Phase 3 certificates.
"""
import pytest

from reflect.ring import QI, ZS2, ZZ
from reflect.vinberg import vinberg

# Mcleod p = 3 facet counts: n + (number of rows of his Table 2 present in dim n)
MCLEOD_P3 = {2: 3, 3: 4, 4: 6, 5: 7, 6: 8, 7: 9, 8: 10, 9: 12, 10: 14, 11: 15, 12: 18, 13: 22}
MCLEOD_P3_IDEAL = {3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 2, 10: 3, 11: 5, 12: 6, 13: 13}
# Mark p = 5 facet counts (n = 8: 12 nodes in Mark's Figure 1; others cross-checked by AlVin)
MARK_P5 = {2: 4, 3: 6, 4: 7, 5: 8, 6: 10, 7: 11, 8: 12}
# Bugaenko d = 1 + sqrt2 facet counts (cross-checked with AlVin root lists)
BUGAENKO = {2: 3, 3: 5, 4: 7, 5: 10, 6: 34}


@pytest.mark.parametrize("n", sorted(MCLEOD_P3))
def test_mcleod_p3(n):
    r = vinberg(ZZ, QI(3), n, delta_max=200)
    assert r.status == "reflective"
    assert len(r.roots) == MCLEOD_P3[n]
    assert r.termination["ideal_vertices"] == MCLEOD_P3_IDEAL.get(n, 0)
    assert r.termination["finite_vertices"] + r.termination["ideal_vertices"] > 0


@pytest.mark.parametrize("n", sorted(MARK_P5))
def test_mark_p5(n):
    r = vinberg(ZZ, QI(5), n, delta_max=200)
    assert r.status == "reflective"
    assert len(r.roots) == MARK_P5[n]


def test_mark_p5_n9_does_not_terminate():
    r = vinberg(ZZ, QI(5), 9, delta_max=400)
    assert r.status == "undetermined"


@pytest.mark.parametrize("n", sorted(BUGAENKO))
def test_bugaenko(n):
    r = vinberg(ZS2, QI(1, 1), n, delta_max=2000)
    assert r.status == "reflective"
    assert len(r.roots) == BUGAENKO[n]
    assert r.termination["ideal_vertices"] == 0


def test_bugaenko_n7_does_not_terminate():
    r = vinberg(ZS2, QI(1, 1), 7, delta_max=3000, max_roots=120)
    assert r.status == "undetermined"


@pytest.mark.parametrize("d,nmax", [(7, 3), (11, 4), (13, 2), (17, 3), (19, 2), (23, 2)])
def test_thesis_table_3_1(d, nmax):
    """Mcleod thesis Table 3.1 (and Mark): f_d^n reflective exactly for 2 <= n <= nmax."""
    for n in range(2, nmax + 1):
        assert vinberg(ZZ, QI(d), n, delta_max=20000).status == "reflective", n
    assert vinberg(ZZ, QI(d), nmax + 1, delta_max=400, max_roots=80).status == "undetermined"
