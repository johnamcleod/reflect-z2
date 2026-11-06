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


def test_guglielmetti_5_4r2_is_p7a():
    """Calibration (c): Guglielmetti's d = 5+4 sqrt2 equals p7a times eps^2 (same lattice);
    n = 5: 20 facets, f-vector (96, 240, 228, 102, 20, 1) via CoxIter; n = 6 not reflective."""
    import json
    from reflect.ring import QI
    from reflect.run import run_path
    eps = QI.parse("1+r2")
    assert QI.parse("-1+2*r2") * eps * eps == QI.parse("5+4*r2")
    r5 = json.loads(run_path("Zr2", QI.parse("-1+2*r2"), 5).read_text())
    assert r5["status"] == "reflective" and r5["n_roots"] == 20
    assert r5["coxiter_check"]["f_vector"] == [96, 240, 228, 102, 20, 1]
    assert r5["coxiter_check"]["cocompact"] == "yes"
    r6 = json.loads(run_path("Zr2", QI.parse("-1+2*r2"), 6).read_text())
    assert r6["status"] != "reflective"


@pytest.mark.parametrize("D", ["7+5*r2", "11+8*r2", "5+4*r2"])
def test_unit_square_rescaling_gives_same_chamber(D):
    """<-D,1,1,1> with D = d u^2 is isometric to <-d,1,1,1> via x0 -> u x0, so Vinberg
    must return the same chamber.  (AlVin at the pinned commit does NOT for D = 7+5r2 and
    11+8r2: it misses the root (r2-1, 1+r2, 0, 0); see OPEN_QUESTIONS OQ6.)"""
    import json
    from reflect.facet import same_line
    from reflect.linalg import sqrt_K
    from reflect.params import canonical
    from reflect.ring import QI, ZS2
    from reflect.run import run_path
    from reflect.vinberg import vinberg
    Dq = QI.parse(D)
    d = canonical(Dq)
    u = sqrt_K(Dq / d)
    u = u if u.sign() > 0 else -u
    r = vinberg(ZS2, Dq, 3, delta_max=1e5, max_roots=100, time_limit=300)
    ours = [tuple(QI.parse(x) for x in q["k"])
            for q in json.loads(run_path("Zr2", d, 3).read_text())["roots"]]
    mapped = [(u * k[0],) + tuple(k[1:]) for k in (x.k for x in r.roots)]
    assert r.status == "reflective" and len(mapped) == len(ours)
    assert all(any(same_line(m, o) for o in ours) for m in mapped)
