"""Checks behind the n -> n+1 lemma (OPEN_QUESTIONS.md OQ1): at a norm-1 root of an
integral Z[sqrt2]-lattice, no Coxeter angle pi/3, pi/6, pi/8 can occur; and every
certified dimension range is upward-closed (consistency with the lemma)."""
import json
import pathlib

from reflect.ring import EPS, QI, ZS2, cos2

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_no_odd_angles_at_norm_one_root():
    # cos^2 = t^2/m with t in O, m totally positive, 2t/m in O.  Write t = m u / 2.
    # For k in {3, 6, 8} we need m u^2 = 4 cos^2(pi/k) and then t = m u / 2 = 2cos^2/u.
    for k in (3, 6, 8):
        c4 = cos2(k) * 4               # 1, 3, 2 + sqrt2
        for j in range(-6, 7):         # u ranges over all units up to sign and squares
            for s in (1, -1):
                u = s * EPS ** j
                t = c4 / (2 * u)
                assert not ZS2.contains(t), (k, u, t)
        # u non-unit: m = c4 / u^2 must be integral, so u^2 | c4; c4 has no square factors
        # except units (1, 3, 2 + sqrt2 = sqrt2 * EPS are squarefree), so u is a unit.
        for c in (QI(1), QI(3), QI(2, 1)):
            _, fs = ZS2.factor(c)
            assert all(e == 1 for _, e in fs)


def test_certified_ranges_are_upward_closed():
    certs = {p.name for p in (ROOT / "data/certs").glob("*.json")}
    by_label = {}
    for name in certs:
        label, n = name[:-5].rsplit("_n", 1)
        by_label.setdefault(label, set()).add(int(n))
    for label, ns in by_label.items():
        n0 = min(ns)
        assert ns == set(range(n0, 11)), (label, sorted(ns))
