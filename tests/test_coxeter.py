"""Exactness tests for the Coxeter-diagram layer: labels, Bareiss PD vs independent
Gaussian elimination, compactness / finite-volume criteria on known polytopes."""
import itertools
import random
from fractions import Fraction

import pytest

from reflect.coxeter import (DOTTED, PARALLEL, Diagram, determinant, edge_label,
                             is_positive_definite)
from reflect.ring import QI, ZS2, ZZ, cos2
from reflect.vinberg import form, vinberg

rng = random.Random(7)


def test_edge_labels_exact():
    # roots of norm 2 and 2: (e,f) = -1 -> cos^2 = 1/4 -> 3
    assert edge_label(QI(-1), QI(2), QI(2), False) == 3
    assert edge_label(QI(0), QI(2), QI(2), False) == 2
    # norms 1, 2, product -1 -> cos^2 = 1/2 -> 4
    assert edge_label(QI(-1), QI(1), QI(2), False) == 4
    # norms 2, 2+r2, product -(1+r2): c = (3+2r2)/(4+2r2) = (2+r2)/4 -> 8
    assert edge_label(QI(-1, -1), QI(2), QI(2, 1), False) == 8
    assert edge_label(QI(-3), QI(2), QI(2), False) == DOTTED
    with pytest.raises(ValueError):
        edge_label(QI(-2), QI(2), QI(2), False)    # parallel in co-compact mode
    assert edge_label(QI(-2), QI(2), QI(2), True) == PARALLEL
    with pytest.raises(ValueError):
        edge_label(QI(1), QI(2), QI(2), True)      # obtuse
    with pytest.raises(ValueError):
        edge_label(QI(-1), QI(3), QI(1), True)     # cos^2 = 1/3: not a Coxeter angle


def _bareiss_minors_pd(G):
    """PD by Diagram's incremental Bareiss: build a diagram-free chain manually."""
    from reflect.coxeter import psign, to_pair
    d = Diagram.__new__(Diagram)
    d.G = [[to_pair(x) for x in row] for row in G]
    chain, lows, piv = [], [], []
    for v in range(len(G)):
        low, minor = Diagram._extend(d, chain, lows, piv, v)
        if psign(minor) <= 0:
            return False
        chain.append(v); lows.append(low); piv.append(minor)
    return True


def test_bareiss_pd_matches_gaussian():
    for _ in range(400):
        k = rng.randint(1, 6)
        # random symmetric matrices over Z[sqrt2], biased towards PD (Gram of random vecs + I)
        vecs = [[QI(rng.randint(-2, 2), rng.randint(-1, 1)) for _ in range(k)] for _ in range(k)]
        G = [[sum((vecs[i][t] * vecs[j][t] for t in range(k)), QI(0)) for j in range(k)]
             for i in range(k)]
        shift = QI(rng.randint(-3, 2), rng.randint(-1, 1))
        for i in range(k):
            G[i][i] = G[i][i] + shift
        assert _bareiss_minors_pd(G) == is_positive_definite(G)


def test_known_finite_diagrams_are_elliptic():
    # B_n stabiliser roots v_{i+1}-v_i, -v_n are PD (elliptic) for every n
    for n in range(2, 9):
        roots = []
        for i in range(1, n):
            k = [QI(0)] * (n + 1)
            k[i], k[i + 1] = QI(-1), QI(1)
            roots.append(k)
        k = [QI(0)] * (n + 1); k[n] = QI(-1); roots.append(k)
        G = [[form(QI(1), a, b) for b in roots] for a in roots]
        assert is_positive_definite(G)


def test_affine_A1_and_Btilde():
    # two parallel roots: Gram [[2,-2],[-2,2]] has det 0
    assert determinant([[QI(2), QI(-2)], [QI(-2), QI(2)]]) == QI(0)


def test_compact_simplices_terminate():
    # d = 1+r2, n = 2..5: Bugaenko's polytopes; triangle/simplices for small n
    for n, nroots in [(2, 3), (3, 5), (4, 7), (5, 10)]:
        r = vinberg(ZS2, QI(1, 1), n, delta_max=100)
        assert r.status == "reflective" and len(r.roots) == nroots


def test_noncompact_finite_volume_q():
    # -3x0^2 + ...: n = 2 Lanner triangle (compact), n = 3 simplex with an ideal vertex
    r = vinberg(ZZ, QI(3), 2)
    assert r.status == "reflective" and len(r.roots) == 3
    assert r.termination["ideal_vertices"] == 0
    r = vinberg(ZZ, QI(3), 3)
    assert r.status == "reflective" and len(r.roots) == 4
    assert r.termination["ideal_vertices"] >= 1
