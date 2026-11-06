"""Phase 3: non-reflectivity certificates = infinite-order integral isometries g with
g(P) = P, where P is the fundamental chamber of the reflection subgroup W(L).

Construction (methods.md, Phase 3).  Let S0 = {e_1..e_n} be the stabiliser simple roots;
their mirrors meet at the vertex v0 of P (elliptic of rank n, type B_n).  For any other
vertex S' of P (an elliptic rank-n set of facets found by Vinberg) and any bijection
pi: S0 -> S' preserving the Gram matrix, there is at most one g in O'(V) with
g(e_i) = pi(e_i) and g(v0) a positive multiple of the vertex vector v' of S'.  If g is
integral, then g(P) = P: g maps the cone of P at v0 (cut out by S0) onto the cone of P
at v' (cut out by S'), so g(P) and P are W-chambers sharing interior points.
"""
from __future__ import annotations

import itertools
import json
import math
import pathlib
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .coxeter import Diagram, is_positive_definite
from .linalg import (charpoly, eye, from_columns, inverse, kernel, matmul, matvec,
                     poly_eval, sqrt_K, transpose)
from .ring import ONE, ZERO, QI, ZS2
from .vinberg import form

ROOT = pathlib.Path(__file__).resolve().parent.parent
CERTS = ROOT / "data" / "certs"


def gram_F(d: QI, n: int):
    return [[(-d if i == 0 else ONE) if i == j else ZERO for j in range(n + 1)]
            for i in range(n + 1)]


def vertex_vector(d: QI, S: Sequence[Sequence[QI]]) -> List[QI]:
    """Time-like vector orthogonal to the roots in S (rank n), future-pointing."""
    rows = [[(-d * s[0]) if i == 0 else s[i] for i in range(len(s))] for s in S]
    ker = kernel(rows)
    assert len(ker) == 1, "S does not have rank n"
    v = ker[0]
    if (-d * v[0]).sign() > 0:  # (v, v0) = -d v_0 must be < 0
        v = [-x for x in v]
    assert form(d, v, v).sign() < 0, "vertex vector is not time-like"
    return v


def euler_phi(k: int) -> int:
    r, m, p = k, k, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def finite_order_bound(dim: int) -> int:
    """Largest possible order of a finite-order element of GL_dim(O_K), [K:Q] = 2:
    eigenvalues are roots of unity of degree <= 2 dim over Q, so phi(order of each) <=
    2 dim, and the order is the lcm.  We return max k with phi(k) <= 2 dim (bounds each
    eigenvalue order); the element order is bounded by the lcm of such k, but g^k = I is
    tested for every k up to `lcm_bound` below."""
    return max(k for k in range(1, 20 * dim + 10) if euler_phi(k) <= 2 * dim)


def is_finite_order(g, dim: int) -> Optional[int]:
    """Exact: return the order of g if g^k = I for some k <= B, else None.  B = lcm of all
    k with phi(k) <= 2 dim is impractical; instead we use: g has finite order iff g^N = I
    with N = lcm{k : phi(k) <= 2 dim}, computed by repeated squaring."""
    ks = [k for k in range(1, 20 * dim + 10) if euler_phi(k) <= 2 * dim]
    N = 1
    for k in ks:
        N = N * k // math.gcd(N, k)
    P = _matpow(g, N)
    if P != eye(dim):
        return None
    # find the exact order: smallest divisor of N
    order = N
    for p in _prime_factors(N):
        while order % p == 0 and _matpow(g, order // p) == eye(dim):
            order //= p
    return order


def _prime_factors(N):
    out, p = [], 2
    while p * p <= N:
        if N % p == 0:
            out.append(p)
            while N % p == 0:
                N //= p
        p += 1
    if N > 1:
        out.append(N)
    return out


def _matpow(g, k):
    R = eye(len(g))
    B = g
    while k:
        if k & 1:
            R = matmul(R, B)
        B = matmul(B, B)
        k >>= 1
    return R


def loxodromic_bracket(g) -> Optional[Tuple[Fraction, Fraction, List[QI]]]:
    """Exact certificate of a real eigenvalue > 1 at id: rationals 1 < a < b with
    chi(a) chi(b) < 0 (signs exact in K at id).  Returns (a, b, chi) or None."""
    import numpy as np
    chi = charpoly(g)
    roots = np.roots([c.f() for c in reversed(chi)])
    for r in sorted(roots, key=lambda z: -abs(z)):
        if abs(r.imag) > 1e-9 * max(1, abs(r)) or r.real <= 1 + 1e-9:
            continue
        lam = r.real
        # simplest rationals 1 < a < lam < b with an exact sign change of chi at id
        for den in (2, 10, 100, 10 ** 3, 10 ** 4, 10 ** 6, 10 ** 9, 10 ** 12):
            a = Fraction(math.floor(lam * den), den)
            b = Fraction(math.floor(lam * den) + 1, den)
            if a <= 1:
                continue
            sa = poly_eval(chi, QI(a)).sign()
            sb = poly_eval(chi, QI(b)).sign()
            if sa * sb < 0:
                return a, b, chi
    return None


def _bijections(G0, Gs, norms0, norms1):
    """All bijections pi (as lists) with Gs[pi[i]][pi[j]] == G0[i][j]."""
    n = len(G0)
    out = []

    def rec(i, used, pi):
        if i == n:
            out.append(list(pi))
            return
        for c in range(n):
            if c in used or norms1[c] != norms0[i]:
                continue
            if all(Gs[c][pi[j]] == G0[i][j] for j in range(i)):
                pi.append(c); used.add(c)
                rec(i + 1, used, pi)
                pi.pop(); used.discard(c)

    rec(0, set(), [])
    return out


def isometry_from_vertices(d: QI, n: int, S0, v0, S1, v1) -> Optional[List[List[QI]]]:
    """The unique g in O'(V) with g(S0[i]) = S1[i] and g(v0) = lambda v1 (lambda > 0)."""
    lam2 = form(d, v0, v0) / form(d, v1, v1)
    lam = sqrt_K(lam2)
    if lam is None:
        return None
    M0 = from_columns(list(S0) + [v0])
    M1 = from_columns(list(S1) + [[lam * x for x in v1]])
    M0i = inverse(M0)
    return matmul(M1, M0i)


def check_isometry(d: QI, n: int, g) -> Dict:
    F = gram_F(d, n)
    integral = all(ZS2.contains(x) for row in g for x in row)
    preserves = matmul(transpose(g), matmul(F, g)) == F
    future = g[0][0].sign() > 0         # (g v0, v0) = -d g_00 < 0
    identity = g == eye(n + 1)
    return {"integral": integral, "preserves_form": preserves, "future": future,
            "identity": identity}


def search(d: QI, n: int, roots: List[Sequence[QI]], max_vertices: int = 100000,
           log: bool = False) -> Dict:
    """Search symmetries of P mapping the v0-vertex to other found vertices.

    Returns {"certificate": ... or None, "symmetries": [...finite-order ones...]}."""
    diag = Diagram(n, allow_parallel=False, track_parabolic=False)
    for i, k in enumerate(roots):
        row = [form(d, roots[j], k) for j in range(i)] + [form(d, k, k)]
        diag.add(row, row[-1])
    S0 = list(range(n))
    G0 = [[diag.Gq[i][j] for j in S0] for i in S0]
    norms0 = [diag.norms[i] for i in S0]
    v0 = [ONE] + [ZERO] * n
    finite_syms = []
    tried = 0
    for S in sorted(diag.ell_top, key=lambda s: sorted(s)):
        if S == frozenset(S0):
            continue
        Sl = sorted(S)
        norms1 = [diag.norms[i] for i in Sl]
        if sorted(norms1, key=lambda m: m.key()) != sorted(norms0, key=lambda m: m.key()):
            continue
        Gs = [[diag.Gq[i][j] for j in Sl] for i in Sl]
        v1 = None
        for pi in _bijections(G0, Gs, norms0, norms1):
            tried += 1
            if v1 is None:
                v1 = vertex_vector(d, [roots[i] for i in Sl])
            img = [Sl[c] for c in pi]
            g = isometry_from_vertices(d, n, [roots[i] for i in S0], v0,
                                       [roots[i] for i in img], v1)
            if g is None:
                continue
            chk = check_isometry(d, n, g)
            if not (chk["integral"] and chk["preserves_form"] and chk["future"]) or chk["identity"]:
                continue
            br = loxodromic_bracket(g)
            sym = {"g": g, "S0": S0, "S1": img, "v1": v1}
            if br is not None:
                return {"certificate": _make_cert(d, n, roots, [sym], br), "tried": tried}
            finite_syms.append(sym)
            if log:
                print("  finite-order symmetry ->", img)
        if tried > max_vertices:
            break
    # products of two finite-order symmetries
    for s1, s2 in itertools.combinations(finite_syms, 2):
        g = matmul(s1["g"], s2["g"])
        br = loxodromic_bracket(g)
        if br is not None:
            return {"certificate": _make_cert(d, n, roots, [s1, s2], br, product=True),
                    "tried": tried}
    return {"certificate": None, "tried": tried,
            "finite_symmetries": [[str(x) for row in s["g"] for x in row] for s in finite_syms]}


def _mat_json(g):
    return [[str(x) for x in row] for row in g]


def _make_cert(d, n, roots, syms, br, product=False) -> Dict:
    a, b, chi = br
    g = syms[0]["g"] if not product else matmul(syms[0]["g"], syms[1]["g"])
    factors = []
    for s in syms:
        factors.append({
            "g": _mat_json(s["g"]),
            "vertex_from": [[str(x) for x in roots[i]] for i in s["S0"]],
            "vertex_to": [[str(x) for x in roots[i]] for i in s["S1"]],
            "vertex_to_indices": s["S1"],
        })
    return {"d": str(d), "n": n, "ring": "Z[sqrt2]",
            "gram": "diag(-d,1,...,1)",
            "g": _mat_json(g),
            "factors": factors,
            "infinite_order": {"method": "real eigenvalue > 1 at id",
                               "charpoly": [str(c) for c in chi],
                               "bracket": [str(a), str(b)]},
            "facet_roots_needed": sorted({tuple(str(x) for x in roots[i])
                                          for s in syms for i in s["S1"]}),
            "n_roots_used": len(roots)}
