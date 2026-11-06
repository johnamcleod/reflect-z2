#!/usr/bin/env python3
"""Standalone verifier for non-reflectivity certificates (data/certs/*.json).

Usage:  python3 verify_certificate.py data/certs/*.json [--facets]

Self-contained exact Z[sqrt2] arithmetic (independent of the search code).  For a
certificate with factors g_1 (, g_2) and g = g_1 (g_2) it checks:

  1. every entry of every g_i lies in Z[sqrt2];
  2. g_i^T F g_i = F with F = diag(-d, 1, ..., 1);
  3. g_i^{-1} (computed here) is integral and g_i g_i^{-1} = I;
  4. g_i preserves the future cone: (g_i v0, v0) < 0 at id;
  5. g_i maps the listed vertex roots S_from (in order) exactly onto S_to;
  6. every vector in S_from, S_to is a root of L (primitive, totally positive norm m,
     m | 2k_i, m | 2 d k_0) and S_from, S_to are elliptic of rank n (Gram positive
     definite, exact leading minors at id);
  7. g = product of the factors, and the characteristic polynomial of g recomputed
     here equals the stored one and changes sign on the stored interval (a, b) with
     a > 1 (exact signs at id): g has a real eigenvalue > 1, hence infinite order;
  8. (--facets) the roots of every S_to are facets of P: re-runs the deterministic
     Vinberg algorithm (reflect.vinberg) up to their distance and finds them.
     S_from is always the stabiliser simple system of v0, i.e. facets of P by step 1 of
     Vinberg's algorithm.

Given 1-8, g_i(P) = P for each factor (methods.md, Phase 3 step 4), so g is an
infinite-order element of Sym(P) and L is not reflective.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction

# ---------------- minimal exact Q(sqrt2) ---------------------------------------------


def parse(s):
    s = s.replace(" ", "")
    if "r2" not in s:
        return (Fraction(s), Fraction(0))
    idx = None
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1] not in "*/":
            idx = i
            break
    if idx is None:
        a, b = "0", s
    else:
        a, b = s[:idx], s[idx:]
    b = b.replace("*r2", "").replace("r2", "")
    b = {"": "1", "+": "1", "-": "-1"}.get(b, b)
    return (Fraction(a), Fraction(b))


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def mul(x, y): return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def neg(x): return (-x[0], -x[1])
def is_int(x): return x[0].denominator == 1 and x[1].denominator == 1
def is_zero(x): return x[0] == 0 and x[1] == 0


Z0, Z1 = (Fraction(0), Fraction(0)), (Fraction(1), Fraction(0))


def div(x, y):
    n = y[0] * y[0] - 2 * y[1] * y[1]
    num = mul(x, (y[0], -y[1]))
    return (num[0] / n, num[1] / n)


def sign(x):
    a, b = x
    if a >= 0 and b >= 0:
        return 0 if a == 0 and b == 0 else 1
    if a <= 0 and b <= 0:
        return -1
    t = a * a - 2 * b * b
    return (1 if t > 0 else -1) if a > 0 else (-1 if t > 0 else 1)


def conj(x): return (x[0], -x[1])
def tpos(x): return sign(x) > 0 and sign(conj(x)) > 0
def divides(m, x): return is_int(div(x, m))


def gcd(x, y):
    def nabs(z): return abs(z[0] * z[0] - 2 * z[1] * z[1])
    while not is_zero(y):
        q = div(x, y)
        q = (Fraction(round(q[0])), Fraction(round(q[1])))
        r = sub(x, mul(q, y))
        assert nabs(r) < nabs(y)
        x, y = y, r
    return x


def is_unit(x):
    return is_int(x) and abs(x[0] * x[0] - 2 * x[1] * x[1]) == 1


def matmul(A, B):
    return [[sum_(mul(A[i][k], B[k][j]) for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def sum_(it):
    s = Z0
    for x in it:
        s = add(s, x)
    return s


def transpose(A): return [list(r) for r in zip(*A)]
def eye(n): return [[Z1 if i == j else Z0 for j in range(n)] for i in range(n)]


def inverse(A):
    n = len(A)
    M = [list(A[i]) + [Z1 if i == j else Z0 for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if not is_zero(M[r][c]))
        M[c], M[p] = M[p], M[c]
        inv = div(Z1, M[c][c])
        M[c] = [mul(x, inv) for x in M[c]]
        for r in range(n):
            if r != c and not is_zero(M[r][c]):
                f = M[r][c]
                M[r] = [sub(x, mul(f, y)) for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def pd(G):
    k = len(G)
    A = [list(r) for r in G]
    for i in range(k):
        if sign(A[i][i]) <= 0:
            return False
        for r in range(i + 1, k):
            f = div(A[r][i], A[i][i])
            for c in range(i, k):
                A[r][c] = sub(A[r][c], mul(f, A[i][c]))
    return True


def charpoly(A):
    """det(xI - A) by Faddeev-LeVerrier; coefficients c_0..c_n."""
    n = len(A)
    c = [Z0] * (n + 1)
    c[n] = Z1
    M = [[Z0] * n for _ in range(n)]
    for k in range(1, n + 1):
        AM = matmul(A, M)
        M = [[add(AM[i][j], c[n - k + 1] if i == j else Z0) for j in range(n)] for i in range(n)]
        AMk = matmul(A, M)
        tr = sum_(AMk[i][i] for i in range(n))
        c[n - k] = mul(tr, (Fraction(-1, k), Fraction(0)))
    return c


def peval(c, x):
    r = Z0
    for a in reversed(c):
        r = add(mul(r, x), a)
    return r


# ---------------- the checks ---------------------------------------------------------


def verify(cert, check_facets=False):
    d = parse(cert["d"])
    n = cert["n"]
    N = n + 1
    F = [[(neg(d) if i == 0 else Z1) if i == j else Z0 for j in range(N)] for i in range(N)]

    def form(x, y):
        return sum_([mul(mul(neg(d), x[0]), y[0])] + [mul(a, b) for a, b in zip(x[1:], y[1:])])

    def is_root(k):
        if not all(is_int(x) for x in k):
            return False
        g = Z0
        for x in k:
            g = gcd(g, x) if not is_zero(g) else x
        if not is_unit(g):
            return False
        m = form(k, k)
        return (tpos(m) and all(divides(m, mul((Fraction(2), Fraction(0)), x)) for x in k[1:])
                and divides(m, mul(mul((Fraction(2), Fraction(0)), d), k[0])))

    gs = []
    for i, fac in enumerate(cert["factors"]):
        g = [[parse(x) for x in row] for row in fac["g"]]
        tag = f"factor {i}"
        assert all(is_int(x) for row in g for x in row), f"{tag}: non-integral entry"
        assert matmul(transpose(g), matmul(F, g)) == F, f"{tag}: does not preserve the form"
        gi = inverse(g)
        assert all(is_int(x) for row in gi for x in row), f"{tag}: inverse not integral"
        assert matmul(g, gi) == eye(N), f"{tag}: g g^-1 != I"
        assert sign(g[0][0]) > 0, f"{tag}: does not preserve the future cone"
        S_from = [[parse(x) for x in r] for r in fac["vertex_from"]]
        S_to = [[parse(x) for x in r] for r in fac["vertex_to"]]
        assert len(S_from) == n and len(S_to) == n
        for s, t in zip(S_from, S_to):
            gs_ = [sum_(mul(g[r][c], s[c]) for c in range(N)) for r in range(N)]
            assert gs_ == t, f"{tag}: g does not map the vertex roots as stated"
        # S_from must be the stabiliser simple roots v_{i+1}-v_i, -v_n
        stab = []
        for j in range(1, n):
            k = [Z0] * N
            k[j], k[j + 1] = neg(Z1), Z1
            stab.append(k)
        k = [Z0] * N
        k[n] = neg(Z1)
        stab.append(k)
        assert S_from == stab, f"{tag}: S_from is not the stabiliser simple system"
        for S in (S_from, S_to):
            assert all(is_root(s) for s in S), f"{tag}: vertex vector is not a root"
            G = [[form(a, b) for b in S] for a in S]
            assert pd(G), f"{tag}: vertex roots not elliptic of rank n"
            assert all(sign(G[a][b]) <= 0 for a in range(n) for b in range(n) if a != b), \
                f"{tag}: vertex roots not pairwise non-obtuse"
        gs.append(g)
    g = gs[0]
    for h in gs[1:]:
        g = matmul(g, h)
    assert g == [[parse(x) for x in row] for row in cert["g"]], "g != product of factors"
    chi = charpoly(g)
    assert chi == [parse(c) for c in cert["infinite_order"]["charpoly"]], "charpoly mismatch"
    a, b = (Fraction(x) for x in cert["infinite_order"]["bracket"])
    assert 1 < a < b
    sa, sb = sign(peval(chi, (a, Fraction(0)))), sign(peval(chi, (b, Fraction(0))))
    assert sa * sb < 0, "no sign change: eigenvalue > 1 not certified"
    if check_facets:
        _check_facets(cert)
    return True


def _check_facets(cert):
    import pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from reflect.ring import QI, ZS2
    from reflect.vinberg import vinberg
    d = QI.parse(cert["d"])
    n = cert["n"]
    need = {tuple(QI.parse(x) for x in r) for fac in cert["factors"] for r in fac["vertex_to"]}
    # distance of the farthest needed root: delta = d k0^2 / m
    def delta(k):
        m = -d * k[0] * k[0] + sum((x * x for x in k[1:]), QI(0))
        return (d * k[0] * k[0] / m).f()
    dmax = max(delta(k) for k in need) * (1 + 1e-9) + 1e-9
    res = vinberg(ZS2, d, n, delta_max=dmax, max_roots=10 ** 6, time_limit=10 ** 6)
    found = {r.k for r in res.roots}
    missing = need - found
    assert not missing, f"vertex roots not produced by Vinberg's algorithm: {missing}"


def main(argv):
    facets = "--facets" in argv
    paths = [a for a in argv if not a.startswith("--")]
    bad = 0
    for p in paths:
        cert = json.load(open(p))
        try:
            verify(cert, facets)
            print(f"OK    {p}  (d={cert['d']}, n={cert['n']})")
        except (AssertionError, ValueError, KeyError, StopIteration) as e:
            bad += 1
            print(f"FAIL  {p}: {e}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
