"""Vinberg's algorithm for an arbitrary rank-3 lattice over Z[sqrt2] with Gram matrix G
(signature (2,1) at id, positive definite at sigma), e.g. an entry of Mark's table.

Roots: primitive r with m = Q(r) totally positive and 2 G r in m O^3.  Then 2r/m lies in
the dual L^#, and delta L^# <= L (delta = det G), so m | 2 delta: the norms are the totally
positive divisors of 2 delta modulo unit squares (one representative m per class; each
mirror then carries exactly one root +-r of norm exactly m).

Control vector v0 in L, timelike (Q(v0) < 0 at id), chosen by the caller; if v0 is in the
interior of a chamber no root is orthogonal to it.  Candidates with distance
D(r) = (r,v0)^2 / (m |Q(v0)|) <= Delta satisfy, at id, H(r) = Q(r) - 2 (r,v0)^2/Q(v0)
= m (1 + 2 D) <= m (1 + 2 Delta) (H is the positive majorant w.r.t. v0) and, at sigma,
sigma Q(r) = sigma(m).  So they are lattice points of Z^6 (coordinates a_j + b_j sqrt2)
in the ellipsoid H_id + sigma Q_sigma <= m_id (1 + 2 Delta) + m_sigma.  We enumerate
it by Fincke-Pohst on an exactly LLL-reduced basis (rational arithmetic with a 40-digit
rational sqrt2; floats only inside the enumeration, with a safety margin), and every
candidate is checked exactly.  Floats never accept or reject a root.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .coxeter import Diagram
from .params import root_norms
from .ring import QI, ZS2

S2 = Fraction(math.isqrt(2 * 10 ** 80), 10 ** 40)   # sqrt2 to ~40 digits (rational)


def ipG(G, x, y) -> QI:
    return sum((x[i] * G[i][j] * y[j] for i in range(3) for j in range(3)), QI(0))


def _emb(x: QI, s) -> Fraction:
    return Fraction(x.a) + Fraction(x.b) * s


def _six_grams(G, v0):
    """Grams (A_id, A_sigma) on Z^6 (y = (a0,b0,a1,b1,a2,b2), x_j = a_j + b_j sqrt2) of
    H_id(x) and Q_sigma(x), with sqrt2 replaced by the rational S2."""
    Q0 = ipG(G, v0, v0)
    Gv = [sum((G[i][j] * v0[j] for j in range(3)), QI(0)) for i in range(3)]
    Gid = [[_emb(G[i][j], S2) for j in range(3)] for i in range(3)]
    Gsg = [[_emb(G[i][j], -S2) for j in range(3)] for i in range(3)]
    gv = [_emb(x, S2) for x in Gv]
    q0 = _emb(Q0, S2)
    Hid = [[Gid[i][j] - 2 * gv[i] * gv[j] / q0 for j in range(3)] for i in range(3)]
    # E_id: y -> x_id ; column for (a_j) is e_j, for (b_j) is S2 e_j ; sigma: -S2
    def col(k, s):
        j, t = divmod(k, 2)
        v = [Fraction(0)] * 3
        v[j] = Fraction(1) if t == 0 else s
        return v
    Ai = [[Fraction(0)] * 6 for _ in range(6)]
    As = [[Fraction(0)] * 6 for _ in range(6)]
    for p in range(6):
        for q in range(6):
            ci, cj = col(p, S2), col(q, S2)
            si, sj = col(p, -S2), col(q, -S2)
            Ai[p][q] = sum(ci[i] * Hid[i][j] * cj[j] for i in range(3) for j in range(3))
            As[p][q] = sum(si[i] * Gsg[i][j] * sj[j] for i in range(3) for j in range(3))
    return Ai, As


def _lll(A: List[List[Fraction]], delta=Fraction(3, 4)):
    """Exact LLL on a positive definite Gram matrix A (Fractions).  Returns (U, A') with
    A' = U^T A U, U unimodular integer (columns = new basis in old coordinates)."""
    n = len(A)
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    B = [row[:] for row in A]

    def gso():
        mu = [[Fraction(0)] * n for _ in range(n)]
        Bs = [Fraction(0)] * n
        for i in range(n):
            for j in range(i):
                s = B[i][j] - sum(mu[j][k] * mu[i][k] * Bs[k] for k in range(j))
                mu[i][j] = s / Bs[j]
            Bs[i] = B[i][i] - sum(mu[i][k] ** 2 * Bs[k] for k in range(i))
        return mu, Bs

    def add(i, j, c):  # b_i += c b_j
        for r in range(n):
            U[r][i] += c * U[r][j]
        for r in range(n):
            B[i][r] += c * B[j][r]
        for r in range(n):
            B[r][i] = B[i][r]
        B[i][i] = B[i][i]  # already updated via row; recompute diagonal exactly below

    def recompute():
        # B = U^T A U exactly
        for i in range(n):
            for j in range(n):
                B[i][j] = sum(U[p][i] * A[p][q] * U[q][j] for p in range(n) for q in range(n))

    k = 1
    recompute()
    while k < n:
        mu, Bs = gso()
        for j in range(k - 1, -1, -1):
            c = round(mu[k][j])
            if c:
                for r in range(n):
                    U[r][k] -= c * U[r][j]
                recompute()
                mu, Bs = gso()
        if Bs[k] >= (delta - mu[k][k - 1] ** 2) * Bs[k - 1]:
            k += 1
        else:
            for r in range(n):
                U[r][k], U[r][k - 1] = U[r][k - 1], U[r][k]
            recompute()
            k = max(k - 1, 1)
    return U, B


def _fincke_pohst(Ar: np.ndarray, C: float):
    """All integer z with z^T Ar z <= C (floats; caller re-checks exactly)."""
    n = Ar.shape[0]
    L = np.linalg.cholesky(Ar)            # Ar = L L^T
    R = L.T                               # upper triangular, Ar = R^T R
    q = np.zeros((n, n))
    for i in range(n):
        q[i][i] = R[i][i] ** 2
        for j in range(i + 1, n):
            q[i][j] = R[i][j] / R[i][i]
    z = [0] * n
    out = []

    def rec(i, rem):
        c = -sum(q[i][j] * z[j] for j in range(i + 1, n))
        r = math.sqrt(max(rem, 0) / q[i][i])
        lo, hi = math.ceil(c - r - 1e-9), math.floor(c + r + 1e-9)
        for t in range(lo, hi + 1):
            z[i] = t
            used = q[i][i] * (t - c) ** 2
            if used > rem + 1e-7 * (1 + C):
                continue
            if i == 0:
                out.append(tuple(z))
            else:
                rec(i - 1, rem - used)
        z[i] = 0

    rec(n - 1, C * (1 + 1e-9) + 1e-9)
    return out


class Lattice3:
    def __init__(self, G, v0):
        self.G = G
        self.v0 = list(v0)
        self.Q0 = ipG(G, v0, v0)
        assert self.Q0.sign() < 0, "v0 must be timelike at id"
        from .coxeter import determinant
        self.delta = determinant(G)
        self.norms = root_norms(-self.delta)    # tot. pos. divisors of 2 delta mod unit^2
        self.Ai, self.As = _six_grams(G, v0)
        self._Gid = np.array([[G[i][j].f() for j in range(3)] for i in range(3)])
        self._Gsg = np.array([[G[i][j].fs() for j in range(3)] for i in range(3)])

    def _reduced(self, Ci: float, Cs: float):
        """LLL basis for H_id/Ci + Q_sigma/Cs (region H_id <= Ci, Q_sigma <= Cs lies in
        {form <= 2})."""
        ci, cs = Fraction(Ci).limit_denominator(10 ** 12), Fraction(Cs).limit_denominator(10 ** 12)
        A = [[self.Ai[p][q] / ci + self.As[p][q] / cs for q in range(6)] for p in range(6)]
        U, Ared = _lll(A)
        return U, np.array([[float(x) for x in row] for row in Ared])

    def _vec(self, U, z) -> List[QI]:
        y = [sum(U[p][k] * z[k] for k in range(6)) for p in range(6)]
        return [QI(y[2 * j], y[2 * j + 1]) for j in range(3)]

    def is_root(self, r, m) -> bool:
        if ipG(self.G, r, r) != m or not ZS2.is_unit(ZS2.gcd_list(r)):
            return False
        Gr = [sum((self.G[i][j] * r[j] for j in range(3)), QI(0)) for i in range(3)]
        return all(ZS2.divides(m, 2 * x) for x in Gr)

    def roots_upto(self, Delta: float) -> List[Tuple[QI, QI, Tuple[QI, ...]]]:
        """(exact distance, m, r) for all roots with (r,v0) < 0 and distance <= Delta."""
        out = []
        for m in self.norms:
            Ci = m.f() * (1 + 2 * Delta) * (1 + 1e-9)
            Cs = m.fs() * (1 + 1e-9)
            U, Ar = self._reduced(Ci, Cs)
            pts = _fincke_pohst(Ar, 2.0)
            if not pts:
                continue
            # float prefilter Q ~ m at both embeddings (generous tolerance scaled by the
            # absolute-value form); survivors are re-checked exactly below
            Z = np.array(pts, dtype=float)
            Y = Z @ np.array(U, dtype=float).T
            r2 = math.sqrt(2)
            keep = np.ones(len(pts), dtype=bool)
            for s, Gf, mm in ((r2, self._Gid, m.f()), (-r2, self._Gsg, m.fs())):
                X = Y[:, 0::2] + s * Y[:, 1::2]
                q = np.einsum("ni,ij,nj->n", X, Gf, X)
                qa = np.einsum("ni,ij,nj->n", np.abs(X), np.abs(Gf), np.abs(X))
                keep &= np.abs(q - mm) <= 1e-6 * (1 + qa)
            for idx in np.nonzero(keep)[0]:
                z = pts[idx]
                if not any(z):
                    continue
                r = self._vec(U, z)
                if ipG(self.G, r, r) != m:
                    continue
                t = ipG(self.G, r, self.v0)
                if t.sign() >= 0:     # keep the root with (r, v0) < 0 (the other is -r)
                    continue
                if not self.is_root(r, m):
                    continue
                dist = t * t / (m * (-self.Q0))
                if dist.f() <= Delta * (1 + 1e-12):
                    out.append((dist, m, tuple(r)))
        return out


def vinberg3(G, v0, Delta0: float = 4.0, Delta_max: float = 1e5,
             max_roots: int = 200) -> Dict:
    """Vinberg's algorithm (v0 must lie in the interior of a chamber: asserted, since a
    root orthogonal to v0 would be found with distance 0)."""
    L = Lattice3(G, v0)
    acc: List[Tuple[QI, ...]] = []
    diag = Diagram(2, allow_parallel=False, track_parabolic=False)
    done_upto = -1.0
    Delta = Delta0
    status = "undetermined"
    while Delta <= Delta_max and len(acc) < max_roots:
        cands = [c for c in L.roots_upto(Delta) if c[0].f() > done_upto]
        cands.sort(key=lambda c: c[0].f())
        # batches of exactly equal distance
        i = 0
        while i < len(cands):
            j = i
            while j < len(cands) and cands[j][0] == cands[i][0]:
                j += 1
            batch = cands[i:j]
            assert not batch[0][0].is_zero(), "v0 lies on a mirror"
            surv = [r for _, m, r in batch
                    if all(ipG(G, a, r).sign() <= 0 for a in acc)]
            for a in range(len(surv)):
                for b in range(a + 1, len(surv)):
                    assert ipG(G, surv[a], surv[b]).sign() <= 0, "same-batch obtuse pair"
            for r in surv:
                row = [ipG(G, a, r) for a in acc] + [ipG(G, r, r)]
                diag.add(row, row[-1])
                acc.append(r)
            if surv and len(acc) >= 3 and diag.check_termination()["terminated"]:
                return {"status": "reflective", "roots": acc, "Delta": Delta}
            i = j
        done_upto = Delta
        Delta *= 2
    return {"status": status, "roots": acc, "Delta": Delta / 2}
