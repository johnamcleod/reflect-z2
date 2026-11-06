"""Exact linear algebra over K = Q(sqrt2) with QI entries."""
from __future__ import annotations

from fractions import Fraction
from typing import List, Optional, Sequence

from .ring import ONE, ZERO, QI

Mat = List[List[QI]]


def mat(rows) -> Mat:
    return [[QI.of(x) for x in r] for r in rows]


def eye(n: int) -> Mat:
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def matmul(A: Mat, B: Mat) -> Mat:
    n, m, p = len(A), len(B), len(B[0])
    out = []
    for i in range(n):
        row = []
        Ai = A[i]
        for j in range(p):
            s = ZERO
            for k in range(m):
                if not Ai[k].is_zero() and not B[k][j].is_zero():
                    s = s + Ai[k] * B[k][j]
            row.append(s)
        out.append(row)
    return out


def matvec(A: Mat, v: Sequence[QI]) -> List[QI]:
    return [sum((a * b for a, b in zip(row, v)), ZERO) for row in A]


def transpose(A: Mat) -> Mat:
    return [list(r) for r in zip(*A)]


def from_columns(cols: Sequence[Sequence[QI]]) -> Mat:
    return transpose([list(c) for c in cols])


def inverse(A: Mat) -> Optional[Mat]:
    n = len(A)
    M = [list(A[i]) + [ONE if i == j else ZERO for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if not M[r][c].is_zero()), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        inv = ONE / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and not M[r][c].is_zero():
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def kernel(A: Mat) -> List[List[QI]]:
    """Basis of the right kernel of A (rows x cols) over K."""
    rows, cols = len(A), len(A[0])
    M = [list(r) for r in A]
    piv_cols = []
    r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if not M[i][c].is_zero()), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = ONE / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(rows):
            if i != r and not M[i][c].is_zero():
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in piv_cols]
    basis = []
    for f in free:
        v = [ZERO] * cols
        v[f] = ONE
        for i, pc in enumerate(piv_cols):
            v[pc] = -M[i][f]
        basis.append(v)
    return basis


def rank(A: Mat) -> int:
    return len(A[0]) - len(kernel(A))


def charpoly(A: Mat) -> List[QI]:
    """Coefficients c_0..c_n of det(xI - A) = sum c_k x^k (Faddeev-LeVerrier, exact)."""
    n = len(A)
    coeffs = [ZERO] * (n + 1)
    coeffs[n] = ONE
    M = [[ZERO] * n for _ in range(n)]
    I = eye(n)
    c = ONE
    for k in range(1, n + 1):
        # M_k = A M_{k-1} + c_{n-k+1} I
        AM = matmul(A, M) if k > 1 else [[ZERO] * n for _ in range(n)]
        M = [[AM[i][j] + (coeffs[n - k + 1] if i == j else ZERO) for j in range(n)]
             for i in range(n)]
        AMk = matmul(A, M)
        tr = sum((AMk[i][i] for i in range(n)), ZERO)
        coeffs[n - k] = tr * QI(Fraction(-1, k))
    return coeffs


def poly_eval(coeffs: Sequence[QI], x: QI) -> QI:
    r = ZERO
    for c in reversed(coeffs):
        r = r * x + c
    return r


def sqrt_K(x: QI) -> Optional[QI]:
    """Exact square root in K if it exists (positive at id), else None."""
    x = QI.of(x)
    if x.is_zero():
        return ZERO
    a, b = Fraction(x.a), Fraction(x.b)
    N = a * a - 2 * b * b
    if N < 0:
        return None
    rn = _frac_sqrt(N)
    if rn is None:
        return None
    for s in (rn, -rn):
        p2 = (a + s) / 2
        if p2 < 0:
            continue
        p = _frac_sqrt(p2)
        if p is None:
            continue
        if p == 0:
            q2 = a / 2
            q = _frac_sqrt(q2) if q2 >= 0 else None
            if q is None:
                continue
            cand = QI(0, q)
        else:
            cand = QI(p, b / (2 * p))
        for c in (cand, -cand):
            if c * c == x and c.sign() > 0:
                return c
    return None


def _frac_sqrt(q: Fraction) -> Optional[Fraction]:
    import math
    q = Fraction(q)
    if q < 0:
        return None
    n, d = q.numerator, q.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn * rn == n and rd * rd == d:
        return Fraction(rn, rd)
    return None
