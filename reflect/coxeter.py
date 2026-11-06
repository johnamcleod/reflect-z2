"""Coxeter diagrams of root systems in a hyperbolic O_K-lattice, decided exactly.

Ellipticity of a set S of roots is decided by exact positive-definiteness of the raw
Gram submatrix G_S = ((e_i, e_j)) (entries in O_K).  G_S is congruent to the normalised
Gram matrix via diag(m_i^(-1/2)), so it has the same inertia.  Leading principal minors
are computed fraction-free (Bareiss) in O_K, incrementally along a DFS, and their signs
are taken exactly at id.

Elements of O_K in the hot loops are plain int pairs (a, b) = a + b sqrt2.
"""
from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from .ring import QI, _sgn, cos2

Pair = Tuple[int, int]


# --- fast O_K arithmetic on int pairs ------------------------------------------------
def pmul(x: Pair, y: Pair) -> Pair:
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def psub(x: Pair, y: Pair) -> Pair:
    return (x[0] - y[0], x[1] - y[1])


def pdiv_exact(x: Pair, y: Pair) -> Pair:
    """x / y, asserting the quotient lies in O_K (Bareiss divisions are exact)."""
    n = y[0] * y[0] - 2 * y[1] * y[1]
    a = x[0] * y[0] - 2 * x[1] * y[1]
    b = -x[0] * y[1] + x[1] * y[0]
    qa, ra = divmod(a, n)
    qb, rb = divmod(b, n)
    if ra or rb:
        raise ArithmeticError(f"inexact Bareiss division {x}/{y}")
    return (qa, qb)


def psign(x: Pair) -> int:
    return _sgn(x[0], x[1])


def to_pair(x: QI) -> Pair:
    assert x.is_integral(), x
    return (int(x.a), int(x.b))


# --- labels ---------------------------------------------------------------------------
DOTTED = "dotted"      # ultraparallel mirrors
PARALLEL = "inf"       # parallel mirrors (meet at infinity)


def edge_label(g: QI, mi: QI, mj: QI, allow_parallel: bool) -> object:
    """Coxeter label of a pair of roots with (e_i, e_j) = g <= 0 and norms mi, mj.
    Returns 2 (orthogonal), 3, 4, 6, 8, PARALLEL or DOTTED.  Raises on anything else."""
    if g.sign() > 0:
        raise ValueError(f"obtuse pair: (e_i,e_j) = {g} > 0")
    c = (g * g) / (mi * mj)
    for k in (2, 3, 4, 6, 8):
        if c == cos2(k):
            return k
    if c == QI(1):
        if not allow_parallel:
            raise ValueError("parallel mirrors (cos^2 = 1) cannot occur for a co-compact group")
        return PARALLEL
    if c.sign() > 0 and (c - 1).sign() > 0:
        return DOTTED
    raise ValueError(f"non-Coxeter angle: cos^2 = {c}")


class Diagram:
    """Roots with exact Gram matrix and labels; supports incremental enumeration of
    elliptic subdiagrams (and, for the non-compact case, affine components)."""

    def __init__(self, n: int, allow_parallel: bool, track_parabolic: bool):
        self.n = n                      # hyperbolic dimension; rank of a vertex = n
        self.allow_parallel = allow_parallel
        self.track_parabolic = track_parabolic
        self.G: List[List[Pair]] = []   # raw Gram matrix, int pairs
        self.Gq: List[List[QI]] = []
        self.norms: List[QI] = []
        self.labels: Dict[Tuple[int, int], object] = {}
        self.compatible: List[Set[int]] = []  # j with label not dotted
        self.adjacent: List[Set[int]] = []    # j with label != 2
        # enumeration results
        self.ell_count = defaultdict(int)     # size -> number of elliptic sets
        self.ell_top: Set[FrozenSet[int]] = set()      # elliptic sets of size n
        self.ell_sub: Dict[FrozenSet[int], int] = {}   # size n-1 -> #containing vertices
        self.affine: Set[FrozenSet[int]] = set()       # connected affine subdiagrams

    def __len__(self):
        return len(self.G)

    def add(self, gram_row: Sequence[QI], norm: QI):
        """Add a root given its products with all existing roots and itself."""
        i = len(self.G)
        assert len(gram_row) == i + 1 and gram_row[i] == norm
        row = [to_pair(x) for x in gram_row]
        for j in range(i):
            self.G[j].append(row[j])
            self.Gq[j].append(gram_row[j])
        self.G.append(row)
        self.Gq.append(list(gram_row))
        self.norms.append(norm)
        self.compatible.append(set())
        self.adjacent.append(set())
        for j in range(i):
            lab = edge_label(gram_row[j], self.norms[j], norm, self.allow_parallel)
            self.labels[(j, i)] = self.labels[(i, j)] = lab
            # 'compatible' = not ultraparallel.  Parallel pairs (non-compact case only)
            # give a zero 2x2 minor, so the DFS records them as affine A~1 and never
            # extends through them.
            if lab != DOTTED:
                self.compatible[i].add(j)
                self.compatible[j].add(i)
            if lab != 2:
                self.adjacent[i].add(j)
                self.adjacent[j].add(i)
        self._enumerate_containing(i)

    # ---- incremental Bareiss along a chain ----------------------------------------
    def _extend(self, chain: List[int], lows: List[List[Pair]], piv: List[Pair], v: int
                ) -> Tuple[List[Pair], Pair]:
        """Given chain s_1..s_k with stored lows[t][j] = M^(j)_{s_t, s_{j+1}} (j <= t) and
        pivots piv[j] = leading minor j+1 (piv[-1]... ), return (low_v, minor_{k+1})."""
        k = len(chain)
        G = self.G
        cur = [G[v][s] for s in chain]
        cvv = G[v][v]
        low_v = []
        prev = (1, 0)
        for j in range(k):
            a = cur[j]
            low_v.append(a)
            pj = piv[j]
            for t in range(j + 1, k):
                cur[t] = pdiv_exact(psub(pmul(pj, cur[t]), pmul(a, lows[t][j])), prev)
            cvv = pdiv_exact(psub(pmul(pj, cvv), pmul(a, a)), prev)
            prev = pj
        low_v.append(cvv)
        return low_v, cvv

    def _enumerate_containing(self, v: int):
        """Enumerate all elliptic subsets (size <= n) containing v and otherwise older
        vertices; record sizes n and n-1 and (optionally) affine components."""
        n = self.n
        comp = self.compatible
        low_v, minor = self._extend([], [], [], v)
        if psign(minor) <= 0:
            raise ValueError("root of non-positive norm")
        chain, lows, piv = [v], [low_v], [minor]
        self._record(chain)
        cands = sorted(comp[v])
        self._dfs(chain, lows, piv, cands)

    def _dfs(self, chain, lows, piv, cands):
        n = self.n
        for idx, w in enumerate(cands):
            low_w, minor = self._extend(chain, lows, piv, w)
            s = psign(minor)
            if s > 0:
                if len(chain) + 1 <= n:
                    chain.append(w)
                    lows.append(low_w)
                    piv.append(minor)
                    self._record(chain)
                    if len(chain) < n:
                        nxt = [u for u in cands[idx + 1:] if u in self.compatible[w]]
                        self._dfs(chain, lows, piv, nxt)
                    chain.pop(); lows.pop(); piv.pop()
            elif s == 0 and self.track_parabolic:
                self._record_affine(chain, w)

    def _record(self, chain):
        k = len(chain)
        self.ell_count[k] += 1
        if k == self.n:
            S = frozenset(chain)
            if S in self.ell_top:
                return
            self.ell_top.add(S)
            for x in S:
                sub = S - {x}
                self.ell_sub[sub] = self.ell_sub.get(sub, 0) + 1
        elif k == self.n - 1:
            S = frozenset(chain)
            self.ell_sub.setdefault(S, 0)

    def _record_affine(self, chain, w):
        T = set(chain) | {w}
        # connected component of w inside T
        comp, stack = {w}, [w]
        while stack:
            x = stack.pop()
            for y in self.adjacent[x]:
                if y in T and y not in comp:
                    comp.add(y)
                    stack.append(y)
        self.affine.add(frozenset(comp))

    # ---- criteria -----------------------------------------------------------------
    def parabolic_sets(self) -> List[Tuple[FrozenSet[int], ...]]:
        """Parabolic subdiagrams of rank n-1: unions of mutually orthogonal affine
        components with total rank sum(|C|-1) = n-1."""
        comps = sorted(self.affine, key=lambda c: sorted(c))
        target = self.n - 1
        out = []

        def orth(c1, c2):
            return all(self.labels[(x, y)] == 2 for x in c1 for y in c2)

        def rec(start, chosen, rank):
            if rank == target:
                out.append(tuple(chosen))
                return
            for i in range(start, len(comps)):
                c = comps[i]
                r = len(c) - 1
                if rank + r > target:
                    continue
                if any(c & x for x in chosen):
                    continue
                if all(orth(c, x) for x in chosen):
                    chosen.append(c)
                    rec(i + 1, chosen, rank + r)
                    chosen.pop()

        rec(0, [], 0)
        return out

    def edge_vertex_counts(self) -> Tuple[Dict[FrozenSet[int], int], int, int]:
        """For each elliptic (n-1)-set (an edge of P), the number of vertices on it
        (finite: elliptic n-sets; ideal: parabolic rank n-1 sets).  Returns
        (counts, #finite vertices, #ideal vertices)."""
        counts = dict(self.ell_sub)
        n_ideal = 0
        if self.track_parabolic:
            for P in self.parabolic_sets():
                n_ideal += 1
                for pick in itertools.product(*[sorted(c) for c in P]):
                    S = frozenset().union(*P) - set(pick)
                    if S in counts:
                        counts[S] += 1
                    else:
                        # S must be elliptic (remove one vertex from each component)
                        raise AssertionError("edge of an ideal vertex not enumerated")
        return counts, len(self.ell_top), n_ideal

    def check_termination(self) -> Dict:
        counts, nfin, nideal = self.edge_vertex_counts()
        bad = [S for S, c in counts.items() if c != 2]
        over = [S for S, c in counts.items() if c > 2]
        if over:
            raise AssertionError(f"edge with {counts[over[0]]} vertices: {sorted(over[0])}")
        ok = (nfin + nideal) > 0 and not bad
        return {"terminated": ok, "finite_vertices": nfin, "ideal_vertices": nideal,
                "edges": len(counts), "open_edges": len(bad)}

    # ---- export ---------------------------------------------------------------------
    def coxeter_matrix_rows(self) -> List[Tuple[int, int, object]]:
        out = []
        for i in range(len(self.G)):
            for j in range(i + 1, len(self.G)):
                lab = self.labels[(i, j)]
                if lab != 2:
                    out.append((i, j, lab))
        return out

    def to_coxiter(self) -> str:
        """CoxIter input: '<#vertices> <dim>' then 'i j w' (1-based), w = 0 for an
        ultraparallel (dotted) pair and 1 for parallel (infinite) -- see
        tests/test_coxiter_bridge.py for the check of this convention."""
        lines = [f"{len(self.G)} {self.n}"]
        for i, j, lab in self.coxeter_matrix_rows():
            w = {DOTTED: 1, PARALLEL: 0}.get(lab, lab)
            lines.append(f"{i + 1} {j + 1} {w}")
        return "\n".join(lines) + "\n"


def is_positive_definite(G: Sequence[Sequence[QI]]) -> bool:
    """Exact PD test of a symmetric matrix over K by leading principal minors at id
    (Gaussian elimination in K).  Independent of the Bareiss code above."""
    k = len(G)
    A = [[QI.of(x) for x in row] for row in G]
    for i in range(k):
        if A[i][i].sign() <= 0:
            return False
        for r in range(i + 1, k):
            f = A[r][i] / A[i][i]
            for c in range(i, k):
                A[r][c] = A[r][c] - f * A[i][c]
    return True


def determinant(G: Sequence[Sequence[QI]]) -> QI:
    k = len(G)
    A = [[QI.of(x) for x in row] for row in G]
    det = QI(1)
    for i in range(k):
        p = next((r for r in range(i, k) if not A[r][i].is_zero()), None)
        if p is None:
            return QI(0)
        if p != i:
            A[i], A[p] = A[p], A[i]
            det = -det
        det = det * A[i][i]
        for r in range(i + 1, k):
            f = A[r][i] / A[i][i]
            for c in range(i, k):
                A[r][c] = A[r][c] - f * A[i][c]
    return det
