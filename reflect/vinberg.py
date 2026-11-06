"""Vinberg's algorithm (Bugaenko's version for totally real rings) for the lattice
L = O^(n+1) with Gram matrix diag(-d, 1, ..., 1), control vector v0.

Generic over the Ring interface (ZZ for the calibration over Q, ZS2 for Z[sqrt2]).
All acceptance / termination decisions are exact; floats only order the batches.
"""
from __future__ import annotations

import itertools
import json
import math
import time
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

from .coxeter import Diagram
from .params import root_norms as zs2_root_norms
from .ring import ONE, ZERO, QI, Ring, ZS2, ZZ

Vec = Tuple[QI, ...]


def form(d: QI, x: Sequence[QI], y: Sequence[QI]) -> QI:
    s = -d * x[0] * y[0]
    for a, b in zip(x[1:], y[1:]):
        s = s + a * b
    return s


def root_norms(ring: Ring, d: QI) -> List[QI]:
    if ring is ZS2:
        return zs2_root_norms(d)
    D = int(2 * d.a)
    return [QI(m) for m in range(1, D + 1) if D % m == 0]


def is_root(ring: Ring, d: QI, k: Sequence[QI]) -> bool:
    """Exact root test: primitive, totally positive norm m, m | 2 k_i (i>=1), m | 2 d k_0."""
    if not all(ring.contains(x) for x in k):
        return False
    if not ring.is_unit(ring.gcd_list(k)):
        return False
    m = form(d, k, k)
    if not ring.is_totally_positive(m):
        return False
    return all(ring.divides(m, 2 * x) for x in k[1:]) and ring.divides(m, 2 * d * k[0])


class SumOfSquares:
    """Solutions of y_1^2 + ... + y_s^2 = R in the ring with y_1 >= ... >= y_s >= 0 at id."""

    def __init__(self, ring: Ring):
        self.ring = ring
        self._sq = {}
        self._memo = {}
        self.deadline = None  # optional wall-clock deadline (time.time()); raises TimeoutError

    def squares_upto(self, R: QI) -> List[QI]:
        key = R.key()
        if key not in self._sq:
            self._sq[key] = list(self.ring.squares_bounded(R))  # descending at id
        return self._sq[key]

    def solve(self, R: QI, s: int, upper: Optional[QI] = None) -> List[Tuple[QI, ...]]:
        key = (R.key(), s, None if upper is None else upper.key())
        if key in self._memo:
            return self._memo[key]
        ring = self.ring
        out: List[Tuple[QI, ...]] = []
        if s == 0:
            if R.is_zero():
                out.append(())
        elif R.is_zero():
            out.append(tuple([ZERO] * s))
        elif ring.is_totally_nonneg(R):
            if self.deadline is not None and time.time() > self.deadline:
                raise TimeoutError("time_limit")
            for y in self.squares_upto(R):
                if upper is not None and y.gt(upper):
                    continue
                y2 = y * y
                # remaining s-1 entries are <= y at id, so their squares are <= y^2 at id
                if (R - y2).gt(y2 * (s - 1)):
                    break  # y only decreases from here on
                rest = R - y2
                if s == 1:
                    if rest.is_zero():
                        out.append((y,))
                    continue
                for tail in self.solve(rest, s - 1, y):
                    out.append((y,) + tail)
        self._memo[key] = out
        return out


@dataclass
class Root:
    k: Vec
    m: QI
    delta: QI           # d k0^2 / m  (sinh^2 of the distance to v0, times 1)
    batch: int

    def to_json(self):
        return {"k": [str(x) for x in self.k], "m": str(self.m), "delta": str(self.delta),
                "delta_float": self.delta.f(), "batch": self.batch}


@dataclass
class RunResult:
    ring: str
    d: str
    n: int
    status: str                         # reflective / undetermined
    roots: List[Root]
    stabiliser: Dict
    termination: Dict
    batches: int
    max_delta: float
    seconds: float
    caps: Dict
    labels: List[Tuple[int, int, object]] = field(default_factory=list)
    coxiter: str = ""

    def to_json(self):
        return {"ring": self.ring, "d": self.d, "n": self.n, "status": self.status,
                "n_roots": len(self.roots), "roots": [r.to_json() for r in self.roots],
                "stabiliser": self.stabiliser, "termination": self.termination,
                "batches": self.batches, "max_delta_processed": self.max_delta,
                "seconds": self.seconds, "caps": self.caps,
                "coxeter_edges": [[i, j, lab] for i, j, lab in self.labels],
                "coxiter_graph": self.coxiter}


def _scale_factor(ring: Ring, m: QI) -> QI:
    """c with: m | 2x  <=>  c | x, normalised positive at id."""
    g = ring.gcd(m, QI(2))
    c = m / g
    assert ring.contains(c)
    return c if c.sign() > 0 else -c


def stabiliser_check(ring: Ring, d: QI, n: int, norms: List[QI], sos: SumOfSquares) -> Dict:
    """All roots with k0 = 0.  Every such root is a signed permutation of a vector with
    entries y_1 >= ... >= y_n >= 0 at id; list those 'sorted' representatives and check
    they are exactly (1,0,..,0) [m=1] and (1,1,0,..,0) [m=2], i.e. the root system of the
    stabiliser of v0 is B_n with simple roots v_{i+1}-v_i, -v_n."""
    found = []
    for m in norms:
        c = _scale_factor(ring, m)
        R = m / (c * c)
        if not ring.contains(R):
            continue
        for y in sos.solve(R, n):
            k = (ZERO,) + tuple(c * t for t in y)
            if is_root(ring, d, k):
                found.append((str(m), [str(x) for x in k[1:]]))
    expected = [("1", ["1"] + ["0"] * (n - 1))]
    if n >= 2:
        expected.append(("2", ["1", "1"] + ["0"] * (n - 2)))
    ok = sorted(found) == sorted(expected)
    return {"type_Bn": ok, "sorted_roots": found}


def initial_roots(n: int) -> List[Vec]:
    out = []
    for i in range(1, n):
        k = [ZERO] * (n + 1)
        k[i], k[i + 1] = QI(-1), QI(1)
        out.append(tuple(k))
    k = [ZERO] * (n + 1)
    k[n] = QI(-1)
    out.append(tuple(k))
    return out


def k0_candidates(ring: Ring, d: QI, m: QI, delta_max: Fraction) -> List[QI]:
    """k0 > 0 at id with d k0^2/m <= delta_max (id) and sigma(m) + sigma(d) sigma(k0)^2 >= 0,
    and m | 2 d k0."""
    out = []
    A = math.sqrt(float(delta_max) * m.f() / d.f()) + 1e-9
    if ring is ZS2:
        B = math.sqrt(m.fs() / -d.fs()) + 1e-9
        box = ring.elements_in_box(A + 1, B + 1)
    else:
        box = (QI(a) for a in range(1, int(A) + 2))
    for k0 in box:
        if k0.sign() <= 0:
            continue
        if (d * k0 * k0 / m - QI(delta_max)).sign() > 0:
            continue
        if ring is ZS2 and (m + d * k0 * k0).sign_sigma() < 0:
            continue
        if not ring.divides(m, 2 * d * k0):
            continue
        out.append(k0)
    return out


def vinberg(ring: Ring, d: QI, n: int, *, delta_max: float = 200.0, max_roots: int = 200,
            max_batches: int = 100000, time_limit: float = 600.0,
            log: bool = False) -> RunResult:
    t0 = time.time()
    d = QI.of(d)
    if ring is ZS2:
        assert d.sign() > 0 and d.sign_sigma() < 0, "d must be admissible"
    else:
        assert d.b == 0 and d.sign() > 0
    compact = ring is ZS2
    norms = root_norms(ring, d)
    sos = SumOfSquares(ring)
    stab = stabiliser_check(ring, d, n, norms, sos)
    if not stab["type_Bn"]:
        raise NotImplementedError(f"stabiliser of v0 is not B_n: {stab}")

    diag = Diagram(n, allow_parallel=not compact, track_parabolic=not compact)
    accepted: List[Root] = []

    def add_root(k: Vec, m: QI, delta: QI, batch: int):
        row = [form(d, r.k, k) for r in accepted] + [m]
        diag.add(row, m)
        accepted.append(Root(k, m, delta, batch))

    for k in initial_roots(n):
        add_root(k, form(d, k, k), ZERO, 0)

    # batches: all (m, k0) with delta <= delta_max, ordered by delta (exact ties)
    dm = Fraction(delta_max).limit_denominator(10 ** 6)
    pairs = []
    for m in norms:
        for k0 in k0_candidates(ring, d, m, dm):
            pairs.append((d * k0 * k0 / m, m, k0))
    pairs.sort(key=lambda p: p[0].f())
    batches: List[List[Tuple[QI, QI, QI]]] = []
    for p in pairs:
        if batches and batches[-1][0][0] == p[0]:
            batches[-1].append(p)
        else:
            batches.append([p])

    status, term = "undetermined", diag.check_termination()
    nb = 0
    max_delta = 0.0
    cap_hit = "delta_max"
    for bi, batch in enumerate(batches, start=1):
        if nb >= max_batches:
            cap_hit = "max_batches"
            break
        if time.time() - t0 > time_limit:
            cap_hit = "time_limit"
            break
        if len(accepted) >= max_roots:
            cap_hit = "max_roots"
            break
        nb += 1
        delta = batch[0][0]
        max_delta = delta.f()
        cands: List[Tuple[Vec, QI]] = []
        for _, m, k0 in batch:
            c = _scale_factor(ring, m)
            r = m + d * k0 * k0
            R = r / (c * c)
            if not ring.contains(R):
                continue
            for y in sos.solve(R, n):
                k = (k0,) + tuple(c * t for t in y)
                if not ring.is_unit(ring.gcd_list(k)):
                    continue
                assert form(d, k, k) == m
                cands.append((k, m))
        # acceptance against all earlier batches
        survivors = [(k, m) for k, m in cands
                     if all(form(d, r.k, k).sign() <= 0 for r in accepted)]
        # same-distance survivors must be mutually non-obtuse (Vinberg); assert it
        for (k1, _), (k2, _) in itertools.combinations(survivors, 2):
            if form(d, k1, k2).sign() > 0:
                raise AssertionError(f"same-batch roots with positive product: {k1} {k2}")
        for k, m in survivors:
            assert is_root(ring, d, k)
            add_root(k, m, delta, bi)
        if survivors:
            term = diag.check_termination()
            if log:
                print(f"  batch {bi} delta={delta.f():.4f} +{len(survivors)} -> "
                      f"{len(accepted)} roots; {term}")
            if term["terminated"]:
                status = "reflective"
                cap_hit = None
                break
    else:
        cap_hit = "delta_max"

    return RunResult(ring=ring.name, d=str(d), n=n, status=status, roots=accepted,
                     stabiliser=stab, termination=term, batches=nb, max_delta=max_delta,
                     seconds=time.time() - t0,
                     caps={"delta_max": delta_max, "max_roots": max_roots,
                           "max_batches": max_batches, "time_limit": time_limit,
                           "stopped_by": cap_hit},
                     labels=diag.coxeter_matrix_rows(), coxiter=diag.to_coxiter())
