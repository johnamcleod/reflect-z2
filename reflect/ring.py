"""Exact arithmetic in K = Q(sqrt2) and O_K = Z[sqrt2] (and Z / Q as the b = 0 case).

Elements are QI(a, b) = a + b*sqrt2 with a, b ints or Fractions.  Every sign decision
is exact (compare a^2 with 2 b^2); floats are produced only by `.f()` / `.fs()` for
sorting and display.

Two ring objects implement the interface used by the Vinberg code:
  ZZ  -- the rational integers (one real embedding), elements QI(a, 0)
  ZS2 -- Z[sqrt2] with embeddings id: sqrt2 -> +sqrt2 and sigma: sqrt2 -> -sqrt2
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache
from typing import Iterator, List, Optional, Tuple

SQRT2 = math.sqrt(2.0)


def _sgn(a, b) -> int:
    """Exact sign of a + b*sqrt2 for rational a, b."""
    if a >= 0 and b >= 0:
        return 0 if (a == 0 and b == 0) else 1
    if a <= 0 and b <= 0:
        return -1
    t = a * a - 2 * b * b  # never 0 here: sqrt2 is irrational and (a,b) != 0
    if a > 0:
        return 1 if t > 0 else -1
    return -1 if t > 0 else 1


def _norm_frac(x):
    if isinstance(x, Fraction) and x.denominator == 1:
        return int(x.numerator)
    return x


class QI:
    """a + b*sqrt2, a, b rational (int or Fraction)."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = _norm_frac(a)
        self.b = _norm_frac(b)

    # ---- construction helpers -------------------------------------------------
    @staticmethod
    def of(x) -> "QI":
        if isinstance(x, QI):
            return x
        if isinstance(x, (int, Fraction)):
            return QI(x, 0)
        if isinstance(x, (tuple, list)) and len(x) == 2:
            return QI(x[0], x[1])
        if isinstance(x, str):
            return QI.parse(x)
        raise TypeError(f"cannot convert {x!r} to QI")

    @staticmethod
    def parse(s: str) -> "QI":
        """Parse 'a', 'b*r2', 'a+b*r2', 'a-b*r2' (r2 = sqrt2), with rational a, b."""
        s = s.replace(" ", "").replace("sqrt2", "r2").replace("√2", "r2")
        if "r2" not in s:
            return QI(Fraction(s), 0)
        # split at the last +/- that is not at position 0 and not after '/'
        idx = None
        for i in range(len(s) - 1, 0, -1):
            if s[i] in "+-" and s[i - 1] not in "*/eE":
                idx = i
                break
        if idx is None:
            a_str, b_str = "0", s
        else:
            a_str, b_str = s[:idx], s[idx:]
            if "r2" in a_str:
                a_str, b_str = b_str, a_str
        b_str = b_str.replace("*r2", "").replace("r2", "")
        if b_str in ("", "+"):
            b_str = "1"
        elif b_str == "-":
            b_str = "-1"
        return QI(Fraction(a_str), Fraction(b_str))

    # ---- ring operations -------------------------------------------------------
    def __add__(self, o):
        o = QI.of(o)
        return QI(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.a, -self.b)

    def __sub__(self, o):
        o = QI.of(o)
        return QI(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return QI.of(o) - self

    def __mul__(self, o):
        o = QI.of(o)
        return QI(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __pow__(self, k: int):
        if k < 0:
            return (QI(1) / self) ** (-k)
        r, base = QI(1), self
        while k:
            if k & 1:
                r = r * base
            base = base * base
            k >>= 1
        return r

    def conj(self) -> "QI":
        """Galois conjugate sigma(x)."""
        return QI(self.a, -self.b)

    def norm(self):
        """N(x) = x * sigma(x) = a^2 - 2 b^2 (rational)."""
        return _norm_frac(self.a * self.a - 2 * self.b * self.b)

    def __truediv__(self, o):
        o = QI.of(o)
        n = o.norm()
        if n == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt2)")
        num = self * o.conj()
        return QI(Fraction(num.a) / n, Fraction(num.b) / n)

    def __rtruediv__(self, o):
        return QI.of(o) / self

    # ---- predicates ------------------------------------------------------------
    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_integral(self) -> bool:
        return isinstance(self.a, int) and isinstance(self.b, int)

    def is_rational(self) -> bool:
        return self.b == 0

    def sign(self) -> int:
        """Sign at the identity embedding."""
        return _sgn(self.a, self.b)

    def sign_sigma(self) -> int:
        """Sign at sigma."""
        return _sgn(self.a, -self.b)

    def is_totally_positive(self) -> bool:
        return self.sign() > 0 and self.sign_sigma() > 0

    # exact comparisons at id
    def gt(self, o) -> bool:
        return (self - o).sign() > 0

    def ge(self, o) -> bool:
        return (self - o).sign() >= 0

    def lt(self, o) -> bool:
        return (self - o).sign() < 0

    def le(self, o) -> bool:
        return (self - o).sign() <= 0

    # floats -- display / ordering heuristics only
    def f(self) -> float:
        return float(self.a) + float(self.b) * SQRT2

    def fs(self) -> float:
        return float(self.a) - float(self.b) * SQRT2

    # ---- misc ------------------------------------------------------------------
    def __eq__(self, o):
        try:
            o = QI.of(o)
        except TypeError:
            return NotImplemented
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __bool__(self):
        return not self.is_zero()

    def key(self) -> Tuple:
        return (self.a, self.b)

    def __repr__(self):
        return f"QI({self.a}, {self.b})"

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        bs = "" if self.b == 1 else "-" if self.b == -1 else f"{self.b}*"
        bterm = f"{bs}r2" if bs in ("", "-") else f"{bs}r2"
        if self.a == 0:
            return bterm
        sign = "+" if (self.b > 0) else ""
        return f"{self.a}{sign}{bterm}"

    def to_json(self):
        return [str(self.a), str(self.b)]

    @staticmethod
    def from_json(v) -> "QI":
        return QI(Fraction(v[0]), Fraction(v[1]))


ZERO, ONE = QI(0), QI(1)
R2 = QI(0, 1)
EPS = QI(1, 1)  # fundamental unit 1 + sqrt2, N(EPS) = -1


def cos2(k: int) -> QI:
    """Exact cos^2(pi/k) for the angles that can occur over Q(sqrt2): k in {2,3,4,6,8}."""
    table = {2: QI(0), 3: QI(Fraction(1, 4)), 4: QI(Fraction(1, 2)), 6: QI(Fraction(3, 4)),
             8: QI(Fraction(1, 2), Fraction(1, 4))}
    return table[k]


def isqrt_floor_frac(x: Fraction) -> int:
    """floor(sqrt(x)) for rational x >= 0, exact."""
    if x < 0:
        raise ValueError
    n, d = x.numerator, x.denominator
    r = math.isqrt(n * d) // d
    while (r + 1) * (r + 1) * d <= n:
        r += 1
    while r * r * d > n:
        r -= 1
    return r


class Ring:
    """Interface over which the Vinberg code is written.  Two instances: ZZ and ZS2."""

    name: str
    real_quadratic: bool

    def __init__(self, name: str, real_quadratic: bool):
        self.name = name
        self.real_quadratic = real_quadratic

    def __repr__(self):
        return self.name

    # embeddings ---------------------------------------------------------------
    def n_embeddings(self) -> int:
        return 2 if self.real_quadratic else 1

    def is_totally_positive(self, x: QI) -> bool:
        x = QI.of(x)
        if self.real_quadratic:
            return x.is_totally_positive()
        return x.sign() > 0

    def is_totally_nonneg(self, x: QI) -> bool:
        x = QI.of(x)
        if self.real_quadratic:
            return x.sign() >= 0 and x.sign_sigma() >= 0
        return x.sign() >= 0

    # membership / divisibility --------------------------------------------------
    def contains(self, x: QI) -> bool:
        x = QI.of(x)
        if not x.is_integral():
            return False
        return self.real_quadratic or x.b == 0

    def divides(self, x: QI, y: QI) -> bool:
        """x | y in the ring (x != 0)."""
        x, y = QI.of(x), QI.of(y)
        if x.is_zero():
            return y.is_zero()
        return self.contains(y / x)

    def norm_abs(self, x: QI) -> int:
        """Euclidean function: |N(x)| for Z[sqrt2], |x| for Z."""
        x = QI.of(x)
        return abs(x.norm()) if self.real_quadratic else abs(x.a)

    def divmod(self, y: QI, x: QI) -> Tuple[QI, QI]:
        """Euclidean division y = q x + r with norm_abs(r) < norm_abs(x)."""
        y, x = QI.of(y), QI.of(x)
        t = y / x
        qa = round(Fraction(t.a))
        qb = round(Fraction(t.b)) if self.real_quadratic else 0
        q = QI(qa, qb)
        r = y - q * x
        assert self.norm_abs(r) < self.norm_abs(x), "Euclidean step failed"
        return q, r

    def gcd(self, x: QI, y: QI) -> QI:
        x, y = QI.of(x), QI.of(y)
        while not y.is_zero():
            _, r = self.divmod(x, y)
            x, y = y, r
        return x

    def gcd_list(self, xs) -> QI:
        g = ZERO
        for x in xs:
            g = self.gcd(g, QI.of(x))
            if self.is_unit(g):
                return ONE
        return g

    # units ---------------------------------------------------------------------
    def is_unit(self, x: QI) -> bool:
        x = QI.of(x)
        return self.contains(x) and self.norm_abs(x) == 1

    def fundamental_unit(self) -> Optional[QI]:
        return EPS if self.real_quadratic else None

    def unit_log(self, u: QI) -> Tuple[int, int]:
        """Write a unit u = s * EPS^k, return (s, k).  For ZZ, k = 0."""
        u = QI.of(u)
        assert self.is_unit(u), u
        if not self.real_quadratic:
            return (1 if u.a > 0 else -1), 0
        s = 1 if u.sign() > 0 else -1
        v = u * s
        k = 0
        # v > 0; reduce towards 1 using |v| at id
        while v.gt(ONE):
            v = v / EPS
            k += 1
        while v.lt(ONE):
            v = v * EPS
            k -= 1
        assert v == ONE
        return s, k

    def associates(self, x: QI, kmin: int = -3, kmax: int = 3) -> List[QI]:
        x = QI.of(x)
        if not self.real_quadratic:
            return [x, -x]
        return [s * x * EPS ** k for s in (1, -1) for k in range(kmin, kmax + 1)]

    def totally_positive_units_mod_squares_is_trivial(self) -> bool:
        """Exact proof-by-computation that the totally positive units are the squares.

        Every unit is +-EPS^k.  sigma(EPS) = 1 - sqrt2 < 0 and EPS > 0, so EPS^k has
        sign pattern (+, (-1)^k); -EPS^k has (-, -(-1)^k).  Totally positive iff the
        sign is + and k is even, i.e. iff the unit is a square (EPS^(k/2))^2.
        We check the sign patterns of +-EPS^0, +-EPS^1 exactly (the pattern is
        multiplicative so this covers every k).
        """
        if not self.real_quadratic:
            return True
        assert EPS.sign() == 1 and EPS.sign_sigma() == -1
        assert EPS.norm() == -1
        pats = {}
        for s in (1, -1):
            for k in (0, 1):
                u = s * EPS ** k
                pats[(s, k)] = (u.sign(), u.sign_sigma())
        return pats == {(1, 0): (1, 1), (1, 1): (1, -1), (-1, 0): (-1, -1), (-1, 1): (-1, 1)}

    def reduce_mod_unit_squares(self, x: QI) -> QI:
        """Canonical representative of the orbit x * EPS^(2Z).

        Minimise S(y) = max(y^2, sigma(y)^2) (compared exactly at id); ties are broken
        by the larger value of |y| at id.  The sign of x is kept.  S(x EPS^(2k)) is a max
        of an increasing and a decreasing geometric sequence in k, hence unimodal, so an
        exact local walk finds the minimum.
        """
        x = QI.of(x)
        if not self.real_quadratic or x.is_zero():
            return x
        e2 = EPS * EPS

        def S(y: QI) -> QI:
            p, q = y * y, y.conj() * y.conj()
            return p if p.ge(q) else q

        y = x
        for step in (e2, ONE / e2):
            while S(y * step).lt(S(y)):
                y = y * step
        # y is a minimiser; collect ties among neighbours
        cands = [y] + [c for c in (y * e2, y / e2) if S(c) == S(y)]
        best = cands[0]
        for c in cands[1:]:
            if abs_gt(c, best):
                best = c
        return best

    # enumeration -----------------------------------------------------------------
    def squares_bounded(self, r: QI, upper: Optional[QI] = None, lower_zero: bool = True
                        ) -> Iterator[QI]:
        """All x in the ring with r - x^2 totally >= 0, 0 <= x at id (if lower_zero),
        and x <= upper at id (if given).  r is an element of the ring."""
        r = QI.of(r)
        if not self.is_totally_nonneg(r):
            return
        if not self.real_quadratic:
            hi = isqrt_floor_frac(Fraction(r.a))
            if upper is not None:
                hi = min(hi, int(math.floor(Fraction(upper.a))))
            lo = 0 if lower_zero else -hi
            for a in range(hi, lo - 1, -1):
                yield QI(a, 0)
            return
        # |x| <= sqrt(r) =: A, |sx| <= sqrt(sr) =: B.  x = a + b r2:
        # a = (x + sx)/2, b = (x - sx)/(2 r2)  =>  |a| <= (A+B)/2, |b| <= (A+B)/(2 r2)
        A = math.sqrt(max(r.f(), 0.0))
        B = math.sqrt(max(r.fs(), 0.0))
        amax = int((A + B) / 2 + 2)
        bmax = int((A + B) / (2 * SQRT2) + 2)
        out = []
        for b in range(-bmax, bmax + 1):
            for a in range(-amax, amax + 1):
                x = QI(a, b)
                if lower_zero and x.sign() < 0:
                    continue
                if upper is not None and x.gt(upper):
                    continue
                if self.is_totally_nonneg(r - x * x):
                    out.append(x)
        out.sort(key=lambda z: -z.f())
        yield from out

    def elements_in_box(self, A: float, B: float) -> Iterator[QI]:
        """Elements x with |x| <= A, |sigma x| <= B (float box, for heuristics/tests)."""
        if not self.real_quadratic:
            for a in range(-int(A), int(A) + 1):
                yield QI(a)
            return
        amax = int((A + B) / 2 + 1)
        bmax = int((A + B) / (2 * SQRT2) + 1)
        for a in range(-amax, amax + 1):
            for b in range(-bmax, bmax + 1):
                x = QI(a, b)
                if abs(x.f()) <= A + 1e-12 and abs(x.fs()) <= B + 1e-12:
                    yield x

    # factorisation -------------------------------------------------------------
    def factor(self, x: QI) -> Tuple[QI, List[Tuple[QI, int]]]:
        """Factor a nonzero ring element as unit * prod p_i^e_i (primes up to associates,
        each normalised by reduce_mod_unit_squares and positive at id)."""
        x = QI.of(x)
        assert self.contains(x) and not x.is_zero()
        if not self.real_quadratic:
            n = abs(x.a)
            fs = []
            p = 2
            while p * p <= n:
                e = 0
                while n % p == 0:
                    n //= p
                    e += 1
                if e:
                    fs.append((QI(p), e))
                p += 1
            if n > 1:
                fs.append((QI(n), 1))
            u = QI(1 if x.a > 0 else -1)
            return u, fs
        N = abs(x.norm())
        primes = []
        for p, _ in _int_factor(N):
            for pi in primes_above(p):
                if pi not in primes:
                    primes.append(pi)
        fs = []
        y = x
        for pi in primes:
            e = 0
            while self.divides(pi, y):
                y = y / pi
                e += 1
            if e:
                fs.append((pi, e))
        assert self.is_unit(y), (x, y)
        return y, fs


def abs_gt(p: QI, q: QI) -> bool:
    """|p| > |q| at id, exact."""
    return (p * p - q * q).sign() > 0


def abs_ge(p: QI, q: QI) -> bool:
    return (p * p - q * q).sign() >= 0


def _int_factor(n: int) -> List[Tuple[int, int]]:
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


@lru_cache(maxsize=None)
def primes_above(p: int) -> Tuple[QI, ...]:
    """Generators (positive at id, minimal size) of the prime ideals of Z[sqrt2] above p."""
    if p == 2:
        return (R2,)
    r = p % 8
    if r in (3, 5):
        return (QI(p),)
    assert r in (1, 7), p
    # solve a^2 - 2 b^2 = +-p
    gens = []
    bmax = int(math.isqrt(p)) + 2
    for b in range(1, bmax):
        for s in (1, -1):
            a2 = s * p + 2 * b * b
            if a2 > 0:
                a = math.isqrt(a2)
                if a * a == a2:
                    gens.append(QI(a, b))
                    gens.append(QI(a, -b))
        if gens:
            break
    assert gens, p
    # the two primes are (a + b r2) and (a - b r2); they are not associates
    pi1, pi2 = gens[0], gens[1]
    assert not ZS2.divides(pi1, pi2), "split primes should be distinct"
    return (pi1, pi2)


ZZ = Ring("Z", real_quadratic=False)
ZS2 = Ring("Z[sqrt2]", real_quadratic=True)
