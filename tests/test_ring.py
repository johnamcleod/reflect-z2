import itertools
import random
from fractions import Fraction

import pytest

from reflect.ring import (EPS, ONE, R2, ZERO, QI, ZS2, ZZ, _sgn, cos2, primes_above)

rng = random.Random(12345)


def rand_qi(lo=-30, hi=30):
    return QI(rng.randint(lo, hi), rng.randint(lo, hi))


def test_ring_axioms():
    for _ in range(300):
        x, y, z = rand_qi(), rand_qi(), rand_qi()
        assert (x + y) + z == x + (y + z)
        assert (x * y) * z == x * (y * z)
        assert x * (y + z) == x * y + x * z
        assert x * y == y * x
        assert x - x == ZERO
        assert x * ONE == x
        assert (x * y).conj() == x.conj() * y.conj()
        assert (x * y).norm() == x.norm() * y.norm()
        if not y.is_zero():
            assert (x / y) * y == x


def test_sqrt2_squared():
    assert R2 * R2 == QI(2)
    assert EPS.norm() == -1
    assert EPS.conj() == QI(1, -1)


def test_exact_sign_matches_float_far_from_zero():
    for _ in range(2000):
        x = rand_qi(-1000, 1000)
        if abs(x.f()) > 1e-6:
            assert x.sign() == (1 if x.f() > 0 else -1)
        if abs(x.fs()) > 1e-6:
            assert x.sign_sigma() == (1 if x.fs() > 0 else -1)


def test_exact_sign_near_zero():
    # Pell-type elements a + b r2 with a^2 - 2b^2 = +-1: |x| tiny at one embedding
    u = EPS ** 30          # huge at id, tiny at sigma
    tiny = u.conj()        # = (1 - r2)^30 > 0, about 1e-12
    assert tiny.sign() == 1
    assert (-tiny).sign() == -1
    assert (EPS ** 31).conj().sign() == -1  # (1-r2)^31 < 0
    # rational components
    assert _sgn(Fraction(1, 3), Fraction(-1, 5)) == 1  # 1/3 - r2/5 = 0.333-0.283
    assert _sgn(Fraction(-7, 5), 1) == 1                 # -1.4 + 1.414
    assert _sgn(Fraction(-17, 12), 1) == -1              # -1.41666 + 1.41421


def test_divisibility_and_gcd():
    for _ in range(200):
        x, y = rand_qi(), rand_qi()
        if x.is_zero() or y.is_zero():
            continue
        g = ZS2.gcd(x, y)
        assert ZS2.divides(g, x) and ZS2.divides(g, y)
        # g generates (x, y): Bezout-free check via divisibility of every common divisor
        for d in (R2, QI(3), QI(3, 1), QI(5, 2)):
            if ZS2.divides(d, x) and ZS2.divides(d, y):
                assert ZS2.divides(d, g)
    assert ZS2.divides(R2, QI(2))
    assert not ZS2.divides(QI(2), R2)
    assert ZS2.divides(QI(3, 1), QI(7))
    assert ZZ.divides(QI(3), QI(12)) and not ZZ.divides(QI(5), QI(12))


def test_units():
    for k in range(-5, 6):
        u = EPS ** k
        assert ZS2.is_unit(u)
        assert ZS2.unit_log(u) == (1, k)
        assert ZS2.unit_log(-u) == (-1, k)
    assert not ZS2.is_unit(R2)


def test_totally_positive_units_are_squares():
    assert ZS2.totally_positive_units_mod_squares_is_trivial()
    # brute-force confirmation on a range of exponents
    for s, k in itertools.product((1, -1), range(-8, 9)):
        u = s * EPS ** k
        assert u.is_totally_positive() == (s == 1 and k % 2 == 0)


def test_reduce_mod_unit_squares():
    assert ZS2.reduce_mod_unit_squares(EPS) == EPS          # tie with r2 - 1, keep 1+r2
    assert ZS2.reduce_mod_unit_squares(EPS ** 7) == EPS
    assert ZS2.reduce_mod_unit_squares(QI(-1, 1)) == EPS    # r2-1 = EPS^-1 ~ EPS (mod EPS^2)
    for _ in range(100):
        x = rand_qi(1, 20)
        y = x * EPS ** (2 * rng.randint(-4, 4))
        assert ZS2.reduce_mod_unit_squares(x) == ZS2.reduce_mod_unit_squares(y)


def test_primes_above():
    assert primes_above(2) == (R2,)
    assert primes_above(3) == (QI(3),)
    assert primes_above(5) == (QI(5),)
    for p in (7, 17, 23, 31, 41, 47):
        p1, p2 = primes_above(p)
        assert abs(p1.norm()) == p and abs(p2.norm()) == p
        assert not ZS2.divides(p1, p2)
        assert ZS2.divides(p1 * p2, QI(p))


def test_factor():
    for x in (QI(14, 7), QI(6), QI(2), QI(4, 6), QI(-9, 1)):
        u, fs = ZS2.factor(x)
        y = u
        for p, e in fs:
            y = y * p ** e
        assert y == x and ZS2.is_unit(u)


def test_squares_bounded():
    r = QI(5, 2)  # 5 + 2 r2 ~ 7.83, sigma ~ 2.17
    xs = list(ZS2.squares_bounded(r))
    for x in xs:
        assert ZS2.is_totally_nonneg(r - x * x) and x.sign() >= 0
    # brute force
    brute = [x for x in ZS2.elements_in_box(4, 4)
             if x.sign() >= 0 and ZS2.is_totally_nonneg(r - x * x)]
    assert set(xs) == set(brute)
    assert QI(1, 1) in xs  # (1+r2)^2 = 3 + 2 r2 <= 5 + 2r2, sigma: 3-2r2 <= 5-2r2


def test_cos2_values():
    assert cos2(8) == QI(Fraction(1, 2), Fraction(1, 4))
    assert (cos2(8) * 4) == QI(2, 1)


def test_parse_roundtrip():
    for _ in range(100):
        x = QI(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
        assert QI.parse(str(x)) == x
        assert QI.from_json(x.to_json()) == x
