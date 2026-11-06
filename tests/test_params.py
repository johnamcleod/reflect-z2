import itertools

from reflect.params import (admissible, admissible_classes, enumerate_params, root_norms)
from reflect.ring import EPS, QI, ZS2


PARAMS = enumerate_params(50)


def test_parameter_count_and_kinds():
    kinds = [p.kind for p in PARAMS]
    assert kinds.count("ramified") == 1
    assert kinds.count("inert") == 2            # p = 3, 5 (norms 9, 25)
    assert kinds.count("split") == 12           # 2 primes above each of 7,17,23,31,41,47
    assert kinds.count("unit") == 1
    assert sorted({p.p for p in PARAMS if p.kind == "split"}) == [7, 17, 23, 31, 41, 47]


def test_generators_admissible_and_generate_ideal():
    for p in PARAMS:
        assert admissible(p.d)
        assert p.d.is_integral()
        assert abs(p.d.norm()) == p.ideal_norm
        # d is an associate of the ideal generator
        assert ZS2.divides(p.ideal_gen, p.d) and ZS2.divides(p.d, p.ideal_gen)


def test_exactly_one_admissible_class_per_ideal():
    for p in PARAMS:
        assert admissible_classes(p.d) == [p.d]


def test_split_primes_distinct():
    by_p = {}
    for p in PARAMS:
        if p.kind == "split":
            by_p.setdefault(p.p, []).append(p.d)
    for ds in by_p.values():
        a, b = ds
        assert not ZS2.divides(a, b)                  # different ideals
        assert a != b and a != b * EPS ** 2


def test_unit_case():
    unit = [p for p in PARAMS if p.kind == "unit"][0]
    assert unit.d == QI(1, 1)


def test_root_norms():
    for p in PARAMS:
        ms = root_norms(p.d)
        assert ms == p.root_norms
        assert QI(1) in ms and QI(2) in ms and QI(2, 1) in ms
        for m in ms:
            assert m.is_totally_positive()
            assert ZS2.divides(m, 2 * p.d)
        # pairwise distinct modulo unit squares
        for m1, m2 in itertools.combinations(ms, 2):
            q = m1 / m2
            assert not (ZS2.is_unit(q) and q.is_totally_positive())
        # completeness: number of classes = number of divisors of 2d up to units
        _, fs = ZS2.factor(2 * p.d)
        ndiv = 1
        for _, e in fs:
            ndiv *= e + 1
        assert len(ms) == ndiv
