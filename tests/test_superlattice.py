"""7b machinery: HNF over Z[sqrt2], decomposition h = R_a R_b, stability of L' = L + R_a L,
and equality of the pruned / unpruned general-lattice Vinberg enumeration."""
import json
import pathlib

import pytest

from reflect.linalg import eye, inverse, matmul, rank
from reflect.ring import ONE, ZERO, QI, ZS2
from reflect.run import run_path
from reflect.superlattice import (_is_root_L, decompose, hnf_basis, in_lattice,
                                  pure_translations, superlattice_basis, vertex_symmetries,
                                  vinberg_general)
from reflect.vinberg import form

D17 = QI.parse("-1+3*r2")


def test_hnf_spans_same_module():
    gens = [[QI(2), QI(0), QI(0)], [QI(0), QI(1, 1), QI(0)], [QI(0), QI(0), QI(3)],
            [QI(1), QI(1), QI(1)], [QI(0, 1), QI(2), QI(1, -1)]]
    B = hnf_basis(gens, 3)
    Binv = inverse(B)
    for g in gens:  # every generator lies in the lattice spanned by B
        assert in_lattice(Binv, g) is not None
    # and every basis column is an O-combination of generators (B lies in the module):
    # here the generators contain e1..e3 up to units/2,3 plus (1,1,1) -> the module is O^3
    assert ZS2.is_unit(B[0][0] * B[1][1] * B[2][2])


@pytest.fixture(scope="module")
def p17a_n3():
    p = run_path("Zr2", D17, 3)
    if not p.exists():
        pytest.skip("grid run missing")
    roots = [tuple(QI.parse(x) for x in r["k"]) for r in json.loads(p.read_text())["roots"]][:60]
    syms = vertex_symmetries(D17, 3, roots)
    trans = pure_translations(syms, 3)
    assert trans
    return roots, syms, trans


def test_decomposition_and_stability(p17a_n3):
    _, _, trans = p17a_n3
    n = 3
    h = trans[0]
    I = eye(n + 1)
    assert rank([[h[i][j] - I[i][j] for j in range(n + 1)] for i in range(n + 1)]) == 2
    decs = decompose(D17, n, h, [[ONE] + [ZERO] * n])
    assert decs
    dec = decs[0]
    assert matmul(dec["Ra"], dec["Rb"]) == h
    assert not _is_root_L(D17, dec["a"]) and not _is_root_L(D17, dec["b"])
    B, D, info = superlattice_basis(D17, n, dec["Ra"])
    Binv = inverse(B)
    for e in [[ONE if i == j else ZERO for i in range(n + 1)] for j in range(n + 1)]:
        assert in_lattice(Binv, e) is not None          # L <= L'
    for M in (dec["Ra"], dec["Rb"], h):
        MB = matmul(M, B)
        for j in range(n + 1):
            assert in_lattice(Binv, [MB[i][j] for i in range(n + 1)]) is not None


def test_pruned_enumeration_equals_unpruned(p17a_n3):
    _, _, trans = p17a_n3
    n = 3
    dec = decompose(D17, n, trans[0], [[ONE] + [ZERO] * n])[0]
    B, D, _ = superlattice_basis(D17, n, dec["Ra"])
    a = vinberg_general(D17, n, B, D, delta_max=300, max_roots=500)
    b = vinberg_general(D17, n, B, D, delta_max=300, max_roots=500, prune=False)
    assert a["roots"] == b["roots"] and a["n_roots"] > n
