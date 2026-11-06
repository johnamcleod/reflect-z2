"""7b: does an infinite-order symmetry of P come from a reflective super-lattice?

For a pure translation h in Sym(P) (rank(h - 1) = 2, loxodromic) and any space-like
a in U = im(h - 1), R_a h R_a = h^-1 and R_b := R_a h is again a reflection.  Then
L' = L + R_a(L) is stable under R_a, R_b and h (because h in O(L) and R_a h = R_b,
R_a h^-1 = h R_a), so h = R_a R_b lies in W(L').  a is never a root of L (otherwise
h would lie in W(L) and fix P pointwise).  We then decide whether L' is reflective
with the general-lattice Vinberg algorithm below.
"""
from __future__ import annotations

import itertools
import json
import math
import pathlib
import time
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .certificates import (_bijections, check_isometry, finite_order_bound, gram_F,
                           isometry_from_vertices, loxodromic_bracket, vertex_vector)
from .coxeter import Diagram
from .linalg import (eye, from_columns, inverse, kernel, matmul, matvec, rank, transpose)
from .ring import ONE, ZERO, QI, ZS2, xgcd
from .vinberg import SumOfSquares, form

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "superlattices"


# ---------------------------------------------------------------- symmetries ----------
def vertex_symmetries(d: QI, n: int, roots, limit: int = 400):
    """All integral g in O'(L), g != 1, mapping the v0-vertex onto a found vertex
    (hence g(P) = P), over the found roots."""
    diag = Diagram(n, allow_parallel=False, track_parabolic=False)
    for i, k in enumerate(roots):
        row = [form(d, roots[j], k) for j in range(i)] + [form(d, k, k)]
        diag.add(row, row[-1])
    S0 = list(range(n))
    G0 = [[diag.Gq[i][j] for j in S0] for i in S0]
    norms0 = [diag.norms[i] for i in S0]
    v0 = [ONE] + [ZERO] * n
    out = []
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
            if v1 is None:
                v1 = vertex_vector(d, [roots[i] for i in Sl])
            img = [Sl[c] for c in pi]
            g = isometry_from_vertices(d, n, [roots[i] for i in S0], v0,
                                       [roots[i] for i in img], v1)
            if g is None:
                continue
            chk = check_isometry(d, n, g)
            if chk["integral"] and chk["preserves_form"] and chk["future"] and not chk["identity"]:
                out.append({"g": g, "to": img})
                if len(out) >= limit:
                    return out
    return out


def _key(g):
    return tuple(x.key() for row in g for x in row)


def pure_translations(syms, n: int, max_len: int = 2, max_out: int = 20, max_gens: int = 10):
    """Words of length <= max_len in the (max_gens nearest) symmetries that are pure
    translations, i.e. rank(h - 1) = 2 and loxodromic (also tested for h^2)."""
    gens = []
    seen = set()
    # nearest first: cosh dist(v0, g v0) is proportional to g[0][0]
    syms = sorted(syms, key=lambda s: s["g"][0][0].f())[:max_gens]
    for s in syms:
        for g in (s["g"], inverse(s["g"])):
            k = _key(g)
            if k not in seen:
                seen.add(k)
                gens.append(g)
    I = eye(n + 1)
    found, fkeys = [], set()
    frontier = [I]
    allk = {_key(I)}
    for length in range(1, max_len + 1):
        new = []
        for w in frontier:
            for g in gens:
                h = matmul(w, g)
                k = _key(h)
                if k in allk:
                    continue
                allk.add(k)
                new.append(h)
                for p in (1, 2):
                    hp = h if p == 1 else matmul(h, h)
                    M = [[hp[i][j] - I[i][j] for j in range(n + 1)] for i in range(n + 1)]
                    if rank(M) == 2 and loxodromic_bracket(hp) is not None:
                        kk = _key(hp)
                        if kk not in fkeys and _key(inverse(hp)) not in fkeys:
                            fkeys.add(kk)
                            found.append(hp)
                            if len(found) >= max_out:
                                return found
        frontier = new
        if len(allk) > 20000:
            break
    return found


# ---------------------------------------------------------------- reflections -------
def reflection(d: QI, a: Sequence[QI]):
    n1 = len(a)
    aa = form(d, a, a)
    F = gram_F(d, n1 - 1)
    Fa = matvec(F, a)
    return [[(ONE if i == j else ZERO) - QI(2) * a[i] * Fa[j] / aa for j in range(n1)]
            for i in range(n1)]


def decompose(d: QI, n: int, h, candidates_x) -> List[Dict]:
    """Decompositions h = R_a R_b with a = x - h x for x in candidates_x."""
    I = eye(n + 1)
    out = []
    for x in candidates_x:
        hx = matvec(h, x)
        a = [p - q for p, q in zip(x, hx)]
        if all(t.is_zero() for t in a):
            continue
        if form(d, a, a).sign() <= 0:
            continue
        Ra = reflection(d, a)
        Rb = matmul(Ra, h)
        M = [[Rb[i][j] - I[i][j] for j in range(n + 1)] for i in range(n + 1)]
        assert matmul(Rb, Rb) == I and rank(M) == 1, "R_a h is not a reflection"
        # b spans the image of Rb - 1
        b = next(col for col in transpose(M) if any(not t.is_zero() for t in col))
        out.append({"a": a, "b": b, "Ra": Ra, "Rb": Rb})
    return out


# ---------------------------------------------------------------- HNF over O_K ------
def _content_den(vs) -> QI:
    """A D in O_K with D * v integral for all v (product of denominators, crude)."""
    den = 1
    for v in vs:
        for x in v:
            for c in (Fraction(x.a), Fraction(x.b)):
                den = den * c.denominator // math.gcd(den, c.denominator)
    return QI(den)


def hnf_basis(gens: List[List[QI]], dim: int) -> List[List[QI]]:
    """Basis (columns) of the O_K-module generated by integral column vectors gens,
    assumed of full rank dim; lower-triangular Hermite form."""
    cols = [list(g) for g in gens]
    basis = []
    for r in range(dim):
        # combine all columns with nonzero entry in row r into one
        piv = None
        rest = []
        for c in cols:
            if c[r].is_zero():
                rest.append(c)
                continue
            if piv is None:
                piv = c
                continue
            g, s, t = xgcd(ZS2, piv[r], c[r])
            u, v = piv[r] / g, c[r] / g
            newp = [s * x + t * y for x, y in zip(piv, c)]
            other = [v * x - u * y for x, y in zip(piv, c)]   # entry r becomes 0
            assert other[r].is_zero()
            piv = newp
            rest.append(other)
        assert piv is not None, "generators do not span"
        basis.append(piv)
        cols = [c for c in rest if any(not x.is_zero() for x in c)]
    return from_columns(basis)


def superlattice_basis(d: QI, n: int, Ra) -> Tuple[List[List[QI]], QI, Dict]:
    """L' = L + R_a L.  Returns (basis matrix B with columns in V coords, D with D B
    integral, info)."""
    E = eye(n + 1)
    gens = [list(c) for c in transpose(E)] + [list(c) for c in transpose(Ra)]
    D = _content_den(gens)
    gens_int = [[D * x for x in g] for g in gens]
    B_int = hnf_basis(gens_int, n + 1)
    B = [[x / D for x in row] for row in B_int]
    # index [L' : L] = N(det B^-1) ... record det
    from .coxeter import determinant
    det = determinant(B)
    return B, D, {"det_B": str(det), "index_norm": str(abs((ONE / det).norm()))}


def residue_set(B, D: QI):
    """Image of D*L' in (O/D)^(n+1) (D a rational integer, L = O^(n+1) <= L' <= D^-1 L).
    y in L is in D*L' iff its componentwise residue mod D lies in this set."""
    Dn = int(D.a)
    assert D.b == 0 and Dn > 0
    cols = [[D * B[i][j] for i in range(len(B))] for j in range(len(B))]
    gens = []
    for c in cols:
        for mult in (QI(1), QI(0, 1)):
            v = [mult * x for x in c]
            assert all(ZS2.contains(x) for x in v)
            gens.append(tuple((int(x.a) % Dn, int(x.b) % Dn) for x in v))
    zero = tuple((0, 0) for _ in cols)
    S = {zero}
    frontier = [zero]
    while frontier:
        nxt = []
        for s in frontier:
            for g in gens:
                u = tuple(((a + c) % Dn, (b + e) % Dn) for (a, b), (c, e) in zip(s, g))
                if u not in S:
                    S.add(u)
                    nxt.append(u)
        frontier = nxt
    return Dn, S


def in_lattice(Binv, x) -> Optional[List[QI]]:
    c = matvec(Binv, x)
    return c if all(ZS2.contains(t) for t in c) else None


# ---------------------------------------------------------------- general Vinberg ---
def _signed_perms(y: Sequence[QI]):
    seen = set()
    for p in set(itertools.permutations(range(len(y)))):
        base = [y[i] for i in p]
        nz = [i for i, t in enumerate(base) if not t.is_zero()]
        for signs in itertools.product((1, -1), repeat=len(nz)):
            v = list(base)
            for i, s in zip(nz, signs):
                if s < 0:
                    v[i] = -v[i]
            k = tuple(t.key() for t in v)
            if k not in seen:
                seen.add(k)
                yield v


def vinberg_general(d: QI, n: int, B, D: QI, *, delta_max: float = 2000.0,
                    max_roots: int = 120, time_limit: float = 900.0,
                    max_perms: int = 3_000_000, prune: bool = True) -> Dict:
    """Vinberg's algorithm for the lattice L' with basis B (columns, V coordinates) in
    the space with form diag(-d,1,...,1), control vector v0.  Exact."""
    t0 = time.time()
    F = gram_F(d, n)
    Fp = matmul(transpose(B), matmul(F, B))           # Gram of L'
    den = _content_den(Fp)
    c = den * den                                      # totally positive scaling
    Fpc = [[c * x for x in row] for row in Fp]
    assert all(ZS2.contains(x) for row in Fpc for x in row)
    # exponent of L'^#/L': lcm of denominators of Fpc^-1
    inv = inverse(Fpc)
    e = ONE
    for row in inv:
        for q in row:
            if q.is_zero():
                continue
            # denominator ideal of q: smallest t with t q integral
            num_den = _content_den([[q]])
            t = num_den
            g = ZS2.gcd(t, t * q) if not (t * q).is_zero() else t
            t = t / g if ZS2.contains(t / g) else t
            # make lcm(e, t)
            e = e * t / ZS2.gcd(e, t)
    from .params import totally_positive_associate
    _, fs = ZS2.factor(2 * e)
    norms = set()
    for ex in itertools.product(*[range(k + 1) for _, k in fs]):
        x = ONE
        for (p, _), k in zip(fs, ex):
            x = x * p ** k
        norms.add(totally_positive_associate(x))
    norms = sorted(norms, key=lambda z: z.f())
    Binv = inverse(B)
    Dn, RES = residue_set(B, D)

    def qform(x, y):  # scaled form c * (x, y)
        return c * form(d, x, y)

    def is_root_Lp(x):
        co = in_lattice(Binv, x)
        if co is None:
            return False, None
        if not ZS2.is_unit(ZS2.gcd_list(co)):
            return False, None
        m = qform(x, x)
        if not m.is_totally_positive():
            return False, None
        for j in range(n + 1):
            bj = [B[i][j] for i in range(n + 1)]
            if not ZS2.contains(2 * qform(x, bj) / m):
                return False, None
        return True, m

    sos = SumOfSquares(ZS2)
    nperm = 0

    def scale_i(m):
        """Roots x of L' with c(x,x) = m satisfy m | 2c x_i (pair with v_i in L <= L'),
        so y_i = D x_i is a multiple of s = D m / gcd(m, 2c)."""
        if not prune:
            return D
        g = ZS2.gcd(m, 2 * c)
        s = D * m / g
        return s if s.sign() > 0 else -s

    def scale_0(m):
        if not prune:
            return ONE
        g = ZS2.gcd(m, 2 * c * d)
        s = D * m / g
        return s if s.sign() > 0 else -s

    def candidates(m, y0):
        """x in L' with c (x,x) = m and x_0 = y0 / D."""
        nonlocal nperm
        R0 = D * D * m / c                      # (y,y) for y = D x
        if not ZS2.contains(R0):
            return
        if not ZS2.divides(scale_0(m), y0):
            return
        R = R0 + d * y0 * y0
        s = scale_i(m)
        Rz = R / (s * s)
        if not ZS2.contains(Rz):
            return
        r0 = (int(y0.a) % Dn, int(y0.b) % Dn)
        for zs in sos.solve(Rz, n):
            ys = [s * z for z in zs]
            for v in _signed_perms(ys):
                nperm += 1
                if nperm > max_perms:
                    raise TimeoutError("signed-permutation budget exceeded")
                res = (r0,) + tuple((int(t.a) % Dn, int(t.b) % Dn) for t in v)
                if res not in RES:
                    continue
                x = [y0 / D] + [t / D for t in v]
                if any(form(d, x, s_).sign() > 0 for s_ in simple_now):
                    continue
                ok, mm = is_root_Lp(x)
                if ok and mm == m:
                    yield x

    # stabiliser of v0: roots with x_0 = 0, chamber containing a generic rho
    simple_now: List = []
    stab = []
    for m in norms:
        stab.extend(candidates(m, ZERO))
    rho = [ZERO] + [QI(Fraction((n + 1 - i) * 1000 + i * i, 997)) for i in range(1, n + 1)]
    for r in stab:
        assert not form(d, r, rho).is_zero(), "rho not generic"
    neg = [r for r in stab if form(d, r, rho).sign() < 0]
    neg.sort(key=lambda r: (form(d, r, rho) * form(d, r, rho) / qform(r, r)).f())
    simple = []
    for r in neg:
        if all(qform(r, s).sign() <= 0 for s in simple):
            simple.append(r)
    diag = Diagram(n, allow_parallel=False, track_parabolic=False)
    acc = []

    def add(x):
        row = [qform(a, x) for a in acc] + [qform(x, x)]
        diag.add(row, row[-1])
        acc.append(x)

    for s in simple:
        add(s)
    simple_now.extend(simple)
    # k0 candidates: x_0 = y0/D with y0 > 0
    pairs = []
    for m in norms:
        A = math.sqrt(delta_max * (m / c).f() / d.f()) * D.f() + 2
        Bs = math.sqrt(max((D * D * m / c).fs(), 0) / -d.fs()) + 2
        for y0 in ZS2.elements_in_box(A, Bs):
            if y0.sign() <= 0:
                continue
            delta = d * y0 * y0 * c / (D * D * m)
            if delta.f() > delta_max:
                continue
            if (D * D * m / c + d * y0 * y0).sign_sigma() < 0:
                continue
            pairs.append((delta, m, y0))
    pairs.sort(key=lambda p: p[0].f())
    batches = []
    for p in pairs:
        if batches and batches[-1][0][0] == p[0]:
            batches[-1].append(p)
        else:
            batches.append([p])
    status, term = "undetermined", diag.check_termination() if diag.ell_top else None
    stop = "delta_max"
    try:
        for batch in batches:
            if time.time() - t0 > time_limit:
                stop = "time_limit"
                break
            if len(acc) >= max_roots:
                stop = "max_roots"
                break
            surv = []
            for _, m, y0 in batch:
                for x in candidates(m, y0):
                    if all(qform(a, x).sign() <= 0 for a in acc):
                        surv.append(x)
            for x in surv:
                add(x)
            if surv:
                term = diag.check_termination()
                if term["terminated"]:
                    status, stop = "reflective", None
                    break
    except TimeoutError as ex:
        stop = str(ex)
    return {"status": status, "n_roots": len(acc), "stop": stop,
            "stabiliser_simple": [[str(t) for t in s] for s in simple],
            "roots": [[str(t) for t in x] for x in acc], "termination": term,
            "norms": [str(m) for m in norms], "scale": str(c),
            "coxiter_graph": diag.to_coxiter(), "seconds": time.time() - t0}


# ---------------------------------------------------------------- driver -----------
def analyse(label: str, d: QI, n: int, roots, log=print, vinberg_kw=None) -> Dict:
    syms = vertex_symmetries(d, n, roots)
    trans = pure_translations(syms, n)
    res = {"label": label, "d": str(d), "n": n, "n_symmetries": len(syms),
           "n_pure_translations": len(trans), "decompositions": []}
    if not trans:
        res["outcome"] = "no decomposition (no pure translation among short words)"
        return res
    # candidate points x: v0 and the found vertices' vectors
    v0 = [ONE] + [ZERO] * n
    xs = [v0]
    for s in syms[:20]:
        xs.append(matvec(s["g"], v0))
    best = None
    for h in trans[:5]:
        for dec in decompose(d, n, h, xs)[:6]:
            assert not _is_root_L(d, dec["a"]) and not _is_root_L(d, dec["b"])
            B, D, info = superlattice_basis(d, n, dec["Ra"])
            entry = {"h": [[str(x) for x in r] for r in h],
                     "a": [str(x) for x in dec["a"]], "b": [str(x) for x in dec["b"]],
                     "basis": [[str(x) for x in r] for r in B], **info}
            Binv = inverse(B)
            # sanity: R_a, R_b, h preserve L'
            for M in (dec["Ra"], dec["Rb"], h):
                MB = matmul(M, B)
                assert all(in_lattice(Binv, [MB[i][j] for i in range(n + 1)]) is not None
                           for j in range(n + 1)), "L' not stable"
            vg = vinberg_general(d, n, B, D, **(vinberg_kw or {}))
            entry["vinberg"] = {k: v for k, v in vg.items() if k not in ("roots", "coxiter_graph")}
            entry["vinberg_roots"] = vg["roots"]
            entry["coxiter_graph"] = vg["coxiter_graph"]
            res["decompositions"].append(entry)
            log(f"  {label} n={n}: L' index-norm {info['index_norm']}: {vg['status']} "
                f"({vg['n_roots']} roots, stop={vg['stop']}, {vg['seconds']:.1f}s)")
            if vg["status"] == "reflective":
                best = entry
                break
        if best:
            break
    res["outcome"] = ("reflective super-lattice found" if best else
                      "decompositions found; no reflective L' within caps")
    return res


def _is_root_L(d, a):
    from .vinberg import is_root
    g = ZS2.gcd_list(a)
    if g.is_zero():
        return False
    # scale a to a primitive integral vector
    den = _content_den([a])
    v = [den * x for x in a]
    g = ZS2.gcd_list(v)
    v = [x / g for x in v]
    return is_root(ZS2, d, v)
