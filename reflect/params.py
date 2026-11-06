"""Phase 1: admissible prime parameters d of Z[sqrt2] with N(p) <= bound, and the
possible root norms for each (totally positive divisors of 2d modulo unit squares)."""
from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, field
from typing import List

from .ring import EPS, ONE, R2, QI, ZS2, primes_above


def admissible(d: QI) -> bool:
    return d.sign() > 0 and d.sign_sigma() < 0


def make_admissible(x: QI) -> QI:
    """Multiply x by a unit so that the sign pattern becomes (+, -)."""
    for u in (ONE, EPS, -ONE, -EPS):
        y = x * u
        if admissible(y):
            return y
    raise AssertionError("unreachable: units realise every sign pattern")


def canonical(d: QI) -> QI:
    """Canonical admissible representative of d * EPS^(2Z) (admissibility is preserved
    by multiplication by EPS^2, which is totally positive)."""
    assert admissible(d)
    c = ZS2.reduce_mod_unit_squares(d)
    assert admissible(c)
    return c


def admissible_classes(gen: QI, krange=range(-6, 7)) -> List[QI]:
    """All admissible associates u*gen (u = +-EPS^k, k in krange), grouped modulo unit
    squares; returns the list of canonical representatives (expected length 1)."""
    reps = set()
    for s, k in itertools.product((1, -1), krange):
        y = s * gen * EPS ** k
        if admissible(y):
            reps.add(ZS2.reduce_mod_unit_squares(y))
    return sorted(reps, key=lambda z: (z.f(), z.fs()))


def totally_positive_associate(x: QI) -> QI:
    for u in (ONE, EPS, -ONE, -EPS):
        y = x * u
        if y.is_totally_positive():
            return ZS2.reduce_mod_unit_squares(y)
    raise AssertionError("unreachable")


def root_norms(d: QI) -> List[QI]:
    """Totally positive divisors of 2d modulo EPS^(2Z), canonical representatives,
    sorted by value at id."""
    u, fs = ZS2.factor(2 * d)
    out = set()
    for exps in itertools.product(*[range(e + 1) for _, e in fs]):
        x = ONE
        for (p, _), k in zip(fs, exps):
            x = x * p ** k
        out.add(totally_positive_associate(x))
    return sorted(out, key=lambda z: z.f())


@dataclass
class Param:
    label: str          # e.g. "p7a", "inert3", "ram2", "unit"
    kind: str           # split / inert / ramified / unit
    p: int              # rational prime below (0 for the unit case)
    ideal_gen: QI       # a generator of the ideal (before sign adjustment)
    ideal_norm: int     # N(ideal) (1 for the unit case)
    d: QI               # canonical admissible generator
    root_norms: List[QI] = field(default_factory=list)

    def to_json(self):
        return {"label": self.label, "kind": self.kind, "p": self.p,
                "ideal_gen": str(self.ideal_gen), "ideal_norm": self.ideal_norm,
                "d": str(self.d), "d_float": [self.d.f(), self.d.fs()],
                "N_d": self.d.norm(),
                "root_norms": [str(m) for m in self.root_norms]}


def enumerate_params(bound: int = 50) -> List[Param]:
    out: List[Param] = []
    for p in range(2, bound + 1):
        if any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
            continue
        gens = primes_above(p)
        if p == 2:
            kind, normp = "ramified", 2
        elif p % 8 in (3, 5):
            kind, normp = "inert", p * p
        else:
            kind, normp = "split", p
        if normp > bound:
            continue
        for i, g in enumerate(gens):
            assert abs(g.norm()) == normp
            classes = admissible_classes(g)
            assert len(classes) == 1, (p, g, classes)
            d = canonical(make_admissible(g))
            assert d == classes[0]
            label = {"ramified": "ram2", "inert": f"inert{p}"}.get(kind, f"p{p}{'ab'[i]}")
            out.append(Param(label, kind, p, g, normp, d, root_norms(d)))
    d = canonical(EPS)
    out.append(Param("unit", "unit", 0, ONE, 1, d, root_norms(d)))
    return out


def main():
    import pathlib
    ps = enumerate_params()
    path = pathlib.Path(__file__).resolve().parent.parent / "data" / "params.json"
    path.write_text(json.dumps([p.to_json() for p in ps], indent=1))
    for p in ps:
        print(f"{p.label:8s} {p.kind:9s} N(id)={p.ideal_norm:3d}  d={str(p.d):10s} "
              f"(id {p.d.f():8.4f}, sigma {p.d.fs():8.4f})  N(d)={p.d.norm():4d}  "
              f"m in {{{', '.join(str(m) for m in p.root_norms)}}}")


if __name__ == "__main__":
    main()
