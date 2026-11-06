"""Phase 1 (brief 2): candidate parameters d from Mark's table.

    python -m reflect.candidates            -> data/candidates.json

C = { canonical(-delta(entry)) }.  By Mark's theorem every squarefree admissible d with
<-d,1,1> reflective lies in C; each d in C is then decided directly (reflect.decide2).
"""
from __future__ import annotations

import json
from collections import defaultdict
from typing import Dict, List

from .mark import ROOT, load
from .params import admissible, canonical, enumerate_params, root_norms
from .ring import QI, ZS2

OUT = ROOT / "data" / "candidates.json"


def prime_kind(pi: QI) -> str:
    n = abs(pi.norm())
    if n == 2:
        return "ramified"
    r = int(round(n ** 0.5))
    if r * r == n:          # norm q^2: an inert rational prime q (q = 3, 5 mod 8)
        assert r % 8 in (3, 5)
        return "inert"
    assert n % 8 in (1, 7)  # norm p: a prime above a split p
    return "split"


def label_for(d: QI, grid: Dict[str, str]) -> str:
    if str(d) in grid:
        return grid[str(d)]
    a, b = int(d.a), int(d.b)
    f = lambda x: f"m{-x}" if x < 0 else str(x)
    return f"d{f(a)}_{f(b)}"


def build() -> List[Dict]:
    grid = {str(p.d): p.label for p in enumerate_params()}
    by_d = defaultdict(list)
    for e in load():
        d = canonical(-e["delta"])
        assert admissible(d), (e["id"], d)
        by_d[str(d)].append(e)
    out = []
    for ds, es in by_d.items():
        d = QI.parse(ds)
        unit, fac = ZS2.factor(d)
        assert all(k == 1 for _, k in fac)
        primes = [{"pi": str(p), "norm": abs(p.norm()), "kind": prime_kind(p)} for p, _ in fac]
        out.append({
            "label": label_for(d, grid), "d": ds, "d_float": [d.f(), d.fs()], "N_d": d.norm(),
            "primes": primes, "in_grid": ds in grid,
            "root_norms": [str(m) for m in root_norms(d)],
            "entries": [e["id"] for e in es],
            "shapes": sorted({e["shape"] for e in es}),
        })
    out.sort(key=lambda c: (abs(c["N_d"]), c["d_float"][0]))
    return out


def main():
    c = build()
    OUT.write_text(json.dumps(c, indent=1))
    for x in c:
        print(f"{x['label']:10s} d={x['d']:10s} N={x['N_d']:5d}  "
              f"{'·'.join(str(p['norm']) + p['kind'][0] for p in x['primes']) or 'unit':14s} "
              f"entries={len(x['entries']):3d} grid={x['in_grid']}")
    print(len(c), "candidates")


if __name__ == "__main__":
    main()
