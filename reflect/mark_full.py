"""Full chambers of all entries of Mark's table, by our own rank-3 Vinberg (vinberg3).

    python -m reflect.mark_full [--jobs 4]   -> data/mark_full.json

For each entry the control vector is an interior point of Mark's polygon, built from the
listed roots: the sum of >= 3 listed vertices, or (two vertices) their sum pushed off the
common wall.  Checks per entry: Vinberg terminates; Mark's listed roots are all simple
roots of the chamber found (up to positive unit multiples); the number of sides equals
the side count of the shape string.  This also verifies that the listed roots are simple
roots of a chamber (OPEN_QUESTIONS OQ5).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

from .facet import same_line
from .linalg import kernel
from .mark import ROOT, load, shape_base
from .ring import QI
from .vinberg3 import ipG, vinberg3

OUT = ROOT / "data" / "mark_full.json"


def _vertex(G, a, b, others):
    rows = [[sum((G[i][j] * a[i] for i in range(3)), QI(0)) for j in range(3)],
            [sum((G[i][j] * b[i] for i in range(3)), QI(0)) for j in range(3)]]
    k = kernel(rows)
    assert len(k) == 1
    p = k[0]
    assert ipG(G, p, p).sign() < 0, "vertex not timelike"
    for o in others:
        s = ipG(G, p, o).sign()
        if s:
            return p if s < 0 else [-x for x in p]
    raise AssertionError("cannot orient vertex")


def _integral(x):
    from math import lcm
    D = 1
    for t in x:
        D = lcm(D, Fraction(t.a).denominator, Fraction(t.b).denominator)
    return [t * D for t in x]


def interior_point(e, eps):
    G, R = e["G"], e["R"]
    base, k = shape_base(e["shape"])
    h = len(R)
    pairs = [(i, i + 1) for i in range(h - 1)] + ([(h - 1, 0)] if k == 1 else [])
    verts = [_vertex(G, R[a], R[b], [R[c] for c in range(h) if c not in (a, b)])
             for a, b in pairs]
    x = [sum((v[i] for v in verts), QI(0)) for i in range(3)]
    zero = [r for r in R if ipG(G, x, r).is_zero()]
    x = [x[i] - eps * sum((r[i] for r in zero), QI(0)) for i in range(3)]
    assert all(ipG(G, x, r).sign() < 0 for r in R) and ipG(G, x, x).sign() < 0
    return _integral(x)


def full_chamber(e):
    t = time.time()
    G, R = e["G"], e["R"]
    base, k = shape_base(e["shape"])
    n_sides = len(base) * k
    last = None
    for eps in (QI(Fraction(1, 100)), QI(Fraction(1, 10000)), QI(Fraction(1, 10 ** 6))):
        v0 = interior_point(e, eps)
        res = vinberg3(G, v0, Delta_max=1e7, max_roots=60)
        roots = res["roots"]
        listed_found = all(any(same_line(r, s) for s in roots) for r in R)
        last = {"id": e["id"], "status": res["status"], "n_sides": len(roots),
                "expected_sides": n_sides, "listed_found": listed_found,
                "v0": [str(x) for x in v0], "eps": str(eps),
                "roots": [[str(x) for x in r] for r in roots],
                "seconds": round(time.time() - t, 2)}
        if res["status"] == "reflective" and listed_found:
            break
    last["ok"] = (last["status"] == "reflective" and last["listed_found"]
                  and last["n_sides"] == n_sides)
    return last


def main(argv):
    jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 4
    ids = [int(a) for a in argv if a.isdigit() and "--jobs" not in argv[max(0, argv.index(a) - 1):argv.index(a)]]
    E = load()
    todo = [E[i] for i in ids] if ids else E
    done = json.loads(OUT.read_text()) if OUT.exists() else {}
    todo = [e for e in todo if str(e["id"]) not in done or not done[str(e["id"])]["ok"]]
    with Pool(jobs) as pool:
        for r in pool.imap_unordered(full_chamber, todo):
            done[str(r["id"])] = r
            OUT.write_text(json.dumps(done, indent=0))
            if not r["ok"]:
                print(f"entry {r['id']}: {r['status']} sides {r['n_sides']}/{r['expected_sides']} "
                      f"listed_found={r['listed_found']}", flush=True)
    bad = [k for k, v in done.items() if not v["ok"]]
    print(f"{len(done)} entries, {len(bad)} not ok: {bad[:20]}")


if __name__ == "__main__":
    main(sys.argv[1:])
