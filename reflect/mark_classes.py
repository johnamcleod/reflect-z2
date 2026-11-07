"""Isometry classes of the entries of Mark's table (OPEN_QUESTIONS OQ4: 522 entries vs the
432 lattices of the theorem).

    python -m reflect.mark_classes      -> data/mark_classes.json

Invariants (determinant class, angle cycle up to dihedral symmetry, multiset of root-norm
classes) split the entries into groups; isometric entries lie in the same group.  Within
a group, entries are compared by chamber matching: an isometry phi: L_A -> L_B, composed
with an element of W(L_B), maps A's chamber onto B's, so A's listed simple roots go to
u_i * (simple roots of B) with units u_i > 0.  We enumerate these injections against the
FULL simple-root cycle of B, solve phi on three independent roots, and accept phi iff
phi, phi^-1 are integral and phi^T G_B phi = G_A (an explicit isometry certificate).

Full cycles: entries with shape 'abc...' list all sides; for '(..)^2' entries the full
cycle is computed with reflect.vinberg3 when needed (cached in data/mark_full.json).
A negative answer (no injection gives an isometry) relies on B's cycle being the complete
set of simple roots of a chamber: for computed cycles this is our Vinberg run; for fully
listed entries it is Mark's list, which matches our own Vinberg wherever we ran it.
"""
from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from typing import Dict, List, Optional

from .linalg import from_columns, inverse, matmul, rank, sqrt_K, transpose
from .mark import ROOT, ip, load, shape_base
from .mark_full import OUT as FULL, full_chamber
from .params import canonical
from .ring import QI, ZS2

OUT = ROOT / "data" / "mark_classes.json"


def canon_shape(s):
    b, k = shape_base(s)
    full = b * k
    c = []
    for t in (full, full[::-1]):
        c += [t[i:] + t[:i] for i in range(len(t))]
    return min(c)


def invariant(e):
    b, k = shape_base(e["shape"])
    norms = sorted(str(ZS2.reduce_mod_unit_squares(m)) for m in e["M"]) * k
    return (str(canonical(-e["delta"])), canon_shape(e["shape"]), tuple(sorted(norms)))


def _unit_scale(M: QI, m: QI) -> Optional[QI]:
    r = M / m
    if not ZS2.contains(r) or not ZS2.is_unit(r):
        return None
    u = sqrt_K(r)
    if u is None:
        return None
    return u if u.sign() > 0 else -u


def isometry(A, B_G, B_roots) -> Optional[List[List[QI]]]:
    """phi: L_A -> L_B (columns in B's basis) mapping A's listed roots to unit multiples
    of B's full root cycle, or None."""
    GA, RA, MA = A["G"], A["R"], A["M"]
    h = len(RA)
    nB = [ip(B_G, x, x) for x in B_roots]
    GBB = [[ip(B_G, a, b) for b in B_roots] for a in B_roots]
    GAA = [[ip(GA, RA[i], RA[j]) for j in range(h)] for i in range(h)]
    opts = [[(j, u) for j in range(len(B_roots)) for u in [_unit_scale(MA[i], nB[j])]
             if u is not None] for i in range(h)]
    idx = None
    for a in range(h):
        for b in range(a + 1, h):
            for c in range(b + 1, h):
                if idx is None and rank(from_columns([RA[a], RA[b], RA[c]])) == 3:
                    idx = (a, b, c)
    Rinv = inverse(from_columns([RA[i] for i in idx]))

    def rec(i, used, pi):
        if i == h:
            img = [[pi[t][1] * x for x in B_roots[pi[t][0]]] for t in range(h)]
            phi = matmul(from_columns([img[t] for t in idx]), Rinv)
            if not all(ZS2.contains(x) for row in phi for x in row):
                return None
            if matmul(transpose(phi), matmul(B_G, phi)) != GA:
                return None
            pinv = inverse(phi)
            if pinv is None or not all(ZS2.contains(x) for row in pinv for x in row):
                return None
            return phi
        for j, u in opts[i]:
            if j in used:
                continue
            if all(pi[a][1] * u * GBB[pi[a][0]][j] == GAA[a][i] for a in range(i)):
                pi.append((j, u))
                used.add(j)
                res = rec(i + 1, used, pi)
                if res is not None:
                    return res
                used.discard(j)
                pi.pop()
        return None
    return rec(0, set(), [])


def main():
    E = load()
    full = json.loads(FULL.read_text()) if FULL.exists() else {}

    def cycle(e):
        if shape_base(e["shape"])[1] == 1:
            return e["R"]
        key = str(e["id"])
        if key not in full or not full[key]["ok"]:
            full[key] = full_chamber(e)
            FULL.write_text(json.dumps(full, indent=0))
            assert full[key]["ok"], f"full chamber failed for entry {e['id']}"
        return [[QI.parse(x) for x in r] for r in full[key]["roots"]]

    groups = defaultdict(list)
    for e in E:
        groups[invariant(e)].append(e["id"])
    parent = {e["id"]: e["id"] for e in E}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    certs = []
    for key, ids in groups.items():
        if len(ids) == 1:
            continue
        # fully listed members first (no Vinberg needed), then half-listed ones
        ids = sorted(ids, key=lambda i: (shape_base(E[i]["shape"])[1], i))
        reps: List[int] = []
        for i in ids:
            hit = None
            for r in reps:
                phi = isometry(E[i], E[r]["G"], cycle(E[r]))
                if phi is not None:
                    hit = (r, phi)
                    break
            if hit:
                parent[find(i)] = find(hit[0])
                certs.append({"from": i, "to": hit[0],
                              "phi": [[str(x) for x in row] for row in hit[1]]})
            else:
                reps.append(i)
        print(f"group {key[0]:10s} {key[1]:14s} size {len(ids):2d} -> {len(reps)} classes",
              flush=True)
    classes = defaultdict(list)
    for e in E:
        classes[find(e["id"])].append(e["id"])
    OUT.write_text(json.dumps({"n_entries": len(E), "n_classes": len(classes),
                               "classes": sorted(classes.values()),
                               "isometries": certs}, indent=0))
    print(f"{len(E)} entries -> {len(classes)} isometry classes")


if __name__ == "__main__":
    main()
