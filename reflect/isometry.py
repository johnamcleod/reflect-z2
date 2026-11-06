"""Calibration (d): match Mark's table entries with <-d,1,1> (n = 2 Vinberg polygons).

    python -m reflect.isometry       -> data/isometries/<label>.json, summary on stdout

For a reflective <-d,1,1> we know ALL simple roots of our chamber P.  If phi: L -> <-d,1,1>
is an isometry, then composing with an element of W(<-d,1,1>) we may assume phi maps
Mark's chamber onto P, so phi maps Mark's listed (consecutive) simple roots injectively
onto our simple roots, preserving norms and inner products.  We enumerate those injections,
solve for phi on 3 independent roots and keep phi iff phi, phi^-1 are integral and
phi^T F phi = G.  A match is an explicit isometry certificate.  'No match for an entry'
proves non-isometry only if Mark's listed roots are simple roots of a chamber of L
(Mark's claim; checked here only for consistency: they are roots, pairwise non-obtuse,
with the printed angles).
"""
from __future__ import annotations

import json
from typing import Dict, List, Optional

from .candidates import OUT as CAND
from .linalg import from_columns, inverse, matmul, rank, sqrt_K, transpose
from .mark import ROOT, ip, load
from .ring import QI, ZS2
from .run import run_path
from .vinberg import form

OUT = ROOT / "data" / "isometries"


def _unit_scale(M: QI, m: QI) -> Optional[QI]:
    """The unit u > 0 (at id) with u^2 m = M, if any.  A simple root is determined by its
    mirror only up to such a factor (u s is primitive with norm u^2 m)."""
    r = M / m
    if not ZS2.contains(r) or not ZS2.is_unit(r):
        return None
    u = sqrt_K(r)
    if u is None:
        return None
    return u if u.sign() > 0 else -u


def _injections(GM, RM, MM, d, ours):
    """Injective maps i -> (j, u_i) with r_i |-> u_i s_j preserving all inner products."""
    n_ours = [form(d, x, x) for x in ours]
    h = len(RM)
    G_ours = [[form(d, a, b) for b in ours] for a in ours]
    G_mark = [[ip(GM, RM[i], RM[j]) for j in range(h)] for i in range(h)]
    opts = [[(j, u) for j in range(len(ours))
             for u in [_unit_scale(MM[i], n_ours[j])] if u is not None] for i in range(h)]

    def rec(i, used, pi):
        if i == h:
            yield list(pi)
            return
        for j, u in opts[i]:
            if j in used:
                continue
            if all(pi[a][1] * u * G_ours[pi[a][0]][j] == G_mark[a][i] for a in range(i)):
                pi.append((j, u))
                used.add(j)
                yield from rec(i + 1, used, pi)
                used.discard(j)
                pi.pop()
    yield from rec(0, set(), [])


def isometry_for(entry, d: QI, ours) -> Optional[Dict]:
    GM, RM, MM = entry["G"], entry["R"], entry["M"]
    F = [[(-d if i == j == 0 else QI(1) if i == j else QI(0)) for j in range(3)] for i in range(3)]
    # three independent listed roots
    idx = None
    for a in range(len(RM)):
        for b in range(a + 1, len(RM)):
            for c in range(b + 1, len(RM)):
                if rank(from_columns([RM[a], RM[b], RM[c]])) == 3:
                    idx = (a, b, c)
                    break
            if idx:
                break
        if idx:
            break
    if idx is None:
        return None
    Rinv = inverse(from_columns([RM[i] for i in idx]))
    for pi in _injections(GM, RM, MM, d, ours):
        img = [[pi[i][1] * x for x in ours[pi[i][0]]] for i in range(len(RM))]
        S = from_columns([img[i] for i in idx])
        phi = matmul(S, Rinv)
        if not all(ZS2.contains(x) for row in phi for x in row):
            continue
        if matmul(transpose(phi), matmul(F, phi)) != GM:
            continue
        pinv = inverse(phi)
        if pinv is None or not all(ZS2.contains(x) for row in pinv for x in row):
            continue
        # every listed root goes where pi says
        if any([sum((phi[r][c] * RM[i][c] for c in range(3)), QI(0)) for r in range(3)]
               != img[i] for i in range(len(RM))):
            continue
        return {"entry": entry["id"], "shape": entry["shape"],
                "map": [[j, str(u)] for j, u in pi],
                "phi": [[str(x) for x in row] for row in phi]}
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    entries = load()
    summary = []
    for c in json.loads(CAND.read_text()):
        d = QI.parse(c["d"])
        p = run_path("Zr2", d, 2)
        if not p.exists():
            summary.append((c["label"], "no run"))
            continue
        run = json.loads(p.read_text())
        if run["status"] != "reflective":
            summary.append((c["label"], f"<-d,1,1> {run['status']}: no entry can match"))
            continue
        ours = [tuple(QI.parse(x) for x in r["k"]) for r in run["roots"]]
        matches = [m for m in (isometry_for(entries[i], d, ours) for i in c["entries"]) if m]
        (OUT / f"{c['label']}.json").write_text(json.dumps(
            {"label": c["label"], "d": c["d"], "n_sides_ours": len(ours),
             "entries_checked": c["entries"], "matches": matches}, indent=1))
        summary.append((c["label"], f"{len(matches)} of {len(c['entries'])} entries isometric"
                        + ("" if matches else "  <-- NO MATCH")))
    for s in summary:
        print(f"{s[0]:10s} {s[1]}")


if __name__ == "__main__":
    main()
