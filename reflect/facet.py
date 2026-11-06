"""Calibration (e): the facet u^perp, u = -v_{n+1}, of the (d, n+1) polytope reproduces the
(d, n) polytope (docs/monotonicity.md, step 4).

For each facet root e meeting u (elliptic pair, angle pi/k): k odd -> error; otherwise the
wall of F through that ridge is the primitive projection of e onto u^perp, i.e. e with its
last coordinate set to 0, made primitive.  The resulting set must equal the (d, n) roots
up to positive unit multiples.
"""
from __future__ import annotations

import json
from typing import List, Optional, Sequence, Tuple

from .coxeter import DOTTED, edge_label
from .linalg import sqrt_K
from .ring import QI, ZS2
from .run import run_path
from .vinberg import form, is_root


def _primitive(v: Sequence[QI]) -> Tuple[QI, ...]:
    g = ZS2.gcd_list(list(v))
    return tuple(x / g for x in v)


def same_line(a, b) -> bool:
    """b = u a for a unit u > 0 at id."""
    i = next(j for j, x in enumerate(a) if not x.is_zero())
    if b[i].is_zero():
        return False
    u = b[i] / a[i]
    return ZS2.contains(u) and ZS2.is_unit(u) and u.sign() > 0 and all(
        u * x == y for x, y in zip(a, b))


def facet_roots(d: QI, roots_up: List[Sequence[QI]]) -> List[Tuple[QI, ...]]:
    N = len(roots_up[0])
    u = tuple([QI(0)] * (N - 1) + [QI(-1)])
    out = []
    for e in roots_up:
        e = tuple(e)
        if e == u:
            continue
        g = form(d, e, u)
        lab = edge_label(g, form(d, e, e), QI(1), False)
        if lab == DOTTED:
            continue
        if lab in (3,):                    # odd angle at the facet: impossible (lemma)
            raise AssertionError(f"odd angle pi/{lab} at the facet u^perp: {e}")
        w = _primitive(e[:-1])       # projection onto u^perp, dropped last coordinate
        if not any(same_line(w, x) for x in out):
            out.append(w)
    return out


def check_pair(d: QI, n: int) -> Optional[str]:
    """None if (d,n+1) facet reproduces (d,n); else a message.  Both runs must exist."""
    up = json.loads(run_path("Zr2", d, n + 1).read_text())
    lo = json.loads(run_path("Zr2", d, n).read_text())
    assert up["status"] == "reflective" and lo["status"] == "reflective"
    R_up = [tuple(QI.parse(x) for x in r["k"]) for r in up["roots"]]
    R_lo = [tuple(QI.parse(x) for x in r["k"]) for r in lo["roots"]]
    F = facet_roots(d, R_up)
    for w in F:
        if not is_root(ZS2, d, w):
            return f"projected wall {w} is not a root of L_n"
    miss = [r for r in R_lo if not any(same_line(r, w) for w in F)]
    extra = [w for w in F if not any(same_line(r, w) for r in R_lo)]
    if miss or extra:
        return f"{len(miss)} roots of (d,n) missing, {len(extra)} extra walls"
    return None
