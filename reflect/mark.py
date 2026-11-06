"""Mark's table (arXiv:1512.01133, Appendix B): loading and exact consistency checks.

data/mark_table.json is produced by tools/parse_mark.py from the TeX source.
"""
from __future__ import annotations

import json
import pathlib
import re
from typing import Dict, List

from .coxeter import determinant, is_positive_definite
from .ring import QI, ZS2

ROOT = pathlib.Path(__file__).resolve().parent.parent
TABLE = ROOT / "data" / "mark_table.json"

COS2 = {2: QI(0), 3: QI(1) / 4, 4: QI(1) / 2, 6: QI(3) / 4, 8: QI(2, 1) / 4}


def q(p) -> QI:
    return QI(p[0], p[1])


def load() -> List[Dict]:
    out = []
    for e in json.loads(TABLE.read_text())["entries"]:
        e = dict(e)
        e["G"] = [[q(x) for x in r] for r in e["gram"]]
        e["R"] = [[q(x) for x in r] for r in e["roots"]]
        e["M"] = [q(x) for x in e["norms"]]
        e["delta"] = q(e["det"])
        out.append(e)
    return out


def ip(G, x, y) -> QI:
    return sum((x[i] * G[i][j] * y[j] for i in range(3) for j in range(3)), QI(0))


def sigma_matrix(G):
    return [[x.conj() for x in r] for r in G]


def angle_label(G, a, b) -> str:
    g = ip(G, a, b)
    c = g * g / (ip(G, a, a) * ip(G, b, b))
    for k, v in COS2.items():
        if c == v:
            return str(k)
    return "dotted" if (c - 1).sign() > 0 else "?"


def shape_base(shape: str):
    m = re.fullmatch(r"\((\d+)\)\^(\d+)", shape)
    return (m.group(1), int(m.group(2))) if m else (shape, 1)


def check_entry(e) -> Dict[str, bool]:
    """All exact checks for one entry; returns {check: passed}."""
    G, R, M, delta = e["G"], e["R"], e["M"], e["delta"]
    res = {}
    res["symmetric"] = all(G[i][j] == G[j][i] for i in range(3) for j in range(3))
    # signature (2,1) at id: -G is not PD, and det < 0 with a PD 2x2 leading block is
    # not guaranteed; use: sigma(G) PD and det(G) < 0 at id and some positive vector.
    det = determinant(G)
    res["sigma_positive_definite"] = is_positive_definite(sigma_matrix(G))
    res["det_negative_at_id"] = det.sign() < 0
    # with det < 0 at id the signature is (2,1) or (0,3); exclude negative definite
    res["signature_2_1"] = res["det_negative_at_id"] and not is_positive_definite(
        [[-x for x in r] for r in G])
    res["det_equals_delta"] = det == delta
    res["det_is_delta_times_unit"] = ZS2.is_unit(det / delta) if not delta.is_zero() else False
    res["det_norm"] = delta.norm() == e["det_norm"]
    _, fac = ZS2.factor(delta)
    res["squarefree"] = all(k == 1 for _, k in fac)
    res["unscaled"] = ZS2.is_unit(ZS2.gcd_list([x for r in G for x in r]))
    ok_roots = True
    for r, m in zip(R, M):
        if ip(G, r, r) != m or not ZS2.is_totally_positive(m):
            ok_roots = False
            break
        if not ZS2.is_unit(ZS2.gcd_list(r)):
            ok_roots = False
            break
        Gr = [sum((G[j][i] * r[i] for i in range(3)), QI(0)) for j in range(3)]
        if not all(ZS2.divides(m, 2 * x) for x in Gr):
            ok_roots = False
            break
    res["roots_are_roots"] = ok_roots
    res["roots_non_obtuse"] = all(ip(G, R[i], R[j]).sign() <= 0
                                  for i in range(len(R)) for j in range(i + 1, len(R)))
    base, k = shape_base(e["shape"])
    n = len(R)
    labs = "".join(angle_label(G, R[i], R[(i + 1) % n]) for i in range(n if k == 1 else n - 1))
    res["shape_matches"] = len(base) == n and labs in base + base
    return res
