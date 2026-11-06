"""Subprocess bridges to AlVin and CoxIter (independent cross-checks)."""
from __future__ import annotations

import pathlib
import re
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

from .ring import QI, ZS2, ZZ, Ring

ROOT = pathlib.Path(__file__).resolve().parent.parent
ALVIN = ROOT / "tools/src/AlVin/build/alvin"
COXITER = ROOT / "tools/src/CoxIter/build/coxiter"


def tools_available() -> bool:
    return ALVIN.exists() and COXITER.exists()


def _alvin_coeff(x: QI) -> str:
    """AlVin syntax for a + b sqrt2: 'a+bT' (T = sqrt2)."""
    a, b = int(x.a), int(x.b)
    if b == 0:
        return str(a)
    s = "" if a == 0 else str(a)
    if b > 0 and s:
        s += "+"
    s += ("" if b == 1 else "-" if b == -1 else str(b)) + "T"
    return s


def _parse_alvin_entry(s: str) -> QI:
    s = s.replace(" ", "")
    m = re.fullmatch(r"([+-]?\d+)?(?:([+-]?\d*)\*?T\(2\))?", s)
    if not m:
        raise ValueError(f"cannot parse AlVin coefficient {s!r}")
    a = int(m.group(1)) if m.group(1) else 0
    if m.group(2) is None:
        b = 0
    else:
        g = m.group(2)
        b = 1 if g in ("", "+") else -1 if g == "-" else int(g)
    return QI(a, b)


def run_alvin(ring: Ring, d: QI, n: int, timeout: float = 600, extra: List[str] = ()
              ) -> Dict:
    """Run AlVin on -d x0^2 + x1^2 + ... + xn^2.  Returns roots (as QI tuples), whether
    it ended, and the CoxIter graph it wrote."""
    coeffs = [_alvin_coeff(-d)] + ["1"] * n
    args = [str(ALVIN)]
    if ring is ZS2:
        args.append("-k=Q[sqrt 2]")
    args += ["-qf", ",".join(coeffs), *extra]
    with tempfile.TemporaryDirectory() as td:
        (pathlib.Path(td) / "output").mkdir()
        try:
            p = subprocess.run(args, cwd=td, capture_output=True, text=True, timeout=timeout)
            out, timed_out = p.stdout, False
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
            timed_out = True
        graphs = list((pathlib.Path(td) / "output").glob("*.coxiter"))
        graph = graphs[0].read_text() if graphs else None
    roots = []
    for line in out.splitlines():
        mm = re.match(r"\s*e(\d+)\s*=\s*\((.*)\)\s*$", line)
        if mm:
            roots.append(tuple(_parse_alvin_entry(t) for t in mm.group(2).split(",")))
    header = re.search(r"Quadratic form \((\d+),(\d+)\):\s*(.*)", out)
    return {"roots": roots, "ended": "Algorithm ended" in out and not timed_out,
            "timed_out": timed_out, "graph": graph, "stdout": out,
            "form": header.group(3).strip() if header else None,
            "signature": (int(header.group(1)), int(header.group(2))) if header else None}


def run_coxiter(graph: str, timeout: float = 600, extra: List[str] = ()) -> Dict:
    """Run CoxIter on a graph (stdin) with compactness, f-vector, Euler characteristic."""
    p = subprocess.run([str(COXITER), "-c", "-fv", *extra], input=graph,
                       capture_output=True, text=True, timeout=timeout)
    out = p.stdout

    def grab(pat):
        m = re.search(pat, out)
        return m.group(1).strip() if m else None

    fv = grab(r"f-vector: \(([^)]*)\)")
    return {"cocompact": grab(r"Cocompact: (\w+)"),
            "finite_covolume": grab(r"Finite covolume: (\w+)"),
            "f_vector": [int(x) for x in fv.split(",")] if fv else None,
            "vertices_at_infinity": grab(r"Number of vertices at infinity: (\d+)"),
            "euler": grab(r"Euler characteristic: ([-\d/]+)"),
            "dimension": grab(r"Dimension: (\d+)"),
            "stdout": out}
