"""Phase 3 driver: search a certificate for every undetermined run in data/runs.

    python -m reflect.certify            # all undetermined Z[sqrt2] runs
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

from .certificates import CERTS, search
from .params import enumerate_params
from .ring import QI, ZS2
from .run import RUNS, run_path
from .vinberg import vinberg


def cert_path(label: str, n: int) -> pathlib.Path:
    return CERTS / f"{label}_n{n}.json"


def certify_one(label: str, d: QI, n: int, extend_roots=(400, 800), log=True):
    path = run_path("Zr2", d, n)
    run = json.loads(path.read_text())
    if run["status"] == "reflective":
        return None
    roots = [tuple(QI.parse(x) for x in r["k"]) for r in run["roots"]]
    attempts = []
    for N in sorted({min(len(roots), x) for x in (20, 40, 80, 160, len(roots))}):
        t = time.time()
        res = search(d, n, roots[:N])
        attempts.append({"prefix": N, "tried": res["tried"], "s": round(time.time() - t, 2)})
        if res["certificate"]:
            cert = res["certificate"]
            cert["label"] = label
            cert["attempts"] = attempts
            CERTS.mkdir(parents=True, exist_ok=True)
            cert_path(label, n).write_text(json.dumps(cert, indent=1))
            return cert
    # extend the Vinberg run and retry
    for mr in extend_roots:
        r = vinberg(ZS2, d, n, delta_max=10 ** 5, max_roots=mr, time_limit=1800)
        if r.status == "reflective":
            raise RuntimeError(f"{label} n={n}: became reflective with more roots -- rerun grid")
        roots = [x.k for x in r.roots]
        t = time.time()
        res = search(d, n, roots)
        attempts.append({"prefix": len(roots), "tried": res["tried"], "s": round(time.time() - t, 2)})
        if res["certificate"]:
            cert = res["certificate"]
            cert["label"] = label
            cert["attempts"] = attempts
            cert_path(label, n).write_text(json.dumps(cert, indent=1))
            return cert
    return {"failed": True, "attempts": attempts}


def main(argv):
    only = set(argv)
    for p in enumerate_params():
        if only and p.label not in only:
            continue
        for n in range(2, 11):
            path = run_path("Zr2", p.d, n)
            if not path.exists():
                print(f"{p.label} n={n}: no run", flush=True)
                continue
            if json.loads(path.read_text())["status"] == "reflective":
                continue
            if cert_path(p.label, n).exists():
                print(f"{p.label:7s} n={n:2d}: certificate exists", flush=True)
                continue
            t = time.time()
            c = certify_one(p.label, p.d, n)
            if c and not c.get("failed"):
                print(f"{p.label:7s} n={n:2d}: CERTIFIED  factors={len(c['factors'])} "
                      f"eigenvalue in {tuple(c['infinite_order']['bracket'])} "
                      f"{time.time() - t:.1f}s", flush=True)
            else:
                print(f"{p.label:7s} n={n:2d}: no certificate found {c}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
