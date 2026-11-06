"""Run the 7b super-lattice analysis for each certified non-reflective case.

    python -m reflect.super_driver [labels...] [--all-n]
By default only the smallest non-reflective dimension n0 of each parameter is analysed.
"""
from __future__ import annotations

import json
import sys
import time

from .certify import cert_path
from .params import enumerate_params
from .ring import QI
from .run import run_path
from .superlattice import OUT, analyse


def main(argv):
    all_n = "--all-n" in argv
    labels = {a for a in argv if not a.startswith("--")}
    OUT.mkdir(parents=True, exist_ok=True)
    for p in enumerate_params():
        if labels and p.label not in labels:
            continue
        ns = [n for n in range(2, 11) if cert_path(p.label, n).exists()]
        if not ns:
            continue
        for n in (ns if all_n else ns[:1]):
            out = OUT / f"{p.label}_n{n}.json"
            if out.exists():
                print(f"{p.label} n={n}: exists", flush=True)
                continue
            run = json.loads(run_path("Zr2", p.d, n).read_text())
            roots = [tuple(QI.parse(x) for x in r["k"]) for r in run["roots"]][:80]
            t = time.time()
            res = analyse(p.label, p.d, n, roots,
                          log=lambda s: print(s, flush=True),
                          vinberg_kw=dict(delta_max=3000, max_roots=150, time_limit=400))
            res["seconds"] = time.time() - t
            out.write_text(json.dumps(res, indent=1))
            print(f"{p.label:7s} n={n}: {res['outcome']} (symmetries {res['n_symmetries']}, "
                  f"pure translations {res['n_pure_translations']}) {res['seconds']:.0f}s",
                  flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
