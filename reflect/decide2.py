"""Decide n = 2 for every candidate d (data/candidates.json).

    python -m reflect.decide2 [--jobs 4] [labels...]

Runs Vinberg at n = 2 (stored in data/runs/), and for a non-terminating run searches for
a non-reflectivity certificate (data/certs/<label>_n2.json).
"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

from .candidates import OUT as CAND
from .certify import cert_path, certify_one
from .ring import QI, ZS2
from .run import run_one, run_path

CAPS = dict(delta_max=1e8, max_roots=400, time_limit=1800)


def decide(c):
    label, d = c["label"], QI.parse(c["d"])
    t = time.time()
    p = run_path("Zr2", d, 2)
    run = json.loads(p.read_text()) if p.exists() else None
    if run is None or (run["status"] != "reflective" and not cert_path(label, 2).exists()
                       and run.get("caps", {}).get("delta_max", 0) < CAPS["delta_max"]):
        run = run_one(ZS2, d, 2, label, **CAPS)
    status = run["status"]
    cert = None
    if status != "reflective":
        if cert_path(label, 2).exists():
            cert = cert_path(label, 2).stem
        else:
            try:
                r = certify_one(label, d, 2)
                if r and not r.get("failed"):
                    cert = cert_path(label, 2).stem
            except RuntimeError as e:          # became reflective with more roots
                status = f"recheck: {e}"
    return {"label": label, "d": c["d"], "status": status, "n_roots": run["n_roots"],
            "cert": cert, "seconds": round(time.time() - t, 1)}


def main(argv):
    jobs = 4
    if "--jobs" in argv:
        jobs = int(argv[argv.index("--jobs") + 1])
    labels = {a for a in argv if not a.startswith("--") and not a.isdigit()}
    cands = [c for c in json.loads(CAND.read_text()) if not labels or c["label"] in labels]
    with Pool(jobs) as pool:
        for r in pool.imap_unordered(decide, cands):
            print(f"{r['label']:10s} d={r['d']:11s} {r['status']:13s} roots={r['n_roots']:3d} "
                  f"cert={r['cert']} {r['seconds']}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
