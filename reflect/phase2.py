"""Brief 2, Phases 2-3: for each new d reflective at n = 2, go up in n until the first n0
that does not terminate within the caps, then search a certificate at n0.  By the
monotonicity lemma every n > n0 then follows from n0.

    python -m reflect.phase2 [--jobs 4] [labels...]     -> data/phase2.json (+ runs, certs)
"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

from .candidates import OUT as CAND
from .certify import cert_path, certify_one
from .mark import ROOT
from .ring import QI, ZS2
from .run import run_one, run_path

OUT = ROOT / "data" / "phase2.json"
CAPS = dict(delta_max=2e4, max_roots=400, time_limit=1800)


def climb(c):
    label, d = c["label"], QI.parse(c["d"])
    t = time.time()
    faces = {}
    n0 = None
    cert = None
    for n in range(3, 11):
        p = run_path("Zr2", d, n)
        run = json.loads(p.read_text()) if p.exists() else run_one(ZS2, d, n, label, **CAPS)
        if run["status"] == "reflective":
            faces[n] = run["n_roots"]
            continue
        n0 = n
        if cert_path(label, n).exists():
            cert = cert_path(label, n).stem
        else:
            try:
                r = certify_one(label, d, n)
                cert = cert_path(label, n).stem if r and not r.get("failed") else None
            except RuntimeError as e:   # terminated with more roots: record, do not force
                cert = f"recheck: {e}"
        break
    return {"label": label, "d": c["d"], "faces": faces, "n0": n0, "cert_n0": cert,
            "seconds": round(time.time() - t, 1)}


def main(argv):
    jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 4
    labels = {a for a in argv if not a.startswith("--") and not a.isdigit()}
    todo = []
    for c in json.loads(CAND.read_text()):
        if c["in_grid"] or (labels and c["label"] not in labels):
            continue
        p = run_path("Zr2", QI.parse(c["d"]), 2)
        if p.exists() and json.loads(p.read_text())["status"] == "reflective":
            todo.append(c)
    done = json.loads(OUT.read_text()) if OUT.exists() else {}
    with Pool(jobs) as pool:
        for r in pool.imap_unordered(climb, todo):
            done[r["label"]] = r
            OUT.write_text(json.dumps(done, indent=1))
            print(f"{r['label']:10s} d={r['d']:11s} reflective n={sorted(r['faces'])} "
                  f"faces={r['faces']} n0={r['n0']} cert={r['cert_n0']} {r['seconds']}s",
                  flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
