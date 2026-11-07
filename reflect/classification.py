"""Brief 2: the classification of squarefree admissible d over Z[sqrt2], assembled from
data/ only (candidates, n = 2 runs, Phase 2 climbs, certificates, isometry calibration).

    python -m reflect.classification      -> data/classification.json

For each candidate d (from Mark's table) and for the extra conjugate-partner d checked
directly (EXTRA), records: the reflective range {2..N(d)} (N(d) = 1 means none), the first
non-reflective n0, the certificate at n0, faces per reflective n, and the isometry
calibration against Mark's entries.  By the monotonicity lemma every n > n0 is
non-reflective.  d outside the candidate list are non-reflective in every n *provided*
Mark's table is complete; see OPEN_QUESTIONS OQ7 (p17a is reflective but missing).
"""
from __future__ import annotations

import json

from .candidates import OUT as CAND
from .certify import cert_path
from .mark import ROOT
from .ring import QI
from .run import run_path

OUT = ROOT / "data" / "classification.json"
ISO = ROOT / "data" / "isometries"
NS = range(2, 11)
# conjugate partners -sigma(d) of reflective candidates that are not candidates themselves
EXTRA = [{"label": "dm1_12", "d": "-1+12*r2", "why": "partner -sigma(d) of d1_12"}]


def _run(d, n):
    p = run_path("Zr2", d, n)
    return json.loads(p.read_text()) if p.exists() else None


def record(c) -> dict:
    label, d = c["label"], QI.parse(c["d"])
    faces, n0, cert = {}, None, None
    for n in NS:
        r = _run(d, n)
        if r and r["status"] == "reflective":
            faces[n] = r["n_roots"]
            continue
        n0 = n
        if cert_path(label, n).exists():
            cert = cert_path(label, n).stem
        break
    top = max(faces) if faces else 1
    assert sorted(faces) == list(range(2, top + 1)), (label, faces)  # interval (lemma)
    iso = ISO / f"{label}.json"
    match = None
    if iso.exists():
        j = json.loads(iso.read_text())
        match = {"entries_checked": len(j["entries_checked"]), "isometric": len(j["matches"])}
    return {"label": label, "d": c["d"], "N_d": d.norm(), "in_grid": c.get("in_grid", False),
            "mark_entries": len(c.get("entries", [])), "max_reflective_n": top,
            "faces": faces, "n0": n0, "cert_n0": cert, "mark_match": match,
            "note": c.get("why")}


def build():
    C = json.loads(CAND.read_text())
    recs = [record(c) for c in C] + [record(c) for c in EXTRA]
    return recs


def main():
    recs = build()
    OUT.write_text(json.dumps(recs, indent=1))
    refl = [r for r in recs if r["max_reflective_n"] >= 2]
    print(f"{len(recs)} d checked, {len(refl)} reflective at n = 2")
    for r in refl:
        print(f"{r['label']:8s} {r['d']:10s} N(d)={r['max_reflective_n']} faces={r['faces']} "
              f"n0={r['n0']} cert={r['cert_n0']} match={r['mark_match']}")


if __name__ == "__main__":
    main()
