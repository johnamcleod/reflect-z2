"""Run Vinberg for one (ring, d, n) and store the artefact in data/runs/.

    python -m reflect.run Z 5 2..9            # calibration over Q
    python -m reflect.run Zr2 1+r2 2..7       # over Z[sqrt2]
    python -m reflect.run Zr2 all 2..10       # the whole parameter grid
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

from .bridges import run_alvin, run_coxiter, tools_available
from .params import enumerate_params
from .ring import QI, ZS2, ZZ
from .vinberg import vinberg

ROOT = pathlib.Path(__file__).resolve().parent.parent
RUNS = ROOT / "data" / "runs"


def run_path(ring_name: str, d: QI, n: int) -> pathlib.Path:
    tag = re.sub(r"[^0-9a-zA-Z+\-]", "_", str(d)).replace("*", "")
    return RUNS / f"{ring_name}_{tag}_n{n}.json"


def run_one(ring, d: QI, n: int, label: str = "", crosscheck: bool = True, **caps) -> dict:
    res = vinberg(ring, d, n, **caps)
    out = res.to_json()
    out["label"] = label
    if crosscheck and tools_available():
        if res.status == "reflective":
            cx = run_coxiter(res.coxiter)
            out["coxiter_check"] = {k: v for k, v in cx.items() if k != "stdout"}
    rn = "Z" if ring is ZZ else "Zr2"
    path = run_path(rn, d, n)
    RUNS.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    return out


def parse_range(s: str):
    if ".." in s:
        a, b = s.split("..")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ring", choices=["Z", "Zr2"])
    ap.add_argument("d")
    ap.add_argument("ns")
    ap.add_argument("--delta-max", type=float, default=2000.0)
    ap.add_argument("--max-roots", type=int, default=150)
    ap.add_argument("--time-limit", type=float, default=900.0)
    ap.add_argument("--stop-after-first-fail", action="store_true",
                    help="skip higher n once a run is undetermined (saves time)")
    a = ap.parse_args(argv)
    ring = ZZ if a.ring == "Z" else ZS2
    if a.d == "all":
        targets = [(p.d, p.label) for p in enumerate_params()]
    else:
        targets = [(QI.parse(a.d), "")]
    for d, label in targets:
        for n in parse_range(a.ns):
            out = run_one(ring, d, n, label, delta_max=a.delta_max, max_roots=a.max_roots,
                          time_limit=a.time_limit)
            cx = out.get("coxiter_check", {})
            print(f"{label or a.ring:8s} d={d} n={n:2d}: {out['status']:12s} roots={out['n_roots']:3d} "
                  f"batches={out['batches']:5d} delta<={out['max_delta_processed']:.1f} "
                  f"stop={out['caps']['stopped_by']} {out['seconds']:.1f}s "
                  f"coxiter: cocompact={cx.get('cocompact')} fv={cx.get('finite_covolume')} "
                  f"ideal={cx.get('vertices_at_infinity')} chi={cx.get('euler')}", flush=True)
            if a.stop_after_first_fail and out["status"] != "reflective":
                break


if __name__ == "__main__":
    main()
