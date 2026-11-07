"""Brief 2 classification (reflect/classification.py, report section) and OQ7."""
import json
import pathlib

import pytest

from reflect.candidates import OUT as CAND
from reflect.certify import cert_path
from reflect.classification import EXTRA, build
from reflect.facet import same_line
from reflect.mark_full import interior_point
from reflect.ring import QI
from reflect.run import run_path
from reflect.vinberg3 import vinberg3

ROOT = pathlib.Path(__file__).resolve().parent.parent
RECS = build()


def test_every_candidate_and_extra_is_recorded():
    labels = {c["label"] for c in json.loads(CAND.read_text())} | {e["label"] for e in EXTRA}
    assert {r["label"] for r in RECS} == labels


@pytest.mark.parametrize("r", RECS, ids=lambda r: r["label"])
def test_decided_with_certificate_at_n0(r):
    # reflective range is an interval {2..N} (lemma) and n0 = N+1 carries a certificate
    top = r["max_reflective_n"]
    assert sorted(r["faces"]) == list(range(2, top + 1))
    assert r["n0"] == top + 1
    assert r["cert_n0"] and cert_path(r["label"], r["n0"]).exists()


def test_reflective_candidates_match_mark_except_p17a():
    refl = [r for r in RECS if r["max_reflective_n"] >= 2]
    missing = [r["label"] for r in refl if not r["mark_match"] or r["mark_match"]["isometric"] == 0]
    assert missing == ["p17a"]


def test_oq7_p17a_hexagon_by_independent_vinberg():
    d = QI.parse("-1+3*r2")
    run = json.loads(run_path("Zr2", d, 2).read_text())
    assert run["status"] == "reflective" and run["n_roots"] == 6
    G = [[-d, QI(0), QI(0)], [QI(0), QI(1), QI(0)], [QI(0), QI(0), QI(1)]]
    R = [[QI.parse(x) for x in t["k"]] for t in run["roots"]]
    e = {"G": G, "R": [R[i] for i in (5, 0, 1, 3, 2, 4)], "shape": "842842", "id": -1}
    res = vinberg3(G, interior_point(e, QI(0)), Delta_max=1e7, max_roots=60)
    assert res["status"] == "reflective" and len(res["roots"]) == 6
    assert all(any(same_line(a, b) for b in res["roots"]) for a in R)


def test_oq7_only_842_hexagon_in_table_is_p17b():
    from reflect.mark import load
    from reflect.mark_classes import canon_shape
    hits = [e["id"] for e in load() if canon_shape(e["shape"]) == canon_shape("(842)^2")]
    assert hits == [173]
    j = json.loads((ROOT / "data/isometries/p17b.json").read_text())
    assert [m["entry"] for m in j["matches"]] == [173]
    j = json.loads((ROOT / "data/isometries/p17a.json").read_text())
    assert j["matches"] == [] and 173 not in j["entries_checked"]
