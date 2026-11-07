"""Write the LaTeX tables of the technical report (tech_report/tables/*.tex) from data/.

    python3 tools/make_paper_tables.py

Every number in the tables comes from the stored runs, certificates and isometries, so the
report can be regenerated after any change to data/.
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from reflect.candidates import OUT as CAND          # noqa: E402
from reflect.classification import EXTRA, build    # noqa: E402
from reflect.literature import CELLS, N0, NEW, source  # noqa: E402
from reflect.ring import QI                         # noqa: E402
from reflect.run import run_path                    # noqa: E402

OUT = ROOT / "tech_report" / "tables"
CERTS = ROOT / "data" / "certs"
CITE = {"Mark18": "Mar18", "Bug90": "Bug90", "Bog": "Bog20", "Gug17": "Gug17"}


def tex_num(s: str) -> str:
    """'-1+3*r2' -> '-1+3\\sqrt2' (math mode content)."""
    q = QI.parse(s)
    a, b = q.a, q.b
    if b == 0:
        return f"{a}"
    bs = "" if abs(b) == 1 else f"{abs(b)}"
    if a == 0:
        return f"{'-' if b < 0 else ''}{bs}\\sqrt2"
    return f"{a}{'-' if b < 0 else '+'}{bs}\\sqrt2"


def factor_tex(c) -> str:
    if c.get("primes") is None:          # extra d (not a candidate): factor it here
        from reflect.ring import ZS2
        _, fac = ZS2.factor(QI.parse(c["d"]))
        c = {"primes": [{"pi": str(q)} for q, _ in fac]}
    if not c["primes"]:
        return r"\text{unit}"
    return r"\cdot".join(f"({tex_num(p['pi'])})" for p in c["primes"])


def polygon_cycle(run) -> str:
    """Angle labels k (angle pi/k) around the n = 2 polygon, as a string like 842842."""
    N = run["n_roots"]
    lab = {}
    for i, j, k in run["coxeter_edges"]:
        lab[(i, j)] = lab[(j, i)] = k
    adj = {i: [j for j in range(N) if j != i and lab.get((i, j), 2) != "dotted"]
           for i in range(N)}
    assert all(len(v) == 2 for v in adj.values()), adj
    cyc, prev, cur = [0], None, 0
    while True:
        nxt = [j for j in adj[cur] if j != prev][0]
        if nxt == 0:
            break
        cyc.append(nxt)
        prev, cur = cur, nxt
    ks = [str(lab.get((cyc[i], cyc[(i + 1) % N]), 2)) for i in range(N)]
    rots = [ks[i:] + ks[:i] for i in range(N)]
    rots += [r[::-1] for r in rots]
    return "".join(min(rots))


def run_of(d, n):
    p = run_path("Zr2", QI.parse(d), n)
    return json.loads(p.read_text()) if p.exists() else None


def reflective_table(recs, cands):
    rows = []
    refl = sorted([r for r in recs if r["max_reflective_n"] >= 2],
                  key=lambda r: (-r["max_reflective_n"], abs(r["N_d"]), r["label"]))
    for r in refl:
        c = cands.get(r["label"], {})
        iso = r["mark_match"]["isometric"] if r["mark_match"] else 0
        srcs = []
        for n in sorted(r["faces"]):
            s = source(r["label"], n, iso)
            srcs.append(f"\\cite{{{CITE[s[0]]}}}" if s else r"\textbf{new}")
        # compress identical consecutive sources
        faces = ", ".join(str(r["faces"][n]) for n in sorted(r["faces"]))
        n0 = f"{r['n0']}" + (f"~\\cite{{{CITE[N0[r['label']][0]]}}}" if r["label"] in N0 else "")
        rows.append(f"\\texttt{{{r['label'].replace('_', chr(92) + '_')}}} & ${tex_num(r['d'])}$ & "
                    f"${r['N_d']}$ & ${factor_tex(c)}$ & $2$--${r['max_reflective_n']}$ & {faces} & "
                    f"{n0} & {', '.join(srcs)} \\\\")
    (OUT / "reflective.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")
    return refl


def invariants_table(refl):
    rows = []
    for r in refl:
        for n in sorted(r["faces"]):
            run = run_of(r["d"], n)
            cx = run.get("coxiter_check") or {}
            fv = cx.get("f_vector")
            fvs = "(" + ", ".join(str(x) for x in fv) + ")" if fv else "--"
            cyc = polygon_cycle(run) if n == 2 else ""
            new = r"$\star$" if (r["label"], n) in NEW else ""
            rows.append(f"\\texttt{{{r['label'].replace('_', chr(92) + '_')}}} & {n} & "
                        f"{run['n_roots']} & {fvs} & ${cx.get('euler', '--')}$ & "
                        f"\\texttt{{{cyc}}} & {new} \\\\")
            assert cx.get("cocompact") == "yes", (r["label"], n)
    (OUT / "invariants.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")


def candidates_table(recs, cands):
    rows = []
    for r in sorted(recs, key=lambda r: (abs(r["N_d"]), r["label"])):
        c = cands.get(r["label"], {"entries": [], "primes": None, "d": r["d"]})
        run = run_of(r["d"], 2)
        if r["max_reflective_n"] >= 2:
            st = f"R ({run['n_roots']})"
        else:
            st = f"N ({run['n_roots']} roots)"
        m = r["mark_match"]
        iso = f"{m['isometric']}/{m['entries_checked']}" if m else "--"
        extra = r"$^\dagger$" if r["label"] in {e["label"] for e in EXTRA} else ""
        rows.append(f"\\texttt{{{r['label'].replace('_', chr(92) + '_')}}}{extra} & "
                    f"${tex_num(r['d'])}$ & ${r['N_d']}$ & ${factor_tex(c)}$ & "
                    f"{len(c['entries'])} & {st} & {iso} \\\\")
    (OUT / "candidates.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")


def certificates_table(recs):
    rows = []
    for r in sorted(recs, key=lambda r: (abs(r["N_d"]), r["label"])):
        p = CERTS / f"{r['cert_n0']}.json"
        j = json.loads(p.read_text())
        io = j["infinite_order"]
        br = io.get("bracket")
        brs = f"$[{br[0]}, {br[1]}]$" if br else "--"
        nf = len(j.get("factors", [])) or 1
        rows.append(f"\\texttt{{{r['label'].replace('_', chr(92) + '_')}}} & {r['n0']} & "
                    f"{j.get('n_roots_used', '--')} & {nf} & {brs} \\\\")
    (OUT / "certificates.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")


def numbers(recs, cands):
    refl = [r for r in recs if r["max_reflective_n"] >= 2]
    certs = sorted(CERTS.glob("*.json"))
    macros = {
        "NCand": len(cands), "NRefl": len(refl), "NNon": len(recs) - len(refl),
        "NCerts": len(certs), "NNew": len(NEW),
        "NEntries": len(__import__("reflect.mark", fromlist=["load"]).load()),
    }
    (OUT / "numbers.tex").write_text(
        "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cands = {c["label"]: c for c in json.loads(CAND.read_text())}
    recs = build()
    refl = reflective_table(recs, cands)
    invariants_table(refl)
    candidates_table(recs, cands)
    certificates_table(recs)
    numbers(recs, cands)
    print("wrote", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
