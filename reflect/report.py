"""Phase 4: results.csv, diagrams/ and report.md, regenerated from data/ only.

    python -m reflect.report
"""
from __future__ import annotations

import csv
import json
import pathlib
import subprocess

from .params import enumerate_params
from .ring import QI
from .run import run_path
from .vinberg import form

ROOT = pathlib.Path(__file__).resolve().parent.parent
DIAG = ROOT / "diagrams"
CERTS = ROOT / "data" / "certs"
SUPER = ROOT / "data" / "superlattices"
NS = range(2, 11)
Q_PATTERN = {1: 19, 2: 14, 3: 13, 5: 8, 7: 3, 11: 4, 17: 3}


def load_run(d: QI, n: int):
    p = run_path("Zr2", d, n)
    return json.loads(p.read_text()) if p.exists() else None


def cell(label: str, d: QI, n: int):
    r = load_run(d, n)
    cert = CERTS / f"{label}_n{n}.json"
    if r and r["status"] == "reflective":
        return "reflective", r["n_roots"], r
    if cert.exists():
        return "non-reflective", cert.stem, r
    lower = [m for m in range(2, n) if (CERTS / f"{label}_n{m}.json").exists()]
    if lower:  # n -> n+1 lemma (methods.md, "Extension to all n")
        return "non-reflective-lemma", f"{label}_n{min(lower)}", r
    return "undetermined", None, r


# ---- diagrams -------------------------------------------------------------------------
def write_diagram(label: str, d: QI, n: int, run: dict):
    DIAG.mkdir(exist_ok=True)
    stem = DIAG / f"{label}_n{n}"
    roots = [[QI.parse(x) for x in r["k"]] for r in run["roots"]]
    gram = [[str(form(d, a, b)) for b in roots] for a in roots]
    edges = run["coxeter_edges"]
    (stem.with_suffix(".json")).write_text(json.dumps({
        "label": label, "d": str(d), "n": n, "form": f"-({d}) x0^2 + x1^2 + ... + x{n}^2",
        "roots": [r["k"] for r in run["roots"]], "norms": [r["m"] for r in run["roots"]],
        "gram": gram, "coxeter_edges": edges,
        "coxiter": run.get("coxiter_check")}, indent=1))
    # TikZ (circular layout; Coxeter conventions: no edge = pi/2, plain = pi/3,
    # labelled 4, 6, 8; dashed = ultraparallel)
    import math
    N = len(roots)
    lines = ["\\begin{tikzpicture}[scale=2.2,v/.style={circle,draw,fill=white,inner sep=1.2pt,font=\\tiny}]"]
    for i in range(N):
        a = 2 * math.pi * i / N
        lines.append(f"  \\node[v] ({i}) at ({math.cos(a):.3f},{math.sin(a):.3f}) {{{i + 1}}};")
    for i, j, lab in edges:
        if lab == "dotted":
            lines.append(f"  \\draw[dashed] ({i}) -- ({j});")
        elif lab == 3:
            lines.append(f"  \\draw ({i}) -- ({j});")
        else:
            lines.append(f"  \\draw ({i}) -- node[fill=white,inner sep=0.5pt,font=\\tiny] {{{lab}}} ({j});")
    lines.append("\\end{tikzpicture}")
    stem.with_suffix(".tex").write_text("\n".join(lines) + "\n")
    # Graphviz -> PNG
    dot = ["graph G {", "  node [shape=circle, fontsize=10, width=0.3, fixedsize=true];",
           '  overlap=false; splines=true;']
    for i in range(N):
        dot.append(f'  {i} [label="{i + 1}"];')
    for i, j, lab in edges:
        if lab == "dotted":
            dot.append(f"  {i} -- {j} [style=dashed, color=gray50];")
        elif lab == 3:
            dot.append(f"  {i} -- {j};")
        else:
            dot.append(f'  {i} -- {j} [label="{lab}", fontsize=9];')
    dot.append("}")
    stem.with_suffix(".dot").write_text("\n".join(dot) + "\n")
    try:
        subprocess.run(["neato", "-Tpng", "-o", str(stem.with_suffix(".png")),
                        str(stem.with_suffix(".dot"))], check=True, capture_output=True)
    except Exception as e:  # noqa: BLE001
        print("graphviz failed for", stem, e)


# ---- report -----------------------------------------------------------------------------
def main():
    params = enumerate_params()
    rows = []
    summary = []
    for p in params:
        row = {"label": p.label, "kind": p.kind, "p": p.p, "ideal_norm": p.ideal_norm,
               "d": str(p.d), "N(d)": p.d.norm()}
        maxn = 1
        cells = {}
        for n in NS:
            status, info, run = cell(p.label, p.d, n)
            if status == "reflective":
                row[f"n{n}"] = f"reflective ({info} faces)"
                write_diagram(p.label, p.d, n, run)
                if maxn == n - 1:
                    maxn = n
            elif status == "non-reflective":
                row[f"n{n}"] = f"non-reflective (cert {info})"
            elif status == "non-reflective-lemma":
                row[f"n{n}"] = f"non-reflective (lemma from {info})"
            else:
                row[f"n{n}"] = "undetermined"
            cells[n] = (status, info, run)
        row["max_reflective_n"] = maxn if maxn >= 2 else "none"
        rows.append(row)
        summary.append((p, maxn, cells))
    with open(ROOT / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    write_report(summary)


def _sup_line(label, n):
    f = SUPER / f"{label}_n{n}.json"
    if not f.exists():
        return None
    s = json.loads(f.read_text())
    out = s["outcome"]
    best = [e for e in s.get("decompositions", []) if e["vinberg"]["status"] == "reflective"]
    if best:
        e = best[0]
        out += f": [L':L] has norm {e['index_norm']}, {e['vinberg']['n_roots']} faces"
    return out


def write_report(summary):
    L = []
    L.append("# Reflectivity of −d·x₀² + x₁² + … + xₙ² over ℤ[√2]\n")
    L.append("Generated by `python -m reflect.report` from `data/` (do not edit by hand).\n")
    L.append("## Summary table\n")
    L.append("R k = reflective with k faces (Coxeter diagram in `diagrams/`); "
             "N = non-reflective, with a certificate in `data/certs/` that "
             "`verify_certificate.py` re-checks; N* = non-reflective by the n → n+1 lemma from a "
             "certified lower n (no direct certificate found); U = undetermined.\n")
    hdr = "| param | type | N(𝔭) | d | " + " | ".join(f"n={n}" for n in NS) + " | max n |"
    L.append(hdr)
    L.append("|" + "---|" * (5 + len(NS)))
    for p, maxn, cells in summary:
        cs = []
        for n in NS:
            st, info, _ = cells[n]
            cs.append(f"R {info}" if st == "reflective" else "N" if st == "non-reflective"
                      else "N*" if st == "non-reflective-lemma" else "U")
        L.append(f"| {p.label} | {p.kind} | {p.ideal_norm} | {p.d} | " + " | ".join(cs) +
                 f" | {maxn if maxn >= 2 else '–'} |")
    L.append("")
    L.append("## Maximal reflective dimension against |N(d)|\n")
    L.append("| |N(d)| | params (type) | d ≈ | σ(d) ≈ | max reflective n |")
    L.append("|---|---|---|---|---|")
    for p, maxn, _ in sorted(summary, key=lambda t: (abs(t[0].d.norm()), t[0].label)):
        L.append(f"| {abs(p.d.norm())} | {p.label} ({p.kind}) | {p.d.f():.3f} | {p.d.fs():.3f} "
                 f"| {maxn if maxn >= 2 else 'none'} |")
    L.append("")
    L.append("### The two primes above a split p\n")
    L.append("For split p = 𝔭𝔭′ the two canonical generators differ by which embedding is "
             "small: the 'a' generator has the smaller value d at the time-like embedding id.\n")
    L.append("| p | a: d, max n | b: d, max n |")
    L.append("|---|---|---|")
    by = {p.label: (p, m) for p, m, _ in summary}
    for q in sorted({p.p for p, _, _ in summary if p.kind == "split"}):
        a, b = by.get(f"p{q}a"), by.get(f"p{q}b")
        if a and b:
            f = lambda t: f"{t[0].d.f():.2f}, {t[1] if t[1] >= 2 else 'none'}"
            L.append(f"| {q} | {f(a)} | {f(b)} |")
    L.append("")
    L.append("Over ℚ the analogous pattern is d = 1, 2, 3, 5, 7, 11, 17 → n = "
             + ", ".join(str(v) for v in Q_PATTERN.values()) + ".\n")
    L.append("## Reflective super-lattice check (non-reflective cases)\n")
    any_sup = False
    for p, maxn, cells in summary:
        for n in NS:
            line = _sup_line(p.label, n)
            if line:
                any_sup = True
                L.append(f"- {p.label}, n = {n}: {line}")
    if not any_sup:
        L.append("- (not yet run)")
    L.append("")
    (ROOT / "report.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
