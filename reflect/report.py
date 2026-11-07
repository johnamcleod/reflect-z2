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


def _observations(summary):
    by = {p.label: (p, m) for p, m, _ in summary}
    pairs, a_ge, a_gt, b_gt = [], 0, 0, []
    for q in sorted({p.p for p, _, _ in summary if p.kind == "split"}):
        a, b = by.get(f"p{q}a"), by.get(f"p{q}b")
        if a and b:
            ma, mb = a[1], b[1]
            pairs.append(q)
            a_ge += ma >= mb
            a_gt += ma > mb
            if mb > ma:
                b_gt.append(q)
    top = sorted(summary, key=lambda t: -t[1])
    none = [p.label for p, m, _ in summary if m < 2]
    big = [f"{p.label} ({m})" for p, m, _ in top if m >= 3]
    out = [
        f"- **The maximal reflective n falls quickly with |N(d)|.** Only {', '.join(big)} reach n ≥ 3. "
        f"{', '.join(none) or 'None'} are non-reflective already at n = 2.",
        f"- **The two primes above a split p behave differently.** The generator that is small at the "
        f"time-like embedding id ('a') does at least as well in {a_ge} of {len(pairs)} pairs, "
        f"and strictly better in {a_gt} of them. "
        + (f"The exception is p = {', '.join(map(str, b_gt))}. " if b_gt else "")
        + "So the size of |σ(d)| matters as well as |N(d)|.",
        "- **Split, inert or ramified makes no visible difference beyond the norm.** The ramified prime "
        "(√2, norm 2) and the inert prime 3 (norm 9) fit the decline with |N(d)|.",
        f"- **The maxima are far below the ℚ pattern.** The largest here is "
        f"{top[0][1]} (for {top[0][0].label}), against 19 for d = 1 over ℚ. This is expected: our "
        "groups are co-compact, whereas the ℚ groups have cusps.",
        "- **Non-reflectivity extends to every n.** It passes from n to n+1 "
        "(methods.md, 'Extension to all n'). Every non-reflective cell here therefore holds for all larger n as well.",
    ]
    return "\n".join(out) + "\n"


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
    L.append("**Scope of the super-lattice check.** These results say *not found*, not *ruled out*.")
    L.append("")
    L.append("- **Which translations were tried.** Only pure translations h that are words of length ≤ 2 in at most 10 vertex symmetries.")
    L.append("- **Which super-lattices were built.** Up to 3 distinct L′ = L + R_a(L) per case.")
    L.append("- **Caps on each Vinberg run on L′.** δ ≤ 3000, about 150 roots, 400 s.")
    L.append("- **Effectively untested candidates.** Some L′ of large index stopped at the time limit before any root was found.")
    L.append("- **Whole cases that timed out.** A case marked 'exceeded time budget' was cut off by a 1500 s hard limit.")
    L.append("")
    L.append("## Observations\n")
    L.append(_observations(summary))
    L.append("")
    L += classification_section()
    (ROOT / "report.md").write_text("\n".join(L) + "\n")


def classification_section():
    """Brief 2 (additional_notes.md §7): all squarefree admissible d."""
    from .classification import build
    from .literature import N0, NEW, REFS, source
    recs = build()
    for r in recs:                      # diagrams for every reflective cell
        for n in r["faces"]:
            write_diagram(r["label"], QI.parse(r["d"]), n, load_run(QI.parse(r["d"]), n))
    with open(ROOT / "results_classification.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "d", "N(d)", "in_grid", "mark_entries", "max_reflective_n",
                    "faces", "n0", "cert_n0", "mark_isometric"])
        for r in recs:
            w.writerow([r["label"], r["d"], r["N_d"], r["in_grid"], r["mark_entries"],
                        r["max_reflective_n"] if r["max_reflective_n"] >= 2 else "none",
                        " ".join(f"n{n}:{k}" for n, k in sorted(r["faces"].items())),
                        r["n0"], r["cert_n0"],
                        r["mark_match"]["isometric"] if r["mark_match"] else ""])
    refl = sorted([r for r in recs if r["max_reflective_n"] >= 2],
                  key=lambda r: (-r["max_reflective_n"], abs(r["N_d"]), r["label"]))
    non = [r for r in recs if r["max_reflective_n"] < 2]
    L = ["## Classification over ℤ[√2] (brief 2)\n",
         "**Method.** The method rests on three facts:\n",
         "- **Monotonicity lemma** (docs/monotonicity.md): reflectivity at n implies reflectivity "
         "at every n′ < n. So the reflective n for each d form an interval {2, …, N(d)}.",
         "- **Mark's table** (arXiv:1512.01133): ⟨−d,1,1⟩ is strongly squarefree, so if it is "
         "reflective it should appear in the table. The 522 table entries give 70 candidate "
         "values of d (data/candidates.json).",
         "- **Direct decision.** Every candidate is decided by our own runs: Vinberg up to N(d), "
         "then a certificate at n₀ = N(d)+1, which `verify_certificate.py` re-checks. Every "
         "n > n₀ is then non-reflective by the lemma.\n",
         f"**Result.** {len(refl)} values of d are reflective, and {len(non)} are non-reflective "
         "in every n ≥ 2. The non-reflective ones are the remaining candidates plus the extra "
         "conjugate partner dm1_12; each is certified at n = 2.\n",
         "**Sources.** Each reflective cell (d, n) is annotated with the earliest paper we "
         "found that contains it ([Mark18], [Bug90], [Bog], [Gug17]; list below). "
         "**🆕 new** marks a cell that appears in none of them; "
         "the attributions are checked as described in `reflect/literature.py`.\n",
         "| param | d | norm of d | reflective n | faces | n₀ (cert) | Mark entries isometric "
         "| sources per reflective n |",
         "|---|---|---|---|---|---|---|---|"]
    new_cells = []
    for r in refl:
        m = r["mark_match"]
        iso = m["isometric"] if m else 0
        mm = f"{iso} of {m['entries_checked']}" if m else "–"
        if m and iso == 0:
            mm += " **(missing, OQ7)**"
        srcs = []
        for n in sorted(r["faces"]):
            s = source(r["label"], n, iso)
            if s:
                srcs.append(f"n={n}: [{s[0]}]")
            else:
                srcs.append(f"**n={n}: 🆕 new**")
                new_cells.append((r, n))
        n0 = f"{r['n0']} ({r['cert_n0']})"
        if r["label"] in N0:
            n0 += f" [{N0[r['label']][0]}]"
        L.append(f"| {r['label']} | {r['d']} | {r['N_d']} | 2–{r['max_reflective_n']} | "
                 + ", ".join(f"{k}" for _, k in sorted(r["faces"].items()))
                 + f" | {n0} | {mm} | " + "; ".join(srcs) + " |")
    L.append("")
    L.append("### New reflective cases\n")
    L.append("These cells are reflective by our runs, with the Coxeter diagram in `diagrams/` "
             "and a CoxIter check, and appear in none of the sources below:\n")
    for r, n in new_cells:
        L.append(f"- **{r['label']}, d = {r['d']}, n = {n}**: {r['faces'][n]} faces; "
                 f"{NEW.get((r['label'], n), 'UNEXPLAINED')}.")
    L.append("")
    L.append("**Non-reflectivity above the reflective range.** The first non-reflective n₀ was "
             "already known only for unit (n₀ = 7, [Bug90]) and p7a (n₀ = 6, consistent with "
             "the range n ≤ 5 in [Gug17]). For the 43 candidates non-reflective at n = 2, "
             "non-reflectivity also follows from [Mark18] (no isometric entry). For all other "
             "d the certificate at n₀ is new.\n")
    L.append("**References.**\n")
    for k, v in REFS.items():
        L.append(f"- [{k}] {v}")
    L.append("")
    L.append("**Non-reflective for all n** (certificate at n = 2): "
             + ", ".join(f"{r['label']} ({r['d']})" for r in non) + ".\n")
    L.append("**Caveats.**\n")
    L.append("- **Mark's table is incomplete (OQ7).** ⟨−(−1+3√2),1,1⟩ (p17a) is reflective: its "
             "chamber is the hexagon (π/8, π/4, π/2)², confirmed independently by "
             "`reflect.vinberg3`. Yet none of the 8 table entries with its determinant is "
             "isometric to it. The table's only hexagon of this shape (entry 173) is p17b. "
             "So the step \"d outside the candidate list ⇒ non-reflective\" is conditional on "
             "the table having no further omissions.")
    L.append("- **Conjugate-partner check.** For every reflective d we checked its conjugate "
             "partner −σ(d) directly. All but one were already candidates. The exception, "
             "−1+12√2 (partner of the 24-gon d = 1+12√2), is certified non-reflective at n = 2. "
             "Reflectivity is not symmetric between partners (for example, p31b is reflective "
             "and p31a is not), so this check cannot rule out further omissions.")
    L.append("- **Duplicates in the table (OQ4).** The 522 entries fall into at least 144 isometry "
             "classes. Explicit isometries merge at least 144 entries, so there are at most 378 "
             "classes, not 432. That check is incomplete and does not affect the result.")
    L.append("")
    return L


if __name__ == "__main__":
    main()
