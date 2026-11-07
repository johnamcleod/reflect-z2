"""Where each reflective cell (d, n) of the classification already appears in the literature.

A cell is attributed to a source only after an explicit check (see OPEN_QUESTIONS OQ6/OQ7
and methods.md):
  - n = 2: an exact isometry <-d,1,1> ~ an entry of Mark's table (data/isometries/);
  - n = 3: Bogachev's rows are diagonal <-D,1,1,1> with D = d * (unit square), checked
    via canonical(); his non-diagonal rows L(6), L(7), L(10), L(11) have determinants
    -(2+3r2) (L(6), L(7)), -(3+4r2), -(11+12r2) up to unit squares -- none equals -d for inert3 or dm2_3;
  - Bugaenko and Guglielmetti: as quoted in Guglielmetti's thesis (Bugaenko 1990 itself
    was not available, OPEN_QUESTIONS OQ1).
Cells not listed here are new as far as these sources go.
"""
from __future__ import annotations

REFS = {
    "Mark18": "A. Mark, The classification of rank 3 reflective hyperbolic lattices over "
              "ℤ[√2], Math. Proc. Camb. Phil. Soc. 164 (2018) 221–257, arXiv:1512.01133",
    "Bug90": "V. O. Bugaenko, Selecta Math. Sov. 9 (1990) (not read: cited as [Bug90] in "
             "Guglielmetti's thesis and in Belolipetsky's survey arXiv:1506.03111)",
    "Bog": "N. Bogachev, From geometry to arithmetic of compact hyperbolic Coxeter "
           "polytopes, arXiv:2003.11944 (candidate list L(1)–L(11))",
    "Gug17": "R. Guglielmetti, PhD thesis, Université de Fribourg, 2017, table of forms "
             "⟨−α,1,…,1⟩ for quadratic Pisot α and the list of reflective forms over ℚ[√2]",
}

# (label, n) -> (source key, detail)
CELLS = {
    ("unit", 2): ("Bug90", "n ≤ 6; also Mark18"), ("unit", 3): ("Bug90", "n ≤ 6; Bogachev L(1)"),
    ("unit", 4): ("Bug90", "n ≤ 6"), ("unit", 5): ("Bug90", "n ≤ 6"),
    ("unit", 6): ("Bug90", "n ≤ 6"),
    ("ram2", 3): ("Bog", "L(5)"),
    ("p7a", 3): ("Bog", "L(3) = ⟨−5−4√2,1,1,1⟩"),
    ("p7a", 4): ("Gug17", "⟨1−2Θ₂,1,1,1,1⟩ and ⟨−5−4Θ₂,…⟩"),
    ("p7a", 5): ("Gug17", "⟨−5−4Θ₂,1,…,1⟩, n ≤ 5"),
    ("p7b", 3): ("Bog", "L(2) (and L(4) ≅ L(2))"),
    ("p23a", 3): ("Bog", "L(8) = ⟨−7−6√2,1,1,1⟩, 7+6√2 = (−3+4√2)(1+√2)²"),
}
# Reflective cells found in none of the sources, with the reason (tests check that this is
# exactly the set of reflective cells without a source).
NEW = {
    ("p17a", 2): "no isometric entry in Mark's table; its only (842)² hexagon is p17b (OQ7)",
    ("ram2", 4): "Bogachev covers n = 3 only; no source for ⟨−√2,1,…,1⟩ in higher n",
    ("ram2", 5): "as n = 4",
    ("p7b", 4): "Bogachev covers n = 3 only (L(2)); not in Guglielmetti's list",
    ("inert3", 3): "no row of Bogachev's L(1)–L(11) has determinant class −(3+3√2)",
    ("dm2_3", 3): "no row of Bogachev's L(1)–L(11) has determinant class −(−2+3√2)",
}
# n = 2 cells come from Mark's table via data/isometries/ (all but p17a, OQ7).
# Upper ends (first non-reflective n0) already in the literature:
N0 = {"unit": ("Bug90", "reflective ⇔ n ≤ 6"), "p7a": ("Gug17", "table entry n ≤ 5")}


def source(label: str, n: int, mark_isometric: int) -> tuple[str, str] | None:
    if n == 2 and (label, n) not in CELLS:
        return ("Mark18", f"{mark_isometric} isometric entries") if mark_isometric else None
    return CELLS.get((label, n))
