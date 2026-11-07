#!/bin/bash
# Build the arXiv submission bundle: tech_report/arxiv_submission.tar.gz
#   - main.tex + tables/ + figures/ (pdflatex, no BibTeX: the bibliography is inline)
#   - anc/: ancillary files (verifier, code, certificates, chamber data)
# The bundle is test-compiled in a clean directory before it is packed.
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root
python3 tools/make_paper_tables.py

OUT=tech_report/arxiv
rm -rf "${OUT:?}"
mkdir -p "$OUT/tables" "$OUT/figures" "$OUT/anc/data" "$OUT/anc/tools"
cp tech_report/main.tex "$OUT/"
cp tech_report/tables/*.tex "$OUT/tables/"
cp tech_report/figures/*.tex "$OUT/figures/"

# ancillary files
cp verify_certificate.py "$OUT/anc/"
rsync -a --exclude '__pycache__' reflect "$OUT/anc/"
cp -R data/certs data/isometries "$OUT/anc/data/"
cp data/classification.json data/candidates.json data/mark_full.json data/phase2.json "$OUT/anc/data/"
mkdir -p "$OUT/anc/chambers"
cp diagrams/*.json "$OUT/anc/chambers/"
cp tools/parse_mark.py tools/make_paper_tables.py "$OUT/anc/tools/"
cat > "$OUT/anc/README.txt" <<'EOF'
Ancillary files for "Reflectivity of the quadratic forms -d x0^2 + x1^2 + ... + xn^2
over Z[sqrt2]: a computer-assisted classification".

verify_certificate.py   independent verifier of the non-reflectivity certificates
                         (own exact arithmetic; the facet check re-runs reflect.vinberg)
                         usage: python3 verify_certificate.py data/certs/*.json
reflect/                 the Python package (exact arithmetic in Q(sqrt2), Vinberg's
                         algorithm, rank-3 Vinberg, certificates, Mark-table checks)
data/certs/              all certificates, <label>_n<n>.json
data/isometries/         explicit isometries <-d,1,1> -> entries of Mark's table
data/classification.json per-d summary (reflective range, n0, certificate, matches)
data/candidates.json     the 70 candidates d with the ids of the table entries
data/mark_full.json      full chambers of selected table entries (rank-3 Vinberg)
chambers/                every reflective case: roots, norms, Gram matrix, Coxeter edges,
                         CoxIter invariants
tools/parse_mark.py      parser for the TeX source of arXiv:1512.01133 (the table itself
                         is not redistributed; regenerate it from the arXiv source)

Element notation: "a+b*r2" means a + b*sqrt(2).
EOF

# test-compile in a clean copy
TMP=$(mktemp -d)
cp -R "$OUT"/main.tex "$OUT"/tables "$OUT"/figures "$TMP/"
( cd "$TMP" && for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null; done )
cp "$TMP/main.pdf" tech_report/main.pdf
rm -rf "${TMP:?}"

tar -C "$OUT" -czf tech_report/arxiv_submission.tar.gz main.tex tables figures anc
du -h tech_report/arxiv_submission.tar.gz
echo "OK: tech_report/arxiv_submission.tar.gz"
