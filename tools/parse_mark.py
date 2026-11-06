#!/usr/bin/env python3
"""Parse the appendix tables of Mark, arXiv:1512.01133 (TeX source) into JSON.

    python3 -I tools/parse_mark.py <dir with the unpacked arXiv source> data/mark_table.json

Pure text parsing, no evaluation of the input.  Elements a + b sqrt2 are stored as
[a, b] integer pairs.  Consistency checks are done separately (tests/test_mark_table.py).
"""
import json
import re
import sys
from pathlib import Path

FILES = ["triangles", "quadrilaterals", "pentagons", "hexagons", "septagons", "octagons",
         "decagons", "undecagons", "dodecagons", "tetradecagons", "hexadecagons",
         "octadecagons", "icosagons", "24-gons"]

NUM = r"-?\d+"


def parse_elem(s):
    """'-24 - 17w', '-3w', '1 - w', '7' -> [a, b]."""
    s = s.replace(" ", "")
    terms = re.findall(r"[+-]?[^+-]+", s)
    if not s or "".join(terms) != s:
        raise ValueError(f"cannot parse element {s!r}")
    a = b = 0
    for t in terms:
        if t.endswith("w"):
            c = t[:-1]
            b += 1 if c in ("", "+") else -1 if c == "-" else int(c)
        else:
            if not re.fullmatch(r"[+-]?\d+", t):
                raise ValueError(f"cannot parse element {s!r}")
            a += int(t)
    return [a, b]


def rows_of(block, start_pat, ncols):
    out = []
    for line in block.splitlines():
        line = line.strip()
        if not re.match(start_pat, line):
            continue
        body = re.sub(r"\\\\\*?\s*$", "", line)
        toks = [t.strip() for t in body.split("&")]
        if toks and toks[0].startswith("r_"):
            toks = toks[1:]
        if len(toks) != ncols or not all(re.fullmatch(NUM, t) for t in toks):
            raise ValueError(f"bad row {line!r}")
        out.append([int(t) for t in toks])
    return out


def parse_block(block, src, idx):
    det = re.search(r"det = \$(.*?)\$\s*&\s*det norm = \$(.*?)\$", block)
    shape = re.search(r"\\multicolumn\{2\}\{c\}\{\$(.*?)\$\}", block)
    qf = re.search(r"quadratic form.*?\\begin\{array\}\{[^}]*\}(.*?)\\end\{array\}", block, re.S)
    rl = re.search(r"root list.*?\\begin\{array\}\{[^}]*\}(.*?)\\end\{array\}", block, re.S)
    if not (det and shape and qf and rl):
        raise ValueError(f"{src} entry {idx}: missing field")
    gram_rows = rows_of(qf.group(1), r"^-?\d", 6)
    root_rows = rows_of(rl.group(1), r"^r_\{?\d+\}?", 8)
    if len(gram_rows) != 3:
        raise ValueError(f"{src} entry {idx}: {len(gram_rows)} Gram rows")
    sh = shape.group(1).replace("\\left(", "(").replace("\\right)", ")").replace(" ", "")
    return {
        "file": src, "index_in_file": idx,
        "det": parse_elem(det.group(1)), "det_norm": int(det.group(2).replace(" ", "")),
        "shape": sh,
        "gram": [[r[2 * j:2 * j + 2] for j in range(3)] for r in gram_rows],
        "roots": [[r[2 * j:2 * j + 2] for j in range(3)] for r in root_rows],
        "norms": [r[6:8] for r in root_rows],
    }


def main(src_dir, out):
    entries = []
    for f in FILES:
        text = (Path(src_dir) / f"{f}.tex").read_text()
        text = "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("%"))
        blocks = text.split("\\hline\\hline")[1:]
        for i, b in enumerate(blocks):
            e = parse_block(b, f, i)
            e["id"] = len(entries)
            entries.append(e)
    Path(out).write_text(json.dumps({"source": "Mark, arXiv:1512.01133v1, Appendix B (TeX source)",
                                     "n_entries": len(entries), "entries": entries}, indent=0))
    print(len(entries), "entries")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
