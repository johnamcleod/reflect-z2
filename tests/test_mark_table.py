"""Mark's appendix (arXiv:1512.01133), parsed by tools/parse_mark.py: every entry passes
exact consistency checks (reflect.mark.check_entry)."""
import pytest

from reflect.mark import load, check_entry

ENTRIES = load()


def test_entry_count():
    # the paper's abstract says 432 lattices; the TeX appendix has 522 entries, several of
    # them isometric to each other (see data/isometries and OPEN_QUESTIONS OQ4)
    assert len(ENTRIES) == 522


@pytest.mark.parametrize("e", ENTRIES, ids=lambda e: f"{e['id']}_{e['shape']}")
def test_entry_consistent(e):
    res = check_entry(e)
    assert all(res.values()), {k: v for k, v in res.items() if not v}
