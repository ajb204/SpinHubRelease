from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_retired_magma_and_slice4d_are_quarantined_not_exposed_by_notebook():
    text = (ROOT / 'app/notebook.py').read_text()
    assert 'def AddMagmaTab' not in text
    assert 'def AddTabFour4D' not in text
    assert (ROOT / 'legacy/magma').is_dir()
    assert (ROOT / 'legacy/slice4d').is_dir()
