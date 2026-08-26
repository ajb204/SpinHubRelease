from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_active_notebook_does_not_connect_legacy_usta_tab():
    source = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'STDFrame' not in source
    assert 'AddTabSTD' not in source
    assert '("uSTA",' not in source
