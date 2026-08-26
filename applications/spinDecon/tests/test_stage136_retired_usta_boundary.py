from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def test_active_notebook_contains_no_usta_tab_route():
    source = (ROOT / 'app/notebook.py').read_text()
    assert 'AdduSTA' not in source
    assert 'uSTA_sims' not in source
