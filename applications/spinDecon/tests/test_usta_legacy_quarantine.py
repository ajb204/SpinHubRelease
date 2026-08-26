from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]




def test_superseded_std_source_is_quarantined():
    legacy = ROOT / 'legacy' / 'usta' / 'std_frame_legacy.py'
    assert legacy.is_file()
    assert 'class STDFrame' in legacy.read_text()
