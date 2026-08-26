from pathlib import Path


def test_fid_phase_rule_includes_2d_second_trace():
    source = (Path(__file__).parents[1] / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
    assert '(ndim == 2 and sel == 2)' in source
    assert '(ndim == 3 and sel in (3, 4))' in source
