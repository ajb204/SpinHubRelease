from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_project_summary_uses_canonical_peak_workspace_imports():
    source = (ROOT / 'gui' / 'reporting' / 'project_summary.py').read_text()
    assert 'from spinDecon.gui.workspaces.peaks import peakFrame' in source
    assert 'from spinDecon.gui.workspaces.peaks import peakFitFrame' in source
    assert 'from .Frames.peakFrame' not in source
    assert 'from .Frames.peakFitFrame' not in source


def test_peak_workspace_exports_full_peak_list_boundary():
    source = (ROOT / 'gui' / 'workspaces' / 'peaks.py').read_text()
    assert 'PeakListFrame' in source
    assert 'conn_data' not in source
