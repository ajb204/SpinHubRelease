from types import SimpleNamespace
from spinDecon.analysis.pseudo_service import PseudoAxisService


class _Store:
    def __init__(self):
        self.analysis = {}
        self.metadata = {}
    def mark_pseudo_series_reviewed(self, **details):
        self.analysis['pseudo_series_reviewed'] = details or True


class _Box:
    def GetValue(self):
        return 'peaks.tab'


def test_pseudo_service_owns_fitting_review_and_full_peak_resolution(tmp_path):
    calls = []
    legacy = SimpleNamespace(
        store=_Store(), data_store=None, fullPeakBox=_Box(),
        get_fuda_dir=lambda: str(tmp_path / 'fit'),
        _resolve_spec_file=lambda value: str(tmp_path / value),
        OnButtonSave=lambda *args: calls.append('save'),
        _notify_analysis_changed=lambda: calls.append('notify'),
    )
    service = PseudoAxisService(legacy)
    assert service.fit_dir() == str(tmp_path / 'fit')
    assert service.full_peak_file() == str(tmp_path / 'peaks.tab')
    assert service.mark_series_reviewed()
    assert service.series_reviewed()
    assert calls == ['save', 'notify']


def test_canonical_pseudo_and_oned_workspace_imports_are_used_by_notebook():
    from pathlib import Path
    root = Path(__file__).parents[1]
    source = (root / 'app' / 'notebook.py').read_text()
    assert 'from ..gui.workspaces.oned import OneDFrame' in source
    assert 'from ..gui.workspaces.pseudo2d import Pseudo2D' in source
