"""Regression coverage for removal of deprecated NMR topology widgets."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECON_FRAME = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
DECON_TAB = (ROOT / 'app' / 'notebook.py').read_text()
WORKFLOW = (ROOT / 'gui' / 'workspaces' / 'workflow.py').read_text()
PROCESS = (ROOT / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
PEAK = (ROOT / 'gui' / 'workspaces' / 'peak_review.py').read_text()
VIEWER = (ROOT / 'domain' / 'dimensions' / 'viewer_contract.py').read_text()


def _block(text, name, next_name=None):
    start = text.index('    def ' + name)
    if next_name:
        end = text.index('    def ' + next_name, start)
    else:
        end = text.find('\n    def ', start + 8)
        if end < 0:
            end = len(text)
    return text[start:end]


def test_spectrum_box_no_longer_constructs_or_lays_out_topology_widgets():
    block = _block(DECON_FRAME, 'spectrum_box(self):', 'pre_read_disabling')
    assert 'wx.RadioBox' not in block
    assert 'pseudoAxis' not in block
    assert 'dimBox' not in block
    assert 'pseudoBox' not in block
    assert 'dimensionRow' not in block
    assert 'buttonProcess' in block and 'buttonRead' in block and 'buttonReOrganise' in block


def test_workflow_routes_topology_through_state_api_not_nmr_widgets():
    block = _block(DECON_TAB, 'apply_workflow_dataset_type(self, spectral_dimensions, pseudo_axis):', 'mark_workflow_series_inspected')
    assert "getattr(tab, 'apply_dataset_type', None)" in block
    assert 'dimBox' not in block
    assert 'pseudoBox' not in block


def test_dataset_type_api_preserves_pseudo4d_guard_and_refreshes_dimensions():
    block = _block(DECON_FRAME, 'apply_dataset_type(self, spectral_dimensions, pseudo_axis, *, show_error=True):', '_active_topology')
    assert 'if dim == 4 and pseudo:' in block
    assert "errorMessage('pseudo4d not yet supported')" in block
    assert 'self.state.sync_from_values(spectral_dimensions=dim, pseudo_axis=pseudo)' in block
    assert 'self.SetDim()' in block
    assert 'self._update_full_peak_controls()' in block


def test_remaining_consumers_do_not_read_removed_widgets():
    for text in (WORKFLOW, PROCESS, PEAK, VIEWER):
        assert '.dimBox' not in text
        assert '.pseudoBox' not in text


def test_deprecated_topology_widget_api_is_absent_from_production_python():
    """The removed controls/handlers must not survive as a hidden compatibility API."""
    forbidden = ('dimBox', 'pseudoBox', 'dimEntry', 'pseudoBoxCheck')
    offenders = {}
    for path in ROOT.rglob('*.py'):
        if 'tests' in path.parts or 'legacy' in path.parts:
            continue
        source = path.read_text(errors='ignore')
        hits = [name for name in forbidden if name in source]
        if hits:
            offenders[str(path.relative_to(ROOT))] = hits
    assert offenders == {}
