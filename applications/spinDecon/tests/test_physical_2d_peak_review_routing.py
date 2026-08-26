from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_review_peaks_routes_physical_2d_to_projection_and_peakframe():
    text = (ROOT / 'app' / 'workflow_controller.py').read_text(encoding='utf-8')
    block = text[text.index("if action_key == 'review_peaks':"):text.index("if action_key == 'fit_spectrum':")]
    assert "mode.spectral_dimensions == 2 and not mode.has_pseudo_axis" in block
    assert "self.AddTabTwo(True, tab)" in block
    assert "tab.OnButtonPeaky(None)" in block
    assert "return self.select_page('Projections')" in block


def test_full_2d_show_synchronises_projection_and_peakframe():
    text = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text(encoding='utf-8')
    start = text.index('    def select_full_peak(')
    block = text[start:text.index('    def get_full3d_peak_selection_context', start)]
    assert "topology.spectral_dim_count == 2 and not topology.has_pseudo_axis" in block
    assert "focus_2d_peak" in block
    assert "self._live_peak_frame()" in block
    assert "focus_peak(selected_peak)" in block
    assert "self._physical_2d_peak_for_full_record(record)" in block


def test_projection_exposes_physical_2d_focus_api():
    text = (ROOT / 'gui' / 'workspaces' / 'projection.py').read_text(encoding='utf-8')
    assert '    def focus_2d_peak(self, peak, width_fraction=0.10):' in text
    assert 'self.axes.set_xlim' in text
    assert 'self.axes.set_ylim' in text


def test_peakframe_external_focus_reuses_persistent_selection_ornament():
    text = (ROOT / 'gui' / 'workspaces' / 'peak_review.py').read_text(encoding='utf-8')
    start = text.index('    def focus_peak(self, peak, width_fraction=0.10):')
    block = text[start:text.index('    def read_peaklist_file', start)]
    assert "self.select = [selected_index]" in block
    assert "self.selection_type = 'single'" in block
    assert 'self._update_selection_artists()' in block


def test_physical_2d_projection_focus_draws_selected_peak_ornament():
    text = (ROOT / 'gui' / 'workspaces' / 'projection.py').read_text(encoding='utf-8')
    start = text.index('    def focus_2d_peak(self, peak, width_fraction=0.10):')
    block = text[start:text.index('    def select_full_peak_from_list', start)]
    assert 'self._selected_2d_peak_xy = (x, y)' in block
    assert "marker='x'" in block
    assert 's=140' in block
    assert '0.3, 0.5, 0.0, 1.0' in block


def test_reference_selection_routes_to_physical_2d_projection_highlight():
    text = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text(encoding='utf-8')
    start = text.index('    def select_reference_peak(self, peak_name):')
    block = text[start:text.index('    def _select_pseudo2d_full_peak_in_projection', start)]
    assert 'topology.spectral_dim_count == 2 and not topology.has_pseudo_axis' in block
    assert 'focus_projection(peak)' in block
