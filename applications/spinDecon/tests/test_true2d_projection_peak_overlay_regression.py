from pathlib import Path


def _source():
    return Path('gui/workspaces/projection.py').read_text()


def test_full_peak_records_are_supported_by_projection_overlay_converter():
    source = _source()
    block = source[source.index('    def _peak_points_for_overlay'):source.index('    def _projection_view')]
    assert 'isinstance(entry, dict)' in block
    assert 'entry.get("coordinates")' in block
    assert 'entry.get("name"' in block


def test_true2d_draw_renders_peak_overlay_from_persistent_toolbar_state():
    source = _source()
    block = source[source.index('    def draw_2d'):source.index('    def load_decon_data')]
    assert 'if self.cb_grid.IsChecked()' in block
    assert 'self._cached_peak_overlay()' in block
    assert 'self._peak_points_for_overlay(full_peaks, swap=True)' in block
    assert 'scatter_xy_points' in block


def test_true2d_peaks_toolbar_redraws_without_resetting_zoom():
    source = _source()
    block = source[source.index('    def on_cb_grid(self, event):'):source.index('    def on_cb_grid_auto', source.index('    def on_cb_grid(self, event):'))]
    assert 'self.draw_2d(keepaxes=True)' in block
    assert 'self.canvas.draw_idle()' in block
