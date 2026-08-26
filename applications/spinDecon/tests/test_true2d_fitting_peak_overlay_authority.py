from pathlib import Path


def test_true2d_projected_overlay_is_published_from_full_peak_list():
    source = Path('gui/workspaces/nmr.py').read_text()
    block = source[source.index('    def _rebuild_projected_peak_lists'):source.index('    def get_projected_peak_overlay')]
    two_d = block[block.index("if int(getattr(self, 'dim', 0) or 0) == 2:"):block.index("if int(getattr(self, 'dim', 0) or 0) != 3:")]
    assert 'self.get_full_peak_payload()' in two_d
    assert "store.projected_peak_lists[('full', display_x, display_y)]" in two_d
    assert "'source': 'full'" in two_d
    assert 'self.get_reference_peaks()' not in two_d
    assert "('reference', display_x, display_y)" not in two_d


def test_true2d_fitting_peaks_toolbar_requests_full_overlay():
    source = Path('gui/workspaces/pseudo3d.py').read_text()
    block = source[source.index('    def _reference_peak_overlay'):source.index('    def _reference_peak_by_name')]
    assert "list_key = 'full' if self._is_physical_2d_adapter() else 'reference'" in block
