from pathlib import Path


def test_peakframe_threshold_sync_updates_minimum_contour_control():
    src = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_review.py').read_text()
    start = src.index('    def sync_main_threshold')
    end = src.index('    def _analysis_spectrum_path', start)
    block = src[start:end]
    assert 'self.thresh = self.tabOne.dmax * float(self.tabOne.threshBox.GetValue())' in block
    assert "self.textbox0.SetValue(minimum)" in block
    assert 'self.ax_reset = 0' in block
    assert 'self.draw_figure()' in block
