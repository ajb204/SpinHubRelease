from pathlib import Path
ROOT = Path(__file__).parents[1]
DECON = (ROOT / "gui" / "workspaces" / "nmr.py").read_text()
PEAK = (ROOT / "gui" / "workspaces" / "peak_review.py").read_text()

def test_peakframe_is_singleton_per_decon_frame():
    assert "def _live_peak_frame(self):" in DECON
    assert "existing = self._live_peak_frame()" in DECON
    assert "existing.Raise()" in DECON

def test_threshold_refreshes_existing_peakframe_without_opening_one():
    start = DECON.index("def OnButtonNoise")
    end = DECON.index("def prepare_workflow", start)
    block = DECON[start:end]
    assert "peak_frame = self._live_peak_frame()" in block
    assert "sync_main_threshold" in block
    assert "OnButtonPeaky" not in block

def test_peakframe_can_sync_threshold_in_place():
    assert "def sync_main_threshold(self, redraw=True):" in PEAK
    assert "self.tabOne.threshBox.GetValue()" in PEAK
    assert "self.ax_reset = 0" in PEAK
