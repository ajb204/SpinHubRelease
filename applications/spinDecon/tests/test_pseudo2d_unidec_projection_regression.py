from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_pseudo2d_fit_peaks_uses_spectral_projection():
    text = (ROOT / 'gui' / 'workspaces' / 'peak_fit.py').read_text()
    assert "elif self.dim == 1:" in text
    assert "get_pseudo2d_projection_data(ensure_file=True)" in text
    assert "array = self.data" in text


def test_pseudo2d_decon_targets_projection_file_and_publishes_results():
    text = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    assert "def get_pseudo2d_projection_data" in text
    assert "nmrglue_project2D_1D" in text
    assert "decset['infile'] = pseudo_spec" in text
    assert "def _load_pseudo2d_projection_decon_outputs" in text
    assert "save_spectrum('pseudo2d_projection_decon'" in text
    assert "self._load_pseudo2d_projection_peaks(self.corrFile)" in text
    assert "if spectral_dim_count == 1 and topology.has_pseudo_axis:" in text


def test_projection_window_can_show_pseudo2d_peaks_and_deconvolution():
    text = (ROOT / 'gui' / 'workspaces' / 'projection.py').read_text()
    assert "def _pseudo2d_projection(self, decon=False)" in text
    assert "pseudo2d_projection_decon" in text
    assert "def _pseudo2d_peak_overlay" in text
    assert "if self.cb_grid.IsChecked():" in text
    assert "if self.cb_calc.IsChecked():" in text
