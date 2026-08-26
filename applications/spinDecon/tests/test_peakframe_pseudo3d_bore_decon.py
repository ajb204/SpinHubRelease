from pathlib import Path

TEXT = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_review.py').read_text()

def test_pseudo3d_gets_one_dimensional_bore():
    assert 'def _has_1d_bore(self):' in TEXT
    assert 'return self.spectral_dim_count == 3 or self._is_pseudo3d_dataset()' in TEXT
    assert "view = self.peak_service.pseudo3d_view('raw')" in TEXT
    assert "trace = cube[:, y1, x1]" in TEXT

def test_pseudo3d_peakframe_decon_reads_projection_dot_decon():
    service = (Path(__file__).parents[1] / 'analysis' / 'peak_service.py').read_text()
    assert "decon_path = (analysis_path or '') + '.decon'" in service
    assert "key = ('peakframe_decon', labels[0], labels[1], 'n')" in service
    assert "store.save_view(key, **view)" in service
