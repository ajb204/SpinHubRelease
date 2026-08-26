from pathlib import Path
SRC = Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py'
TEXT = SRC.read_text()
def test_pseudo2d_fit_uses_original_data_and_protocol_flag():
    assert "decset['pseudo2DFit'] = '1'" in TEXT
    assert "decset['infile'] = self._resolve_input_path(self.infileBox.GetValue())" in TEXT
    assert "self.cb_decon3d.IsChecked() and self.cb_decback.IsChecked()" in TEXT
def test_pseudo2d_shows_only_f1_fixed_radius():
    assert 'show_f1_fit = (topology.spectral_dim_count == 2) or (topology.spectral_dim_count == 1 and topology.has_pseudo_axis)' in TEXT
    assert 'show_f2_fit = (topology.spectral_dim_count == 2)' in TEXT
