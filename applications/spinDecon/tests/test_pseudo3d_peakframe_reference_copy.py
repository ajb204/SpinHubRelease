from pathlib import Path


def test_pseudo3d_peakframe_reference_copy_is_wired_into_completion():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    assert "def _promote_pseudo3d_peakframe_reference_list" in source
    assert "destination = os.path.join(spec_dir, 'test.ft2.2D.list')" in source
    assert "self.referencePeakBox.SetValue(relative)" in source
    assert "self.state.reference_peak_file = relative" in source
    assert "caller == 'peakframe' and topology.spectral_dim_count == 2 and topology.has_pseudo_axis" in source
    assert "peak_path = self._promote_pseudo3d_peakframe_reference_list(peak_path)" in source


def test_peakframe_recon_saves_current_list_before_async_run():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_review.py').read_text()
    recon = source[source.index('    def OnButtonRecon'):source.index('    def OnButtonSaveDecon')]
    assert 'self.SavePeakList(peak_path)' in recon
    assert "caller='peakframe', recon=True" in recon
