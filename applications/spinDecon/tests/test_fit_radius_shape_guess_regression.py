from pathlib import Path

TEXT = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_fit.py').read_text()


def test_radius_guess_uses_current_fit_shape_not_raw_level_scan():
    start = TEXT.index('    def on_guess(self,event=None):')
    end = TEXT.index('    def on_save(self,event=None):', start)
    body = TEXT[start:end]
    assert '_radii_from_current_peak_shape(level=.10)' in body
    assert 'estimate_level_radius(' not in body
    assert 'self.owner.data' not in body


def test_radius_guess_status_escapes_literal_percent():
    start = TEXT.index('    def on_guess(self,event=None):')
    end = TEXT.index('    def on_save(self,event=None):', start)
    body = TEXT[start:end]
    assert "'10%% radius from fitted peak shape: %s'%detail" in body


def test_shape_radius_uses_all_current_pseudovoigt_controls():
    assert 'self.psf_sliders[dim].val' in TEXT
    assert 'self.lorentz_sliders[dim].val' in TEXT
    assert 'self.voigt_sliders[dim].val' in TEXT
    assert 'gaussian=numpy.exp(-4.0*numpy.log(2.0)*(radius/g)**2)' in TEXT
    assert 'lorentzian=1.0/(1.0+4.0*(radius/l)**2)' in TEXT
