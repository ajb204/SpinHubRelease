from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECTION = ROOT / 'gui' / 'workspaces' / 'projection.py'
NMR = ROOT / 'gui' / 'workspaces' / 'nmr.py'


def test_projection_2d_decon_toolbar_redraws_main_contours():
    text = PROJECTION.read_text()
    start = text.index('    def _toolbar_decon(self, active):')
    end = text.index('    def _toolbar_peaks', start)
    body = text[start:end]
    assert 'self.draw_2d(keepaxes=True)' in body
    assert 'self._configure_2d_trace_visibility()' in body


def test_pseudo3d_2d_decon_is_optional_and_toolbar_tracks_data():
    text = (ROOT / 'gui' / 'workspaces' / 'pseudo3d.py').read_text()
    assert 'def _decon_pseudo3d_view(self):' in text
    assert 'decon_available = decon_view is not None' in text
    assert 'self.toolbar.enable_decon(decon_available)' in text
    assert "colors='red'" in text
    assert "colors='green'" in text


def test_pseudo3d_decon_toggle_can_rebuild_after_late_decon_run():
    text = (ROOT / 'gui' / 'workspaces' / 'pseudo3d.py').read_text()
    start = text.index('    def _toolbar_decon(self, active):')
    end = text.index('    def _toolbar_peaks', start)
    body = text[start:end]
    assert 'self._decon_pseudo3d_view()' in body
    assert "not getattr(self, 'calc_artists', [])" in body
    assert 'self.draw_figureGO()' in body


def test_true_2d_pseudo3d_decon_never_falls_back_to_raw():
    text = NMR.read_text()
    start = text.index('    def get_pseudo3d_view(')
    end = text.index('    def get_spectrum_view(', start)
    body = text[start:end]
    assert "if key == 'raw':" in body
    assert "data = spectrum.get('data')" in body


def test_publishing_decon_invalidates_cached_pseudo3d_decon_view():
    text = (ROOT / 'project' / 'data_store.py').read_text()
    assert 'self.views.pop(("pseudo3d", "decon"), None)' in text


def test_pseudo3d_contour_toggle_supports_modern_matplotlib_artist_api():
    text = (ROOT / 'gui' / 'workspaces' / 'pseudo3d.py').read_text()
    assert "collections = getattr(cs, 'collections', None)" in text
    assert "else [cs]" in text
