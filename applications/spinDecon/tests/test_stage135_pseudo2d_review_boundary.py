from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def test_pseudo2d_review_persistence_uses_pseudo_service():
    source = (ROOT / 'gui/workspaces/pseudo2d.py').read_text()
    method = source[source.index('    def on_mark_intensities_correct'):source.index('    def _sync_review_button')]
    assert "self.pseudo_service.mark_series_reviewed(source='pseudo2d_fitting')" in method
    assert 'self.tabOne' not in method
