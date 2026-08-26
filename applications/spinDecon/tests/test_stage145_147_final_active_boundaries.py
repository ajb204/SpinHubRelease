from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_slice2d_scientific_state_is_service_owned():
    text = (ROOT / 'gui/workspaces/slice2d.py').read_text()
    constructor = text[text.index('def __init__(self,parent,tabOne):'):text.index('def _make_modeless_window')]
    assert 'self.tabOne=' not in constructor
    assert 'parent.tabOne' not in constructor
    assert 'SliceService(tabOne)' in constructor


def test_peak_fit_no_longer_stores_legacy_workspace():
    text = (ROOT / 'gui/workspaces/peak_fit.py').read_text()
    constructor = text[text.index('class peakFitFrame'):text.index('def create_main_panel')]
    assert 'self.tabOne = parent' not in constructor
    assert 'self.peak_fit_service.visible_axes' in text


def test_peak_review_uses_peak_service_for_projection_and_persistence():
    text = (ROOT / 'gui/workspaces/peak_review.py').read_text()
    assert 'self.tabOne=parent' not in text
    projection = text[text.index('def _projection_payload'):text.index('def _display_payload')]
    assert 'self.peak_service.projection_payload' in projection
    assert 'self.tabOne' not in projection
    commit = text[text.index('def _commit_projection_peaks'):text.index('def focus_peak')]
    assert 'self.peak_service.commit_projection_peaks' in commit
    assert 'self.tabOne' not in commit
