from pathlib import Path


def test_projection_does_not_own_historical_line_fitting_mode():
    source = Path('gui/workspaces/projection.py').read_text()
    forbidden = (
        'Unidec_line_fitting',
        'self.line_fitting',
        'def peak_fit(',
        'def fuda_thread(',
        'def overlap_thread(',
        'def print_lorentzian(',
    )
    for token in forbidden:
        assert token not in source


def test_peakfit_remains_canonical_nmr_launched_fitting_workspace():
    nmr = Path('gui/workspaces/nmr.py').read_text()
    peak_fit = Path('gui/workspaces/peak_fit.py').read_text()
    assert 'def OnButtonPeakFit(' in nmr
    assert 'decon.gui.workspaces import peak_fit as peakFitFrame' in nmr
    assert 'class peakFitFrame' in peak_fit
    assert 'PeakFitService' in peak_fit


def test_projection_line_fitting_history_is_quarantined():
    root = Path('legacy/projection_line_fitting')
    assert (root / 'README.md').exists()
    assert (root / 'projection_pre_cleanup_snapshot.py').exists()
