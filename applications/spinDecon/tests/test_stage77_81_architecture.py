from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_gui_does_not_import_deconframe_module():
    for base in (ROOT / 'gui', ROOT / 'integrations'):
        for path in base.rglob('*.py'):
            text = path.read_text()
            assert 'Frames.deconFrame' not in text, path
            assert 'from .deconFrame' not in text, path
            assert 'from ..deconFrame' not in text, path






def test_peak_io_helper_is_gui_independent():
    text = (ROOT / 'processing' / 'peak_io.py').read_text()
    assert 'import wx' not in text
    assert 'from spinDecon.domain.peaks import peakEntry' in text
