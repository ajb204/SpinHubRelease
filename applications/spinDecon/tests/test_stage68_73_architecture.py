from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]




def test_canonical_gui_code_no_longer_imports_frames_toolbar_or_widgets():
    text = '\n'.join(p.read_text(errors='ignore') for p in (ROOT / 'gui').rglob('*.py'))
    assert 'decon.Frames.matplotlib_toolbar' not in text
    assert 'decon.Frames.widgets' not in text
