from pathlib import Path

from spinDecon.processing.peak_io import read_peak_list


def test_canonical_peak_io_is_available():
    assert callable(read_peak_list)


def test_active_gui_uses_processing_peak_io():
    root = Path(__file__).resolve().parents[1]
    source = (root / "gui" / "workspaces" / "peak_review.py").read_text()
    assert "decon.processing.peak_io" in source
    assert "decon.misc.peak_io" not in source
