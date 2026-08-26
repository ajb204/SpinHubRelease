from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]




def test_peak_algorithms_have_canonical_analysis_ownership():
    fit = (ROOT / "gui" / "workspaces" / "peak_fit.py").read_text()
    assert "decon.analysis.peak_picker" in fit
    assert "decon.analysis.peak_shape_estimator" in fit
    assert (ROOT / "analysis" / "peak_picker.py").exists()
    assert (ROOT / "analysis" / "peak_shape_estimator.py").exists()
