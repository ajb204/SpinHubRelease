from types import SimpleNamespace
from spinDecon.analysis.services import attach_analysis_services


def test_analysis_service_composition_attaches_complete_service_set():
    context = SimpleNamespace()
    workspace = object()
    result = attach_analysis_services(context, workspace)
    assert result is context
    assert context.legacy_nmr_workspace is workspace
    for name in ("full3d", "one_d", "projection", "peaks", "full_peaks", "peak_fit", "slices", "pseudo", "diffusion"):
        assert getattr(context, name) is not None


def test_notebook_delegates_service_wiring_to_composition_root():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "app" / "notebook.py").read_text()
    assert "attach_analysis_services(self.app_context, self.nmr_workspace)" in source
    assert "Full3DService(self.tabOne)" not in source
    assert "PeakService(self.tabOne)" not in source
