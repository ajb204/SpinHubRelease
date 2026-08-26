"""Regression checks for mechanical errors introduced during dimension refactoring."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MIGRATED_VIEWERS = (
    "gui/workspaces/oned.py",
    "gui/workspaces/projection.py",
    "gui/workspaces/pseudo2d.py",
    "gui/workspaces/pseudo3d.py",
    "gui/workspaces/slice1d.py",
    "gui/workspaces/slice2d.py",
)


def test_viewer_dimension_contract_imports_only_contract_symbols():
    forbidden = {"string", "copy", "math", "numpy", "os", "re", "sys"}
    for relative in MIGRATED_VIEWERS:
        text = (ROOT / relative).read_text()
        for line in text.splitlines():
            if line.startswith("from spinDecon.domain.dimensions.viewer_contract import"):
                imported = {name.strip() for name in line.split("import", 1)[1].split(",")}
                assert not (imported & forbidden), (relative, imported & forbidden)


def test_process_frame_has_no_double_self_dimension_access():
    text = (ROOT / "gui/dialogs/processing/process.py").read_text()
    assert "self.self.has_pseudo_axis" not in text
    assert "inst.self.has_pseudo_axis" not in text


def test_peakframe_pseudo3d_selection_uses_canonical_topology():
    text = (ROOT / "gui/workspaces/peak_review.py").read_text()
    start = text.index("def _projection_payload")
    end = text.index("def _display_payload", start)
    projection = text[start:end]
    assert "getattr(self.tabOne, 'pseudo'" not in projection
    assert "self.peak_service.projection_payload" in projection
    assert "self.tabOne" not in projection


def test_projection_distinguishes_true_2d_from_two_spectral_plus_pseudo():
    text = (ROOT / "gui/workspaces/projection.py").read_text()
    start = text.index("def create_main_panel")
    end = text.index("def create_status_bar", start)
    panel = text[start:end]
    assert "if self._is_true_2d_spectrum():" in panel
    assert "if self._is_3p_projection_case():" in panel
    assert "if(self.spectral_dim_count==3):\n            if self._is_3p_projection_case()" not in panel
    assert "Projection requires cached raw 3p projection view" in panel


def test_projection_3p_decon_path_is_reachable_with_two_spectral_dimensions():
    text = (ROOT / "gui/workspaces/projection.py").read_text()
    start = text.index("def load_decon_data")
    end = text.index("def background_save", start)
    block = text[start:end]
    assert "self._is_3p_projection_case() or self.spectral_dim_count == 3" in block


def test_peakframe_decon_is_explicitly_a_2d_analysis_target():
    text = (ROOT / "gui/workspaces/peak_review.py").read_text()
    start = text.index("def OnButtonDecon")
    end = text.index("def OnButtonRecon", start)
    block = text[start:end]
    assert "target = self._analysis_spectrum_path()" in block
    assert "caller='peakframe'" in block
    assert "input_override=target" in block
    assert "dimension_override=2" in block


def test_peakframe_recon_uses_same_2d_target_and_its_2d_peak_list():
    text = (ROOT / "gui/workspaces/peak_review.py").read_text()
    start = text.index("def OnButtonRecon")
    end = text.index("def OnButtonSaveDecon", start)
    block = text[start:end]
    assert "peak_path = self._peak_list_path()" in block
    assert "self.SavePeakList(peak_path)" in block
    assert "recon=True" in block
    assert "input_override=target" in block
    assert "dimension_override=2" in block
    assert "peak_list_override=peak_path" in block


def test_main_decon_keeps_2d_bore_hybrid_mode():
    text = (ROOT / "gui/workspaces/nmr.py").read_text()
    start = text.index("def OnButtonDecon")
    end = text.index("#Analyse results from decon", start)
    block = text[start:end]
    assert "dec3d=self.cb_decon3d.IsChecked()*1. #bore mode?" in block
    assert "self.OnButtonReadPeak(None)" in block
    assert "self._resolve_spec_file(self.referencePeakBox.GetValue())" in block
    assert "decset['dec3d']=1" in block


def test_recon_ignores_bore_checkbox_and_accepts_explicit_peakframe_list():
    text = (ROOT / "gui/workspaces/nmr.py").read_text()
    start = text.index("def OnButtonDecon")
    end = text.index("#Analyse results from decon", start)
    block = text[start:end]
    assert "if recon:" in block
    assert "dec3d=0" in block
    assert "peak_list_override" in block
    assert "Protocol3P applies only to a main-spectrum job" in block
