from pathlib import Path


def test_physical_3d_workflow_pick_forces_use_2d_peaklist():
    text = (Path(__file__).parents[1] / 'app' / 'notebook.py').read_text()
    start = text.index("        if action_key == 'peak_pick':")
    end = text.index("        if action_key == 'review_peaks':", start)
    block = text[start:end]
    assert "not mode.has_pseudo_axis and mode.spectral_dimensions >= 3" in block
    assert "use_2d = getattr(tab, 'cb_decon3d', None)" in block
    assert "use_2d.SetValue(True)" in block
    assert block.index("use_2d.SetValue(True)") < block.index("tab.OnButtonDecon(None)")
