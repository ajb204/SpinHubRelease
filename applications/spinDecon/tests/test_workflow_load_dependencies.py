from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_review_materialises_all_inputs_in_dependency_order():
    text = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def ensure_workflow_review_inputs_loaded(self):')
    end = text.index('    def ensure_reference_peak_list_loaded(self):', start)
    block = text[start:end]
    assert block.index('ensure_workflow_spectrum_loaded()') < block.index('ensure_reference_peak_list_loaded()')
    assert block.index('ensure_reference_peak_list_loaded()') < block.index('ensure_full_peak_list_loaded()')
    assert block.index('ensure_full_peak_list_loaded()') < block.index('ensure_deconvolution_loaded()')


def test_physical_3d_pick_materialises_reference_list():
    text = (ROOT / 'app' / 'notebook.py').read_text()
    start = text.index("        if action_key == 'peak_pick':")
    end = text.index("        if action_key == 'review_peaks':", start)
    block = text[start:end]
    assert 'mode.spectral_dimensions >= 3' in block
    assert "ensure_reference_peak_list_loaded" in block
    assert block.index('ensure_reference_peak_list_loaded') < block.index('tab.OnButtonDecon(None)')


def test_peak_shape_and_reference_actions_load_spectrum_first():
    text = (ROOT / 'app' / 'notebook.py').read_text()
    peak_shape = text[text.index("        if action_key == 'peak_shape':"):text.index("        if action_key == 'reference_peaks':")]
    reference = text[text.index("        if action_key == 'reference_peaks':"):text.index("        if action_key == 'peak_pick':")]
    assert 'if not ensure_loaded(): return False' in peak_shape
    assert 'if not ensure_loaded(): return False' in reference


def test_manual_decon_load_obeys_spectrum_reference_dependency():
    text = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonAnalyse(self,event):')
    end = text.index('    def ', start + 8)
    block = text[start:end]
    assert block.index('ensure_workflow_spectrum_loaded()') < block.index('ensure_reference_peak_list_loaded()')
    assert block.index('ensure_reference_peak_list_loaded()') < block.index('_load_decon_outputs')
