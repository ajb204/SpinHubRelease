from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_diffusion_service_is_gui_independent_and_wired_to_context():
    service = (ROOT / 'analysis' / 'diffusion_service.py').read_text()
    context = (ROOT / 'app' / 'context.py').read_text()
    composition = (ROOT / 'analysis' / 'services.py').read_text()
    notebook = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'import wx' not in service
    assert 'diffusion: Any = None' in context
    assert 'context.diffusion = DiffusionService(legacy_workspace)' in composition
    assert 'attach_analysis_services(self.app_context, self.nmr_workspace)' in notebook
