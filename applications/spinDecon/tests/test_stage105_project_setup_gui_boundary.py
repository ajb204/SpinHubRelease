from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_project_setup_dialog_has_gui_owner():
    canonical = (ROOT / 'gui' / 'dialogs' / 'project_setup.py').read_text()
    assert 'class ProjectSetupDialog(wx.Dialog)' in canonical
    assert 'def run_project_setup' in canonical
    assert 'from spinDecon.project.service import ProjectService' in canonical


def test_application_uses_gui_project_setup_boundary():
    app = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'from spinDecon.gui.dialogs.project_setup import run_project_setup' in app
    assert 'from ..project.setup import run_project_setup' not in app
