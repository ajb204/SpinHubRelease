from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_workflow_package_no_longer_owns_wx_action_routing():
    canonical = (ROOT / 'app' / 'workflow_controller.py').read_text()
    compat = (ROOT / 'workflow' / 'actions.py').read_text()
    notebook = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'class WorkflowController' in canonical
    assert 'import wx' in canonical
    assert 'from spinDecon.app import workflow_controller as _canonical' in compat
    assert 'from .workflow_controller import WorkflowController' in notebook


def test_workflow_domain_modules_are_wx_free():
    for name in ('model.py', 'status.py', 'registry.py'):
        source = (ROOT / 'workflow' / name).read_text()
        assert 'import wx' not in source
        assert 'from wx' not in source
