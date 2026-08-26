from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_workflow_overview_is_a_gui_workspace_without_active_compat_alias():
    canonical = (ROOT / 'gui' / 'workspaces' / 'workflow.py').read_text()
    assert 'class WorkflowOverviewPanel' in canonical
    assert not (ROOT / 'workflow' / 'overview.py').exists()
    assert (ROOT / 'legacy' / 'compatibility' / 'workflow' / 'overview.py').exists()
