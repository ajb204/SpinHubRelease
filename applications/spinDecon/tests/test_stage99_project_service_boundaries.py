from pathlib import Path


def test_project_service_uses_canonical_application_and_workflow_boundaries():
    source = (Path(__file__).parents[1] / "project" / "service.py").read_text()
    assert "from spinDecon.app.launcher import open_project_gui" in source
    assert "from spinDecon.workflow.registry import WORKFLOW_BY_KEY, WORKFLOW_REGISTRY" in source
    assert "decon.decon_tab" not in source
    assert "decon.workflow_registry" not in source


def test_workflow_status_prefers_explicit_legacy_workspace_bridge():
    from types import SimpleNamespace
    from spinDecon.workflow.status import _legacy_workspace
    explicit = object()
    old = object()
    assert _legacy_workspace(SimpleNamespace(legacy_nmr_workspace=explicit, tabOne=old)) is explicit
    assert _legacy_workspace(SimpleNamespace(tabOne=old)) is old
    assert _legacy_workspace(None) is None
