from types import SimpleNamespace
from spinDecon.app.context import ApplicationContext
from spinDecon.analysis.services import attach_analysis_services
from spinDecon.workflow.status import _legacy_workspace


def test_context_exposes_canonical_nmr_workspace_and_keeps_legacy_alias():
    workspace = SimpleNamespace()
    context = ApplicationContext()
    attach_analysis_services(context, workspace)
    assert context.nmr_workspace is workspace
    assert context.legacy_nmr_workspace is workspace
    assert _legacy_workspace(context) is workspace
