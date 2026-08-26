from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]




def test_active_consumers_use_canonical_pdf_viewer_import():
    for rel in ('gui/workspaces/nmr.py',):
        text = (ROOT / rel).read_text()
        assert 'from spinDecon.gui.dialogs.pdf_viewer import PDFViewer' in text
        assert 'from spinDecon.pdfViewer import PDFViewer' not in text


def test_canonical_workflow_package_has_self_consistent_imports():
    model = (ROOT / 'workflow/model.py').read_text()
    status = (ROOT / 'workflow/status.py').read_text()
    overview = (ROOT / 'gui/workspaces/workflow.py').read_text()
    app = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'from ..domain.analysis_mode import AnalysisMode, WorkflowKind' in model
    assert 'from .model import StageRequirement, WorkflowPlan' in status
    assert 'from spinDecon.workflow.model import StageRequirement, build_workflow_plan' in overview
    assert 'from spinDecon.workflow.status import StageStatus, evaluate_workflow, recommended_action' in overview
    assert 'from spinDecon.gui.workspaces.workflow import WorkflowOverviewPanel' in app
    assert 'from ..workflow.registry import WORKFLOW_BY_KEY' in app


def test_active_tree_has_no_frames_imports():
    offenders = []
    for path in ROOT.rglob('*.py'):
        rel = path.relative_to(ROOT)
        if rel.parts[0] in {'Frames', 'legacy', 'tests'}:
            continue
        text = path.read_text(errors='ignore')
        if 'from spinDecon.Frames' in text or 'from .Frames' in text or 'import spinDecon.Frames' in text:
            offenders.append(str(rel))
    assert offenders == []
