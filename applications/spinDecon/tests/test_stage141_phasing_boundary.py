from pathlib import Path
from spinDecon.analysis.phasing_service import PhasingService
from spinDecon.app.context import ApplicationContext


def test_phasing_service_is_gui_independent():
    source = Path('analysis/phasing_service.py').read_text()
    assert 'import wx' not in source


def test_application_context_exposes_phasing_service_slot():
    assert hasattr(ApplicationContext(), 'phasing')


def test_phasing_workspace_no_longer_stores_tabone():
    source = Path('gui/workspaces/phasing.py').read_text()
    assert 'self.tabOne=' not in source
    assert 'parent.tabOne' not in source
