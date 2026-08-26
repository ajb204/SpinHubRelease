from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _assert_no_wx(package):
    for path in (ROOT / package).rglob('*.py'):
        source = path.read_text()
        assert 'import wx' not in source, path
        assert 'from wx' not in source, path


def test_project_package_is_gui_independent():
    _assert_no_wx('project')


def test_workflow_model_package_is_gui_independent():
    _assert_no_wx('workflow')
