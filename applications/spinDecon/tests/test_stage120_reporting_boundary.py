from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_project_summary_is_owned_by_gui_reporting():
    compat = (ROOT / 'project' / 'summary.py').read_text()
    canonical = (ROOT / 'gui' / 'reporting' / 'project_summary.py').read_text()
    nmr = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    assert 'from spinDecon.gui.reporting import project_summary as _canonical' in compat
    assert 'def generate_project_summary' in canonical
    assert 'from spinDecon.gui.reporting.project_summary import' in nmr
    assert 'from spinDecon.project.summary import' not in nmr
