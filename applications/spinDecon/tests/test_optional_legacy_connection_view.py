"""Regression coverage for headless/modern workflows without legacy tabFour."""
from pathlib import Path


def _source():
    return (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()


def test_getconn_does_not_require_legacy_tabfour():
    source = _source()
    assert "tab_four = getattr(self.parent, 'tabFour', None)" in source
    assert "self.parent.tabFour.update_conn_data()" not in source


def test_decon_finish_does_not_require_legacy_tabfour():
    source = _source()
    assert "self.parent.tabFour.conn_data = self.conn_data" not in source
    assert "Analysis complete." in source
    assert "self.calcy.complete_decon_progress(True)" in source
