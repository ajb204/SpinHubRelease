from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_project_service_has_no_wx_ownership():
    source = (ROOT / 'project' / 'service.py').read_text()
    assert 'import wx' not in source
    assert 'wx.App(' not in source
    assert 'wx.GetApp(' not in source
    assert 'from spinDecon.app.launcher import open_project_gui' in source


def test_app_launcher_owns_wx_application_lifecycle():
    source = (ROOT / 'app' / 'launcher.py').read_text()
    assert 'import wx' in source
    assert 'wx.GetApp()' in source
    assert 'wx.App(False)' in source
    assert 'wx_app.MainLoop()' in source
