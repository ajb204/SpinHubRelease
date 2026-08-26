from pathlib import Path


def test_nmrpipe_processing_module_imports_no_wx_at_module_scope():
    source = Path(__file__).parents[1].joinpath('processing', 'nmrpipe_scripts.py').read_text()
    assert '\nimport wx\n' not in source
    assert 'ShellOutputFrame' not in source
    assert 'wx.CallAfter' not in source


def test_gui_execution_owns_wx_processing_orchestration():
    source = Path(__file__).parents[1].joinpath('gui', 'dialogs', 'processing', 'execution.py').read_text()
    assert 'import wx' in source
    assert 'ShellOutputFrame' in source
    assert 'wx.CallAfter' in source


def test_nmr_workspace_uses_gui_nmrpipe_adapter():
    source = Path(__file__).parents[1].joinpath('gui', 'workspaces', 'nmr.py').read_text()
    assert 'gui.dialogs.processing.nmrpipe_adapter import nmrPipe' in source
