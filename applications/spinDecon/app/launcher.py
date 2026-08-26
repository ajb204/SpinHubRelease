"""GUI launch boundary for opening an existing Decon project.

Keeping wx application ownership here allows :mod:`decon.project.service` to
remain usable by non-GUI project browsers and automation.
"""
from __future__ import annotations


def open_project_gui(parameter_file, *, state, workflow=None, show=True):
    import wx
    from spinDecon.app.notebook import MyApp

    wx_app = wx.GetApp()
    owns_app = wx_app is None
    if owns_app:
        wx_app = wx.App(False)

    frame = MyApp(str(parameter_file), showFlg=show, state=state)
    if workflow is not None:
        frame.notebook.open_workflow(workflow)
    if owns_app:
        wx_app.MainLoop()
    return frame
