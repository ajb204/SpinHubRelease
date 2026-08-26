#!/usr/bin/python
"""Standalone UniDecNMR window.

This frame hosts the existing deconvolution panel so the main NMR tab can stay
lightweight while the deconvolution controls live in their own window.
"""

from __future__ import annotations

import wx

from spinDecon.gui.workspaces.nmr import NMRWorkspace


class UniDecNMRFrame(wx.Frame):
    def __init__(self, controller, title: str = "UniDecNMR"):
        self.controller = controller
        parent = controller.GetTopLevelParent() if hasattr(controller, "GetTopLevelParent") else None
        super().__init__(parent, title=title, style=wx.DEFAULT_FRAME_STYLE | wx.RESIZE_BORDER)

        self.panel = NMRWorkspace(self, controller.deconParFile, state=controller.state, include_decon_box=True)
        # Keep the logical notebook parent so the existing callbacks continue to
        # operate on the main application tabs.
        self.panel.parent = controller.parent
        self.panel.unidecWindow = self

        self.controller.unidecWindow = self

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Fit()

        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _on_close(self, event):
        try:
            if getattr(self.controller, "unidecWindow", None) is self:
                self.controller.unidecWindow = None
        except Exception:
            pass
        try:
            self.panel.unidecWindow = None
        except Exception:
            pass
        event.Skip()
