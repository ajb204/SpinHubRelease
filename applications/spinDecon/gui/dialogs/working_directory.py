"""Small startup dialog for choosing spinDecon's project working directory.

Keep this module deliberately lightweight: it is shown before matplotlib,
scipy, the notebook UI, or processing modules are imported so an application
launch gives immediate visual feedback even on a cold Python start.
"""
from __future__ import annotations

from pathlib import Path
import wx


class WorkingDirectoryDialog(wx.Dialog):
    """Ask the user which directory should become the process working dir."""

    def __init__(self, parent, initial_directory: str | Path | None = None):
        super().__init__(parent, title="Welcome to spinDecon!")
        initial = Path(initial_directory or Path.home()).expanduser().resolve(strict=False)

        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)

        title = wx.StaticText(panel, label="Welcome to spinDecon!")
        font = title.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        font.SetPointSize(font.GetPointSize() + 2)
        title.SetFont(font)
        outer.Add(title, 0, wx.ALL, 14)

        message = wx.StaticText(
            panel,
            label="Choose the working directory for this spinDecon session.",
        )
        outer.Add(message, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 14)

        row = wx.BoxSizer(wx.HORIZONTAL)
        self.path_box = wx.TextCtrl(panel, value=str(initial), size=(480, -1))
        browse = wx.Button(panel, label="Choose...")
        browse.Bind(wx.EVT_BUTTON, self._choose_directory)
        row.Add(self.path_box, 1, wx.EXPAND | wx.RIGHT, 8)
        row.Add(browse, 0)
        outer.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 14)

        buttons = wx.StdDialogButtonSizer()
        open_button = wx.Button(panel, wx.ID_OK, label="Open")
        cancel_button = wx.Button(panel, wx.ID_CANCEL)
        buttons.AddButton(open_button)
        buttons.AddButton(cancel_button)
        buttons.Realize()
        outer.Add(buttons, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 14)

        panel.SetSizer(outer)
        dialog_sizer = wx.BoxSizer(wx.VERTICAL)
        dialog_sizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(dialog_sizer)
        self.Bind(wx.EVT_BUTTON, self._on_open, id=wx.ID_OK)
        self.Fit()
        self.CentreOnScreen()
        self.path_box.SetFocus()
        self.path_box.SetInsertionPointEnd()

    def _choose_directory(self, _event):
        current = self.path_box.GetValue().strip()
        default = Path(current).expanduser() if current else Path.home()
        if not default.is_dir():
            default = default.parent if default.parent.is_dir() else Path.home()
        dlg = wx.DirDialog(
            self,
            message="Choose spinDecon working directory",
            defaultPath=str(default),
            style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST,
        )
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.path_box.SetValue(dlg.GetPath())
        finally:
            dlg.Destroy()

    def _on_open(self, _event):
        text = self.path_box.GetValue().strip()
        path = Path(text).expanduser().resolve(strict=False) if text else None
        if path is None or not path.is_dir():
            wx.MessageBox(
                "Please choose an existing working directory.",
                "spinDecon",
                wx.OK | wx.ICON_WARNING,
                parent=self,
            )
            return
        self.EndModal(wx.ID_OK)

    def selected_directory(self) -> Path:
        return Path(self.path_box.GetValue().strip()).expanduser().resolve(strict=False)


def choose_working_directory(parent=None, initial_directory: str | Path | None = None) -> Path | None:
    """Return the chosen directory, or ``None`` when startup is cancelled."""
    dlg = WorkingDirectoryDialog(parent, initial_directory=initial_directory)
    try:
        if dlg.ShowModal() != wx.ID_OK:
            return None
        return dlg.selected_directory()
    finally:
        dlg.Destroy()
