"""First-run project setup shown before the main Decon frame is constructed."""
from __future__ import annotations
from pathlib import Path
import os
import wx

from spinDecon.project.service import ProjectService, DEFAULT_PARAMETER_NAME
from spinDecon.processing.vpar_decon import inspect_acquisition, find_child_acquisitions


class ProjectSetupDialog(wx.Dialog):
    def __init__(self, parent, directory: str | Path):
        super().__init__(parent, title="Welcome to spinDecon!")
        self.working_dir = Path(directory).expanduser().resolve(strict=False)
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(panel, label="Welcome to spinDecon!")
        font = title.GetFont(); font.SetWeight(wx.FONTWEIGHT_BOLD); title.SetFont(font)
        outer.Add(title, 0, wx.ALL, 12)
        outer.Add(wx.StaticText(panel, label=f"Working directory: {self.working_dir}"), 0, wx.LEFT|wx.RIGHT|wx.BOTTOM, 12)
        outer.Add(wx.StaticText(panel, label="Please setup paths"), 0, wx.LEFT|wx.RIGHT|wx.BOTTOM, 12)

        grid = wx.FlexGridSizer(rows=2, cols=3, vgap=7, hgap=7); grid.AddGrowableCol(1, 1)
        self.raw_box = wx.TextCtrl(panel, value="", size=(380, -1), style=wx.TE_PROCESS_ENTER)
        self.spec_box = wx.TextCtrl(panel, value="./spec", size=(380, -1))
        for label, box, handler in (("Raw data:", self.raw_box, self._choose_raw), ("OutPath:", self.spec_box, self._choose_spec)):
            grid.Add(wx.StaticText(panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL); grid.Add(box, 1, wx.EXPAND)
            button=wx.Button(panel,label="...",size=(40,24)); button.Bind(wx.EVT_BUTTON,handler); grid.Add(button,0)
        outer.Add(grid, 0, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 12)

        summary_box = wx.StaticBoxSizer(wx.VERTICAL, panel, "Selected acquisition")
        # wxWidgets requires controls managed by a StaticBoxSizer to be
        # children of the StaticBox itself (not merely siblings on the panel).
        self.summary_text = wx.StaticText(summary_box.GetStaticBox(), label="Select a raw-data folder to inspect the acquisition.")
        self.summary_text.Wrap(520)
        summary_box.Add(self.summary_text, 0, wx.EXPAND|wx.ALL, 8)
        outer.Add(summary_box, 0, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 12)

        self.dim_box = wx.RadioBox(panel, label="Spectral dimensions", choices=["1D","2D","3D","4D"], majorDimension=4, style=wx.RA_SPECIFY_COLS)
        self.dim_box.SetSelection(0)
        self.pseudo_box = wx.CheckBox(panel, label="pseudoaxis")
        self.raw_box.Bind(wx.EVT_KILL_FOCUS, self._raw_path_changed)
        self.raw_box.Bind(wx.EVT_TEXT_ENTER, self._raw_path_changed)
        outer.Add(self.dim_box, 0, wx.LEFT|wx.RIGHT|wx.BOTTOM, 12); outer.Add(self.pseudo_box, 0, wx.LEFT|wx.RIGHT|wx.BOTTOM, 12)
        # Create the standard buttons with ``panel`` as their parent.  Using
        # Dialog.CreateStdDialogButtonSizer() creates buttons parented to the
        # dialog itself; putting that sizer inside the panel sizer triggers a
        # wxWidgets parent/sizer assertion on wxPython 4.3.
        buttons = wx.StdDialogButtonSizer()
        ok_button = wx.Button(panel, wx.ID_OK)
        cancel_button = wx.Button(panel, wx.ID_CANCEL)
        buttons.AddButton(ok_button)
        buttons.AddButton(cancel_button)
        buttons.Realize()
        outer.Add(buttons, 0, wx.EXPAND | wx.ALL, 12)
        panel.SetSizer(outer)

        # The dialog owns the panel; the panel owns every control managed by
        # ``outer``.  This keeps the wxWidgets containing-window hierarchy
        # consistent on macOS as well as the other platforms.
        dialog_sizer = wx.BoxSizer(wx.VERTICAL)
        dialog_sizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(dialog_sizer)
        self.Bind(wx.EVT_BUTTON, self._on_ok, id=wx.ID_OK)
        self.Fit()
        self.CentreOnParent()

    def _relative_to_project(self, chosen):
        chosen=os.path.abspath(chosen); root=str(self.working_dir)
        try: common=os.path.commonpath([chosen,root])
        except ValueError: common=''
        if common == root:
            value=os.path.relpath(chosen,root)
            return './' if value == '.' else './'+value.replace(os.sep,'/')
        return chosen

    def _choose_raw(self, event):
        dlg=wx.DirDialog(self,message="Choose raw data directory",defaultPath=self.raw_box.GetValue() or str(self.working_dir))
        if dlg.ShowModal()==wx.ID_OK:
            self.raw_box.SetValue(self._relative_to_project(dlg.GetPath()))
            self._inspect_raw_path()
        dlg.Destroy()

    def _raw_path_changed(self, event):
        self._inspect_raw_path()
        event.Skip()

    def _resolved_raw_path(self):
        value=self.raw_box.GetValue().strip()
        if not value: return None
        path=Path(value).expanduser()
        if not path.is_absolute(): path=self.working_dir/path
        return path.resolve(strict=False)

    def _inspect_raw_path(self):
        path=self._resolved_raw_path()
        if path is None:
            self.summary_text.SetLabel("Select a raw-data folder to inspect the acquisition."); return
        info=inspect_acquisition(path)
        if info is None:
            children=find_child_acquisitions(path)
            if children:
                rows=[f"No raw acquisition is present directly in this folder. {len(children)} acquisition(s) were found one level below:"]
                for child in children[:6]:
                    vendor='Bruker' if child.vendor=='bruk' else 'Varian'
                    extra=' + pseudo axis' if child.pseudo_axis else ''
                    seq=f" - {child.sequence}" if child.sequence else ''
                    rows.append(f"  {child.path.name}: {vendor}, {child.dimension}D{extra}{seq}")
                rows.append("Select one of these experiment folders as Raw data.")
                self.summary_text.SetLabel('\n'.join(rows))
            else:
                self.summary_text.SetLabel("No Bruker or Varian raw acquisition was detected in this folder.")
            self.summary_text.Wrap(520); self.Layout(); self.Fit(); return
        vendor='Bruker' if info.vendor=='bruk' else 'Varian'
        heading=f"{vendor} - detected {info.dimension}D" + (" + pseudo axis" if info.pseudo_axis else "")
        rows=[heading]
        if info.sequence: rows.append(f"Pulse sequence: {info.sequence}")
        if info.nuclei: rows.append("Nuclei: " + " / ".join(info.nuclei))
        if info.acquisition_time: rows.append(f"Acquisition time: {info.acquisition_time}")
        if info.pseudo_axis:
            desc=f"{info.pseudo_axis_size} points" if info.pseudo_axis_size else "detected"
            if info.pseudo_axis_columns: desc += " (" + ", ".join(info.pseudo_axis_columns) + ")"
            rows.append("Pseudo axis: " + desc)
        rows.append("Detection is best-effort; you can override the controls below.")
        self.summary_text.SetLabel('\n'.join(rows)); self.summary_text.Wrap(520)
        if info.dimension and 1 <= info.dimension <= 4: self.dim_box.SetSelection(info.dimension-1)
        self.pseudo_box.SetValue(bool(info.pseudo_axis))
        self.Layout(); self.Fit()

    def _choose_spec(self, event):
        current=self.spec_box.GetValue().strip() or './spec'; default=current if os.path.isabs(current) else self.working_dir/current
        dlg=wx.DirDialog(self,message="Choose output directory",defaultPath=os.path.abspath(default))
        if dlg.ShowModal()==wx.ID_OK: self.spec_box.SetValue(self._relative_to_project(dlg.GetPath()))
        dlg.Destroy()

    def _on_ok(self,event):
        if not self.raw_box.GetValue().strip():
            wx.MessageBox("Please select the raw data directory.","Project setup",wx.OK|wx.ICON_WARNING,parent=self); return
        self.EndModal(wx.ID_OK)

    def values(self):
        raw=self.raw_box.GetValue().strip(); raw_path=Path(raw).expanduser()
        if not raw_path.is_absolute(): raw_path=self.working_dir/raw_path
        return raw_path.resolve(strict=False), self.spec_box.GetValue().strip() or './spec', self.dim_box.GetSelection()+1, bool(self.pseudo_box.GetValue())


def run_project_setup(parent, *, service: ProjectService | None=None, directory: str | Path | None=None, parameter_name: str=DEFAULT_PARAMETER_NAME):
    service=service or ProjectService(); directory=Path(directory or Path.cwd()).resolve(strict=False)
    dlg=ProjectSetupDialog(parent,directory)
    try:
        if dlg.ShowModal()!=wx.ID_OK: return None
        raw,spec,dim,pseudo=dlg.values()
        state=service.create_initial_parameter_file(directory,raw,dimension=dim,pseudo_axis=pseudo,spec_path=spec,parameter_name=parameter_name)
        os.chdir(state.working_dir)
        return state
    finally: dlg.Destroy()
