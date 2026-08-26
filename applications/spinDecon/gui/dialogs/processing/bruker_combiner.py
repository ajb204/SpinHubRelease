"""GUI for combining numbered Bruker 1D experiments into one pseudo-2D SER."""
from pathlib import Path
import wx

from spinDecon.processing.bruker_combiner import (
    combine_bruker_experiments, discover_numbered_experiments, inspect_combination,
)


class CombineBrukerFrame(wx.Dialog):
    def __init__(self, parent, raw_dir):
        super().__init__(parent, title='Combine Bruker experiments', size=(820, 620),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.proc = parent
        self.raw_dir = Path(raw_dir)
        self.experiments = []

        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(wx.StaticText(panel, label=f'Fid path: {self.raw_dir}'), 0, wx.ALL | wx.EXPAND, 8)

        range_row = wx.BoxSizer(wx.HORIZONTAL)
        range_row.Add(wx.StaticText(panel, label='Start folder:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.start = wx.TextCtrl(panel, size=(90, -1))
        range_row.Add(self.start, 0, wx.RIGHT, 12)
        range_row.Add(wx.StaticText(panel, label='Finish folder:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.finish = wx.TextCtrl(panel, size=(90, -1))
        range_row.Add(self.finish, 0, wx.RIGHT, 12)
        self.discover_btn = wx.Button(panel, label='Discover')
        range_row.Add(self.discover_btn, 0)
        outer.Add(range_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        outer.Add(wx.StaticText(panel, label='Select experiments to combine:'), 0, wx.LEFT | wx.RIGHT | wx.TOP, 8)
        self.listbox = wx.CheckListBox(panel)
        outer.Add(self.listbox, 1, wx.ALL | wx.EXPAND, 8)

        outer.Add(wx.StaticText(panel, label='Inspection:'), 0, wx.LEFT | wx.RIGHT, 8)
        self.details = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.HSCROLL, size=(-1, 160))
        outer.Add(self.details, 0, wx.ALL | wx.EXPAND, 8)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.combine_btn = wx.Button(panel, label='Combine')
        cancel = wx.Button(panel, wx.ID_CANCEL, label='Cancel')
        buttons.AddStretchSpacer(1)
        buttons.Add(self.combine_btn, 0, wx.RIGHT, 8)
        buttons.Add(cancel, 0)
        outer.Add(buttons, 0, wx.ALL | wx.EXPAND, 8)
        panel.SetSizer(outer)

        self.discover_btn.Bind(wx.EVT_BUTTON, self.on_discover)
        self.listbox.Bind(wx.EVT_CHECKLISTBOX, self.on_selection_changed)
        self.combine_btn.Bind(wx.EVT_BUTTON, self.on_combine)
        self.on_discover(None)

    def _range(self):
        def parse(ctrl):
            text = ctrl.GetValue().strip()
            return int(text) if text else None
        return parse(self.start), parse(self.finish)

    def on_discover(self, event):
        try:
            start, finish = self._range()
            if start is not None and finish is not None and finish < start:
                raise ValueError('Finish folder must be greater than or equal to start folder.')
            self.experiments = discover_numbered_experiments(self.raw_dir, start, finish)
        except Exception as exc:
            wx.MessageBox(str(exc), 'Cannot discover experiments', wx.OK | wx.ICON_ERROR, self)
            return
        self.listbox.Clear()
        for exp in self.experiments:
            self.listbox.Append(f'{exp.number}    {exp.raw_kind}    {exp.raw_bytes} bytes    record {exp.record_bytes} bytes')
        for i in range(len(self.experiments)):
            self.listbox.Check(i, True)
        self._refresh_inspection()

    def _selected(self):
        return [exp for i, exp in enumerate(self.experiments) if self.listbox.IsChecked(i)]

    def on_selection_changed(self, event):
        self._refresh_inspection()
        event.Skip()

    def _refresh_inspection(self):
        selected = self._selected()
        info = inspect_combination(selected)
        lines = [f'Selected rows: {len(selected)}']
        if selected:
            lines.append('Experiments: ' + ', '.join(str(e.number) for e in selected))
        if info.numeric_varying_parameters:
            lines.append('')
            lines.append('Numeric pseudo-axis parameters detected:')
            for key, values in info.numeric_varying_parameters.items():
                preview = ', '.join(values[:8]) + (' ...' if len(values) > 8 else '')
                lines.append(f'  {key}: {preview}')
        if info.varying_parameters:
            other = [k for k in info.varying_parameters if k not in info.numeric_varying_parameters]
            if other:
                lines.append('')
                lines.append('Other varying scalar parameters (manifest only): ' + ', '.join(other))
        if info.warnings:
            lines.append('')
            lines.extend('WARNING: ' + x for x in info.warnings)
        if info.errors:
            lines.append('')
            lines.extend('ERROR: ' + x for x in info.errors)
        self.details.SetValue('\n'.join(lines))
        self.combine_btn.Enable(bool(selected) and not info.errors)

    def on_combine(self, event):
        selected = self._selected()
        if not selected:
            return
        existing = any((self.raw_dir / name).exists() for name in ('ser', 'acqus', 'acqu2s', 'pulseprogram', 'combine_manifest.json'))
        overwrite = False
        if existing:
            answer = wx.MessageBox(
                'Combined raw-data files already exist in the fid path. Replace them?\n\n'
                'Numbered source folders will not be modified.',
                'Replace combined data?', wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING, self)
            if answer != wx.YES:
                return
            overwrite = True
        try:
            combine_bruker_experiments(self.raw_dir, selected, overwrite=overwrite)
        except Exception as exc:
            wx.MessageBox(str(exc), 'Combination failed', wx.OK | wx.ICON_ERROR, self)
            return
        wx.MessageBox(
            f'Combined {len(selected)} experiments into {self.raw_dir / "ser"}.\n'
            'Synthetic acqus/acqu2s, pseudo-axis list files and combine_manifest.json were created.',
            'Combination complete', wx.OK | wx.ICON_INFORMATION, self)
        self.EndModal(wx.ID_OK)
