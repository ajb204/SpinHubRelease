import os
import logging
import re
from pathlib import Path

import wx

from spinDecon.gui.dialogs.errors import errorMessage
from spinDecon.gui.dialogs.shell_output import run_command_with_output
from spinDecon.processing.vpar_decon import vpar,GetParBruk
from spinDecon.processing.dimension_contract import legacy_vpar_dimension



class ConversionFrame(wx.Frame):
    """Standalone conversion window for Bruker/Varian/GE fid->nmrPipe script generation."""

    def __init__(self, parent):
        super().__init__(parent, title='Conversion', size=(720, 540))
        self.proc = parent
        self.script_path = None
        self.script_frame = None
        self.status_bar = None
        self._hover_default_status = 'Ready'
        self.inst = None

        self.dim_choices = getattr(parent, 'dim_choices', ['H1', 'N15', 'C13', 'F19', 'P31'])
        self.refList = ['Water', 'Auto', 'Manual']

        self._build_ui()
        self._load_from_file()
        # Seed only missing keys: existing live edits are authoritative.
        state = getattr(self.proc, 'state', None)
        if state is not None:
            state.seed_gui_settings(self.collect_updates(update_state=False))
        self.apply_live_settings()
        self._bind_live_state_controls()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _bind_live_state_controls(self):
        """Mirror conversion widget edits into shared live GUI state."""
        def changed(event):
            try:
                self.collect_updates()
            except Exception:
                pass
            try:
                event.Skip()
            except Exception:
                pass

        names = ['RealName', 'xcenBox', 'nusFil', 'refBox']
        names += [f'label{i}' for i in range(1, 5)]
        names += [f'cb_rk{i}' for i in range(1, 4)]
        seen = set()
        for name in names:
            ctrl = getattr(self, name, None)
            if ctrl is None or id(ctrl) in seen:
                continue
            seen.add(id(ctrl))
            for evt in (wx.EVT_TEXT, wx.EVT_CHOICE, wx.EVT_CHECKBOX):
                try:
                    ctrl.Bind(evt, changed)
                except Exception:
                    pass

    def apply_live_settings(self):
        state = getattr(self.proc, 'state', None)
        live = getattr(state, 'gui_settings', {}) if state is not None else {}
        for idx in range(1, min(getattr(self, 'label_count', 0), 4) + 1):
            key = f'label{idx}'; ctrl = getattr(self, key, None)
            if key in live and ctrl is not None:
                try: ctrl.SetValue(str(live[key]))
                except Exception: pass
        for key, ctrl_name in [('RealName','RealName'), ('xcen','xcenBox'), ('nusFil','nusFil')]:
            ctrl = getattr(self, ctrl_name, None)
            if key in live and ctrl is not None:
                try: ctrl.SetValue(str(live[key]))
                except Exception: pass
        if self._has_pseudo_axis() and hasattr(self, 'RealName'):
            try:
                if not self._clean_label_text(self.RealName.GetValue()):
                    self.RealName.SetValue('pseudo')
            except Exception:
                pass
        if 'refBox' in live and hasattr(self, 'refBox'):
            try: self.refBox.SetSelection(int(live['refBox']))
            except Exception: pass
        for idx in range(1, 4):
            key=f'rk{idx}'; ctrl=getattr(self, f'cb_rk{idx}', None)
            if key in live and ctrl is not None:
                try: ctrl.SetValue(str(live[key]).lower() in ('y','yes','1','true'))
                except Exception: pass

    # ------------------------------------------------------------------
    # UI


    def _set_status_text(self, text):
        try:
            if getattr(self, 'status_bar', None) is not None:
                self.status_bar.SetStatusText(str(text or ''))
        except Exception:
            pass

    def _hover_lookup_text(self, widget):
        hover_map = getattr(self, '_hover_help_map', {})
        current = widget
        while current is not None:
            try:
                text = hover_map.get(current)
            except Exception:
                text = None
            if text:
                return text
            try:
                current = current.GetParent()
            except Exception:
                current = None
        return self._hover_default_status

    def _on_hover_event(self, event):
        try:
            widget = event.GetEventObject()
        except Exception:
            widget = None
        self._set_status_text(self._hover_lookup_text(widget))
        try:
            event.Skip()
        except Exception:
            pass

    def _on_hover_leave_frame(self, event):
        self._set_status_text(self._hover_default_status)
        try:
            event.Skip()
        except Exception:
            pass

    def _install_hover_help(self, widget, text):
        if widget is None:
            return
        try:
            if not hasattr(self, '_hover_help_map'):
                self._hover_help_map = {}
            self._hover_help_map[widget] = text
        except Exception:
            pass
        try:
            widget.SetToolTip(text)
        except Exception:
            pass
        try:
            widget.Bind(wx.EVT_ENTER_WINDOW, self._on_hover_event)
            widget.Bind(wx.EVT_MOTION, self._on_hover_event)
        except Exception:
            pass

    def _install_hover_map(self, mapping):
        for item in mapping:
            if not item:
                continue
            widget, text = item
            self._install_hover_help(widget, text)


    def _hover_widget_name(self, widget):
        for name, value in self.__dict__.items():
            if value is widget:
                return name
        return ''

    def _guess_hover_text(self, widget, name=''):
        name_l = str(name or '').lower()
        label = ''
        try:
            label = str(widget.GetLabel()).strip()
        except Exception:
            label = ''
        label_l = label.lower()

        label_map = {
            'p0': 'Zero-order phase correction.',
            'p1': 'First-order phase correction.',
            'lp': 'Linear prediction.',
            'poly': 'Polynomial baseline correction.',
            'f1180': 'Apply first-point correction.',
            'flip': 'Signal polarity / flip handling.',
            'window': 'Select the apodization window function.',
            'op1': 'Primary parameter for the selected window function.',
            'op2': 'Secondary parameter for the selected window function.',
            'first pt': 'First-point scaling factor used by the window.',
            'direct dimension:': 'Settings for the direct dimension.',
            'dimensions:': 'Configure nucleus labels for each dimension.',
            'referencing:': 'Referencing settings.',
            'sparse sampling:': 'Sparse-sampling schedule settings.',
            'status': 'Status summary.',
            'dimension:': 'Current data dimensionality.',
            'format:': 'Detected input data format.',
            'converted': 'Converted output status indicator.',
            'processed': 'Processed output status indicator.',
            'script target': 'Choose the target script mode.',
            'smile maxiter:': 'Maximum number of SMILE iterations.',
            'outpath:': 'Working output directory.',
            'specpath:': 'Processed spectrum output directory.',
            'fid select:': 'Choose which FID trace or slice is used for the preview.',
            'x ranges (ppm):': 'Set the displayed ppm window.',
            'min:': 'Lower ppm limit of the display window.',
            'max:': 'Upper ppm limit of the display window.',
        }
        if label_l in label_map:
            return label_map[label_l]

        class_name = type(widget).__name__.lower()

        if 'showscript' in name_l:
            return 'Preview the generated script.'
        if 'guessbtn' in name_l:
            return 'Guess the processing settings from the current file.'
        if 'runbtn' in name_l:
            return 'Run the current processing or conversion step.'
        if 'savebtn' in name_l:
            return 'Save the current settings.'
        if 'closebtn' in name_l:
            return 'Close this window.'
        if 'buttonconv' in name_l:
            return 'Open the conversion settings.'
        if 'buttonproc' in name_l:
            return 'Open the processing settings.'
        if 'buttonstore' in name_l:
            return 'Choose the output directory.'
        if 'opendirfilebtn' in name_l:
            return 'Choose the output directory.'
        if 'openspecdirbtn' in name_l:
            return 'Choose the spectrum output folder.'
        if 'nusbrowse' in name_l:
            return 'Browse for a sparse-sampling schedule.'
        if 'sld_0_auto_btn' in name_l or 'sld_1_auto_btn' in name_l:
            return 'Estimate the phase automatically from the current signal.'
        if 'sld_0_mode_btn' in name_l:
            return 'Switch between coarse and fine P0 control.'
        if 'sld_1_mode_btn' in name_l:
            return 'Switch between coarse and fine P1 control.'
        if 'scriptmodebox' in name_l:
            return 'Choose whether the generated script targets standard processing or SMILE.'
        if 'maxiterbox' in name_l:
            return 'Maximum number of SMILE iterations.'
        if 'refbox' in name_l:
            return 'Choose the referencing mode.'
        if 'xcenbox' in name_l:
            return 'Enter the manual reference position in ppm.'
        if 'realname' in name_l:
            return 'Enter the real-axis label used for pseudo-dimension data.'
        if 'nusfil' in name_l:
            return 'Path to the sparse-sampling schedule file.'
        if 'path_label' in name_l:
            return 'Current processing script path.'
        if 'scriptbox' in name_l:
            return 'Editable processing script text.'
        if 'convertedautobtn' in name_l:
            return 'Automatically refresh the converted status.'
        if 'processedautobtn' in name_l:
            return 'Automatically refresh the processed status.'
        if 'convertedlamp' in name_l:
            return 'Converted output status light.'
        if 'processedlamp' in name_l:
            return 'Processed output status light.'
        if 'cb_show_fid' in name_l:
            return 'Display the raw FID instead of the phased spectrum.'
        if 'fidselect' in name_l:
            return 'Choose which FID trace or slice is used for the preview.'
        if 'xminbox' in name_l:
            return 'Lower ppm limit of the display window.'
        if 'xmaxbox' in name_l:
            return 'Upper ppm limit of the display window.'

        if 'check' in class_name:
            if 'base' in name_l and 'lin' in name_l:
                return 'Enable linear baseline correction.'
            if 'base' in name_l and 'sol' in name_l:
                return 'Enable digital solvent suppression.'
            if 'base' in name_l and 'pol' in name_l:
                return 'Apply polynomial baseline correction.'
            if 'f1180' in name_l:
                return 'Apply first-point correction.'
            if 'lp' in name_l:
                return 'Enable linear prediction.'
            if 'rk' in name_l:
                return 'Enable the Rance-Kay option for this indirect dimension.'
            if 'show_fid' in name_l:
                return 'Display the raw FID instead of the phased spectrum.'
            return 'Toggle this option.'

        if 'togglebutton' in class_name:
            if 'auto' in label_l:
                return 'Automatically refresh the related status value.'
            if label_l == 'c':
                return 'Switch between coarse and fine control.'
            return 'Toggle this option.'

        if 'radiobox' in class_name:
            return 'Choose one of the available modes.'

        if 'slider' in class_name:
            if '0' in name_l:
                return 'Adjust the zero-order phase.'
            if '1' in name_l:
                return 'Adjust the first-order phase.'
            return 'Adjust the current value.'

        if 'combobox' in class_name:
            if 'windowbox' in name_l:
                return 'Choose the apodization window function.'
            if 'cb_ft' in name_l:
                return 'Choose the flip mode for this dimension.'
            if 'refbox' in name_l:
                return 'Choose the referencing mode.'
            return 'Choose a value from the list.'

        if 'textctrl' in class_name:
            if 'p0' in name_l:
                return 'Enter the zero-order phase correction.'
            if 'p1' in name_l:
                return 'Enter the first-order phase correction.'
            if 'firstpoint' in name_l:
                return 'Enter the first-point scaling factor.'
            if 'win2val' in name_l:
                return 'Enter the primary window parameter.'
            if 'win3val' in name_l:
                return 'Enter the secondary window parameter.'
            if 'label' in name_l:
                return 'Enter the nucleus label for this dimension.'
            if 'maxiterbox' in name_l:
                return 'Maximum number of SMILE iterations.'
            if 'xcenbox' in name_l:
                return 'Enter the manual reference position in ppm.'
            if 'nusfil' in name_l:
                return 'Path to the sparse-sampling schedule file.'
            if 'scriptbox' in name_l:
                return 'Editable processing script text.'
            if 'path' in name_l:
                return 'Current file or folder path.'
            return 'Enter a value.'

        if 'button' in class_name:
            if label_l == '...':
                if 'spec' in name_l:
                    return 'Choose the spectrum output folder.'
                if 'dir' in name_l:
                    return 'Choose the output directory.'
                return 'Browse for a file or folder.'
            if label_l == 'autophase':
                return 'Estimate the phase automatically from the current signal.'
            if label_l == 'load':
                return 'Load the sparse-sampling schedule file.'
            if label_l == 'advanced':
                if 'conv' in name_l:
                    return 'Open advanced conversion settings.'
                if 'proc' in name_l:
                    return 'Open advanced processing settings.'
                return 'Open advanced settings.'
            if label_l == 'datastore':
                return 'Choose the output directory.'
            if label_l == 'guess':
                return 'Guess the processing settings from the current file.'
            return f'{label or "Button"} action.'

        if 'statictext' in class_name:
            if label_l == 'p0:':
                return 'Zero-order phase correction.'
            if label_l == 'p1:':
                return 'First-order phase correction.'
            if label_l == 'f1180':
                return 'Apply first-point correction.'
            if label_l == 'lp':
                return 'Linear prediction.'
            if label_l == 'poly':
                return 'Polynomial baseline correction.'
            if label_l == 'flip':
                return 'Signal polarity / flip handling.'
            if label_l == 'window':
                return 'Apodization window function.'
            if label_l == 'op1':
                return 'Primary parameter for the selected window function.'
            if label_l == 'op2':
                return 'Secondary parameter for the selected window function.'
            if label_l == 'first pt':
                return 'First-point scaling factor.'
            if label_l.startswith('dim '):
                return f'Processing settings for {label_l}.'
            if label_l.startswith('label '):
                return 'Enter the nucleus label for this dimension.'
            if label_l == 'direct dimension:':
                return 'Settings for the direct dimension.'
            if label_l == 'referencing:':
                return 'Referencing settings.'
            if label_l == 'sparse sampling:':
                return 'Sparse-sampling schedule settings.'
            if label_l == 'status':
                return 'Status summary.'
        return ''

    def _install_default_hover_help(self):
        try:
            items = list(self.__dict__.items())
        except Exception:
            items = []
        for name, widget in items:
            if widget is None:
                continue
            if name.startswith('_'):
                continue
            if name in {'panel', 'status_bar', 'statusBox', 'statusSizer', 'statusPanel', 'canvas', 'toolbar'}:
                continue
            try:
                if not isinstance(widget, wx.Window):
                    continue
            except Exception:
                continue
            try:
                if hasattr(self, '_hover_help_map') and widget in self._hover_help_map:
                    continue
            except Exception:
                pass
            text = self._guess_hover_text(widget, name)
            if not text:
                continue
            self._install_hover_help(widget, text)

    def _build_ui(self):
        panel = wx.Panel(self)
        self.panel = panel
        panel.SetBackgroundColour(wx.Colour(255, 255, 255))

        root = wx.BoxSizer(wx.VERTICAL)

        self.label_count = self._dimension_count()
        self.pseudo_axis = self._has_pseudo_axis()

        self.label_controls = []
        self.rk_controls = []
        self.rk_rows = []

        def add_label_row(parent_sizer, index, with_rk=False):
            row = wx.BoxSizer(wx.HORIZONTAL)
            box_parent = parent_sizer.GetStaticBox()
            label_text = wx.StaticText(box_parent, label=f'Label {index}')
            label_ctrl = wx.TextCtrl(box_parent, size=(100, 22))
            setattr(self, f'label{index}_text', label_text)
            setattr(self, f'label{index}', label_ctrl)
            self.label_controls.append(label_ctrl)
            row.Add(label_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.Add(label_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)

            if with_rk:
                rk_index = index - 1
                rk_text = wx.StaticText(box_parent, label=f'Rance-Kay {rk_index}')
                rk_ctrl = wx.CheckBox(box_parent, -1, '')
                setattr(self, f'rk{rk_index}_text', rk_text)
                setattr(self, f'cb_rk{rk_index}', rk_ctrl)
                self.rk_controls.append(rk_ctrl)
                self.rk_rows.append((rk_text, rk_ctrl))
                row.Add(rk_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
                row.Add(rk_ctrl, 0, wx.ALIGN_CENTER_VERTICAL)

            parent_sizer.Add(row, 0, wx.BOTTOM, 6)

        # ---- Dimensions ---------------------------------------------------------
        dimensions_box = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Dimensions:'), wx.VERTICAL)
        add_label_row(dimensions_box, 1, with_rk=False)
        for idx in range(2, min(self.label_count, 4) + 1):
            add_label_row(dimensions_box, idx, with_rk=True)

        if self.pseudo_axis:
            self.real_axis_row = wx.BoxSizer(wx.HORIZONTAL)
            self.real_axis_text = wx.StaticText(dimensions_box.GetStaticBox(), label='Real axis')
            self.RealName = wx.TextCtrl(dimensions_box.GetStaticBox(), size=(100, 22))
            self.real_axis_row.Add(self.real_axis_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            self.real_axis_row.Add(self.RealName, 0, wx.ALIGN_CENTER_VERTICAL)
            dimensions_box.Add(self.real_axis_row, 0, wx.BOTTOM, 0)

        # ---- Referencing --------------------------------------------------------
        referencing_static = wx.StaticBox(panel, -1, 'Referencing:')
        referencing_box = wx.StaticBoxSizer(referencing_static, wx.VERTICAL)
        self.refBox = wx.RadioBox(
            referencing_static,
            label='',
            choices=self.refList,
            majorDimension=1,
            style=wx.RA_SPECIFY_ROWS,
        )
        self.Bind(wx.EVT_RADIOBOX, self._on_reference_change, self.refBox)

        self.xcenBox = wx.TextCtrl(referencing_static, size=(80, 22))
        referencing_box.Add(self.refBox, 0, wx.EXPAND | wx.ALL, 8)

        manual_row = wx.BoxSizer(wx.HORIZONTAL)
        manual_row.Add(wx.StaticText(referencing_static, label='Manual ppm'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        manual_row.Add(self.xcenBox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.xcarStatus = wx.StaticText(referencing_static, label='xCar ? ppm')
        manual_row.Add(self.xcarStatus, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 12)
        self.xcenBox.Bind(wx.EVT_TEXT, self._on_manual_ppm_change)
        referencing_box.Add(manual_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        # ---- Sparse sampling ----------------------------------------------------
        sparse_static = wx.StaticBox(panel, -1, 'Sparse sampling:')
        sparse_box = wx.StaticBoxSizer(sparse_static, wx.VERTICAL)
        self.nusLabel = wx.StaticText(sparse_static, label='Schedule')
        self.nusFil = wx.TextCtrl(sparse_static, size=(380, 22))
        self.nusBrowse = wx.Button(sparse_static, label='Load')
        self.nusBrowse.Bind(wx.EVT_BUTTON, self._browse_nus)

        schedule_row = wx.BoxSizer(wx.HORIZONTAL)
        schedule_row.Add(self.nusLabel, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        schedule_row.Add(self.nusFil, 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        schedule_row.Add(self.nusBrowse, 0, wx.ALIGN_CENTER_VERTICAL)
        sparse_box.Add(schedule_row, 0, wx.EXPAND | wx.ALL, 8)

        # ---- Layout -------------------------------------------------------------
        top_row = wx.BoxSizer(wx.HORIZONTAL)
        top_row.Add(dimensions_box, 1, wx.EXPAND | wx.RIGHT, 10)
        top_row.Add(referencing_box, 0, wx.EXPAND)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.showScriptBtn = wx.Button(panel, label='Show Script')
        self.runBtn = wx.Button(panel, label='Run')
        self.closeBtn = wx.Button(panel, label='Close')
        btn_row.Add(self.showScriptBtn, 0, wx.RIGHT, 6)
        btn_row.Add(self.runBtn, 0, wx.RIGHT, 6)
        btn_row.Add(self.closeBtn, 0)

        self.showScriptBtn.Bind(wx.EVT_BUTTON, self.on_show_script)
        self.runBtn.Bind(wx.EVT_BUTTON, self.on_run)
        self.closeBtn.Bind(wx.EVT_BUTTON, self._on_close)

        self._install_hover_map([
            (getattr(self, 'label1_text', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label1', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label2_text', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label2', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label3_text', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label3', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label4_text', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'label4', None), 'Enter the nucleus label for this dimension.'),
            (getattr(self, 'RealName', None), 'Enter the real-axis label used for pseudo-dimension data.'),
            (self.refBox, 'Choose the referencing mode.'),
            (self.xcenBox, 'Enter the manual reference position in ppm.'),
            (self.xcarStatus, 'Calculated reference position; updates as the reference mode changes.'),
            (self.nusLabel, 'Path to the sparse-sampling schedule file.'),
            (self.nusFil, 'Path to the sparse-sampling schedule file.'),
            (self.nusBrowse, 'Browse for a sparse-sampling schedule.'),
            (self.showScriptBtn, 'Preview the generated conversion script.'),
            (self.runBtn, 'Run the conversion.'),
            (self.closeBtn, 'Close the conversion window.'),
        ])

        root.Add(top_row, 0, wx.ALL | wx.EXPAND, 10)
        root.Add(sparse_box, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 10)
        root.Add(btn_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.ALIGN_RIGHT, 10)

        panel.SetSizer(root)
        panel.Layout()

        try:
            self.status_bar = self.CreateStatusBar(1)
        except Exception:
            self.status_bar = None

        panel.Fit()
        best = panel.GetBestSize()
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_hover_leave_frame)
        self.SetMinSize(best)
        self.SetClientSize(best)
        self.Layout()
        try:
            self.Fit()
        except Exception:
            pass
        try:
            self._set_status_bar_text(self._reference_status_text())
            self._install_default_hover_help()
        except Exception:
            pass
    # ------------------------------------------------------------------
    # Helpers
    def _controls_for_dim(self):
        spectral = max(1, int(self.proc.spectral_dim_count))
        return spectral, max(0, spectral - 1)

    def _combo_items(self, ctrl):
        items = []
        self._stdout_log(f'combo item reader widget={type(ctrl).__name__}')

        readers = []
        for attr in ('GetStrings', 'GetItems'):
            if hasattr(ctrl, attr):
                readers.append(attr)
        if hasattr(ctrl, 'GetCount') and hasattr(ctrl, 'GetString'):
            readers.append('GetCount/GetString')

        for reader in readers:
            try:
                if reader == 'GetCount/GetString':
                    count = int(ctrl.GetCount())
                    items = [ctrl.GetString(i) for i in range(count)]
                else:
                    items = list(getattr(ctrl, reader)())
                self._stdout_log(f'combo items read via {reader}: {items!r}')
                return items
            except Exception as exc:
                self._stdout_log(f'combo item read via {reader} failed: {exc!r}')
                pass

        try:
            current = ctrl.GetValue()
        except Exception as exc:
            self._stdout_log(f'combo value fallback read failed: {exc!r}')
            return []
        if current:
            items = [current]
        self._stdout_log(f'combo items fallback using current value only: {items!r}')
        return items

    def _replace_combo_control(self, attr_name, old_ctrl, new_choices, value):
        parent = old_ctrl.GetParent()
        sizer = old_ctrl.GetContainingSizer()
        try:
            style = old_ctrl.GetWindowStyle()
        except Exception:
            style = 0
        try:
            size = old_ctrl.GetSize()
        except Exception:
            size = wx.DefaultSize
        pass

        new_ctrl = wx.ComboBox(parent, size=size, choices=new_choices, style=style)
        try:
            new_ctrl.SetValue(value)
        except Exception:
            pass

        if sizer is not None:
            replaced = False
            try:
                replaced = bool(sizer.Replace(old_ctrl, new_ctrl))
                self._stdout_log(f'combo sizer Replace returned {replaced}')
            except Exception as exc:
                self._stdout_log(f'combo sizer Replace failed: {exc!r}')
            if not replaced:
                try:
                    item = sizer.GetItem(old_ctrl)
                    children = sizer.GetChildren()
                    idx = children.index(item) if item in children else -1
                    flag = item.GetFlag() if item is not None else 0
                    border = item.GetBorder() if item is not None else 0
                    proportion = item.GetProportion() if item is not None else 0
                    sizer.Detach(old_ctrl)
                    if idx >= 0:
                        sizer.Insert(idx, new_ctrl, proportion, flag, border)
                    else:
                        sizer.Add(new_ctrl, proportion, flag, border)
                    replaced = True
                    self._stdout_log('combo sizer fallback reinsert succeeded')
                except Exception as exc2:
                    self._stdout_log(f'combo sizer fallback reinsert failed: {exc2!r}')

        try:
            old_ctrl.Destroy()
        except Exception:
            pass
        setattr(self, attr_name, new_ctrl)
        if hasattr(self, '_hover_help_map'):
            try:
                hover_text = self._hover_help_map.get(old_ctrl)
            except Exception:
                hover_text = None
            if hover_text:
                self._install_hover_help(new_ctrl, hover_text)
        if hasattr(self, 'label_controls') and attr_name.startswith('label'):
            try:
                idx = int(attr_name.replace('label', '')) - 1
                if 0 <= idx < len(self.label_controls):
                    self.label_controls[idx] = new_ctrl
            except Exception:
                pass
        try:
            parent.Layout()
            parent.Fit()
            self.Layout()
        except Exception:
            pass
        return new_ctrl

    def _safe_set_combo(self, attr_name, ctrl, value):
        if value is None:
            return ctrl
        value = str(value).strip()
        if not value:
            return ctrl
        self._stdout_log(
            f'combo update start attr={attr_name}, widget={type(ctrl).__name__}, '
            f'has_SetItems={hasattr(ctrl, "SetItems")}, has_Append={hasattr(ctrl, "Append")}, '
            f'has_SetStringSelection={hasattr(ctrl, "SetStringSelection")}, incoming={value!r}'
        )
        pass

        choices = self._combo_items(ctrl)
        self._stdout_log(f'combo before update choices={choices!r}, incoming={value!r}')
        pass

        try:
            if value not in choices:
                new_choices = list(choices) + [value]
                self._stdout_log(f'combo missing {value!r}; new choices will be {new_choices!r}')
                applied = False
                try:
                    if hasattr(ctrl, 'SetItems'):
                        ctrl.SetItems(new_choices)
                        applied = True
                        self._stdout_log('combo updated via SetItems')
                except Exception as exc:
                    self._stdout_log(f'combo SetItems failed: {exc!r}')
                    pass
                if not applied:
                    try:
                        if hasattr(ctrl, 'Clear'):
                            ctrl.Clear()
                        for item in new_choices:
                            ctrl.Append(item)
                        applied = True
                        self._stdout_log('combo updated via Clear/Append')
                    except Exception as exc:
                        self._stdout_log(f'combo Clear/Append failed: {exc!r}')
                        pass
                if not applied:
                    ctrl = self._replace_combo_control(attr_name, ctrl, new_choices, value)
                    self._stdout_log('combo replaced with a new control instance')
        except Exception:
            pass

        selected = False
        try:
            if hasattr(ctrl, 'SetStringSelection'):
                selected = bool(ctrl.SetStringSelection(value))
                self._stdout_log(f'SetStringSelection({value!r}) -> {selected}')
        except Exception as exc:
            self._stdout_log(f'SetStringSelection({value!r}) failed: {exc!r}')
            selected = False
        if not selected:
            try:
                ctrl.SetValue(value)
                self._stdout_log(f'SetValue({value!r}) applied')
            except Exception as exc:
                self._stdout_log(f'SetValue({value!r}) failed: {exc!r}')
        try:
            if hasattr(ctrl, 'FindString') and hasattr(ctrl, 'SetSelection'):
                idx = ctrl.FindString(value)
                self._stdout_log(f'FindString({value!r}) -> {idx}')
                if idx != -1:
                    ctrl.SetSelection(idx)
                    self._stdout_log(f'SetSelection({idx}) applied')
        except Exception as exc:
            self._stdout_log(f'SetSelection via FindString failed: {exc!r}')
        try:
            ctrl.Refresh()
            ctrl.Update()
        except Exception as exc:
            self._stdout_log(f'combo refresh/update failed: {exc!r}')
        after = self._combo_items(ctrl)
        try:
            current = ctrl.GetValue()
        except Exception as exc:
            current = f'<error reading value: {exc!r}>'
        self._stdout_log(f'combo after update choices={after!r}, selected={current!r}')
        pass
        return ctrl

    def _safe_set_text(self, ctrl, value):
        if value is None:
            return
        try:
            ctrl.SetValue(str(value))
        except Exception:
            pass

    def _clean_label_text(self, value):
        text = '' if value is None else str(value)
        text = text.strip().replace(' ', '').replace('<', '').replace('>', '')
        if text in ('', '0', 'None'):
            return ''
        return text

    def _default_label_text(self, idx):
        # Do not invent Varian nucleus labels when procpar is ambiguous.
        if self._spectrometer_type() == 'var':
            return ''
        return f'H1_{idx}'

    def _stdout_log(self, message):
        pass

    def _data_path_candidates(self):
        candidates = []
        for candidate in (
            self.proc._raw_output_dir(),
        ):
            if candidate and candidate not in candidates:
                candidates.append(candidate)
        self._stdout_log(f'label fallback path candidates: {candidates}')
        #logging.debug('Conversion label fallback path candidates: %s', candidates)
        return candidates

    def _bruker_label_default(self, idx):
        if not self._is_bruker():
            #msg = f'label{idx} Bruker fallback skipped because data type is not bruker'
            #print(msg)
            #print('Conversion label fallback skipped for label%d because data type is not bruker' % idx)
            return ''
        dim = self._dimension_count()
        if idx < 1 or idx > min(dim, 4):
            msg = f'label{idx} Bruker fallback skipped because dimension={dim}'
            #print(msg)
            #print('Conversion label fallback skipped for label%d because dimension=%s', idx, dim)
            return ''
        acqu_files = {
            1: 'acqus',
            2: 'acqu2s',
            3: 'acqu3s',
            4: 'acqu4s',
        }
        root = self._root()
        target_name = acqu_files.get(idx, 'acqus')
        #print(f'Attempting Bruker fallback for label{idx}: target file {target_name} (dim={dim})')
        for base_path in self._data_path_candidates():
            acqu_path = os.path.join(base_path, target_name)
            exists = os.path.exists(acqu_path)
            #print(f'Checking {acqu_path}: exists={exists}')
            #print('Conversion label fallback: looking for label%d in %s' %( idx, acqu_path))
            if not exists:
                #print(f'{acqu_path} missing; skipping')
                continue

            try:
                #print(f'Reading NUC1 from {acqu_path}')
                raw=GetParBruk(acqu_path,('','NUC1'),verb='n')[0]
                #raw = root.Parse(acqu_path, 'NUC1', default='')
                #print(f'NUC1 from {acqu_path}: {raw!r}')
            except Exception:
                #print(f'Parse failed for {acqu_path}')
                #print('Conversion label fallback: parse failed for %s' % acqu_path)
                raw = ''
            cleaned = self._clean_label_text(raw)
            #print(f'label{idx} raw={raw!r} cleaned={cleaned!r} from {acqu_path}')
            #print('Conversion label fallback: label%d raw=%r cleaned=%r from %s'%( idx, raw, cleaned, acqu_path))
            if cleaned:
                #print(f'label{idx} resolved to {cleaned!r}')
                return cleaned
        #print(f'No NUC1 value found for label{idx} in any candidate path')
        #print('Conversion label fallback: no NUC1 found for label%d', idx)
        return ''

    def _dimension_count(self):
        """Number of spectral label controls.

        Conversion has one control per physical axis: spectral labels are
        represented by labelN controls and the optional real physical axis by
        RealName.  The two counts must therefore not be conflated.
        """
        return max(1, int(getattr(self.proc, 'spectral_dim_count', getattr(self.proc, 'dim', 1))))

    def _has_pseudo_axis(self):
        return bool(getattr(self.proc, 'has_pseudo_axis', False))

    def _physical_dimension_count(self):
        return int(getattr(self.proc, 'physical_dim_count', self._dimension_count() + int(self._has_pseudo_axis())))

    def _legacy_backend_dim(self):
        """Temporary adapter for the pre-topology vpar backend only.

        No GUI window consumes this encoding.  It can be deleted when vpar is
        migrated to DatasetTopology.
        """
        return legacy_vpar_dimension(self.proc.topology)

    def _update_widget_visibility(self):
        is_bruker = self._is_bruker()
        for row in getattr(self, 'rk_rows', []):
            for widget in row:
                try:
                    widget.Show(not is_bruker)
                except Exception:
                    pass
        try:
            text = self._reference_status_text()
            if hasattr(self, 'xcarStatus'):
                self.xcarStatus.SetLabel(text)
            self._set_status_bar_text(text)
        except Exception:
            pass
        try:
            if hasattr(self, 'panel'):
                self.panel.Layout()
            self.Layout()
            self.Fit()
        except Exception:
            pass

    def _root(self):
        return getattr(self.proc, 'parent', self.proc)

    def _outpath_dir(self):
        root = self._root()
        for attr in ('outPathBox', 'dirCtrl'):
            ctrl = getattr(root, attr, None)
            if ctrl is not None:
                try:
                    value = str(ctrl.GetValue()).strip()
                except Exception:
                    value = ''
                if value:
                    return value
        try:
            return str(self.proc._raw_output_dir()).strip()
        except Exception:
            return ''

    def _schedule_base_dir(self):
        base = self._outpath_dir()
        if base:
            return base
        return ''

    def _normalize_nus_schedule(self, value, base_dir=None):
        value = str(value or '').strip()
        if not value:
            return ''
        base_dir = base_dir if base_dir is not None else self._schedule_base_dir()
        if not base_dir:
            normalized = os.path.normpath(value)
            return '' if normalized in ('.', '') else normalized
        try:
            abs_base = os.path.abspath(base_dir)
        except Exception:
            abs_base = base_dir
        try:
            abs_value = os.path.abspath(value) if os.path.isabs(value) else os.path.abspath(os.path.join(abs_base, value))
        except Exception:
            abs_value = value
        try:
            rel = os.path.relpath(abs_value, abs_base)
            if rel and not rel.startswith('..'):
                return '' if rel in ('.', '') else os.path.normpath(rel)
        except Exception:
            pass
        if os.path.isabs(value):
            return abs_value
        normalized = os.path.normpath(value)
        return '' if normalized in ('.', '') else normalized

    def _default_nus_schedule(self):
        base_dir = self._schedule_base_dir()
        root = self._root()
        try:
            saved = str(root.Parse(root.deconParFile, 'nusFil', default='')).strip()
        except Exception:
            saved = ''
        if saved:
            return self._normalize_nus_schedule(saved, base_dir=base_dir)

        if base_dir:
            # Conversion writes the expanded NUS sampling table as `schedule`.
            # Prefer that canonical project-local file; retain `nuslist` as a
            # compatibility fallback for older projects/imported datasets.
            for filename in ('schedule', 'nuslist'):
                candidate = os.path.join(base_dir, filename)
                if os.path.exists(candidate):
                    return self._normalize_nus_schedule(candidate, base_dir=base_dir)
        return ''

    def _spectrometer_type(self):
        # ProcessFrame owns spectrometer detection.  Its raw-data path does not
        # change while this window exists, so conversion code only consumes the
        # cached result rather than probing the same directory again.
        return getattr(self.proc, 'tp', None)

    def _is_bruker(self):
        return self._spectrometer_type() == 'bruk'

    def _rk_count(self):
        return max(0, int(self.proc.spectral_dim_count) - 1)

    def _rk_values_for_build(self):
        if self._is_bruker():
            return [False] * self._rk_count()
        return [ctrl.IsChecked() for ctrl in self._rk_controls()]

    def _label_value(self, idx, default=None):
        ctrl = getattr(self, f'label{idx}', None)
        if ctrl is not None:
            try:
                value = self._clean_label_text(ctrl.GetValue())
                if value:
                    return value
            except Exception:
                pass
        root = self._root()
        try:
            value = self._clean_label_text(root.ParseAllStr(root.deconParFile, f'label{idx}'))
            if value:
                return value
        except Exception:
            pass
        if self._is_bruker():
            bruker_value = self._bruker_label_default(idx)
            if bruker_value:
                return bruker_value
        if default is not None:
            return default
        return self._default_label_text(idx)

    def _reference_labels(self):
        labels = [self._label_value(i) for i in range(1, self._dimension_count() + 1)]
        if self._has_pseudo_axis():
            try:
                real_name = str(self.RealName.GetValue()).strip()
            except Exception:
                real_name = ''
            if not real_name or real_name in ('0', 'None'):
                real_name = str(self._root().Parse(self._root().deconParFile, 'RealName', default='')).strip()
            if real_name not in ('', '0', 'None'):
                labels.append(real_name)
        return labels

    def _current_reference_mode(self):
        try:
            return self.refList[self.refBox.GetSelection()]
        except Exception:
            return 'Water'

    def _reference_vpar(self, convert=False):
        self._ensure_store()
        labels = self._reference_labels()
        rk = self._rk_values_for_build()
        o1p = self._selected_reference()
        outdir = self.proc._spec_output_dir()
        if not outdir:
            outdir = './spec'
        fid_path = str(self.proc._raw_output_dir()).strip()
        try:
            inst = getattr(self.proc, 'vpar', None)
            if inst is None:
                inst = vpar()
                self.proc.vpar = inst
            inst.Setup(
                self.proc,
                outdir,
                self.proc.dim,
                labels,
                rk,
                nuslist=self.nusFil.GetValue().strip(),
                o1p=o1p,
                FidPath=fid_path,
            )
            if convert:
                inst.Convert()
            return inst
        except Exception:
            pass
            raise

    def _reference_status_text(self):
        mode = self._current_reference_mode()
        if mode == 'Manual':
            try:
                manual = float(self.xcenBox.GetValue())
            except Exception:
                return 'xCar ? ppm'
            return f'xCar {manual:g} ppm'

        try:
            inst = self._reference_vpar(convert=True)
        except Exception:
            pass
            return 'xCar ? ppm'

        xcar = getattr(inst, 'waterppm', None)
        if xcar is None:
            xcar = getattr(inst, 'waterppmTOF', None)
        if xcar is None:
            return 'xCar ? ppm'
        try:
            xcar = float(xcar)
        except Exception:
            pass

        temp = getattr(inst, 'temp', None)
        if temp is not None and str(temp).strip() not in ('', 'None'):
            try:
                return f'xCar {xcar:g} ppm ({float(temp):g} K)'
            except Exception:
                return f'xCar {xcar} ppm ({temp} K)'
        return f'xCar {xcar:g} ppm'

    def _set_status_bar_text(self, text):
        try:
            if getattr(self, 'status_bar', None) is not None:
                self.SetStatusText(text or '')
        except Exception:
            pass

    def _widget_visible(self, name):
        ctrl = getattr(self, name, None)
        if ctrl is None:
            return False
        for method in ('IsShownOnScreen', 'IsShown'):
            try:
                fn = getattr(ctrl, method, None)
                if callable(fn):
                    return bool(fn())
            except Exception:
                pass
        return True

    def _visible_text(self, name, default=''):
        ctrl = getattr(self, name, None)
        if ctrl is None or not self._widget_visible(name):
            return default
        try:
            return str(ctrl.GetValue()).strip()
        except Exception:
            return default

    def _visible_checked(self, name, default=False):
        ctrl = getattr(self, name, None)
        if ctrl is None or not self._widget_visible(name):
            return default
        try:
            return bool(ctrl.IsChecked())
        except Exception:
            return default

    def _visible_selection(self, name, default=0):
        ctrl = getattr(self, name, None)
        if ctrl is None or not self._widget_visible(name):
            return default
        try:
            return int(ctrl.GetSelection())
        except Exception:
            return default

    def _refresh_reference_status(self):
        text = self._reference_status_text()
        try:
            if hasattr(self, 'xcarStatus'):
                self.xcarStatus.SetLabel(text)
        except Exception:
            pass
        self._hover_default_status = text
        self._set_status_bar_text(text)
        try:
            if hasattr(self, 'panel'):
                self.panel.Layout()
            self.Layout()
        except Exception:
            pass

    def _on_manual_ppm_change(self, event):
        if self._current_reference_mode() == 'Manual':
            self._refresh_reference_status()
        self._hover_default_status = self._reference_status_text()
        if event is not None:
            event.Skip()

    def _parse_bool(self, root, key, default=False):
        try:
            raw = str(root.Parse(root.deconParFile, key, default='')).strip().lower()
        except Exception:
            return default
        if raw in ('y', 'yes', 'true', '1', 't'):
            return True
        if raw in ('n', 'no', 'false', '0', 'f'):
            return False
        return default

    def _load_from_file(self):
        root = self._root()

        for idx in range(1, min(self.label_count, 4) + 1):
            #print("   " ,idx)
            #print(self._is_bruker())
            ctrl = getattr(self, f'label{idx}', None)
            if ctrl is None:
                continue
            try:
                raw = root.ParseAllStr(root.deconParFile, f'label{idx}')
            except Exception:
                pass
                raw = ''
            raw = self._clean_label_text(raw)
            #print(raw,not raw)
            #print(f'label{idx} from system file {root.deconParFile}: {raw!r}')
            #print('Conversion label load: label%d from system file %s -> %r' % (idx, root.deconParFile, raw))
            if not raw and self._is_bruker():
                raw = self._bruker_label_default(idx)
            if not raw:
                raw = self._default_label_text(idx)
            self._safe_set_text(ctrl, raw)

        for idx in range(1, 4):
            ctrl = getattr(self, f'cb_rk{idx}', None)
            if ctrl is None:
                continue
            ctrl.SetValue(self._parse_bool(root, f'rk{idx}', False))
        if hasattr(self, 'RealName'):
            try:
                real_name = str(root.Parse(root.deconParFile, 'RealName', default='')).strip()
            except Exception:
                real_name = ''
            if real_name in ('', '0', 'None'):
                real_name = 'pseudo'
            self._safe_set_text(self.RealName, real_name)
        self._safe_set_text(self.xcenBox, root.Parse(root.deconParFile, 'xcen', default=''))
        self._safe_set_text(self.nusFil, self._default_nus_schedule())

        try:
            self.refBox.SetSelection(int(root.Parse(root.deconParFile, 'refBox', default=0)))
        except Exception:
            pass

        self._on_reference_change(None)
        self._update_widget_visibility()
        self._refresh_reference_status()

    def _label_controls(self):
        controls = []
        for name in ('label1', 'label2', 'label3', 'label4'):
            ctrl = getattr(self, name, None)
            if ctrl is not None:
                controls.append(ctrl)
        return controls

    def _rk_controls(self):
        controls = []
        for name in ('cb_rk1', 'cb_rk2', 'cb_rk3'):
            ctrl = getattr(self, name, None)
            if ctrl is not None:
                controls.append(ctrl)
        return controls

    def _load_from_parent(self):
        self._load_from_file()

    def _copy_to_parent(self):
        """Legacy no-op: conversion now owns its own widgets/state."""
        return

    def _selected_reference(self):
        mode = self.refList[self.refBox.GetSelection()]
        if mode == 'Manual':
            try:
                return float(self.xcenBox.GetValue())
            except Exception:
                raise ValueError('Manual referencing requires a numeric ppm value.')
        return mode

    def _on_reference_change(self, event):
        manual = self.refList[self.refBox.GetSelection()] == 'Manual'
        self.xcenBox.Enable(manual)
        self._refresh_reference_status()
        if event is not None:
            event.Skip()

    def _format_summary(self, inst):
        lines = []
        if hasattr(inst, 'tp'):
            lines.append(f'Source type: {inst.tp}')
        if hasattr(inst, 'seqfil'):
            lines.append(f'Pulse program: {inst.seqfil}')
        if hasattr(inst, 'np'):
            lines.append(f'Direct points: {getattr(inst, "np", "?")}')
        if hasattr(inst, 'np2'):
            lines.append(f'Direct half points: {getattr(inst, "np2", "?")}')
        for key in ('ni', 'ni2', 'ni3', 'nz'):
            if hasattr(inst, key):
                lines.append(f'{key}: {getattr(inst, key)}')
        for key, label in (
            ('aqseq', 'aqseq'),
            ('bytor', 'byte order'),
            ('DTYPA', 'DTYPA'),
            ('DECIM', 'DECIM'),
            ('DSPFVS', 'DSPFVS'),
            ('GRPDLY', 'GRPDLY'),
        ):
            if hasattr(inst, key):
                lines.append(f'{label}: {getattr(inst, key)}')
        if hasattr(inst, 'brukerFmt'):
            lines.append(f'Bruker format: {getattr(inst, "brukerFmt")}')
        if hasattr(inst, 'waterppm'):
            lines.append(f'Reference ppm: {getattr(inst, "waterppm")}')
        if hasattr(inst, 'f1ppm'):
            lines.append(f'Indirect 1 carrier: {getattr(inst, "f1ppm")}')
        if hasattr(inst, 'f2ppm'):
            lines.append(f'Indirect 2 carrier: {getattr(inst, "f2ppm")}')
        if hasattr(inst, 'f3ppm'):
            lines.append(f'Indirect 3 carrier: {getattr(inst, "f3ppm")}')
        if hasattr(self, 'nusFil') and self.nusFil.GetValue().strip():
            lines.append(f'NUS schedule: {self.nusFil.GetValue().strip()}')
        return '\n'.join(lines) if lines else 'No conversion metadata available.'

    def _script_target(self):
        try:
            outdir = self.proc._spec_output_dir()
        except Exception:
            outdir = ''
        if not outdir:
            outdir = './spec'
        return os.path.join(outdir, 'fid.test.com')

    def _ensure_store(self):
        raw_dir = str(self.proc._raw_output_dir()).strip()
        if not os.path.exists(raw_dir):
            raise ValueError('Datastore does not exist. Check the project paths first.')

    def _source_metadata_status(self):
        raw_dir = str(self.proc._raw_output_dir()).strip()
        procpar = os.path.join(raw_dir, 'procpar')
        acqus = os.path.join(raw_dir, 'acqus')
        if os.path.isfile(procpar):
            return 'var', procpar
        if os.path.isfile(acqus):
            return 'bruk', acqus
        return None, None

    def _missing_label_indices(self):
        missing = []
        for idx in range(1, self._dimension_count() + 1):
            ctrl = getattr(self, f'label{idx}', None)
            value = self._clean_label_text(ctrl.GetValue()) if ctrl is not None else ''
            if not value:
                missing.append(idx)
        if self._has_pseudo_axis() and hasattr(self, 'RealName'):
            if not self._clean_label_text(self.RealName.GetValue()):
                missing.append('pseudo')
        return missing

    def _preflight_conversion(self):
        """Validate user-fixable inputs before vpar/script generation."""
        self.collect_updates()
        self._ensure_store()
        vendor, _metadata_path = self._source_metadata_status()
        if vendor is None:
            raise ValueError(
                'Source metadata is missing.\n\n'
                'DECON found the raw-data folder, but could not find the acquisition '
                'metadata required for automatic conversion.\n\n'
                'Expected:\n  Varian: procpar\n  Bruker: acqus\n\n'
                'Automatic conversion cannot safely determine the acquisition parameters '
                'without this file. If you know the required NMRPipe parameters, you can '
                'create or edit the conversion script manually.'
            )

        missing = self._missing_label_indices()
        if missing:
            dims = ', '.join(str(x) for x in missing)
            varian_note = (
                '\n\nThis is common with Varian datasets because nucleus labels cannot always '
                'be inferred reliably from procpar.'
            ) if vendor == 'var' else ''
            raise ValueError(
                'Axis labels are required.\n\n'
                f'Label(s) {dims} could not be determined automatically.{varian_note}\n\n'
                'Please enter the missing axis label values in the Advanced tab before '
                'using automatic conversion.'
            )

        mode = self._current_reference_mode()
        if mode == 'Manual':
            self._selected_reference()
            return
        try:
            inst = self._reference_vpar(convert=True)
            xcar = getattr(inst, 'waterppm', None)
            if xcar is None:
                xcar = getattr(inst, 'waterppmTOF', None)
            if xcar is None:
                raise ValueError('No reference position was calculated.')
        except Exception as exc:
            if mode == 'Water':
                raise ValueError(
                    'Water referencing could not be determined.\n\n'
                    'The acquisition temperature or other metadata needed to calculate '
                    'the water reference could not be read reliably.\n\n'
                    'Choose Manual referencing and enter the reference position in ppm, '
                    'or correct the source acquisition metadata.'
                ) from exc
            raise ValueError(
                'Automatic referencing could not be determined.\n\n'
                'Choose Manual referencing and enter the reference position in ppm, '
                'or correct the source acquisition metadata.'
            ) from exc

    def _build_vpar(self):
        self._preflight_conversion()
        # Make the live state agree with the widgets before conversion script
        # inference. vpar below continues to consume the widgets directly.
        self.collect_updates()
        self._ensure_store()
        labels = [getattr(self, f'label{i}').GetValue().strip()
                  for i in range(1, self._dimension_count() + 1)]
        if self._has_pseudo_axis():
            real_name = self.RealName.GetValue().strip()
            if real_name:
                labels.append(real_name)

        rk = self._rk_values_for_build()

        o1p = self._selected_reference()
        outdir = self.proc._spec_output_dir()
        if not outdir:
            outdir = './spec'
        fid_path = str(self.proc._raw_output_dir()).strip()
        inst = getattr(self.proc, 'vpar', None)
        if inst is None:
            inst = vpar()
            self.proc.vpar = inst
        inst.Setup(
            self.proc,
            outdir,
            self._legacy_backend_dim(),
            labels,
            rk,
            nuslist=self.nusFil.GetValue().strip(),
            o1p=o1p,
            FidPath=fid_path,
        )
        if not getattr(inst, 'initialized', False):
            reason = getattr(inst, 'setup_error', None)
            if reason == 'missing_source_metadata':
                raise ValueError('Source metadata is missing (expected procpar or acqus).')
            if reason == 'missing_raw_data':
                raise ValueError('No raw fid/ser data file could be found in the source folder.')
            raise ValueError('Conversion setup could not be initialized from the source data.')
        return inst

    def _generate_guess_script(self):
        inst = self._build_vpar()
        script_path = inst.BuildConversionScript()
        if script_path == -1 or getattr(inst, 'abort', 0) == 1:
            raise RuntimeError('Conversion inference failed for the current settings.')
        self.inst = inst
        self.script_path = script_path
        return script_path

    def _ensure_script_ready(self):
        script_path = self._script_target()
        if not os.path.exists(script_path):
            return self._generate_guess_script()
        self.script_path = script_path
        return script_path

    def _save_script_file(self, script_path=None):
        script_path = script_path or self._ensure_script_ready()
        if not os.path.exists(script_path):
            raise ValueError('No script exists to save.')
        text = Path(script_path).read_text()
        Path(script_path).write_text(text)
        self.script_path = script_path
        return script_path

    def _run_script_path(self, script_path, on_finish=None):
        def _finish(*_args, **_kwargs):
            try:
                self.proc.UpdateLampLights()
            except Exception:
                pass
            # Continue slice conversion/preview processing in the same output
            # window. ProcessSlice invokes the final UI callback when the whole
            # workflow, rather than just fid.test.com, has completed.
            if self.inst is not None and hasattr(self.inst, 'ProcessSlice'):
                try:
                    self.inst.ProcessSlice(output_frame=output_frame, on_finish=_workflow_done)
                    return
                except Exception as exc:
                    # Slice preparation is part of conversion, not an optional
                    # cosmetic step.  Previously this exception was swallowed,
                    # which made a failed Bruker/nmrglue extraction look like a
                    # successful conversion while slice.fid was simply absent.
                    try:
                        output_frame.append_text('\nPreview slice preparation failed: %s\n' % exc)
                        output_frame.finish_workflow(False)
                        output_frame.set_status('Failed')
                    except Exception:
                        pass
                    return
            _workflow_done()

        def _workflow_done(*_args, **_kwargs):
            try:
                output_frame.start_step('Refresh display')
            except Exception:
                pass
            try:
                self.proc.UpdateLampLights()
            except Exception:
                pass
            try:
                # A freshly converted multidimensional preview should open on
                # the most informative trace, rather than always FID 1.
                self.proc.auto_select_strongest_fid()
            except Exception:
                pass
            try:
                # New data should be immediately visible at its natural
                # vertical scale (the same behaviour as the toolbar Draw).
                self.proc.draw_figure(reset_y=True)
            except Exception:
                pass
            try:
                output_frame.append_text('\nConversion workflow complete.\n')
                output_frame.finish_workflow(True)
                output_frame.set_status('Complete')
            except Exception:
                pass
            if on_finish is not None:
                try:
                    on_finish()
                except Exception:
                    pass

        # Conversion has only a few meaningful milestones.  Keep these coarse
        # so a routine user can see progress without having to understand the
        # NMRPipe transcript.
        needs_extract = False
        try:
            needs_extract = self.inst._slice_dim_count() > 1 and not (self.inst.dim == '3p' and self.inst.tp == 'bruk')
        except Exception:
            pass
        conversion_steps = ['Convert raw data to NMRPipe']
        if needs_extract:
            conversion_steps.append('Prepare preview slice')
        conversion_steps.extend(['Process preview spectrum', 'Refresh display'])
        output_frame = run_command_with_output(
            ['csh', script_path],
            parent=self,
            title='Conversion Output',
            on_finish=_finish,
            final=False,
            label='Convert raw data to NMRPipe',
            workflow_steps=conversion_steps,
            workflow_step=0,
        )
        return script_path

    def _browse_nus(self, event):
        cwd = self._schedule_base_dir() or os.getcwd()
        try:
            cwd = os.path.abspath(cwd)
        except Exception:
            pass
        dlg = wx.FileDialog(
            self,
            message='Choose NUS schedule',
            defaultDir=cwd,
            defaultFile='',
            wildcard='All files (*.*)|*.*',
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
        )
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.nusFil.SetValue(self._normalize_nus_schedule(dlg.GetPath(), base_dir=self._schedule_base_dir()))
        finally:
            dlg.Destroy()

    # ------------------------------------------------------------------
    # Actions
    def on_show_script(self, event):
        try:
            script_path = self._ensure_script_ready()
            if getattr(self, 'script_frame', None) is not None:
                try:
                    self.script_frame.Raise()
                    self.script_frame.SetFocus()
                    return
                except Exception:
                    self.script_frame = None
            self.script_frame = ConversionScriptFrame(self, script_path=script_path)
            self.script_frame.Show(True)
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def collect_updates(self, update_state=True):
        """Return current conversion widget values without writing to disk."""
        root = self._root()
        updates = {}
        for idx in range(1, min(self.label_count, 4) + 1):
            name = f'label{idx}'
            ctrl = getattr(self, name, None)
            if ctrl is not None:
                updates[name] = self._clean_label_text(ctrl.GetValue())
        for idx in range(1, 4):
            name = f'cb_rk{idx}'
            if hasattr(self, name) and self._widget_visible(name):
                updates[f'rk{idx}'] = self.Parent.IntToBool(self._visible_checked(name, self._parse_bool(root, f'rk{idx}', False)))
        if hasattr(self, 'RealName') and self._widget_visible('RealName'):
            updates['RealName'] = self._clean_label_text(self.RealName.GetValue())
        if hasattr(self, 'xcenBox') and self._widget_visible('xcenBox'):
            updates['xcen'] = self._visible_text('xcenBox', '')
        if hasattr(self, 'nusFil') and self._widget_visible('nusFil'):
            updates['nusFil'] = self._visible_text('nusFil', '')
        if hasattr(self, 'refBox') and self._widget_visible('refBox'):
            updates['refBox'] = str(self._visible_selection('refBox', 0))
        state = getattr(self.proc, 'state', None)
        if state is not None and update_state:
            state.update_gui_settings(updates)
            setter = getattr(self.proc, 'set_dimension_labels', None)
            if callable(setter):
                setter([updates.get(f'label{i}', '') for i in range(1, min(self.label_count, 4) + 1)], refresh=True)
        return updates

    def on_run(self, event):
        try:
            self.proc.save_current_gui_state(reason='conversion-run')
            script_path = self._generate_guess_script()
            self._save_script_file(script_path)
            self._run_script_path(script_path)
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def _on_close(self, event):
        try:
            self.collect_updates()
        except Exception:
            pass
        try:
            if getattr(self, 'script_frame', None) is not None:
                try:
                    self.script_frame.Destroy()
                except Exception:
                    pass
                self.script_frame = None
            if hasattr(self.proc, 'conv_frame') and self.proc.conv_frame is self:
                self.proc.conv_frame = None
        except Exception:
            pass
        self.Destroy()


class ConversionScriptFrame(wx.Frame):
    """Editable fid.test.com window that does not touch GUI save/load state."""

    def __init__(self, conv, script_path=None):
        super().__init__(conv, title='fid.test.com', size=(980, 480))
        self.conv = conv
        self.script_path = script_path or self.conv._script_target()
        self.inst = None
        self.panel = None
        self.status_box = None
        self.extra_box = None
        self.status_text = None
        self.extra_text = None
        self.path_label = None
        self.scriptBox = None
        self._build_ui()
        self._load_script()
        self._refresh_metadata()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _status_lines(self):
        return [
            'Byte order: inferred from acqus',
            'Data type: inferred from acqus',
            'Digital filter: inferred from acqus',
            'Bruker format: inferred from acqus',
        ]

    def _extra_lines(self, inst=None):
        inst = inst or getattr(self.conv, 'inst', None)
        if inst is None:
            try:
                inst = self.conv._build_vpar()
            except Exception:
                inst = None
        if inst is None:
            return ['No extra information available.']

        lines = []

        def add(label, value):
            if value is None:
                return
            text = str(value).strip()
            if text and text not in ('0', 'None'):
                lines.append(f'{label}: {text}')

        add('Pulse sequence', getattr(inst, 'seqfil', None))
        add('Spectrometer type', getattr(inst, 'tp', None))
        add('Acquisition order', getattr(inst, 'aqseq', None))
        add('Data path', getattr(inst, 'FidPath', None))
        add('Output path', getattr(inst, 'outdir', None))
        add('Dimension', getattr(inst, 'dim', None))
        add('Labels', ', '.join([str(x) for x in getattr(inst, 'labb', []) if str(x).strip() not in ('', '0', 'None')]))
        add('Rance-Kay', ', '.join(['yes' if x else 'no' for x in getattr(inst, 'rk', [])]))
        add('Reference mode', getattr(inst, 'o1p', None))
        add('Water ppm', getattr(inst, 'waterppm', None))
        add('Water ppm (TOF)', getattr(inst, 'waterppmTOF', None))
        if hasattr(inst, 'temp'):
            add('Temperature', f'{getattr(inst, "temp")} K')
        add('NUS schedule', getattr(inst, 'nuslist', None))
        add('Bruker format', getattr(inst, 'brukerFmt', None))
        for key in ('bytor', 'DTYPA', 'DECIM', 'DSPFVS', 'GRPDLY', 'BYTORDA'):
            add(key, getattr(inst, key, None))
        for key in ('np', 'np2', 'ni', 'ni2', 'ni3', 'nz', 'sw', 'sw1', 'sw2', 'sw3'):
            add(key.upper(), getattr(inst, key, None))
        for key in ('SFO1', 'SFO2', 'SFO3', 'SFO4', 'BF1', 'BF2', 'BF3', 'BF4', 'O1', 'O2', 'O3', 'O4'):
            add(key, getattr(inst, key, None))
        for key in ('NUC1', 'NUC2', 'NUC3', 'NUC4'):
            add(key, getattr(inst, key, None))
        for key in ('f1ppm', 'f2ppm', 'f3ppm'):
            add(key, getattr(inst, key, None))
        return lines or ['No extra information available.']

    def _refresh_metadata(self):
        if self.status_text is not None:
            self.status_text.SetValue('\n'.join(self._status_lines()))
        if self.extra_text is not None:
            self.extra_text.SetValue('\n'.join(self._extra_lines()))

    def _build_ui(self):
        self.panel = wx.Panel(self)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        root = wx.BoxSizer(wx.VERTICAL)

        top_row = wx.BoxSizer(wx.HORIZONTAL)
        status_static = wx.StaticBox(self.panel, -1, 'Status')
        extra_static = wx.StaticBox(self.panel, -1, 'Extra information')
        self.status_box = wx.StaticBoxSizer(status_static, wx.VERTICAL)
        self.extra_box = wx.StaticBoxSizer(extra_static, wx.VERTICAL)
        self.status_text = wx.TextCtrl(
            self.panel,
            style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP,
            size=(280, 72),
        )
        self.extra_text = wx.TextCtrl(
            self.panel,
            style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP,
            size=(520, 72),
        )
        self.status_box.Add(self.status_text, 1, wx.EXPAND | wx.ALL, 4)
        self.extra_box.Add(self.extra_text, 1, wx.EXPAND | wx.ALL, 4)
        top_row.Add(self.status_box, 0, wx.RIGHT | wx.EXPAND, 12)
        top_row.Add(self.extra_box, 1, wx.EXPAND)
        root.Add(top_row, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)

        self.path_label = wx.StaticText(self.panel, label=os.path.basename(self.script_path))
        root.Add(self.path_label, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)

        self.scriptBox = wx.TextCtrl(self.panel, style=wx.TE_MULTILINE | wx.TE_RICH | wx.TE_DONTWRAP)
        self.scriptBox.SetMinSize((920, 240))

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.guessBtn = wx.Button(self.panel, label='Guess')
        self.saveBtn = wx.Button(self.panel, label='Save')
        self.runBtn = wx.Button(self.panel, label='Run')
        self.closeBtn = wx.Button(self.panel, label='Close')
        btn_row.Add(self.guessBtn, 0, wx.RIGHT, 6)
        btn_row.Add(self.saveBtn, 0, wx.RIGHT, 6)
        btn_row.Add(self.runBtn, 0, wx.RIGHT, 6)
        btn_row.Add(self.closeBtn, 0)

        self.guessBtn.Bind(wx.EVT_BUTTON, self.on_guess)
        self.saveBtn.Bind(wx.EVT_BUTTON, self.on_save)
        self.runBtn.Bind(wx.EVT_BUTTON, self.on_run)
        self.closeBtn.Bind(wx.EVT_BUTTON, self._on_close)

        root.Add(self.scriptBox, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 10)
        root.Add(btn_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.ALIGN_RIGHT, 10)

        self.panel.SetSizer(root)
        self.panel.Layout()
        self.panel.Fit()
        best = self.panel.GetBestSize()
        self.SetMinSize(best)
        self.SetClientSize(best)
        self.Layout()

    def _load_script(self):
        if self.script_path and os.path.exists(self.script_path):
            try:
                self.scriptBox.SetValue(Path(self.script_path).read_text())
                return
            except Exception:
                pass
        try:
            script_path = self.conv._generate_guess_script()
            self.inst = getattr(self.conv, 'inst', None)
            self.script_path = script_path
            self.path_label.SetLabel(os.path.basename(script_path))
            self.scriptBox.SetValue(Path(script_path).read_text())
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def on_guess(self, event):
        try:
            self.conv.proc.save_current_gui_state(reason='conversion-script-guess')
            script_path = self.conv._generate_guess_script()
            self.inst = getattr(self.conv, 'inst', None)
            self.script_path = script_path
            self.path_label.SetLabel(os.path.basename(script_path))
            self.scriptBox.SetValue(Path(script_path).read_text())
            self._refresh_metadata()
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def on_save(self, event):
        try:
            if not self.script_path:
                self.script_path = self.conv._script_target()
            Path(self.script_path).write_text(self.scriptBox.GetValue())
            self.conv.script_path = self.script_path
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def on_run(self, event):
        try:
            self.on_save(event)
            self.conv.proc.save_current_gui_state(reason='conversion-script-run')
            self.conv.script_path = self.script_path

            def _finish(*_args, **_kwargs):
                try:
                    self.conv.proc.UpdateLampLights()
                except Exception:
                    pass
                try:
                    self.conv.proc.draw_figure()
                except Exception:
                    pass

            self.conv._run_script_path(self.script_path, on_finish=_finish)
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def _on_close(self, event):
        try:
            if hasattr(self.conv, 'script_frame') and self.conv.script_frame is self:
                self.conv.script_frame = None
        except Exception:
            pass
        self.Destroy()

