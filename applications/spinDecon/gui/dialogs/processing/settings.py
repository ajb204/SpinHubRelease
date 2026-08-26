import re
import os
from pathlib import Path

import wx

from spinDecon.gui.dialogs.errors import errorMessage
from spinDecon.gui.dialogs.shell_output import ShellOutputFrame


def _shell_ascii(text):
    """Force executable processing scripts to ASCII shell punctuation.

    This is intentionally duplicated at the GUI save boundary: the script editor
    can load/cache text independently of the NMRPipe writer, and Run saves the
    editor contents immediately before execution.
    """
    return (str(text).replace("\u201c", '"').replace("\u201d", '"')
                     .replace("\u2018", "'").replace("\u2019", "'")
                     .replace("\u2014", "--").replace("\u2013", "-")
                     .replace("\u2026", "..."))


class ProcessingFrame(wx.Frame):
    """Standalone processing window for baseline, LP, and save/run controls."""

    def __init__(self, parent):
        super().__init__(parent, title='Processing', size=(1160, 880))
        self.proc = parent
        self.spectral_dim_count = int(parent.spectral_dim_count)
        self.physical_dim_count = int(parent.physical_dim_count)
        self.has_pseudo_axis = bool(parent.has_pseudo_axis)
        self.topology = parent.topology
        self.dim = self.spectral_dim_count  # compatibility alias: spectral only
        self.spectral_dim_count_choices = getattr(parent, 'dim_choices', ['H1', 'N15', 'C13', 'F19', 'P31'])
        self.ftlisty = ['Auto', 'Neg', 'Alt', 'AltNeg', 'Real']
        self.ftdic = {name: idx for idx, name in enumerate(self.ftlisty)}
        self.ftdic['y'] = 1
        self.ftdic['n'] = 0
        self.script_mode = 'process'
        self.script_path = None
        self.script_frame = None
        self._autoapodise_progress_frame = None
        self._script_cache = {'n': '', 'y': ''}
        self._last_script_lp = 'n'
        self.status_bar = None
        self._hover_default_status = 'Ready'

        try:
            self.status_bar = self.CreateStatusBar(1)
            self.status_bar.SetStatusText('Ready')
        except Exception:
            self.status_bar = None

        self._build_ui()
        self._load_from_file()
        # Hydrate shared state from disk-loaded controls without replacing any
        # newer values already edited in another window.
        state = getattr(self.proc, 'state', None)
        if state is not None:
            state.seed_gui_settings(self._processing_updates(update_state=False))
        self.apply_live_settings()
        self._bind_live_state_controls()
        self.proc._bind_processing_controls(self)
        self.p0.Bind(wx.EVT_TEXT, self._on_direct_phase_change)
        self.p1.Bind(wx.EVT_TEXT, self._on_direct_phase_change)
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _bind_live_state_controls(self):
        """Keep ProjectState.gui_settings current as processing widgets change.

        The parameter file remains a commit boundary; this handler only snapshots
        widgets into shared memory so other windows and script generators see the
        newest values before an explicit/automatic commit.
        """
        def changed(event):
            try:
                self._processing_updates()
            except Exception:
                pass
            try:
                event.Skip()
            except Exception:
                pass

        seen = set()
        for name, ctrl in vars(self).items():
            if id(ctrl) in seen:
                continue
            if not (name.startswith(('p0', 'p1', 'win2Val', 'win3Val', 'firstPoint', 'windowBox', 'cb_', 'lp', 'bl'))
                    or name in ('maxIterBox', 'ncpusBox', 'scriptModeBox')):
                continue
            seen.add(id(ctrl))
            for evt in (wx.EVT_TEXT, wx.EVT_CHOICE, wx.EVT_CHECKBOX):
                try:
                    ctrl.Bind(evt, changed)
                except Exception:
                    pass

    def apply_live_settings(self):
        """Overlay unsaved shared settings after the normal file load."""
        state = getattr(self.proc, 'state', None)
        live = getattr(state, 'gui_settings', {}) if state is not None else {}
        if 'ProcTarg' in live:
            self._apply_script_target(live['ProcTarg'])
        for key, value in live.items():
            ctrl = getattr(self, key, None)
            if ctrl is not None:
                try:
                    ctrl.SetValue(str(value))
                    continue
                except Exception:
                    pass
            if key.startswith('window'):
                ctrl = getattr(self, 'windowBox' + key[len('window'):], None)
                try: ctrl.SetSelection(int(value))
                except Exception: pass
            elif key.startswith('flip'):
                ctrl = getattr(self, 'cb_ft' + key[len('flip'):], None)
                try: ctrl.SetValue(str(value))
                except Exception: pass
            elif key in ('lin', 'poly', 'sol'):
                name = {'lin':'cb_baseLin','poly':'cb_basepol','sol':'cb_baseSol'}[key]
                ctrl = getattr(self, name, None)
                try: ctrl.SetValue(str(value).lower() in ('y','yes','1','true'))
                except Exception: pass
            elif key.startswith(('f', 'lp', 'bl')):
                # Dimension-specific boolean keys.
                if key.startswith('f') and key.endswith('180'): name = 'cb_' + key
                elif key.startswith('lp'): name = 'cb_' + key
                elif key.startswith('bl'): name = 'cb_basepol' + key[2:]
                else: name = ''
                ctrl = getattr(self, name, None)
                try: ctrl.SetValue(str(value).lower() in ('y','yes','1','true'))
                except Exception: pass

    # ------------------------------------------------------------------

    def _set_status_text(self, text):
        try:
            if getattr(self, 'status_bar', None) is not None:
                self.status_bar.SetStatusText(str(text or ''))
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
        panel.SetBackgroundColour(wx.Colour(255, 255, 255))

        root = wx.BoxSizer(wx.VERTICAL)

        # Settings.
        self.dataLbl2 = wx.StaticBox(panel, -1, 'Settings')
        self.dataSizer2 = wx.StaticBoxSizer(self.dataLbl2, wx.VERTICAL)

        direct_row = wx.BoxSizer(wx.HORIZONTAL)
        self.baseLab = wx.StaticText(self.dataLbl2, label='Direct dimension:')
        self.cb_baseLin = wx.CheckBox(self.dataLbl2, -1, 'Linear baseline', style=wx.ALIGN_RIGHT)
        self.cb_baseSol = wx.CheckBox(self.dataLbl2, -1, 'Digital solvent suppress', style=wx.ALIGN_RIGHT)
        direct_row.Add(self.baseLab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        direct_row.Add(self.cb_baseLin, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 16)
        direct_row.Add(self.cb_baseSol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.dataSizer2.Add(direct_row, 0, wx.ALL | wx.EXPAND, 8)

        self.sizer2 = wx.GridBagSizer(10, 10)
        self.f1180Lab = wx.StaticText(self.dataLbl2, label='f1180')
        self.lpLab = wx.StaticText(self.dataLbl2, label='LP')
        self.polyLab = wx.StaticText(self.dataLbl2, label='Poly')
        self.p0Lab = wx.StaticText(self.dataLbl2, label='P0')
        self.p1Lab = wx.StaticText(self.dataLbl2, label='P1')
        self.flipLab = wx.StaticText(self.dataLbl2, label='Flip')
        self.windowLab = wx.StaticText(self.dataLbl2, label='Window')
        self.windowOp1 = wx.StaticText(self.dataLbl2, label='Op1')
        self.windowOp2 = wx.StaticText(self.dataLbl2, label='Op2')
        self.firstPointFactor = wx.StaticText(self.dataLbl2, label='First Pt')
        self.lab0 = wx.StaticText(self.dataLbl2, label='dim 1', size=(50, 22))

        # Dim 1 processing controls mirrored from the original process window.
        self.p0 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
        self.p1 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
        self.cb_basepol = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
        self.cb_ft0 = wx.ComboBox(self.dataLbl2, -1, choices=self.ftlisty, style=wx.CB_READONLY, size=(-1, 22))
        self.windowBox0 = wx.ComboBox(self.dataLbl2, -1, choices=['GM', 'SP', 'EM'], style=wx.CB_READONLY, size=(-1, 22))
        self.windowBox0.Bind(wx.EVT_COMBOBOX, self.set_win)
        self.firstPoint0 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
        self.win2Val0 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
        self.win3Val0 = wx.TextCtrl(self.dataLbl2, size=(40, 22))


        def add_dim_headers(r):
            self.sizer2.Add(self.f1180Lab, (r, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.lpLab, (r, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.polyLab, (r, 3), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.p0Lab, (r, 4), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.p1Lab, (r, 5), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.flipLab, (r, 6), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.windowLab, (r, 7), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.windowOp1, (r, 8), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.windowOp2, (r, 9), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.firstPointFactor, (r, 10), flag=wx.ALIGN_CENTER_VERTICAL)

        add_dim_headers(0)
        self.sizer2.Add(self.lab0, (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.Add(self.cb_basepol, (1, 3), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.Add(self.p0, (1, 4))
        self.sizer2.Add(self.p1, (1, 5))
        self.sizer2.Add(self.cb_ft0, (1, 6))
        self.sizer2.Add(self.windowBox0, (1, 7))
        self.sizer2.Add(self.win2Val0, (1, 8))
        self.sizer2.Add(self.win3Val0, (1, 9))
        self.sizer2.Add(self.firstPoint0, (1, 10))


        current_row = 2
        if isinstance(self.spectral_dim_count, int):
            if self.spectral_dim_count >= 2:
                self.lab1 = wx.StaticText(self.dataLbl2, label='dim 2', size=(50, 22))
                self.cb_f1180 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_lp1 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_basepol1 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.p0_1 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.p1_1 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.cb_ft1 = wx.ComboBox(self.dataLbl2, -1, choices=self.ftlisty, style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox1 = wx.ComboBox(self.dataLbl2, -1, choices=['GM', 'SP', 'EM'], style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox1.Bind(wx.EVT_COMBOBOX, self.set_win)
                self.firstPoint1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win2Val1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win3Val1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))

                self.sizer2.Add(self.lab1, (current_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
                self.sizer2.Add(self.cb_f1180, (current_row, 1))
                self.sizer2.Add(self.cb_lp1, (current_row, 2))
                self.sizer2.Add(self.cb_basepol1, (current_row, 3))
                self.sizer2.Add(self.p0_1, (current_row, 4))
                self.sizer2.Add(self.p1_1, (current_row, 5))
                self.sizer2.Add(self.cb_ft1, (current_row, 6))
                self.sizer2.Add(self.windowBox1, (current_row, 7))
                self.sizer2.Add(self.win2Val1, (current_row, 8))
                self.sizer2.Add(self.win3Val1, (current_row, 9))
                self.sizer2.Add(self.firstPoint1, (current_row, 10))
                current_row += 1

            if self.spectral_dim_count >= 3:
                self.lab2 = wx.StaticText(self.dataLbl2, label='dim 3', size=(50, 22))
                self.cb_f2180 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_lp2 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_basepol2 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.p0_2 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.p1_2 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.cb_ft2 = wx.ComboBox(self.dataLbl2, -1, choices=self.ftlisty, style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox2 = wx.ComboBox(self.dataLbl2, -1, choices=['GM', 'SP', 'EM'], style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox2.Bind(wx.EVT_COMBOBOX, self.set_win)
                self.firstPoint2 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win2Val2 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win3Val2 = wx.TextCtrl(self.dataLbl2, size=(40, 22))

                self.sizer2.Add(self.lab2, (current_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
                self.sizer2.Add(self.cb_f2180, (current_row, 1))
                self.sizer2.Add(self.cb_lp2, (current_row, 2))
                self.sizer2.Add(self.cb_basepol2, (current_row, 3))
                self.sizer2.Add(self.p0_2, (current_row, 4))
                self.sizer2.Add(self.p1_2, (current_row, 5))
                self.sizer2.Add(self.cb_ft2, (current_row, 6))
                self.sizer2.Add(self.windowBox2, (current_row, 7))
                self.sizer2.Add(self.win2Val2, (current_row, 8))
                self.sizer2.Add(self.win3Val2, (current_row, 9))
                self.sizer2.Add(self.firstPoint2, (current_row, 10))
                current_row += 1

            if self.spectral_dim_count == 4:
                self.lab3 = wx.StaticText(self.dataLbl2, label='dim 4', size=(50, 22))
                self.cb_f3180 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_lp3 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.cb_basepol3 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
                self.p0_3 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.p1_3 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
                self.cb_ft3 = wx.ComboBox(self.dataLbl2, -1, choices=self.ftlisty, style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox3 = wx.ComboBox(self.dataLbl2, -1, choices=['GM', 'SP', 'EM'], style=wx.CB_READONLY, size=(-1, 22))
                self.windowBox3.Bind(wx.EVT_COMBOBOX, self.set_win)
                self.firstPoint3 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win2Val3 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
                self.win3Val3 = wx.TextCtrl(self.dataLbl2, size=(40, 22))

                self.sizer2.Add(self.lab3, (current_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
                self.sizer2.Add(self.cb_f3180, (current_row, 1))
                self.sizer2.Add(self.cb_lp3, (current_row, 2))
                self.sizer2.Add(self.cb_basepol3, (current_row, 3))
                self.sizer2.Add(self.p0_3, (current_row, 4))
                self.sizer2.Add(self.p1_3, (current_row, 5))
                self.sizer2.Add(self.cb_ft3, (current_row, 6))
                self.sizer2.Add(self.windowBox3, (current_row, 7))
                self.sizer2.Add(self.win2Val3, (current_row, 8))
                self.sizer2.Add(self.win3Val3, (current_row, 9))
                self.sizer2.Add(self.firstPoint3, (current_row, 10))
                current_row += 1

        elif self.spectral_dim_count == 2 and self.has_pseudo_axis:
            self.lab1 = wx.StaticText(self.dataLbl2, label='dim 2', size=(50, 22))
            self.cb_f1180 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
            self.cb_lp1 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
            self.cb_basepol1 = wx.CheckBox(self.dataLbl2, -1, '', style=wx.ALIGN_RIGHT)
            self.p0_1 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
            self.p1_1 = wx.TextCtrl(self.dataLbl2, size=(50, 22))
            self.cb_ft1 = wx.ComboBox(self.dataLbl2, -1, choices=self.ftlisty, style=wx.CB_READONLY, size=(-1, 22))
            self.windowBox1 = wx.ComboBox(self.dataLbl2, -1, choices=['GM', 'SP', 'EM'], style=wx.CB_READONLY, size=(-1, 22))
            self.windowBox1.Bind(wx.EVT_COMBOBOX, self.set_win)
            self.firstPoint1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
            self.win2Val1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))
            self.win3Val1 = wx.TextCtrl(self.dataLbl2, size=(40, 22))

            self.sizer2.Add(self.lab1, (current_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(self.cb_f1180, (current_row, 1))
            self.sizer2.Add(self.cb_lp1, (current_row, 2))
            self.sizer2.Add(self.cb_basepol1, (current_row, 3))
            self.sizer2.Add(self.p0_1, (current_row, 4))
            self.sizer2.Add(self.p1_1, (current_row, 5))
            self.sizer2.Add(self.cb_ft1, (current_row, 6))
            self.sizer2.Add(self.windowBox1, (current_row, 7))
            self.sizer2.Add(self.win2Val1, (current_row, 8))
            self.sizer2.Add(self.win3Val1, (current_row, 9))
            self.sizer2.Add(self.firstPoint1, (current_row, 10))

        self.dataSizer2.Add(self.sizer2, 1, wx.ALL | wx.EXPAND, 8)

        dim_rows = []
        for idx in range(1, 5):
            lab = getattr(self, f'lab{idx-1}', None)
            if lab is None:
                continue
            dim_rows.append((idx, lab))
        self._install_hover_map([
            (getattr(self, 'lab1', None), 'Processing settings for dimension 2.'),
            (getattr(self, 'cb_f1180', None), 'Apply first-point correction to dimension 2.'),
            (getattr(self, 'cb_lp1', None), 'Enable linear prediction for dimension 2.'),
            (getattr(self, 'cb_basepol1', None), 'Apply polynomial baseline correction to dimension 2.'),
            (getattr(self, 'p0_1', None), 'Enter the phase correction for dimension 2.'),
            (getattr(self, 'p1_1', None), 'Enter the phase correction for dimension 2.'),
            (getattr(self, 'cb_ft1', None), 'Choose the flip mode for dimension 2.'),
            (getattr(self, 'windowBox1', None), 'Choose GM, SP, or EM apodization for dimension 2.'),
            (getattr(self, 'firstPoint1', None), 'First-point scaling factor for dimension 2.'),
            (getattr(self, 'win2Val1', None), 'Window parameter for dimension 2.'),
            (getattr(self, 'win3Val1', None), 'Window parameter for dimension 2.'),
            (getattr(self, 'lab2', None), 'Processing settings for dimension 3.'),
            (getattr(self, 'cb_f2180', None), 'Apply first-point correction to dimension 3.'),
            (getattr(self, 'cb_lp2', None), 'Enable linear prediction for dimension 3.'),
            (getattr(self, 'cb_basepol2', None), 'Apply polynomial baseline correction to dimension 3.'),
            (getattr(self, 'p0_2', None), 'Enter the phase correction for dimension 3.'),
            (getattr(self, 'p1_2', None), 'Enter the phase correction for dimension 3.'),
            (getattr(self, 'cb_ft2', None), 'Choose the flip mode for dimension 3.'),
            (getattr(self, 'windowBox2', None), 'Choose GM, SP, or EM apodization for dimension 3.'),
            (getattr(self, 'firstPoint2', None), 'First-point scaling factor for dimension 3.'),
            (getattr(self, 'win2Val2', None), 'Window parameter for dimension 3.'),
            (getattr(self, 'win3Val2', None), 'Window parameter for dimension 3.'),
            (getattr(self, 'lab3', None), 'Processing settings for dimension 4.'),
            (getattr(self, 'cb_f3180', None), 'Apply first-point correction to dimension 4.'),
            (getattr(self, 'cb_lp3', None), 'Enable linear prediction for dimension 4.'),
            (getattr(self, 'cb_basepol3', None), 'Apply polynomial baseline correction to dimension 4.'),
            (getattr(self, 'p0_3', None), 'Enter the phase correction for dimension 4.'),
            (getattr(self, 'p1_3', None), 'Enter the phase correction for dimension 4.'),
            (getattr(self, 'cb_ft3', None), 'Choose the flip mode for dimension 4.'),
            (getattr(self, 'windowBox3', None), 'Choose GM, SP, or EM apodization for dimension 4.'),
            (getattr(self, 'firstPoint3', None), 'First-point scaling factor for dimension 4.'),
            (getattr(self, 'win2Val3', None), 'Window parameter for dimension 4.'),
            (getattr(self, 'win3Val3', None), 'Window parameter for dimension 4.'),
        ])

        # Buttons.
        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.showScriptBtn = wx.Button(panel, label='Show Script')
        self.runBtn = wx.Button(panel, label='Run')
        self.closeBtn = wx.Button(panel, label='Close')
        #self.autoApodiseBtn = wx.Button(panel, label='AutoApodise')

        #for btn in (self.autoApodiseBtn,self.showScriptBtn, self.runBtn, self.saveBtn, self.closeBtn):
        for btn in (self.showScriptBtn, self.runBtn, self.closeBtn):
            btn_row.Add(btn, 0, wx.RIGHT, 8)

        self.showScriptBtn.Bind(wx.EVT_BUTTON, self.on_show_script)
        self.runBtn.Bind(wx.EVT_BUTTON, self.on_run)
        self.closeBtn.Bind(wx.EVT_BUTTON, self._on_close)
        #self.autoApodiseBtn.Bind(wx.EVT_BUTTON, self.on_auto_apodise)
        
        root.Add(self.dataSizer2, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 10)

        mode_static = wx.StaticBox(panel, -1, 'Script target')
        mode_box = wx.StaticBoxSizer(mode_static, wx.VERTICAL)
        mode_row = wx.BoxSizer(wx.HORIZONTAL)
        self.scriptModeBox = wx.RadioBox(
            mode_static,
            label='',
            choices=['Process', 'SMILE', 'MDDNMR'],
            majorDimension=1,
            style=wx.RA_SPECIFY_ROWS,
        )
        self.scriptModeBox.Bind(wx.EVT_RADIOBOX, self._on_script_target_change)
        self.maxIterLab = wx.StaticText(mode_static, label='SMILE MaxIter:')
        self.maxIterBox = wx.TextCtrl(mode_static, size=(60, 22))
        self.mddMethodBox = wx.ComboBox(mode_static, value='CS', choices=['CS', 'FT'], style=wx.CB_READONLY, size=(70, 22))
        self.mddAlgorithmBox = wx.ComboBox(mode_static, value='IST', choices=['IST', 'IRLS'], style=wx.CB_READONLY, size=(75, 22))
        self.mddIterBox = wx.TextCtrl(mode_static, value='100', size=(55, 22))
        self.mddVEBox = wx.CheckBox(mode_static, label='Virtual Echo')
        self.mddVEBox.SetValue(True)
        self.ncpusLab = wx.StaticText(mode_static, label='CPUs:')
        self.ncpusBox = wx.TextCtrl(mode_static, size=(60, 22))
        mode_row.Add(self.scriptModeBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        mode_row.Add(self.maxIterLab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        mode_row.Add(self.maxIterBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        mode_row.Add(wx.StaticText(mode_static, label='MDD method:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        mode_row.Add(self.mddMethodBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        mode_row.Add(wx.StaticText(mode_static, label='Algorithm:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        mode_row.Add(self.mddAlgorithmBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        mode_row.Add(wx.StaticText(mode_static, label='Iterations:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        mode_row.Add(self.mddIterBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        mode_row.Add(self.mddVEBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        mode_row.Add(self.ncpusLab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        mode_row.Add(self.ncpusBox, 0, wx.ALIGN_CENTER_VERTICAL)
        mode_box.Add(mode_row, 0, wx.ALL | wx.EXPAND, 6)

        self._install_hover_map([
            (self.baseLab, 'Settings for the direct-dimension processing chain.'),
            (self.cb_baseLin, 'Apply linear baseline correction to the direct dimension.'),
            (self.cb_baseSol, 'Apply digital solvent suppression to the direct dimension.'),
            (self.f1180Lab, 'Apply first-point correction.'),
            (self.lpLab, 'Enable linear prediction.'),
            (self.polyLab, 'Apply polynomial baseline correction.'),
            (self.p0Lab, 'Zero-order phase correction.'),
            (self.p1Lab, 'First-order phase correction.'),
            (self.flipLab, 'Signal polarity / flip handling.'),
            (self.windowLab, 'Select the apodization window function.'),
            (self.windowOp1, 'Primary parameter for the selected window function.'),
            (self.windowOp2, 'Secondary parameter for the selected window function.'),
            (self.firstPointFactor, 'First-point scaling factor used by the window.'),
            (self.lab0, 'Processing settings for dimension 1.'),
            (self.p0, 'Enter the phase correction for dimension 1.'),
            (self.p1, 'Enter the phase correction for dimension 1.'),
            (self.cb_basepol, 'Apply polynomial baseline correction to dimension 1.'),
            (self.cb_ft0, 'Choose the flip mode for dimension 1.'),
            (self.windowBox0, 'Choose GM, SP, or EM apodization for dimension 1.'),
            (self.firstPoint0, 'First-point scaling factor for dimension 1.'),
            (self.win2Val0, 'Window parameter for dimension 1.'),
            (self.win3Val0, 'Window parameter for dimension 1.'),
            (self.showScriptBtn, 'Preview the generated processing script.'),
            (self.runBtn, 'Run the current processing script.'),
            (self.closeBtn, 'Close the processing window.'),
            (self.scriptModeBox, 'Choose whether the generated script targets standard processing or SMILE.'),
            (self.maxIterLab, 'Maximum number of SMILE iterations.'),
            (self.maxIterBox, 'Maximum number of SMILE iterations.'),
            (self.ncpusLab, 'Number of CPU threads used by SMILE processing.'),
            (self.ncpusBox, 'Number of CPU threads used by SMILE processing.'),
        ])
        self._install_default_hover_help()

        root.Add(mode_box, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 10)
        root.Add(btn_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.ALIGN_RIGHT, 10)

        panel.SetSizer(root)
        panel.Layout()
        try:
            panel.Fit()
        except Exception:
            pass
        try:
            self.Fit()
        except Exception:
            pass
        best = panel.GetBestSize()
        frame_w = best.x + 20
        frame_h = best.y + 40
        self.SetMinSize((frame_w, frame_h))
        self.SetClientSize((frame_w, frame_h))
        self.Layout()
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_hover_leave_frame)
        self._refresh_status_text()

    # ------------------------------------------------------------------
    def _safe_set_text(self, ctrl, value):
        if value is None:
            return
        try:
            ctrl.SetValue(str(value))
        except Exception:
            pass

    def _safe_set_combo(self, ctrl, value):
        if value is None:
            return
        try:
            ctrl.SetSelection(int(value))
            return
        except Exception:
            pass
        try:
            ctrl.SetValue(str(value))
        except Exception:
            pass

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

    def _refresh_status_text(self):
        try:
            target = self._script_target_value()
        except Exception:
            target = 'Process'
        text = f'Script target: {target}'
        self._hover_default_status = text
        self._set_status_text(text)

    def _update_reconstruction_controls(self):
        try:
            selection = self.scriptModeBox.GetSelection()
        except Exception:
            selection = 0
        for ctrl in (getattr(self, 'maxIterLab', None), getattr(self, 'maxIterBox', None)):
            if ctrl is not None:
                ctrl.Enable(selection == 1)
        for name in ('mddMethodBox', 'mddAlgorithmBox', 'mddIterBox', 'mddVEBox'):
            ctrl = getattr(self, name, None)
            if ctrl is not None:
                ctrl.Enable(selection == 2)

    def _current_lp_flag(self):
        try:
            return 'y' if self.scriptModeBox.GetSelection() == 1 else ('m' if self.scriptModeBox.GetSelection() == 2 else 'n')
        except Exception:
            return 'n'

    def _script_target_value(self):
        try:
            return str(self.scriptModeBox.GetStringSelection()).strip() or 'Process'
        except Exception:
            return 'Process'

    def _apply_script_target(self, value):
        try:
            text = str(value).strip().lower()
        except Exception:
            text = ''
        selection = 2 if text in {'mddnmr', 'mdd', '2', 'm'} else (1 if text in {'smile process', 'smile', '1', 'y', 'yes', 'lp'} else 0)
        try:
            self.scriptModeBox.SetSelection(selection)
        except Exception:
            pass
        self._update_reconstruction_controls()
        self._refresh_status_text()
        return selection

    def _script_target_path(self, lp=None):
        if lp is None:
            lp = self._current_lp_flag()
        filename = 'nmrprocMDDNMR.com' if lp == 'm' else ('nmrprocLP.com' if lp == 'y' else 'nmrproc.test.com')
        try:
            base = self.proc._spec_output_dir()
        except Exception:
            base = ''
        if not base:
            base = './spec'
        return os.path.join(base, filename)

    def _pipefile_for_mode(self):
        if self.spectral_dim_count in (1,):
            return 'test.ft'
        if self.spectral_dim_count in (2,):
            return 'test.ft2'
        if self.spectral_dim_count == 3:
            return 'test.ft3'
        return 'test.ft4'

    def _render_script_text(self, lp=None):
        # Snapshot first: generated scripts must reflect the widgets visible now,
        # never an older parameter-file value.
        self._processing_updates()
        lp = self._current_lp_flag() if lp is None else lp
        try:
            text, _pipefile = self.proc.RenderProcessScript(lp=lp)
            return _shell_ascii(text)
        except Exception:
            pass
            raise

    def _load_script_text(self, lp=None):
        lp = self._current_lp_flag() if lp is None else lp
        path = self._script_target_path(lp)
        if os.path.exists(path):
            try:
                return Path(path).read_text(), path
            except Exception:
                pass
        return self._render_script_text(lp), path

    def _save_script_text(self, text, lp=None, outfile=None):
        lp = self._current_lp_flag() if lp is None else lp
        path = outfile or self._script_target_path(lp)
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        text = _shell_ascii(text)
        out_path.write_text(text)
        self.script_path = path
        return path

    def _generate_and_save_script(self, lp=None, outfile=None):
        lp = self._current_lp_flag() if lp is None else lp
        text = self._render_script_text(lp)
        path = self._save_script_text(text, lp=lp, outfile=outfile)
        return text, path

    def RenderProcessScript(self, lp='n'):
        return self._render_script_text(lp)

    def WriteProcessScript(self, lp='n', outfile=None):
        lp = self._current_lp_flag() if lp is None else lp
        return self._generate_and_save_script(lp, outfile=outfile)

    def renderprocessscript(self, lp='n'):
        return self.RenderProcessScript(lp=lp)

    def writeprocessscript(self, lp='n', outfile=None):
        return self.WriteProcessScript(lp=lp, outfile=outfile)

    def _run_generated_script(self, lp=None, on_finish=None, title='Processing Output'):
        lp = self._current_lp_flag() if lp is None else lp
        text, script_path = self._generate_and_save_script(lp)
        try:
            self.proc.RunProcessScript(script_path, lp=lp, on_finish=on_finish, title=title)
        except Exception:
            pass
            raise
        return script_path

    def RunProcessScript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        """Proxy the main processing runner from the editor window."""
        runner = getattr(self.proc, 'RunProcessScript', None)
        if runner is None:
            runner = getattr(self.proc, 'runprocessscript', None)
        if runner is None:
            raise AttributeError('Process frame has no RunProcessScript method')
        return runner(script_path, lp=lp, on_finish=on_finish, title=title)

    def runprocessscript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        return self.RunProcessScript(script_path, lp=lp, on_finish=on_finish, title=title)

    def _on_script_target_change(self, event):
        self._update_reconstruction_controls()
        self._refresh_status_text()
        frame = getattr(self, 'script_frame', None)
        if frame is not None:
            try:
                frame.refresh_from_target()
            except Exception:
                pass

    def _on_direct_phase_change(self, event):
        if getattr(self.proc, '_syncing_direct_phase_controls', False):
            if event is not None:
                event.Skip()
            return
        try:
            self.proc._sync_direct_phase_from_controls(refresh_plot=True, sync_sliders=True)
        except Exception:
            pass
        if event is not None:
            event.Skip()

    def _load_from_file(self):
        root = self.proc.parent
        self._apply_script_target(root.Parse(root.deconParFile, 'ProcTarg', default='Process'))
        self._safe_set_text(self.maxIterBox, root.ParseInt(root.deconParFile, 'maxIterSMILE', default=0))
        self._safe_set_combo(self.mddMethodBox, root.Parse(root.deconParFile, 'mddMethod', default='CS'))
        self._safe_set_combo(self.mddAlgorithmBox, root.Parse(root.deconParFile, 'mddAlgorithm', default='IST'))
        self._safe_set_text(self.mddIterBox, root.ParseInt(root.deconParFile, 'mddIterations', default=100))
        self.mddVEBox.SetValue(root.Parse(root.deconParFile, 'mddVirtualEcho', default='y') == 'y')
        self._safe_set_text(self.ncpusBox, root.ParseInt(root.deconParFile, 'ncpus', default=max(1, int(getattr(self.proc, 'ncpus', 1)))))
        self._sync_ncpus_to_process_frame()
        self._safe_set_text(self.p0, root.ParseFlt(root.deconParFile, 'p0'))
        self._safe_set_text(self.p1, root.ParseFlt(root.deconParFile, 'p1'))
        self._safe_set_text(self.win2Val0, root.ParseFlt(root.deconParFile, 'win2Val0', default=20))
        self._safe_set_text(self.win3Val0, root.ParseFlt(root.deconParFile, 'win3Val0', default=2))
        self._safe_set_text(self.firstPoint0, root.ParseFlt(root.deconParFile, 'firstPoint0', default=0.5))
        try:
            self.cb_ft0.SetValue(root.Parse(root.deconParFile, 'flip0', default='Auto'))
        except Exception:
            pass
        self._safe_set_combo(self.windowBox0, root.ParseInt(root.deconParFile, 'window0', default=0))
        self._load_dimension_labels(root)
        self._update_window_labels()
        self._refresh_status_text()

        if root.Parse(root.deconParFile, 'lin', default='n') == 'y':
            self.cb_baseLin.SetValue(True)
        if root.Parse(root.deconParFile, 'poly', default='n') == 'y':
            self.cb_basepol.SetValue(True)
        if root.Parse(root.deconParFile, 'sol', default='n') == 'y':
            self.cb_baseSol.SetValue(True)

        if isinstance(self.spectral_dim_count, int):
            if self.spectral_dim_count >= 2:
                self._safe_set_text(self.p0_1, root.ParseFlt(root.deconParFile, 'p0_1'))
                self._safe_set_text(self.p1_1, root.ParseFlt(root.deconParFile, 'p1_1'))
                self._safe_set_text(self.win2Val1, root.ParseFlt(root.deconParFile, 'win2Val1', default=20))
                self._safe_set_text(self.win3Val1, root.ParseFlt(root.deconParFile, 'win3Val1', default=2))
                self._safe_set_text(self.firstPoint1, root.ParseFlt(root.deconParFile, 'firstPoint1', default=0.5))
                try:
                    self.cb_ft1.SetValue(root.Parse(root.deconParFile, 'flip1', default='Auto'))
                except Exception:
                    pass
                self._safe_set_combo(self.windowBox1, root.ParseInt(root.deconParFile, 'window1', default=0))
                if root.Parse(root.deconParFile, 'f1180', default='n') == 'y':
                    self.cb_f1180.SetValue(True)
                if root.Parse(root.deconParFile, 'lp1', default='n') == 'y':
                    self.cb_lp1.SetValue(True)
                if root.Parse(root.deconParFile, 'bl1', default='n') == 'y':
                    self.cb_basepol1.SetValue(True)
            if self.spectral_dim_count >= 3:
                self._safe_set_text(self.p0_2, root.ParseFlt(root.deconParFile, 'p0_2'))
                self._safe_set_text(self.p1_2, root.ParseFlt(root.deconParFile, 'p1_2'))
                self._safe_set_text(self.win2Val2, root.ParseFlt(root.deconParFile, 'win2Val2', default=20))
                self._safe_set_text(self.win3Val2, root.ParseFlt(root.deconParFile, 'win3Val2', default=2))
                self._safe_set_text(self.firstPoint2, root.ParseFlt(root.deconParFile, 'firstPoint2', default=0.5))
                try:
                    self.cb_ft2.SetValue(root.Parse(root.deconParFile, 'flip2', default='Auto'))
                except Exception:
                    pass
                self._safe_set_combo(self.windowBox2, root.ParseInt(root.deconParFile, 'window2', default=0))
                if root.Parse(root.deconParFile, 'f2180', default='n') == 'y':
                    self.cb_f2180.SetValue(True)
                if root.Parse(root.deconParFile, 'lp2', default='n') == 'y':
                    self.cb_lp2.SetValue(True)
                if root.Parse(root.deconParFile, 'bl2', default='n') == 'y':
                    self.cb_basepol2.SetValue(True)
            if self.spectral_dim_count == 4:
                self._safe_set_text(self.p0_3, root.ParseFlt(root.deconParFile, 'p0_3'))
                self._safe_set_text(self.p1_3, root.ParseFlt(root.deconParFile, 'p1_3'))
                self._safe_set_text(self.win2Val3, root.ParseFlt(root.deconParFile, 'win2Val3', default=20))
                self._safe_set_text(self.win3Val3, root.ParseFlt(root.deconParFile, 'win3Val3', default=2))
                self._safe_set_text(self.firstPoint3, root.ParseFlt(root.deconParFile, 'firstPoint3', default=0.5))
                try:
                    self.cb_ft3.SetValue(root.Parse(root.deconParFile, 'flip3', default='Auto'))
                except Exception:
                    pass
                self._safe_set_combo(self.windowBox3, root.ParseInt(root.deconParFile, 'window3', default=0))
                if root.Parse(root.deconParFile, 'f3180', default='n') == 'y':
                    self.cb_f3180.SetValue(True)
                if root.Parse(root.deconParFile, 'lp3', default='n') == 'y':
                    self.cb_lp3.SetValue(True)
                if root.Parse(root.deconParFile, 'bl3', default='n') == 'y':
                    self.cb_basepol3.SetValue(True)
        elif self.spectral_dim_count == 2 and self.has_pseudo_axis:
            self._safe_set_text(self.p0_1, root.ParseFlt(root.deconParFile, 'p0_1'))
            self._safe_set_text(self.p1_1, root.ParseFlt(root.deconParFile, 'p1_1'))
            self._safe_set_text(self.win2Val1, root.ParseFlt(root.deconParFile, 'win2Val1', default=20))
            self._safe_set_text(self.win3Val1, root.ParseFlt(root.deconParFile, 'win3Val1', default=2))
            self._safe_set_text(self.firstPoint1, root.ParseFlt(root.deconParFile, 'firstPoint1', default=0.5))
            try:
                self.cb_ft1.SetValue(root.Parse(root.deconParFile, 'flip1', default='Auto'))
            except Exception:
                pass
            self._safe_set_combo(self.windowBox1, root.ParseInt(root.deconParFile, 'window1', default=0))
            if root.Parse(root.deconParFile, 'f1180', default='n') == 'y':
                self.cb_f1180.SetValue(True)
            if root.Parse(root.deconParFile, 'lp1', default='n') == 'y':
                self.cb_lp1.SetValue(True)
            if root.Parse(root.deconParFile, 'bl1', default='n') == 'y':
                self.cb_basepol1.SetValue(True)

        self._update_window_labels()
        self._refresh_status_text()

    def reload_from_file(self):
        """Reload widget values from the current system save file.

        This is used by the auto-run path so the hidden processing window picks
        up the same saved defaults as the visible manual window.
        """
        self._load_from_file()

    def _ncpus_value(self, default=1):
        ctrl = getattr(self, 'ncpusBox', None)
        if ctrl is not None:
            try:
                return max(1, int(float(ctrl.GetValue())))
            except Exception:
                pass
        try:
            return max(1, int(float(getattr(self.proc, 'ncpus', default))))
        except Exception:
            return max(1, int(default))

    def _sync_ncpus_to_process_frame(self):
        value = self._ncpus_value(default=1)
        try:
            self.proc.ncpus = value
        except Exception:
            pass
        try:
            parent = getattr(self.proc, 'parent', None)
            core_box = getattr(parent, 'coreBox', None)
            if core_box is not None:
                core_box.SetValue(str(value))
        except Exception:
            pass
        return value

    def _load_dimension_labels(self, root):
        getter = getattr(self.proc, 'get_dimension_labels', None)
        labels = list(getter()) if callable(getter) else []
        labels += [f'dim {i}' for i in range(len(labels) + 1, 5)]
        for ctrl_name, label in zip(('lab0', 'lab1', 'lab2', 'lab3'), labels):
            ctrl = getattr(self, ctrl_name, None)
            if ctrl is not None:
                try:
                    ctrl.SetLabel(label)
                except Exception:
                    pass

    def _update_window_labels(self):
        vals = [self.windowBox0.GetValue()]
        if hasattr(self, 'windowBox1'):
            vals.append(self.windowBox1.GetValue())
        if hasattr(self, 'windowBox2'):
            vals.append(self.windowBox2.GetValue())
        if hasattr(self, 'windowBox3'):
            vals.append(self.windowBox3.GetValue())

        if len(vals) > 0 and all(v == vals[0] for v in vals):
            if vals[0] == 'GM':
                self.windowOp1.SetLabel('Gauss')
                self.windowOp2.SetLabel('Anti-Lor')
            elif vals[0] == 'SP':
                self.windowOp1.SetLabel('Power')
                self.windowOp2.SetLabel('Offset')
            elif vals[0] == 'EM':
                self.windowOp1.SetLabel('LB')
                self.windowOp2.SetLabel('N/A')
            else:
                self.windowOp1.SetLabel('Mixed')
                self.windowOp2.SetLabel('Mixed')
        else:
            self.windowOp1.SetLabel('Mixed')
            self.windowOp2.SetLabel('Mixed')

    def set_win(self, event):
        self._update_window_labels()
        self._refresh_status_text()

    # ------------------------------------------------------------------
    # Actions
    def on_show_script(self, event):
        try:
            self._sync_ncpus_to_process_frame()
            script_path = self._script_target_path(self._current_lp_flag())
            if getattr(self, 'script_frame', None) is not None:
                try:
                    self.script_frame.Raise()
                    self.script_frame.SetFocus()
                    return
                except Exception:
                    self.script_frame = None
            self.script_frame = ProcessingScriptFrame(self, script_path=script_path)
            self.script_frame.Show(True)
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def on_run(self, event, on_finish=None):
        self.proc.save_current_gui_state(reason='processing-run')
        self._sync_ncpus_to_process_frame()
        self._refresh_status_text()
        try:
            self._run_generated_script(on_finish=on_finish)
        except Exception:
            errorMessage('Could not run processing script.')

    def _processing_updates(self, update_state=True):
        """Return the settings owned by this window in parameter-file form.

        Keep this list symmetric with _load_from_file(): every editable setting
        loaded by ProcessingFrame is written back by ProcessingFrame.  This
        avoids relying on ProcessFrame's deliberately narrower Save action.
        """
        def text(name, default=''):
            ctrl = getattr(self, name, None)
            if ctrl is None:
                return default
            try:
                return str(ctrl.GetValue()).strip()
            except Exception:
                return default

        def checked(name):
            ctrl = getattr(self, name, None)
            try:
                return 'y' if ctrl is not None and ctrl.GetValue() else 'n'
            except Exception:
                return 'n'

        def selection(name, default=0):
            ctrl = getattr(self, name, None)
            try:
                value = int(ctrl.GetSelection())
                return str(value if value >= 0 else default)
            except Exception:
                return str(default)

        updates = {
            'ProcTarg': self._script_target_value(),
            'maxIterSMILE': text('maxIterBox', '0'),
            'mddMethod': text('mddMethodBox', 'CS'),
            'mddAlgorithm': text('mddAlgorithmBox', 'IST'),
            'mddIterations': text('mddIterBox', '100'),
            'mddVirtualEcho': checked('mddVEBox'),
            'ncpus': text('ncpusBox', '1'),
            'p0': text('p0', '0'),
            'p1': text('p1', '0'),
            'flip0': text('cb_ft0', 'Auto'),
            'window0': selection('windowBox0'),
            'win2Val0': text('win2Val0', '20'),
            'win3Val0': text('win3Val0', '2'),
            'firstPoint0': text('firstPoint0', '0.5'),
            'lin': checked('cb_baseLin'),
            'poly': checked('cb_basepol'),
            'sol': checked('cb_baseSol'),
        }

        # Indirect dimensions use the historical 1-based suffixes.
        for idx in range(1, 4):
            if not hasattr(self, f'windowBox{idx}'):
                continue
            updates.update({
                f'p0_{idx}': text(f'p0_{idx}', '0'),
                f'p1_{idx}': text(f'p1_{idx}', '0'),
                f'flip{idx}': text(f'cb_ft{idx}', 'Auto'),
                f'window{idx}': selection(f'windowBox{idx}'),
                f'win2Val{idx}': text(f'win2Val{idx}', '20'),
                f'win3Val{idx}': text(f'win3Val{idx}', '2'),
                f'firstPoint{idx}': text(f'firstPoint{idx}', '0.5'),
                f'f{idx}180': checked(f'cb_f{idx}180'),
                f'lp{idx}': checked(f'cb_lp{idx}'),
                f'bl{idx}': checked(f'cb_basepol{idx}'),
            })
        state = getattr(self.proc, 'state', None)
        if state is not None and update_state:
            state.update_gui_settings(updates)
        return updates

    collect_updates = _processing_updates

    def _on_close(self, event):
        try:
            self._processing_updates()
        except Exception:
            pass
        try:
            if getattr(self, 'script_frame', None) is not None:
                try:
                    self.script_frame.Destroy()
                except Exception:
                    pass
                self.script_frame = None
            self.proc.processing_frame = None
            self.proc._unbind_processing_controls()
        except Exception:
            pass
        try:
            if getattr(self, '_autoapodise_progress_frame', None) is not None:
                try:
                    self._autoapodise_progress_frame.Destroy()
                except Exception:
                    pass
                self._autoapodise_progress_frame = None
        except Exception:
            pass
        self.Destroy()


    def on_auto_apodise(self, event):
        #print("que?")
        progress = None
        try:
            progress = ShellOutputFrame(self, title='AutoApodise Progress')
            progress.set_status('Starting automatic apodization...')
            progress.append_text('Starting automatic apodization...\n')
            progress.Show()
            self._autoapodise_progress_frame = progress
        except Exception:
            import traceback
            pass
            progress = None
        try:
            handler = getattr(self.proc, 'OnAutoApodise', None)
            if handler is None:
                handler = getattr(self.proc, 'AutoApodise', None)
            if handler is None:
                raise AttributeError('Process frame has no automatic apodization method')
            try:
                handler(event=event, progress_frame=progress)
            except TypeError:
                handler(event)
        except Exception as exc:
            import traceback
            pass
            pass
            errorMessage(str(exc))
            if progress is not None:
                try:
                    progress.append_text('Automatic apodization failed to start.\n')
                    progress.set_status('Failed')
                except Exception:
                    import traceback
                    pass
        finally:
            self._autoapodise_progress_frame = progress



        
class ProcessingScriptFrame(wx.Frame):
    """Editable processing script window for nmrproc.test.com / nmrprocLP.com."""

    def __init__(self, proc, script_path=None):
        super().__init__(proc, title='nmrproc.test.com', size=(980, 520))
        self.proc = proc
        self.panel = None
        self.path_label = None
        self.scriptBox = None
        self.script_path = script_path or self._script_target_path(self._current_lp_flag())

        self._build_ui()
        self._load_script()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _update_reconstruction_controls(self):
        try:
            selection = self.scriptModeBox.GetSelection()
        except Exception:
            selection = 0
        for ctrl in (getattr(self, 'maxIterLab', None), getattr(self, 'maxIterBox', None)):
            if ctrl is not None:
                ctrl.Enable(selection == 1)
        for name in ('mddMethodBox', 'mddAlgorithmBox', 'mddIterBox', 'mddVEBox'):
            ctrl = getattr(self, name, None)
            if ctrl is not None:
                ctrl.Enable(selection == 2)

    def _current_lp_flag(self):
        try:
            return self.proc._current_lp_flag()
        except Exception:
            return 'n'

    def _script_target_path(self, lp=None):
        try:
            return self.proc._script_target_path(lp)
        except Exception:
            filename = 'nmrprocLP.com' if (lp or self._current_lp_flag()) == 'y' else 'nmrproc.test.com'
            try:
                base = self.proc._spec_output_dir()
            except Exception:
                base = ''
            if not base:
                base = './spec'
            return os.path.join(base, filename)

    def _current_path_label(self):
        try:
            return self.script_path or self._script_target_path(self._current_lp_flag())
        except Exception:
            return 'nmrproc.test.com'

    def _build_ui(self):
        self.panel = wx.Panel(self)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        root = wx.BoxSizer(wx.VERTICAL)

        self.path_label = wx.StaticText(self.panel, label=self._current_path_label())
        root.Add(self.path_label, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)

        self.scriptBox = wx.TextCtrl(self.panel, style=wx.TE_MULTILINE | wx.TE_RICH | wx.TE_DONTWRAP)
        self.scriptBox.SetMinSize((920, 360))
        root.Add(self.scriptBox, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 10)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.guessBtn = wx.Button(self.panel, label='Guess')
        self.saveBtn = wx.Button(self.panel, label='Save')
        self.runBtn = wx.Button(self.panel, label='Run')
        self.closeBtn = wx.Button(self.panel, label='Close')
        for btn in (self.guessBtn, self.saveBtn, self.runBtn, self.closeBtn):
            btn_row.Add(btn, 0, wx.RIGHT, 6)

        self.guessBtn.Bind(wx.EVT_BUTTON, self.on_guess)
        self.saveBtn.Bind(wx.EVT_BUTTON, self.on_save)
        self.runBtn.Bind(wx.EVT_BUTTON, self.on_run)
        self.closeBtn.Bind(wx.EVT_BUTTON, self._on_close)
        
        root.Add(btn_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.ALIGN_RIGHT, 10)

        self.panel.SetSizer(root)
        self.panel.Layout()
        try:
            self.panel.Fit()
        except Exception:
            pass
        best = self.panel.GetBestSize()
        self.SetMinSize(best)
        self.SetClientSize(best)
        self.Layout()

    def _load_script(self):
        lp = self._current_lp_flag()
        path = self.script_path or self._script_target_path(lp)
        self.script_path = path
        try:
            if path and os.path.exists(path):
                self.path_label.SetLabel(os.path.basename(path))
                self.scriptBox.SetValue(Path(path).read_text())
                return
        except Exception:
            pass
        try:
            text, saved_path = self.proc._generate_and_save_script(lp, outfile=path)
            self.script_path = saved_path
            self.path_label.SetLabel(os.path.basename(saved_path))
            self.scriptBox.SetValue(text)
        except Exception:
            try:
                text = self.proc._render_script_text(lp)
                self.path_label.SetLabel(os.path.basename(path))
                self.scriptBox.SetValue(text)
            except Exception as exc:
                pass
                errorMessage(str(exc))

    def refresh_from_target(self):
        self._load_script()

    def on_guess(self, event):
        try:
            self.proc.proc.save_current_gui_state(reason='processing-script-guess')
            lp = self._current_lp_flag()
            target_path = self.script_path or self._script_target_path(lp)
            text, saved_path = self.proc._generate_and_save_script(lp, outfile=target_path)
            self.script_path = saved_path
            self.path_label.SetLabel(os.path.basename(saved_path))
            self.scriptBox.SetValue(text)
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def on_save(self, event):
        try:
            if not self.script_path:
                self.script_path = self._script_target_path(self._current_lp_flag())
            out_path = Path(self.script_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            clean_text = _shell_ascii(self.scriptBox.GetValue())
            # Keep editor and disk identical so Run cannot reintroduce smart punctuation.
            if clean_text != self.scriptBox.GetValue():
                self.scriptBox.SetValue(clean_text)
            out_path.write_text(clean_text)
            self.path_label.SetLabel(os.path.basename(self.script_path))
            self.proc.script_path = self.script_path
            return self.script_path
        except Exception as exc:
            pass
            errorMessage(str(exc))
            return None

    def on_run(self, event):
        try:
            script_path = self.on_save(event)
            self.proc.proc.save_current_gui_state(reason='processing-script-run')
            if not script_path:
                return
            runner = getattr(self.proc, 'RunProcessScript', None)
            if runner is None:
                runner = getattr(self.proc, 'runprocessscript', None)
            if runner is None:
                raise AttributeError('Process frame has no RunProcessScript method')
            runner(script_path, lp=self._current_lp_flag())
        except Exception as exc:
            pass
            errorMessage(str(exc))

    def _on_close(self, event):
        try:
            if getattr(self.proc, 'script_frame', None) is self:
                self.proc.script_frame = None
        except Exception:
            pass
        try:
            if getattr(self, '_autoapodise_progress_frame', None) is not None:
                try:
                    self._autoapodise_progress_frame.Destroy()
                except Exception:
                    pass
                self._autoapodise_progress_frame = None
        except Exception:
            pass
        self.Destroy()
