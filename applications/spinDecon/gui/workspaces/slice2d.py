#!/usr/bin/python
"""
Academic Use Licence

These licence terms apply to all licences granted by THE CHANCELLOR, MASTERS AND SCHOLARS OF THE UNIVERSITY OF OXFORD whose administrative offices are at University Offices, Wellington Square, Oxford OX1 2JD, United Kingdom (the "University") for use of UniDecNMR ("the Software") downloaded from the following website: https://github.com/charliebuchanan/UniDecNMR ("the Website")
By downloading the Software through the Source, you (the "Licensee") are confirming that you agree that your use of the Software is subject to these licence terms.

PLEASE READ THESE LICENCE TERMS CAREFULLY BEFORE DOWNLOADING THE SOFTWARE THROUGH THIS WEBSITE.  IF YOU DO NOT AGREE TO THESE LICENCE TERMS YOU SHOULD NOT DOWNLOAD THE SOFTWARE.

THE SOFTWARE IS INTENDED FOR USE BY ACADEMICS CARRYING OUT RESEARCH AND NOT FOR USE BY CONSUMERS OR COMMERCIAL BUSINESSES.

1.	Academic Use Licence
1.1	The Licensee is granted a limited non-exclusive and non-transferable royalty free licence to download and use the Software provided that the Licensee will:
(a)	limit their use of the Software to their own internal academic non-commercial research which is undertaken for the purposes of education or other scholarly use; 
(b)	not use the Software for or on behalf of any third party or to provide a service or integrate all or part of the Software into a product for sale or license to third parties;
(c)	use the Software in accordance with the prevailing instructions and guidance for use given on the Website and comply with procedures on the Website for user identification, authentication and access;
(d)	comply with all applicable laws and regulations with respect to their use of the Software; and 
(e)	ensure that the Copyright Notice "Copyright (c) 2022, University of Oxford" appears prominently wherever the Software is reproduced and on any documents or other material created using the Software.
1.2	The Licensee may only reproduce, modify, transmit or transfer the Software where:
(a)	such reproduction, modification, transmission or transfer is for academic, research or other scholarly use;
(b)	the conditions of this Licence are imposed upon the receiver of the Software or any modified Software;
(c)	all original and modified Source Code is included in any transmitted software program; and
(d)	the Licensee grants the University an irrevocable, indefinite, royalty free, non-exclusive unlimited licence to use and sub-licence any modified Source Code as part of the Software.

1.3	The University reserves the right at any time and without liability or prior notice to the Licensee to revise, modify and replace the functionality and performance of the access to and operation of the Software.
1.4	The Licensee acknowledges and agrees that the University owns all intellectual property rights in the Software.  The Licensee shall not have any right, title or interest in the Software.
1.5	This Licence will terminate immediately and the Licensee will no longer have any right to use the Software or exercise any of the rights granted to the Licensee upon any breach of the conditions in Section 1 of this Licence.

2.	Indemnity and Liability 
2.1	The Licensee shall defend, indemnify and hold harmless the University against any claims, actions, proceedings, losses, damages, expenses and costs (including without limitation court costs and reasonable legal fees) arising out of or in connection with the Licensee's possession or use of the Software, or any breach of these terms by the Licensee. 
2.2	The Software is provided on an 'as is' basis and the Licensee uses the Software at their own risk. No representations, conditions, warranties or other terms of any kind are given in respect of the the Software and all statutory warranties and conditions are excluded to the fullest extent permitted by law. Without affecting the generality of the previous sentences, the University gives no implied or express warranty and makes no representation that the Software or any part of the Software: (a) will enable specific results to be obtained; or (b) meets a particular specification or is comprehensive within its field or that it is error free or will operate without interruption; or (c) is suitable for any particular, or the Licensee's specific purposes. 
2.3	Except in relation to fraud, death or personal injury, the University's liability to the Licensee for any use of the Software, in negligence or arising in any other way out of the subject matter of these licence terms, will not extend to any incidental or consequential damages or losses, or any loss of profits, loss of revenue, loss of data, loss of contracts or opportunity, whether direct or indirect.
2.4	The Licensee hereby irrevocably undertakes to the University not to make any claim against any employee, student, researcher or other individual engaged by the University, being a claim which seeks to enforce against any of them any liability whatsoever in connection with these licence terms or their subject-matter. 

3.	General 
3.1	Severability - If any provision (or part of a provision) of these licence terms is found by any court or administrative body of competent jurisdiction to be invalid, unenforceable or illegal, the other provisions shall remain in force.
3.2	Entire Agreement - These licence terms constitute the whole agreement between the parties and supersede any previous arrangement, understanding or agreement between them relating to the Software. 
3.3	Law and Jurisdiction - These licence terms and any disputes or claims arising out of or in connection with them shall be governed by, and construed in accordance with, the law of England. The Licensee irrevocably submits to the exclusive jurisdiction of the English courts for any dispute or claim that arises out of or in connection with these licence terms.

If you are interested in using the Software commercially, please contact Oxford University Innovation Limited to negotiate a licence. Contact details are enquiries@innovation.ox.ac.uk 

"""
import wx
from spinDecon.domain.dimensions.viewer_contract import spectral_dim_count
from spinDecon.gui.context import context_for, project_for
import string
import copy
import math
import numpy
import os
import re
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar

class _ToolbarToggleState:
    """Small wx.CheckBox-compatible state holder for toolbar-owned toggles."""
    def __init__(self, value=False):
        self._value = bool(value)
        self._enabled = True
    def GetValue(self): return self._value
    def IsChecked(self): return self._value
    def SetValue(self, value): self._value = bool(value)
    def Enable(self, enabled=True): self._enabled = bool(enabled)
    def IsEnabled(self): return self._enabled

import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from wx.lib.mixins.listctrl import ColumnSorterMixin
from spinDecon.project.parameter_store import parse_float as ParseFlt
from spinDecon.gui.workspaces.slice1d import _scientific_unicode
import logging
from spinDecon.gui.dialogs.errors import errorMessage
from spinDecon.gui.widgets.common import PersistentStateButton


##########################################################################
# 2D plotting of NMR slices
#

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

class SliceFrame2D(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent,tabOne):
        wx.Panel.__init__(self, parent=parent)


        self.parent=parent
        self.app_context = context_for(tabOne, parent)
        self.slice_service = (getattr(self.app_context, "slices", None) if self.app_context is not None else None) or SliceService(tabOne)
        self.state = project_for(tabOne, parent)
        self.sym = self.slice_service.symmetry_enabled
        # Copy in the previous variables through the slice boundary.
        self.index_data=self.index(self.slice_service.reference_peaks())
        self.thresh = self.slice_service.threshold()
        #self.offset=copy.deepcopy(tabOne.offset)
        self.offset=0.0                                #
        # Peak data remain in decon_tab/DataStore; this viewer resolves them on demand.

        self.spectrumfile = self.slice_service.spectrum_path
        #get 2d strips from 3d data
        self.GetSlice2d()             #slice up the 2D spectrum


        #
        self.ax_reset1=1       #for keeping the zoom
        self.ax_reset2=1
        self.ax_resetCC=1       #for keeping the zoom
        self.ax_resetHC=1
        self.inc=0            #for incrementing the slices
        self.inc2=0
        self.selection=[]
        # Full-peak editing state.  This is GUI/history state only: the authoritative
        # peak records always live in the authoritative Full Peak List.
        self.full_tool_mode = None
        self.full_selected_name = None
        # Pane that initiated the current Slice2D selection.  Peak identity
        # remains global/authoritative; this is presentation state only.
        self.full_selection_pane = None
        self.full_undo_stack = []
        self.full_redo_stack = []
        self._full_plane = None

        self.create_main_panel()
        #print('22')
        self.draw_figure()
        #print('44')
        self.Show(True)
        self.Fit()

    @property
    def peak(self):
        return self.slice_service.reference_peaks()

    def _make_modeless_window(self, title):
        """Create a modeless tool window owned by the 2D-slices panel."""
        frame = wx.Frame(self.GetTopLevelParent(), title=title,
                         style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(frame)
        frame.Bind(wx.EVT_CLOSE, lambda evt, f=frame: (f.Hide(), evt.Veto()))
        return frame, panel

    def _close_tool_window(self, frame):
        frame.Hide()

    def _show_tool_window(self, frame):
        # Keep the Full 3D peak-tools palette docked visually to the right of
        # the application's main top-level window.  Other modeless settings
        # windows retain their existing placement behaviour.
        if frame is getattr(self, 'fullToolsFrame', None):
            owner = wx.GetTopLevelParent(self)
            if owner is not None:
                pos = owner.GetScreenPosition()
                size = owner.GetSize()
                x = pos.x + size.width
                y = pos.y
                # Keep the palette on the current display when the owner is
                # close to a screen edge.  Prefer the requested right-hand
                # position, falling back inside the display only as needed.
                display_index = wx.Display.GetFromWindow(owner)
                if display_index != wx.NOT_FOUND:
                    area = wx.Display(display_index).GetClientArea()
                    tool_size = frame.GetSize()
                    x = min(x, area.GetRight() - tool_size.width + 1)
                    x = max(x, area.GetLeft())
                    y = min(y, area.GetBottom() - tool_size.height + 1)
                    y = max(y, area.GetTop())
                frame.SetPosition((x, y))
        if not frame.IsShown():
            frame.Show()
        frame.Raise()

    def _contour_group(self, parent, prefix, min_value=None):
        box = wx.StaticBox(parent, label='Contours:')
        sz = wx.StaticBoxSizer(box, wx.HORIZONTAL)
        labels = []
        ctrls = []
        for label, width in [('Min:', 82), ('Fac:', 45), ('Num:', 45)]:
            lab = wx.StaticText(box, label=label)
            ctrl = wx.TextCtrl(box, size=(width, 22), style=wx.TE_PROCESS_ENTER)
            labels.append(lab); ctrls.append(ctrl)
            sz.Add(lab, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 4)
            sz.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 3)
        ctrls[0].SetValue(str(self.thresh if min_value is None else min_value))
        ctrls[1].SetValue('1.2')
        ctrls[2].SetValue('15')
        setattr(self, prefix + '_contour_box', box)
        return sz, ctrls

    def contour_boxes(self):
        # Main contours now belong to the independent Left/Right tool windows.
        # This method creates only the projection-contour modeless window.
        self.projContourFrame, panel = self._make_modeless_window('Projection Contours')
        contourSizer, ctrls = self._contour_group(panel, 'projection')
        self.textbox_minP, self.textbox_maxP, self.textbox_lvlP = ctrls
        for ctrl in ctrls:
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_projection_setting_changed)
        close = wx.Button(panel, label='Close', size=(-1, 24))
        close.Bind(wx.EVT_BUTTON, lambda evt: self._close_tool_window(self.projContourFrame))
        root = wx.BoxSizer(wx.VERTICAL)
        root.Add(contourSizer, 0, wx.EXPAND | wx.ALL, 5)
        root.Add(close, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        panel.SetSizer(root)
        self.projContourFrame.SetClientSize(root.CalcMin())

    def left_box(self):
        self.leftFrame, panel = self._make_modeless_window('Left')
        flags = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        self.width1Lab = wx.StaticText(panel, label='Width (ppm):')
        self.width1Box = wx.TextCtrl(panel, size=(48, 22), style=wx.TE_PROCESS_ENTER)
        self.width1Box.SetValue(str(self.slice_service.peak_shape_width(3)))
        self.listy = self.slice_service.peak_names()
        self.Nbutton = wx.Button(panel, label='Next', size=(-1, 22))
        self.Pbutton = wx.Button(panel, label='Previous', size=(-1, 22))
        self.upDownLab = wx.StaticText(panel, label='Slice:')
        self.sliceSpin1 = wx.SpinCtrl(panel, value='0', min=-99999, max=99999, initial=0, size=(72, 24))
        self.ComboBox1 = wx.ComboBox(panel, size=(90, 22), choices=self.listy, style=wx.CB_READONLY)
        if self.ComboBox1.GetCount(): self.ComboBox1.SetSelection(0)
        self.ComboBox1.Bind(wx.EVT_COMBOBOX, self.on_peak_combo_changed)
        self.Nbutton.Bind(wx.EVT_BUTTON, self.on_N_button)
        self.Pbutton.Bind(wx.EVT_BUTTON, self.on_P_button)
        self.sliceSpin1.Bind(wx.EVT_SPINCTRL, self.on_slice_spin_left)
        row1 = wx.BoxSizer(wx.HORIZONTAL)
        for w in (self.ComboBox1, self.Pbutton, self.Nbutton): row1.Add(w, 0, flags, 3)
        row2 = wx.BoxSizer(wx.HORIZONTAL); row2.Add(self.width1Lab,0,flags,3); row2.Add(self.width1Box,0,flags,3)
        row3 = wx.BoxSizer(wx.HORIZONTAL); row3.Add(self.upDownLab,0,flags,3); row3.Add(self.sliceSpin1,0,flags,3)
        contourSizer, ctrls = self._contour_group(panel, 'left')
        self.textbox0, self.textbox1, self.textbox2 = ctrls  # compatibility aliases: Left/main contour values
        self.width1Box.Bind(wx.EVT_TEXT_ENTER, self.on_left_right_setting_changed)
        for ctrl in ctrls:
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_left_right_setting_changed)
        close = wx.Button(panel, label='Close', size=(-1,24)); close.Bind(wx.EVT_BUTTON, lambda evt: self._close_tool_window(self.leftFrame))
        root=wx.BoxSizer(wx.VERTICAL)
        root.Add(row1,0,wx.ALL,2); root.Add(row2,0,wx.LEFT|wx.RIGHT,2); root.Add(row3,0,wx.LEFT|wx.RIGHT|wx.BOTTOM,2)
        root.Add(contourSizer,0,wx.EXPAND|wx.ALL,4); root.Add(close,0,wx.ALIGN_RIGHT|wx.ALL,5)
        panel.SetSizer(root); self.leftFrame.SetClientSize(root.CalcMin())

    def right_box(self):
        self.rightFrame, panel = self._make_modeless_window('Right')
        flags = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        self.width2Lab = wx.StaticText(panel, label='Width (ppm):')
        self.width2Box = wx.TextCtrl(panel, size=(48,22), style=wx.TE_PROCESS_ENTER); self.width2Box.SetValue(str(self.slice_service.peak_shape_width(2)))
        self.ComboBox2 = wx.ComboBox(panel, size=(90,22), choices=self.listy, style=wx.CB_READONLY)
        if self.ComboBox2.GetCount(): self.ComboBox2.SetSelection(0)
        self.Nbutton2=wx.Button(panel,label='Next',size=(-1,22)); self.Pbutton2=wx.Button(panel,label='Previous',size=(-1,22))
        self.sliceSpin2 = wx.SpinCtrl(panel, value='0', min=-99999, max=99999, initial=0, size=(72,24))
        self.swapbutton=wx.Button(panel,label='Swap',size=(-1,22)); self.upDownLab2=wx.StaticText(panel,label='Slice:')
        self.ComboBox2.Bind(wx.EVT_COMBOBOX, self.on_peak_combo_changed)
        self.swapbutton.Bind(wx.EVT_BUTTON, self.on_swap_button)
        self.Nbutton2.Bind(wx.EVT_BUTTON, self.on_N_button2)
        self.Pbutton2.Bind(wx.EVT_BUTTON, self.on_P_button2)
        self.sliceSpin2.Bind(wx.EVT_SPINCTRL, self.on_slice_spin_right)
        row1=wx.BoxSizer(wx.HORIZONTAL)
        for w in (self.ComboBox2,self.Pbutton2,self.Nbutton2): row1.Add(w,0,flags,3)
        row2=wx.BoxSizer(wx.HORIZONTAL); row2.Add(self.width2Lab,0,flags,3); row2.Add(self.width2Box,0,flags,3); row2.Add(self.swapbutton,0,flags,3)
        row3=wx.BoxSizer(wx.HORIZONTAL); row3.Add(self.upDownLab2,0,flags,3); row3.Add(self.sliceSpin2,0,flags,3)
        contourSizer, ctrls = self._contour_group(panel, 'right')
        self.textbox0_right, self.textbox1_right, self.textbox2_right = ctrls
        self.width2Box.Bind(wx.EVT_TEXT_ENTER, self.on_left_right_setting_changed)
        for ctrl in ctrls:
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_left_right_setting_changed)
        close=wx.Button(panel,label='Close',size=(-1,24)); close.Bind(wx.EVT_BUTTON,lambda evt:self._close_tool_window(self.rightFrame))
        root=wx.BoxSizer(wx.VERTICAL); root.Add(row1,0,wx.ALL,2); root.Add(row2,0,wx.LEFT|wx.RIGHT,2); root.Add(row3,0,wx.LEFT|wx.RIGHT|wx.BOTTOM,2)
        root.Add(contourSizer,0,wx.EXPAND|wx.ALL,4); root.Add(close,0,wx.ALIGN_RIGHT|wx.ALL,5)
        panel.SetSizer(root); self.rightFrame.SetClientSize(root.CalcMin())

    def control_box(self):
        # These controls remain owned/referenced by this panel, but are placed
        # directly on the compact Matplotlib toolbar row.
        # Orthogonal, Deconvolution and 1D are native check tools on the
        # Matplotlib toolbar.  Keep these lightweight state objects because
        # the mature plotting code reads them through GetValue().
        self.cb_flip=_ToolbarToggleState(False); self.cb_decon=_ToolbarToggleState(False)
        self.cb_1d=_ToolbarToggleState(False); self.cb_grid_auto=wx.CheckBox(self,label='Labels')
        # State-only compatibility control: peak visibility is now owned by the Matplotlib toolbar.
        self.cb_grid_auto.Hide()
        self.Bind(wx.EVT_CHECKBOX,self.on_cb_grid_auto,self.cb_grid_auto)
        self.cb_grid_auto.SetValue(1)
        self.leftToolButton=wx.Button(self,label='Left',size=(-1,24)); self.rightToolButton=wx.Button(self,label='Right',size=(-1,24))
        # Compact main-toolbar arrows navigate the focused/reference peak, not
        # the Z-slice offset.  The modeless Left/Right SpinCtrls remain the
        # dedicated inc/inc2 controls for stepping through neighbouring planes.
        self.quickSliceSpin1 = wx.SpinButton(self, size=(22,24), style=wx.SP_VERTICAL)
        self.quickSliceSpin2 = wx.SpinButton(self, size=(22,24), style=wx.SP_VERTICAL)
        self.quickSliceSpin1.SetRange(-99999, 99999)
        self.quickSliceSpin2.SetRange(-99999, 99999)
        self.quickSliceSpin1.Bind(wx.EVT_SPIN_UP, self.on_quick_peak_up_left)
        self.quickSliceSpin1.Bind(wx.EVT_SPIN_DOWN, self.on_quick_peak_down_left)
        self.quickSliceSpin2.Bind(wx.EVT_SPIN_UP, self.on_quick_peak_up_right)
        self.quickSliceSpin2.Bind(wx.EVT_SPIN_DOWN, self.on_quick_peak_down_right)
        self.projContourButton=wx.Button(self,label='Projection Contours',size=(-1,24)); self.peaksToolButton=wx.Button(self,label='Peaks',size=(-1,24))
        # Projection contours moved to the Matplotlib contour tool.  Keep the object only for legacy references.
        self.projContourButton.Hide()
        self.leftToolButton.Bind(wx.EVT_BUTTON,lambda evt:self._show_tool_window(self.leftFrame))
        self.rightToolButton.Bind(wx.EVT_BUTTON,lambda evt:self._show_tool_window(self.rightFrame))
        self.projContourButton.Bind(wx.EVT_BUTTON,lambda evt:self._show_tool_window(self.projContourFrame))
        # UX: the toolbar Peaks button opens the authoritative Full Peak List.
        # Slice2D no longer owns peak-list editing; the Full Peak List is authoritative.
        self.peaksToolButton.Bind(wx.EVT_BUTTON, self.on_full_peak_list)

    def full_peak_tools_box(self):
        """Modeless editor for the main-frame Full 3D peak list.

        Deliberately independent of every legacy Slice2D-local peak structure.
        """
        self.fullToolsFrame, panel = self._make_modeless_window('Full 3D Peak Tools')
        self.fullUndoButton = wx.Button(panel, label='Undo', size=(-1, 24))
        self.fullRedoButton = wx.Button(panel, label='Redo', size=(-1, 24))
        self.fullSelectButton = PersistentStateButton(panel, label='Select', size=(-1, 24))
        self.fullMoveButton = PersistentStateButton(panel, label='Move', size=(-1, 24))
        self.fullRemoveButton = wx.Button(panel, label='Remove', size=(-1, 24))
        self.fullMaxButton = wx.Button(panel, label='Maximise', size=(-1, 24))
        self.fullAddButton = PersistentStateButton(panel, label='Add', size=(-1, 24))
        self.fullShowButton = PersistentStateButton(panel, label='Show', size=(-1, 24))
        self.fullUndoButton.Bind(wx.EVT_BUTTON, self.on_full_undo)
        self.fullRedoButton.Bind(wx.EVT_BUTTON, self.on_full_redo)
        self.fullSelectButton.Bind(wx.EVT_BUTTON, self.on_full_select)
        self.fullMoveButton.Bind(wx.EVT_BUTTON, self.on_full_move)
        self.fullRemoveButton.Bind(wx.EVT_BUTTON, self.on_full_remove)
        self.fullMaxButton.Bind(wx.EVT_BUTTON, self.on_full_maximise)
        self.fullAddButton.Bind(wx.EVT_BUTTON, self.on_full_add)
        self.fullShowButton.Bind(wx.EVT_BUTTON, self.on_full_show)
        close = wx.Button(panel, label='Close', size=(-1, 24))
        close.Bind(wx.EVT_BUTTON, self.on_full_tools_close_button)
        # Full Tools is special: closing it must also release the persistent
        # Tools button.  Replace the generic modeless-window close handler.
        self.fullToolsFrame.Unbind(wx.EVT_CLOSE)
        self.fullToolsFrame.Bind(wx.EVT_CLOSE, self.on_full_tools_close)
        root = wx.BoxSizer(wx.VERTICAL)
        # A narrow vertical tool palette leaves the Slice2D plot unobscured.
        # Give every command the same width for a clean, predictable layout.
        buttons = (self.fullUndoButton, self.fullRedoButton,
                   self.fullSelectButton, self.fullMoveButton,
                   self.fullRemoveButton, self.fullMaxButton,
                   self.fullAddButton, self.fullShowButton, close)
        for widget in buttons:
            widget.SetMinSize((92, 24))
            root.Add(widget, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 4)
        root.AddSpacer(4)
        panel.SetSizer(root); self.fullToolsFrame.SetClientSize(root.CalcMin())
        self._update_full_tool_controls()

    def _full_payload(self):
        return self.slice_service.full_peak_payload()

    def _full_records(self):
        return list(self._full_payload().get('records') or [])

    def _current_reference_name(self):
        i = self.ComboBox1.GetSelection()
        return str(self.peak[i].name) if 0 <= i < len(self.peak) else ''

    def _record_belongs_to_reference(self, record, reference_name=None):
        ref = self._current_reference_name() if reference_name is None else str(reference_name)
        return re.match(r'^%s_(\d+)$' % re.escape(ref), str(record.get('name', ''))) is not None

    def _full_record(self, name=None):
        wanted = self.full_selected_name if name is None else name
        for record in self._full_records():
            if str(record.get('name', '')) == str(wanted):
                return record
        return None

    def _set_full_status(self, text):
        # Full-peak tool instructions/results belong in the application's main
        # status bar, not in the small modeless tool palette.
        frame = wx.GetTopLevelParent(self)
        statusbar = getattr(frame, 'statusbar', None) if frame is not None else None
        if statusbar is not None:
            statusbar.SetStatusText(str(text))

    def _toolbar_tools(self, active):
        if active:
            self._show_tool_window(self.fullToolsFrame)
            self.toolbar.set_tools_active(True)
        else:
            self._hide_full_tools()

    def on_full_tools_toggle(self, event=None):
        if self.fullToolsFrame.IsShown():
            self._hide_full_tools()
        else:
            self._show_tool_window(self.fullToolsFrame)
            if hasattr(self, 'toolbar'):
                self.toolbar.set_tools_active(True)

    def _hide_full_tools(self):
        self.fullToolsFrame.Hide()
        if hasattr(self, 'toolbar'):
            self.toolbar.set_tools_active(False)

    def on_full_tools_close_button(self, event=None):
        self._hide_full_tools()

    def on_full_tools_close(self, event):
        self._hide_full_tools()
        event.Veto()

    def _update_full_tool_controls(self):
        """Derive Full Tools button state from selection and mouse mode."""
        if not hasattr(self, 'fullSelectButton'):
            return
        has_selection = self._full_record() is not None
        modal = self.full_tool_mode not in (None, 'select')

        # As in peakFrame, a completed single selection leaves Select blue.
        self.fullSelectButton.SetActive(self.full_tool_mode == 'select' or
                                        (has_selection and self.full_tool_mode is None))
        self.fullMoveButton.SetActive(self.full_tool_mode == 'move')
        self.fullAddButton.SetActive(self.full_tool_mode == 'add')
        self.fullShowButton.SetActive(self.full_tool_mode == 'show')

        self.fullSelectButton.Enable(not modal)
        self.fullMoveButton.Enable(((not modal) or self.full_tool_mode == 'move') and has_selection)
        self.fullAddButton.Enable(((not modal) or self.full_tool_mode == 'add') and not has_selection)
        self.fullShowButton.Enable(((not modal) or self.full_tool_mode == 'show') and not has_selection)
        self.fullRemoveButton.Enable((not modal) and has_selection)
        self.fullMaxButton.Enable((not modal) and has_selection)
        for button in (self.fullSelectButton, self.fullMoveButton, self.fullAddButton, self.fullShowButton):
            button.Refresh()
        self._update_full_history_buttons()

    def _set_full_tool_mode(self, mode=None):
        if mode == 'move' and self._full_record() is None:
            self._set_full_status('Move: select one Full peak first')
            self._update_full_tool_controls()
            return
        if mode in ('add', 'show') and self._full_record() is not None:
            self._set_full_status('%s: deselect the current Full peak first' % mode.capitalize())
            self._update_full_tool_controls()
            return
        self.full_tool_mode = mode
        self._update_full_tool_controls()
        if mode:
            self._set_full_status('%s: click the upper 2D slice' % mode.capitalize())

    def _clear_full_selection(self, redraw=True):
        self.full_selected_name = None
        self.full_selection_pane = None
        try:
            self.slice_service.clear_peak_selection(redraw_full3d=True)
        except Exception:
            pass
        self._update_full_tool_controls()
        if redraw:
            self.draw_figure(redraw_projections=False)

    def on_full_select(self, event=None):
        active = bool(getattr(self.fullSelectButton, 'IsActive', lambda: False)())
        if active:
            self.full_tool_mode = None
            self._clear_full_selection(redraw=True)
            self._set_full_status('Select: off')
        else:
            self._set_full_tool_mode('select')
            self._set_full_status('Select: click the upper 2D slice to select the nearest Full peak')

    def on_full_move(self, event=None):
        active = bool(getattr(self.fullMoveButton, 'IsActive', lambda: False)())
        if active:
            self._set_full_tool_mode(None)
            self._set_full_status('Move: off')
        else:
            self._set_full_tool_mode('move')

    def on_full_add(self, event=None):
        active = bool(getattr(self.fullAddButton, 'IsActive', lambda: False)())
        if active:
            self._set_full_tool_mode(None)
            self._set_full_status('Add: off')
        else:
            self._set_full_tool_mode('add')

    def on_full_show(self, event=None):
        active = bool(getattr(self.fullShowButton, 'IsActive', lambda: False)())
        if active:
            self._set_full_tool_mode(None)
            self._set_full_status('Show: off')
        else:
            self._set_full_tool_mode('show')
            self._set_full_status('Show: click a contextual Full peak in the upper 2D slice')

    def _full_snapshot(self):
        payload = self._full_payload()
        return copy.deepcopy({'records': payload.get('records') or [], 'rows': payload.get('rows') or []})

    def _push_full_undo(self):
        self.full_undo_stack.append(self._full_snapshot()); self.full_redo_stack = []
        self._update_full_history_buttons()

    def _update_full_history_buttons(self):
        if hasattr(self, 'fullUndoButton'):
            self.fullUndoButton.Enable(bool(self.full_undo_stack))
            self.fullRedoButton.Enable(bool(self.full_redo_stack))

    def _normalise_full_records(self, records):
        dim = self.slice_service.dimension
        rows = []
        for i, record in enumerate(records):
            record['row_index'] = i
            coords = tuple(float(v) for v in record.get('coordinates', ()))
            fields = list(record.get('fields') or [])
            need = dim + 2
            while len(fields) < need: fields.append('0')
            fields[0] = str(record.get('name', fields[0] if fields else ''))
            for j, value in enumerate(coords[:dim]): fields[j + 1] = str(value)
            intensity = record.get('intensity')
            if intensity is not None: fields[dim + 1] = str(float(intensity))
            record['fields'] = fields
            rows.append(list(fields))
        return rows

    def _commit_full_records(self, records):
        records = copy.deepcopy(list(records))
        rows = self._normalise_full_records(records)
        old = self._full_payload()
        self.slice_service.save_full_peak_records(records, rows, source_path=old.get('source_path'))
        # Keep the main-window selection metadata authoritative too.  This
        # routes selection to every open cross-view, including Full 3D.
        if self.full_selected_name and self._full_record(self.full_selected_name) is not None:
            try:
                self.slice_service.select_full_peak(self.full_selected_name, source_view=self, source_pane=self.full_selection_pane)
            except Exception:
                pass
        else:
            try:
                self.slice_service.clear_peak_selection(redraw_full3d=True)
            except Exception:
                pass
        self.draw_figure(redraw_projections=False)
        # Do not create Full 3D merely because a peak changed, but redraw it if
        # it is already open.  _full3d_viewer also respects PageExists().
        self.slice_service.redraw_full3d_if_open()
        self._update_full_tool_controls()

    def _restore_full_snapshot(self, snapshot):
        self._commit_full_records(snapshot.get('records') or [])

    def on_full_undo(self, event=None):
        if not self.full_undo_stack: return
        self.full_redo_stack.append(self._full_snapshot())
        snap = self.full_undo_stack.pop(); self._restore_full_snapshot(snap)
        if self._full_record() is None: self.full_selected_name = None
        self._update_full_tool_controls(); self._set_full_status('Undo')

    def on_full_redo(self, event=None):
        if not self.full_redo_stack: return
        self.full_undo_stack.append(self._full_snapshot())
        snap = self.full_redo_stack.pop(); self._restore_full_snapshot(snap)
        if self._full_record() is None: self.full_selected_name = None
        self._update_full_tool_controls(); self._set_full_status('Redo')

    def _current_full_plane_context(self):
        """Return display/fixed labels for the upper Slice2D plane."""
        if self.slice_service.dimension != 3 or len(self.slice_service.labels) < 3:
            return None
        sele = self.ComboBox1.GetSelection()
        if sele < 0: return None
        # ReSlice2d(orth=0): data[:, fixed-dim1, dim2-window].T is plotted
        # with X=labb[0], Y=labb[2].  labb[1] is the fixed plane coordinate.
        fixed_index = int(self.slice_service.peak_indices[sele][0]) + int(self.inc)
        scale = numpy.asarray(self.slice_service.axis(1))
        fixed_index = max(0, min(fixed_index, len(scale) - 1))
        labels = self.slice_service.labels
        return {'x_label': str(labels[0]), 'y_label': str(labels[2]),
                'fixed_label': str(labels[1]), 'fixed_value': float(scale[fixed_index])}

    def _update_full_record_axes(self, record, x, y, fixed_value=None):
        ctx = self._current_full_plane_context()
        if ctx is None: return
        axes = dict(record.get('axis_values') or {})
        axes[ctx['x_label']] = float(x); axes[ctx['y_label']] = float(y)
        if fixed_value is not None: axes[ctx['fixed_label']] = float(fixed_value)
        record['axis_values'] = axes
        labels = list(self.slice_service.labels); dim = 3
        record['coordinates'] = tuple(float(axes[str(labels[dim - 1 - i])]) for i in range(dim))

    def _full_candidates(self):
        return [r for r in self._full_records() if self._record_belongs_to_reference(r)]

    def _select_full_at(self, x, y, pane='top'):
        ctx = self._full_plane_context_for_pane(pane)
        ref = self._full_reference_for_pane(pane)
        candidates = [r for r in self._full_records() if self._record_belongs_to_reference(r, ref)]
        if ctx is None or not candidates: self._set_full_status('No Full peaks for this reference'); return
        axes = self.axes1 if pane == 'top' else self.axes2
        xr = abs(axes.get_xlim()[1] - axes.get_xlim()[0]) or 1.0
        yr = abs(axes.get_ylim()[1] - axes.get_ylim()[0]) or 1.0
        usable = [r for r in candidates if ctx['x_label'] in r.get('axis_values', {}) and ctx['y_label'] in r.get('axis_values', {})]
        if not usable: return
        record = min(usable, key=lambda r: ((float(r['axis_values'][ctx['x_label']])-x)/xr)**2 + ((float(r['axis_values'][ctx['y_label']])-y)/yr)**2)
        self.full_selected_name = str(record.get('name'))
        self.full_selection_pane = pane
        try: self.slice_service.select_full_peak(self.full_selected_name, source_view=self, source_pane=self.full_selection_pane)
        except Exception: pass
        # Single Select is intentionally persistent, matching peakFrame: each
        # subsequent click replaces the selection until Select is pressed off.
        self.full_tool_mode = 'select'
        self._set_full_status('Selected: %s' % self.full_selected_name)
        self._update_full_tool_controls()
        self.draw_figure(redraw_projections=False)

    def _move_full_at(self, x, y):
        record = self._full_record()
        if record is None: return
        self._push_full_undo(); records = self._full_records()
        for r in records:
            if str(r.get('name')) == self.full_selected_name:
                self._update_full_record_axes(r, x, y); break
        self.full_tool_mode = None
        self._commit_full_records(records)
        self._set_full_status('Moved: %s' % self.full_selected_name)

    def _next_full_name(self, reference_name):
        rx = re.compile(r'^%s_(\d+)$' % re.escape(str(reference_name)))
        nums = [int(m.group(1)) for r in self._full_records() for m in [rx.match(str(r.get('name','')))] if m]
        return '%s_%d' % (reference_name, (max(nums) + 1) if nums else 1)

    def _add_full_at(self, x, y):
        ctx = self._current_full_plane_context(); ref = self._current_reference_name()
        if ctx is None or not ref: return
        self._push_full_undo(); records = self._full_records(); name = self._next_full_name(ref)
        axes = {ctx['x_label']: float(x), ctx['y_label']: float(y), ctx['fixed_label']: float(ctx['fixed_value'])}
        labels = list(self.slice_service.labels)
        coords = tuple(float(axes[str(labels[2-i])]) for i in range(3))
        intensity = self._full_intensity_from_axes(axes)
        fields = [name] + [str(v) for v in coords] + [str(float(intensity))]
        records.append({'name': name, 'coordinates': coords, 'axis_values': axes, 'intensity': float(intensity),
                        'row_index': len(records), 'fields': fields, 'analysis': {}})
        self.full_selected_name = name
        self.full_tool_mode = None
        self._commit_full_records(records)
        self._set_full_status('Added: %s' % name)

    def _full_intensity_from_axes(self, axes):
        indices = []
        for i, label in enumerate(self.slice_service.labels[:3]):
            scale = numpy.asarray(self.slice_service.axis(i))
            indices.append(int(numpy.abs(scale - float(axes[str(label)])).argmin()))
        try: return self.slice_service.sample(indices)
        except Exception: return 0.0

    def on_full_remove(self, event=None):
        if self._full_record() is None: self._set_full_status('Remove: select a Full peak first'); return
        name = self.full_selected_name; self._push_full_undo()
        # Clear the local selection before committing so the authoritative
        # selection metadata and all viewers are cleared in the same commit.
        self.full_selected_name = None
        self.full_selection_pane = None
        self.full_tool_mode = None
        self._commit_full_records([r for r in self._full_records() if str(r.get('name')) != name])
        self._set_full_status('Removed: %s' % name)

    def on_full_maximise(self, event=None):
        record = self._full_record(); plane = self._full_plane
        if record is None: self._set_full_status('Maximise: select a Full peak first'); return
        if not plane: self._set_full_status('Maximise: redraw the 2D slice first'); return
        ctx = self._current_full_plane_context(); axes = record.get('axis_values', {})
        if ctx is None or ctx['x_label'] not in axes or ctx['y_label'] not in axes: return
        X, Y, Z = plane
        xvec = numpy.asarray(X[0, :] if numpy.asarray(X).ndim == 2 else X)
        yvec = numpy.asarray(Y[:, 0] if numpy.asarray(Y).ndim == 2 else Y)
        # Some mesh orientations vary on the opposite dimensions; choose the vectors that actually vary.
        if numpy.asarray(X).ndim == 2 and numpy.ptp(xvec) == 0: xvec = numpy.asarray(X[:, 0])
        if numpy.asarray(Y).ndim == 2 and numpy.ptp(yvec) == 0: yvec = numpy.asarray(Y[0, :])
        arr = numpy.asarray(Z)
        # Robustly locate the closest mesh cell directly; then hill-climb abs intensity in this 2D plane only.
        dist = (numpy.asarray(X)-float(axes[ctx['x_label']]))**2 + (numpy.asarray(Y)-float(axes[ctx['y_label']]))**2
        iy, ix = numpy.unravel_index(int(numpy.nanargmin(dist)), dist.shape)
        for _ in range(30):
            y0,y1=max(0,iy-1),min(arr.shape[0],iy+2); x0,x1=max(0,ix-1),min(arr.shape[1],ix+2)
            sub=numpy.abs(arr[y0:y1,x0:x1]); dy,dx=numpy.unravel_index(int(numpy.nanargmax(sub)),sub.shape)
            ni,nj=y0+dy,x0+dx
            if (ni,nj)==(iy,ix): break
            iy,ix=ni,nj
        newx=float(numpy.asarray(X)[iy,ix]); newy=float(numpy.asarray(Y)[iy,ix])
        self._push_full_undo(); records=self._full_records()
        for r in records:
            if str(r.get('name')) == self.full_selected_name:
                self._update_full_record_axes(r,newx,newy); r['intensity']=float(arr[iy,ix]); break
        self._commit_full_records(records); self._set_full_status('Maximised in current 2D plane: %s' % self.full_selected_name)

    def _handle_full_tool_click(self, event):
        if self.full_tool_mode is None or event.xdata is None or event.ydata is None:
            return False
        if event.inaxes is self.axes1:
            pane = 'top'
        elif event.inaxes is self.axes2 and not self.cb_1d.GetValue():
            pane = 'bottom'
        else:
            return False
        if self.full_tool_mode == 'select':
            self._select_full_at(float(event.xdata), float(event.ydata), pane=pane)
        elif self.full_tool_mode == 'move':
            # Editing remains tied to the upper plane for now; selection itself
            # is pane-aware.  This preserves the existing coordinate semantics.
            if pane != 'top': return False
            self._move_full_at(float(event.xdata), float(event.ydata))
        elif self.full_tool_mode == 'add':
            if pane != 'top': return False
            self._add_full_at(float(event.xdata), float(event.ydata))
        elif self.full_tool_mode == 'show':
            if pane != 'top': return False
            self._show_neighbour_at(float(event.xdata), float(event.ydata))
        return True

    def _lower_is_dependent(self):
        return bool(self.cb_decon.GetValue() or self.cb_flip.GetValue() or self.cb_1d.GetValue())

    def _full_reference_for_pane(self, pane):
        if pane == 'top' or self._lower_is_dependent():
            idx = self.ComboBox1.GetSelection()
        else:
            idx = self.ComboBox2.GetSelection()
        return str(self.peak[idx].name) if 0 <= idx < len(self.peak) else ''

    def _full_plane_context_for_pane(self, pane):
        if self.slice_service.dimension != 3 or len(self.slice_service.labels) < 3:
            return None
        labels = [str(v) for v in self.slice_service.labels[:3]]
        if pane == 'top':
            return {'x_label': labels[0], 'y_label': labels[2]}
        if self.cb_1d.GetValue():
            return None
        return {'x_label': labels[0], 'y_label': labels[1] if self.cb_flip.GetValue() else labels[2]}

    def _selection_visible_in_pane(self, pane):
        if not self.full_selected_name:
            return False
        # Derived lower views (Decon/Orth/1D) deliberately mirror the upper
        # slice.  In ordinary 2D/2D mode each pane owns its selection styling.
        if self._lower_is_dependent():
            return True
        return self.full_selection_pane == pane

    def _full_peak_plot_label(self, record_or_name):
        """Return the within-slice number used to annotate a Full peak.

        Canonical names remain the authoritative identifiers (for example
        ``180_2``); only plot presentation is shortened to the generated
        numeric suffix (``2``).
        """
        if isinstance(record_or_name, dict):
            name = str(record_or_name.get('name', ''))
        else:
            name = str(record_or_name or '')
        match = re.search(r'_(\d+)$', name)
        return match.group(1) if match else name

    def _full_slice_context_for_pane(self, pane, orth=False):
        """Return the fixed 3D axis and current discrete slice for a pane.

        This mirrors ReSlice2d exactly: normal planes fix dimension 1, while
        orthogonal planes fix dimension 2.  Keeping the calculation in index
        space makes Slice2D use the same notion of slice distance as Full3D.
        """
        if self.slice_service.dimension != 3 or len(self.slice_service.labels) < 3:
            return None
        labels = [str(v) for v in self.slice_service.labels[:3]]
        if pane == 'top':
            peak_index = self.ComboBox1.GetSelection()
            increment = int(self.inc)
        elif orth:
            # ReSlice2d lower orthogonal view always uses the left reference;
            # a deconvolved dependent view follows inc, raw orthogonal uses inc2.
            peak_index = self.ComboBox1.GetSelection()
            increment = int(self.inc if self.cb_decon.GetValue() else self.inc2)
        elif self.cb_decon.GetValue():
            peak_index = self.ComboBox1.GetSelection()
            increment = int(self.inc)
        else:
            peak_index = self.ComboBox2.GetSelection()
            increment = int(self.inc2)
        if peak_index < 0 or peak_index >= len(self.slice_service.peak_indices):
            return None
        if orth:
            scale = numpy.asarray(self.slice_service.axis(2))
            current = int(self.slice_service.peak_indices[peak_index][1]) + increment
            slice_label = labels[2]
        else:
            scale = numpy.asarray(self.slice_service.axis(1))
            current = int(self.slice_service.peak_indices[peak_index][0]) + increment
            slice_label = labels[1]
        if len(scale) == 0:
            return None
        current = max(0, min(current, len(scale) - 1))
        return {'slice_label': slice_label, 'slice_scale': scale, 'slice_index': current}

    def _focused_reference_for_pane(self, pane='top'):
        """Return the canonical reference currently focused by a Slice2D pane."""
        if pane == 'top' or self._lower_is_dependent():
            combo = self.ComboBox1
        else:
            combo = self.ComboBox2
        index = combo.GetSelection()
        return str(combo.GetString(index)) if 0 <= index < combo.GetCount() else ''

    def _displayed_full_peak_candidates(self, pane='top', orth=False, contextual_only=False):
        """Return Full records actually eligible for the Slice2D overlay.

        The returned delta is measured in discrete planes, exactly as used for
        ornament greyscale.  The boolean context flag distinguishes peaks belonging
        to the pane's focused canonical reference from other displayed peaks.  Show
        uses this same helper so hit-testing can never navigate to an ornament that
        is not represented by the current display rules.
        """
        if self.slice_service.dimension != 3 or len(self.slice_service.labels) < 3:
            return []
        labels = [str(v) for v in self.slice_service.labels[:3]]
        ctx = {'x_label': labels[0], 'y_label': labels[1] if orth else labels[2]}
        slice_ctx = self._full_slice_context_for_pane(pane, orth=orth)
        if slice_ctx is None:
            return []
        slice_label = slice_ctx['slice_label']
        scale = slice_ctx['slice_scale']
        current = slice_ctx['slice_index']
        result = []
        focused_reference = self._focused_reference_for_pane(pane)
        for record in self._full_records():
            values = record.get('axis_values', {})
            if (ctx['x_label'] not in values or ctx['y_label'] not in values or
                    slice_label not in values):
                continue
            peak_slice_index = int(numpy.argmin(numpy.abs(scale - float(values[slice_label]))))
            delta = peak_slice_index - current
            if abs(delta) > 2:
                continue
            is_contextual = not self._record_belongs_to_reference(record, focused_reference)
            if contextual_only and not is_contextual:
                continue
            result.append((record, delta, ctx, is_contextual))
        return result

    def _show_neighbour_at(self, x, y):
        """Focus the Right pane on the parent slice of the nearest contextual peak."""
        # The upper contour is the comparison anchor and is never changed.
        orth = bool(self.cb_flip.GetValue() and self.cb_decon.GetValue())
        candidates = self._displayed_full_peak_candidates('top', orth=orth, contextual_only=True)
        if not candidates:
            return
        xr = abs(self.axes1.get_xlim()[1] - self.axes1.get_xlim()[0]) or 1.0
        yr = abs(self.axes1.get_ylim()[1] - self.axes1.get_ylim()[0]) or 1.0
        record, delta, ctx, is_contextual = min(
            candidates,
            key=lambda item: ((float(item[0]['axis_values'][item[2]['x_label']]) - x) / xr) ** 2 +
                             ((float(item[0]['axis_values'][item[2]['y_label']]) - y) / yr) ** 2)
        canonical = str(record.get('name', ''))
        parent_name = canonical.rsplit('_', 1)[0] if '_' in canonical else ''
        choices = [str(self.ComboBox2.GetString(i)) for i in range(self.ComboBox2.GetCount())]
        try:
            target = choices.index(parent_name)
        except ValueError:
            self._set_full_status('Show: parent slice %s is not available' % parent_name)
            return
        # Show is a comparison/navigation command.  Dependent lower-pane modes
        # cannot display an independently selected Right reference, so return to
        # the normal 2D/2D interpretation without changing the Left selection.
        self.cb_decon.SetValue(False)
        self.cb_flip.SetValue(False)
        self.cb_1d.SetValue(False)
        self._set_combo_selection(self.ComboBox2, target)
        self.inc2 = 0
        self._sync_slice_spin_controls()
        self.ax_reset2 = 1
        self.selection = []
        self.draw_figure(redraw_projections=False)
        self._set_full_status('Show: %s in Right pane (parent slice %s)' % (canonical, parent_name))

    def _draw_full_peak_overlay(self, target_axes=None, reference_name=None, orth=False, show_labels=None, pane='top'):
        """Draw nearby Full peaks using the same slice-distance cue as Full3D.

        Peak membership is intentionally geometric rather than reference-name
        based: every authoritative Full peak within two discrete planes of the
        displayed slice is eligible.  Editing/selection membership remains
        unchanged elsewhere.
        """
        target_axes = self.axes1 if target_axes is None else target_axes
        if self.slice_service.dimension != 3 or len(self.slice_service.labels) < 3:
            return
        labels = [str(v) for v in self.slice_service.labels[:3]]
        ctx = {'x_label': labels[0], 'y_label': labels[1] if orth else labels[2]}
        slice_ctx = self._full_slice_context_for_pane(pane, orth=orth)
        if slice_ctx is None:
            return
        # Full3D only draws ornaments when its peak-overlay control is active.
        # Slice2D's existing call sites already decide when ornaments are drawn;
        # once drawn, labels accompany every nearby marker just as in Full3D.
        for r, delta, ctx, is_contextual in self._displayed_full_peak_candidates(pane, orth=orth):
            values = r.get('axis_values', {})
            canonical_name = str(r.get('name', ''))
            selected = (canonical_name == str(self.full_selected_name) and
                        self._selection_visible_in_pane(pane))
            color = ('r' if selected else
                     (self.slice_service.full_peak_slice_color(delta) or
                      str(min(0.15 + 0.25 * abs(delta), 0.85))))
            x = float(values[ctx['x_label']]); y = float(values[ctx['y_label']])
            target_axes.scatter(x, y, c=color, s=55, marker='x', zorder=5)
            # Focus membership and geometric slice distance are deliberately
            # independent.  A peak from another canonical reference may lie in
            # exactly the same 2D plane (delta == 0): it stays black, but is
            # presented as contextual using the full canonical identifier.
            # Contextual labels are smaller and tucked closer to their ornament.
            if not is_contextual:
                plot_label = self._full_peak_plot_label(r)
                label_fontsize = 11
                label_offset = (5, 5)
            else:
                plot_label = canonical_name
                label_fontsize = 8
                label_offset = (3, 3)
            target_axes.annotate(plot_label, (x, y),
                                 xytext=label_offset, textcoords='offset points',
                                 ha='left', va='bottom', rotation=0,
                                 fontsize=label_fontsize, color=color, zorder=6)

    def _draw_pane_peak_name(self, target_axes, peak_index):
        """Display the currently viewed reference peak in the plot itself."""
        if peak_index < 0 or peak_index >= len(self.slice_service.peaks):
            return
        name = str(getattr(self.slice_service.peak(peak_index), 'name', ''))
        target_axes.text(0.02, 0.98, name, transform=target_axes.transAxes,
                         ha='left', va='top', fontsize=10, fontweight='bold',
                         zorder=20)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.


        self.fig=Figure()
        self.fig.clear()
        from matplotlib.gridspec import GridSpec
        # Use a compact 2 x 2 layout.  The old 2 x 9 GridSpec left a whole
        # unused column between the slices and the projections and retained
        # Matplotlib's generous default outer margins.
        gs1 = GridSpec(2, 2, width_ratios=[5.0, 1.25],
                       left=0.055, right=0.985, bottom=0.065, top=0.975,
                       wspace=0.10, hspace=0.10)
        self.axesCC = self.fig.add_subplot(gs1[0, 1])
        self.axesHC = self.fig.add_subplot(gs1[1, 1])
        self.axes1 = self.fig.add_subplot(gs1[0, 0])
        self.axes2 = self.fig.add_subplot(gs1[1, 0], sharex=self.axes1)

        self.axes1.xaxis.set_visible(False)

        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(100,100))

        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        # Projection contours are static between explicit Draw! requests.
        # draw_event is also emitted after toolbar zoom/pan, so recapturing the
        # backgrounds here keeps projection-line blitting compatible with the
        # standard Matplotlib navigation toolbar.
        self.canvas.mpl_connect('draw_event', self._on_projection_draw_event)
        self._projection_backgrounds = {}
        self._projection_artists = []
        self._projection_static_ready = False
        self._projection_blit_busy = False

        # Build modeless settings windows first so toolbar buttons can reference them.
        self.contour_boxes()
        self.left_box()
        self.right_box()
        self.full_peak_tools_box()
        self.control_box()

        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view, peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon, orth_callback=self._toolbar_orth, one_d_callback=self._toolbar_1d, contour_callback=self._toolbar_contours, tools_callback=self._toolbar_tools, peaks_active=True, coordinates=False)
        # Keep all application widgets inside the native Matplotlib wx.ToolBar.
        # The toolbar's own wx.TB_BOTTOM border then runs continuously over the
        # complete row instead of ending after Home/Pan/Zoom/Save.
        # Keep Matplotlib's native tools together, followed by our controls and
        # finally the live coordinate readout.
        self.toolbar.AddSeparator()
        # Keep the application controls visually grouped: drawing/display
        # options, focused Left/Right reference navigation, then projection
        # and peak tools.
        for widget in (self.leftToolButton, self.quickSliceSpin1, self.rightToolButton, self.quickSliceSpin2):
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)
        self.toolbar.AddSeparator()
        for widget in (self.peaksToolButton,):
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)
        self.toolbar.bind_control_status_help(self.leftToolButton, 'Open left reference controls')
        self.toolbar.bind_control_status_help(self.quickSliceSpin1, 'Previous or next left reference peak')
        self.toolbar.bind_control_status_help(self.rightToolButton, 'Open right reference controls')
        self.toolbar.bind_control_status_help(self.quickSliceSpin2, 'Previous or next right reference peak')
        self.toolbar.bind_control_status_help(self.peaksToolButton, 'Open peak controls')
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 2)
        self.SetSizer(self.vbox)
        self.Layout()

        logging.info('Checking projections in shared data store...')
        if self._projection_view_xy() is None or self._projection_view_yz() is None:
            raise RuntimeError('2D Slices requires raw projection views in the shared data store')
        logging.info('Done')


    def _projection_view_xy(self):
        labels = self.slice_service.labels
        return self.slice_service.projection(labels[1], labels[0], decon=False, transpose='n')

    def _projection_view_yz(self):
        labels = self.slice_service.labels
        return self.slice_service.projection(labels[2], labels[1], decon=False, transpose='n')

    @property
    def Xs_xy(self): return self._projection_view_xy()['XX']
    @property
    def Ys_xy(self): return self._projection_view_xy()['YY']
    @property
    def Zs_xy(self): return self._projection_view_xy()['ZZ']
    @property
    def Xs_yz(self): return self._projection_view_yz()['XX']
    @property
    def Ys_yz(self): return self._projection_view_yz()['YY']
    @property
    def Zs_yz(self): return self._projection_view_yz()['ZZ']
    @property
    def Xs_xy_xmin(self): return numpy.min(self.Xs_xy)
    @property
    def Xs_xy_xmax(self): return numpy.max(self.Xs_xy)
    @property
    def Xs_xy_ymin(self): return numpy.min(self.Ys_xy)
    @property
    def Xs_xy_ymax(self): return numpy.max(self.Ys_xy)
    @property
    def Xs_yz_xmin(self): return numpy.min(self.Xs_yz)
    @property
    def Xs_yz_xmax(self): return numpy.max(self.Xs_yz)
    @property
    def Xs_yz_ymin(self): return numpy.min(self.Ys_yz)
    @property
    def Xs_yz_ymax(self): return numpy.max(self.Ys_yz)

    #make an index
    def index(self,array):
        index=[]
        for i in range(len(array)):
            index.append((array[i].name))
        return index

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=string.split(line)
            peak.append(linetosave)
        peakfile.close()
        return peak

    def findnear_index(self,test,array):
        #array = numpy.asarray(array)
        idx = (numpy.abs(array - test)).argmin()
        return idx
        #return array[idx]


    def GetLevels(self,max_level,min_level,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels


    #get 2d strips from 3d data
    def GetSlice2d(self):
        # dic=self.tabOne.dic
        # self.data=self.slice_service.data
        # self.index0=self.slice_service.axis(0)
        # self.index1=self.slice_service.axis(1)
        # self.index2=self.slice_service.axis(2)

        self.DECON=0
        if self.slice_service.deconvolution_enabled:
            #self.datadec=self.slice_service.datadec
            if(self.slice_service.datadec.shape==self.slice_service.data.shape):
                logging.info('Shape of deconvolved, and raw are different')
                self.DECON=0
            else:
                self.DECON=1


    #get 2d strips from 3d data
    def ReSlice2d(self,arr,inc,pkl,peak,width,orth=0):

        if(orth==0):
            #print out 2D slice for each peak correlation
            #print "Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width
            ptC=self.slice_service.peak_indices[pkl][0]
            ptC=ptC+inc
            ptH=self.slice_service.peak_indices[pkl][1]

            ptH_max=self.findnear_index(float(peak[pkl].x)+float(width)/2,self.slice_service.axis(2))#find the nearest point to desired chemical shift in carbon index
            ptH_min=self.findnear_index(float(peak[pkl].x)-float(width)/2,self.slice_service.axis(2))#find the nearest point to desired chemical shift in carbon index
            Xs,Ys=self.slice_service.slice_meshes(orth=False, pt_c=ptC, pt_h_max=ptH_max, pt_h_min=ptH_min)
            Zs=arr[:,ptC,ptH_max:ptH_min].transpose() #extract the relevant 2d slice
        else:
            #print out 2D slice for each peak correlation
            #print "Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width

            #ptC=self.findnear_index(float(peak[pkl][2]),self.index1)#find the nearest point to desired chemical shift in carbon index
            ptC_max=self.findnear_index(float(peak[pkl].y)+float(width)/2.,self.slice_service.axis(1))#find the nearest point to desired chemical shift in carbon index
            ptC_min=self.findnear_index(float(peak[pkl].y)-float(width)/2.,self.slice_service.axis(1))#find the nearest point to desired chemical shift in carbon index
            ptH=self.slice_service.peak_indices[pkl][1]
            #ptH=self.findnear_index(float(peak[pkl][2]),self.index2)#find the nearest point to desired chemical shift in carbon index
            ptH=ptH+inc
            Xs,Ys=self.slice_service.slice_meshes(orth=True, pt_h=ptH, pt_c_max=ptC_max, pt_c_min=ptC_min)
            Zs=arr[:,ptC_max:ptC_min,ptH].transpose() #extract the relevant 2d slice

        #print 'Done!'
        return Xs,Ys,Zs




    def _draw_1d_lower_panel(self, peak_index, show_decon=False, show_labels=False):
        """Draw the centre 1D cut for the selected left/reference peak.

        axes2 intentionally retains sharex=axes1: the trace is the dimension-0
        cut through the centre of the 2D slice shown immediately above it.
        """
        view = self.slice_service.reference_1d_view(peak_index)
        if not view:
            self.axes2.text(0.5, 0.5, '1D slice unavailable',
                            transform=self.axes2.transAxes, ha='center', va='center')
            return

        x = view['x']
        noise_sigma = view.get('noise_sigma')
        if noise_sigma is None or not numpy.isfinite(noise_sigma) or noise_sigma <= 0:
            self.axes2.text(0.5, 0.5, 'S/N unavailable: no valid noise estimate',
                            transform=self.axes2.transAxes, ha='center', va='center')
            return
        raw = numpy.asarray(view['raw']) / noise_sigma
        threshold = float(view['threshold']) / noise_sigma

        self.axes2.set_xlabel(view.get('x_label', self.slice_service.labels[0] + ' (ppm)'), fontsize=8)
        self.axes2.set_ylabel('S/N', fontsize=8)
        self.axes2.format_coord = lambda xpos, snr: 'x=%.4f, S/N=%.3f' % (xpos, snr)
        self.axes2.plot(x, raw, 'r', label='Data', lw=1.0)
        self.axes2.axhline(threshold, color='g', linestyle='--', lw=0.9, label='Noise Threshold')

        if show_decon and view.get('decon') is not None:
            self.axes2.plot(x, numpy.asarray(view['decon']) / noise_sigma, 'b', label='Deconvolved', lw=1.0)

        # Draw the picked-peak intensity stems whenever Peaks is enabled.
        # Peak markers are resolved from the authoritative Full Peak List so
        # interactive and report rendering use the same persistent source.
        markers = list(view.get('markers', []) or [])
        if show_labels and not markers:
            reference_name = str(getattr(view.get('peak'), 'name', ''))
            x_label = str(self.slice_service.labels[0]) if self.slice_service.labels else None
            for record in self._full_records():
                if not self._record_belongs_to_reference(record, reference_name):
                    continue
                values = record.get('axis_values', {})
                if x_label not in values:
                    continue
                try:
                    markers.append({
                        'x': float(values[x_label]),
                        'height': float(record.get('intensity', 0.0)),
                        'label': self._full_peak_plot_label(record),
                    })
                except (TypeError, ValueError):
                    continue
        if show_labels:
            for marker in markers:
                xpos = marker['x']
                self.axes2.plot((xpos, xpos), (0, marker['height'] / noise_sigma), 'k', lw=1.0)

        ymin, ymax = self.axes2.get_ylim()
        offset = -ymin / 2.0
        reference_name = str(getattr(view.get('peak'), 'name', ''))
        x_label = str(self.slice_service.labels[0]) if self.slice_service.labels else None
        for record in self._full_records():
            canonical_name = str(record.get('name', ''))
            selected = (canonical_name == str(self.full_selected_name) and
                        self._selection_visible_in_pane('bottom'))
            if not self._record_belongs_to_reference(record, reference_name):
                continue
            values = record.get('axis_values', {})
            if x_label not in values or not (show_labels or selected):
                continue
            xpos = float(values[x_label])
            self.axes2.text(xpos, -offset, self._full_peak_plot_label(record),
                            fontsize=8, rotation=0, ha='center', va='center',
                            color='r' if selected else 'k', zorder=6)

        self.axes2.yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: '%.1f' % value))
        self.axes2.legend(fontsize=8)

    def on_cb_1d(self, event):
        self.ax_reset2 = 1
        self.draw_figure()

    def draw_figure(self, redraw_projections=False):
        self._sync_slice_spin_controls()
        """Redraw the slice panels and update projection overlays.

        Projection contour plots are deliberately cached.  They are rebuilt
        only for the initial draw or when Draw! explicitly requests it; the
        green/blue guide lines are persistent animated artists updated by
        blitting.
        """

        if(self.ax_resetHC==0):
            x_minHC,x_maxHC=self.axesHC.get_xlim()
            y_minHC,y_maxHC=self.axesHC.get_ylim()
        if(self.ax_resetCC==0):
            x_minCC,x_maxCC=self.axesCC.get_xlim()
            y_minCC,y_maxCC=self.axesCC.get_ylim()
        if(self.ax_reset1==0):
            x_min1,x_max1=self.axes1.get_xlim()
            y_min1,y_max1=self.axes1.get_ylim()
        if(self.ax_reset2==0):
            x_min2,x_max2=self.axes2.get_xlim()
            y_min2,y_max2=self.axes2.get_ylim()

        #colormap=cm.seismic
        colormap=cm.Reds
        #colormap=cm.RdYlBu
        #from matplotlib.gridspec import GridSpec
        #gs1=GridSpec(2,9)

        #self.fig.clear()

        Width1=float(self.width1Box.GetValue())
        Width2=float(self.width2Box.GetValue())
        if(Width1==0):
            Width1=1
            self.width1Box.SetValue(str(1))
        if(Width2==0):
            Width2=1
            self.width2Box.SetValue(str(1))

        max_levelP=float(self.textbox_maxP.GetValue())
        min_levelP=float(self.textbox_minP.GetValue())
        ctr_levelP=int(self.textbox_lvlP.GetValue())
        levelsP=self.GetLevels(max_levelP,min_levelP,ctr_levelP)

        max_level=float(self.textbox1.GetValue())
        min_level=float(self.textbox0.GetValue())
        ctr_level=int(self.textbox2.GetValue())
        levels_left=self.GetLevels(max_level,min_level,ctr_level)
        max_level_right=float(self.textbox1_right.GetValue())
        min_level_right=float(self.textbox0_right.GetValue())
        ctr_level_right=int(self.textbox2_right.GetValue())
        levels_right=self.GetLevels(max_level_right,min_level_right,ctr_level_right)

        sele1=self.ComboBox1.GetSelection()
        sele2=self.ComboBox2.GetSelection()

        dimC=float(self.slice_service.peaks[sele1].ppmJ)
        dimH=float(self.slice_service.peaks[sele1].ppmK)

        label=self.slice_service.peaks[self.ComboBox1.GetSelection()].name

        dimC2=float(self.slice_service.peaks[sele2].ppmJ)
        dimH2=float(self.slice_service.peaks[sele2].ppmK)

        orth_cb=self.cb_flip.GetValue()
        dec_cb=self.cb_decon.GetValue()
        one_d_cb=self.cb_1d.GetValue()
        grid_cb=self.cb_grid_auto.GetValue()


        ###### Projection panels (cached contours + blitted guide lines) ######
        if redraw_projections or not self._projection_static_ready:
            self._draw_projection_static(levelsP, colormap, sele1, dimC, dimH, label)

        self._update_projection_guides(sele1, sele2, dimC2, dimH2)


        ###################################################################
        #Subplot 3 - the slice NUMBER TWO
        #self.axes2 = self.fig.add_subplot(gs1[1,:-2])
        self.axes2.clear()


        if one_d_cb:
            # 1D mode has priority for the lower panel.  It uses ComboBox1
            # (Left) and shares axes1's X axis by construction.
            self._enable_right_slice_spin(False)
            self._enable_right_slice_spin(False)
            self._draw_1d_lower_panel(sele1, show_decon=dec_cb, show_labels=grid_cb)
            self._draw_pane_peak_name(self.axes2, sele1)
            # Do not set X limits here: sharex keeps the lower trace locked to
            # the exact ppm window of the 2D slice above.
            self.ax_reset2 = 0
        else:
            if(orth_cb):
                self.axes2.set_xlabel(self.slice_service.labels[0]+' (ppm)',fontsize=8)
                self.axes2.set_ylabel(self.slice_service.labels[1]+' (ppm)',fontsize=8)
            else:
                self.axes2.set_xlabel(self.slice_service.labels[0]+' (ppm)',fontsize=8)
                self.axes2.set_ylabel(self.slice_service.labels[2]+' (ppm)',fontsize=8)

            if(orth_cb and dec_cb):
                self._enable_right_slice_spin(False)
                self._enable_right_slice_spin(False)
                Xs,Ys,Zs=self.ReSlice2d(self.slice_service.datadec,self.inc,sele1,self.slice_service.peaks,Width2,orth=orth_cb)
            elif(orth_cb): #orthoganol
                self._enable_right_slice_spin(True)
                self._enable_right_slice_spin(True)
                Xs,Ys,Zs=self.ReSlice2d(self.slice_service.data,self.inc2,sele1,self.slice_service.peaks,Width2,orth=orth_cb)
            elif(dec_cb): #decon
                self._enable_right_slice_spin(False)
                self._enable_right_slice_spin(False)
                Xs,Ys,Zs=self.ReSlice2d(self.slice_service.datadec,self.inc,sele1,self.slice_service.peaks,Width1,orth=orth_cb)
            else: #other combo box
                self._enable_right_slice_spin(True)
                self._enable_right_slice_spin(True)
                Xs,Ys,Zs=self.ReSlice2d(self.slice_service.data,self.inc2,sele2,self.slice_service.peaks,Width1,orth=orth_cb)

            self.axes2.contour(Xs, Ys, Zs, levels_right, cmap=colormap, norm=colors.Normalize(vmin=-numpy.max(levels_right), vmax=numpy.max(levels_right))) #plot pdb network
            lower_index = sele1 if (orth_cb or dec_cb) else sele2
            lower_name = str(self.slice_service.peaks[lower_index].name)
            self._draw_full_peak_overlay(self.axes2, lower_name, orth=orth_cb, show_labels=grid_cb, pane='bottom')
            self._draw_pane_peak_name(self.axes2, lower_index)



            #do the main peak labels
            if(orth_cb and dec_cb):
                    self.axes2.scatter(dimC,dimC,c='g',s=100)
            if(orth_cb):#if orthoganol
                    self.axes2.scatter(dimC,dimC,c='g',s=100)
            elif(dec_cb):
                    self.axes2.scatter(dimH,dimH,c='g',s=100)
            else:
                    self.axes2.scatter(dimC2,dimH2,c='g',s=100)


            y_max2a=Ys[0][0]
            y_min2a=Ys[(len(Ys))-1][0]
            x_max2a=Xs[0][0]
            x_min2a=Xs[0][(len(Xs[0]))-1]

            if(grid_cb):#horizontal line
                xl=(x_max2a,x_min2a)
                if(self.cb_flip.GetValue()): #if we want orthogonal, get value from first tick box
                    hl=(dimC,dimC)
                elif(self.cb_decon.GetValue()): #otherwise, if want to see decon
                    hl=(dimH,dimH)
                else: #else go into the other combo-box.
                    hl=(dimH2,dimH2)
                self.axes2.plot(xl,hl,'blue', lw = 0.5)

                # plot slice-position guide
                if(orth_cb): #if orthoganol
                    xl=(dimC,dimC2)
                    hl=(dimC,dimC)
                    self.axes2.plot(xl,hl,'cyan') #horizontal
                    yd=(y_min2a,dimC)
                    xd=(dimC,dimC)
                    self.axes2.plot(xd,yd,'cyan') #vertical 1
                    yd=(y_min2a,dimC)
                    xd=(dimC2,dimC2)
                    self.axes2.plot(xd,yd,'cyan') #vertical 2
                elif(dec_cb):
                    xl=(dimC,dimC2)
                    hl=(dimH,dimH)
                    self.axes2.plot(xl,hl,'cyan') #horizontal
                    yd=(y_min2a,dimH)
                    xd=(dimC2,dimC2)
                    self.axes2.plot(xd,yd,'cyan') #vertical 1
                    yd=(y_min2a,dimH2)
                    xd=(dimC,dimC)
                    self.axes2.plot(xd,yd,'cyan') #vertical 2
                else:
                    xl=(dimC,dimC2)
                    hl=(dimH2,dimH2)
                    self.axes2.plot(xl,hl,'cyan') #horizontal
                    yd=(y_max2a,dimH2)
                    xd=(dimC2,dimC2)
                    self.axes2.plot(xd,yd,'cyan') #vertical 1
                    yd=(y_max2a,dimH2)
                    xd=(dimC,dimC)
                    self.axes2.plot(xd,yd,'cyan') #vertical 2


            if(self.ax_reset2==1):
                y_max2=float(numpy.nanmax(Ys))
                y_min2=float(numpy.nanmin(Ys))
                x_max2=float(numpy.nanmax(Xs))
                x_min2=float(numpy.nanmin(Xs))
                xmaxi = max(x_max2,x_min2)
                xmini = min(x_max2,x_min2)
                ymaxi = max(y_max2,y_min2)
                ymini = min(y_max2,y_min2)
                self.axes2.set_xlim(xmaxi,xmini)
                self.axes2.set_ylim(ymaxi,ymini)
                self.ax_reset2=0
            else:
                xmaxi = max(x_max2,x_min2)
                xmini = min(x_max2,x_min2)
                ymaxi = max(y_max2,y_min2)
                ymini = min(y_max2,y_min2)
                self.axes2.set_xlim(xmaxi,xmini)
                self.axes2.set_ylim(ymaxi,ymini)



        ##############################################################333
        #Subplot 3 - the slice
        #self.axes1 = self.fig.add_subplot(gs1[0,:-2])

        self.axes1.clear()
        if(orth_cb and dec_cb):
            self.axes1.set_xlabel(self.slice_service.labels[0]+' (ppm)',fontsize=8)
            self.axes1.set_ylabel(self.slice_service.labels[1]+' (ppm)',fontsize=8)
            Xs,Ys,Zs=self.ReSlice2d(self.slice_service.data,self.inc,sele1,self.slice_service.peaks,Width2,orth=orth_cb)
        else:
            self.axes1.set_xlabel(self.slice_service.labels[0]+' (ppm)',fontsize=8)
            self.axes1.set_ylabel(self.slice_service.labels[2]+' (ppm)',fontsize=8)

            Xs,Ys,Zs=self.ReSlice2d(self.slice_service.data,self.inc,sele1,self.slice_service.peaks,Width1)


        self.axes1.contour(Xs, Ys, Zs, levels_left, cmap=colormap, norm=colors.Normalize(vmin=-numpy.max(levels_left), vmax=numpy.max(levels_left))) #plot pdb network
        self._full_plane = (numpy.asarray(Xs), numpy.asarray(Ys), numpy.asarray(Zs))
        self._draw_full_peak_overlay(self.axes1, str(self.slice_service.peaks[sele1].name),
                                     orth=(orth_cb and dec_cb), show_labels=grid_cb, pane='top')
        self._draw_pane_peak_name(self.axes1, sele1)
        #do the main peak labels
        if(orth_cb and dec_cb):
            #print dimC,dimC2,dimH,dimH2
            self.axes1.scatter(dimC,dimC2,c='g',s=100)
        else:
            self.axes1.scatter(dimC,dimH,c='g',s=100)

        y_max2a=Ys[0][0]
        y_min2a=Ys[(len(Ys))-1][0]
        x_max2a=Xs[0][0]
        x_min2a=Xs[0][(len(Xs[0]))-1]
        #print x_max2a,x_min2a,y_max2a,y_min2a
        if(grid_cb):  # if authoritative peak overlays are enabled
            xl=(x_min2a,x_max2a)
            if(orth_cb and dec_cb):
                hl=(dimC,dimC)
            else:
                hl=(dimH,dimH)
            self.axes1.plot(xl,hl,'green', lw=0.5) #horizontal


            if(orth_cb and dec_cb):
                xl=(dimC2,dimC)
                hl=(dimC2,dimC2)
                self.axes1.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimC)
                xd=(dimC,dimC)
                self.axes1.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimC)
                xd=(dimC2,dimC2)
                self.axes1.plot(xd,yd,'cyan') #vertical 2
            else:
                xl=(dimC,dimC2)
                hl=(dimH,dimH)
                self.axes1.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimH)
                xd=(dimC,dimC)
                self.axes1.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimH)
                xd=(dimC2,dimC2)
                self.axes1.plot(xd,yd,'cyan') #vertical 2



        if(self.ax_reset1==1):
            # Derive limits from the newly extracted coordinate grids.  Do not
            # reuse corner values from the previous slice: depending on axis
            # orientation/transposition those can leave a stale Y range after a
            # ComboBox peak change.
            y_max1=float(numpy.nanmax(Ys))
            y_min1=float(numpy.nanmin(Ys))
            x_max1=float(numpy.nanmax(Xs))
            x_min1=float(numpy.nanmin(Xs))
            xmaxi = max(x_max1,x_min1)
            xmini = min(x_max1,x_min1)
            ymaxi = max(y_max1,y_min1)
            ymini = min(y_max1,y_min1)
            self.axes1.set_xlim(xmaxi,xmini)
            self.axes1.set_ylim(ymaxi,ymini)
            self.ax_reset1=0
        else:
            xmaxi = max(x_max1,x_min1)
            xmini = min(x_max1,x_min1)
            ymaxi = max(y_max1,y_min1)
            ymini = min(y_max1,y_min1)
            self.axes1.set_xlim(xmaxi,xmini)
            self.axes1.set_ylim(ymaxi,ymini)

        # Draw the lower-left X-axis title inside the rectangular plot border,
        # matching the in-axes label treatment used by the projection plots.
        # Place it at the bottom-left rather than below the axis.
        axes2_x_label = self.axes2.get_xlabel()
        self.axes2.set_xlabel('')
        self.axes2.text(0.02, 0.03, axes2_x_label, transform=self.axes2.transAxes,
                        ha='left', va='bottom', fontsize=8)

        #self.fig.tight_layout()
        self.canvas.draw()

    def _draw_projection_static(self, levelsP, colormap, sele1, dimC, dimH, label):
        """Completely rebuild the two projection plots and guide artists."""
        # Preserve a toolbar-selected ROI unless an explicit axis reset was
        # requested (Draw!/peak reset retains the historical reset behaviour).
        old_cc = (self.axesCC.get_xlim(), self.axesCC.get_ylim()) if self.ax_resetCC == 0 else None
        old_hc = (self.axesHC.get_xlim(), self.axesHC.get_ylim()) if self.ax_resetHC == 0 else None

        self.axesCC.clear()
        self.axesHC.clear()
        norm = colors.Normalize(vmin=-numpy.max(levelsP), vmax=numpy.max(levelsP))
        self.axesCC.contour(self.Ys_xy, self.Xs_xy, self.Zs_xy, levelsP, cmap=colormap, norm=norm)
        # Keep projection x ticks on the bottom: the compact layout leaves
        # deliberately little room above the top-right projection.  Replace
        # the conventional x-axis title with a compact nucleus label inside
        # the lower-right corner of the axes.
        cc_x_label = str(self.slice_service.labels[0]).replace('(ppm)', '').strip()
        self.axesCC.set_xlabel('')
        self.axesCC.xaxis.set_ticks_position('bottom')
        self.axesCC.xaxis.set_label_position('bottom')
        self.axesCC.tick_params(axis='x', which='both', top=False, labeltop=False,
                                bottom=True, labelbottom=True)
        self.axesCC.text(0.98, 0.03, cc_x_label, transform=self.axesCC.transAxes,
                         ha='right', va='bottom', fontsize=8)
        self.axesCC.set_ylabel(self.slice_service.labels[1]+' (ppm)', fontsize=8)

        self.axesHC.contour(self.Xs_yz, self.Ys_yz, self.Zs_yz, levelsP, cmap=colormap, norm=norm)
        self.axesHC.text(dimH, dimC, label, fontsize=8)
        self.axesHC.scatter(dimH, dimC, c='k', s=50, marker='x', zorder=2)
        for pk in self.slice_service.peaks:
            self.axesHC.scatter(pk.ppmK, pk.ppmJ, c='k', s=10, marker='x', zorder=2)
        hc_x_label = str(self.slice_service.labels[2]).replace('(ppm)', '').strip()
        self.axesHC.set_xlabel('')
        self.axesHC.xaxis.set_ticks_position('bottom')
        self.axesHC.xaxis.set_label_position('bottom')
        self.axesHC.tick_params(axis='x', which='both', top=False, labeltop=False,
                                bottom=True, labelbottom=True)
        self.axesHC.text(0.98, 0.03, hc_x_label, transform=self.axesHC.transAxes,
                         ha='right', va='bottom', fontsize=8)
        self.axesHC.set_ylabel(self.slice_service.labels[1]+' (ppm)', fontsize=8)

        if self.ax_resetCC == 1:
            x_maxCC = self.Ys_xy[0][0]
            x_minCC = self.Ys_xy[len(self.Ys_xy)-1][0]
            y_maxCC = self.Xs_xy[0][0]
            y_minCC = self.Xs_xy[0][len(self.Xs_xy[0])-1]
            self.axesCC.set_xlim(x_maxCC, x_minCC)
            self.axesCC.set_ylim(y_maxCC, y_minCC)
            self.ax_resetCC = 0
        elif old_cc is not None:
            self.axesCC.set_xlim(*old_cc[0]); self.axesCC.set_ylim(*old_cc[1])

        if self.ax_resetHC == 1:
            y_maxHC = self.Ys_yz[0][0]
            y_minHC = self.Ys_yz[len(self.Ys_yz)-1][0]
            x_maxHC = self.Xs_yz[0][0]
            x_minHC = self.Xs_yz[0][len(self.Xs_yz[0])-1]
            self.axesHC.set_xlim(x_maxHC, x_minHC)
            self.axesHC.set_ylim(y_maxHC, y_minHC)
            self.ax_resetHC = 0
        elif old_hc is not None:
            self.axesHC.set_xlim(*old_hc[0]); self.axesHC.set_ylim(*old_hc[1])

        # Animated artists are excluded from normal full draws and are painted
        # over the cached projection backgrounds by _blit_projection_guides.
        cc_g, = self.axesCC.plot([], [], c='g', animated=True)
        cc_b, = self.axesCC.plot([], [], c='b', animated=True)
        hc_hg, = self.axesHC.plot([], [], c='g', animated=True)
        hc_hb, = self.axesHC.plot([], [], c='b', animated=True)
        hc_vg, = self.axesHC.plot([], [], c='g', animated=True)
        hc_vb, = self.axesHC.plot([], [], c='b', animated=True)
        self._projection_artists = [cc_g, cc_b, hc_hg, hc_hb, hc_vg, hc_vb]
        self._projection_static_ready = True
        self._projection_backgrounds = {}

    def _update_projection_guides(self, sele1, sele2, dimC2, dimH2):
        """Update only the green/blue projection guide-line coordinates."""
        if len(self._projection_artists) != 6:
            return
        yval = self.slice_service.axis(1)[self.slice_service.peaks[sele1].indexJ + self.inc]
        cc_ex = (self.Xs_xy_ymin, self.Xs_xy_ymax)
        hc_hex = (self.Xs_yz_xmin, self.Xs_yz_xmax)
        hc_vex = (self.Xs_yz_ymin, self.Xs_yz_ymax)
        data = [
            (cc_ex, (yval, yval)), (cc_ex, (dimC2, dimC2)),
            (hc_hex, (yval, yval)), (hc_hex, (dimC2, dimC2)),
            ((self.slice_service.peaks[sele1].ppmK, self.slice_service.peaks[sele1].ppmK), hc_vex),
            ((dimH2, dimH2), hc_vex),
        ]
        for artist, (xs, ys) in zip(self._projection_artists, data):
            artist.set_data(xs, ys)
        self._blit_projection_guides()

    def _on_projection_draw_event(self, event):
        """Refresh blit backgrounds after full draw, resize, zoom or pan."""
        if self._projection_blit_busy or not self._projection_static_ready:
            return
        try:
            self._projection_backgrounds = {
                self.axesCC: self.canvas.copy_from_bbox(self.axesCC.bbox),
                self.axesHC: self.canvas.copy_from_bbox(self.axesHC.bbox),
            }
            self._blit_projection_guides()
        except Exception:
            self._projection_backgrounds = {}

    def _blit_projection_guides(self):
        """Paint only projection guide lines over cached static plots."""
        if len(self._projection_artists) != 6 or not self._projection_backgrounds:
            return
        self._projection_blit_busy = True
        try:
            for ax in (self.axesCC, self.axesHC):
                background = self._projection_backgrounds.get(ax)
                if background is None:
                    continue
                self.canvas.restore_region(background)
                for artist in self._projection_artists:
                    if artist.axes is ax and artist.get_visible():
                        ax.draw_artist(artist)
                self.canvas.blit(ax.bbox)
        finally:
            self._projection_blit_busy = False

    def on_cb_grid(self, event):
        self.draw_figure()

    def on_cb_decon(self, event):
        """Toggle the matching deconvolved slice; never launch analysis here."""
        if self.cb_decon.GetValue() and self.slice_service.deconvolution_enabled == 0:
            self.cb_decon.SetValue(False)
            if hasattr(self, 'toolbar'):
                self.toolbar.set_decon_active(False)
            return
        self.inc = 0
        self.inc2 = 0
        self.ax_reset1 = 1
        self.ax_reset2 = 1
        self.draw_figure()

    def on_cb_grid_auto(self, event):
        self.draw_figure()

    def _reset_slice_axes_for_peak_change(self):
        """Reset slice offsets and limits after either peak identifier changes.

        A new reference peak can have a different centre in both displayed
        dimensions.  Keeping any of the old limits makes the contour data change
        while one axis can still describe the previous peak, which is especially
        confusing when navigating from the diagnostics table.
        """
        self.inc = 0
        self.inc2 = 0
        self.ax_reset1 = 1
        self.ax_reset2 = 1
        self.ax_resetCC = 1
        self.ax_resetHC = 1

    def _apply_settings_redraw(self, projection=False):
        """Validate display settings and redraw without changing slice/peak state."""
        try:
            float(self.width1Box.GetValue()); float(self.width2Box.GetValue())
            float(self.textbox0.GetValue()); float(self.textbox1.GetValue()); int(self.textbox2.GetValue())
            float(self.textbox0_right.GetValue()); float(self.textbox1_right.GetValue()); int(self.textbox2_right.GetValue())
            float(self.textbox_minP.GetValue()); float(self.textbox_maxP.GetValue()); int(self.textbox_lvlP.GetValue())
        except (ValueError, TypeError):
            return
        self.selection = []
        self.draw_figure(redraw_projections=bool(projection))

    def on_left_right_setting_changed(self, event):
        obj = event.GetEventObject() if event is not None else None
        value = obj.GetValue() if obj is not None and hasattr(obj, 'GetValue') else None
        self._apply_settings_redraw(False)
        if event is not None:
            event.Skip()

    def on_projection_setting_changed(self, event):
        self._apply_settings_redraw(True)
        if event is not None:
            event.Skip()

    def _set_combo_selection(self, combo, index):
        """Set both the logical and native displayed value of a read-only combo."""
        if index < 0 or index >= combo.GetCount():
            return False
        text = combo.GetString(index)
        combo.SetSelection(index)
        # Some wx native backends can leave the visible text stale after a
        # programmatic SetSelection in a modeless window. SetStringSelection
        # updates the native text control as well as the selection index.
        combo.SetStringSelection(text)
        combo.Refresh()
        combo.Update()
        # Re-assert after the current button event/draw has unwound. This is
        # harmless on platforms that update immediately and fixes stale native
        # combo text on those that defer painting.
        wx.CallAfter(self._refresh_combo_display, combo, index)
        return True

    def _refresh_combo_display(self, combo, index):
        if combo and index >= 0 and index < combo.GetCount():
            combo.SetStringSelection(combo.GetString(index))
            combo.Refresh()
            combo.Update()

    def _sync_slice_spin_controls(self):
        """Keep the modeless Z-slice selectors aligned with inc/inc2."""
        # The compact toolbar SpinButtons intentionally have no displayed value:
        # they navigate ComboBox1/ComboBox2 and therefore must not mirror inc.
        for name, value in (('sliceSpin1', self.inc), ('sliceSpin2', self.inc2)):
            control = getattr(self, name, None)
            if control is not None and control.GetValue() != int(value):
                control.SetValue(int(value))

    def _enable_right_slice_spin(self, enabled):
        """Enable Right-plane stepping and focused-reference navigation as appropriate."""
        control = getattr(self, 'sliceSpin2', None)
        if control is not None:
            control.Enable(bool(enabled))
        # Keep the main toolbar's Right reference navigator available.  It
        # changes ComboBox2 rather than inc2 and is also a target of Show.
        quick = getattr(self, 'quickSliceSpin2', None)
        if quick is not None:
            quick.Enable(True)

    def _set_slice_from_spin(self, side, requested):
        """Apply an absolute spin value only when the corresponding data slice exists."""
        current = int(self.inc if side == 'left' else self.inc2)
        requested = int(requested)
        delta = requested - current
        target = self._slice_target_in_bounds(side, delta)
        if target is None:
            self._sync_slice_spin_controls()
            return False
        if side == 'left':
            self.inc = target
        else:
            self.inc2 = target
        self.selection = []
        self._sync_slice_spin_controls()
        self.draw_figure()
        return True

    def on_slice_spin_left(self, event):
        self._set_slice_from_spin('left', event.GetEventObject().GetValue())

    def on_slice_spin_right(self, event):
        self._set_slice_from_spin('right', event.GetEventObject().GetValue())

    def _slice_target_in_bounds(self, side, delta):
        """Return the proposed slice offset if it indexes the 3-D data safely."""
        arr = self.slice_service.data
        orth = bool(self.cb_flip.GetValue())
        decon = bool(self.cb_decon.GetValue())
        if side == 'left':
            peak_index = self.ComboBox1.GetSelection()
            current = self.inc
            # axes1 uses an orthogonal (H-index) slice only for Orth+Decon.
            use_h_axis = orth and decon
        else:
            current = self.inc2
            if orth:
                peak_index = self.ComboBox1.GetSelection()
                use_h_axis = True
            else:
                peak_index = self.ComboBox2.GetSelection()
                use_h_axis = False
        if peak_index == wx.NOT_FOUND:
            return None
        base = self.slice_service.peak_indices[peak_index][1 if use_h_axis else 0]
        axis_len = arr.shape[2 if use_h_axis else 1]
        target = current + delta
        absolute = base + target
        if absolute < 0 or absolute >= axis_len:
            return None
        return target

    def _move_peak_selection(self, combo, delta, reset_left=False, reset_right=False):
        """Move a peak combo by one item without ever creating an invalid selection."""
        count = combo.GetCount()
        current = combo.GetSelection()
        if count <= 0:
            return False
        if current == wx.NOT_FOUND:
            current = 0 if delta >= 0 else count - 1
        target = current + delta
        if target < 0 or target >= count:
            return False
        if reset_left:
            self.ax_reset1 = 1
            self.inc = 0
            if self.cb_flip.GetValue() or self.cb_decon.GetValue():
                self.ax_reset2 = 1
        if reset_right:
            self.ax_reset1 = 1
            self.ax_reset2 = 1
            self.inc2 = 0
        result = self._set_combo_selection(combo, target)
        # Programmatic selection does not emit EVT_COMBOBOX, so redraw here.
        self.selection = []
        self.draw_figure()
        return True

    def on_peak_combo_changed(self, event):
        """Immediately redraw when a Left/Right peak identifier is selected."""
        self._reset_slice_axes_for_peak_change()
        self.selection = []
        self.draw_figure()
        if event is not None:
            event.Skip()

    def _toolbar_orth(self, active):
        """Toggle the orthogonal display and redraw the affected slice panes."""
        self.cb_flip.SetValue(bool(active))
        self.inc = 0
        self.inc2 = 0
        self.ax_reset1 = 1
        self.ax_reset2 = 1
        # Projection contours remain cached; draw_figure only rebuilds the
        # slice panes and the existing animated projection guides are blitted.
        self.draw_figure()

    def _toolbar_1d(self, active):
        """Toggle the lower 1D trace with an immediate display update."""
        self.cb_1d.SetValue(bool(active))
        self.ax_reset2 = 1
        # Keep projection contours cached and use their existing blit path.
        self.draw_figure()

    def _toolbar_decon(self, active):
        self.cb_decon.SetValue(bool(active))
        self.on_cb_decon(None)

    def _toolbar_peaks(self, active):
        self.cb_grid_auto.SetValue(bool(active))
        self.on_cb_grid_auto(None)

    def _toolbar_contours(self):
        self._show_tool_window(self.projContourFrame)

    def redraw_view(self):
        # Explicit redraw is the projection refresh point: rebuild both contour
        # projections and then capture fresh blit backgrounds for their guides.
        self._reset_slice_axes_for_peak_change()
        self.selection = []
        self.draw_figure(redraw_projections=True)

    def on_draw_button(self, event):
        self.redraw_view()
        if event is not None and hasattr(event, 'Skip'):
            event.Skip()

    def on_N_button(self, event):
        self._move_peak_selection(self.ComboBox1, +1, reset_left=True)

    def on_P_button(self, event):
        self._move_peak_selection(self.ComboBox1, -1, reset_left=True)

    def on_quick_peak_up_left(self, event):
        """Toolbar up arrow: advance the focused Left/reference peak."""
        self._move_peak_selection(self.ComboBox1, +1, reset_left=True)

    def on_quick_peak_down_left(self, event):
        """Toolbar down arrow: move to the previous Left/reference peak."""
        self._move_peak_selection(self.ComboBox1, -1, reset_left=True)

    def on_quick_peak_up_right(self, event):
        """Toolbar up arrow: advance the focused Right/reference peak."""
        self._move_peak_selection(self.ComboBox2, +1, reset_right=True)

    def on_quick_peak_down_right(self, event):
        """Toolbar down arrow: move to the previous Right/reference peak."""
        self._move_peak_selection(self.ComboBox2, -1, reset_right=True)


    def on_swap_button(self, event):
        self.ax_reset1=1
        self.ax_reset2=1
        tmp=self.ComboBox1.GetSelection()
        self._set_combo_selection(self.ComboBox1, self.ComboBox2.GetSelection())
        self._set_combo_selection(self.ComboBox2, tmp)
        self.selection = []
        self.draw_figure()

    def on_Up_button(self, event):
        target = self._slice_target_in_bounds('left', +1)
        if target is not None:
            self.inc = target
            self._sync_slice_spin_controls()
            self.selection = []
            self.draw_figure()

    def on_Down_button(self, event):
        target = self._slice_target_in_bounds('left', -1)
        if target is not None:
            self.inc = target
            self._sync_slice_spin_controls()
            self.selection = []
            self.draw_figure()

    def on_full_peak_list(self, event):
        """Open the authoritative Full Peak List from the Slice2D toolbar."""
        if self.slice_service is None:
            logging.warning('Full Peak List viewer is not available from Slice2D')
            return None
        return self.slice_service.open_full_peak_list(event)

    def on_N_button2(self, event):
        self._move_peak_selection(self.ComboBox2, +1, reset_right=True)

    def on_P_button2(self, event):
        self._move_peak_selection(self.ComboBox2, -1, reset_right=True)


    def on_Up_button2(self, event):
        target = self._slice_target_in_bounds('right', +1)
        if target is not None:
            self.inc2 = target
            self._sync_slice_spin_controls()
            self.selection = []
            self.draw_figure()

    def on_Down_button2(self, event):
        target = self._slice_target_in_bounds('right', -1)
        if target is not None:
            self.inc2 = target
            self._sync_slice_spin_controls()
            self.selection = []
            self.draw_figure()




    def on_pick(self, event):
        """Handle authoritative Full Peak List interactions only."""
        self._handle_full_tool_click(event)

    def on_text_enter(self, event):
        self.draw_figure()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)

    def on_exit(self, event):
        self.Destroy()

    def on_about(self, event):
        msg = """ A demo using wxPython with matplotlib:

         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw!")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER,
            self.on_flash_status_off,
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
