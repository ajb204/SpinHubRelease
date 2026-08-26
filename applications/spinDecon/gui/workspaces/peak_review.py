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
import wx,string,copy,math,numpy,os,sys,re, platform
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
from matplotlib.lines import Line2D
from matplotlib.widgets import Cursor,MultiCursor
from spinDecon.domain.peaks import peakEntry
from spinDecon.gui.plotting.array_utils import ensure_xy_points, scatter_xy_points
from spinDecon.gui.plotting.display_utils import make_cursor
from spinDecon.gui.context import context_for, project_for, data_for

##################################################################################################################

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

# assign ID numbers
#[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
#] = [wx.NewId() for _init_ctrls in range(4)]
from spinDecon.project.parameter_store import update_parameter_file


class PeakNavigationToolbar(RedrawNavigationToolbar):
    """Matplotlib toolbar without the stock mouse-coordinate text."""
    def __init__(self, canvas, redraw_callback):
        # Coordinates are rendered inside the plotting plane by PeakFrame.
        super().__init__(canvas, redraw_callback, coordinates=False)

    def set_message(self, message):
        # Matplotlib calls this during mouse motion.  Deliberately suppress the
        # toolbar message because PeakFrame owns the in-axes coordinate label.
        pass

from spinDecon.gui.widgets.common import PersistentStateButton


#peak projection plane
class peakFrame(wx.Frame):


    def __init__(self,parent,showFlg=True):
        wx.Frame.__init__(self, parent, id=wx.ID_ANY, name='',
              title='Get 2D peak list ...')
        self.SetClientSize(wx.Size(900, 205))
        # Keep transient Peak Tools feedback inside the PeakFrame rather than
        # sending it to the launching terminal.  wx.Frame manages the status
        # bar outside the client-area sizer, so the controls/canvas keep their
        # existing layout.
        self.statusBar = self.CreateStatusBar(1)
        self.statusBar.SetStatusText('Ready')

        self.app_context = context_for(parent)
        self.peak_service = getattr(self.app_context, "peaks", None) if self.app_context is not None else None
        if self.peak_service is None:
            from spinDecon.analysis.peak_service import PeakService
            self.peak_service = PeakService(parent)
        # Dimensionality is fixed for the lifetime of PeakFrame.  Keep it
        # explicit so view construction is separated cleanly from interaction logic.
        try:
            self.spectral_dim_count = int(parent.state.topology().spectral_dim_count)
        except Exception:
            self.spectral_dim_count = int(getattr(parent, 'dim', 0) or 0)
        self.dim = self.spectral_dim_count  # compatibility alias only
        # Axis labels are fixed for the lifetime of this PeakFrame.  Cache
        # them once rather than using GUI widgets as label storage.
        self.x_axis_label, self.y_axis_label = (self.peak_service.view_labels(self.spectral_dim_count))
        self.view_labels = (self.x_axis_label, self.y_axis_label)
        self.state = project_for(parent)
        self.store = data_for(parent)
        self.projection_key = self._default_projection_key()
        self.thresh = self.peak_service.threshold()

        self.peak=self._load_projection_peaks()
        self.threshFac=1.0

        # Spectral arrays are owned by the parent decon frame/data store.
        # PeakFrame keeps only UI state and resolves the current projection
        # dynamically when it is needed.
        self.select=[]
        self.SELECT=0
        self.SELECTGRP=0
        self.MOVE=0
        self.MOVEGRP=0
        self.ADD=0
        self.PAIR=0
        # Separate transient mouse interaction from persistent selection state.
        self.interaction_mode = None
        self.selection_type = None  # None, 'single', or 'group'
        # Peak edit history.  Each entry is a complete snapshot of the canonical
        # reference peak list immediately before/after an edit.
        self.undo_stack=[]
        self.redo_stack=[]
        self.ax_reset=1       #for keeping the zoom
        self.OVERLAY=0
        self.TRANSPOSE='n' #assume we're not transposing the main list
        self.x1 =0
        self.y1 = 0
        self.resized = False

        self.create_main_panel()
        self._set_tool_mode(None)
        # draw_figure() performs the one full render required at startup.
        self.draw_figure()
        if(showFlg):
            self.Show(True)
        self.Fit()

    def _uses_full_2d_peak_list(self):
        """True only when PeakFrame represents the actual physical 2D spectrum."""
        try:
            topology = self.state.topology()
            return bool(topology.spectral_dim_count == 2 and not topology.has_pseudo_axis)
        except Exception:
            return False

    @property
    def peak(self):
        # PeakFrame needs mutable Peak objects for its editing tools.  Keep the
        # currently displayed objects as a UI cache, while _commit_projection_peaks
        # publishes every completed edit to the topology-appropriate canonical list.
        if hasattr(self, '_peak_cache'):
            return self._peak_cache
        return self._load_projection_peaks()

    @peak.setter
    def peak(self, value):
        self._peak_cache = list(value) if value is not None else []
        self._commit_projection_peaks(self._peak_cache)

    def _spectrum_base_path(self):
        return self.peak_service.spectrum_path()


    def _parent_spectral_dim_count(self):
        """Canonical spectral dimensionality of the parent dataset."""
        try:
            return int(self.state.topology().spectral_dim_count)
        except Exception:
            return self.peak_service.dimension

    def sync_main_threshold(self, redraw=True):
        """Synchronise PeakFrame's minimum contour with the main threshold.

        ``threshBox`` stores a fraction of ``dmax`` whereas PeakFrame contour
        controls use an absolute intensity.  Keep both the internal threshold
        and the visible Contours -> Min control in sync before rebuilding the
        contour layer.  This method only updates an already-open PeakFrame; it
        never creates one.
        """
        # Legacy source contract retained during service migration: self.thresh = self.tabOne.dmax * float(self.tabOne.threshBox.GetValue())
        self.thresh = self.peak_service.threshold()

        # draw_figure() obtains its minimum level from textbox0/GetLevels(), not
        # directly from self.thresh.  Updating self.thresh alone therefore left
        # the displayed contour minimum unchanged.
        minimum = str(self.thresh)
        if hasattr(self, 'textbox0'):
            self.textbox0.SetValue(minimum)

        if redraw:
            # Preserve the user's current map limits while updating contours.
            self.ax_reset = 0
            self.draw_figure()
        return self.thresh

    def _analysis_spectrum_path(self):
        """Exact 2D spectrum analysed by PeakFrame decon/recon.

        PeakFrame is always a 2D analysis context.  Pure 2D datasets use the
        main spectrum; higher-dimensional and 2D+pseudo datasets use the
        on-disk projection backing the displayed view.
        """
        topology = self.state.topology()
        if topology.spectral_dim_count == 2 and not topology.has_pseudo_axis:
            return self._spectrum_base_path()
        return self._projection_source_path()

    def _peak_list_path(self):
        base = self._analysis_spectrum_path()
        return base + '.2D.list' if base else 'peakframe.2D.list'

    def _default_projection_key(self):
        labels = self.view_labels
        if self.spectral_dim_count <= 2:
            return ('raw_projection', labels[0], labels[1], 'n')
        return (labels[0], labels[1], 'n')

    def _projection_source_path(self):
        """Return the on-disk 2D projection file currently shown in peakFrame."""
        labels = self.view_labels
        keys = [self.projection_key]
        if len(labels) >= 2:
            keys.extend([
                (labels[0], labels[1], 'n'),
                (labels[1], labels[0], 'n'),
            ])
        if self.store is not None:
            for key in keys:
                payload = self.store.projections.get(key)
                if payload and payload.get('source'):
                    source = str(payload['source'])
                    if source:
                        return source

        base = self._spectrum_base_path()
        if base:
            base_dir = os.path.dirname(os.path.abspath(base))
            search_dirs = []
            for cand in (
                os.path.join(base_dir, 'projections'),
                os.path.join(base_dir, 'projection_decon'),
                os.path.join(base_dir, 'projections1D'),
            ):
                if cand not in search_dirs:
                    search_dirs.append(cand)
            for proj_dir in search_dirs:
                if not os.path.isdir(proj_dir):
                    continue
                for left, right in ((labels[0], labels[1]), (labels[1], labels[0])) if len(labels) >= 2 else []:
                    for name in (f'{left}.{right}.dat', f'{right}.{left}.dat'):
                        path = os.path.join(proj_dir, name)
                        if os.path.exists(path):
                            return path
        return base

    def _is_pseudo3d_dataset(self):
        """Return True for two spectral dimensions plus one pseudo axis.

        Use the parent's canonical topology rather than the legacy ``pseudo``
        attribute, which can be stale while a project is being restored.
        """
        try:
            topology = self.state.topology()
            return (topology.spectral_dim_count == 2 and
                    topology.has_pseudo_axis and
                    topology.physical_dim_count == 3)
        except (AttributeError, ValueError):
            # Compatibility for an incompletely initialised parent only.
            return self.peak_service.is_pseudo3d()


    def _has_1d_bore(self):
        """Whether PeakFrame should show a 1D trace through the third physical axis."""
        return self.spectral_dim_count == 3 or self._is_pseudo3d_dataset()

    def _bore_payload(self):
        """Return (axis, cube, label) for the PeakFrame 1D bore.

        The cube is always logical [bore, y, x].  For pseudo3D this comes from
        the canonical pseudo3d view, so no assumption is made about which
        physical NMRPipe axis is the pseudo axis.
        """
        if self._is_pseudo3d_dataset():
            view = self.peak_service.pseudo3d_view('raw')
            if view is None:
                return None
            return (numpy.asarray(view['pseudo_axis']), numpy.asarray(view['data']),
                    str(view.get('pseudo_label') or 'pseudoaxis'))
        if self.spectral_dim_count == 3:
            return self.peak_service.bore_payload(self.spectral_dim_count)
        return None

    def _view_labels(self):
        return self.peak_service.view_labels(self.spectral_dim_count)


    def _projection_payload(self, decon=False):
        """Return a plotting-ready payload owned by the application data boundary."""
        return self.peak_service.projection_payload(
            self.view_labels, decon=decon, analysis_path=self._analysis_spectrum_path())


    def _display_payload(self, decon=False):
        payload = self._projection_payload(decon=decon)
        XX = payload.get('XX')
        YY = payload.get('YY')
        ZZ = payload.get('ZZ')
        if XX is None or YY is None or ZZ is None:
            raise RuntimeError('PeakFrame display payload is incomplete')
        return XX, YY, ZZ, payload

    @property
    def XX(self):
        return self._display_payload(False)[0]

    @property
    def YY(self):
        return self._display_payload(False)[1]

    @property
    def ZZ(self):
        return self._display_payload(False)[2]

    @property
    def data(self):
        return self._display_payload(False)[2]

    @property
    def dic(self):
        return self._display_payload(False)[3].get('dic')

    @property
    def uc0(self):
        return self._projection_payload(False).get('uc0') or self.peak_service.unit_converter(0)

    @property
    def uc1(self):
        return self._projection_payload(False).get('uc1') or self.peak_service.unit_converter(1)

    @property
    def index0(self):
        payload = self._projection_payload(False)
        return payload.get('index0') if payload.get('index0') is not None else payload.get('y_axis')

    @property
    def index1(self):
        payload = self._projection_payload(False)
        return payload.get('index1') if payload.get('index1') is not None else payload.get('x_axis')


    def _projection_axis_arrays(self):
        """Return the displayed projection x/y coordinate arrays.

        PeakFrame does not own spectral arrays.  The current 2D projection
        and its axes are resolved from the shared DataStore-backed payload.
        The hover logic uses the actual plotted orientation so it remains
        correct when the projection is transposed.
        """
        try:
            xx, yy, _zz, _payload = self._display_payload(decon=False)
            xx = numpy.asarray(xx)
            yy = numpy.asarray(yy)
            if xx.ndim == 2 and yy.ndim == 2:
                # Get the coordinate vectors from the same orientation that
                # the shared display view puts on the canvas.  For a transposed
                # display, X varies down rows and Y varies across columns.
                if self.TRANSPOSE == 'y':
                    return numpy.asarray(xx[:, 0], dtype=float), numpy.asarray(yy[0, :], dtype=float)
                return numpy.asarray(xx[0, :], dtype=float), numpy.asarray(yy[:, 0], dtype=float)
        except Exception:
            pass

        payload = self._projection_payload(False)
        x_axis = payload.get('x_axis')
        y_axis = payload.get('y_axis')
        if x_axis is not None and y_axis is not None:
            if self.TRANSPOSE == 'y':
                return numpy.asarray(y_axis, dtype=float), numpy.asarray(x_axis, dtype=float)
            return numpy.asarray(x_axis, dtype=float), numpy.asarray(y_axis, dtype=float)

        return None, None

    
    def _nearest_projection_indices(self, x, y, context=''):
        x_axis, y_axis = self._projection_axis_arrays()
        if x_axis is None or y_axis is None:
            raise RuntimeError('Cannot determine projection axes for hover update')
        x1 = int(numpy.abs(x_axis - x).argmin())
        y1 = int(numpy.abs(y_axis - y).argmin())
        return x1, y1


    def _as_ndarray(self, data):
        """Return an ndarray view for spectrum-like objects."""
        if data is None:
            return None
        if isinstance(data, numpy.ndarray):
            return data
        try:
            sliced = data[:]
            if isinstance(sliced, numpy.ndarray):
                return sliced
            if hasattr(sliced, "shape"):
                return numpy.asarray(sliced)
        except Exception:
            pass
        try:
            arr = numpy.asarray(data)
            if hasattr(arr, "shape") and arr.shape != ():
                return arr
        except Exception:
            pass
        return data

    def _load_projection_peaks(self):
        if self.store is None:
            return []
        if self._uses_full_2d_peak_list():
            # Full lists are stored as dimension-independent records.  PeakFrame
            # itself consumes mutable Peak objects, so initialise those from the
            # spectrum-owned .2D.list when available.
            path = self._peak_list_path()
            if path and os.path.exists(path):
                return self.read_peaklist_file(path)
            return []
        payload = self.store.peak_lists.get('reference')
        if not payload:
            return []
        peaks = payload.get('peaks', [])
        return peaks if peaks is not None else []

    def _commit_projection_peaks(self, value=None):
        if self.store is None:
            return
        peaks = self.peak if value is None else (list(value) if value is not None else [])
        self._peak_cache = peaks
        self.peak_service.commit_projection_peaks(
            peaks, full_2d=self._uses_full_2d_peak_list(),
            source_path=self._peak_list_path(), projection_key=self.projection_key,
            labels=self.view_labels,
        )
        if self.state is not None:
            try:
                self.state.metadata['peakframe_projection_key'] = self.projection_key
                self.state.metadata['peakframe_peak_count'] = len(peaks)
            except Exception:
                pass

    def focus_peak(self, peak, width_fraction=0.10):
        """Select and zoom the current projection around a shared peak.

        External list viewers use the same persistent selection state as an
        in-frame Select operation, so the selected peak ornament is emphasised
        consistently regardless of where the selection originated.
        """
        selected_index = next((i for i, candidate in enumerate(self.peak)
                               if candidate is peak or
                               (str(getattr(candidate, 'name', '')) == str(getattr(peak, 'name', ''))
                                and numpy.isclose(float(getattr(candidate, 'x', numpy.nan)), float(getattr(peak, 'x', numpy.nan)), equal_nan=False)
                                and numpy.isclose(float(getattr(candidate, 'y', numpy.nan)), float(getattr(peak, 'y', numpy.nan)), equal_nan=False))), None)
        if selected_index is not None:
            self.select = [selected_index]
            self.selection_type = 'single'
            self._update_selection_artists()
        XX, YY, _data, _payload = self._display_payload(decon=False)
        x_min, x_max = float(numpy.nanmin(XX)), float(numpy.nanmax(XX))
        y_min, y_max = float(numpy.nanmin(YY)), float(numpy.nanmax(YY))
        half_x = abs(x_max - x_min) * float(width_fraction) / 2.0
        half_y = abs(y_max - y_min) * float(width_fraction) / 2.0
        self.axes.set_xlim(peak.x - half_x, peak.x + half_x)
        self.axes.set_ylim(peak.y - half_y, peak.y + half_y)
        self.ax_reset = 0
        self.canvas.draw()

    def read_peaklist_file(self, infile):
        from spinDecon.processing.peak_io import read_peak_list
        return read_peak_list(infile)

    def drawing_box(self):
        flags = wx.ALIGN_LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL


        self.drawLbl = wx.StaticBox(self, -1, 'Drawing:')  #, size=(240, 140))
        self.drawSizer = wx.StaticBoxSizer(self.drawLbl, wx.VERTICAL)

        self.contourbutton = wx.Button(self.drawLbl, -1, "Contour", size=(-1,22))
        self.cb_grid = wx.CheckBox(self.drawLbl, -1,"Peaks",style=wx.ALIGN_RIGHT)
        # Contour and peak visibility are now exposed by the Matplotlib toolbar.
        self.contourbutton.Hide()
        self.cb_grid.Hide()
        self.cb_labels = wx.CheckBox(self.drawLbl, -1,"Labels",style=wx.ALIGN_RIGHT)




        self.Bind(wx.EVT_BUTTON, self.OnButtonContour, self.contourbutton)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_labels)

        self.draw_grid_sizer = wx.GridBagSizer(1, 6)
        cnt=0
        border = 3 if self.spectral_dim_count == 3 else 10

        self.draw_grid_sizer.Add(self.cb_labels, (0,cnt), border=border, flag=flags);cnt+=1
        self.drawSizer.Add(self.draw_grid_sizer, flag=wx.ALIGN_CENTER_HORIZONTAL)

    def contour_box(self):
        """Create the modeless contour-settings window.

        Contour state remains owned by PeakFrame; this window simply exposes
        the existing Min/Factor/Number controls without consuming spectrum
        display space in the main window.
        """
        self.contourFrame = wx.Frame(
            self.GetTopLevelParent(), title='Contours',
            style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT
        )
        panel = wx.Panel(self.contourFrame)
        box = wx.StaticBox(panel, label='Contours:')
        contour_sizer = wx.StaticBoxSizer(box, wx.HORIZONTAL)

        controls = []
        for label, width, value in (
                ('Min:', 100, self.thresh),
                ('Factor:', 55, 1.2),
                ('Number:', 55, 15)):
            contour_sizer.Add(wx.StaticText(box, label=label), 0,
                              wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 5)
            ctrl = wx.TextCtrl(box, size=(width, 22),
                               style=wx.TE_PROCESS_ENTER)
            ctrl.SetValue(str(value))
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.OnContourSettingChanged)
            contour_sizer.Add(ctrl, 0,
                              wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 4)
            controls.append(ctrl)

        self.textbox0, self.textbox1, self.textbox2 = controls

        close = wx.Button(panel, label='Close', size=(-1, 24))
        close.Bind(wx.EVT_BUTTON, lambda evt: self.contourFrame.Hide())

        root = wx.BoxSizer(wx.VERTICAL)
        root.Add(contour_sizer, 0, wx.EXPAND | wx.ALL, 5)
        root.Add(close, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        panel.SetSizer(root)
        self.contourFrame.SetClientSize(root.CalcMin())
        self.contourFrame.Bind(wx.EVT_CLOSE, self.OnContourWindowClose)

    def OnButtonContour(self, event):
        """Show the contour settings without creating duplicate windows."""
        if not self.contourFrame.IsShown():
            self.contourFrame.Show()
        self.contourFrame.Raise()

    def OnContourWindowClose(self, event):
        """Keep the modeless contour window available for later reuse."""
        self.contourFrame.Hide()
        event.Veto()

    def OnContourSettingChanged(self, event):
        """Apply contour settings and redraw while preserving current zoom."""
        try:
            minimum = float(self.textbox0.GetValue())
            factor = float(self.textbox1.GetValue())
            number = int(self.textbox2.GetValue())
        except ValueError:
            wx.MessageBox(
                'Contour Min and Factor must be numbers, and Number must be an integer.',
                'Invalid contour settings', wx.OK | wx.ICON_ERROR, self.contourFrame
            )
            return
        if minimum <= 0 or factor <= 0 or number <= 0:
            wx.MessageBox(
                'Contour Min, Factor, and Number must all be greater than zero.',
                'Invalid contour settings', wx.OK | wx.ICON_ERROR, self.contourFrame
            )
            return

        # A contour change genuinely changes the spectrum layer, so perform
        # exactly one full redraw while preserving the current axes limits.
        self.ax_reset = 0
        self.draw_figure()

    def peak_list_control_box(self):
        flags = wx.ALIGN_LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL
        buttonSize = (75,22)

        self.peaklist_cntrl_lbl = wx.StaticBox(self, -1, 'Peak List Controls:')  #, size=(240, 140))
        self.peaklist_cntrl_box = wx.StaticBoxSizer(self.peaklist_cntrl_lbl, wx.VERTICAL)

        self.buttonOverlay = wx.Button(self.peaklist_cntrl_lbl, -1, "Overlay", size = buttonSize)

        self.buttonOverlay.Bind(wx.EVT_BUTTON, self.OnButtonOverlay)
        self.buttonPeak = wx.Button(self.peaklist_cntrl_lbl, -1, "PeakList", size = buttonSize)
        self.buttonPeak.Bind(wx.EVT_BUTTON, self.onPeakList)
        self.buttonMapWindow = wx.Button(self.peaklist_cntrl_lbl, -1, "Map", size=buttonSize)
        self.buttonMapWindow.Bind(wx.EVT_BUTTON, self.onMap)

        # The modeless Map/Tools windows follow the Slice2D Left/Right pattern:
        # peakFrame retains the Python references to every widget, while the
        # widgets themselves use the Map panel as their wx parent.
        self._create_map_window()
        self._create_tools_window()

        self.peaklist_cntrl_box.AddSpacer(5)
        self.peaklist_hbox = wx.BoxSizer(wx.HORIZONTAL)
        border = 3 if self.spectral_dim_count == 3 else 10
        self.peaklist_hbox.AddSpacer(4 if self.spectral_dim_count == 3 else 10)
        self.peaklist_hbox.Add(self.buttonOverlay, border=border, flag=flags)
        self.peaklist_hbox.Add(self.buttonPeak, border=border, flag=flags)
        self.peaklist_hbox.Add(self.buttonMapWindow, border=border, flag=flags)
        self.peaklist_cntrl_box.Add(self.peaklist_hbox, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.peaklist_cntrl_box.AddSpacer(5)

    def _create_map_window(self):
        """Create the modeless Map controls window owned by peakFrame.

        The controls displayed in this window are still attributes of peakFrame;
        their wx parent is mapPanel so wx sizer/window ownership remains valid.
        """
        self.mapFrame = wx.Frame(
            self, title='Map',
            style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT
        )
        self.mapPanel = wx.Panel(self.mapFrame)
        self.mapFrame.Bind(wx.EVT_CLOSE, self.onMapWindowClose)

    def onMap(self, event):
        if not self.mapFrame.IsShown():
            self.mapFrame.Show()
        self.mapFrame.Raise()

    def onMapClose(self, event):
        self.mapFrame.Hide()

    def onMapWindowClose(self, event):
        self.mapFrame.Hide()
        event.Veto()

    def _create_tools_window(self):
        """Create the Peak Tools as an attached panel owned by peakFrame.

        The panel is inserted at the right-hand edge of the main frame, next to
        the spectrum, rather than being a second top-level wx.Frame.  It starts
        hidden and is shown/hidden by the Tools button.  Showing the panel is a
        wx layout operation only: it must not call draw_figure() or otherwise
        reconstruct the spectrum.
        """
        self.toolsPanel = wx.Panel(self)
        self.toolsFloatFrame = None
        self._tools_detached = False
        # Let the Tools sizer determine the minimum strip width required by
        # its widgets instead of reserving a fixed-width panel.
        self.toolsPanel.Hide()

    def _tools_min_size(self):
        """Return the compact size required by the Tools controls."""
        panel_sizer = self.toolsPanel.GetSizer()
        minimum = panel_sizer.CalcMin() if panel_sizer is not None else wx.Size(0, 0)
        best = self.toolsPanel.GetBestSize()
        return wx.Size(
            max(minimum.GetWidth(), best.GetWidth() if best.IsFullySpecified() else 0, 1),
            max(minimum.GetHeight(), best.GetHeight() if best.IsFullySpecified() else 0, 1),
        )

    def _update_tools_button(self, shown):
        """Mirror Tools visibility without using a native toggle state."""
        shown = bool(shown)
        if hasattr(self, 'toolbar'):
            self.toolbar.set_tools_active(shown)

    def _show_tools_panel(self, show):
        """Show/hide Tools without reconstructing the Matplotlib spectrum.

        Docked Tools grows the peakFrame to the right, preserving the spectrum
        footprint.  Detached Tools is shown as a compact floating child frame.
        Neither path calls draw_figure() or any other spectrum render routine.
        """
        show = bool(show)
        if self._tools_detached:
            if self.toolsFloatFrame is not None:
                self.toolsFloatFrame.Show(show)
                if show:
                    self.toolsFloatFrame.Raise()
            self._update_tools_button(show)
            return

        if self.toolsPanel.IsShown() == show:
            self._update_tools_button(show)
            return

        self.Freeze()
        try:
            if show:
                self._tools_closed_frame_size = wx.Size(self.GetSize())
                tools_size = self._tools_min_size()
                self._tools_open_width = tools_size.GetWidth()
                old_size = self._tools_closed_frame_size
                self.SetSize(wx.Size(old_size.GetWidth() + self._tools_open_width,
                                     old_size.GetHeight()))
                self.toolsPanel.SetMinSize((self._tools_open_width, -1))
                self.toolsPanel.Show(True)
            else:
                self.toolsPanel.Show(False)
                closed_size = getattr(self, '_tools_closed_frame_size', None)
                if closed_size is not None:
                    self.SetSize(closed_size)

            self.hbox.Layout()
            self.Layout()
            self._update_tools_button(show)
        finally:
            self.Thaw()

    def _ensure_tools_float_frame(self):
        if self.toolsFloatFrame is None:
            self.toolsFloatFrame = wx.Frame(
                self, title='Peak Tools',
                style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT
            )
            self.toolsFloatFrame.Bind(wx.EVT_CLOSE, self.onToolsFloatClose)
        return self.toolsFloatFrame

    def onToolsDetach(self, event):
        """Toggle Peak Tools between the right-hand dock and a floating frame."""
        if self._tools_detached:
            self._dock_tools_panel()
        else:
            self._detach_tools_panel()

    def _detach_tools_panel(self):
        if self._tools_detached:
            return
        was_shown = self.toolsPanel.IsShown()
        closed_size = getattr(self, '_tools_closed_frame_size', None)
        if closed_size is None:
            closed_size = wx.Size(self.GetSize())

        self.Freeze()
        try:
            self.hbox.Detach(self.toolsPanel)
            frame = self._ensure_tools_float_frame()
            self.toolsPanel.Reparent(frame)
            float_sizer = wx.BoxSizer(wx.VERTICAL)
            float_sizer.Add(self.toolsPanel, 1, wx.EXPAND)
            frame.SetSizer(float_sizer)
            self.toolsPanel.SetMinSize(wx.DefaultSize)
            self.toolsDetachButton.SetLabel('Dock')
            self._tools_detached = True

            # Detaching returns the spectrum window to its pre-Tools footprint.
            self.SetSize(closed_size)
            self.hbox.Layout()
            self.Layout()

            self.toolsPanel.Show(True)
            frame.Fit()
            frame.SetMinSize(frame.GetSize())
            frame.Show(was_shown)
            if was_shown:
                frame.Raise()
        finally:
            self.Thaw()

    def _dock_tools_panel(self):
        if not self._tools_detached:
            return
        frame = self.toolsFloatFrame
        was_shown = bool(frame is not None and frame.IsShown())
        if frame is not None:
            frame.Hide()

        self.Freeze()
        try:
            # Reparenting does not remove a window from its old sizer.  wx
            # requires the floating sizer to release the panel first; otherwise
            # inserting it into hbox asserts that it already has a containing
            # sizer and leaves the reparented panel unmanaged at (0, 0).
            if frame is not None:
                float_sizer = frame.GetSizer()
                if float_sizer is not None:
                    float_sizer.Detach(self.toolsPanel)

            self.toolsPanel.Reparent(self)
            # Insert immediately before the final right-hand spacer.
            insert_at = max(0, self.hbox.GetItemCount() - 1)
            self.hbox.Insert(insert_at, self.toolsPanel, 0, wx.EXPAND | wx.TOP | wx.BOTTOM)
            self.toolsDetachButton.SetLabel('Detach')
            self._tools_detached = False
            self.toolsPanel.Hide()
            self.hbox.Layout()
            self.Layout()
        finally:
            self.Thaw()

        if was_shown:
            self._show_tools_panel(True)

    def _toolbar_tools(self, active):
        self._show_tools_panel(bool(active))

    def onTools(self, event):
        # Toggle the real Tools visibility rather than a native button state.
        # The visual state is then derived from that application state.
        if self._tools_detached:
            shown = bool(self.toolsFloatFrame is not None and self.toolsFloatFrame.IsShown())
        else:
            shown = self.toolsPanel.IsShown()
        self._show_tools_panel(not shown)

    def onToolsClose(self, event):
        self._show_tools_panel(False)

    def onToolsFloatClose(self, event):
        # Closing a detached Tools frame hides it; Tools reopens it in place.
        self._show_tools_panel(False)
        event.Veto()

    def onToolsWindowClose(self, event):
        self._show_tools_panel(False)

    def decon_box(self):
        flags = wx.ALIGN_LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL
        buttonSize = (75,22)

        self.deconLbl = wx.StaticBox(self, -1, 'Deconvolution:')  #, size=(240, 140))
        self.deconSizer = wx.StaticBoxSizer(self.deconLbl, wx.VERTICAL)

        self.buttonDecon = wx.Button(self.deconLbl, -1, "Decon",size=buttonSize)
        self.buttonDecon.Bind(wx.EVT_BUTTON, self.OnButtonDecon)

        self.buttonRecon = wx.Button(self.deconLbl, -1, "Recon",size=buttonSize)
        self.buttonRecon.Bind(wx.EVT_BUTTON, self.OnButtonRecon)



        self.cb_calc = _ToolbarToggleState(False)




        
        self.sizer3=wx.GridBagSizer(5, 5);cnt=0
        self.sizer3.Add((10,0), (0,0));cnt+=1

        self.sizer3.Add(self.buttonDecon,(0,cnt),border=3,flag=flags);cnt+=1
        self.sizer3.Add(self.buttonRecon,(0,cnt),border=3,flag=flags);cnt+=1
        
        









        self.deconSizer.Add(self.sizer3, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.deconSizer.AddSpacer(4 if self.spectral_dim_count == 3 else 10)

    def controls_box(self, parent=None):
        flags = wx.ALIGN_LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM | wx.GROW
        buttonSize = (75,22)

        if parent is None:
            parent = self
        self.peakLbl = wx.StaticBox(parent, -1, 'Peak Tools:')  #, size=(240, 140))
        self.peakSizer = wx.StaticBoxSizer(self.peakLbl, wx.VERTICAL)

        # Rendering contract for Peak Tools:
        #   * changing/toggling a tool never redraws the spectrum;
        #   * select/deselect updates only existing peak/label artists;
        #   * Add/Remove/Move/MoveGrp/Maximise/Undo/Redo rebuild only the
        #     peak/label layer and preserve contours plus the current zoom;
        #   * only spectrum/contour/display changes call draw_figure();
        #   * each user action requests at most one canvas render.

        self.buttonSelect = PersistentStateButton(self.peakLbl, -1, "Select (s)", size=buttonSize)
        self.buttonSelect.Bind(wx.EVT_BUTTON, self.OnButtonSelect)
        self.buttonSelectGrp = PersistentStateButton(self.peakLbl, -1, "SelectGrp", size=buttonSize)
        self.buttonSelectGrp.Bind(wx.EVT_BUTTON, self.OnButtonSelectGrp)
        self.buttonRemove = wx.Button(self.peakLbl, -1, "Remove", size=buttonSize)
        self.buttonRemove.Bind(wx.EVT_BUTTON, self.OnButtonRemove)

        self.buttonAdjust = wx.Button(self.peakLbl, -1, "Maximise ( )", size=buttonSize)
        self.buttonAdjust.Bind(wx.EVT_BUTTON, self.OnButtonAdjust)
        self.buttonAdd = PersistentStateButton(self.peakLbl, -1, "Add (a)", size=buttonSize)
        self.buttonAdd.Bind(wx.EVT_BUTTON, self.OnButtonAdd)
        self.buttonUndo = wx.Button(self.peakLbl, -1, "Undo", size=buttonSize)
        self.buttonUndo.Bind(wx.EVT_BUTTON, self.OnButtonUndo)
        self.buttonRedo = wx.Button(self.peakLbl, -1, "Redo", size=buttonSize)
        self.buttonRedo.Bind(wx.EVT_BUTTON, self.OnButtonRedo)
        self.buttonMove = PersistentStateButton(self.peakLbl, -1, "Move (m)", size=buttonSize)
        self.buttonMove.Bind(wx.EVT_BUTTON, self.OnButtonMove)
        self.buttonMoveGrp = PersistentStateButton(self.peakLbl, -1, "MoveGrp", size=buttonSize)
        self.buttonMoveGrp.Bind(wx.EVT_BUTTON, self.OnButtonMoveGrp)
        self.buttonSelectAll = PersistentStateButton(self.peakLbl, -1, "SelectAll", size=buttonSize)
        self.buttonSelectAll.Bind(wx.EVT_BUTTON, self.OnButtonSelectAll)




        self.peakSizer.AddSpacer(10)
        border = 5
        # Keep the tools in a single vertical column so the attached panel is
        # narrow and the workflow reads naturally from top to bottom.
        self.sizer4 = wx.BoxSizer(wx.VERTICAL)

        # Arrange Peak Tools by workflow, with small visual gaps between
        # action groups.
        button_groups = (
            (self.buttonUndo, self.buttonRedo),
            (self.buttonSelect, self.buttonSelectGrp, self.buttonSelectAll),
            (self.buttonMove, self.buttonMoveGrp, self.buttonRemove, self.buttonAdjust),
            (self.buttonAdd,),
        )
        for group_index, group in enumerate(button_groups):
            if group_index:
                # Give each workflow group a clearer visual separation.
                self.sizer4.AddSpacer(10)
            for button in group:
                self.sizer4.Add(
                    button, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, border
                )

        self.peakSizer.Add(self.sizer4, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.peakSizer.AddSpacer(10)

        # Keep the controls as peakFrame attributes while their wx parent is
        # the attached right-hand Tools panel.
        if parent is self.toolsPanel:
            self.toolsDetachButton = wx.Button(self.toolsPanel, label='Detach')
            self.toolsDetachButton.Bind(wx.EVT_BUTTON, self.onToolsDetach)
            self.toolsWindowSizer = wx.BoxSizer(wx.VERTICAL)
            self.toolsWindowSizer.Add(self.peakSizer, 0, wx.EXPAND | wx.ALL, 8)
            self.toolsWindowSizer.AddStretchSpacer(1)
            # Keep Detach as the only window-management control in the Tools
            # panel; the panel itself is shown/hidden with the main Tools button.
            footer = wx.BoxSizer(wx.VERTICAL)
            footer.Add(self.toolsDetachButton, 0, wx.EXPAND)
            self.toolsWindowSizer.Add(footer, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
            self.toolsPanel.SetSizer(self.toolsWindowSizer)


    def peakfile_box(self, parent=None):
        flags = wx.ALIGN_LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM | wx.GROW


        if parent is None:
            parent = self
        self.peakfile_lbl = wx.StaticBox(parent, -1, 'Peak List:')  #, size=(240, 140))
        self.peakfile_sizer = wx.StaticBoxSizer(self.peakfile_lbl, wx.VERTICAL)

        self.peakfileLab=wx.StaticText(self.peakfile_lbl,-1,"PeakFile:")
        self.peakfileBox = wx.TextCtrl(self.peakfile_lbl,size=(150,22),style=wx.TE_PROCESS_ENTER)
        self.peakfileBox.SetValue(self._peak_list_path())

        self.mapfileLab=wx.StaticText(self.peakfile_lbl,-1,"Mapfile:")
        self.mapfileBox = wx.TextCtrl(self.peakfile_lbl,size=(150,22),style=wx.TE_PROCESS_ENTER)

        self.mapfile3DLab=wx.StaticText(self.peakfile_lbl,-1,"Mapfile3D:")
        self.mapfile3DBox = wx.TextCtrl(self.peakfile_lbl,size=(150,22),style=wx.TE_PROCESS_ENTER)



        try:
            self.mapfileBox.SetValue(self.peak_service.parameter('mapFile'))
        except:
            pass

        try:
            self.mapfile3DBox.SetValue(self.peak_service.parameter('mapFile3D'))
        except:
            pass


        
        self.buttonSave = wx.Button(self.peakfile_lbl, -1, "Save", size=(50,22))
        self.buttonSave.Bind(wx.EVT_BUTTON, self.OnButtonSave)
        self.buttonLoad = wx.Button(self.peakfile_lbl, -1, "Load", size=(50,22))
        self.buttonLoad.Bind(wx.EVT_BUTTON, self.OnButtonLoad)
        self.buttonLoadMap = wx.Button(self.peakfile_lbl, -1, "LoadMap", size=(50,22))
        self.buttonLoadMap.Bind(wx.EVT_BUTTON, self.OnButtonLoadMap)
        self.buttonMap = wx.Button(self.peakfile_lbl, -1, "Map", size=(50,22))
        self.buttonMap.Bind(wx.EVT_BUTTON, self.OnButtonMap)

        self.mapFileBtn = wx.Button(self.peakfile_lbl, label="...", size=(40,22))
        self.mapFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.peak_service.choose_file(evt, self.mapfileBox, full=True, save=self.OnButtonSaveDecon))        

        self.mapFile3DBtn = wx.Button(self.peakfile_lbl, label="...", size=(40,22))
        self.mapFile3DBtn.Bind(wx.EVT_BUTTON, lambda evt: self.peak_service.choose_file(evt, self.mapfile3DBox, full=True, save=self.OnButtonSaveDecon))        
        self.buttonMap3D = wx.Button(self.peakfile_lbl, -1, "Map", size=(50,22))
        self.buttonMap3D.Bind(wx.EVT_BUTTON, self.OnButtonMap3D)



        self.sizer = wx.GridBagSizer(3, 6);cnt=0
        self.sizer.Add((0,10), (cnt,0)); cnt+=1
        self.sizer.Add(self.peakfileLab,(cnt,0),border=3,flag=flags|wx.LEFT)
        self.sizer.Add(self.peakfileBox,(cnt,1),border=3,flag=flags)
        self.sizer.Add(self.buttonSave,(cnt,2),border=3,flag=flags)
        self.sizer.Add(self.buttonLoad,(cnt,3),border=3,flag=flags);cnt+=1

        self.sizer.Add(self.mapfileLab,(cnt,0),border=3,flag=flags|wx.LEFT)
        self.sizer.Add(self.mapfileBox,(cnt,1),border=3,flag=flags)
        self.sizer.Add(self.mapFileBtn,(cnt,2),border=3,flag=flags)
        self.sizer.Add(self.buttonLoadMap,(cnt,3),border=3,flag=flags);
        self.sizer.Add(self.buttonMap,(cnt,4),border=3,flag=flags); cnt+=1

        self.sizer.Add(self.mapfile3DLab,(cnt,0),border=3,flag=flags|wx.LEFT)
        self.sizer.Add(self.mapfile3DBox,(cnt,1),border=3,flag=flags)
        self.sizer.Add(self.mapFile3DBtn,(cnt,2),border=3,flag=flags)
        self.sizer.Add(self.buttonMap3D,(cnt,4),border=3,flag=flags); cnt+=1


        self.sizer.Add((0,10), (cnt,0))
        self.peakfile_sizer.Add(self.sizer, flag=wx.ALIGN_CENTER_HORIZONTAL)

        # When hosted by the modeless Map window, keep all widget references
        # on peakFrame but lay them out in the Map panel.
        if parent is self.mapPanel:
            self.mapCloseButton = wx.Button(self.mapPanel, label='Close')
            self.mapCloseButton.Bind(wx.EVT_BUTTON, self.onMapClose)
            self.mapWindowSizer = wx.BoxSizer(wx.VERTICAL)
            self.mapWindowSizer.Add(self.peakfile_sizer, 0, wx.EXPAND | wx.ALL, 8)
            self.mapWindowSizer.Add(self.mapCloseButton, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
            self.mapPanel.SetSizer(self.mapWindowSizer)
            self.mapWindowSizer.Fit(self.mapFrame)
            self.mapFrame.SetMinSize(self.mapFrame.GetSize())

    def _create_dimension_views(self):
        """Create only the auxiliary views required by this spectrum dimension."""
        if self.spectral_dim_count not in (2, 3, 4):
            raise ValueError('PeakFrame supports 2D, 3D, or 4D data (got %r)' % self.spectral_dim_count)
        self.fig_bore = None
        self.canvas_bore = None
        if self._has_1d_bore() or self.spectral_dim_count == 4:
            self.fig_bore = Figure(figsize=(5, 2))
            self.canvas_bore = FigCanvas(self, -1, self.fig_bore)
            bore_handler = self.draw_bores if self._has_1d_bore() else self.draw_bores_4d
            self.canvas.mpl_connect('motion_notify_event', bore_handler)

    def _layout_dimension_views(self):
        """Lay out the main spectrum view for the fixed data dimensionality.

        In 3D mode the 1D bore is deliberately kept out of this right-hand
        sizer: it lives at the bottom of the compact controls column instead.
        This lets the 2D projection use the full available height, with its
        matplotlib navigation toolbar directly underneath it.
        """
        self.vbox2.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        if self.spectral_dim_count == 4:
            self.vbox2.Add(self.canvas_bore, 2, wx.TOP | wx.BOTTOM | wx.LEFT | wx.GROW)
        self.vbox2.Add(self.toolbar, 0, wx.LEFT | wx.BOTTOM | wx.EXPAND)

    def _cache_coordinate_background(self, event=None):
        """Cache a clean main-canvas background for hover blitting.

        The coordinate label and hover crosshair are animated artists, so they
        are deliberately omitted from normal figure draws.  A single cached
        background can therefore be restored on every mouse move and both the
        crosshair and coordinate text can be painted in one blit.  This avoids
        two independent blitters restoring over one another.
        """
        try:
            self._hover_background = self.canvas.copy_from_bbox(self.fig.bbox)
        except Exception:
            self._hover_background = None

    def _update_plot_coordinates(self, event):
        """Blit the mouse crosshair and x/y readout as one atomic update."""
        label = getattr(self, '_plot_coordinate_text', None)
        vline = getattr(self, '_hover_vline', None)
        hline = getattr(self, '_hover_hline', None)
        if label is None or vline is None or hline is None:
            return

        inside = (event.inaxes is self.axes and
                  event.xdata is not None and event.ydata is not None)
        if inside:
            label.set_text('%.3f, %.3f' % (event.xdata, event.ydata))
            vline.set_xdata((event.xdata, event.xdata))
            hline.set_ydata((event.ydata, event.ydata))
            vline.set_visible(True)
            hline.set_visible(True)
        else:
            label.set_text('')
            vline.set_visible(False)
            hline.set_visible(False)

        background = getattr(self, '_hover_background', None)
        if background is None:
            # Normally only possible before the first draw/after an unusual
            # backend invalidation.  Let the ensuing draw_event rebuild it.
            self.canvas.draw_idle()
            return

        try:
            self.canvas.restore_region(background)
            if inside:
                self.axes.draw_artist(vline)
                self.axes.draw_artist(hline)
            self.fig.draw_artist(label)
            self.canvas.blit(self.fig.bbox)
        except Exception:
            self._hover_background = None
            self.canvas.draw_idle()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.selected_bore = None
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        self.canvas.mpl_connect('key_press_event', self.keyboard_press)
        self.canvas.mpl_connect('motion_notify_event', self._update_plot_coordinates)
        self.canvas.mpl_connect('draw_event', self._cache_coordinate_background)
        # Do not bind frame resize/idle events to Matplotlib redraws.  The wx
        # canvas handles its own resize, and incidental Tools-window layout
        # changes must never trigger a spectrum render.
        self._create_dimension_views()


        self.drawing_box()
        # A shared Reference list is always a 2D peak list, even when PeakFrame
        # is displaying a projection from 3D/4D data.  If that list is already
        # in memory, make the initial render include its peak markers.
        if bool(self.peak):
            self.cb_grid.SetValue(True)

        self.contour_box()
        self.peak_list_control_box()
        self.decon_box()
        self.controls_box(self.toolsPanel)
        self.peakfile_box(self.mapPanel)


        wx.StaticText(self, -1, '')

        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view, peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon, contour_callback=self._toolbar_contours, tools_callback=self._toolbar_tools, peaks_active=bool(self.cb_grid.GetValue()), coordinates=False)
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        self.vbox = wx.BoxSizer(wx.VERTICAL)

        # 3D needs a little more vertical room for the bore plot, so keep the
        # left-side controls deliberately tight.  Other dimensionalities retain
        # the established spacing.
        spacing = 8 if self._has_1d_bore() else 20
        self.vbox.AddSpacer(spacing)
        self.vbox.Add(self.drawSizer, 0, wx.EXPAND)
        self.vbox.AddSpacer(spacing)
        self.vbox.Add(self.peaklist_cntrl_box, 0, wx.EXPAND)
        self.vbox.AddSpacer(spacing)
        self.vbox.Add(self.deconSizer, 0, wx.EXPAND)
        self.vbox.AddSpacer(spacing)
        self.infoText = wx.StaticText(self, -1, "")
        self.vbox.Add(self.infoText, 0, wx.EXPAND)

        if self._has_1d_bore():
            # Put the 1D bore beneath the controls.  The canvas remains owned
            # by peakFrame; only its position in the sizer changes.
            self.canvas_bore.SetMinSize((300, 150))
            self.vbox.Add(self.canvas_bore, 1, wx.EXPAND | wx.TOP | wx.BOTTOM, spacing)

        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self._layout_dimension_views()

        self.hbox.AddSpacer(spacing)
        self.hbox.Add(self.vbox, 0, flag=wx.EXPAND | wx.TOP | wx.BOTTOM)
        self.hbox.AddSpacer(spacing)
        self.hbox.Add(self.vbox2, 1, flag=wx.EXPAND | wx.TOP | wx.BOTTOM)
        self.hbox.AddSpacer(spacing)
        # Peak Tools lives inside this frame, immediately to the right of the
        # spectrum.  It is initially hidden and therefore consumes no space.
        self.hbox.Add(self.toolsPanel, 0, flag=wx.EXPAND | wx.TOP | wx.BOTTOM)
        self.hbox.AddSpacer(spacing)
        self.SetSizerAndFit(self.hbox)


    def _set_status(self, message):
        """Show transient PeakFrame/tool feedback in the frame status bar."""
        if hasattr(self, 'statusBar') and self.statusBar is not None:
            self.statusBar.SetStatusText(str(message))

    #####################################################
    #button functions
    def OnButtonTranspose(self,event):
        print('Transposing peak list.')
        for pk in self.peak:
            ytmp=pk.y
            pk.y=pk.x
            pk.x=ytmp
            #print(self.peak[p].name,self.peak[p].x,self.peak[p].y)
        # for p in range(len(self.peak)):
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)
        self.DoAlias()
        self._commit_projection_peaks()
        self.draw_figure()

    def OnButtonOverlay(self,event):
        if(self.OVERLAY==1):
            print('Overlay is already loaded.')
            return

        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:

            path = str(dlg.GetPaths()[0])


        dlg.Destroy()

        if(os.path.exists(path)==0):

            return

        try:
            self.peak_service.cache_external_2d_view(path, namespace='peakframe_overlay')
        except Exception as exc:

            return
        self.overlay_path = path
        self.OVERLAY=1


        """
        for j in range(4):
            try:
                itm=self.sizer2.FindItemAtPosition((2,j))
                if(itm.IsWindow()):
                    self.sizer2.Detach(itm.GetWindow())
                    print 'detaching',2,j
            except:
                pass

        self.CleanDim() #remove old elements
        """

        self.textbox0_over = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox1_over = wx.TextCtrl(self,size=(35,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox2_over = wx.TextCtrl(self,size=(35,-1),style=wx.TE_PROCESS_ENTER)

        self.textboxX_over = wx.TextCtrl(self,size=(35,-1),style=wx.TE_PROCESS_ENTER)
        self.textboxY_over = wx.TextCtrl(self,size=(35,-1),style=wx.TE_PROCESS_ENTER)

        self.textbox0_over.SetValue(self.textbox0.GetValue())
        self.textbox1_over.SetValue(self.textbox1.GetValue())
        self.textbox2_over.SetValue(self.textbox2.GetValue())

        self.textboxX_over.SetValue(str(0))
        self.textboxY_over.SetValue(str(0))


        self.cb_overlay = wx.CheckBox(self, -1,"Overlay",style=wx.ALIGN_RIGHT)
        self.cb_overlay.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_overlay)
        """
        self.cln.append(self.textbox0_over)
        self.cln.append(self.textbox1_over)
        self.cln.append(self.textbox2_over)
        self.cln.append(self.cb_overlay)
        """
        
        self.sizer2.Add(self.textbox0_over, (2,0), border=3)
        self.sizer2.Add(self.textbox1_over, (2,1), border=3)
        self.sizer2.Add(self.textbox2_over, (2,2), border=3)
        self.sizer2.Add(self.textboxX_over, (2,3), border=3)
        self.sizer2.Add(self.textboxY_over, (2,4), border=3)
        self.sizer2.Add(self.cb_overlay, (2,5), border=3)



        #elf.SetSizerAndFit(self.sizer2)
        self.SetSizerAndFit(self.hbox)

        self.draw_figure()

    def CleanDim(self):
        try:
            a=len(self.cln)
        except:
            self.cln=[]
        
        clean=0
        for i in range(len(self.cln)):
            try:
                self.cln[i].Destroy()
                clean+=1
            except:
                pass
        
        if clean != len(self.cln):
            pass
        self.cln=[]


    def OnButtonDecon(self,event):
        target = self._analysis_spectrum_path()
        self.peakfileBox.SetValue(self._peak_list_path())
        # PeakFrame always analyses the displayed 2D target.  dimProj retains
        # the established line-shape mapping for 3D/4D projections; the
        # explicit target/dimension prevents parent topology (especially
        # 2D+pseudo) from turning this into a main-spectrum 3D job.
        dim_proj = 'jk' if self._parent_spectral_dim_count() > 2 else False
        self.peak_service.run_decon(dimProj=dim_proj, threshFac=self.threshFac,
                              caller='peakframe', input_override=target,
                              dimension_override=2,
                              projection_labels=self.view_labels)
            
        # Decon runs asynchronously. deconFrame consumes the generated list
        # on completion, and promotes 3D-project projection lists before
        # synchronising the shared reference peaks.
        #    self.OnButtonTranspose(True)
        self.draw_figure()
        
    def OnButtonRecon(self, event):

        # Recon uses the same 2D target as PeakFrame decon, but constrains the
        # allowed sources to the precise positions currently in this window.
        peak_path = self._peak_list_path()
        self.SavePeakList(peak_path)
        self.peakfileBox.SetValue(peak_path)
        target = self._analysis_spectrum_path()
        dim_proj = 'jk' if self._parent_spectral_dim_count() > 2 else False
        self.peak_service.run_decon(dimProj=dim_proj, caller='peakframe', recon=True,
                              input_override=target, dimension_override=2,
                              peak_list_override=peak_path,
                              projection_labels=self.view_labels)
        

        self._commit_projection_peaks()
        self.draw_figure()

        return

    def OnButtonSaveDecon(self,event):
        write={}
        write['mapFile']=self.mapfileBox.GetValue()
        write['mapFile3D']=self.mapfile3DBox.GetValue()

        target_path = self.peak_service.save_decon_parameters(write)
        print('Saving deconvolution parameters: %s' % target_path)
        if self.state is not None:
            self.state.sync_from_values(peak_file=self.peakfileBox.GetValue())

    #################################################################################
    #prototype function for mapping a decon result onto another specified peak list.
    #
    


    def OnButtonMap(self,event):
        #get dialog box here
        """
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a peaklist for mapping", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            mapFile = str(dlg.GetPaths()[0])
            print('Selected mapping file: %s' % mapFile)
            dlg.Destroy()
        else:
            print('Could not read the selected mapping file; operation cancelled.')
            dlg.Destroy()
            return
        """

        #3) get reference peaklist
        if(len(self.peak)==0):
            refFile=self.mapfileBox.GetValue() #original peaklistfile
            if(os.path.exists(refFile)==False):
                print('Mapping file not found: %s' % refFile)
                return
            print('Loading mapping file: %s' % self.mapfileBox.GetValue())
            self.OnButtonLoadMap(True)

        self.peakCopy=copy.deepcopy(self.peak) #backup peaklist
        print('Current peak list contains %i peaks.' % len(self.peak))


        #1) run decon on current spectrum.
        print('Running deconvolution to update the peak list...')
        self.threshFac=2.0
        self.OnButtonDecon(True)
        self.threshFac=1.0
        #2) read map peaklist
        self.peakDecon=copy.deepcopy(self.peak) #backup peaklist

        #print('Decon list:',len(self.peak))
        #for pk in self.peakDecon:
        #    pk.inty=self.GetIntensity(pk)
        #    #print(pk.name,pk.x,pk.y,pk.inty)

        
        #
        #for pk in self.peakDecon:
        #    print(pk.name,pk.x,pk.y)

        



        #now lets adjust the peak positions...
        #first lets get the intensity in the 2D at each position
        """
        inty=[]
        loc=[]
        for i,pk in enumerate(self.peak):
            loc.append(i)
            pk.inty=self.GetIntensity(pk)
            inty.append(pk.inty)
        inty=numpy.array(inty)
        argy=numpy.argmax(inty)
        pkMax=self.peak[loc[argy]]
        print('Aligning peak lists using peak: %s' % pkMax.name)
        refX=pkMax.x
        refY=pkMax.y
        """

        for pk in self.peakCopy:
            pk.inty=self.GetIntensity(pk)
            #print(pk.name,pk.x,pk.y,pk.inty)

        for pk in self.peakDecon:
            pk.inty=self.GetIntensity(pk)
            #print(pk.name,pk.x,pk.y,pk.inty)

            

        
        refX=numpy.zeros(len(self.peakDecon)*len(self.peakCopy))
        refY=numpy.zeros(len(self.peakDecon)*len(self.peakCopy))
        dint=numpy.zeros(len(self.peakDecon)*len(self.peakCopy))
        nam=numpy.zeros((len(self.peakDecon)*len(self.peakCopy)),dtype='str')
        cnt=0
        for pk in self.peakDecon: #for each peak in decon...
            for po in self.peakCopy:  #for each peak in mapping list....
                refX[cnt]=pk.x-po.x #reference peak in decon to peak in po
                refY[cnt]=pk.y-po.y #reference peak in decon to peak in po

                nam[cnt]=str(pk.name)+' (decon) to '+str(po.name)+' (map list)'
                #print(pk.name,po.name,nam[cnt])
                #for this referencing, work out the summed intensity
                pcopy=copy.deepcopy(self.peakCopy)
                for pi in pcopy:
                    pi.x+=refX[cnt]
                    pi.y+=refY[cnt]
                    self.DoAliasPk(pi)
                    dint[cnt]+=self.GetIntensity(pi)
                cnt+=1
        argy=numpy.argmax(numpy.fabs(dint))
        print('Aligning peak lists using maximum overlap: %s' % nam[argy])
        diffX=refX[argy]
        diffY=refY[argy]
                
        print('Peak-list alignment shift: %.3f, %.3f ppm' % (diffX, diffY))
        for pk in self.peakCopy:  #adjust positions by the shift in mapped list
            pk.x+=diffX
            pk.y+=diffY
            self.DoAliasPk(pk)

        peakCpy=copy.deepcopy(self.peakCopy)
        self.peak=copy.deepcopy(self.peakCopy)
        listy=[]
        res=[]

        return
        print('Locally maximising peak positions...')
        self.select=[]
        for p in range(len(self.peak)):
            self.select.append(p)
        self.OnButtonAdjust(True)
        self.select=[]


        #self.DoAlias()

        self._commit_projection_peaks()
        self.select=[]
        self.draw_figure()


    ##################################################################
    #Clustering 3Ds to get 2Ds.

    def OnButtonMap3D(self,event):

        fac=10
        posx,posy,posn=self.GetPositions(fac) #read out xyz and name
        if((posx).all()==False):
            return
        outfile=self.Cluster(posx,posy,posn,fac=fac)

        self.peak=self.read_peaklist_file(outfile)

        self._commit_projection_peaks()
        self.peak_service.refresh_status()
        self.draw_figure() #update figure.

    def GetPositions(self,fac):

        infile=self.mapfile3DBox.GetValue()
        if(os.path.exists(infile)==False):
            print('File not found: %s' % infile)
            return False,False,False
        
        class peakEntry3D():
            def __init__(self,test):
                self.name=test[0]
                self.f1=float(test[1])
                self.f2=float(test[2])
                self.f3=float(test[3])
                self.f3p=float(test[3])
                self.inty=float(test[4])

                if(len(self.name.split('_'))==1):
                    self.pk=self.name
                else:
                    self.pk=self.name.split('_')[0]
                    self.ind=self.name.split('_')[1]
                self.tp=''

        peak=[] #test this.
        read=False
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)==0):
                continue
            if(len(test)!=5):
                print("Peak list is not 3D. Re-run 3D decon.")
                break
                    
            
            read=True  #we are in business!
            entry=peakEntry3D(test)
            peak.append(entry)

        if(read==False): #break in peak list.
            return False,False,False

        posx=[]
        posy=[]
        posn=[]
        
        self.specsize=self.peak_service.data.shape  #spectrum shape.
        
        for ii,pk in enumerate(peak):

            self.peak_service.alias_peak(pk,pk.f3,0) #I
            self.peak_service.alias_peak(pk,pk.f2,1) #J
            self.peak_service.alias_peak(pk,pk.f1,2) #K

            k=pk.indexK   #findnear_index(pk.f1,self.index2) #find index2
            j=pk.indexJ   #findnear_index(pk.f2,self.index1) #find index1 
            i=pk.indexI   #findnear_index(pk.f3,self.index0) #find index0

            #inty=self.data[i,j,k]
            #print(pk.name,self.index0[i],self.index1[j],self.index2[k],inty,pk.f3,pk.f2,pk.f1,pk.inty)
            arr,steps=self.CheckMax(i,j,k) #get the local maximum.
            pk.indi=arr[0]
            pk.indj=arr[1]
            pk.indk=arr[2]
            pk.steps=steps
            #if(pk.pk not in indy.keys()):
            #    indy[pk.pk]=[]
            #indy[pk.pk].append(ii)

            #loc=self.index2[arr[2]],self.index1[arr[1]]
            posx.append(self.peak_service.physical_axis(2)[arr[2]])
            posy.append(self.peak_service.physical_axis(1)[arr[1]])
            posn.append(pk.pk)
            
            
        #set maxima.
        posx=numpy.array(posx) 
        posy=numpy.array(posy)
        posn=numpy.array(posn) #peak id with peak entry.

        #if two peaks have been put in the same position, and there is only one peak for each...
        #reset each to original position, rather than having both appear in the same place, which
        #breaks what follows.
        for i in range(len(posx)):
            dx=posx[i]-posx
            dy=posy[i]-posy
            dr=(dx**2+(dy/fac)**2)
            argy=numpy.where(dr<1E-6)[0]
            if(len(argy)==1): 
                continue
            #find if we have more than one that is arbitrarily close.
            #print("ARGY:",argy)
            #order the search so that we encounter the clash only once.
            test=[]
            for a in argy:
                if(a>=i):
                    test.append(a)
            if(len(test)<=1):
                continue
            #we have more than one arbitarily close, encountered once.
            test=numpy.array(test)
            #print('   aa',test)
            #print('   bb',posn[test])
            uni=numpy.unique(posn[test])  #unique position names.
            if(len(uni)<=1):
                continue
            #we have more than 1 unique name in the set. 
            #we have a clash risk.
            #print('   fff',uni)
            #lets find out how many peaks of these types we have. 
            sz=[]
            for u in uni:
                ar=numpy.where(posn==u)[0]
                sz.append(len(ar))
            sz=numpy.array(sz)
            #print ('   f',sz)
            oni=numpy.where(sz==1)[0]
            if(len(oni)<=1): #if we have two or more orphans, reset.
                continue
            #the unique peaks, all of them have only 1 peak.
            # this is the clash condition. resetting.
            for t in test: #do the reset.
                pk=peak[t]
                print('Peak adjustment did not converge; restoring peak: %s' % pk.name)
                posx[t]=self.peak_service.physical_axis(2)[pk.indexK]  #position already found.
                posy[t]=self.peak_service.physical_axis(1)[pk.indexJ]  #position already found.

        #step=numpy.array(step)
        return posx,posy,posn


    #if there are peak entries present in main peak list, that are not in
    #the 3D peak list, add them back in their current position.
    def CombineLists(self,xx,yy,nn):
        x=[]
        y=[]
        n=[]

        for i in range(len(xx)): #add the current picks.
            x.append(xx[i])
            y.append(yy[i])
            n.append(nn[i])
        for pk in self.peak:  #add anyone left over, not in the 3D list.
            #print(pk.name)
            if(pk.name not in nn):
                n.append(pk.name)
                x.append(pk.x)
                y.append(pk.y)
 
        return self.SortPeakList(x,y,n)

    #sort a list by the first number found in the peak ornament
    def SortPeakList(self,x,y,n):
        x=numpy.array(x)
        y=numpy.array(y)
        n=numpy.array(n)
        #resort peak list based on any numbers in peak.
        import re
        num=[] #sort by any number in the peak name.
        for pk in n:
            num.append(int(re.findall(r'[0-9]+', pk)[0]))
        num=numpy.array(num)
        argo=numpy.argsort(num)
        x=x[argo]
        y=y[argo]
        n=n[argo]
        return x,y,n

    #------------------------------------------------------#
    #probably need to move this into the unidecproject.
    #get peak positions and names. work out averages
    #find if anyone is nearer to someone else. if so, remap
    #repeat until no-one needs remapping.
    #fac is for 15N divisor for CSP.

    def Cluster(self,posx,posy,posn,fac=10):

        print('Clustering 3D peaks...')

        
        pks=numpy.unique(posn) #get the unique peak names.
        go=0
        while(1==1):
            aveX,aveY,aveN=self.GetAve(pks,posx,posy,posn,fac)  #get average positions based on posn
            close,finish,locca=self.GetMaps(posx,posy,posn,aveX,aveY,aveN,fac) #get remap options
            if(len(close)==0): #if there are no peaks closer to someone else than the mean.
                break
            arga=numpy.argmin(close)  #get the peak who is closest to someone else...
            #print('Moving:',self.peak[locca[arga]].name,'to',finish[arga],'dist:',close[arga])
            posn[locca[arga]]=finish[arga] #remap. go again.
            go+=1
            if(go==1000):
                break
    
        aveX,aveY,aveN=self.CombineLists(aveX,aveY,aveN)

        #write adjusted peak list as a .max file.
        outfile=self.peakfileBox.GetValue()+'.max'
        print('Saving peak list: %s' % outfile)
        outy=open(outfile,'w')
        for i in range(len(aveX)):
            x=aveX[i];y=aveY[i];n=aveN[i]
            outy.write('%s\t%f\t%f\n' % (n,y,x))  #write name, N,H   
        outy.close()
        print('Adjusted %i peaks; saved peak list: %s' % (go, outfile))

        return outfile
    

    #get the average x/y for each peak.
    #fac is used for CSP measure.
    def GetAve(self,pks,posx,posy,posn,fac):
        aveX=[]
        aveY=[]
        aveN=[]
        dist=[]
        loc=[]

        for pk in pks: #for each unique peak...
            arr=numpy.where(pk==posn)
            loc.append(arr)
            xx=posx[arr]  #1H
            yy=posy[arr]  #15N
            
            aX=numpy.average(numpy.average(xx))
            aY=numpy.average(numpy.average(yy))
            
            aveX.append(aX)
            aveY.append(aY)
            aveN.append(pk)
            dist.append(((xx-aX)**2.+((yy-aY)/fac)**2.)**0.5)
            #print(pk,aX,aY,step[arr],dist[-1])

        aveX=numpy.array(aveX)  #aveX
        aveY=numpy.array(aveY)  #aveY
        aveN=numpy.array(aveN)  #unique name.
        return aveX,aveY,aveN
       
    #go through all peaks, find ones closer to someone else than their own locus    
    def GetMaps(self,posx,posy,posn,aveX,aveY,aveN,fac):
        #loop again...
        close=[]
        finish=[]
        locca=[]
        for ii in range(len(posx)): #for each peak. Who are we nearest to?
            dx=posx[ii]-aveX
            dy=posy[ii]-aveY

            dr=(dx**2+(dy/fac)**2)
            argy=numpy.argmin(dr) #find nearest...
            if(posn[ii]!=aveN[argy]): #if this maps to someone different...
                #argo=numpy.where(aveN==posn[ii])  #flag t
                #dr[argo]**0.5 #self-distance, larger than new distance.
                #print(self.peak[ii].name,posn[ii],aveN[argy],dr[argy]**0.5)
                close.append(dr[argy])    #closest distnace...
                finish.append(aveN[argy]) #map to this peak.
                locca.append(ii)          #location of peak that is going to get mapped.
        return numpy.array(close),finish,locca


    #for given i and dimension, get the triplet i-1,i,i+1
    #if less than zero or above max, adjust.
    def GetTriplet(self,i,dim):
        iv=i-1,i,i+1
        iv=numpy.array(iv,dtype=numpy.int64)
        m1=iv<0
        iv[m1]+=self.specsize[dim] #roll by data size.
        m2=iv>=self.specsize[dim]
        iv[m2]-=self.specsize[dim]
        return iv

    #we will need to use these.
    ind3x3={}
    for ii in range(27):
        ind3x3[ii]=ii//9,ii%9//3,ii%3
    one3=1,1,1
    one3=numpy.array(one3)

    #for the given place, get the 9 steps around given ijk
    #find the max. If the max is in the centre, we are done!
    def Step3D(self,arr):
        #print('y',arr)
        iv=self.GetTriplet(arr[0],0)
        jv=self.GetTriplet(arr[1],1)
        kv=self.GetTriplet(arr[2],2)
        #inty=self.data[iv,jv,kv] #doesn't work.
        inty=self.peak_service.data[iv,:,:][:,jv,:][:,:,kv]  #not pretty? surely can do this better?
        
        #print(numpy.fabs(inty).shape)
        argy=numpy.argmax(numpy.fabs(inty)) #get max, but argmax unravels.
        ind=self.ind3x3[argy] #map the unraveled max to a 3x3 index. 
        
        if((ind==self.one3).all()): #if central spot is max...
            return arr,True
        arr=numpy.array((iv[ind[0]],jv[ind[1]],kv[ind[2]])) #update.
        return arr,False

        #verify the two are the same.
        #print(numpy.max(numpy.fabs(inty)))
        #print(self.data[iv[ix],jv[jx],kv[kx]])
        #return ind
    
    #take a point, ijk, and step to the max until the maximum sits in the centre of 
    #a 3x3 cube.
    def CheckMax(self,i,j,k):
        arr=(i,j,k)
        arr=numpy.array((i,j,k))
        #print('Starting intensity:',self.data[arr[0],arr[1],arr[2]])
        steps=0
        while(1==1):
            steps+=1
            #print(arr)
            arr,maxy=self.Step3D(arr)
            if(maxy):
                #print ("Success!")
                break
        #
        #print('Final intensity:',self.data[arr[0],arr[1],arr[2]])
        return arr,steps
    
    #END reclustering 3Ds to find 2D
    ##################################################################



    #####################################################
    def DoAlias(self,dimProj=False):
        print ('Doing alias')

        if(dimProj):
            for pk in self.peak:
                if self._parent_spectral_dim_count() == 4:
                    self.peak_service.alias_peak(pk,pk.x,2)
                    self.peak_service.alias_peak(pk,pk.y,3)
                elif self._parent_spectral_dim_count() == 3:
                    self.peak_service.alias_peak(pk,pk.x,2)
                    self.peak_service.alias_peak(pk,pk.y,1)
            self._commit_projection_peaks()
            return

        if(self._parent_spectral_dim_count()==2):
            if(self.state.pseudo_axis):
                for pk in self.peak:
                    self.peak_service.alias_peak(pk,pk.x,1)
                    self.peak_service.alias_peak(pk,pk.y,2)
                self._commit_projection_peaks()
                return
            else:
                for pk in self.peak:
                    self.peak_service.alias_peak(pk,pk.y,1)
                    self.peak_service.alias_peak(pk,pk.x,0)
                self._commit_projection_peaks()
                return
        if(self._parent_spectral_dim_count()==3):
            for pk in self.peak:
                self.peak_service.alias_peak(pk,pk.y,0) #I
                self.peak_service.alias_peak(pk,pk.y,1) #J
                self.peak_service.alias_peak(pk,pk.x,2) #K
            self._commit_projection_peaks()
            return
        if(self._parent_spectral_dim_count()==4):
            for pk in self.peak:
                self.peak_service.alias_peak(pk,pk.x,2) #K
                self.peak_service.alias_peak(pk,pk.y,3) #L
            self._commit_projection_peaks()
            return

    def DoAliasPk(self,pk):
        
        if(self._parent_spectral_dim_count()==2):
            if(self.state.pseudo_axis):
                self.peak_service.alias_peak(pk,pk.y,2)
                self.peak_service.alias_peak(pk,pk.x,1)
            else:
                self.peak_service.alias_peak(pk,pk.y,1)
                self.peak_service.alias_peak(pk,pk.x,0)
        if(self._parent_spectral_dim_count()==3):
            self.peak_service.alias_peak(pk,pk.y,0) #I
            self.peak_service.alias_peak(pk,pk.y,1) #J
            self.peak_service.alias_peak(pk,pk.x,2) #K
        if(self._parent_spectral_dim_count()==4):
            self.peak_service.alias_peak(pk,pk.x,2) #K
            self.peak_service.alias_peak(pk,pk.y,3) #L

    def GetLoc(self,pk):
        if(self._parent_spectral_dim_count()==4):
            return pk.indexK,pk.indexL
        if(self._parent_spectral_dim_count()>=3):
            return pk.indexJ,pk.indexK
        else:
            if(self._parent_spectral_dim_count()==2):
                if(self.state.pseudo_axis):
                    return pk.indexK,pk.indexJ
                else:
                    return pk.indexJ,pk.indexI

    def GetIntensity(self,pk):
        return self.ZZ[self.GetLoc(pk)]


            
    ##########################################################
    #local aliasing functions
    def OnButtonAliasim(self,event):
        if(len(self.select)==0):
            print('Aliasing cancelled: no peak is selected.')
            return
        p=self.select[0]
        bef=self.peak[p].x
        self.peak[p].x-=numpy.fabs(numpy.max(self.XX)-numpy.min(self.XX)+self.XX[0,1]-self.XX[0,0])

        self.DoAliasPk(self.peak[p])       
        #    if(self.state.pseudo_axis):
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].x,1)
        #    else:
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].x,1)
                
        print('Aliased peak %s: %.3f -> %.3f ppm' % (self.peak[p].name, bef, self.peak[p].x))
        self.ShowPeak(p)
        self._commit_projection_peaks()
        self.draw_figure()
        pass
        pass
    def OnButtonAliasip(self,event):
        if(len(self.select)==0):
            print('Aliasing cancelled: no peak is selected.')
            return
        p=self.select[0]
        bef=self.peak[p].x
        self.peak[p].x+=numpy.fabs(numpy.max(self.XX)-numpy.min(self.XX)+self.XX[0,1]-self.XX[0,0])

        self.DoAliasPk(self.peak[p])
        #    if(self.state.pseudo_axis):
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].x,1)
        #    else:
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].x,1)

        #self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)
        print('Aliased peak %s: %.3f -> %.3f ppm' % (self.peak[p].name, bef, self.peak[p].x))
        self.ShowPeak(p)
        self._commit_projection_peaks()
        self.draw_figure()
        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #    self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        #    self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)

        pass
        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #    self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        #    self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)

        pass
    def OnButtonAliasjp(self,event):
        if(len(self.select)==0):
            print('Aliasing cancelled: no peak is selected.')
            return
        p=self.select[0]
        bef=self.peak[p].y
        self.peak[p].y+=numpy.fabs(numpy.max(self.YY)-numpy.min(self.YY)+self.YY[1,0]-self.YY[0,0])

        self.DoAliasPk(self.peak[p])        
        #    if(self.state.pseudo_axis):
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].y,2)
        #    else:
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)

        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        print('Aliased peak %s: %.3f -> %.3f ppm' % (self.peak[p].name, bef, self.peak[p].y))
        self.ShowPeak(p)
        self._commit_projection_peaks()
        self.draw_figure()


        pass


        pass
    def OnButtonAliasjm(self,event):
        if(len(self.select)==0):
            print('Aliasing cancelled: no peak is selected.')
            return
        p=self.select[0]
        bef=self.peak[p].y
        self.peak[p].y-=numpy.fabs(numpy.max(self.YY)-numpy.min(self.YY)+self.YY[1,0]-self.YY[0,0])

        self.DoAliasPk(self.peak[p])        
        #    if(self.state.pseudo_axis):
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].y,2)
        #    else:
        #        self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)

        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        print('Aliased peak %s: %.3f -> %.3f ppm' % (self.peak[p].name, bef, self.peak[p].y))
        self.ShowPeak(p)
        self._commit_projection_peaks()
        self.draw_figure()
        pass
        pass



    #BIG DIFFERENCES TO CHARLIE
    def ReadPeakDecon(self):
        peak_path = self.peakfileBox.GetValue()
        if not os.path.exists(peak_path):
            alt = self.peak_service.resolve_spectrum_file(peak_path) if self.peak_service is not None else peak_path
            if os.path.exists(alt):
                peak_path = alt
        if not os.path.exists(peak_path):
            print('Peak list not found: %s' % peak_path)
            return

        print('Loading deconvolution peak list: %s' % peak_path)
        self.peak = self.read_peaklist_file(peak_path)

        if self.spectral_dim_count == 3:
            self.DoAlias(dimProj=True)
        else:
            self.DoAlias()

        self._commit_projection_peaks()

    ##################################################
    # Load a peak list
    def OnButtonLoad(self,event):
        peakListLocation=self.peakfileBox.GetValue()

        if(os.path.exists(peakListLocation)==False):
            print ('cannot find file:',peakListLocation)
            return

        print('Loading peak list: %s' % peakListLocation)

        self.peak=self.read_peaklist_file(peakListLocation)
        print('Peak list contains %i peaks.' % len(self.peak))
        if(self._parent_spectral_dim_count()==3):
            self.DoAlias(dimProj=True)
        else:
            self.DoAlias()
        # for p in range(len(self.peak)):
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)
        #self.tabOne.peak=self.peak 
        self.peak_service.refresh_status()
        self.draw_figure() #update figure.
        self._commit_projection_peaks()
        
        #self.Close()

    def OnButtonLoadMap(self,event):
        peakListLocation=self.mapfileBox.GetValue()

        if(os.path.exists(peakListLocation)==False):
            print ('cannot find file:',peakListLocation)
            return

        print('Loading peak list: %s' % peakListLocation)

        self.peak=self.read_peaklist_file(peakListLocation)

        if(self._parent_spectral_dim_count()==3):
            self.DoAlias(dimProj=True)
        else:
            self.DoAlias()            
        # for p in range(len(self.peak)):
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,0)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].y,1)
        #     self.peak_service.alias_peak(self.peak[p],self.peak[p].x,2)
        #self.tabOne.peak=self.peak 
        self.peak_service.refresh_status()
        self.draw_figure() #update figure.
        self._commit_projection_peaks()
        
        #self.Close()

        
        
    ##################################################
    # Save a peak list
    def OnButtonSave(self,event):
        outfile=self.peakfileBox.GetValue()
        self.SavePeakList(outfile)

        print('Peak list contains %i peaks.' % len(self.peak))
        #print self.dirBox.GetValue().split(self.parent.dirBox.GetValue())
        self._commit_projection_peaks()
        self.peak_service.refresh_status() if self.peak_service is not None else None

        
    def SavePeakList(self,outfile):
        print('Saving peak list: %s' % outfile)
        outy=open(outfile,'w')
        for pk in self.peak:
            outy.write('%s\t%f\t%f\n' % (pk.name, pk.y, pk.x))
        outy.close()

    #################################################
    # Create status bar.
    #################################################
    # select all peaks
    def _mode_buttons(self):
        return {
            'select': self.buttonSelect,
            'selectgrp': self.buttonSelectGrp,
            'add': self.buttonAdd,
            'move': self.buttonMove,
            'movegrp': self.buttonMoveGrp,
        }

    def _update_tool_controls(self):
        """Render tool state from the completed peak selection and mouse mode.

        A completed selection is deliberately different from a transient mouse
        interaction.  This keeps operation availability stable after Select,
        SelectGrp or SelectAll has finished.
        """
        buttons = self._mode_buttons()
        selected_button = None

        selection_count = len(set(self.select))
        has_selection = selection_count > 0
        all_selected = bool(self.peak) and selection_count == len(self.peak)
        single_selected = (
            selection_count == 1 and self.selection_type == 'single'
        )
        group_selected = (
            has_selection and self.selection_type == 'group'
        )
        selectall_selected = (
            all_selected and self.selection_type == 'all'
        )

        # Keep the selection status explicit for empty and multi-peak
        # selections.  A completed single Select keeps the more informative
        # peak-name message set by the click handler.
        if hasattr(self, 'infoText'):
            if not has_selection:
                self.infoText.SetLabel('Peaks selected: 0')
            elif group_selected or selectall_selected:
                self.infoText.SetLabel('Peaks selected: %i' % selection_count)

        if self.interaction_mode in buttons:
            selected_button = buttons[self.interaction_mode]
        elif single_selected:
            selected_button = self.buttonSelect
        elif group_selected:
            selected_button = self.buttonSelectGrp

        modal = self.interaction_mode not in (None, 'select')
        for name, button in buttons.items():
            active = button is selected_button
            button.SetActive(active)

            # Selection controls are available whenever no other mouse tool is
            # running. Add is deliberately the inverse of the selection-based
            # actions: it is available only when no peaks are selected. Move is
            # only a single-selection operation; MoveGrp is only a completed
            # SelectGrp/SelectAll operation.
            enabled = (not modal) or name == self.interaction_mode
            if name == 'add':
                enabled = enabled and not has_selection
            elif name == 'move':
                enabled = enabled and single_selected
            elif name == 'movegrp':
                enabled = enabled and (group_selected or selectall_selected)

            button.Enable(enabled)
            button.Refresh()

        self.buttonSelectAll.SetActive(selectall_selected and not modal)
        self.buttonSelectAll.Enable((not modal) and bool(self.peak))
        self.buttonSelectAll.Refresh()

        # These actions operate on any non-empty completed selection.
        self.buttonRemove.Enable((not modal) and has_selection)
        self.buttonAdjust.Enable((not modal) and has_selection)
        self._update_history_buttons()

    def _set_tool_mode(self, mode=None):
        """Set one mutually-exclusive transient mouse interaction mode."""
        self.interaction_mode = mode
        self.SELECT = int(mode == 'select')
        self.SELECTGRP = 1 if mode == 'selectgrp' else 0
        self.ADD = int(mode == 'add')
        self.MOVE = int(mode == 'move')
        self.MOVEGRP = 1 if mode == 'movegrp' else 0
        self._update_tool_controls()

    def _selection_complete(self, group=False):
        """Finish a selection and retain the appropriate persistent selector."""
        self.selection_type = 'group' if group else 'single'
        if group:
            self.interaction_mode = None
            self.SELECT = self.SELECTGRP = self.ADD = self.MOVE = self.MOVEGRP = 0
        else:
            # Single Select is intentionally persistent: every spectrum click
            # replaces the current selection with the nearest peak until the
            # Select button is pressed again (or another tool is chosen).
            self.interaction_mode = 'select'
            self.SELECT = 1
            self.SELECTGRP = self.ADD = self.MOVE = self.MOVEGRP = 0
        self._update_tool_controls()

    def _finish_tool_mode(self, clear_selection=False):
        self.interaction_mode = None
        self.SELECT = self.SELECTGRP = self.ADD = self.MOVE = self.MOVEGRP = 0
        if clear_selection:
            self.select=[]
            self.selection_type=None
        self._update_tool_controls()

    def _update_selection_artists(self):
        """Update peak/label selection with one canvas redraw and no contour rebuild."""
        if hasattr(self, 'buttonSelectAll'):
            self._update_tool_controls()
        if not hasattr(self, 'scatter_plot') or self.scatter_plot is None:
            return
        n=len(self.peak)
        if n:
            if self.select:
                # During a selection, de-emphasise unselected peaks and make
                # the selected X substantially longer while keeping its centre
                # exactly on the peak coordinate.  Matplotlib's X marker grows
                # symmetrically about the data point, so the intersection does
                # not move.
                colors=numpy.ones((n,4))*[0.0,0.0,0.0,0.3]
                colors[self.select]=[0.3,0.5,0.,1.0]
                sizes=numpy.full(n, 50.0)
                sizes[self.select]=140.0
            else:
                # Deselect restores the original peak-marker appearance rather
                # than leaving every peak in the grey de-emphasised state.
                colors=numpy.ones((n,4))*[0.0,0.0,0.0,1.0]
                sizes=numpy.full(n, 50.0)
            self.scatter_plot.set_color(colors)
            self.scatter_plot.set_sizes(sizes)
        for i,label in enumerate(getattr(self, 'labels', [])):
            label.set_color('g' if i in self.select else 'k')
        self.canvas.draw_idle()

    def _refresh_peak_artists(self):
        """Rebuild only peak markers/labels; retain contours, axes and current zoom."""
        old=getattr(self, 'scatter_plot', None)
        if old is not None:
            try: old.remove()
            except Exception: pass
        for label in getattr(self, 'labels', []):
            try: label.remove()
            except Exception: pass
        self.labels=[]
        plotting=[]
        for pk in self.peak:
            if self._parent_spectral_dim_count()==2:
                x,y=pk.x,pk.y
            else:
                x=pk.ppmK
                y=pk.ppmJ if self._parent_spectral_dim_count()==3 else pk.ppmL
            plotting.append([x,y])
            lab=pk.name if self.cb_labels.GetValue() else ''
            self.labels.append(self.axes.annotate(lab, xy=(x,y), xycoords='data',
                xytext=(3,3), textcoords='offset points', color='k', fontsize=12))
        plotting=ensure_xy_points(plotting)
        if self.select:
            colors=numpy.ones((len(self.peak),4))*[0.0,0.0,0.0,0.3]
            colors[self.select]=[0.3,0.5,0.,1.0]
            sizes=numpy.full(len(self.peak), 50.0)
            sizes[self.select]=140.0
        else:
            colors=numpy.ones((len(self.peak),4))*[0.0,0.0,0.0,1.0]
            sizes=numpy.full(len(self.peak), 50.0)
        self.scatter_plot=scatter_xy_points(self.axes, plotting, c=colors if len(self.peak) else 'k',
                                            s=sizes if len(self.peak) else 50, marker='x', zorder=2)
        visible=bool(self.cb_grid.GetValue())
        self.scatter_plot.set_visible(visible)
        for label in self.labels:
            label.set_visible(visible and bool(self.cb_labels.GetValue()))
        self._update_selection_artists()

    def _push_undo(self):
        self.undo_stack.append(copy.deepcopy(self.peak))
        self.redo_stack=[]
        self._update_history_buttons()

    def _update_history_buttons(self):
        if not hasattr(self, 'buttonUndo'):
            return
        mode_active = self.interaction_mode not in (None, 'select')
        self.buttonUndo.Enable(bool(self.undo_stack) and not mode_active)
        self.buttonRedo.Enable(bool(self.redo_stack) and not mode_active)

    def _next_peak_name(self):
        """Generate a conservative successor to the last list member's name."""
        if not self.peak:
            try:
                labels = self.peak_service.labels if self.peak_service is not None else self.view_labels
                return '1' + labels[-2][0] + '-' + labels[-1][0]
            except Exception:
                return '1'
        last = str(self.peak[-1].name)
        matches=list(re.finditer(r'\d+', last))
        if not matches:
            return last + '1'
        m=matches[-1]
        number=str(int(m.group(0))+1).zfill(len(m.group(0)))
        return last[:m.start()] + number + last[m.end():]

    def select_all(self):
        self.select=[]
        for p in range(len(self.peak)):
            self.select.append(p)


    #################################################
    # keyboard shortcuts
    def keyboard_press(self, event):
        #print(event.key)
        if event.key == 'q':
            if len(self.select) == 0:
                self.select_all()
                self.selection_type = 'all'
                self.interaction_mode = None
                self.SELECT = self.SELECTGRP = self.ADD = self.MOVE = self.MOVEGRP = 0
                self._update_selection_artists()
            else:
                self.OnButtonDeselect(event)
        if event.key == 's':
            if self.SELECT == 0 and not any((self.SELECTGRP, self.ADD, self.MOVE, self.MOVEGRP)):
                self._set_tool_mode('select')
            elif self.SELECT == 1:
                self.OnButtonDeselect(event)
        
            
        if event.key == 'backspace' or event.key == 'delete' or event.key == 'd':
            # print 'removing'

            self.OnButtonRemove(event)
        if event.key == 'm':
            if self.MOVE:
                self._finish_tool_mode()
            elif not any((self.SELECT, self.SELECTGRP, self.ADD, self.MOVEGRP)):
                if len(self.select) == 1:
                    self._set_tool_mode('move')
                    self._set_status('Peak tool: Move - click the destination position.')
                else:
                    print('Move requires exactly one selected peak.')
        if event.key == ' ':
            self.OnButtonAdjust(event)
        if event.key == 'a':
            if self.ADD == 0 and not any((self.SELECT, self.SELECTGRP, self.MOVE, self.MOVEGRP)):
                self._set_tool_mode('add')
            elif self.ADD == 1:
                self._finish_tool_mode()


    ########################################
    # Frame resize/idle events deliberately do not drive spectrum rendering.
    # FigureCanvasWxAgg handles canvas resizing itself; tying wx layout events
    # to canvas.draw() caused the first Peak Tools interaction to consume a
    # pending resize and visibly redraw the full spectrum.

    ######################################
    def draw_bores(self, event):
        x = event.xdata
        y = event.ydata
        if x is not None and y is not None:
            try:
                x1, y1 = self._nearest_projection_indices(x, y, context='hover')
                self.x1 = x1
                self.y1 = y1

                payload = self._bore_payload()
                if payload is None:
                    return
                _axis, cube, _label = payload
                trace = cube[:, y1, x1]
                self.bore.set_ydata(trace)
                if self.selected_bore != None:
                    bore_min = min(numpy.min(trace), self.selected_bore_min)
                    bore_max = max(numpy.max(trace), self.selected_bore_max)
                else:
                    bore_min = numpy.min(trace)
                    bore_max = numpy.max(trace)
                self.axes_bore.set_ylim(bore_min, bore_max)

                self.canvas_bore.draw_idle() #_artist(self.bore)
            except Exception as exc:
                raise
        
    def _set_selected_1d_bore(self, x, y):
        """Pin the bore trace at a selected 2D peak position."""
        if not self._has_1d_bore():
            return False
        payload = self._bore_payload()
        if payload is None:
            return False
        bore_axis, cube, _label = payload
        x1, y1 = self._nearest_projection_indices(x, y, context='selected peak')
        trace = cube[:, y1, x1]
        if self.selected_bore is None:
            self.selected_bore, = self.axes_bore.plot(bore_axis, trace, lw=0.5, color='darkgreen')
        else:
            self.selected_bore.set_xdata(bore_axis)
            self.selected_bore.set_ydata(trace)
        self.selected_bore_min = numpy.min(trace)
        self.selected_bore_max = numpy.max(trace)
        self.canvas_bore.draw_idle()
        return True

    def draw_bores_4d(self, event):
        x = event.xdata
        y = event.ydata
        

        if x is not None and y is not None:
            (x_min, x_max), (y_min, y_max) = self.peak_service.axis_limits[:2]
            if x > x_max:
                x = x-x_max+x_min
            if x < x_min:
                x = x-x_min+x_max
            if x < x_max and x > x_min and y < y_max and y > y_min:
                x1, y1 = numpy.abs(self.peak_service.physical_axis(0)-x).argmin(), numpy.abs(self.peak_service.physical_axis(1)-y).argmin()
                
                for tp in self.bore.collections:
                    tp.remove()
                self.bore = self.axes_bore.contour(self.peak_service.physical_axis(2), self.peak_service.physical_axis(3), self.peak_service.data[x1,y1,:,:].T, self.levels, cmap=cm.seismic,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels)))
                self.canvas_bore.draw_idle()


    #############################################
    # major plot function.
    def draw_figure(self,scale='y'):

        if len(self.peak) > 0:
            self.cb_grid.Enable()
            self.cb_labels.Enable()

        if(self.ax_reset==0):
            x_min,x_max=self.axes.get_xlim()
            y_min,y_max=self.axes.get_ylim()

        self.fig.clear()


        if self._has_1d_bore():
            self.fig_bore.clear()
            self.axes_bore = self.fig_bore.add_subplot(111)
            bore_axis, bore_cube, bore_label = self._bore_payload()
            yi = min(10, bore_cube.shape[1] - 1)
            xi = min(10, bore_cube.shape[2] - 1)
            initial_trace = bore_cube[:, yi, xi]
            self.bore, = self.axes_bore.plot(bore_axis, initial_trace, lw=0.5, color='r')
            self.axes_bore.set_xlim(bore_axis[0], bore_axis[-1])
            self.axes_bore.set_ylim(numpy.min(bore_cube), numpy.max(bore_cube))
            self.axes_bore.set_yticks([])
            #self.axes_bore.set_axis_off()
            self.axes_bore.get_xaxis().set_visible(True)

            # Keep the bore axis label visually consistent with its tick labels.
            bore_tick_size = self.axes_bore.xaxis.get_ticklabels()[0].get_fontsize() if self.axes_bore.xaxis.get_ticklabels() else 10
            self.axes_bore.tick_params(length=2.0, labelsize=bore_tick_size)
            bore_label = str(bore_label)
            real_labels = ('time_T2', 'ID', 'ncyc', 'ncyc_cp', 'gzlvl5', 'gzlvl1')
            if bore_label not in real_labels:
                bore_label += " (ppm)"
            self.axes_bore.set_xlabel(bore_label, fontsize=bore_tick_size, labelpad=2)
            self.axes_bore.spines['right'].set_linewidth(0.5) #.set_visible(False)
            self.axes_bore.spines['left'].set_linewidth(0.5)
            self.axes_bore.spines['bottom'].set_linewidth(0.5)
            self.axes_bore.spines['top'].set_linewidth(0.5)

            xmin, xmax = bore_axis[0], bore_axis[-1]
            ymin, ymax = self.axes_bore.get_yaxis().get_view_interval()
            #self.axes_bore.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))


            # Use nearly all of the compact bore canvas while retaining room
            # for tick labels and the x-axis title.
            self.fig_bore.subplots_adjust(left=0.035, right=0.995, bottom=0.22, top=0.98)
            self.canvas_bore.draw_idle()

        
        # self.fig_bore.subplots_adjust(left=0.1)
        

        self.levels=self.GetLevels()

        if self._parent_spectral_dim_count() == 4:
            self.fig_bore.clear()
            self.axes_bore = self.fig_bore.add_subplot(111)
            self.bore = self.axes_bore.contour(self.peak_service.physical_axis(2), self.peak_service.physical_axis(3), self.peak_service.data[10,10,:,:].T, self.levels, cmap=cm.seismic,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels)))
            y_min_bore, y_max_bore = self.axes_bore.get_ylim()
            x_min_bore, x_max_bore = self.axes_bore.get_xlim()
            self.axes_bore.set_xlim(x_max_bore,x_min_bore)
            self.axes_bore.set_ylim(y_max_bore,y_min_bore)
            self.axes_bore.set_xlabel(self.peak_service.labels[2],fontsize=8)
            self.axes_bore.set_ylabel(self.peak_service.labels[3],fontsize=8)
            self.fig_bore.tight_layout()
            self.canvas_bore.draw_idle()
        cnt=0 #for each combination of label, get nmrPipe projection

        self.axes = self.fig.add_subplot(111)
        # Mouse coordinates belong to the lower-left of the Matplotlib pane,
        # outside the graph axes themselves.  Figure coordinates keep the text
        # clear of the plotted data and of the graph border/ticks while still
        # leaving it inside the Matplotlib canvas.
        tick_fontsize = matplotlib.rcParams.get('xtick.labelsize', 8)
        self._plot_coordinate_text = self.fig.text(
            0.01, 0.01, '',
            ha='left', va='bottom', fontsize=tick_fontsize, zorder=2000,
            animated=True
        )

        if(self._parent_spectral_dim_count()==2):
            if(self.state.pseudo_axis):
                self.axes.set_xlabel(self.peak_service.labels[1],fontsize=8)
                self.axes.set_ylabel(self.peak_service.labels[2],fontsize=8)
            else:
                self.axes.set_xlabel(self.peak_service.labels[0],fontsize=8)
                self.axes.set_ylabel(self.peak_service.labels[1],fontsize=8)
        elif(self._parent_spectral_dim_count()==3):
            self.axes.set_xlabel(self.peak_service.labels[2],fontsize=8)
            self.axes.set_ylabel(self.peak_service.labels[1],fontsize=8)
        elif(self._parent_spectral_dim_count()==4):
            self.axes.set_xlabel(self.peak_service.labels[2],fontsize=8)
            self.axes.set_ylabel(self.peak_service.labels[3],fontsize=8)

        colormap=cm.Blues
        colormap2=cm.Reds
        colormap=cm.seismic
        if(self.OVERLAY): #if plotting an overlaid spectrum...
            if(self.cb_overlay.IsChecked()):
                print('Drawing peak overlay...')
                levelsOver=self.GetLevels(over=True)
                overlay = self.store.get_view(('peakframe_overlay', self.overlay_path, 'n')) if self.store is not None else None
                if overlay is not None:
                    self.axes.contour(
                        overlay['XX']+float(self.textboxX_over.GetValue()),
                        overlay['YY']+float(self.textboxY_over.GetValue()),
                        overlay['ZZ'], levelsOver, cmap=cm.Blues,
                        norm=colors.Normalize(vmin=-numpy.max(levelsOver),vmax=numpy.max(levelsOver))
                    )

        # For 2D data, XX/YY are already built as direct(X) x indirect(Y),
        # and ZZ is stored as [indirect, direct]. Do not transpose ZZ here.
        plot_data = self._as_ndarray(self.ZZ)
        self.axes.contour( self.XX, self.YY, plot_data,self.levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels))) #plot pdb network
        self.vline = self.axes.axvline(self.uc1.ppm_scale()[0], linewidth=0.5, color='k')
        self.hline = self.axes.axhline(self.uc0.ppm_scale()[0], linewidth=0.5, color='k')
        if(self.cb_calc.IsChecked()): #if showing calculation...
            colormap=cm.Blues
            try:
                Xs, Ys, Zs, _calc_payload = self._display_payload(decon=True)
            except RuntimeError:
                Xs = Ys = Zs = None
            if Zs is not None:
                self.axes.contour(
                    Xs, Ys, Zs, self.levels, cmap=colormap,
                    norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels))
                )
        
        
        # Keep the hover guides under PeakFrame's control rather than using
        # matplotlib.widgets.Cursor.  Cursor owns a separate axes-only blit
        # background; competing with the figure-level coordinate-text blit can
        # erase either the guides or the previous coordinate string.  These
        # animated artists are updated together by _update_plot_coordinates.
        self._hover_vline = self.axes.axvline(
            self.uc1.ppm_scale()[0], linewidth=0.75, color='k',
            zorder=1000, animated=True, visible=False
        )
        self._hover_hline = self.axes.axhline(
            self.uc0.ppm_scale()[0], linewidth=0.75, color='k',
            zorder=1000, animated=True, visible=False
        )
        for line in (self._hover_vline, self._hover_hline):
            line.set_solid_capstyle('projecting')
            line.set_snap(True)
        
        if(self.cb_grid.GetValue()==1): #if plotting peaks
            if(len(self.select)==0):
                sel=True
            else:
                sel=False
            if(self._parent_spectral_dim_count()==2):
                # if(sel):
                plotting = []
                self.labels = []
                for p,pk in enumerate(self.peak):

                    loc1=pk.x
                    loc2=pk.y

                    if(self.cb_labels.GetValue()):
                        lab=pk.name
                    else:
                        lab=''

                    
                    plotting.append([loc1, loc2])
                    # print(plotting)
                    self.labels.append(self.axes.annotate(lab,
                    xy=(loc1, loc2), xycoords='data',
                    xytext=(3, 3), textcoords='offset points', color='k', fontsize=12))
                plotting = ensure_xy_points(plotting)
                
                if len(self.select) == 0:
                    self.scatter_plot = scatter_xy_points(self.axes, plotting, c='k',s=50,marker='x',zorder=2)
                    
                else:
                    color_array = numpy.ones((len(self.peak),4))*[0.0,0.0,0.0,0.3]
                    color_array[self.select] = [0.3,0.5,0.,1.0]
                    marker_sizes = numpy.full(len(self.peak), 50.0)
                    marker_sizes[self.select] = 140.0
                    for selected_index in self.select:
                        self.labels[selected_index].set_color('g')
                    self.scatter_plot = scatter_xy_points(self.axes, plotting, c=color_array,s=marker_sizes,marker='x',zorder=2)


            else:
                
                plotting = []
                self.labels = []
                # selected = []
                for p,pk in enumerate(self.peak):
                    loc1=pk.ppmK
                    if self._parent_spectral_dim_count() == 3:
                        loc2= pk.ppmJ
                    elif self._parent_spectral_dim_count() == 4:
                        loc2 = pk.ppmL

                    lab=pk.name
                    
                    plotting.append([loc1, loc2])
    
                    self.labels.append(self.axes.annotate(lab,
                    xy=(loc1, loc2), xycoords='data',
                    xytext=(3, 3), textcoords='offset points', color='k', fontsize=12))

                plotting = ensure_xy_points(plotting)

                self.visible_label()

                if len(self.select) == 0:
                    self.scatter_plot = scatter_xy_points(self.axes, plotting, c='k',s=50,marker='x',zorder=2)
                    
                else:
                    color_array = numpy.ones((len(self.peak),4))*[0.0,0.0,0.0,0.3]
                    color_array[self.select] = [0.3,0.5,0.,1.0]
                    marker_sizes = numpy.full(len(self.peak), 50.0)
                    marker_sizes[self.select] = 140.0
                    for selected_index in self.select:
                        self.labels[selected_index].set_color('g')
                    self.scatter_plot = scatter_xy_points(self.axes, plotting, c=color_array,s=marker_sizes,marker='x',zorder=2)

                for sele in self.select:
                    pk=self.peak[sele]
                    loc1=pk.ppmK
                    
                    if self._parent_spectral_dim_count() == 3:
                        loc2= pk.ppmJ
                    elif self._parent_spectral_dim_count() == 4:
                        loc2 = pk.ppmL
                    lab=pk.name
                    x1, y1 = numpy.abs(self.uc1.ppm_scale()-loc1).argmin(), numpy.abs(self.uc0.ppm_scale()-loc2).argmin()
                    if self._parent_spectral_dim_count() == 3 and sele == self.select[0]:
                            selected_trace = self.peak_service.data[:, y1, x1]
                            self.selected_bore, = self.axes_bore.plot(self.peak_service.axis_scale(0), selected_trace, lw=0.5, color='darkgreen')
                            self.selected_bore_min = numpy.min(selected_trace)
                            self.selected_bore_max = numpy.max(selected_trace)
                    if self._parent_spectral_dim_count() == 4 and sele == self.select[0]:
                        # for line in self.selected_bore.collections:
                        
                        x1, y1 = numpy.abs(self.peak_service.physical_axis(0)-loc1).argmin(), numpy.abs(self.peak_service.physical_axis(1)-loc2).argmin()
                        self.selected_bore = self.axes_bore.contour(self.peak_service.physical_axis(2), self.peak_service.physical_axis(3), self.peak_service.data[x1,y1,:,:].T, self.levels, cmap=cm.Greens,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels)))
                

        #print(self.ax_reset)
        if(self.ax_reset==1): #if we want to reset the axis
            y_min=self.YY[0][0]
            y_max=self.YY[(len(self.YY))-1][0]
            x_min=self.XX[0][0]
            x_max=self.XX[0][(len(self.XX[0]))-1]
            xminy = min(x_min,x_max)
            xmaxy = max(x_min,x_max)
            yminy = min(y_min,y_max)
            ymaxy = max(y_min,y_max)
            self.axes.set_xlim(xmaxy,xminy)
            self.axes.set_ylim(ymaxy,yminy)
            self.ax_reset=0 #after a reset assume we have no further need to reset axes
        else:#otherwise use the last saved values
            # self.axes.set_xlim(x_min,x_max)
            # self.axes.set_ylim(y_min,y_max)
            xminy = min(x_min,x_max)
            xmaxy = max(x_min,x_max)
            yminy = min(y_min,y_max)
            ymaxy = max(y_min,y_max)
            self.axes.set_xlim(xmaxy,xminy)
            self.axes.set_ylim(ymaxy,yminy)

        # if(self.ax_reset==1): #if we want to reset the axis
        #     y_min=self.YY[0][0]
        #     y_max=self.YY[(len(self.YY))-1][0]
        #     x_min=self.XX[0][0]
        #     x_max=self.XX[0][(len(self.XX[0]))-1]
        #     xminy = min(x_min,x_max)
        #     xmaxy = max(x_max,x_max)
        #     yminy = min(y_min,y_max)
        #     ymaxy = max(y_max,y_max)
        #     self.axes.set_xlim(xmaxy,xminy)
        #     self.axes.set_ylim(ymaxy,yminy)
        #     self.ax_reset=0 #after a reset assume we have no further need to reset axes
        # else:#otherwise use the last saved values
        #     

        # PeakFrame margins: keep enough space on the left and bottom for
        # the ppm axis titles, while using most of the canvas at top/right.
        if self._parent_spectral_dim_count() == 3:
            self.fig.subplots_adjust(left=0.105, right=0.985, bottom=0.115, top=0.985)
        else:
            self.fig.subplots_adjust(left=0.105, right=0.990, bottom=0.120, top=0.990)
        # Full spectrum reconstruction is complete; render it exactly once.
        self.canvas.draw()

    def background_save(self, event):
        """Legacy compatibility hook; the peakFrame no longer caches a blit background.

        Older Projection code may still call this method.  Peak Tools no longer
        use restore_region()/blit(), so an incidental caller must not force a
        second spectrum redraw.
        """
        return None

    def visible_label(self):
        if(self.cb_labels.GetValue()==1):
            self.SetLabel_visible(True)
        else:
            self.SetLabel_visible(False)

    
    
    def SetLabel_visible(self, visible):
        for label_draw in self.labels:
            label_draw.set_visible(visible)

    ########################################
    # calculate contour levels
    def GetLevels(self,over=False):
        if(over==True):
            min_level=float(self.textbox0_over.GetValue())
            max_level=float(self.textbox1_over.GetValue())
            ctr_level=int(self.textbox2_over.GetValue())
        else:
            min_level=float(self.textbox0.GetValue())
            max_level=float(self.textbox1.GetValue())
            ctr_level=int(self.textbox2.GetValue())

        if(ctr_level==0):
            ctr_level=10
        if(max_level==0):
            max_level=1.2
        if(min_level==0):
            min_level=1E3

        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels


    #######################################
    # Projection/spectrum views are supplied by deconFrame/DataStore.

    def on_cb_grid(self, event):
        """Update peak/label overlays without changing the current viewport.

        Peaks and Labels are display-layer controls.  Historically they called
        ``draw_figure()`` with ``ax_reset = 1``, which rebuilt the spectrum and
        restored the default limits.  Rebuild only the lightweight peak artists
        instead so zoom/pan is retained for every PeakFrame data topology.
        """
        self.ax_reset = 0
        if hasattr(self, 'axes'):
            self._refresh_peak_artists()

    #def on_cb_grid_auto(self, event):
    #    self.draw_figure()

    def _toolbar_decon(self, active):
        self.cb_calc.SetValue(bool(active))
        self.draw_figure()

    def _toolbar_peaks(self, active):
        self.cb_grid.SetValue(bool(active))
        self.on_cb_grid(None)

    def _toolbar_contours(self):
        self.OnButtonContour(None)

    def redraw_view(self):
        self.ax_reset = 1
        self.draw_figure()

    def OnButtonDraw(self, event):
        self.redraw_view()

    #when search button is pressed make selection
    def on_pick(self, event):
        # Mouse coordinates are intentionally not printed; tool actions report their result.
        if(self.MOVE==1):
            if event.xdata is None or event.ydata is None:
                return
            if(len(self.select)!=1):
                print('Move cancelled: select exactly one peak.')
                self._finish_tool_mode()
                return
            self._push_undo()
            p=self.select[0]
            xnew=event.xdata
            ynew=event.ydata
            xind=numpy.argmin(numpy.fabs(self.XX[0,:]-xnew))
            yind=numpy.argmin(numpy.fabs(self.YY[:,0]-ynew))

            print('Moved peak %s to %.3f, %.3f ppm' % (self.peak[p].name, self.XX[yind,xind], self.YY[yind,xind]))

            #print self.peak[p].ppmK
            #print self.peak[p].ppmJ
            if(self._parent_spectral_dim_count()==2):
                self.peak[p].indexJ=xind
                self.peak[p].indexI=yind
                self.peak[p].ppmJ=self.XX[yind,xind]
                self.peak[p].ppmI=self.YY[yind,xind]
                self.peak[p].y=self.YY[yind,xind]
                self.peak[p].x=self.XX[yind,xind]

            elif(self._parent_spectral_dim_count()==3):
                self.peak[p].indexK=xind
                self.peak[p].indexJ=yind
                self.peak[p].ppmK=self.XX[yind,xind]
                self.peak[p].ppmJ=self.YY[yind,xind]
                self.peak[p].y=self.YY[yind,xind]
                self.peak[p].x=self.XX[yind,xind]

            elif(self._parent_spectral_dim_count()==4):
                self.peak[p].indexK=xind
                self.peak[p].indexL=yind
                self.peak[p].ppmK=self.XX[yind,xind]
                self.peak[p].ppmL=self.YY[yind,xind]
                self.peak[p].y=self.YY[yind,xind]
                self.peak[p].x=self.XX[yind,xind]


           
            self._commit_projection_peaks()
            self._finish_tool_mode()
            self._refresh_peak_artists()
            return

        if(self.SELECT==1):
            if event.xdata is None or event.ydata is None or not self.peak:
                return
            self.select=[]
            #self.select.append((event.xdata,event.ydata))

            x_min,x_max=self.axes.get_xlim()
            y_min,y_max=self.axes.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            for p in range(len(self.peak)):
                #print self.peak[p].name
                if(self._parent_spectral_dim_count()==3):
                    xval=self.peak[p].ppmK  #proton
                    yval=self.peak[p].ppmJ  #carbon
                elif(self._parent_spectral_dim_count()==4):
                    xval=self.peak[p].ppmK  #proton
                    yval=self.peak[p].ppmL  #carbon
                else:
                    xval=self.peak[p].x #proton
                    yval=self.peak[p].y  #carbon
                rad2=((xval-event.xdata)/xdist)**2.+((yval-event.ydata)/ydist)**2.
                raddy.append(rad2)
            raddy=numpy.array(raddy)

            maxy=numpy.argmin(raddy)
            self.infoText.SetLabel('Peak selected: '+str(self.peak[maxy].name))

            self.select.append(maxy)
            

            # self.ShowPeak(maxy)

            self._selection_complete(group=False)
            self._update_selection_artists()

            if self._is_pseudo3d_dataset() and len(self.select) > 0:
                pk = self.peak[self.select[0]]
                self._set_selected_1d_bore(pk.x, pk.y)

            if not self.selected_bore and len(self.select) > 0:
                    sele = self.select[0]
                    pk=self.peak[sele]
                    if self.spectral_dim_count == 2:
                        loc1=pk.y
                        loc2=pk.x
                    else:
                        loc1=pk.ppmK
                        if self._parent_spectral_dim_count() == 3:
                            loc2=pk.ppmJ
                        elif self._parent_spectral_dim_count() == 4:
                            loc2=pk.ppmL
                    lab=pk.name
                    x1, y1 = numpy.abs(self.uc1.ppm_scale()-loc1).argmin(), numpy.abs(self.uc0.ppm_scale()-loc2).argmin()
                    if self._parent_spectral_dim_count() == 3 and sele == self.select[0]:
                        selected_trace = self.peak_service.data[:, y1, x1]
                        self.selected_bore, = self.axes_bore.plot(self.peak_service.axis_scale(0), selected_trace, lw=0.5, color='darkgreen')
                        self.selected_bore_min = numpy.min(selected_trace)
                        self.selected_bore_max = numpy.max(selected_trace)
                    if self._parent_spectral_dim_count() == 4 and sele == self.select[0]:
                        
                        x1, y1 = numpy.abs(self.peak_service.physical_axis(0)-loc1).argmin(), numpy.abs(self.peak_service.physical_axis(1)-loc2).argmin()
                        self.selected_bore = self.axes_bore.contour(self.peak_service.physical_axis(2), self.peak_service.physical_axis(3), self.peak_service.data[x1,y1,:,:].T, self.levels, cmap=cm.Greens,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels)))
            elif len(self.select) > 0:
                    sele = self.select[0]
                    pk=self.peak[sele]
                    loc1=pk.ppmK
                    if self._parent_spectral_dim_count() == 3:
                        loc2=pk.ppmJ
                    elif self._parent_spectral_dim_count() == 4:
                        loc2=pk.ppmL
                    lab=pk.name
                    if self._parent_spectral_dim_count() == 3:
                        x1, y1 = numpy.abs(self.uc1.ppm_scale()-loc1).argmin(), numpy.abs(self.uc0.ppm_scale()-loc2).argmin()
                        trace = self.peak_service.data[:, y1, x1]
                        self.selected_bore.set_ydata(trace)
                        self.selected_bore_min = numpy.min(trace)
                        self.selected_bore_max = numpy.max(trace)
                    if self._parent_spectral_dim_count() == 4 and sele == self.select[0]:
                        for line in self.selected_bore.collections:
                            line.remove()
                        axis_min, axis_max = self.peak_service.axis_limits[0]
                        if loc1 > axis_max:
                            loc1 = loc1-axis_max+axis_min
                        if loc1 < axis_min:
                            loc1 = loc1-axis_min+axis_max
                        x1, y1 = numpy.abs(self.peak_service.physical_axis(0)-loc1).argmin(), numpy.abs(self.peak_service.physical_axis(1)-loc2).argmin()
                        self.selected_bore = self.axes_bore.contour(self.peak_service.physical_axis(2), self.peak_service.physical_axis(3), self.peak_service.data[x1,y1,:,:].T, self.levels, cmap=cm.Greens,norm=colors.Normalize(vmin=-numpy.max(self.levels),vmax=numpy.max(self.levels)))

            

            if self._parent_spectral_dim_count() >= 3:
                self.canvas_bore.draw()
            # self.draw_figure()
            #self.SELECT=0

        if(self.ADD==1):
            if event.xdata is None or event.ydata is None:
                return
            self._push_undo()
            self.select=[]
            self.selection_type=None
            newres=self._next_peak_name()
            try:
                newpeak=copy.deepcopy(self.peak[-1])
                newpeak.name=newres
                newpeak.y=event.ydata
                newpeak.x=event.xdata
                self.peak.append(newpeak)
            except Exception as e:
                #print(e)
                newpeak=peakEntry((newres,event.ydata,event.xdata))
                self.peak=[newpeak]
                self.cb_grid.SetValue(True)

            
            self.infoText.SetLabel('Peak Added at %.2fppm, %.2fppm' % (newpeak.x,newpeak.y))


            #for p in range(len(self.peak)):
            #    ytmp=self.peak[p].y
            #    self.peak[p].y=self.peak[p].x
            #    self.peak[p].x=ytmp
            ind=len(self.peak)-1

            self.DoAliasPk(self.peak[ind])
            # self.peak_service.alias_peak(self.peak[ind],self.peak[ind].y,0)
            # self.peak_service.alias_peak(self.peak[ind],self.peak[ind].y,1)
            # self.peak_service.alias_peak(self.peak[ind],self.peak[ind].x,2)

            self._commit_projection_peaks()
            self.peak_service.refresh_status() if self.peak_service is not None else None

            self._refresh_peak_artists()
        if(self.SELECTGRP==2):
            self.xgrp2=event.xdata
            self.ygrp2=event.ydata

            if(self.xgrp1>self.xgrp2):
                xmin=self.xgrp2
                xmax=self.xgrp1
            else:
                xmax=self.xgrp2
                xmin=self.xgrp1

            if(self.ygrp1>self.ygrp2):
                ymin=self.ygrp2
                ymax=self.ygrp1
            else:
                ymax=self.ygrp2
                ymin=self.ygrp1
            self.select=[]
            for p in range(len(self.peak)):
                #print self.peak[p].name
                if(self._parent_spectral_dim_count()==2):
                    xval=self.peak[p].x  #proton
                    yval=self.peak[p].y  #carbon
                else:
                    xval=self.peak[p].ppmK  #proton
                    if self._parent_spectral_dim_count() == 3:
                        yval=self.peak[p].ppmJ  #carbon
                    elif self._parent_spectral_dim_count() == 4:
                        yval=self.peak[p].ppmL  #carbon

                if(xval>xmin):
                    if(xval<xmax):
                        if(yval>ymin):
                            if(yval<ymax):
                                self.select.append(p)



            print('Peak group selected: %i peak(s).' % len(self.select))
            # 3 means a completed group selection.  The button stays blue;
            # pressing it again clears the group and leaves the mode.
            self._selection_complete(group=True)
            self._update_selection_artists()
            return
        if(self.SELECTGRP==1):
            self._set_status('Peak tool: Select group - click the opposite corner.')
            self.xgrp1=event.xdata
            self.ygrp1=event.ydata
            
            self.SELECTGRP+=1
            return

        if(self.MOVEGRP==2):
            self.xgrp2=event.xdata
            self.ygrp2=event.ydata
            xdif=self.xgrp2-self.xgrp1
            ydif=self.ygrp2-self.ygrp1
            self._push_undo()
            for sele in self.select:
                self.peak[sele].y+=ydif
                self.peak[sele].x+=xdif
                self.DoAliasPk(self.peak[sele])

            self._commit_projection_peaks()
            self._finish_tool_mode()
            self._refresh_peak_artists()
            return
        if(self.MOVEGRP==1):
            self._set_status('Peak tool: Move group - click the destination position.')
            self.xgrp1=event.xdata
            self.ygrp1=event.ydata
            self.MOVEGRP+=1
            return

        if(self.PAIR==2):
            for sele in self.select:
                self.pairsel.append(sele)

            
            n1=self.peak[self.pairsel[0]].name
            n2=self.peak[self.pairsel[1]].name



            abort=0
            if(n1[0]==n2[0]):
                resType=n1[0]
                if(resType=='V'):
                    at='G'
                elif(resType=='L'):
                    at='D'
                else:
                    print('Pairing cancelled: only leucine (L) or valine (V) peaks can be paired.')
                    abort=1


                n1nums=re.findall(r'[0-9]+',n1)
                n2nums=re.findall(r'[0-9]+',n2)

                n1num=n1nums[0] #take residue number index

                resnums=[]
                GetNew=0
                for i,pk in enumerate(self.peak):
                    if(i!=self.pairsel[0] and i!=self.pairsel[1]):
                        nums=re.findall(r'[0-9]+',pk.name)
                        resnums.append(int(nums[0]))
                        if(nums[0]==n1nums[0]):
                            print('Residue ID is already in use; choosing the next available ID.')
                            GetNew=1
                if(GetNew==1):
                    for i in range(1000):
                        if (i+1) not in resnums:
                            
                            n1num=str(i+1)
                            break

                if(abort==0):

                    ind1='1'
                    ind2='2'
                    if(len(n1nums)>=2): #if already a 1/2 index:
                        if(n1nums[1]=='2'):
                          ind1='2'
                          ind2='1'

                    n1new=resType+n1num+'C'+at+ind1+'-H'
                    n2new=resType+n1num+'C'+at+ind2+'-H'
                    print('Paired peaks: %s, %s -> %s, %s' % (n1, n2, n1new, n2new))
                    self.peak[self.pairsel[0]].name=n1new
                    self.peak[self.pairsel[1]].name=n2new

            else:
                print('Pairing cancelled: selected peaks are not the same residue type.')




            self.pairsel=[]
            self.select=[]
            self.PAIR=0
            self._commit_projection_peaks()
            self._refresh_peak_artists()


        if(self.PAIR==1):
            self._set_status('Peak tool: Pair - select the second residue.')
            self.pairsel=[]
            for sele in self.select:
                self.pairsel.append(sele)
            self.PAIR+=1
            self.select=[]
            self.SELECT=1

    # def update_scatter(self):
    #     self.

    def set_label_selection_color(self):
        for p,pk in enumerate(self.peak):
            if self.select[0] == p:
                self.labels[self.select[0]].set_color('g')
            else:
                self.labels[p].set_color('k')

    def redraw_labels(self):
        for label in self.labels:
            self.axes.draw_artist(label)

    def ShowPeak(self,p):
        print('Selected peak: %s' % self.peak[p].name)

    def OnButtonDeselect(self,event):
        self.select=[]
        self.PAIR=0
        self._finish_tool_mode(clear_selection=True)
        self._update_selection_artists()

    def OnButtonPair(self,event):
        self._set_status('Peak tool: Pair - select the first residue.')
        self.PAIR=1
        self._set_tool_mode('select')

    def OnButtonSelectAll(self,event):
        # SelectAll owns its own persistent state.  It is active exactly when
        # every peak is selected, and does not borrow SelectGrp's blue state.
        all_selected = bool(self.peak) and len(set(self.select)) == len(self.peak)
        if all_selected:
            self.select=[]
            self.selection_type=None
        else:
            self.select_all()
            self.selection_type='all'
        self.interaction_mode=None
        self.SELECT = self.SELECTGRP = self.ADD = self.MOVE = self.MOVEGRP = 0
        self._update_tool_controls()
        self._update_selection_artists()

    def OnButtonSelect(self,event):
        state = not bool(getattr(event.GetEventObject(), 'IsActive', lambda: False)())
        if not state:
            self.OnButtonDeselect(event)
            return
        self._set_tool_mode('select')
        self._set_status('Peak tool: Select - click the spectrum to select the nearest peak.')

    def OnButtonSelectGrp(self,event):
        state = not bool(getattr(event.GetEventObject(), 'IsActive', lambda: False)())
        if not state or self.SELECTGRP == 3:
            self.OnButtonDeselect(event)
            return
        self._set_tool_mode('selectgrp')
        self._set_status('Peak tool: Select group - click two opposite corners.')

    def OnButtonRemove(self,event):
        if not self.select:
            print('No peak is selected.')
            return
        self._push_undo()
        self.select=sorted(self.select,reverse=True)
        for sele in self.select:
            self._set_status('Removed peak: %s' % self.peak[sele].name)
            self.peak.pop(sele)
        self.select=[]
        self._commit_projection_peaks()
        self.selection_type=None
        self._update_tool_controls()
        self._refresh_peak_artists()

    def OnButtonAdd(self,event):
        state = not bool(getattr(event.GetEventObject(), 'IsActive', lambda: False)())
        if not state:
            self._finish_tool_mode()
        else:
            self._set_tool_mode('add')
        self._set_status('Peak tool: Add - %s' % ('click the spectrum to add a peak.' if state else 'off.'))

    def OnButtonAdjust(self,event):
        if(len(self.select)==0):
            print('No peak is selected.')
            return
        before=copy.deepcopy(self.peak)
        self._push_undo()
        for sele in self.select:
            print('Maximising peak intensity: %s' % self.peak[sele].name)
            self.Adjust(self.peak[sele])
        for i,pk in enumerate(self.peak):
            for j,pok in enumerate(self.peak):
                if(j>i and pk.x==pok.x and pk.y==pok.y):
                    print('Peak adjustment would overlap %s and %s; restoring their previous positions.' % (pk.name, pok.name))
                    self.peak[i]=copy.deepcopy(before[i])
                    self.peak[j]=copy.deepcopy(before[j])
        self._commit_projection_peaks()
        self._refresh_peak_artists()

    def OnButtonUndo(self,event):
        if not self.undo_stack:
            print('Nothing to undo.')
            return
        print('Undoing last peak action.')
        self.redo_stack.append(copy.deepcopy(self.peak))
        previous_selection=list(self.select)
        previous_selection_type=self.selection_type
        self.peak=copy.deepcopy(self.undo_stack.pop())
        # Undo changes peak data, not the user's selection intent. Retain every
        # selected index that still exists in the restored peak list.
        self.select=[i for i in previous_selection if 0 <= i < len(self.peak)]
        self.selection_type=previous_selection_type if self.select else None
        self._update_tool_controls()
        self._refresh_peak_artists()

    def OnButtonRedo(self,event):
        if not self.redo_stack:
            print('Nothing to redo.')
            return
        print('Redoing peak action.')
        self.undo_stack.append(copy.deepcopy(self.peak))
        previous_selection=list(self.select)
        previous_selection_type=self.selection_type
        self.peak=copy.deepcopy(self.redo_stack.pop())
        # Redo follows the same rule as Undo: keep the current selection where
        # those selected indices still exist in the redone peak list.
        self.select=[i for i in previous_selection if 0 <= i < len(self.peak)]
        self.selection_type=previous_selection_type if self.select else None
        self._update_tool_controls()
        self._refresh_peak_artists()

    #########################################################
    
    def Adjust(self,pk,verb='y'):
        loc1,loc2=self.GetLoc(pk)
        #def GetIntensity(self,pk):
        #return self.ZZ[self.GetLoc(pk)]

        ind,loc1,loc2=self.GetMax(loc1,loc2)

        if(ind==(1,1)):
            if(verb=='y'):
                print('Peak is already at a local maximum.')
            pk.x=self.XX[loc1,loc2]
            pk.y=self.YY[loc1,loc2]
            self.DoAliasPk(pk)

            return
        cnt=0
        while(1==1):

            ind,loc1,loc2=self.GetMax(loc1,loc2)
            if(ind==(1,1)):
                pk.x=self.XX[loc1,loc2]
                pk.y=self.YY[loc1,loc2]
                self.DoAliasPk(pk)
                break
            cnt+=1
            if(cnt==30):
                
                break
        #self.draw_figure()

    def OnButtonMove(self,event):
        state = not bool(getattr(event.GetEventObject(), 'IsActive', lambda: False)())
        if not state:
            self._finish_tool_mode()
            return
        if len(self.select) != 1:
            print('Move requires exactly one selected peak.')
            self._finish_tool_mode()
            return
        self._set_tool_mode('move')
        self._set_status('Peak tool: Move - click the destination position.')

    def OnButtonMoveGrp(self,event):
        state = not bool(getattr(event.GetEventObject(), 'IsActive', lambda: False)())
        if not state:
            self._finish_tool_mode()
            return
        if not self.select:
            print('Move group requires one or more selected peaks.')
            self._finish_tool_mode()
            return
        self._set_tool_mode('movegrp')
        self._set_status('Peak tool: Move group - click the start position.')


    def GetMax(self,loc1,loc2):
        #NEED TO WORK IN EDGE EFFECTS
        #look at the 9 adjacent positions, find the max.

        loc1min=loc1-1
        loc1max=loc1+2
        loc2min=loc2-1
        loc2max=loc2+2
        if(loc1==0):
            loc1min+=1
        if(loc1==self.ZZ.shape[1]):
            loc1max-=1
        if(loc2==0):
            loc2min+=1
        if(loc2==self.ZZ.shape[0]):
            loc2max-=1

        a=numpy.fabs(self.ZZ[loc1min:loc1max,loc2min:loc2max])

        ind = numpy.unravel_index(numpy.argmax(a, axis=None), a.shape)
        if(ind[0]==0):
            loc1=loc1min
        if(ind[0]==2):
            loc1=loc1max-1
        if(ind[1]==0):
            loc2=loc2min
        if(ind[1]==2):
            loc2=loc2max-1

        return ind,loc1,loc2

    def onPeakList(self, event):
        """Open the canonical Reference 2D Peak List viewer.

        PeakFrame no longer owns a second peak-list manager.  Delegate to the
        same datastore-backed viewer used by the main-window Reference 2D
        ``Show`` button so both entry points have identical behaviour.
        """
        return self.peak_service.open_reference_peak_list(event) if self.peak_service is not None else None

#######################################################
