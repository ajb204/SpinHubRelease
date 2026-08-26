#!/usr/bin/python
"""General raw-data viewer for one spectral dimension plus one pseudo axis."""
import os
import traceback
import wx
from spinDecon.gui.context import context_for, project_for, data_for
import numpy as np
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas

from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar
from spinDecon.domain.pseudo_axis import PseudoAxisTable, pseudo_axis_path, load_saved_column, save_selected_column
from spinDecon.domain.dimensions.viewer_contract import topology_for


def _p2d_lifecycle_debug(message):
    return


class Pseudo2D(wx.Panel):
    """Contour pseudo-2D inspector with a mouse-selected 1D trace."""
    def __init__(self, parent, tabOne):
        super().__init__(parent=parent, id=wx.ID_ANY)
        self.parent, self.tabOne = parent, tabOne
        self.app_context = context_for(tabOne, parent)
        self.state = project_for(tabOne, parent)
        self.topology = topology_for(tabOne)
        self.pseudo_service = self.app_context.pseudo if self.app_context is not None else None
        if self.pseudo_service is None:
            from spinDecon.analysis.pseudo_service import PseudoAxisService
            self.pseudo_service = PseudoAxisService(tabOne)
        self.data = np.asarray(self.pseudo_service.data)
        if self.data.ndim != 2:
            self.data = np.squeeze(self.data)
        if self.data.ndim != 2:
            raise RuntimeError('Pseudo2D requires a two-dimensional processed spectrum')
        self.x = self._spectral_axis()
        # Normalise orientation to [pseudo row, direct spectral point].
        if self.data.shape[1] != len(self.x) and self.data.shape[0] == len(self.x):
            self.data = self.data.T
        self.axis_table = PseudoAxisTable.load(pseudo_axis_path(self))
        self.axis_column = self.axis_table.default_column(load_saved_column(self, 'pseudo2DDisplayAxis'))
        self.row = 0
        self._motion_pending = None
        self._slice_blit_background = None
        self._slice_blit_busy = False
        self._pseudo2d_debug_motion_count = 0
        self._pseudo2d_debug_last_row = None
        self.selected_fit_ppm = None
        self.selected_fit_slices = []
        self.selected_fit_name = ''
        self._build_ui()
        self._set_axis(self.axis_column)
        self.draw_figure()

    def _spectral_axis(self):
        return self.pseudo_service.spectral_axis(self.data.shape)

    def _build_ui(self):
        self.fig = Figure(constrained_layout=False)
        # Keep the interactive contour on the left and reserve the right-hand
        # side for the complete selected-peak surface (ppm x pseudo axis x intensity).
        gs = self.fig.add_gridspec(1, 2, width_ratios=(1.25, 1.0), wspace=0.28)
        self.axes = self.fig.add_subplot(gs[0, 0])
        self.peak3d_axes = self.fig.add_subplot(gs[0, 1], projection='3d')
        self.trace_axes = self.axes.twinx()
        self.trace_axes.patch.set_visible(False)
        self.trace_axes.set_navigate(False)
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view,
                                               contour_callback=self.on_contours,
                                               coordinates=True)
        self.toolbar.Realize()
        self.axisChoice = wx.ComboBox(self.toolbar, choices=self.axis_table.data_columns,
                                      style=wx.CB_READONLY, size=(150, -1))
        if self.axis_column in self.axis_table.data_columns:
            self.axisChoice.SetStringSelection(self.axis_column)
        self.axisChoice.Bind(wx.EVT_COMBOBOX, self.on_axis_choice)
        self.toolbar.AddSeparator()
        self.toolbar.AddControl(wx.StaticText(self.toolbar, label='Pseudo axis:'))
        self.toolbar.AddControl(self.axisChoice)
        self.toolbar.AddSeparator()
        self.fittingToolButton = wx.Button(self.toolbar, -1, 'Fitting', size=(-1, 22))
        self.fittingToolButton.Bind(wx.EVT_BUTTON, self.show_fitting_window)
        self.toolbar.AddControl(self.fittingToolButton)
        self.analysisToolButton = wx.Button(self.toolbar, -1, 'Analysis', size=(-1, 22))
        self.analysisToolButton.Bind(wx.EVT_BUTTON, self.show_analysis_selector)
        self.toolbar.AddControl(self.analysisToolButton)
        try:
            self.toolbar.bind_control_status_help(self.fittingToolButton, 'Inspect pseudo2D restrained-fit results')
            self.toolbar.bind_control_status_help(self.analysisToolButton, 'Open pseudo2D analysis controls')
        except AttributeError:
            pass
        self.toolbar.Realize()
        s = wx.BoxSizer(wx.VERTICAL); s.Add(self.canvas, 1, wx.EXPAND); s.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(s)
        self._make_contour_frame()
        self._motion_cid = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self._motion_suspended_for_fitting = False
        self.Bind(wx.EVT_WINDOW_DESTROY, self._on_owner_destroy_debug)
        self.canvas.Bind(wx.EVT_WINDOW_DESTROY, self._on_canvas_destroy_debug)
        pass
        self.canvas.mpl_connect('draw_event', self._on_draw)
        self.canvas.mpl_connect('resize_event', self._invalidate_slice_blit)
        self.canvas.mpl_connect('figure_leave_event', self._on_contour_pointer_leave)
        self.axes.callbacks.connect('xlim_changed', self._on_xlim_changed)

    def _make_contour_frame(self):
        self.contourFrame = wx.Frame(self, title='Contours', style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        p = wx.Panel(self.contourFrame); s = wx.BoxSizer(wx.HORIZONTAL)
        default = max(float(np.nanmax(np.abs(self.data))) * 0.03, np.finfo(float).eps)
        self.contourMin = wx.TextCtrl(p, value=str(default), size=(100,22), style=wx.TE_PROCESS_ENTER)
        self.contourFactor = wx.TextCtrl(p, value='1.2', size=(60,22), style=wx.TE_PROCESS_ENTER)
        self.contourNumber = wx.TextCtrl(p, value='15', size=(60,22), style=wx.TE_PROCESS_ENTER)
        for label, ctrl in [('Min:', self.contourMin), ('Factor:', self.contourFactor), ('Number:', self.contourNumber)]:
            s.Add(wx.StaticText(p, label=label), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 4); s.Add(ctrl, 0, wx.ALL, 4)
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.redraw_view)
        p.SetSizerAndFit(s); fs=wx.BoxSizer(wx.VERTICAL); fs.Add(p,1,wx.EXPAND); self.contourFrame.SetSizerAndFit(fs)
        self.contourFrame.Bind(wx.EVT_CLOSE, lambda e: self.contourFrame.Hide())

    def _set_axis(self, column):
        self.axis_column = column
        try:
            self.y = np.asarray(self.axis_table.numeric_values(column), dtype=float)
            self.ylabels = None
        except Exception:
            self.y = np.arange(len(self.axis_table.rows), dtype=float)
            self.ylabels = [str(r.get(column, '')) for r in self.axis_table.rows]
        n = min(len(self.y), self.data.shape[0]); self.y = self.y[:n]; self.data = self.data[:n]

    def _levels(self):
        try: lo=float(self.contourMin.GetValue()); fac=float(self.contourFactor.GetValue()); n=max(1,int(self.contourNumber.GetValue()))
        except Exception: lo=max(np.nanmax(np.abs(self.data))*0.03, np.finfo(float).eps); fac=1.2; n=15
        return lo * np.power(fac, np.arange(n, dtype=float))

    def draw_figure(self, keepaxes=False):
        oldx = self.axes.get_xlim() if keepaxes and self.axes.has_data() else None
        self.axes.clear(); self.trace_axes.clear()
        levels=self._levels(); X,Y=np.meshgrid(self.x,self.y)
        self.axes.contour(X,Y,self.data,levels=levels,linewidths=0.6)
        self.axes.contour(X,Y,-self.data,levels=levels,linewidths=0.6)
        labels = self.pseudo_service.labels
        self.axes.set_xlabel(((labels[0] if labels else '') or 'Direct') + ' (ppm)')
        self.axes.set_ylabel(self.axis_column or 'Pseudo axis')
        if self.ylabels is not None and len(self.ylabels) <= 30:
            self.axes.set_yticks(self.y); self.axes.set_yticklabels(self.ylabels)
        if oldx is not None:
            lo, hi = sorted(oldx); self.axes.set_xlim(hi, lo)
        elif len(self.x):
            self.axes.set_xlim(float(np.nanmax(self.x)), float(np.nanmin(self.x)))
        self.trace_line, = self.trace_axes.plot(self.x, self.data[self.row], linewidth=1.0, animated=True)
        # The slice trace is a transient mouse-inspection aid.  Keep it hidden
        # until the pointer is actually over the contour axes.
        self.trace_line.set_visible(False)
        self.trace_axes.set_ylabel('Intensity')
        self.selection_line = self.axes.axhline(self.y[self.row], linewidth=0.8, animated=True)
        self.selection_line.set_visible(False)
        self._slice_trace_visible = False
        self.fitting_peak_line = None
        if self.selected_fit_ppm is not None:
            self.fitting_peak_line = self.axes.axvline(self.selected_fit_ppm, linewidth=1.0)
        self._draw_selected_peak_3d()
        self._slice_blit_background = None
        self._rescale_trace_axis(); self.canvas.draw_idle()

    def set_fitting_peak(self, ppm, slices, peak_name=''):
        """Show a fitting selection on the contour and as stacked 1D traces in 3D."""
        self.selected_fit_ppm = float(ppm)
        self.selected_fit_slices = list(slices or [])
        self.selected_fit_name = str(peak_name or '')
        if self.fitting_peak_line is not None:
            try: self.fitting_peak_line.remove()
            except Exception: pass
        self.fitting_peak_line = self.axes.axvline(self.selected_fit_ppm, linewidth=1.0)
        self._draw_selected_peak_3d()
        self._invalidate_slice_blit()
        self.canvas.draw_idle()

    def _draw_selected_peak_3d(self):
        ax = self.peak3d_axes
        ax.clear()
        slices = self.selected_fit_slices
        if not slices:
            ax.text2D(.5, .5, 'Select a peak in Fitting', ha='center', va='center', transform=ax.transAxes)
            ax.set_xlabel('ppm'); ax.set_ylabel(self.axis_column or 'Slice'); ax.set_zlabel('Intensity')
            return
        arrays = [np.asarray(a, dtype=float) for a in slices if np.asarray(a).ndim == 2 and np.asarray(a).shape[1] >= 2]
        if not arrays:
            return
        npts = min(a.shape[0] for a in arrays)
        nslices = min(len(arrays), len(self.y))
        arrays = arrays[:nslices]
        ppm = arrays[0][:npts, 0]
        z = np.vstack([a[:npts, 1] for a in arrays])
        yy = np.asarray(self.y[:nslices], dtype=float)
        # Display the pseudo2D selection as a stack of ordinary 1D NMR
        # traces.  Each spectrum remains at its real pseudo-axis Y value;
        # intensity is represented only by Z (no filled/interpolated surface).
        for yval, trace in zip(yy, z):
            ax.plot(ppm, np.full(ppm.shape, yval, dtype=float), trace, linewidth=0.8)
        ax.set_xlabel('ppm'); ax.set_ylabel(self.axis_column or 'Slice'); ax.set_zlabel('Intensity')
        if self.selected_fit_name: ax.set_title(self.selected_fit_name)
        if len(ppm): ax.set_xlim(float(np.nanmax(ppm)), float(np.nanmin(ppm)))

    def _rescale_trace_axis(self):
        xmin,xmax=self.axes.get_xlim(); lo,hi=sorted((xmin,xmax)); mask=(self.x>=lo)&(self.x<=hi)
        vals=self.data[:,mask] if np.any(mask) else self.data
        mn=float(np.nanmin(vals)); mx=float(np.nanmax(vals)); pad=(mx-mn)*0.05 or max(abs(mx),1.0)*0.05
        self.trace_axes.set_ylim(mn-pad,mx+pad)

    def _on_xlim_changed(self, axes):
        if hasattr(self, 'trace_axes'):
            self._invalidate_slice_blit()
            self._rescale_trace_axis()
            self.canvas.draw_idle()

    def _invalidate_slice_blit(self, event=None):
        self._slice_blit_background = None

    def _on_draw(self, event=None):
        """Cache the static contour/background image for fast slice motion.

        The trace axis is a ``twinx`` axis and therefore sits on top of the
        contour axis for Matplotlib hit testing.  The moving artists are
        marked animated so they are excluded from this cached background.
        """
        if self._slice_blit_busy or not hasattr(self, 'trace_line'):
            return
        try:
            self._slice_blit_background = self.canvas.copy_from_bbox(self.fig.bbox)
            pass
            self._blit_slice(debug_reason='draw_event')
        except Exception:
            self._slice_blit_background = None

    def _set_slice_trace_visible(self, visible):
        """Show the mouse-selected 1D slice only while over the contour."""
        visible = bool(visible)
        if getattr(self, '_slice_trace_visible', False) == visible:
            return False
        self._slice_trace_visible = visible
        if hasattr(self, 'trace_line'):
            self.trace_line.set_visible(visible)
        if hasattr(self, 'selection_line'):
            self.selection_line.set_visible(visible)
        return True

    def _on_contour_pointer_leave(self, event=None):
        if self._set_slice_trace_visible(False):
            if not self._blit_slice(debug_reason='leave'):
                self.canvas.draw_idle()

    def _blit_slice(self, debug_reason=None):
        if self._slice_blit_background is None:
            if debug_reason:
                pass
            return False
        try:
            self._slice_blit_busy = True
            self.canvas.restore_region(self._slice_blit_background)
            self.axes.draw_artist(self.selection_line)
            self.trace_axes.draw_artist(self.trace_line)
            self.canvas.blit(self.fig.bbox)
            if debug_reason:
                pass
            return True
        except Exception as exc:
            pass
            self._slice_blit_background = None
            return False
        finally:
            self._slice_blit_busy = False

    def on_motion(self, event):
        """Select and display the pseudo2D row underneath the mouse.

        Row selection belongs to the contour axis, not the overlaid intensity
        (``twinx``) axis.  Test the mouse against the contour axes' pixel bbox
        and transform its display Y coordinate back through ``axes.transData``.
        This remains correct regardless of which overlapping axes Matplotlib
        reports in ``event.inaxes``.
        """
        self._pseudo2d_debug_motion_count += 1
        n = self._pseudo2d_debug_motion_count
        inaxes_name = ('contour' if event.inaxes is self.axes else
                       'trace' if event.inaxes is self.trace_axes else
                       repr(event.inaxes))
        # Print the first few events and then periodically, so we can tell
        # whether WX/Matplotlib is delivering motion events at all without
        # flooding the terminal during normal mouse movement.
        verbose = n <= 12 or (n % 50 == 0)
        if verbose:
            pass

        if not len(self.y) or event.x is None or event.y is None:
            if verbose:
                pass
            return
        inside = bool(self.axes.bbox.contains(event.x, event.y))
        if verbose:
            pass
        if not inside:
            if self._set_slice_trace_visible(False):
                if not self._blit_slice(debug_reason='outside'):
                    self.canvas.draw_idle()
            return
        try:
            contour_xy = self.axes.transData.inverted().transform((event.x, event.y))
            pseudo_y = float(contour_xy[1])
        except (TypeError, ValueError, OverflowError) as exc:
            pass
            return
        if not np.isfinite(pseudo_y):
            pass
            return

        row = int(np.nanargmin(np.abs(self.y - pseudo_y)))
        visibility_changed = self._set_slice_trace_visible(True)
        if verbose or row != self.row:
            pass
        if row == self.row:
            if visibility_changed:
                if not self._blit_slice(debug_reason='enter'):
                    self.canvas.draw_idle()
            return

        old_row = self.row
        self.row = row
        trace = np.asarray(self.data[row, :])
        # Load the complete 1D spectral slice for this pseudo-axis row.
        self.trace_line.set_data(self.x, trace)
        self.selection_line.set_ydata((self.y[row], self.y[row]))
        pass

        # Only the trace and row marker move.  Keep the contours cached and
        # blit these animated artists for responsive mouse tracking.
        if not self._blit_slice(debug_reason='motion'):
            pass
            self.canvas.draw_idle()

    def on_axis_choice(self, event):
        self._set_axis(self.axisChoice.GetValue()); save_selected_column(self, self.axis_column, 'pseudo2DDisplayAxis'); self.row=min(self.row,len(self.y)-1); self.draw_figure(keepaxes=True)

    def on_contours(self, event=None): self.contourFrame.Show(); self.contourFrame.Raise()
    def _suspend_motion_for_fitting(self):
        pass
        cid = getattr(self, '_motion_cid', None)
        if cid is not None:
            try:
                self.canvas.mpl_disconnect(cid)
                pass
            except Exception as exc:
                pass
            self._motion_cid = None
        self._motion_suspended_for_fitting = True
        self._invalidate_slice_blit()

    def _resume_motion_after_fitting(self):
        pass
        if not getattr(self, '_motion_suspended_for_fitting', False):
            return
        try:
            if self and getattr(self, 'canvas', None) and self._motion_cid is None:
                self._motion_cid = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
                pass
        except Exception as exc:
            pass
        self._motion_suspended_for_fitting = False

    def _on_owner_destroy_debug(self, event):
        pass
        event.Skip()

    def _on_canvas_destroy_debug(self, event):
        pass
        event.Skip()

    def _report_fitting_frame(self):
        """Create/refresh the Fitting inspector without showing it."""
        frame = getattr(self, 'fittingFrame', None)
        if frame is None or not frame:
            self.fittingFrame = Pseudo2DFittingFrame(self)
            self.fittingFrame.Hide()
        else:
            self.fittingFrame.refresh_results()
        return self.fittingFrame

    def fitting_window_report_data(self):
        return self._report_fitting_frame().fitting_window_report_data()

    def export_fitting_report_figures(self, report_dir, units=None):
        return self._report_fitting_frame().export_fitting_report_figures(report_dir, units)

    def available_downstream_analyses(self):
        """Return analyses available for one spectral dimension plus pseudo axis."""
        return ['Diffusion', 'Decay']

    def selected_downstream_analysis(self):
        return self.pseudo_service.downstream_analysis

    def _notify_analysis_changed(self):
        self.pseudo_service.notify_changed()

    def save_downstream_analysis(self, selection):
        """Persist the selected pseudo2D analysis in the normal project state."""
        selection = str(selection or '').strip()
        if selection not in self.available_downstream_analyses():
            return False
        self.pseudo_service.set_downstream_analysis(selection)
        return True

    def _make_analysis_frame(self):
        frame = wx.Frame(self, title='Analysis',
                         style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(frame)
        row = wx.BoxSizer(wx.HORIZONTAL)
        self.AnalCombo = wx.ComboBox(panel, -1, choices=self.available_downstream_analyses(),
                                     size=(120, 22), style=wx.CB_READONLY)
        saved = self.selected_downstream_analysis()
        if saved and self.AnalCombo.FindString(saved) != wx.NOT_FOUND:
            self.AnalCombo.SetStringSelection(saved)
        elif self.AnalCombo.GetCount():
            self.AnalCombo.SetSelection(0)
        save_button = wx.Button(panel, -1, 'Save', size=(-1, 22))
        open_button = wx.Button(panel, -1, 'Open', size=(-1, 22))
        close_button = wx.Button(panel, -1, 'Close', size=(-1, 22))
        save_button.Bind(wx.EVT_BUTTON, self.OnSaveAnalysisButton)
        open_button.Bind(wx.EVT_BUTTON, self.OnOpenAnalysisButton)
        close_button.Bind(wx.EVT_BUTTON, lambda evt: frame.Hide())
        for widget in (self.AnalCombo, save_button, open_button, close_button):
            row.Add(widget, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        panel.SetSizerAndFit(row)
        frame.SetSizerAndFit(wx.BoxSizer(wx.VERTICAL))
        frame.GetSizer().Add(panel, 1, wx.EXPAND)
        frame.Fit()
        frame.Bind(wx.EVT_CLOSE, lambda evt: frame.Hide())
        return frame

    def show_analysis_selector(self, event=None):
        frame = getattr(self, 'analysisFrame', None)
        if frame is None or not frame:
            self.analysisFrame = self._make_analysis_frame()
        saved = self.selected_downstream_analysis()
        if saved and self.AnalCombo.FindString(saved) != wx.NOT_FOUND:
            self.AnalCombo.SetStringSelection(saved)
        self.analysisFrame.Show()
        self.analysisFrame.Raise()

    def OnSaveAnalysisButton(self, event=None):
        selection = self.AnalCombo.GetValue()
        if not self.save_downstream_analysis(selection):
            wx.MessageBox('Choose an analysis type first.', 'Analysis', wx.OK | wx.ICON_WARNING)
            return False
        return True

    def OnOpenAnalysisButton(self, event=None):
        selection = self.AnalCombo.GetValue()
        if not self.save_downstream_analysis(selection):
            wx.MessageBox('Choose an analysis type first.', 'Analysis', wx.OK | wx.ICON_WARNING)
            return False
        if selection == 'Diffusion':
            return self.open_diffusion_analysis()
        # uSTA and Decay are deliberately persisted now; their pseudo2D launch
        # paths can be attached here without changing the project-state format.
        wx.MessageBox('%s analysis is not yet connected for pseudo2D.' % selection,
                      'Analysis', wx.OK | wx.ICON_INFORMATION)
        return True

    def open_diffusion_analysis(self):
        """Launch the existing pseudo-dimensional diffusion panel modelessly."""
        from spinDecon.gui.workspaces import pseudo2d_diffusion
        frame = getattr(self, 'diffusionFrame', None)
        if frame is None or not frame:
            frame = wx.Frame(self, title='Diffusion Analysis',
                             style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
            # The legacy diffusion panel consults parent.tabOne as well as its
            # explicit tabOne argument, so expose the project controller on the
            # modeless host frame.
            frame.tabOne = self.tabOne
            panel = pseudo2d_diffusion.Pseudo2DDiffusion(frame, self.tabOne)
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(panel, 1, wx.EXPAND)
            frame.SetSizer(sizer)
            frame.SetSize((1100, 760))
            self.diffusionFrame = frame
        self.diffusionFrame.Show()
        self.diffusionFrame.Raise()
        return True

    def show_fitting_window(self, event=None):
        frame = getattr(self, 'fittingFrame', None)
        pass
        self._suspend_motion_for_fitting()
        try:
            if frame is None or not frame:
                pass
                self.fittingFrame = Pseudo2DFittingFrame(self)
            else:
                pass
                self.fittingFrame.refresh_results()
            pass
            self.fittingFrame.Show()
            pass
            self.fittingFrame.Raise()
            pass
        except Exception as exc:
            pass
            self._resume_motion_after_fitting()
            raise
    def redraw_view(self, event=None):
        """Toolbar Draw: restore the complete pseudo2D X/Y view."""
        # A redraw is deliberately different from contour/axis-choice refresh:
        # discard any pan/zoom limits and reconstruct both contour dimensions
        # from the complete data extents.
        self._invalidate_slice_blit()
        self.draw_figure(keepaxes=False)
    def onFocus(self, event):
        if event is not None: event.Skip()

class Pseudo2DFittingFrame(wx.Frame):
    """Modeless inspector for restrained pseudo2D fit results in ``fit/``.

    This is intentionally a reduced-dimensional analogue of Pseudo3D's
    Fitting palette: peak selection chooses a .dat/.out pair and the right-hand
    plot rasters through the 1D fitted pseudo slices.  Slice labels are ordinal
    by design; numerical pseudo-axis values remain in the result files only.
    """
    def __init__(self, owner):
        super().__init__(owner, title='Fitting',
                         style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        self.owner = owner
        self.tabOne = owner.tabOne  # legacy construction compatibility
        self.pseudo_service = owner.pseudo_service
        self.fit_dir = self._fit_dir()
        self.current_peak = None
        self.slices = []
        self.slice_index = 0
        panel = wx.Panel(self)

        self.peakList = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for col, label in enumerate(('Peak', 'Group', '%err', 'f01(ppm)', 'w1(Hz)', 'g1', 'Phase (deg)', 'wD/wA')):
            self.peakList.InsertColumn(col, label)
        self.peakList.SetColumnWidth(0, 95)
        self.peakList.SetColumnWidth(1, 60)
        for col in range(2, 8): self.peakList.SetColumnWidth(col, 85)
        self.peakList.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_peak_selected)

        self.fig = Figure(constrained_layout=True)
        self.axes = self.fig.add_subplot(111)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.nav = RedrawNavigationToolbar(self.canvas, self.redraw_view, coordinates=True)
        self.nav.Realize()

        self.sliceChoice = wx.ComboBox(panel, style=wx.CB_READONLY, size=(100, -1))
        self.prevButton = wx.Button(panel, label='-', size=(32, -1))
        self.nextButton = wx.Button(panel, label='+', size=(32, -1))
        self.diffBox = wx.CheckBox(panel, label='Difference')
        # Reviewing pseudo2D fits is explicit workflow evidence: viewing the
        # fitting palette alone must not complete the workflow stage.
        self.correctButton = wx.Button(panel, label='Mark intensities as correct')
        self.correctButton.SetToolTip('Accept the displayed pseudo2D fitting results and save this review in the project.')
        self.prevButton.Bind(wx.EVT_BUTTON, lambda evt: self.step_slice(-1))
        self.nextButton.Bind(wx.EVT_BUTTON, lambda evt: self.step_slice(1))
        self.sliceChoice.Bind(wx.EVT_COMBOBOX, self.on_slice_choice)
        self.diffBox.Bind(wx.EVT_CHECKBOX, self.on_display_change)
        self.correctButton.Bind(wx.EVT_BUTTON, self.on_mark_intensities_correct)

        controls = wx.BoxSizer(wx.HORIZONTAL)
        controls.Add(wx.StaticText(panel, label='Slice:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        controls.Add(self.sliceChoice, 0, wx.RIGHT, 4)
        controls.Add(self.prevButton, 0, wx.RIGHT, 2)
        controls.Add(self.nextButton, 0, wx.RIGHT, 12)
        controls.Add(self.diffBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        controls.AddStretchSpacer(1)
        controls.Add(self.correctButton, 0, wx.ALIGN_CENTER_VERTICAL)

        plot = wx.BoxSizer(wx.VERTICAL)
        plot.Add(controls, 0, wx.EXPAND | wx.BOTTOM, 5)
        plot.Add(self.canvas, 1, wx.EXPAND)
        plot.Add(self.nav, 0, wx.EXPAND)
        main = wx.BoxSizer(wx.HORIZONTAL)
        main.Add(self.peakList, 2, wx.EXPAND | wx.ALL, 8)
        main.Add(plot, 5, wx.EXPAND | wx.TOP | wx.RIGHT | wx.BOTTOM, 8)
        panel.SetSizer(main)
        self.SetMinSize((850, 450)); self.SetSize((1100, 600))
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.Bind(wx.EVT_ACTIVATE, self._on_activate_debug)
        self.Bind(wx.EVT_SHOW, self._on_show_debug)
        self.Bind(wx.EVT_WINDOW_DESTROY, self._on_destroy_debug)
        self.canvas.Bind(wx.EVT_WINDOW_DESTROY, self._on_canvas_destroy_debug)
        pass
        self.refresh_results()

    def _fit_dir(self):
        return self.pseudo_service.fit_dir()

    @staticmethod
    def _natural_key(value):
        try: return (0, int(value))
        except (TypeError, ValueError): return (1, str(value))

    def _peak_names(self):
        if not os.path.isdir(self.fit_dir): return []
        names = {os.path.splitext(f)[0] for f in os.listdir(self.fit_dir) if f.endswith('.dat')}
        return sorted(names, key=self._natural_key)

    def fitting_files_complete(self):
        """Whether every Full 1D peak has the .dat/.out pair shown here."""
        payload = self.pseudo_service.full_peak_payload()
        expected = []
        for peak in payload.get('peaks', []) or payload.get('records', []) or []:
            name = getattr(peak, 'name', None)
            if name is None and isinstance(peak, dict): name = peak.get('name') or peak.get('Name')
            if name is not None and str(name).strip(): expected.append(str(name).strip())
        if not expected:
            for row in payload.get('rows', []) or []:
                if row: expected.append(str(row[0]).strip())
        if not expected:
            path = self.pseudo_service.full_peak_file()
            try:
                with open(path) as handle:
                    for line in handle:
                        fields = line.split()
                        if len(fields) < 2: continue
                        try: float(fields[1])
                        except ValueError: continue
                        expected.append(fields[0])
            except OSError:
                return False
        return bool(expected) and all(os.path.isfile(os.path.join(self.fit_dir, name + ext))
                                      for name in expected for ext in ('.dat', '.out'))

    def _parse_out(self, peak):
        vals = {'group':'', 'f01':'', 'w1':'', 'g1':'', 'dw1':'', 'd_over_a':'', 'phase':'', 'width_ratio':'', 'intensity':[]}
        path = os.path.join(self.fit_dir, peak + '.out')
        try:
            with open(path) as handle:
                for line in handle:
                    fields = line.split()
                    if fields and fields[0] == '#': fields = fields[1:]
                    if len(fields) < 2: continue
                    key = fields[0].split('(')[0]
                    if key == 'Overlap_group' and len(fields) >= 2:
                        vals['group'] = fields[1]
                    elif key in ('f01', 'w1', 'g1', 'dw1'):
                        try: vals[key] = '%.3f' % float(fields[1])
                        except ValueError: vals[key] = fields[1]
                    elif key == 'D/A':
                        try: vals['d_over_a'] = float(fields[1])
                        except ValueError: pass
                    elif not line.lstrip().startswith('#') and len(fields) >= 3:
                        try: vals['intensity'].append(float(fields[1]))
                        except ValueError: pass
        except OSError: pass
        try:
            ratio = float(vals['d_over_a'])
            vals['phase'] = '%.2f' % np.degrees(np.arctan(ratio))
        except (TypeError, ValueError):
            vals['phase'] = ''
        try:
            aw, dw = float(vals['w1']), float(vals['dw1'])
            vals['width_ratio'] = '%.2f' % (dw / aw) if aw else ''
        except (TypeError, ValueError):
            vals['width_ratio'] = ''
        return vals

    def _parse_dat(self, peak):
        """Parse C++ Protocol2PFit output into one [x,data,calc] array per slice."""
        path = os.path.join(self.fit_dir, peak + '.dat')
        slices, current = [], []
        try:
            with open(path) as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) == 4 and not fields[0].startswith('#'):
                        try: current.append((float(fields[1]), float(fields[2]), float(fields[3])))
                        except ValueError: continue
                    elif not fields and current:
                        slices.append(np.asarray(current, dtype=float)); current = []
            if current: slices.append(np.asarray(current, dtype=float))
        except OSError:
            return []
        return [arr for arr in slices if arr.size]

    def _percent_error(self, slices):
        if not slices: return ''
        data = np.concatenate([a[:,1] for a in slices]); calc = np.concatenate([a[:,2] for a in slices])
        denom = float(np.nanmax(np.abs(data))) if data.size else 0.0
        return ('%.3f' % (np.sqrt(np.nanmean((data-calc)**2)) / denom * 100.0)) if denom else ''

    def refresh_results(self):
        self._sync_review_button()
        previous = self.current_peak
        self.peakList.DeleteAllItems()
        names = self._peak_names()
        select = 0
        for row, peak in enumerate(names):
            slices = self._parse_dat(peak); vals = self._parse_out(peak)
            idx = self.peakList.InsertItem(row, peak)
            for col, value in enumerate((vals['group'], self._percent_error(slices), vals['f01'], vals['w1'], vals['g1'], vals['phase'], vals['width_ratio']), start=1):
                self.peakList.SetItem(idx, col, value)
            if peak == previous: select = row
        if names:
            self.peakList.Select(select); self.peakList.Focus(select)
        else:
            self.axes.clear(); self.axes.text(.5,.5,'No pseudo2D fit results found',ha='center',va='center',transform=self.axes.transAxes)
            self.canvas.draw_idle()

    def fitting_window_report_data(self):
        """Return the exact rows/groups displayed by the pseudo2D Fitting window."""
        self.refresh_results()
        columns = [self.peakList.GetColumn(col).GetText() for col in range(self.peakList.GetColumnCount())]
        rows = [[self.peakList.GetItem(row, col).GetText() for col in range(len(columns))]
                for row in range(self.peakList.GetItemCount())]
        grouped, singles, order = {}, [], []
        for row in rows:
            peak, group = row[0], row[1]
            if group and group not in ('-', '-1'):
                if group not in grouped:
                    grouped[group] = []; order.append(group)
                grouped[group].append(peak)
            else:
                singles.append({'group': None, 'peaks': [peak]})
        return columns, rows, ([{'group': group, 'peaks': grouped[group]} for group in order] + singles)

    def export_fitting_report_figures(self, report_dir, units=None):
        """Export one pseudo2D overview and one fit plot per slice for each group.

        Report plots deliberately differ from the interactive Fitting window:
        they use the complete x range stored in the .dat files and annotate all
        fitted members of an overlap group.  The pseudo axis is temporarily
        fixed to ``exp_no`` and the fitting palette is never shown.
        """
        report_dir = str(report_dir); os.makedirs(report_dir, exist_ok=True)
        if units is None: _columns, _rows, units = self.fitting_window_report_data()
        old_column = getattr(self.owner, 'axis_column', None)
        old_peak, old_slice = self.current_peak, self.slice_index
        paths = []
        try:
            if 'exp_no' in getattr(self.owner.axis_table, 'headers', []): self.owner._set_axis('exp_no')
            for number, source_unit in enumerate(units, start=1):
                unit = dict(source_unit); peaks = list(unit.get('peaks', []))
                unit['overview_figures'] = []; unit['slice_figures'] = []
                parsed = {peak: self._parse_dat(peak) for peak in peaks}
                fits = {peak: self._parse_out(peak) for peak in peaks}
                # The GUI 3D panel displays one selected fit at a time.  Use the
                # first member as the group's representative, but save only one
                # contour/3D canvas for the group.
                representative = peaks[0] if peaks else None
                if representative:
                    try:
                        self.owner.set_fitting_peak(float(fits[representative].get('f01')), parsed[representative], representative)
                        self.owner.draw_figure(keepaxes=False)
                        overview = 'pseudo2d_fit_%03d_group.pdf' % number
                        self.owner.canvas.print_figure(os.path.join(report_dir, overview), bbox_inches='tight')
                        unit['overview_figures'].append((representative, overview))
                    except (TypeError, ValueError):
                        overview = ''
                else: overview = ''

                nslices = max([len(v) for v in parsed.values()] or [0])
                group_slices = []
                for slice_index in range(nslices):
                    available = [(peak, arrs[slice_index]) for peak, arrs in parsed.items() if slice_index < len(arrs)]
                    if not available: continue
                    fig = Figure(figsize=(2.45, 1.7)); ax = fig.add_subplot(111)
                    # Draw one data trace (normally common to all group members)
                    # and each member's fitted curve, over the complete stored x range.
                    peak0, arr0 = available[0]
                    ax.plot(arr0[:, 0], arr0[:, 1], color='r', linewidth=0.55, label='Data')
                    xmin, xmax = np.inf, -np.inf
                    for peak, arr in available:
                        ax.plot(arr[:, 0], arr[:, 2], color='b', linewidth=0.5, label=str(peak))
                        if len(arr): xmin=min(xmin, float(np.nanmin(arr[:,0]))); xmax=max(xmax, float(np.nanmax(arr[:,0])))
                    for peak in peaks:
                        try: pos=float(fits[peak].get('f01'))
                        except (TypeError, ValueError): continue
                        ax.axvline(pos, linewidth=0.55, linestyle='--')
                        ax.text(pos, 0.97, str(peak), rotation=90, va='top', ha='right', fontsize=5.5, transform=ax.get_xaxis_transform())
                    if np.isfinite(xmin) and np.isfinite(xmax): ax.set_xlim(xmax, xmin)
                    title = ('Group %s' % unit.get('group')) if unit.get('group') is not None else str(peak0)
                    ax.set_title('%s / exp %d' % (title, slice_index + 1), fontsize=7)
                    ax.tick_params(labelsize=6); fig.tight_layout(pad=0.35)
                    name = 'pseudo2d_fit_%03d_slice_%03d.pdf' % (number, slice_index + 1)
                    fig.savefig(os.path.join(report_dir, name), bbox_inches='tight'); group_slices.append(name)
                unit['slice_figures'] = group_slices
                paths.append((overview, unit))
        finally:
            if old_column:
                try: self.owner._set_axis(old_column)
                except Exception: pass
            self.current_peak, self.slice_index = old_peak, old_slice
            self.owner.draw_figure(keepaxes=False)
        return paths

    def on_peak_selected(self, event):
        self.current_peak = self.peakList.GetItemText(event.GetIndex())
        old_slice = self.slice_index
        self.slices = self._parse_dat(self.current_peak)
        self.current_fit = self._parse_out(self.current_peak)
        self.sliceChoice.Clear()
        for i in range(len(self.slices)): self.sliceChoice.Append(str(i + 1))
        self.slice_index = min(old_slice, max(0, len(self.slices)-1))
        if self.slices: self.sliceChoice.SetSelection(self.slice_index)
        try:
            self.owner.set_fitting_peak(float(self.current_fit.get('f01')), self.slices, self.current_peak)
        except (TypeError, ValueError, AttributeError):
            pass
        self.draw_current()

    def on_slice_choice(self, event=None):
        sel = self.sliceChoice.GetSelection()
        if sel >= 0: self.slice_index = sel
        self.draw_current()

    def step_slice(self, delta):
        if not self.slices: return
        self.slice_index = (self.slice_index + delta) % len(self.slices)
        self.sliceChoice.SetSelection(self.slice_index)
        self.draw_current()

    def on_display_change(self, event=None): self.draw_current(keepaxes=True)
    def redraw_view(self, event=None): self.draw_current(keepaxes=False)

    def on_mark_intensities_correct(self, event=None):
        """Persist explicit acceptance of the currently reviewed pseudo2D fits."""
        # Do not allow acceptance of an incomplete fit set.  The same Full 1D
        # peak list / .dat+.out criterion is used by workflow completion for
        # the extraction stage.
        if not self.fitting_files_complete():
            wx.MessageBox(
                'The pseudo2D intensity series cannot be marked as correct because '
                'one or more Full 1D peaks do not have complete fitting results.',
                'Fitting', wx.OK | wx.ICON_WARNING)
            return False

        notebook = getattr(self.owner, 'parent', None)
        pass
        mark = getattr(notebook, 'mark_workflow_series_inspected', None)
        if callable(mark):
            try:
                if mark():
                    self.correctButton.SetLabel('Intensities marked as correct')
                    self.correctButton.Disable()
                    return True
            except Exception as exc:
                wx.MessageBox(
                    'The pseudo2D intensity review could not be saved.\n\n%s' % exc,
                    'Fitting', wx.OK | wx.ICON_ERROR)
                return False

        # Standalone embedding uses the same application service boundary as
        # the notebook workflow; the service owns persistence and notification.
        if not self.pseudo_service.mark_series_reviewed(source='pseudo2d_fitting'):
            return False
        self.correctButton.SetLabel('Intensities marked as correct')
        self.correctButton.Disable()
        return True

    def _sync_review_button(self):
        reviewed = self.pseudo_service.series_reviewed()
        pass
        self.correctButton.SetLabel('Intensities marked as correct' if reviewed else 'Mark intensities as correct')
        self.correctButton.Enable(not reviewed)

    def draw_current(self, keepaxes=False):
        oldx = self.axes.get_xlim() if keepaxes and self.axes.has_data() else None
        oldy = self.axes.get_ylim() if keepaxes and self.axes.has_data() else None
        self.axes.clear()
        if self.slices:
            arr = self.slices[self.slice_index]
            self.axes.plot(arr[:,0], arr[:,1], color='r', linewidth=0.7, label='Data')
            self.axes.plot(arr[:,0], arr[:,2], color='b', linewidth=0.7, label='Calc')
            try:
                peak_ppm = float(self.current_fit.get('f01', ''))
                fit_intensities = self.current_fit.get('intensity', [])
                if self.slice_index < len(fit_intensities):
                    self.axes.vlines(peak_ppm, 0.0, fit_intensities[self.slice_index], color='k', linewidth=1.2, label='Selected peak')
            except (TypeError, ValueError):
                pass
            if self.diffBox.GetValue(): self.axes.plot(arr[:,0], arr[:,1]-arr[:,2], label='Difference')
            self.axes.axhline(0.0, linewidth=0.5)
            self.axes.set_xlabel('F1 (ppm)'); self.axes.set_ylabel('Intensity')
            self.axes.set_title('%s — Slice %d' % (self.current_peak, self.slice_index + 1))
            self.axes.legend(loc='best')
            if oldx is not None:
                lo, hi = sorted(oldx); self.axes.set_xlim(hi, lo)
            elif len(arr):
                self.axes.set_xlim(float(np.nanmax(arr[:,0])), float(np.nanmin(arr[:,0])))
            if oldy is not None: self.axes.set_ylim(oldy)
        self.canvas.draw_idle()

    def _on_activate_debug(self, event):
        pass
        event.Skip()

    def _on_show_debug(self, event):
        pass
        event.Skip()

    def _on_destroy_debug(self, event):
        pass
        event.Skip()

    def _on_canvas_destroy_debug(self, event):
        pass
        event.Skip()

    def on_close(self, event):
        pass
        try:
            self.Hide()
            pass
            wx.CallAfter(self.owner._resume_motion_after_fitting)
            pass
        except Exception as exc:
            pass
            raise

