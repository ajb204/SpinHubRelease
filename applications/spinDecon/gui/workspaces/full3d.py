#!/usr/bin/python
import wx, string, copy, math, numpy, os
from spinDecon.gui.context import context_for, project_for, data_for
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
from matplotlib.widgets import Slider
# from .frameFeatures import drawing_box, contour_box

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
class Full3D(wx.Panel):
    """Interactive orthogonal slice viewer for a 3D spectrum.

    Scientific data and coordinate interpretation are owned by the main
    deconFrame/DataStore.  This panel retains only GUI state and matplotlib
    artists.  Static slice changes use a clean redraw; mouse-driven traces use
    blitting for smooth interaction.
    """

    def __init__(self, parent, tabOne):
        wx.Panel.__init__(self, parent=parent)
        self.tabOne = tabOne
        self.parent = parent
        self.app_context = context_for(tabOne, parent)
        self.state = project_for(tabOne, parent)
        self.store = data_for(tabOne, parent)
        self.full3d = getattr(self.app_context, "full3d", None) if self.app_context is not None else None

        self.bore_dim = 1
        self.n = 0
        self.horizontal = False
        self.vertical = False
        self.trackers_locked = False
        self.last_mouse_x = None
        self.last_mouse_y = None
        self._slider_updating = False
        self._blit_background = None
        self._view = None
        self._drawing_static = False

        try:
            self.thresh = float(self.store.dmax) * float(tabOne.threshBox.GetValue())
        except Exception:
            self.thresh = 1.0

        self.create_main_panel()
        self._set_initial_slice()
        self.draw_slider()
        self.draw_figure(keepaxes=False)
        self._refresh_control_availability()
        self.Show(True)
        self.Fit()

    # ------------------------------------------------------------------ GUI
    def drawing_box(self):
        # Keep toolbar controls borderless, matching the cleaner Slice2D row.
        self.vbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.contourbutton = wx.Button(self, -1, "Contours", size=(-1, 22))
        self.cb_grid = wx.CheckBox(self, -1, "Peaks", style=wx.ALIGN_RIGHT)
        # Legacy state widgets are retained for existing plotting code but no longer shown.
        self.contourbutton.Hide()
        self.cb_grid.Hide()
        self.cb_calc = _ToolbarToggleState(False)
        self.Bind(wx.EVT_BUTTON, self.on_contour_button, self.contourbutton)
        self.Bind(wx.EVT_CHECKBOX, self.on_draw_button, self.cb_grid)

    def contour_box(self):
        # Keep contour state owned by Full3D, but present the editors in a
        # modeless utility window rather than consuming space in the main row.
        self.contourFrame = wx.Frame(self, title='Contours', style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(self.contourFrame)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.text1 = wx.StaticText(panel, -1, 'Min:')
        self.text2 = wx.StaticText(panel, -1, 'Factor:')
        self.text3 = wx.StaticText(panel, -1, 'Number:')
        self.textbox0 = wx.TextCtrl(panel, size=(82, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox1 = wx.TextCtrl(panel, size=(50, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox2 = wx.TextCtrl(panel, size=(50, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox0.SetValue(str(self.thresh))
        self.textbox1.SetValue('1.2')
        self.textbox2.SetValue('15')
        for ctrl in (self.textbox0, self.textbox1, self.textbox2):
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_draw_button)
        for widget in (self.text1, self.textbox0, self.text2, self.textbox1, self.text3, self.textbox2):
            sizer.Add(widget, 0, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.contourCloseButton = wx.Button(panel, -1, "Close", size=(-1, 22))
        self.contourCloseButton.Bind(wx.EVT_BUTTON, self._hide_contour_frame)
        sizer.Add(self.contourCloseButton, 0, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        panel.SetSizerAndFit(sizer)
        frame_sizer = wx.BoxSizer(wx.VERTICAL)
        frame_sizer.Add(panel, 1, wx.EXPAND)
        self.contourFrame.SetSizerAndFit(frame_sizer)
        self.contourFrame.Bind(wx.EVT_CLOSE, self._hide_contour_frame)

    def _hide_contour_frame(self, event):
        self.contourFrame.Hide()
        if event is not None:
            event.Veto()

    def on_contour_button(self, event):
        if not self.contourFrame.IsShown():
            self.contourFrame.Show()
        self.contourFrame.Raise()
        self.textbox0.SetFocus()

    def control_box(self):
        # Slice/orientation controls now live entirely in the native Matplotlib
        # toolbar.  State remains on Full3D; the toolbar is only the controller.
        self.controlSizer = wx.BoxSizer(wx.HORIZONTAL)

    def create_main_panel(self):
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(100, 100))
        self.toolbar = RedrawNavigationToolbar(
            self.canvas, self.redraw_view,
            peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon,
            contour_callback=self._toolbar_contours,
            slice_up_callback=lambda: self.on_next_button(None),
            slice_down_callback=lambda: self.on_prev_button(None),
            transpose_callback=lambda: self.on_help_button(None),
            horizontal_callback=self._toolbar_horizontal,
            vertical_callback=self._toolbar_vertical,
            horizontal_active=self.horizontal, vertical_active=self.vertical,
            coordinates=False)
        self.slider_fig = Figure()
        self.slider_canvas = FigCanvas(self, -1, self.slider_fig)
        self.slider_canvas.SetMinSize(wx.Size(-1, 20))

        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP
        self.drawing_box()
        self.contour_box()
        self.control_box()

        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND | wx.LEFT | wx.TOP | wx.RIGHT, 3)
        # The slice-position strip stays immediately below the spectrum but is
        # kept shallow so the spectrum receives essentially all spare height.
        self.vbox.Add(self.slider_canvas, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 3)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 3)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.cid_click = self.canvas.mpl_connect('button_press_event', self.on_click)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.on_key)
        # Keep the blit cache in step with native Matplotlib navigation.
        # Zoom/pan/home/back/forward redraw the canvas themselves; Full3D must
        # adopt that new canvas as its static background rather than rebuilding
        # the scientific plot (which would overwrite the navigation limits).
        self.cid_draw = self.canvas.mpl_connect('draw_event', self._on_canvas_draw)
        self.Bind(wx.EVT_SIZE, self.OnSize)

    # -------------------------------------------------------------- model API
    def _spec(self):
        return self.full3d.view_spec(self.bore_dim)

    def _slice_view(self):
        return self.full3d.slice_view(self.bore_dim, self.n)

    def _set_initial_slice(self):
        spec = self._spec()
        if spec is None or len(spec['slice_scale']) == 0:
            self.n = 0
        else:
            self.n = min(10, len(spec['slice_scale']) - 1)

    def _has_decon_data(self):
        return getattr(self.store, 'datadec', None) is not None

    def _has_full_peaks(self):
        try:
            return self.full3d.has_full_peaks()
        except Exception:
            return False

    def _refresh_control_availability(self):
        has_decon = self._has_decon_data()
        has_peaks = self._has_full_peaks()
        self.cb_calc.Enable(has_decon)
        self.cb_grid.Enable(has_peaks)
        if not has_decon:
            self.cb_calc.SetValue(False)
        if not has_peaks:
            self.cb_grid.SetValue(False)
        if hasattr(self, 'toolbar'):
            self.toolbar.enable_peaks(has_peaks)
            self.toolbar.set_peaks_active(self.cb_grid.GetValue())
            self.toolbar.enable_decon(has_decon)
            self.toolbar.set_decon_active(self.cb_calc.GetValue())

    # ------------------------------------------------------------- navigation
    def set_slice_index(self, value, sync_slider=True, keepaxes=True):
        """Public slice-navigation API used by deconFrame peak selection."""
        self._set_slice_index(value, sync_slider=sync_slider, keepaxes=keepaxes)

    def _set_slice_index(self, value, sync_slider=True, keepaxes=True):
        spec = self._spec()
        if spec is None or len(spec['slice_scale']) == 0:
            return
        new_n = max(0, min(int(value), len(spec['slice_scale']) - 1))
        self.n = new_n
        if sync_slider and hasattr(self, 'slider'):
            self._slider_updating = True
            try:
                self.slider.set_val(float(spec['slice_scale'][self.n]))
                self.slider_canvas.draw_idle()
            finally:
                self._slider_updating = False
        self.draw_figure(keepaxes=keepaxes)

    def on_prev_button(self, event):
        self._set_slice_index(self.n + 1)

    def on_next_button(self, event):
        self._set_slice_index(self.n - 1)

    def _transpose(self, step):
        """Change orientation, retaining a locked physical 3D position.

        A locked tracker represents a point in the three original spectral
        dimensions.  On transpose the displayed x/y/slice roles change, so
        preserve that physical point rather than preserving an unrelated slice
        index.  Unlocked transpose retains the historical index behaviour.
        """
        old_view = self._slice_view()
        locked_position = None
        if self.trackers_locked and old_view is not None:
            x = self.last_mouse_x
            y = self.last_mouse_y
            if x is not None and y is not None:
                locked_position = {
                    old_view['slice_dim']: float(old_view['slice_value']),
                    old_view['x_dim']: float(x),
                    old_view['y_dim']: float(y),
                }

        self.bore_dim = (self.bore_dim + int(step)) % 3
        spec = self._spec()
        if spec is not None:
            scale = numpy.asarray(spec['slice_scale'])
            if locked_position is not None and len(scale):
                target = locked_position[spec['slice_dim']]
                self.n = int(numpy.argmin(numpy.abs(scale - target)))
                self.last_mouse_x = locked_position[spec['x_dim']]
                self.last_mouse_y = locked_position[spec['y_dim']]
            else:
                self.n = max(0, min(self.n, len(scale) - 1))
        self.draw_slider()
        self.draw_figure(keepaxes=False)

    def on_help_button(self, event):
        self._transpose(1)

    def _transpose_reverse(self):
        self._transpose(-1)

    def draw_slider(self):
        self.slider_fig.clear()
        ax = self.slider_fig.add_subplot(111)
        spec = self._spec()
        if spec is None or len(spec['slice_scale']) == 0:
            self.slider_canvas.draw_idle()
            return
        scale = numpy.asarray(spec['slice_scale'])
        self.n = max(0, min(self.n, len(scale) - 1))
        self.slider = Slider(
            ax, spec['slice_label'], float(numpy.min(scale)), float(numpy.max(scale)),
            valinit=float(scale[self.n]), valfmt='%.2f', color='gray'
        )
        def slider_update(_value):
            if self._slider_updating:
                return
            idx = int(numpy.argmin(numpy.abs(scale - float(self.slider.val))))
            if idx != self.n:
                self.n = idx
                self.draw_figure(keepaxes=True)
        self.slider.on_changed(slider_update)
        self.slider_canvas.draw_idle()

    def update_bore_number(self):
        """Compatibility wrapper retained for legacy callers."""
        spec = self._spec()
        if spec is not None:
            self.n = max(0, min(self.n, len(spec['slice_scale']) - 1))

    def return_current_slice(self):
        view = self._slice_view()
        return None if view is None else view['raw'].T

    def return_current_slice_dec(self):
        view = self._slice_view()
        if view is None or view['decon'] is None:
            return None
        return view['decon'].T

    # ------------------------------------------------------------- trace state
    def _toolbar_horizontal(self, active):
        self.horizontal = bool(active)
        self._refresh_dynamic_display()

    def _toolbar_vertical(self, active):
        self.vertical = bool(active)
        self._refresh_dynamic_display()

    # Compatibility handlers for any legacy callers.
    def on_horiz_button(self, event):
        active = bool(event.GetEventObject().GetValue()) if event is not None else not self.horizontal
        self.toolbar.set_horizontal_active(active)
        self._toolbar_horizontal(active)

    def on_vert_button(self, event):
        active = bool(event.GetEventObject().GetValue()) if event is not None else not self.vertical
        self.toolbar.set_vertical_active(active)
        self._toolbar_vertical(active)

    def _refresh_dynamic_display(self):
        """Refresh trace data without rebuilding or changing the viewport.

        Native Matplotlib navigation redraws only the non-animated artists.
        Recalculate the marginal traces before displaying them so an H/V toggle
        works immediately after Zoom/Home/Back/Forward as well as after a Full3D
        scientific redraw.
        """
        if not hasattr(self, 'line_h'):
            return
        x, y = self.last_mouse_x, self.last_mouse_y
        if x is None or y is None:
            x, y = self._default_tracker_position()
        if x is None or y is None:
            self._configure_trace_visibility()
            self._blit_dynamic()
            return
        self._update_dynamic_from_position(x, y)

    def _mouse_main_coordinates(self, event):
        if not hasattr(self, 'axes') or event.x is None or event.y is None:
            return None, None
        if not self.axes.bbox.contains(event.x, event.y):
            return None, None
        try:
            x, y = self.axes.transData.inverted().transform((event.x, event.y))
            return float(x), float(y)
        except Exception:
            return None, None

    def on_click(self, event):
        if getattr(event, 'button', None) != 1:
            return
        x, y = self._mouse_main_coordinates(event)
        if x is None or y is None:
            return
        self.last_mouse_x = x
        self.last_mouse_y = y
        self.trackers_locked = not self.trackers_locked
        self._update_dynamic_from_position(x, y)

    def on_move(self, event):
        if self.trackers_locked:
            return
        x, y = self._mouse_main_coordinates(event)
        if x is None or y is None:
            return
        self.last_mouse_x = x
        self.last_mouse_y = y
        self._update_dynamic_from_position(x, y)

    def _default_tracker_position(self):
        view = self._view or self._slice_view()
        if view is None:
            return None, None
        x = float(view['x_scale'][len(view['x_scale']) // 2])
        y = float(view['y_scale'][len(view['y_scale']) // 2])
        return x, y

    def _update_dynamic_from_position(self, x, y):
        traces = self.full3d.cross_sections(self.bore_dim, self.n, x, y)
        if traces is None:
            return
        self.last_mouse_x = traces['x_ppm']
        self.last_mouse_y = traces['y_ppm']

        self.line_h.set_data(traces['horizontal_axis'], traces['horizontal'])
        self.line_v.set_data(traces['vertical'], traces['vertical_axis'])
        if traces['horizontal_decon'] is not None:
            self.line_h_dec.set_data(traces['horizontal_axis'], traces['horizontal_decon'])
        if traces['vertical_decon'] is not None:
            self.line_v_dec.set_data(traces['vertical_decon'], traces['vertical_axis'])
        self.cross_h.set_ydata([traces['y_ppm'], traces['y_ppm']])
        self.cross_v.set_xdata([traces['x_ppm'], traces['x_ppm']])

        self.coord_x.set_text('%s: %.2f ppm' % (self._view['x_label'], traces['x_ppm']))
        self.coord_y.set_text('%s: %.2f ppm' % (self._view['y_label'], traces['y_ppm']))
        self.coord_z.set_text('%s: %.2f ppm%s' % (
            self._view['slice_label'], self._view['slice_value'],
            '  [locked]' if self.trackers_locked else ''
        ))
        self._configure_trace_limits(traces)
        self._configure_trace_visibility()
        self._blit_dynamic()

    def _configure_trace_limits(self, traces=None):
        """Keep marginal intensity axes fixed to the complete raw 3D cube."""
        limits = self.full3d.intensity_limits()
        if limits is None:
            return
        lo, hi = limits
        self.axes_h.set_ylim(lo, hi)
        self.axes_v.set_xlim(lo, hi)

    def _configure_trace_visibility(self):
        if not hasattr(self, 'line_h'):
            return
        show_calc = bool(self.cb_calc.IsChecked()) and self._has_decon_data()
        self.axes_h.yaxis.set_visible(bool(self.horizontal))
        self.axes_v.xaxis.set_visible(bool(self.vertical))
        self.line_h.set_visible(bool(self.horizontal))
        self.cross_h.set_visible(bool(self.horizontal))
        self.line_h_dec.set_visible(bool(self.horizontal and show_calc))
        self.line_v.set_visible(bool(self.vertical))
        self.cross_v.set_visible(bool(self.vertical))
        self.line_v_dec.set_visible(bool(self.vertical and show_calc))

    # ------------------------------------------------------------- rendering
    def GetLevels(self):
        try:
            min_level = float(self.textbox0.GetValue())
        except Exception:
            min_level = self.thresh
        try:
            factor = float(self.textbox1.GetValue())
        except Exception:
            factor = 1.2
        try:
            count = int(self.textbox2.GetValue())
        except Exception:
            count = 15
        if count <= 0:
            count = 10
        if factor == 0:
            factor = 1.2
        if min_level == 0:
            min_level = 1e3
        pos = [min_level]
        for _ in range(count - 1):
            pos.append(pos[-1] * factor)
        pos = numpy.asarray(pos)
        return numpy.concatenate((-pos[::-1], pos))

    def _draw_peak_overlay(self):
        if not self.cb_grid.IsChecked():
            return
        for peak in self.full3d.peak_overlay(self.bore_dim, self.n):
            color = peak.get('color', '0.15')
            self.axes.scatter(peak['x'], peak['y'], color=color, marker='x', s=80, zorder=10)
            if peak.get('label'):
                self.axes.annotate(
                    peak['label'], (peak['x'], peak['y']), xytext=(3, 3),
                    textcoords='offset points', color=color, fontsize=10, zorder=10
                )

    def draw_figure(self, keepaxes=True):
        old_xlim = old_ylim = None
        if keepaxes and hasattr(self, 'axes'):
            try:
                old_xlim = self.axes.get_xlim()
                old_ylim = self.axes.get_ylim()
            except Exception:
                pass

        self._refresh_control_availability()
        view = self._slice_view()
        if view is None:
            return
        self._view = view
        self.n = int(view['slice_index'])
        levels = self.GetLevels()

        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        # The toolbar reports coordinates in the current displayed frame.  x and y
        # follow the two visible axes; z is the fixed coordinate selected by the
        # slider.  This deliberately follows the displayed graph after transpose.
        fixed_z = float(view['slice_value'])
        self.axes.format_coord = lambda x, y, z=fixed_z: 'x=%.4f, y=%.4f, z=%.4f' % (x, y, z)
        self.axes_h = self.axes.twinx()   # x=ppm, y=intensity
        self.axes_v = self.axes.twiny()   # x=intensity, y=ppm
        # Use more of the canvas: keep only the margins needed for labels/ticks.
        self.fig.subplots_adjust(left=0.085, bottom=0.10, right=0.955, top=0.975)
        self.axes_h.patch.set_visible(False)
        self.axes_v.patch.set_visible(False)
        self.axes_h.set_navigate(False)
        self.axes_v.set_navigate(False)
        self.axes_h.set_zorder(self.axes.get_zorder() - 1)
        self.axes_v.set_zorder(self.axes.get_zorder() - 1)

        norm = colors.Normalize(vmin=-numpy.max(levels), vmax=numpy.max(levels))
        self.axes.contour(view['x_scale'], view['y_scale'], view['raw'], levels, cmap=cm.seismic, norm=norm)
        if self.cb_calc.IsChecked() and view['decon'] is not None:
            self.axes.contour(view['x_scale'], view['y_scale'], view['decon'], levels, cmap=cm.gray, norm=norm)
        self._draw_peak_overlay()

        # Match the main axis-title font sizes to their corresponding tick labels.
        x_tick_labels = self.axes.get_xticklabels()
        y_tick_labels = self.axes.get_yticklabels()
        x_label_size = x_tick_labels[0].get_fontsize() if x_tick_labels else matplotlib.rcParams['xtick.labelsize']
        y_label_size = y_tick_labels[0].get_fontsize() if y_tick_labels else matplotlib.rcParams['ytick.labelsize']
        self.axes.set_xlabel('%s (ppm)' % view['x_label'], fontsize=x_label_size)
        self.axes.set_ylabel('%s (ppm)' % view['y_label'], fontsize=y_label_size)
        self.axes_h.set_ylabel('Intensity')
        self.axes_v.set_xlabel('Intensity')
        if old_xlim is not None and old_ylim is not None:
            self.axes.set_xlim(old_xlim)
            self.axes.set_ylim(old_ylim)
        else:
            self.axes.set_xlim(float(view['x_scale'][0]), float(view['x_scale'][-1]))
            self.axes.set_ylim(float(view['y_scale'][0]), float(view['y_scale'][-1]))

        # Dynamic artists.  They are excluded from the background and redrawn
        # by blitting on mouse motion.
        x0, y0 = self.last_mouse_x, self.last_mouse_y
        if x0 is None or y0 is None:
            x0, y0 = self._default_tracker_position()
            self.last_mouse_x, self.last_mouse_y = x0, y0
        self.line_h, = self.axes_h.plot([], [], color='r', lw=0.6)
        self.line_h_dec, = self.axes_h.plot([], [], color='b', lw=0.6)
        self.line_v, = self.axes_v.plot([], [], color='r', lw=0.6)
        self.line_v_dec, = self.axes_v.plot([], [], color='b', lw=0.6)
        self.cross_h = self.axes.axhline(y0, color='0.25', lw=0.5)
        self.cross_v = self.axes.axvline(x0, color='0.25', lw=0.5)
        self.coord_x = self.axes.text(0.995, 0.99, '', ha='right', va='top', transform=self.axes.transAxes)
        self.coord_y = self.axes.text(0.995, 0.95, '', ha='right', va='top', transform=self.axes.transAxes)
        self.coord_z = self.axes.text(0.995, 0.91, '', ha='right', va='top', transform=self.axes.transAxes)

        # Dynamic trace/crosshair artists must not be baked into a normal
        # Matplotlib draw.  This lets native navigation redraw the static axes
        # and lets draw_event safely adopt that result as the new blit cache.
        for artist in self._dynamic_artists():
            artist.set_animated(True)

        # The Full3D view has several decorated axes.  Newer Matplotlib can
        # warn when tight_layout cannot satisfy all of their margins, so use
        # explicit stable margins instead of suppressing the warning.
        self.fig.subplots_adjust(left=0.10, right=0.94, bottom=0.11, top=0.94, wspace=0.08, hspace=0.08)
        self._configure_trace_visibility()
        self._drawing_static = True
        try:
            self.canvas.draw()
            self._blit_background = self.canvas.copy_from_bbox(self.fig.bbox)
        finally:
            self._drawing_static = False

        self._update_dynamic_from_position(x0, y0)

    def update_view(self):
        """Compatibility name: slice changes now use a clean redraw."""
        self.draw_figure(keepaxes=True)

    def _dynamic_artists(self):
        return [
            self.line_h, self.line_h_dec, self.line_v, self.line_v_dec,
            self.cross_h, self.cross_v, self.coord_x, self.coord_y, self.coord_z,
        ]

    def _on_canvas_draw(self, event):
        """Adopt a completed native Matplotlib draw as the blit background.

        This is intentionally not a call to draw_figure().  Navigation tools
        own the axes limits; Full3D only refreshes the pixel cache used for the
        horizontal/vertical traces and locked crosshair.
        """
        if self._drawing_static or not hasattr(self, 'axes'):
            return
        try:
            self._blit_background = self.canvas.copy_from_bbox(self.fig.bbox)
            self._blit_dynamic()
        except Exception:
            self._blit_background = None

    def _blit_dynamic(self):
        if self._blit_background is None or not hasattr(self, 'axes'):
            self.canvas.draw_idle()
            return
        try:
            self.canvas.restore_region(self._blit_background)
            for artist in self._dynamic_artists():
                if artist.get_visible():
                    artist.axes.draw_artist(artist)
            self.canvas.blit(self.fig.bbox)
        except Exception:
            # Resizes/backend changes can invalidate a background.  Recover
            # with one clean redraw rather than leaving blitting artefacts.
            self.canvas.draw_idle()

    def save_background(self):
        """Compatibility wrapper: refresh pixels without rebuilding axes."""
        self.canvas.draw_idle()

    def _toolbar_decon(self, active):
        self.cb_calc.SetValue(bool(active))
        self.draw_figure(keepaxes=True)

    def _toolbar_peaks(self, active):
        self.cb_grid.SetValue(bool(active))
        self.on_draw_button(None)

    def _toolbar_contours(self):
        self.on_contour_button(None)

    def redraw_view(self):
        # Explicit pencil redraw clears the transient peak selection and
        # rebuilds Full3D at its complete view.  Other Matplotlib navigation
        # tools are deliberately left untouched.
        try:
            self.full3d.clear_peak_selection(redraw_full3d=False)
        except Exception:
            pass
        self.draw_figure(keepaxes=False)

        # The pencil is the explicit reset action: draw_figure(False) uses the
        # complete ppm extent of the current plane rather than retaining any
        # native Matplotlib zoom/pan limits.

    def on_draw_button(self, event):
        # Checkboxes also use this handler and must not clear peak selection.
        self.draw_figure(keepaxes=True)

    def on_key(self, event):
        if event.key == 'up':
            self.on_next_button(None)
        elif event.key == 'down':
            self.on_prev_button(None)
        elif event.key == 'q':
            self.on_help_button(None)
        elif event.key == 'w':
            self._transpose_reverse()
        elif event.key == 'h':
            self.horizontal = not self.horizontal
            self.toolbar.set_horizontal_active(self.horizontal)
            self._refresh_dynamic_display()
        elif event.key == 'v':
            self.vertical = not self.vertical
            self.toolbar.set_vertical_active(self.vertical)
            self._refresh_dynamic_display()

    def OnSize(self, event):
        event.Skip()
        # A resize invalidates only the pixel cache, not the NMR plot or its
        # viewport.  Let wx/Matplotlib perform its native resize draw; the
        # draw_event handler will capture the resulting background.
        self._blit_background = None
        if self.IsShownOnScreen():
            wx.CallAfter(self.canvas.draw_idle)
