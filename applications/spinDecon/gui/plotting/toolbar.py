"""Project-wide Matplotlib toolbar extensions."""
from pathlib import Path

import wx
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg

from spinDecon.gui.widgets.status_help import _find_status_bar, bind_status_help

_HIDDEN_TOOL_NAMES = {"Home", "Subplots"}


class RedrawNavigationToolbar(NavigationToolbar2WxAgg):
    """Common NMR toolbar with Peaks, optional Contours and Redraw tools.

    Peaks is a native check tool (like Pan/Zoom), so wx supplies the normal
    depressed/darker active state.  Frames keep ownership of plotting state;
    callbacks are deliberately small adapters supplied by each frame.
    """

    _decon_tool = ("Deconvolution", "Show deconvolved spectrum", "decon_overlay", "toggle_decon")
    _orth_tool = ("Orthogonal", "Show orthogonal slice", "orthogonal", "toggle_orth")
    _one_d_tool = ("1D", "Show 1D slice", "one_d", "toggle_one_d")
    _peak_tool = ("Peaks", "Show peak positions", "peaks", "toggle_peaks")
    _contour_tool = ("Contours", "Contour settings", "contours", "open_contours")
    _sliders_tool = ("Sliders", "Show indirect phasing sliders", "sliders", "toggle_sliders")
    _fid_spectrum_tool = ("FID to spectrum", "FID to spectrum / Re-process projections", "fid_to_spectrum", "fid_to_spectrum")
    _tools_tool = ("Tools", "Peak tools", "pickaxe", "toggle_tools")
    _slice_up_tool = ("Slice Up", "Previous slice", "arrow_up", "slice_up")
    _slice_down_tool = ("Slice Down", "Next slice", "arrow_down", "slice_down")
    _transpose_tool = ("Transpose", "Rotate displayed axes by 90 degrees", "rotate_90", "transpose_view")
    _horizontal_tool = ("Horizontal", "Show horizontal trace", "horizontal_trace", "toggle_horizontal")
    _vertical_tool = ("Vertical", "Show vertical trace", "vertical_trace", "toggle_vertical")
    _redraw_tool = ("Redraw", "Redraw", "redraw_pencil", "redraw_view")
    _hidden_tools = _HIDDEN_TOOL_NAMES
    _native_items = tuple(
        item for item in NavigationToolbar2WxAgg.toolitems
        if item is None or item[0] not in _HIDDEN_TOOL_NAMES
    )
    _asset_dir = Path(__file__).with_name("assets")

    def __init__(self, canvas, redraw_callback, *args, peak_callback=None,
                 decon_callback=None, orth_callback=None, one_d_callback=None,
                 contour_callback=None, sliders_callback=None, reprocess_callback=None, tools_callback=None, slice_up_callback=None,
                 slice_down_callback=None, transpose_callback=None, horizontal_callback=None,
                 vertical_callback=None, peaks_active=False, decon_active=False, orth_active=False,
                 one_d_active=False, tools_active=False, horizontal_active=False, vertical_active=False, coordinates=True,
                 style=wx.TB_BOTTOM, **kwargs):
        self._redraw_callback = redraw_callback
        self._peak_callback = peak_callback
        self._decon_callback = decon_callback
        self._orth_callback = orth_callback
        self._one_d_callback = one_d_callback
        self._contour_callback = contour_callback
        self._sliders_callback = sliders_callback
        self._reprocess_callback = reprocess_callback
        self._tools_callback = tools_callback
        self._slice_up_callback = slice_up_callback
        self._slice_down_callback = slice_down_callback
        self._transpose_callback = transpose_callback
        self._horizontal_callback = horizontal_callback
        self._vertical_callback = vertical_callback
        # Project tools are deliberately ordered with Redraw at the far left,
        # followed by the calculated/deconvolved overlay and then Peaks.
        prefix = [self._redraw_tool]
        if decon_callback is not None:
            prefix.append(self._decon_tool)
        if orth_callback is not None:
            prefix.append(self._orth_tool)
        if one_d_callback is not None:
            prefix.append(self._one_d_tool)
        if peak_callback is not None:
            prefix.append(self._peak_tool)
        if contour_callback is not None:
            prefix.append(self._contour_tool)
        if sliders_callback is not None:
            prefix.append(self._sliders_tool)
        if reprocess_callback is not None:
            prefix.append(self._fid_spectrum_tool)
        if tools_callback is not None:
            prefix.append(self._tools_tool)
        if slice_up_callback is not None:
            prefix.append(self._slice_up_tool)
        if slice_down_callback is not None:
            prefix.append(self._slice_down_tool)
        if transpose_callback is not None:
            prefix.append(self._transpose_tool)
        if horizontal_callback is not None:
            prefix.append(self._horizontal_tool)
        if vertical_callback is not None:
            prefix.append(self._vertical_tool)
        self.toolitems = tuple(prefix) + self._native_items

        # This mirrors Matplotlib's wx constructor, with Peaks added to the
        # native check-tool set.  Matplotlib itself hard-codes only Pan/Zoom.
        wx.ToolBar.__init__(self, canvas.GetParent(), -1, style=style)
        if wx.Platform == '__WXMAC__':
            self.SetToolBitmapSize(self.GetToolBitmapSize() * self.GetDPIScaleFactor())
        self.wx_ids = {}
        self._status_help_by_id = {}
        # Ordered native layout used for hover hit-testing on wxMac.  Cocoa's
        # wx.ToolBar implementation exposes neither FindToolForPosition nor
        # per-tool rectangles, so preserve the exact tool/separator sequence.
        self._status_layout = []
        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.AddSeparator()
                self._status_layout.append(None)
                continue
            self.wx_ids[text] = self.AddTool(
                -1, bitmap=self._icon(f"{image_file}.svg"),
                bmpDisabled=wx.NullBitmap, label=text, shortHelp=tooltip_text,
                kind=(wx.ITEM_CHECK if text in ["Deconvolution", "Orthogonal", "1D", "Peaks", "Tools", "Horizontal", "Vertical", "Sliders", "Pan", "Zoom"] else wx.ITEM_NORMAL)
            ).Id
            self.Bind(wx.EVT_TOOL, getattr(self, callback), id=self.wx_ids[text])
            self._status_help_by_id[self.wx_ids[text]] = tooltip_text
            self._status_layout.append(self.wx_ids[text])
        self._coordinates = coordinates
        if coordinates:
            self.AddStretchableSpace()
            self._label_text = wx.StaticText(self, style=wx.ALIGN_RIGHT)
            self.AddControl(self._label_text)
        self.Realize()
        # wxMac reliably delivers mouse motion for the toolbar, but its native
        # FindToolForPosition() returns None for NavigationToolbar2WxAgg tools.
        # Cache each realised tool rectangle instead and do our own hit testing.
        self._hovered_status_tool_id = None
        self._status_tool_rects = self._build_status_tool_rects()
        self.Bind(wx.EVT_MOTION, self._on_toolbar_motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_toolbar_leave)
        NavigationToolbar2.__init__(self, canvas)
        self.set_peaks_active(peaks_active)
        self.set_decon_active(decon_active)
        self.set_orth_active(orth_active)
        self.set_one_d_active(one_d_active)
        self.set_tools_active(tools_active)
        self.set_horizontal_active(horizontal_active)
        self.set_vertical_active(vertical_active)


    def _build_status_tool_rects(self):
        """Return hit rectangles for native toolbar tools.

        wxMac/Cocoa does not expose native tool rectangles (and
        ``FindToolForPosition`` returns ``None``).  Measurements from the
        realised toolbar show that a native item occupies bitmap width + 6
        logical pixels and a separator occupies 15 logical pixels.  These
        values exactly reproduce the realised widths (for example 9 tools +
        2 separators = 372 px with 32 px bitmaps).  Other platforms first use
        their native rectangle APIs and fall back to the same layout model.
        """
        rects = []

        if wx.Platform != '__WXMAC__':
            get_tool_rect = getattr(self, 'GetToolRect', None)
            if get_tool_rect is not None:
                try:
                    for tool_id in self._status_help_by_id:
                        candidate = get_tool_rect(tool_id)
                        if candidate is not None and candidate.width > 0 and candidate.height > 0:
                            rects.append((tool_id, wx.Rect(candidate)))
                    if len(rects) == len(self._status_help_by_id):
                        return rects
                except Exception:
                    rects = []

        bitmap_width = int(self.GetToolBitmapSize().width)
        tool_width = bitmap_width + 6
        separator_width = 15
        height = max(1, int(self.GetClientSize().height))
        x = 0
        for tool_id in self._status_layout:
            if tool_id is None:
                x += separator_width
            else:
                rects.append((tool_id, wx.Rect(x, 0, tool_width, height)))
                x += tool_width
        return rects

    def _tool_id_at_position(self, pos):
        """Return the native toolbar tool ID containing *pos*, if any."""
        # Height can change after the toolbar is inserted into its final sizer,
        # so rebuild if necessary.  Horizontal geometry is deterministic.
        if (not self._status_tool_rects or
                self._status_tool_rects[0][1].height != self.GetClientSize().height):
            self._status_tool_rects = self._build_status_tool_rects()
        for tool_id, rect in self._status_tool_rects:
            if rect.Contains(pos):
                return tool_id
        return None

    def _set_tool_status_help(self, tool_id):
        """Display help for *tool_id*, restoring prior status when absent."""
        frame, bar = _find_status_bar(self)
        if bar is None:
            return

        text = self._status_help_by_id.get(tool_id)
        active = getattr(frame, '_mpl_toolbar_hover_owner', None)
        if text:
            if active is not self:
                try:
                    frame._mpl_toolbar_hover_previous = bar.GetStatusText()
                except Exception:
                    frame._mpl_toolbar_hover_previous = ''
            frame._mpl_toolbar_hover_owner = self
            frame._mpl_toolbar_hover_text = text
            bar.SetStatusText(text)
        elif active is self:
            try:
                if bar.GetStatusText() == getattr(frame, '_mpl_toolbar_hover_text', None):
                    bar.SetStatusText(getattr(frame, '_mpl_toolbar_hover_previous', '') or '')
            finally:
                frame._mpl_toolbar_hover_owner = None
                frame._mpl_toolbar_hover_text = None
                frame._mpl_toolbar_hover_previous = ''

    def _on_toolbar_motion(self, event):
        """Update status help using cached realised toolbar geometry."""
        try:
            tool_id = self._tool_id_at_position(event.GetPosition())
            if tool_id != self._hovered_status_tool_id:
                self._hovered_status_tool_id = tool_id
                self._set_tool_status_help(tool_id)
        except Exception:
            # Status help must never interfere with plotting/navigation.
            pass
        finally:
            event.Skip()

    def _on_toolbar_leave(self, event):
        """Restore the status text when the pointer leaves the toolbar."""
        try:
            self._hovered_status_tool_id = None
            self._set_tool_status_help(None)
        except Exception:
            pass
        finally:
            event.Skip()

    def bind_control_status_help(self, widget, text):
        """Give an AddControl() widget the same status-bar hover behaviour."""
        bind_status_help(self, widget, text)

    @classmethod
    def _icon(cls, name):
        if name in {"redraw_pencil.svg", "decon_overlay.svg", "orthogonal.svg", "one_d.svg", "peaks.svg", "contours.svg", "sliders.svg", "fid_to_spectrum.svg", "pickaxe.svg", "arrow_up.svg", "arrow_down.svg", "rotate_90.svg", "horizontal_trace.svg", "vertical_trace.svg"}:
            svg = (cls._asset_dir / name).read_bytes()
            try:
                dark = wx.SystemSettings.GetAppearance().IsDark()
            except AttributeError:
                bg = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
                fg = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOWTEXT)
                bg_lum = (.299 * bg.red + .587 * bg.green + .114 * bg.blue) / 255
                fg_lum = (.299 * fg.red + .587 * fg.green + .114 * fg.blue) / 255
                dark = fg_lum - bg_lum > .2
            if dark:
                svg = svg.replace(b"#000000", b"#ffffff")
            size = wx.ArtProvider().GetDIPSizeHint(wx.ART_TOOLBAR)
            return wx.BitmapBundle.FromSVG(svg, size)
        return super()._icon(name)


    def set_decon_active(self, active):
        if "Deconvolution" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Deconvolution"], bool(active))

    def enable_decon(self, enabled=True):
        if "Deconvolution" in self.wx_ids:
            self.EnableTool(self.wx_ids["Deconvolution"], bool(enabled))
            if not enabled:
                self.ToggleTool(self.wx_ids["Deconvolution"], False)

    def toggle_decon(self, event=None):
        active = self.GetToolState(self.wx_ids["Deconvolution"])
        if self._decon_callback is not None:
            self._decon_callback(active)

    def set_orth_active(self, active):
        if "Orthogonal" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Orthogonal"], bool(active))

    def toggle_orth(self, event=None):
        active = self.GetToolState(self.wx_ids["Orthogonal"])
        if self._orth_callback is not None:
            self._orth_callback(active)

    def set_one_d_active(self, active):
        if "1D" in self.wx_ids:
            self.ToggleTool(self.wx_ids["1D"], bool(active))

    def toggle_one_d(self, event=None):
        active = self.GetToolState(self.wx_ids["1D"])
        if self._one_d_callback is not None:
            self._one_d_callback(active)

    def set_peaks_active(self, active):
        if "Peaks" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Peaks"], bool(active))

    def enable_peaks(self, enabled=True):
        if "Peaks" in self.wx_ids:
            self.EnableTool(self.wx_ids["Peaks"], bool(enabled))

    def toggle_peaks(self, event=None):
        active = self.GetToolState(self.wx_ids["Peaks"])
        if self._peak_callback is not None:
            self._peak_callback(active)

    def set_tools_active(self, active):
        if "Tools" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Tools"], bool(active))

    def toggle_tools(self, event=None):
        active = self.GetToolState(self.wx_ids["Tools"])
        if self._tools_callback is not None:
            self._tools_callback(active)


    def slice_up(self, event=None):
        if self._slice_up_callback is not None:
            self._slice_up_callback()

    def slice_down(self, event=None):
        if self._slice_down_callback is not None:
            self._slice_down_callback()

    def transpose_view(self, event=None):
        if self._transpose_callback is not None:
            self._transpose_callback()

    def set_horizontal_active(self, active):
        if "Horizontal" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Horizontal"], bool(active))

    def toggle_horizontal(self, event=None):
        active = self.GetToolState(self.wx_ids["Horizontal"])
        if self._horizontal_callback is not None:
            self._horizontal_callback(active)

    def set_vertical_active(self, active):
        if "Vertical" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Vertical"], bool(active))

    def toggle_vertical(self, event=None):
        active = self.GetToolState(self.wx_ids["Vertical"])
        if self._vertical_callback is not None:
            self._vertical_callback(active)

    def open_contours(self, event=None):
        if self._contour_callback is not None:
            self._contour_callback()

    def set_sliders_active(self, active):
        if "Sliders" in self.wx_ids:
            self.ToggleTool(self.wx_ids["Sliders"], bool(active))

    def toggle_sliders(self, event=None):
        active = self.GetToolState(self.wx_ids["Sliders"])
        if self._sliders_callback is not None:
            self._sliders_callback(active)

    def fid_to_spectrum(self, event=None):
        """Run the frame-specific FID-to-spectrum/re-processing action."""
        if self._reprocess_callback is not None:
            self._reprocess_callback(event)

    def redraw_view(self, event=None):
        if self._redraw_callback is not None:
            self._redraw_callback()
