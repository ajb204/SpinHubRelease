from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def test_full_list_refreshes_status_and_pseudo2d_uses_full_count():
    text = (ROOT / "gui" / "workspaces" / "nmr.py").read_text()
    assert "self.Status()" in text[text.index("def load_full_peak_list"):text.index("def refresh_reference_peak_views")]
    assert "_full = self.get_full_peak_payload()" in text

def test_pseudo2d_analyse_does_not_parse_connections():
    text = (ROOT / "gui" / "workspaces" / "nmr.py").read_text()
    block = text[text.index("def OnButtonAnalyse"):text.index("def OnButtonAnalyse", text.index("def OnButtonAnalyse")+1) if text.count("def OnButtonAnalyse") > 1 else len(text)]
    assert "self._load_pseudo2d_projection_peaks(self.corrFile)" in block
    assert "self.GetConn(" not in block
    assert "self._load_decon_peak_list(self.corrFile)" in block

def test_pseudo2d_projection_has_no_contour_tool_and_draw_resets_zoom():
    text = (ROOT / "gui" / "workspaces" / "projection.py").read_text()
    assert "contour_callback=(None if self._is_pseudo2d_projection_case() else self._toolbar_contours)" in text
    redraw = text[text.index("def redraw_view"):text.index("def on_draw_button") ]
    assert "self.draw_figure(keepaxes=False)" in redraw

def test_pseudo2d_selection_updates_marker_and_3d_traces():
    text = (ROOT / "gui" / "workspaces" / "pseudo2d.py").read_text()
    assert "def set_fitting_peak" in text
    assert "self.axes.axvline(self.selected_fit_ppm" in text
    assert "projection='3d'" in text
    assert "for yval, trace in zip(yy, z):" in text
    assert "ax.plot(ppm, np.full(ppm.shape, yval, dtype=float), trace" in text
    assert "ax.plot_surface" not in text
    assert "self.owner.set_fitting_peak" in text
    assert "np.nanmax(ppm)), float(np.nanmin(ppm))" in text

def test_pseudo2d_hover_trace_is_transient_and_tracks_nearest_y_slice():
    text = (ROOT / "gui" / "workspaces" / "pseudo2d.py").read_text()
    assert "figure_leave_event', self._on_contour_pointer_leave" in text
    assert "self.trace_line.set_visible(False)" in text
    assert "def _set_slice_trace_visible" in text
    assert "row = int(np.nanargmin(np.abs(self.y - pseudo_y)))" in text
    assert "self._set_slice_trace_visible(True)" in text
    assert "self._set_slice_trace_visible(False)" in text
    assert "self.trace_line.set_data(self.x, trace)" in text


def test_main_nmr_tab_draws_canonical_pseudo2d_projection_on_read():
    text = (ROOT / "gui" / "workspaces" / "nmr.py").read_text()
    assert "def _draw_main_pseudo2d_projection" in text
    helper = text[text.index("def _draw_main_pseudo2d_projection"):text.index("#Write status menu")]
    assert "self.get_pseudo2d_projection_data(ensure_file=False)" in helper
    assert "self.axes.plot(scale, projected" in helper
    assert "self.axes.set_xlim(float(numpy.nanmax(scale)), float(numpy.nanmin(scale)))" in helper
    read = text[text.index("def OnButtonRead"):text.index("def OnButtonExtract")]
    assert "if spectral_dim_count == 1 and topology.has_pseudo_axis:" in read
    assert "self._draw_main_pseudo2d_projection()" in read

def test_pseudo2d_projection_tools_edit_authoritative_full_1d_list():
    text = (ROOT / "gui" / "workspaces" / "projection.py").read_text()
    assert "tools_callback=(self._toolbar_tools if self._is_pseudo2d_projection_case() else None)" in text
    assert "self.fullToolsFrame, panel = self._make_modeless_window('Tools')" in text
    for label in ("Undo", "Redo", "Select", "Move", "Add", "Maximise", "Remove"):
        assert "label='%s'" % label in text
    assert "dimension=1" in text[text.index("def _commit_full_records"):text.index("def on_full_undo")]
    assert "r['coordinates']=(x,)" in text
    assert "r['coordinates']=(float(xv[i]),)" in text
    assert "refresh_full_peak_list_viewers" in text
    assert "self.tabOne.store.save_peak_list('full'" in text
    assert "pos.x + size.width" in text
