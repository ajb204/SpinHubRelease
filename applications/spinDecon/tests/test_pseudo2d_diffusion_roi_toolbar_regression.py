from pathlib import Path

SRC = Path(__file__).parents[1] / "gui" / "workspaces" / "pseudo2d_diffusion.py"
TOOLBAR = Path(__file__).parents[1] / "gui" / "plotting" / "toolbar.py"


def test_diffusion_uses_project_toolbar_and_embeds_controls():
    text = SRC.read_text()
    assert "RedrawNavigationToolbar" in text
    assert '_redraw_tool = ("Draw!"' in text
    assert "self.savebutton = wx.Button(tb" in text
    assert "tb.AddControl(self.noiseMax)" not in text  # controls are added through the ordered tuple
    assert "self.roi_button = wx.Button(tb" in text
    assert "self.vbox.Add(self.hbox" not in text


def test_project_toolbar_hides_home_and_subplots():
    text = TOOLBAR.read_text()
    assert '_HIDDEN_TOOL_NAMES = {"Home", "Subplots"}' in text


def test_diffusion_roi_manager_is_live_and_removable():
    text = SRC.read_text()
    assert "class DiffusionROIFrame" in text
    assert "self.diffusion.remove_rois" in text
    assert "self.diffusion.highlight_rois" in text
    assert "self._append_roi(a, b)" in text
    assert "self._append_roi(start, end)" in text
    assert "def _rebuild_roi_overlays" in text
    assert "getattr(event, 'dblclick', False)" in text


def test_roi_summary_table_has_fit_statistics_and_merge():
    text = SRC.read_text()
    assert "wx.ListCtrl(panel, style=wx.LC_REPORT)" in text
    assert '"D (fit)"' in text
    assert '"Fit error"' in text
    assert '"D (Gaussian)"' in text
    assert 'label="Merge selected"' in text
    assert "def merge_rois" in text
    assert "def _gradient_fit" in text
    assert "numpy.linalg.lstsq(design, logy, rcond=None)" in text
    assert '"Gaussian error"' in text
    assert "stats.get('gaussian_sigma')" in text


def test_region_gaussian_fit_is_drawn_in_main_histogram_axis():
    text = SRC.read_text()
    assert "def _plot_roi_histograms" in text
    assert "self.axes_err.plot(gaussian['xfit'], gaussian['yfit']" in text
    assert "stats.get('kind') != 'region'" in text


def test_roi_attenuation_fit_does_not_depend_on_preliminary_amplitude_fit():
    text = SRC.read_text()
    build = text[text.index("def _build_roi_data"):text.index("def _append_roi")]
    assert "/ self.asv" not in build
    assert "amps = self.asv" not in build
    assert "self.normalisedDiffFull[:, lo:hi]" in build
    assert "numpy.nanmean(region, axis=1)" in build
    gradient = text[text.index("def _gradient_fit"):text.index("def _gaussian_fit_for_range")]
    assert "mask.sum() < 2" in gradient
    assert "len(xfit) > 2" in gradient
