from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def test_pseudo2d_raw_viewer_and_diffusion_are_separate():
    viewer = (ROOT / 'gui' / 'workspaces' / 'pseudo2d.py').read_text()
    diffusion = (ROOT / 'gui' / 'workspaces' / 'pseudo2d_diffusion.py').read_text()
    tabs = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'class Pseudo2D(wx.Panel)' in viewer
    assert 'PseudoAxisTable.load' in viewer
    assert "contour_callback=self.on_contours" in viewer
    assert "mpl_connect('motion_notify_event'" in viewer
    assert "callbacks.connect('xlim_changed'" in viewer
    assert 'class Pseudo2DDiffusion(wx.Panel)' in diffusion
    assert 'AddTabPseudo2DDiffusion' in tabs

def test_pseudo2d_load_opens_projection_then_inspector():
    source = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('        # Set up tabs from the canonical topology established by makeinp().')
    block = source[start:source.index('        self.READ=1', start)]
    assert 'topology = self._active_topology()' in block
    assert 'if spectral_dim_count == 1:' in block
    assert 'if topology.has_pseudo_axis:' in block
    assert 'self.parent.AddTabTwo(True, self)' in block
    assert 'self.parent.AddTabPseudo2D(True, self)' in block
    assert "self.parent.select_page('Projections')" in block
    assert 'if(self.dim' not in block


def test_pseudo2d_does_not_require_filesystem_projection_cache():
    source = (ROOT / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('        needs_projection_cache = bool(')
    block = source[start:source.index('        if needs_projection_cache:', start)]
    assert 'topology.spectral_dim_count == 1' not in block

def test_projection_has_one_spectral_plus_pseudo_mode():
    source = (ROOT / 'gui' / 'workspaces' / 'projection.py').read_text()
    assert 'def _is_pseudo2d_projection_case' in source
    assert 'numpy.sum(data, axis=0)' in source


def test_pseudo2d_motion_uses_primary_axis_coordinates_and_blitting():
    viewer = (ROOT / 'gui' / 'workspaces' / 'pseudo2d.py').read_text()
    assert "self.axes.bbox.contains(event.x, event.y)" in viewer
    assert "self.axes.transData.inverted().transform((event.x, event.y))" in viewer
    assert "self.trace_line.set_data(self.x, trace)" in viewer
    assert "animated=True" in viewer
    assert "copy_from_bbox(self.fig.bbox)" in viewer
    assert "self.trace_axes.draw_artist(self.trace_line)" in viewer
    assert "self.axes.draw_artist(self.selection_line)" in viewer
    assert "self.canvas.blit(self.fig.bbox)" in viewer


def test_pseudo2d_redraw_resets_both_axes():
    viewer = (ROOT / 'gui' / 'workspaces' / 'pseudo2d.py').read_text()
    start = viewer.index('    def redraw_view(self, event=None):')
    block = viewer[start:viewer.index('    def onFocus', start)]
    assert 'self.draw_figure(keepaxes=False)' in block
    assert 'keepaxes=True' not in block
