from pathlib import Path


def _source():
    return (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'pseudo2d.py').read_text()


def test_pseudo2d_analysis_selector_has_requested_choices_and_save_button():
    source = _source()
    assert "return ['Diffusion', 'Decay']" in source
    assert "'uSTA'" not in source[source.index('def available_downstream_analyses'):source.index('def selected_downstream_analysis')]
    assert "wx.Button(panel, -1, 'Save'" in source
    assert "self.save_downstream_analysis(selection)" in source
    assert 'self.pseudo_service.set_downstream_analysis(selection)' in source


def test_pseudo2d_toolbar_exposes_analysis_selector():
    source = _source()
    assert "wx.Button(self.toolbar, -1, 'Analysis'" in source
    assert 'self.analysisToolButton.Bind(wx.EVT_BUTTON, self.show_analysis_selector)' in source


def test_diffusion_opens_existing_pseudo_dimensional_diffusion_panel_in_new_frame():
    source = _source()
    assert "if selection == 'Diffusion':" in source
    assert 'return self.open_diffusion_analysis()' in source
    assert 'from spinDecon.gui.workspaces import pseudo2d_diffusion' in source
    assert "wx.Frame(self, title='Diffusion Analysis'" in source
    assert 'frame.tabOne = self.tabOne' in source
    assert 'pseudo2d_diffusion.Pseudo2DDiffusion(frame, self.tabOne)' in source
