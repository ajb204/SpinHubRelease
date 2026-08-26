from pathlib import Path
ROOT=Path(__file__).parents[1]

def test_radius_inspector_is_separate_and_opens_to_right():
    text=(ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert 'class FitRadiusFrame(wx.Frame)' in text
    assert 'self.radius_window = FitRadiusFrame(self)' in text
    assert 'self.radius_window.SetPosition(wx.Point(x+w+8,y))' in text
    assert "label='Save'" not in text and "label='Close'" in text
    assert 'self.status=self.CreateStatusBar(2)' in text

def test_radius_view_covers_1d_and_two_spectral_modes_and_refreshes():
    text=(ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert 'self.show_3d_peak_preview = bool(self.dim in (1, 2))' in text
    assert "if self.owner.dim == 1:" in text
    assert "projection='3d'" in text
    assert text.count('self.draw_3d_peak_preview()') >= 2

def test_guess_uses_current_fitted_peak_shape_at_ten_percent():
    text=(ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert '_radii_from_current_peak_shape(level=.10)' in text
    assert 'estimate_level_radius(self.owner.data,self.owner.indexes,self.owner.maxima,level=.10)' not in text
    assert '_shape_radius_at_level' in text

def test_radius_save_uses_normal_parent_save_load_path():
    text=(ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert "self.tabOne.fitRadBox.SetValue('%.6g'%self.radius_f1)" in text
    assert "box.SetValue('%.6g'%value)" in text
    assert "parent_save=getattr(self.owner.tabOne,'OnButtonSave',None)" in text

def test_radius_window_does_not_depend_on_visible_unidec_radius_widgets():
    text=(ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert "getattr(self.tabOne, 'fitF%dBox' % which, None)" in text
    assert "getattr(self.tabOne, 'fitRadBox', None)" in text

def test_2d_radius_surfaces_are_peak_centred_and_intensity_synchronised():
    text = (ROOT/'gui'/'workspaces'/'peak_fit.py').read_text()
    assert "-float(self.owner.indexes[1][i1])" in text
    assert "-float(self.owner.indexes[0][i0])" in text
    assert "patch=patch*(reference/centre_value)" in text
    assert "ax.set_xlim(-self.owner.radius_f2,self.owner.radius_f2)" in text
    assert "ax.set_ylim(-self.owner.radius_f1,self.owner.radius_f1)" in text
    assert "Peak-centred representative surfaces" in text
