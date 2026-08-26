from pathlib import Path
TEXT=(Path(__file__).parents[1]/'gui'/'workspaces'/'peak_fit.py').read_text()

def test_save_buttons_removed_from_fit_windows():
    assert "save=wx.Button(panel,label='Save')" not in TEXT
    assert 'self.saveButton = wx.Button' not in TEXT

def test_both_windows_snapshot_and_prompt_on_close():
    assert 'self._opening_values = self._snapshot_values()' in TEXT
    assert 'self._opening_fit_values = self._snapshot_fit_values()' in TEXT
    assert "'The extraction radius has changed. Save changes?'" in TEXT
    assert "'Peak-fit parameters have changed. Save changes?'" in TEXT
    assert 'wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION' in TEXT

def test_parent_close_coordinates_radius_without_double_prompt():
    assert 'self._closing_pair=True' in TEXT
    assert "not getattr(self.owner, '_closing_pair', False)" in TEXT
    assert 'if radius_changed:' in TEXT and 'radius.on_radius_entered()' in TEXT
    assert 'self._save_fit_changes()' in TEXT

def test_opening_snapshot_is_safe_before_matplotlib_sliders_exist():
    text = TEXT
    method = text.split('def _snapshot_fit_values(self):', 1)[1].split('def _fit_has_unsaved_changes', 1)[0]
    assert "getattr(self, 'psf_sliders', None)" in method
    assert 'sliders_ready' in method
    assert "getattr(self.tabOne, 'sig%dBox' % which)" in method
