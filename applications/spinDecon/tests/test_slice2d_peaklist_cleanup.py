from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SLICE = (ROOT / 'gui' / 'workspaces' / 'slice2d.py').read_text()
PEAK_LIST = (ROOT / 'gui' / 'workspaces' / 'full_peak_list.py').read_text()


def test_slice2d_peaks_button_opens_canonical_full_peak_list():
    assert 'self.peaksToolButton.Bind(wx.EVT_BUTTON, self.on_full_peak_list)' in SLICE
    block = SLICE[SLICE.index('    def on_full_peak_list(self, event):'):]
    assert 'return self.slice_service.open_full_peak_list(event)' in block
    assert 'return self.tabOne.OnButtonFullPeakList(event)' not in block


def test_legacy_slice2d_peaks_window_and_conn_manager_are_removed():
    assert 'def peaks_box(' not in SLICE
    # The historical AssMan/AssManFrame conn_data viewer is fully quarantined.
    assert 'bool=AssMan(self)' not in SLICE
    assert 'class AssMan(' not in SLICE
    assert 'class AssManFrame(' not in SLICE
    assert (ROOT / 'legacy' / 'compatibility' / 'gui' / 'slice2d_assignment_viewer.py').exists()
    assert 'self.NOEbutton' not in SLICE
    assert 'def on_select_button(' not in SLICE
    assert 'def on_deselect_button(' not in SLICE
    assert 'def on_delete_button(' not in SLICE
    assert 'def on_add_button(' not in SLICE
    assert 'def NOETest(' not in SLICE
    assert 'def FishNOE(' not in SLICE
    assert 'def AddNOE(' not in SLICE


def test_unported_noe_and_4d_semantics_are_retained_only_as_future_marker():
    marker = PEAK_LIST[PEAK_LIST.index('# FUTURE FULL-PEAK / NOE MIGRATION MARKER'):]
    assert '_FUTURE_FULL_PEAK_NOE_LEGACY' in marker
    assert 'def NOETest(' in marker
    assert 'def AddNOE(' in marker
    assert 'def OnRemove(' in marker
    assert 'def OnSave(' in marker
    assert 'entry.f4' in marker
    assert 'entry.distScore' in marker
    assert 'entry.distppm' in marker
    assert 'reciprocated' in marker.lower()
