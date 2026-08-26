from pathlib import Path


SOURCE = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_review.py').read_text()


def test_peakframe_peaklist_button_delegates_to_canonical_reference_viewer():
    block = SOURCE[SOURCE.index('    def onPeakList(self, event):'):]
    block = block[:block.index('\n#######################################################')]
    assert 'self.peak_service.open_reference_peak_list(event)' in block
    assert 'self.tabOne.OnButtonReferencePeakList(event)' not in block


def test_peakframe_contains_no_legacy_peak_manager_implementation():
    assert 'class peakMan(' not in SOURCE
    assert 'class peakManFrame(' not in SOURCE
    assert 'class SortedListCtrl(' not in SOURCE
    assert 'ColumnSorterMixin' not in SOURCE
