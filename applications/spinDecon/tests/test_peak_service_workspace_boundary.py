from spinDecon.analysis.peak_service import PeakService


class _UC:
    ppms_scale = [9.0, 8.0]


class _Legacy:
    dim = 2
    labb = ['15N', '1H']
    uc0 = _UC()
    uc1 = _UC()
    uc0min, uc0max = 6.0, 10.0
    uc1min, uc1max = 100.0, 130.0
    def _resolve_spec_file(self, path): return '/spec/' + path
    def OnButtonReferencePeakList(self, event): return ('opened', event)
    def set_reference_peaks(self, peaks): self.peak = list(peaks); return self.peak
    def refresh_reference_peak_views(self): self.refreshed = True


def test_peak_service_exposes_workspace_neutral_axis_and_reference_operations():
    legacy = _Legacy()
    service = PeakService(legacy)
    assert service.axis_limits[:2] == ((6.0, 10.0), (100.0, 130.0))
    assert service.axis_scale(0).tolist() == [9.0, 8.0]
    assert service.resolve_spectrum_file('x.dat') == '/spec/x.dat'
    assert service.open_reference_peak_list('evt') == ('opened', 'evt')
    assert service.set_reference_peaks([1, 2]) == [1, 2]
    assert legacy.refreshed is True
