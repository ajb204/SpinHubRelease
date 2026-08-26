from pathlib import Path
from spinDecon.analysis.full_peak_service import FullPeakListService


class Store:
    def __init__(self): self.saved = None
    def save_peak_list(self, key, **payload): self.saved = (key, payload)

class Legacy:
    def __init__(self):
        self.store = Store()
        self.calls = []
        self.conn_data = ['legacy-only']
    def get_full_peak_payload(self): return {'peaks': ['canonical']}
    def _rebuild_projected_peak_lists(self): self.calls.append('rebuild')
    def _notify_analysis_changed(self): self.calls.append('notify')
    def refresh_full_peak_list_viewers(self): self.calls.append('refresh')
    def Status(self): self.calls.append('status')


def test_full_peak_service_is_authoritative_and_does_not_expose_conn_data():
    legacy = Legacy()
    service = FullPeakListService(legacy)
    assert service.payload()['peaks'] == ['canonical']
    assert not hasattr(service, 'connections')
    assert not hasattr(service, 'conn_data')
    service.save(peaks=['p'], records=['r'], rows=['row'], dimension=3)
    assert legacy.store.saved[0] == 'full'
    assert legacy.conn_data == ['legacy-only']
    assert legacy.calls == ['rebuild', 'notify', 'refresh', 'status']


def test_connection_compatibility_is_not_declared_as_peak_authority():
    text = Path(__file__).parents[1].joinpath('analysis', 'full_peak_service.py').read_text()
    assert "save_peak_list('full'" in text
    assert "save_peak_list('reference'" not in text
