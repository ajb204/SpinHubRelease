from spinDecon.analysis.peak_service import PeakService


class Store:
    def __init__(self):
        self.views = {}
    def get_view(self, key):
        return self.views.get(key)
    def save_view(self, key, **payload):
        self.views[key] = payload


class Legacy:
    dim = 2
    labb = ('N', 'H')
    dmax = 1.0
    store = Store()
    peak = []
    spectrumfile = 'spec.ft2'
    def _is_pseudo3d_topology(self):
        return False
    def get_spectrum_view(self, decon=False, transpose='n'):
        return {'ZZ': [[2]] if decon else [[1]], 'x_axis': [8.0], 'y_axis': [120.0]}
    def get_projection_view(self, *args, **kwargs):
        raise AssertionError('2D data should use spectrum_view')


def test_peak_service_owns_2d_plane_selection():
    service = PeakService(Legacy())
    assert service.view_labels() == ('H', 'N')
    assert service.projection_payload(('H', 'N'))['ZZ'] == [[1]]
    assert service.projection_payload(('H', 'N'), decon=True)['ZZ'] == [[2]]
