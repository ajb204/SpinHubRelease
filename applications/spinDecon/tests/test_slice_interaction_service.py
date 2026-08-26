from spinDecon.analysis.slice_service import SliceService

class Control:
    def GetValue(self): return '0.1'
    def IsChecked(self): return True
class Topology:
    spectral_dim_count = 3
class Legacy:
    dmax=10; threshBox=Control(); cb_grid=Control(); dim=3; labb=('A','B','C')
    DECON=1
    def _active_topology(self): return Topology()
    def Status(self): self.status_called=True

def test_slice_interaction_boundary():
    legacy=Legacy(); service=SliceService(legacy)
    assert service.spectral_dimension == 3
    assert service.symmetry_enabled
    assert not hasattr(service, 'connections')
    assert not hasattr(service, 'load_connections')
    assert not hasattr(service, 'analyse_connections')
    service.refresh_status(); assert legacy.status_called

def test_slice_meshes_and_decon_alias_are_service_owned():
    import numpy
    legacy = Legacy()
    legacy.datadec = numpy.arange(24).reshape(2, 3, 4)
    legacy.XX = numpy.arange(24).reshape(2, 3, 4)
    legacy.YY = legacy.XX + 100
    legacy.ZZ = legacy.XX + 200
    service = SliceService(legacy)
    assert service.datadec is legacy.datadec
    xs, ys = service.slice_meshes(orth=False, pt_c=1, pt_h_max=0, pt_h_min=3)
    numpy.testing.assert_array_equal(xs, legacy.XX[:, 1, 0:3].transpose())
    numpy.testing.assert_array_equal(ys, legacy.ZZ[:, 1, 0:3].transpose())
    xs, ys = service.slice_meshes(orth=True, pt_h=2, pt_c_max=0, pt_c_min=2)
    numpy.testing.assert_array_equal(xs, legacy.XX[:, 0:2, 2].transpose())
    numpy.testing.assert_array_equal(ys, legacy.YY[:, 0:2, 2].transpose())
