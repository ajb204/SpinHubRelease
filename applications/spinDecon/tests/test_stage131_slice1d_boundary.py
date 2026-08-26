from spinDecon.analysis.slice_service import SliceService


class _Control:
    def GetValue(self): return '0.1'


class _Legacy:
    dmax = 10.0
    threshBox = _Control()
    DECON = 1
    labb = ['1H', '15N']
    uc0min, uc0max = 1.0, 9.0
    uc1min, uc1max = 100.0, 130.0


def test_slice_service_exposes_viewer_metadata_without_gui_parent_chain():
    service = SliceService(_Legacy())
    assert service.decon_enabled is True
    assert service.label(0) == '1H'
    assert service.axis_limits(1) == (100.0, 130.0)
