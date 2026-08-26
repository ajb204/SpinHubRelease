from spinDecon.analysis.peak_fit_service import PeakFitService

class Store:
    metadata={}
    def mark_peak_shape_determined(self, **kwargs): self.marked=kwargs
class Legacy:
    peak_shape_fitted=False; uSTA=False
    def __init__(self): self.store=Store()
    def OnButtonSave(self,event): self.saved=event
    def peak_shape_saved(self,was_already_fitted=False): self.completed=was_already_fitted

def test_peak_fit_persistence_boundary():
    legacy=Legacy(); service=PeakFitService(legacy)
    service.save_fit_preferences(7, True)
    assert (legacy.peak_fit_count, legacy.peak_fit_link_widths)==(7,True)
    service.save_project(); assert legacy.saved is None
    service.mark_peak_shape_determined(2,[1,2],[.5,.5],[.5,.5])
    assert legacy.store.marked['dimension']==2
    assert legacy.completed is False

class ParameterLegacy:
    def __init__(self):
        self.received = None

    def update_pseudo3d_parameters(self, values):
        self.received = values
        return "updated"


def test_peak_fit_parameter_update_uses_mapping_contract():
    legacy = ParameterLegacy()
    service = PeakFitService(legacy)
    values = {"3p_radF1": "0.125", "3p_radF2": "0.25"}

    assert service.update_pseudo3d_parameters(values) == "updated"
    assert legacy.received == values
    assert legacy.received is not values
