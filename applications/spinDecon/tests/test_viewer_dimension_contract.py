from types import SimpleNamespace
import pytest
from spinDecon.project.state import ProjectState
from spinDecon.domain.dimensions.viewer_contract import topology_for, spectral_dim_count, physical_dim_count, spectral_physical_indices, pseudo_physical_index

@pytest.mark.parametrize('spectral,pseudo,physical', [
 (1,False,1),(1,True,2),(2,False,2),(2,True,3),
 (3,False,3),(3,True,4),(4,False,4),(4,True,5),
])
def test_viewers_share_canonical_counts(spectral,pseudo,physical):
    state=ProjectState(dimension=spectral,pseudo_axis=pseudo)
    tab=SimpleNamespace(state=state, dim=spectral)
    assert spectral_dim_count(tab)==spectral
    assert physical_dim_count(tab)==physical
    assert len(topology_for(tab).spectral_axes)==spectral

def test_viewer_axis_identity_distinguishes_physical_and_spectral():
    state=ProjectState(dimension=2,pseudo_axis=True)
    tab=SimpleNamespace(state=state,dim=2)
    assert pseudo_physical_index(tab)==0
    assert spectral_physical_indices(tab)==(1,2)
