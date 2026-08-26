import pytest

from spinDecon.domain.topology import DatasetTopology
from spinDecon.domain.dimensions.peak_contract import (
    peak_coordinate_count, peak_list_suffix, spectral_physical_indices,
)

@pytest.mark.parametrize('spectral,pseudo', [
    (1, False), (1, True), (2, False), (2, True),
    (3, False), (3, True), (4, False), (4, True),
])
def test_peak_coordinates_are_always_spectral(spectral, pseudo):
    topology = DatasetTopology.from_counts(spectral, pseudo)
    assert peak_coordinate_count(topology) == spectral
    assert peak_list_suffix(topology) == f'.{spectral}D.list'


def test_pseudo_axis_is_not_a_peak_coordinate():
    topology = DatasetTopology.from_counts(2, True, pseudo_physical_index=1)
    assert topology.physical_dim_count == 3
    assert spectral_physical_indices(topology) == (0, 2)
    assert peak_coordinate_count(topology) == 2


def test_peak_fit_frame_uses_topology_not_legacy_pseudo_flag():
    """Pseudo fitting must not depend on the independently mutable GUI flag."""
    from pathlib import Path
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'peak_fit.py').read_text()
    init = source[source.index('    def __init__'):source.index('    def on_scroll')]
    assert 'topology = self.tabOne._active_topology()' in init
    assert 'topology.has_pseudo_axis' in init
    assert 'self.tabOne.pseudo' not in init
    assert 'self.data.ndim != self.dim' in init
