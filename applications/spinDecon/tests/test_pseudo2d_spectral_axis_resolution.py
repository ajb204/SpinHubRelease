from pathlib import Path


def test_pseudo2d_spectral_axis_uses_canonical_topology_not_real_label():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('    def _spectral_physical_axes(self):')
    end = source.index('    def _spectral_axis_labels(self):', start)
    block = source[start:end]
    assert 'topology = self._active_topology()' in block
    assert 'for axis in topology.spectral_axes:' in block
    assert "_REAL_AXIS_LABELS" not in block
