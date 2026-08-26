from pathlib import Path


def test_projection_overlay_uses_full_peak_list_not_conn_data():
    source = Path('gui/workspaces/projection.py').read_text()
    block = source[source.index('    def _cached_peak_overlay'):source.index('    def _peak_points_for_overlay')]
    assert 'full_peak_payload' in block
    assert 'conn_data' not in block.replace('legacy ``conn_data``', '')
    assert 'self.tabOne=' not in source
