from pathlib import Path


def test_slice1d_view_markers_are_built_from_full_peak_payload():
    source = Path('gui/workspaces/nmr.py').read_text()
    start = source.index('    def get_reference_1d_view(self, peak_index):')
    end = source.index('    def get_full3d_view_spec', start)
    method = source[start:end]
    assert 'self.get_full_peak_payload()' in method
    assert "record.get('axis_values')" in method
    assert "record.get('intensity')" in method
    assert 'conn_data' not in method.replace('conn_data used to provide these markers', '')
