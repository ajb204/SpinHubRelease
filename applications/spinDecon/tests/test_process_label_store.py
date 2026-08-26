from pathlib import Path

from spinDecon.domain.dimensions.labels import canonical_spectral_labels, discover_bruker_labels


def test_duplicate_raw_labels_are_canonicalised():
    assert canonical_spectral_labels(['H', 'H']) == ['H_1', 'H_2']
    assert canonical_spectral_labels(['H', 'N']) == ['H', 'N']


def test_bruker_labels_are_discovered_without_conversion_window(tmp_path):
    (tmp_path / 'acqus').write_text('##$NUC1= <1H>\n')
    (tmp_path / 'acqu2s').write_text('##$NUC1= <1H>\n')
    assert discover_bruker_labels(str(tmp_path), 2) == ['1H', '1H']
    assert canonical_spectral_labels(discover_bruker_labels(str(tmp_path), 2)) == ['1H_1', '1H_2']


def test_process_family_consumes_authoritative_label_api():
    root = Path(__file__).resolve().parents[1]
    process = (root / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
    conversion = (root / 'gui' / 'dialogs' / 'processing' / 'conversion.py').read_text()
    processing = (root / 'gui' / 'dialogs' / 'processing' / 'settings.py').read_text()
    projection = (root / 'gui' / 'dialogs' / 'processing' / 'projections.py').read_text()
    assert 'def get_dimension_labels(self):' in process
    assert 'def set_dimension_labels(self, labels, refresh=True):' in process
    assert "self.parent.ParseAllStr(self._parameter_file_path(), key)" in process
    assert "getattr(self.proc, 'set_dimension_labels', None)" in conversion
    assert "getattr(self.proc, 'get_dimension_labels', None)" in processing
    assert "getattr(self.process_parent, 'get_dimension_labels', None)" in projection
