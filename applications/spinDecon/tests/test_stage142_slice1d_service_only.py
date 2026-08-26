from pathlib import Path


def test_slice1d_does_not_store_or_read_tabone_after_construction():
    source = Path('gui/workspaces/slice1d.py').read_text()
    assert 'self.tabOne' not in source
    assert 'SliceService(tabOne)' in source
