from pathlib import Path


def test_oned_uses_service_for_scientific_data_after_construction():
    source = Path('gui/workspaces/oned.py').read_text()
    assert 'OneDService(tabOne)' in source
    assert 'else tabOne.' not in source
