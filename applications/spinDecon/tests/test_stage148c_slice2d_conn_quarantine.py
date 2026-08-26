from pathlib import Path


def test_active_slice2d_has_no_conn_data_or_magma_api():
    source = Path('gui/workspaces/slice2d.py').read_text()
    assert 'conn_data' not in source
    assert '.connections' not in source
    assert 'MAGMA' not in source
    assert 'NOETest' not in source
    assert 'FishNOE' not in source
    assert 'AddNOE' not in source


def test_slice_service_does_not_expose_legacy_connectivity():
    source = Path('analysis/slice_service.py').read_text()
    assert 'conn_data' not in source
    assert 'def connections' not in source
    assert 'load_connections' not in source
    assert 'analyse_connections' not in source


def test_conn_data_slice2d_code_is_preserved_only_in_legacy():
    legacy = Path('legacy/slice2d/conn_data_workspace_snapshot.py')
    assert legacy.exists()
    source = legacy.read_text()
    assert 'conn_data' in source
    assert 'LEGACY_NOETest' in source
