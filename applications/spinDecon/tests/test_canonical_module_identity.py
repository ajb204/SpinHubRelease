

def test_data_store_root_alias_preserves_identity():
    from spinDecon.project.data_store import DataStore as LegacyDataStore
    from spinDecon.project.data_store import DataStore as CanonicalDataStore
    assert LegacyDataStore is CanonicalDataStore
