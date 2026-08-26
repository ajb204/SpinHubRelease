from pathlib import Path


def test_1d_and_projection_services_do_not_expose_conn_data():
    for path in (Path("analysis/oned_service.py"), Path("analysis/projection_service.py")):
        source = path.read_text()
        assert "conn_data" not in source
        assert "def connections" not in source


def test_current_1d_and_projection_gui_do_not_consume_conn_data():
    for path in (Path("gui/workspaces/oned.py"), Path("gui/workspaces/projection.py")):
        source = path.read_text()
        assert ".connections" not in source
        assert "self.conn_data" not in source


def test_cpmg_and_decay_do_not_contain_legacy_noe_manager():
    for path in (Path("gui/workspaces/cpmg.py"), Path("gui/workspaces/decay.py")):
        source = path.read_text()
        assert "class NOEMan" not in source
        assert "class NOEManFrame" not in source


def test_conn_data_gui_history_is_preserved_in_legacy():
    root = Path("legacy/conn_data")
    assert (root / "README.md").exists()
    assert (root / "oned_pre148d_snapshot.py").exists()
    assert (root / "projection_pre148d_snapshot.py").exists()
    assert (root / "cpmg_pre148d_snapshot.py").exists()
    assert (root / "decay_pre148d_snapshot.py").exists()
