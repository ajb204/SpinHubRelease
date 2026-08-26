from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SLICE = (ROOT / "gui" / "workspaces" / "slice2d.py").read_text()
SERVICE = (ROOT / "analysis" / "slice_service.py").read_text()
PEAK_LIST = (ROOT / "gui" / "workspaces" / "full_peak_list.py").read_text()


def test_slice2d_has_no_local_peak_or_connection_manager_ui():
    assert "def peaks_box(" not in SLICE
    assert "self.NOEbutton" not in SLICE
    assert "def update_conn_data(" not in SLICE
    assert "def on_search_button(" not in SLICE
    assert "def on_save_button(" not in SLICE
    assert "def on_load_button(" not in SLICE


def test_full_peak_list_is_authoritative_and_conn_data_is_explicitly_legacy():
    assert "self.slice_service.open_full_peak_list(event)" in SLICE
    assert "FUTURE FULL-PEAK / NOE MIGRATION MARKER" in PEAK_LIST
    assert "conn_data is legacy and is not an authoritative peak store" in PEAK_LIST
    assert "def open_full_peak_list" in SERVICE
