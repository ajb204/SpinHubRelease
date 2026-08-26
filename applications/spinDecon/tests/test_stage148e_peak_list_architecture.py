from pathlib import Path


def test_active_model_has_no_conn_data_or_conn_entry():
    root = Path(__file__).resolve().parents[1]
    offenders = []
    for package in ("analysis", "app", "domain", "gui", "integrations", "processing", "project", "workflow"):
        for path in (root / package).rglob("*.py"):
            text = path.read_text(errors="ignore")
            code = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
            if "self.conn_data" in code or "connEntry(" in code:
                offenders.append(str(path.relative_to(root)))
    assert offenders == []


def test_datastore_peak_lists_support_multiple_roles():
    from spinDecon.project.data_store import DataStore
    store = DataStore()
    record = {"name": "P1", "coordinates": (1.0, 2.0), "axis_values": {}, "analysis": {}}
    for role in ("reference", "full", "decon"):
        store.save_peak_list(role, role=role, dimensionality=2, records=[record])
    assert set(("reference", "full", "decon")).issubset(store.peak_lists)
    assert all(store.peak_lists[r]["records"][0]["coordinates"] == (1.0, 2.0) for r in ("reference", "full", "decon"))
