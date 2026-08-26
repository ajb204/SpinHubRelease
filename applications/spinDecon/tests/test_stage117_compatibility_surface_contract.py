from pathlib import Path

ROOT_COMPAT = {
    "analysis_mode.py", "data_store.py", "dataset_topology.py", "decon_service.py",
    "decon_tab.py", "dimension_guard.py", "dimension_labels.py", "parameter_store.py",
    "pdfViewer.py", "peak_dimension_contract.py", "peak_picker.py", "peak_shape_estimator.py",
    "project_defaults.py", "project_service.py", "project_setup.py", "project_state.py",
    "project_summary.py", "pseudo_axis_table.py", "shiftXPostFilter.py",
    "viewer_dimension_contract.py", "workflow_model.py", "workflow_overview.py",
    "workflow_registry.py", "workflow_status.py",
}




def test_misc_is_compatibility_only():
    root = Path(__file__).resolve().parents[1]
    for path in (root / "misc").glob("*.py"):
        if path.name == "__init__.py":
            continue
        assert len(path.read_text(errors="ignore").splitlines()) <= 12, path.name
