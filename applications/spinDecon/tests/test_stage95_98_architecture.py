from pathlib import Path

ROOT = Path(__file__).parents[1]






def test_active_packages_do_not_import_root_compatibility_modules():
    compatibility = {
        'analysis_mode', 'dataset_topology', 'dimension_guard', 'dimension_labels',
        'peak_dimension_contract', 'viewer_dimension_contract', 'pseudo_axis_table',
        'workflow_model', 'workflow_status', 'workflow_overview', 'data_store',
        'parameter_store', 'project_defaults', 'project_service', 'project_setup',
        'project_state', 'project_summary', 'decon_service', 'peak_picker',
        'peak_shape_estimator', 'PeakShapeOptimizer', 'shiftXPostFilter', 'pdfViewer',
    }
    offenders = []
    for package in ('app', 'analysis', 'domain', 'gui', 'integrations', 'processing', 'project', 'workflow'):
        for path in (ROOT / package).rglob('*.py'):
            text = path.read_text(errors='ignore')
            for module in compatibility:
                if f'from spinDecon.{module} import' in text:
                    offenders.append(f'{path.relative_to(ROOT)} -> {module}')
    assert offenders == []
