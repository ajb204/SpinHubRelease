"""Architecture guardrails for the post-refactor package boundaries."""
import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]
ACTIVE = ('app', 'analysis', 'domain', 'gui', 'integrations', 'processing', 'project', 'workflow')
ROOT_COMPAT = {
    'analysis_mode', 'data_store', 'dataset_topology', 'decon_service', 'decon_tab',
    'dimension_guard', 'dimension_labels', 'parameter_store', 'peak_dimension_contract',
    'project_defaults', 'project_service', 'project_setup', 'project_state', 'project_summary',
    'pseudo_axis_table', 'viewer_dimension_contract', 'workflow_model', 'workflow_overview',
    'workflow_registry', 'workflow_status',
}


def _imports(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            yield node.module


def test_active_packages_do_not_depend_on_root_compatibility_or_frames():
    violations = []
    for package in ACTIVE:
        for path in (ROOT / package).rglob('*.py'):
            for module in _imports(path):
                if module == 'decon.Frames' or module.startswith('decon.Frames.'):
                    violations.append((path, module))
                if module.startswith('decon.') and module.split('.')[1] in ROOT_COMPAT:
                    violations.append((path, module))
    assert not violations, violations
