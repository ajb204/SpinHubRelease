import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_decon_service_is_gui_independent():
    path = ROOT / 'project' / 'decon_service.py'
    text = path.read_text()
    tree = ast.parse(text)
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
        elif isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
    assert not any(name == 'wx' or name.startswith('decon.gui') or name.startswith('decon.misc.shell_output') for name in imports)
    assert 'launch_with_output' not in text
