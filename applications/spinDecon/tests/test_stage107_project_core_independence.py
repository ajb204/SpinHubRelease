import ast
from pathlib import Path
ROOT=Path(__file__).parents[1]
CORE=('data_store.py','decon_service.py','defaults.py','parameter_store.py','state.py')

def _imports(path):
    tree=ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            yield node.module

def test_project_core_is_gui_independent():
    bad=[]
    for name in CORE:
        path=ROOT/'project'/name
        for module in _imports(path):
            if module=='wx' or module.startswith('wx.') or module.startswith('decon.gui'):
                bad.append((name,module))
    assert bad==[]
