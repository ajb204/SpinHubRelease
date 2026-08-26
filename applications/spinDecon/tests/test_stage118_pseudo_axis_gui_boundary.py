import ast
from pathlib import Path


def test_pseudo_axis_domain_has_no_wx_dependency():
    root = Path(__file__).resolve().parents[1]
    path = root / "domain" / "pseudo_axis.py"
    tree = ast.parse(path.read_text())
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    assert not any(name == "wx" or name.startswith("wx.") for name in imports)
    assert "show_pseudo_axis_table" not in path.read_text()


def test_pseudo_axis_table_presentation_is_gui_owned():
    root = Path(__file__).resolve().parents[1]
    path = root / "gui" / "dialogs" / "pseudo_axis.py"
    assert path.exists()
    assert "def show_pseudo_axis_table" in path.read_text()
