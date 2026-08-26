"""Regression checks for GUI context migration import boundaries."""
import ast
from pathlib import Path


def test_gui_context_imports_only_context_helpers():
    root = Path(__file__).resolve().parents[1]
    allowed = {"context_for", "project_for", "data_for"}
    offenders = []
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "decon.gui.context":
                imported = {alias.name for alias in node.names}
                unexpected = imported - allowed
                if unexpected:
                    offenders.append((str(path.relative_to(root)), sorted(unexpected)))
    assert not offenders, offenders


def test_gui_context_helpers_are_imported_when_used():
    root = Path(__file__).resolve().parents[1]
    helpers = {"context_for", "project_for", "data_for"}
    offenders = []
    context_module = root / "gui" / "context.py"
    for path in root.rglob("*.py"):
        if path == context_module:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        used = {
            node.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Name)
            and isinstance(node.ctx, ast.Load)
            and node.id in helpers
        }
        if not used:
            continue

        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "decon.gui.context":
                imported.update(alias.asname or alias.name for alias in node.names)

        missing = used - imported
        if missing:
            offenders.append((str(path.relative_to(root)), sorted(missing)))

    assert not offenders, offenders
