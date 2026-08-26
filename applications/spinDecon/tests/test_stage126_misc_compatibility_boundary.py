import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]
ACTIVE = ("app", "gui", "analysis", "project", "workflow", "processing", "integrations", "domain")


def test_active_packages_do_not_import_misc_compatibility_modules():
    offenders = []
    for package in ACTIVE:
        for path in (ROOT / package).rglob("*.py"):
            tree = ast.parse(path.read_text(), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom):
                    names = [node.module or ""]
                else:
                    continue
                if any(name == "decon.misc" or name.startswith("decon.misc.") for name in names):
                    offenders.append(str(path.relative_to(ROOT)))
    assert offenders == []
