from pathlib import Path




def test_active_code_does_not_import_indiana_backend():
    root = Path(__file__).resolve().parents[1]
    active_roots = ["app", "analysis", "domain", "gui", "integrations", "processing", "project", "workflow"]
    for dirname in active_roots:
        for path in (root / dirname).rglob("*.py"):
            text = path.read_text(errors="ignore")
            assert "decon.INDIANA" not in text, path
            assert "INDIANA.cellDiff" not in text, path
