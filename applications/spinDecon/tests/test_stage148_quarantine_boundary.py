from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parents[1]
ACTIVE_PACKAGES = ("analysis", "app", "domain", "gui", "integrations", "line_fitting", "processing", "project", "workflow")
REMOVED_NAMESPACES = ("Frames", "misc", "INDIANA", "archive", "tools")
QUARANTINED_CURRENT_BRIDGES = (
    ROOT / "gui" / "workspaces" / "slice4d.py",
    ROOT / "integrations" / "magma",
    ROOT / "integrations" / "usta",
)


def test_obsolete_top_level_namespaces_are_absent():
    assert [name for name in REMOVED_NAMESPACES if (ROOT / name).exists()] == []


def test_retired_usta_indiana_magma_slice4d_are_preserved_under_legacy():
    required = (
        ROOT / "legacy" / "usta",
        ROOT / "legacy" / "indiana",
        ROOT / "legacy" / "magma",
        ROOT / "legacy" / "slice4d",
    )
    assert [str(path) for path in required if not path.exists()] == []
    assert [str(path) for path in QUARANTINED_CURRENT_BRIDGES if path.exists()] == []


def test_compatibility_source_is_quarantined():
    compat = ROOT / "legacy" / "compatibility"
    assert (compat / "root").is_dir()
    assert (compat / "Frames").is_dir()
    assert (compat / "misc").is_dir()


def test_active_source_does_not_import_quarantined_namespaces():
    forbidden = ("decon.legacy", "decon.Frames", "decon.misc", "decon.INDIANA")
    offenders = []
    for package in ACTIVE_PACKAGES:
        for path in (ROOT / package).rglob("*.py"):
            tree = ast.parse(path.read_text(errors="ignore"))
            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.extend(alias.name for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imports.append(node.module)
            hits = [name for name in imports if any(name == token or name.startswith(token + ".") for token in forbidden)]
            if hits:
                offenders.append((str(path.relative_to(ROOT)), hits))
    assert offenders == []


def test_canonical_journey_workspaces_remain_present():
    required = (
        "oned.py", "pseudo2d.py", "slice2d.py", "pseudo3d.py", "full3d.py",
        "nmr.py", "projection.py", "slice1d.py", "peak_review.py", "peak_fit.py",
        "cpmg.py", "decay.py", "phasing.py", "workflow.py",
    )
    base = ROOT / "gui" / "workspaces"
    assert [name for name in required if not (base / name).is_file()] == []
