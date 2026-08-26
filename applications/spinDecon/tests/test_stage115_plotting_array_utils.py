from pathlib import Path

from spinDecon.gui.plotting.array_utils import ensure_xy_points


def test_canonical_plotting_array_helper_is_available():
    assert callable(ensure_xy_points)


def test_active_gui_does_not_import_misc_array_utils():
    root = Path(__file__).resolve().parents[1]
    offenders = []
    for path in (root / "gui").rglob("*.py"):
        if "decon.misc.array_utils" in path.read_text(errors="ignore"):
            offenders.append(str(path.relative_to(root)))
    assert offenders == []
