from pathlib import Path
import ast


def _peakframe_source():
    return Path(__file__).resolve().parents[1].joinpath('gui', 'workspaces', 'peak_review.py').read_text()


def _method_node(name):
    tree = ast.parse(_peakframe_source())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'peakFrame')
    return next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name)


def test_peak_and_label_toggle_is_overlay_only_and_preserves_zoom():
    node = _method_node('on_cb_grid')
    calls = [n for n in ast.walk(node) if isinstance(n, ast.Call)]
    attrs = [n.func.attr for n in calls if isinstance(n.func, ast.Attribute)]
    assert '_refresh_peak_artists' in attrs
    assert 'draw_figure' not in attrs

    assignments = [n for n in ast.walk(node) if isinstance(n, ast.Assign)]
    ax_reset_values = []
    for assignment in assignments:
        for target in assignment.targets:
            if isinstance(target, ast.Attribute) and target.attr == 'ax_reset':
                if isinstance(assignment.value, ast.Constant):
                    ax_reset_values.append(assignment.value.value)
    assert ax_reset_values == [0]


def test_toolbar_peaks_uses_same_zoom_preserving_handler():
    node = _method_node('_toolbar_peaks')
    calls = [n for n in ast.walk(node) if isinstance(n, ast.Call)]
    attrs = [n.func.attr for n in calls if isinstance(n.func, ast.Attribute)]
    assert 'on_cb_grid' in attrs
