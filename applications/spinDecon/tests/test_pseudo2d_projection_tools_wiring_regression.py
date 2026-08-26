from pathlib import Path


def test_pseudo2d_projection_canvas_is_wired_to_peak_tool_click_handler():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'projection.py').read_text()
    assert "self._pseudo2d_tool_click_cid = self.canvas.mpl_connect(" in source
    assert "'button_press_event', self.on_pick)" in source
    assert "self._handle_full_tool_click(event)" in source


def test_pseudo2d_projection_tools_mutate_authoritative_full_1d_list():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'projection.py').read_text()
    assert "save_peak_list('full'" in source
    assert "dimension=1" in source
    assert "self._move_full_at(event.xdata)" in source
    assert "self._add_full_at(event.xdata)" in source
    assert "def on_full_remove" in source
    assert "def on_full_maximise" in source
