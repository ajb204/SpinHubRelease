from pathlib import Path


def test_show_fid_toggle_requests_y_axis_reset():
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "process.py").read_text()
    start = source.index("    def on_show_fid(self, event):")
    end = source.index("    def on_fid_select(self, event):", start)
    handler = source[start:end]
    assert "self.draw_figure(reset_y=True)" in handler
