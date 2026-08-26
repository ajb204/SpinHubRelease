from pathlib import Path


def test_phase_slider_window_is_docked_exactly_below_projection_frame():
    source = (Path(__file__).parents[1] / 'gui' / 'dialogs' / 'processing' / 'projections.py').read_text()
    start = source.index('    def _position_phase_frame(self) -> None:')
    end = source.index('\n    def _hide_phase_frame', start)
    body = source[start:end]

    assert 'rect = self.GetScreenRect()' in body
    assert 'x = int(rect.x)' in body
    assert 'y = int(rect.y + rect.height)' in body
    assert 'self.phaseFrame.SetPosition((x, y))' in body
    # The old screen-edge fallback could move the sliders above or sideways,
    # which violated exact docking to the Projections frame.
    assert 'GetClientArea' not in body
    assert 'pos.y - psize.height' not in body
