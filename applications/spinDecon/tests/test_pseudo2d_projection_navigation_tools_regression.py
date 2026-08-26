from pathlib import Path


def _source():
    return (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'projection.py').read_text()


def test_pseudo2d_tools_cancel_matplotlib_navigation_mode():
    source = _source()
    assert 'def _cancel_projection_navigation(self):' in source
    for handler in ('on_full_select', 'on_full_move', 'on_full_add',
                    'on_full_remove', 'on_full_maximise'):
        start = source.index('    def %s' % handler)
        end = source.find('\n    def ', start + 5)
        body = source[start:end if end != -1 else None]
        assert 'self._cancel_projection_navigation()' in body


def test_projection_redraw_releases_navigation_for_pseudo2d():
    source = _source()
    start = source.index('    def redraw_view(self):')
    end = source.index('    def on_draw_button', start)
    body = source[start:end]
    assert 'self._is_pseudo2d_projection_case()' in body
    assert 'self._cancel_projection_navigation()' in body
