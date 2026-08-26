from pathlib import Path


def test_peak_shape_revision_does_not_auto_launch_pseudo_extraction():
    """Ordinary GUI edits may invalidate workflow evidence but never advance it."""
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('    def peak_shape_saved(')
    end = source.index('    def _mark_pseudo3d_recompute_complete', start)
    body = source[start:end]
    assert "run_workflow_action" not in body
    assert "wx.CallAfter" not in body
    assert "pseudo_intensities_stale = True" in body
    assert "extraction deferred to workflow" in body
