from pathlib import Path

ROOT = Path(__file__).parents[1]




def test_analysis_helpers_have_canonical_owners():
    assert (ROOT / 'analysis/peak_shape_optimizer.py').is_file()
    assert (ROOT / 'analysis/shiftx_post_filter.py').is_file()
    peak_fit = (ROOT / 'gui/workspaces/peak_fit.py').read_text()
    projection = (ROOT / 'gui/workspaces/projection.py').read_text()
    assert 'decon.analysis import peak_shape_optimizer' in peak_fit
    assert 'decon.analysis.shiftx_post_filter' in projection


def test_retired_usta_simulation_prototype_is_quarantined():
    historical = ROOT / 'legacy/usta/uSTA_sims_frame.py'
    assert historical.is_file()
    # It was already entirely commented out and is not an active workspace.
    assert 'class uSTA_sims_frame(wx.Panel):' in historical.read_text()
    shell = (ROOT / 'app' / 'notebook.py').read_text()
    assert 'from .Frames.uSTA.uSTA_sims_frame import uSTA_sims_frame' not in '\n'.join(
        line for line in shell.splitlines() if not line.lstrip().startswith('#')
    )
