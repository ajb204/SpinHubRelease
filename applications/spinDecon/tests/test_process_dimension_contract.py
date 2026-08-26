import pytest

from spinDecon.domain.topology import DatasetTopology
from spinDecon.processing.dimension_contract import (
    conversion_axis_count,
    processing_axis_count,
    legacy_vpar_dimension,
)


@pytest.mark.parametrize(
    "spectral,pseudo,conversion_rows,processing_rows",
    [
        (1, False, 1, 1),
        (1, True, 2, 1),
        (2, False, 2, 2),
        (2, True, 3, 2),
        (3, False, 3, 3),
        (3, True, 4, 3),
        (4, False, 4, 4),
        (4, True, 5, 4),
    ],
)
def test_process_children_use_physical_vs_spectral_contract(
    spectral, pseudo, conversion_rows, processing_rows
):
    topology = DatasetTopology.from_counts(spectral, pseudo)
    assert conversion_axis_count(topology) == conversion_rows
    assert processing_axis_count(topology) == processing_rows


def test_legacy_pseudo_encoding_is_backend_only_adapter():
    assert legacy_vpar_dimension(DatasetTopology.from_counts(1, True)) == "2p"
    assert legacy_vpar_dimension(DatasetTopology.from_counts(2, True)) == "3p"
    assert legacy_vpar_dimension(DatasetTopology.from_counts(2, False)) == 2


def test_pseudo3d_projection_dispatch_precedes_plain_2d():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "process.py").read_text()
    method = source[source.index("    def DoProjections("):source.index("    def RefreshDirectSlice(")]
    pseudo = method.index("if self.has_pseudo_axis and self.spectral_dim_count == 2:")
    plain = method.index("elif self.spectral_dim_count == 2:")
    assert pseudo < plain
    assert "MakeProj3P(pipe_path" in method


def test_pseudo3d_processing_script_does_not_embed_legacy_proj3d_tcl():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "processing" / "nmrpipe_scripts.py").read_text()
    start = source.index("def make_proc_script_3dp(")
    end = source.index("def make_proc_script_3d(", start)
    body = source[start:end]
    assert "make_proj_3dp(frame" not in body
    assert "proj3D.tcl" not in body


def test_process_status_lamps_use_deconmain_visual_style():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "process.py").read_text()
    make_lamp = source[source.index("    def _make_status_lamp("):source.index("    def _dimension_status_text(")]
    set_lamp = source[source.index("    def _set_lamp("):source.index("    def UpdateLampLights(")]
    assert "style=wx.BORDER_SIMPLE" in make_lamp
    assert "wx.Colour(46, 160, 67)" in set_lamp
    assert "wx.Colour(210, 55, 55)" in set_lamp


def test_conversion_status_reuses_physical_fid_candidate_discovery():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "process.py").read_text()
    body = source[source.index("    def _conversion_outputs_exist("):source.index("    def _process_pipefile(")]
    assert "_direct_fid_candidates(show_fid=True)" in body
    assert "os.path.isfile(path)" in body


def test_processed_status_accepts_pseudo3d_spectral_product_and_fallbacks():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "process.py").read_text()
    body = source[source.index("    def _processed_output_exists("):source.index("    def _set_lamp(")]
    assert "pipefile = self._process_pipefile()" in body
    assert "'slice.phased.ft1'" in body
    assert "glob.glob(os.path.join(base, '*.ft*'))" in body
