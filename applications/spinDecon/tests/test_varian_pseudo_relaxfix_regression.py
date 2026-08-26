"""Regression coverage for Varian pseudo-dimensional conversion routing."""

from pathlib import Path


SOURCE_PATH = Path(__file__).parents[1] / 'processing' / 'vpar_decon.py'


def _method_source(name):
    source = SOURCE_PATH.read_text()
    start = source.index(f'    def {name}(')
    next_method = source.find('\n    def ', start + 1)
    return source[start:] if next_method == -1 else source[start:next_method]


def test_pipeparse_relaxfix_is_pseudo3d_only():
    """Varian 2p must not be routed through RelaxFix/fid.final."""
    source = _method_source('PipeParse')

    assert "if(self.dim == '3p'):" in source
    assert "if(self.dim in ('2p', '3p')):" not in source
    assert "if self.dim in ('2p', '3p'):" not in source


def test_varian_pseudo3d_slice_preview_still_uses_relaxfix():
    """The dedicated 3p slice path still requires RelaxFix."""
    source = _method_source('_pipe_parse_varian_3p_slice')

    assert 'RelaxFix.out' in source
    assert 'fid.final' in source
