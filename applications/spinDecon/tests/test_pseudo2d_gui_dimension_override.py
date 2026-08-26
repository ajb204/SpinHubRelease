from pathlib import Path


def _read_block():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('    def OnButtonRead(self,event):')
    end = source.index('    ##############################################', start)
    return source[start:end]


def test_read_uses_project_state_topology_before_makeinp():
    block = _read_block()
    commit = block.index("self.dim = int(getattr(self.state, 'spectral_dimensions', 0) or 1)")
    boundary = block.index("self.makeinp('', resolved)")
    assert commit < boundary
    assert 'dimBox' not in block
    assert 'pseudoBox' not in block


def test_legacy_dimension_migration_remains_at_load_boundary():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('    def makeinp(self,indir,infile):')
    end = source.index('\n    def ', start + 8)
    block = source[start:end]
    assert 'state.canonicalize_loaded_dimensions(' in block
    assert 'self.dim = state.spectral_dimensions' in block
    assert '.dimBox' not in block
    assert '.pseudoBox' not in block
