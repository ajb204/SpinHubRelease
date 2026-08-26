from pathlib import Path


def test_load_project_syncs_paths_but_does_not_clobber_topology_before_read():
    source = (Path(__file__).parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = source.index('    def OnButtonLoadProject(')
    end = source.index('    def OnButtonHousekeeping(', start)
    body = source[start:end]
    sync = body.index('self.state.sync_from_values(')
    read = body.index('self.OnButtonRead(True)')
    assert sync < read
    assert 'working_dir=self.dirBox.GetValue()' in body
    assert 'input_file=self.infileBox.GetValue()' in body
    assert 'spectral_dimensions=' not in body
    assert 'dimension=' not in body
    assert 'pseudo_axis=' not in body
