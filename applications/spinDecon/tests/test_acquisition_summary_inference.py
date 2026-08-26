from spinDecon.processing.vpar_decon import inspect_acquisition, find_child_acquisitions


def _jcamp(path, **values):
    path.write_text('\n'.join(f'##${k}= {v}' for k, v in values.items()) + '\n')


def test_bruker_ignores_single_point_indirect_dimension(tmp_path):
    (tmp_path / 'ser').write_bytes(b'raw')
    _jcamp(tmp_path / 'acqus', PULPROG='<hsqc>', NUC1='<1H>')
    _jcamp(tmp_path / 'acqu2s', TD='128', NUC1='<15N>')
    _jcamp(tmp_path / 'acqu3s', TD='1', NUC1='<13C>')
    info = inspect_acquisition(tmp_path)
    assert info.dimension == 2
    assert info.sequence == 'hsqc'


def test_parent_folder_reports_child_acquisitions(tmp_path):
    child = tmp_path / '1'; child.mkdir()
    (child / 'fid').write_bytes(b'raw')
    _jcamp(child / 'acqus', PULPROG='<zg>', TD='1024')
    found = find_child_acquisitions(tmp_path)
    assert [x.path.name for x in found] == ['1']
