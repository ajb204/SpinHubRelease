from pathlib import Path


def test_processing_scripts_do_not_depend_on_wx_error_dialog():
    source = Path(__file__).parents[1].joinpath('processing', 'nmrpipe_scripts.py').read_text()
    assert 'misc.errors' not in source
    assert 'gui.dialogs.errors' not in source
    assert '_report_processing_error' in source
