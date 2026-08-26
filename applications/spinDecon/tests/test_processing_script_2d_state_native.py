from pathlib import Path


def _source():
    return (Path(__file__).parents[1] / 'processing' / 'nmrpipe_scripts.py').read_text()


def _builder_body():
    text = _source()
    start = text.index('def make_proc_script_2d_state')
    end = text.index('def make_proc_script_2d(', start)
    return text[start:end]


def test_2d_native_builder_reads_processing_values_from_state():
    body = _builder_body()
    assert 'frame.' not in body
    for name in ('p0', 'p1', 'p0_1', 'p1_1', 'windowBox0', 'windowBox1', 'cb_ft0', 'cb_ft1'):
        assert "state.control('%s')" % name in body
    assert "state.value('maxIterBox', 0)" in body
    for name in ('cb_baseSol', 'cb_baseLin', 'cb_basepol', 'cb_basepol1', 'cb_lp1'):
        assert "state.checked('%s')" % name in body
    assert 'process.maxIterBox' not in body
    assert '.GetValue()' not in body


def test_2d_render_and_write_paths_bypass_legacy_context():
    text = _source()
    render = text[text.index('def render_process_script_state'):text.index('def write_process_script_state')]
    write = text[text.index('def write_process_script_state'):text.index('def render_process_script(')]
    assert 'make_proc_script_2d_state(process, state' in render
    assert 'make_proc_script_2d_state(process, state' in write
    assert 'if process.spectral_dim_count <= 4:' in render
    assert 'if process.spectral_dim_count <= 4:' in write


def test_2d_nus_schedule_is_local_not_written_to_process():
    body = _builder_body()
    assert 'nuslist = process._current_nus_schedule()' in body
    assert 'process.nuslist =' not in body
