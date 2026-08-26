from pathlib import Path

from spinDecon.processing.script_context import ProcessingScriptState

ROOT = Path(__file__).resolve().parents[1]


def test_processing_script_state_exposes_plain_checked_and_selection_accessors():
    state = ProcessingScriptState.capture(None, {'cb_baseSol': True, 'p0': '37'})
    assert state.checked('cb_baseSol') is True
    assert state.value('p0') == '37'
    assert state.checked('missing') is False
    assert state.selection('missing', 3) == 3


def test_1d_builder_is_state_native_and_render_bypasses_legacy_context():
    source = (ROOT / 'processing/nmrpipe_scripts.py').read_text()
    body = source[source.index('def make_proc_script_1d_state'):source.index('def make_proc_script_1d(frame')]
    assert "state.control('p0')" in body
    assert "state.control('windowBox0')" in body
    assert "state.checked('cb_baseSol')" in body
    for name in ('p0', 'p1', 'windowBox0', 'win2Val0', 'win3Val0', 'firstPoint0', 'cb_ft0', 'cb_baseSol', 'cb_baseLin', 'cb_basepol'):
        assert f'process.{name}' not in body
        assert f'frame.{name}' not in body
    render = source[source.index('def render_process_script_state'):source.index('def write_process_script_state')]
    assert 'if process.spectral_dim_count == 1:' in render
    assert 'make_proc_script_1d_state(process, state' in render
    assert 'build_processing_script_context(process, state)' not in render
