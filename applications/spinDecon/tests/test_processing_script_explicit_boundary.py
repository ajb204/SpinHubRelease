from pathlib import Path

from spinDecon.processing.script_context import ProcessingScriptState


class Control:
    def __init__(self, value): self.value = value
    def GetValue(self): return self.value


class Processing:
    p0_1 = Control('37')


class ProjectState:
    gui_settings = {'p0_1': '10'}


def test_processing_layer_captures_explicit_plain_state():
    state = ProcessingScriptState.capture(Processing(), ProjectState.gui_settings, ['p0_1'])
    assert state.value('p0_1') == '37'
    assert state.control('p0_1') == '37'


def test_legacy_processing_context_has_been_removed():
    root = Path(__file__).resolve().parents[1]
    source = (root / 'processing' / 'script_context.py').read_text()
    assert 'class FrozenControl' not in source
    assert 'class ProcessingScriptContext' not in source
    assert 'build_processing_script_context' not in source


def test_process_frame_uses_explicit_processing_script_state():
    root = Path(__file__).resolve().parents[1]
    source = (root / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
    assert 'ProcessingScriptContext' not in source
    assert '_processing_script_context' not in source
    assert 'RenderProcessScriptState(self, self._processing_script_state()' in source
    assert 'WriteProcessScriptState(self, self._processing_script_state()' in source


def test_nmrpipe_exposes_explicit_state_api_without_context_adapter():
    root = Path(__file__).resolve().parents[1]
    source = (root / 'processing' / 'nmrpipe_scripts.py').read_text()
    assert 'def render_process_script_state(process, state: ProcessingScriptState' in source
    assert 'def write_process_script_state(process, state: ProcessingScriptState' in source
    assert 'ProcessingScriptContext' not in source
    assert 'build_processing_script_context' not in source
