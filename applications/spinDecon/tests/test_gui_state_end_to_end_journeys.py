"""End-to-end regression journeys for the GUI -> state -> disk -> script contract.

These tests deliberately make persisted, shared-state, and widget values disagree.
The visible/current widget must win at the script boundary.
"""
from pathlib import Path
from types import SimpleNamespace
import pytest

from spinDecon.project.parameter_store import parse_value, update_parameter_file
from spinDecon.project.state import ProjectState
from spinDecon.processing.script_context import ProcessingScriptState


class TextControl:
    def __init__(self, value): self.value = value
    def GetValue(self): return self.value
    def SetValue(self, value): self.value = value


class CheckControl:
    def __init__(self, value): self.value = bool(value)
    def IsChecked(self): return self.value


class ChoiceControl:
    def __init__(self, value): self.value = int(value)
    def GetSelection(self): return self.value


@pytest.mark.parametrize(
    "spectral_dim,pseudo,builder_name,phase_key",
    [
        (1, False, "make_proc_script_1d_state", "p0"),
        (2, False, "make_proc_script_2d_state", "p0_1"),
        (2, True, "make_proc_script_3dp_state", "p0_1"),
        (3, False, "make_proc_script_3d_state", "p0_2"),
        (4, False, "make_proc_script_4d_state", "p0_3"),
    ],
)
def test_current_processing_widget_wins_all_script_routes(
    tmp_path, spectral_dim, pseudo, builder_name, phase_key
):
    """Stale disk/shared values lose to the visible widget before rendering."""
    par = tmp_path / "system.par"
    par.write_text(f"{phase_key} = 10\n")
    state = ProjectState(parameter_file=str(par))
    state.seed_gui_settings({phase_key: parse_value(par, phase_key)})
    state.update_gui_settings({phase_key: "20"})
    processing = SimpleNamespace(**{phase_key: TextControl("37")})
    snapshot = ProcessingScriptState.capture_current(processing, state)
    assert snapshot.value(phase_key) == "37"
    state.update_gui_settings({phase_key: snapshot.value(phase_key)})
    update_parameter_file(par, {phase_key: state.gui_value(phase_key)}, source_path=par)
    assert parse_value(par, phase_key) == "37"

    # Every dimensional dispatch path is explicitly state-native.
    source = (Path(__file__).parents[1] / "processing" / "nmrpipe_scripts.py").read_text()
    assert f"{builder_name}(process, state" in source


def test_processing_snapshot_is_stable_during_script_generation():
    """A later widget edit cannot mutate an in-flight processing script."""
    state = ProjectState()
    control = TextControl("37")
    processing = SimpleNamespace(p0_1=control)
    snapshot = ProcessingScriptState.capture_current(processing, state)
    control.SetValue("99")
    assert snapshot.value("p0_1") == "37"


def test_projection_preview_only_reaches_processing_after_reprocess_promotion(tmp_path):
    par = tmp_path / 'system.par'
    par.write_text('p0_1 = 10\np1_1 = 11\n')
    state = ProjectState(parameter_file=str(par))
    state.seed_gui_settings({'p0_1': '10', 'p1_1': '11'})

    state.update_projection_phase('N15', p0=37, p1=-4)
    # Ordinary save must still see the main processing phase, not preview.
    update_parameter_file(par, state.gui_settings, source_path=par)
    assert parse_value(par, 'p0_1') == '10'

    # Re-process is the explicit promotion boundary.
    preview = state.projection_phase('N15')
    promoted = state.promote_projection_phase({
        'p0_1': f"{preview['p0']:.2f}",
        'p1_1': f"{preview['p1']:.2f}",
    })
    update_parameter_file(par, promoted, source_path=par)
    assert parse_value(par, 'p0_1') == '37.00'
    assert parse_value(par, 'p1_1') == '-4.00'

    snapshot = ProcessingScriptState.capture_current(None, state)
    assert snapshot.value('p0_1') == '37.00'
    assert snapshot.value('p1_1') == '-4.00'


def test_close_reopen_overlay_keeps_unsaved_live_value_until_explicit_reload():
    """Frame hydration may fill gaps but must not overwrite newer live edits."""
    state = ProjectState()
    state.seed_gui_settings({'p0_1': '10', 'label1': 'H1'})
    state.update_gui_settings({'p0_1': '37'})

    # Simulate reopening a frame which first reads the stale parameter file.
    state.seed_gui_settings({'p0_1': '10', 'label1': 'H1', 'label2': 'N15'})
    assert state.gui_value('p0_1') == '37'
    assert state.gui_value('label2') == 'N15'

    # A genuine project reload is a new state boundary and may load disk again.
    reloaded = ProjectState()
    reloaded.seed_gui_settings({'p0_1': '10', 'label1': 'H1', 'label2': 'N15'})
    assert reloaded.gui_value('p0_1') == '10'


def test_conversion_current_widgets_are_collected_before_building_vpar():
    """Conversion snapshots controls before vpar Setup/script inference."""
    source = (Path(__file__).parents[1] / "gui" / "dialogs" / "processing" / "conversion.py").read_text()
    start = source.index("    def _build_vpar(self):")
    end = source.index("    def _generate_guess_script", start)
    body = source[start:end]
    assert body.index("self.collect_updates()") < body.index("inst.Setup(")
    assert "self.label" not in body  # labels are dimension-indexed controls below
    assert "getattr(self, f'label{i}').GetValue()" in body
    assert "nuslist=self.nusFil.GetValue().strip()" in body
    assert "o1p = self._selected_reference()" in body
