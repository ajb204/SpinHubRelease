from pathlib import Path

from spinDecon.project.state import ProjectState
from spinDecon.project.parameter_store import parse_value, update_parameter_file

ROOT = Path(__file__).parents[1]


def _body(path, start_marker, end_marker=None):
    text = (ROOT / path).read_text()
    start = text.index(start_marker)
    if end_marker is None:
        return text[start:]
    end = text.index(end_marker, start + len(start_marker))
    return text[start:end]


def test_project_state_keeps_live_unsaved_gui_values():
    state = ProjectState()
    state.update_gui_settings({'p0_1': '12.50', 'ProcTarg': 'SMILE'})
    assert state.gui_value('p0_1') == '12.50'
    assert state.gui_value('ProcTarg') == 'SMILE'
    assert state.dirty is True


def test_parameter_commit_preserves_unrelated_values(tmp_path):
    par = tmp_path / 'decon.par'
    par.write_text('p0_1 = 0.00\nlabel1 = H1\nunrelated = keep\n')
    update_parameter_file(par, {'p0_1': '17.25', 'label1': 'N15'}, source_path=par)
    assert parse_value(par, 'p0_1') == '17.25'
    assert parse_value(par, 'label1') == 'N15'
    assert parse_value(par, 'unrelated') == 'keep'


def test_process_save_collects_all_open_window_widgets_in_one_write():
    body = _body('gui/dialogs/processing/process.py', '    def save_current_gui_state', '    def OnCombineBruker')
    assert 'sync_current_gui_state()' in body
    assert body.count('update_parameter_file(') == 1
    sync = _body('gui/dialogs/processing/process.py', '    def _collect_current_parameter_updates', '    @staticmethod')
    assert "('processing_frame', 'conv_frame', 'projections_frame')" in sync
    assert 'collector()' in sync


def test_processing_run_commits_widgets_before_generating_script():
    body = _body('gui/dialogs/processing/settings.py', '    def on_run(self, event, on_finish=None):', '    def _processing_updates')
    assert body.index("save_current_gui_state(reason='processing-run')") < body.index('_run_generated_script')
    generated = _body('gui/dialogs/processing/settings.py', '    def _run_generated_script', '    def RunProcessScript')
    assert '_generate_and_save_script' in generated


def test_conversion_run_commits_widgets_before_generating_script():
    body = _body('gui/dialogs/processing/conversion.py', '    def on_run(self, event):', '    def _on_close')
    assert body.index("save_current_gui_state(reason='conversion-run')") < body.index('_generate_guess_script')


def test_auto_processing_does_not_reload_stale_disk_over_current_widgets():
    body = _body('gui/dialogs/processing/process.py', '    def _run_processing_auto', '    def OnButtonProcessingAuto')
    assert "save_current_gui_state(reason='processing-auto')" in body
    assert 'reload_from_file=True' not in body
    assert 'loader()' not in body
    assert 'guess(None)' in body


def test_auto_conversion_commits_before_building_vpar_script():
    body = _body('gui/dialogs/processing/process.py', '    def MakeConvScript')
    assert body.index("save_current_gui_state(reason='conversion-auto')") < body.index('conv._build_vpar()')


def test_projection_reprocess_promotes_phase_to_processing_controls_before_run():
    promote = _body('gui/dialogs/processing/projections.py', '    def _save_projection_phase_values', '    def _refresh_projection_after_reprocess')
    assert "updates['p0_1']" in promote
    assert "updates['p1_1']" in promote
    assert 'ctrl.SetValue(str(value))' in promote
    assert 'update_parameter_file(' not in promote
    run = _body('gui/dialogs/processing/projections.py', '    def _run_silent_reprocess', '    def _on_projection_phase_scroll')
    assert run.index('_save_projection_phase_values()') < run.index('_run_processing_auto')


def test_projection_normal_save_excludes_temporary_phase_deltas():
    body = _body('gui/dialogs/processing/projections.py', '    def collect_updates(self, update_state=True):', '    def OnClose')
    for key in ('cmin', 'cthresh', 'cfac', 'cnum'):
        assert repr(key) in body
    assert "'p0_1'" not in body
    assert "'p1_1'" not in body


def test_processing_edits_are_mirrored_to_live_state_immediately():
    init = _body('gui/dialogs/processing/settings.py', '    def __init__', '    def _bind_live_state_controls')
    assert 'self._bind_live_state_controls()' in init
    bind = _body('gui/dialogs/processing/settings.py', '    def _bind_live_state_controls', '    def apply_live_settings')
    assert 'self._processing_updates()' in bind
    assert 'wx.EVT_TEXT' in bind and 'wx.EVT_CHOICE' in bind and 'wx.EVT_CHECKBOX' in bind


def test_processing_script_render_snapshots_widgets_before_rendering():
    body = _body('gui/dialogs/processing/settings.py', '    def _render_script_text', '    def _load_script_text')
    assert body.index('self._processing_updates()') < body.index('self.proc.RenderProcessScript')


def test_conversion_edits_are_mirrored_to_live_state_immediately():
    init = _body('gui/dialogs/processing/conversion.py', '    def __init__', '    def _bind_live_state_controls')
    assert 'self._bind_live_state_controls()' in init
    bind = _body('gui/dialogs/processing/conversion.py', '    def _bind_live_state_controls', '    def apply_live_settings')
    assert 'self.collect_updates()' in bind
    assert 'wx.EVT_TEXT' in bind and 'wx.EVT_CHOICE' in bind and 'wx.EVT_CHECKBOX' in bind


def test_conversion_script_build_snapshots_widgets_before_vpar_setup():
    body = _body('gui/dialogs/processing/conversion.py', '    def _build_vpar', '    def _generate_guess_script')
    assert body.index('self.collect_updates()') < body.index('inst.Setup(')


def test_projection_preview_is_not_persistent_gui_state_until_promoted():
    state = ProjectState()
    state.update_projection_phase('N15', p0=12.5, p1=-3.25)
    assert state.projection_phase('N15') == {'p0': 12.5, 'p1': -3.25}
    assert 'p0_1' not in state.gui_settings
    assert state.dirty is False

    state.promote_projection_phase({'p0_1': '12.50', 'p1_1': '-3.25'})
    assert state.gui_value('p0_1') == '12.50'
    assert state.gui_value('p1_1') == '-3.25'
    assert state.dirty is True


def test_projection_preview_clear_does_not_erase_promoted_processing_phase():
    state = ProjectState()
    state.update_projection_phase('N15', p0=7.0, p1=8.0)
    state.promote_projection_phase({'p0_1': '7.00', 'p1_1': '8.00'})
    state.clear_projection_phase_preview()
    assert state.projection_phase_preview == {}
    assert state.gui_value('p0_1') == '7.00'
    assert state.gui_value('p1_1') == '8.00'


def test_projection_frame_uses_explicit_preview_and_promotion_state_api():
    state_body = _body('gui/dialogs/processing/projections.py', '    def _projection_phase_state', '    def _projection_phase_key')
    assert 'self.state.projection_phase(' in state_body
    assert 'self.state.update_projection_phase(' in state_body
    assert 'self.state.clear_projection_phase_preview()' in state_body
    promote = _body('gui/dialogs/processing/projections.py', '    def _save_projection_phase_values', '    def _refresh_projection_after_reprocess')
    assert 'self.state.promote_projection_phase(updates)' in promote


def test_seed_gui_settings_never_overwrites_newer_live_edits_or_marks_dirty():
    state = ProjectState()
    state.seed_gui_settings({'label1': 'H1', 'p0_1': '10'})
    assert state.gui_settings == {'label1': 'H1', 'p0_1': '10'}
    assert state.dirty is False
    state.update_gui_settings({'p0_1': '37'})
    state.seed_gui_settings({'p0_1': '10', 'label2': 'N15'})
    assert state.gui_value('p0_1') == '37'
    assert state.gui_value('label2') == 'N15'


def test_processing_frame_hydrates_shared_state_before_applying_live_overlay():
    init = _body('gui/dialogs/processing/settings.py', '    def __init__', '    def _bind_live_state_controls')
    assert init.index('self._load_from_file()') < init.index('state.seed_gui_settings(')
    assert init.index('state.seed_gui_settings(') < init.index('self.apply_live_settings()')
    assert '_processing_updates(update_state=False)' in init


def test_conversion_frame_hydrates_shared_state_before_applying_live_overlay():
    init = _body('gui/dialogs/processing/conversion.py', '    def __init__', '    def _bind_live_state_controls')
    assert init.index('self._load_from_file()') < init.index('state.seed_gui_settings(')
    assert init.index('state.seed_gui_settings(') < init.index('self.apply_live_settings()')
    assert 'collect_updates(update_state=False)' in init


def test_process_conversion_accessors_prefer_live_state_over_stale_disk():
    value_body = _body('gui/dialogs/processing/process.py', '    def _conversion_value', '    def _conversion_checked')
    assert "key in state.gui_settings" in value_body
    assert value_body.index('key in state.gui_settings') < value_body.index("self._parse_param(key")
    checked = _body('gui/dialogs/processing/process.py', '    def _conversion_checked', '    def _sync_conversion_dialog')
    assert "key in state.gui_settings" in checked
    assert checked.index('key in state.gui_settings') < checked.index('self._parse_bool(key')


def test_closed_conversion_manual_reference_uses_live_state_before_disk():
    body = _body('gui/dialogs/processing/process.py', '    def _manual_reference_ppm', '    def _conversion_value')
    assert "'xcen' in state.gui_settings" in body
    assert body.index("'xcen' in state.gui_settings") < body.index("self._parse_param('xcen'")


def test_processing_binding_no_longer_aliases_wx_controls_onto_process_frame():
    body = _body('gui/dialogs/processing/process.py', '    def _bind_processing_controls', '    def _unbind_processing_controls')
    assert "setattr(self, name" not in body
    assert 'self.processing_frame = frame' in body


def test_processing_script_generation_uses_explicit_state_boundary():
    body = _body('gui/dialogs/processing/process.py', '    def RenderProcessScript', '    def BuildProcessScript')
    assert 'RenderProcessScriptState(self, self._processing_script_state()' in body
    assert 'WriteProcessScriptState(self, self._processing_script_state()' in body
    assert 'ProcessingScriptContext' not in (ROOT / 'gui/dialogs/processing/process.py').read_text()

def test_processing_control_schema_does_not_depend_on_removed_process_widget_aliases():
    body = (ROOT / 'processing/script_context.py').read_text()
    assert "'p0_1'" in body and "'p0_2'" in body and "'p0_3'" in body
    assert "hasattr(self, 'cb_f1180')" not in body
    assert 'def _processing_control_names' not in (ROOT / 'gui/dialogs/processing/process.py').read_text()


def test_obsolete_processing_control_snapshot_restore_machinery_is_removed():
    text = (ROOT / 'gui/dialogs/processing/process.py').read_text()
    assert '_snapshot_processing_control_state' not in text
    assert '_restore_processing_control_state' not in text
    assert '_with_restored_processing_state' not in text


def test_direct_phase_reads_processing_owner_or_shared_live_state():
    body = _body('gui/dialogs/processing/process.py', '    def _get_direct_phase_values', '    def _saved_direct_phase_values')
    assert "_processing_live_value('p0'" in body
    assert "_processing_live_value('p1'" in body
    assert 'self.p0.GetValue()' not in body
    assert 'self.p1.GetValue()' not in body


def test_autoapodise_updates_processing_owner_without_process_widget_aliases():
    body = _body('gui/dialogs/processing/process.py', '    def OnAutoApodise', '    def _build_direct_frequency_data')
    assert "_set_processing_widget_value('windowBox0', 'GM')" in body
    assert "_set_processing_widget_value('win3Val0'" in body
    assert "_set_processing_widget_value('win2Val0'" in body
    assert 'self.windowBox0.SetValue' not in body


def test_processing_script_state_is_plain_value_snapshot():
    from spinDecon.processing.script_context import ProcessingScriptState
    class C:
        def __init__(self, value): self.value = value
        def GetValue(self): return self.value
    class Frame: pass
    f = Frame(); f.p0_1 = C('37')
    snap = ProcessingScriptState.capture(f, {'p0_1': '10', 'p1_1': '22'}, ('p0_1', 'p1_1'))
    assert snap.value('p0_1') == '37'
    assert snap.value('p1_1') == '22'
    f.p0_1.value = '99'
    assert snap.value('p0_1') == '37'


def test_process_frame_builds_script_state_before_processing_boundary():
    body = _body('gui/dialogs/processing/process.py', '    def _processing_script_state(self):', '    def RenderProcessScript')
    assert 'ProcessingScriptState.capture_current(frame, state)' in body
    render = _body('gui/dialogs/processing/process.py', '    def RenderProcessScript', '    def BuildProcessScript')
    assert 'RenderProcessScriptState(self, self._processing_script_state()' in render

def test_legacy_script_context_is_removed_and_state_is_plain_values():
    text = (ROOT / 'processing/script_context.py').read_text()
    assert 'class ProcessingScriptContext' not in text
    assert 'class FrozenControl' not in text
    assert 'build_processing_script_context' not in text
    assert 'MappingProxyType(values)' in text


def test_processing_script_schema_is_owned_by_processing_layer_not_process_frame():
    process_text = (ROOT / 'gui/dialogs/processing/process.py').read_text()
    context_text = (ROOT / 'processing/script_context.py').read_text()
    assert 'def _processing_control_names' not in process_text
    assert 'PROCESSING_CONTROL_NAMES = (' in context_text
    for name in ('p0', 'p0_1', 'p0_2', 'p0_3', 'cb_lp1', 'windowBox3'):
        assert repr(name) in context_text


def test_process_frame_captures_current_processing_contract_in_one_call():
    body = _body('gui/dialogs/processing/process.py', '    def _processing_script_state(self):', '    def RenderProcessScript')
    assert 'ProcessingScriptState.capture_current(frame, state)' in body
    assert "getattr(state, 'gui_settings'" not in body
    assert '_processing_control_names' not in body

def test_capture_current_prefers_widget_and_falls_back_to_project_live_state():
    from spinDecon.processing.script_context import ProcessingScriptState
    from spinDecon.project.state import ProjectState
    class C:
        def __init__(self, value): self.value = value
        def GetValue(self): return self.value
    class F: pass
    frame = F(); frame.p0_1 = C('37')
    state = ProjectState(); state.update_gui_settings({'p0_1': '10', 'p1_1': '22'})
    snap = ProcessingScriptState.capture_current(frame, state)
    assert snap.value('p0_1') == '37'
    assert snap.value('p1_1') == '22'
    frame.p0_1.value = '99'
    assert snap.value('p0_1') == '37'


def test_parameter_windows_no_longer_expose_independent_save_buttons():
    processing_ui = _body('gui/dialogs/processing/settings.py', '    def _build_ui', '    def _safe_set_text')
    conversion_ui = _body('gui/dialogs/processing/conversion.py', '    def _build_ui', '    def _controls_for_dim')
    projection_ui = _body('gui/dialogs/processing/projections.py', '    def _build_controls', '    def _bind_live_contour_controls')
    assert "wx.Button(panel, label='Save')" not in processing_ui
    assert "wx.Button(panel, label='Save')" not in conversion_ui
    assert "label='Save'" not in projection_ui


def test_script_editor_save_buttons_remain_available():
    processing = (ROOT / 'gui/dialogs/processing/settings.py').read_text()
    conversion = (ROOT / 'gui/dialogs/processing/conversion.py').read_text()
    # Script editor classes retain their own Save controls because they save
    # script text, not the shared GUI/system parameter state.
    assert processing.count("label='Save'") >= 1
    assert conversion.count("label='Save'") >= 1


def test_projection_contours_are_live_state_and_reopen_prefers_live_over_disk():
    bind = _body('gui/dialogs/processing/projections.py', '    def _bind_live_contour_controls', '    def _set_compact_frame_size')
    assert 'self.collect_updates()' in bind
    assert 'wx.EVT_TEXT' in bind
    defaults = _body('gui/dialogs/processing/projections.py', '    def _set_default_values', '    def _projection_phase_state')
    assert 'state.seed_gui_settings(persisted_contours)' in defaults
    assert 'state.resolved_gui_value(key, value)' in defaults


def test_obsolete_child_parameter_save_handlers_are_removed():
    processing = (ROOT / 'gui' / 'dialogs' / 'processing' / 'settings.py').read_text()
    conversion = (ROOT / 'gui' / 'dialogs' / 'processing' / 'conversion.py').read_text()
    projection = (ROOT / 'gui' / 'dialogs' / 'processing' / 'projections.py').read_text()
    process = (ROOT / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
    assert "reason='processing-save'" not in processing
    assert "reason='conversion-save'" not in conversion
    assert "reason='projection-save'" not in projection
    assert "reason='process-save'" not in process


def test_process_window_no_longer_has_parameter_save_button():
    body = _body('gui/dialogs/processing/process.py', '    def button_box(self):', '    def OnButtonConversion')
    assert "label='Save'" not in body
    assert 'self.buttonSave' not in body


def test_process_close_prompts_only_for_changed_persistable_parameters():
    body = _body('gui/dialogs/processing/process.py', '    def OnClose(self, event):', '    def IsInParse')
    assert 'parameters_changed_since_baseline()' in body
    assert 'Parameters have been changed, would you like to save?' in body
    assert "save_current_gui_state(reason='process-close')" in body
    assert 'discard_parameter_changes()' in body


def test_process_dirty_comparison_is_side_effect_free_and_projection_preview_excluded():
    collect = _body('gui/dialogs/processing/process.py', '    def _collect_current_parameter_updates', '    @staticmethod')
    assert 'saved_gui' in collect and 'saved_metadata' in collect and 'saved_dirty' in collect
    assert 'state.gui_settings.update(saved_gui or {})' in collect
    projection = _body('gui/dialogs/processing/projections.py', '    def collect_updates(self, update_state=True):', '    def OnClose')
    assert 'projection_phase_preview' not in projection


def test_successful_process_save_refreshes_close_baseline():
    body = _body('gui/dialogs/processing/process.py', '    def save_current_gui_state', '    def OnCombineBruker')
    assert 'self._capture_parameter_baseline(updates)' in body
    assert 'update_parameter_file(savefile, updates' in body


def test_discard_restores_live_process_state_without_writing_disk():
    body = _body('gui/dialogs/processing/process.py', '    def discard_parameter_changes', '    def sync_current_gui_state')
    assert 'state.gui_settings.clear()' in body
    assert "state.gui_settings.update(copy.deepcopy(getattr(self, '_gui_settings_baseline', {})))" in body
    assert 'state.projection_phase_preview = {}' in body
    assert 'update_parameter_file' not in body


def test_file_menu_save_routes_to_focused_process_family_after_project_save():
    body = _body('app/notebook.py', '    def OnSaveResults(self, event):', '    def OnLoadResults')
    assert 'active_process = self._focused_process_frame()' in body
    assert body.index('self.notebook.nmr_workspace.OnButtonSave(True)') < body.index("active_process.save_current_gui_state(reason='file-menu-save')")
    focus = _body('app/notebook.py', '    def _focused_process_frame(self):', '    def OnSaveResults')
    assert 'wx.Window.FindFocus()' in focus
    assert "window.__class__.__name__ == 'ProcessFrame'" in focus


def test_process_session_hydration_can_replace_stale_defaults_without_dirtying_state():
    state = ProjectState()
    state.gui_settings.update({'label1': 'H', 'label2': 'H_2'})
    state.dirty = False

    state.hydrate_gui_settings({'label1': 'H', 'label2': 'N'}, overwrite=True)

    assert state.gui_settings['label1'] == 'H'
    assert state.gui_settings['label2'] == 'N'
    assert state.dirty is False


def test_child_hydration_still_preserves_newer_live_labels():
    state = ProjectState()
    state.gui_settings.update({'label1': 'H', 'label2': 'N'})

    state.seed_gui_settings({'label1': 'H', 'label2': 'C'})

    assert state.gui_settings['label2'] == 'N'


def test_process_hydrates_labels_before_getlabs_and_projection_uses_process_labels():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    process_source = (root / 'gui' / 'dialogs' / 'processing' / 'process.py').read_text()
    projection_source = (root / 'gui' / 'dialogs' / 'processing' / 'projections.py').read_text()

    hydrate_pos = process_source.index('self._hydrate_shared_dimension_labels()')
    getlabs_pos = process_source.index('self.GetLabs()', hydrate_pos)
    assert hydrate_pos < getlabs_pos
    assert 'def get_spectral_labels(self):' in process_source
    assert "getter = getattr(self.process_parent, 'get_spectral_labels', None)" in projection_source

def test_process_save_persists_current_nmrpipe_input_file():
    body = _body('gui/dialogs/processing/process.py', '    def collect_updates(self):', '    def _collect_current_parameter_updates')
    assert "'infile': self._current_nmrpipe_input_file()" in body


def test_processing_completion_refreshes_project_lamp_and_workflow():
    body = _body('gui/dialogs/processing/process.py', '    def _update_nmrpipe_file_box', '    def _conversion_outputs_exist')
    assert 'update_project_lamps' in body
    assert 'notify_analysis_changed' in body
