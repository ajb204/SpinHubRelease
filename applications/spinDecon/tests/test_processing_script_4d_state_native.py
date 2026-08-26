from pathlib import Path


def _source():
    return (Path(__file__).parents[1] / "processing" / "nmrpipe_scripts.py").read_text()


def _body(name, next_name):
    text = _source()
    start = text.index("def " + name)
    end = text.index("def " + next_name, start)
    return text[start:end]


def test_4d_builder_reads_all_four_spectral_dimensions_from_state():
    body = _body("make_proc_script_4d_state", "write_direct_phase_script")
    assert "frame." not in body.replace("#FidPath=os.path.join(frame.DataStoreBox.GetValue(),frame.FidPathBox.GetValue())", "").replace("#if(frame.cb_basepol3.IsChecked()):", "")
    for name in ("p0", "p1", "p0_1", "p1_1", "p0_2", "p1_2", "p0_3", "p1_3",
                 "windowBox0", "windowBox1", "windowBox2", "windowBox3",
                 "cb_ft0", "cb_ft1", "cb_ft2", "cb_ft3"):
        assert "state.control('%s')" % name in body
    for name in ("cb_baseSol", "cb_baseLin", "cb_basepol", "cb_basepol1", "cb_basepol2"):
        assert "state.checked('%s')" % name in body


def test_4d_nus_schedule_is_local_and_max_iterations_are_snapshotted():
    body = _body("make_proc_script_4d_state", "write_direct_phase_script")
    assert "nuslist = process._current_nus_schedule()" in body
    assert "process.nuslist =" not in body
    assert "state.value('maxIterBox', 0)" in body
    assert "-sample %s/%s" in body


def test_4d_render_write_paths_bypass_legacy_context():
    text = _source()
    render = text[text.index("def render_process_script_state"):text.index("def write_process_script_state")]
    write = text[text.index("def write_process_script_state"):text.index("def render_process_script(")]
    for body in (render, write):
        assert "make_proc_script_4d_state(process, state" in body
        assert "build_processing_script_context" not in body
        assert "if process.spectral_dim_count <= 4:" in body
