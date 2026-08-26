from spinDecon.project.parameter_store import update_parameter_file, remove_parameter_keys
from spinDecon.project.defaults import UNIDEC_DEFAULTS, is_default_value
from spinDecon.project.service import ProjectService


def test_unidec_canonical_defaults():
    assert UNIDEC_DEFAULTS == {"thresh": 0.08, "fac": 1.4, "conv": 1e-7, "maxiter": 100}
    assert is_default_value("fac", "1.4")
    assert not is_default_value("fac", "2.0")


def test_default_override_can_be_removed(tmp_path):
    par = tmp_path / "spinHub.par"
    update_parameter_file(par, {"dim": 2, "fac": 2.0, "conv": 1e-7})
    remove_parameter_keys(par, {"fac", "conv"})
    text = par.read_text()
    assert "dim = 2" in text
    assert "fac" not in text
    assert "conv" not in text


def test_new_project_does_not_persist_default_unidec_settings(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    state = ProjectService().create_initial_parameter_file(tmp_path, raw, dimension=2, pseudo_axis=True)
    text = (tmp_path / "spinHub.par").read_text()
    for key in ("fac", "conv", "maxiter", "thresh", "ncpus"):
        assert not any(line.startswith(key + " =") for line in text.splitlines())
