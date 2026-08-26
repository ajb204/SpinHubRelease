from pathlib import Path
from spinDecon.project.service import ProjectService


def test_discovery_supports_legacy_and_new_with_legacy_precedence(tmp_path):
    service = ProjectService()
    assert service.discover_parameter_file(tmp_path) is None
    new = tmp_path / 'spinHub.par'; new.write_text('dim = 1\n')
    assert service.discover_parameter_file(tmp_path) == new
    legacy = tmp_path / 'deconParFile'; legacy.write_text('dim = 2\n')
    assert service.discover_parameter_file(tmp_path) == legacy


def test_create_defaults_to_spinhub_and_writes_setup_values(tmp_path):
    raw = tmp_path / 'raw'; raw.mkdir()
    state = ProjectService().create_initial_parameter_file(
        tmp_path, raw, dimension=2, pseudo_axis=True
    )
    par = tmp_path / 'spinHub.par'
    assert par.is_file()
    text = par.read_text()
    assert 'specPath = ./spec' in text
    assert 'dim = 2' in text
    assert 'pseudo = 1' in text
    assert state.spectral_dimensions == 2
    assert state.pseudo_axis is True
    assert (tmp_path / 'spec').is_dir()


def test_create_project_uses_new_default_filename(tmp_path):
    raw = tmp_path / 'raw'; raw.mkdir()
    state = ProjectService().create(raw, tmp_path / 'project', dimension=1)
    assert Path(state.parameter_file).name == 'spinHub.par'
