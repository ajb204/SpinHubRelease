"""Stage 2 tests: ProjectState is the canonical load-time topology boundary."""
from pathlib import Path

import pytest

from spinDecon.project.state import ProjectState

REAL = ("time_T2", "ID", "ncyc", "ncyc_cp", "gzlvl5", "gzlvl1")


def test_dimension_compatibility_field_now_means_spectral_dimensions():
    state = ProjectState(dimension=2, pseudo_axis=True)
    assert state.spectral_dimensions == 2
    assert state.physical_dimensions == 3
    assert state.topology().spectral_dim_count == 2
    assert state.topology().physical_dim_count == 3


def test_legacy_physical_3p_is_migrated_once_at_load_boundary():
    state = ProjectState(dimension=3, pseudo_axis=True)
    migrated = state.canonicalize_loaded_dimensions(
        3, ("time_T2", "15N", "1H"), real_axis_labels=REAL
    )
    assert migrated is True
    assert state.dimension == 2
    assert state.spectral_dimensions == 2
    assert state.physical_dimensions == 3
    assert state.metadata["legacy_dimension_migrated_from"] == 3
    assert state.metadata["dimension_contract"] == "spectral"


def test_canonical_pseudo_state_is_not_reinterpreted():
    state = ProjectState(dimension=2, pseudo_axis=True)
    migrated = state.canonicalize_loaded_dimensions(
        3, ("time_T2", "15N", "1H"), real_axis_labels=REAL
    )
    assert migrated is False
    assert state.dimension == 2


def test_nonpseudo_state_requires_physical_and_spectral_counts_to_match():
    state = ProjectState(dimension=2, pseudo_axis=False)
    state.canonicalize_loaded_dimensions(2, ("15N", "1H"), real_axis_labels=REAL)
    with pytest.raises(ValueError, match="physical dimensions"):
        state.canonicalize_loaded_dimensions(3, ("13C", "15N", "1H"), real_axis_labels=REAL)


def test_ambiguous_pseudo_count_is_not_guessed_without_real_axis_evidence():
    state = ProjectState(dimension=3, pseudo_axis=True)
    with pytest.raises(ValueError, match="physical dimensions"):
        state.canonicalize_loaded_dimensions(
            3, ("13C", "15N", "1H"), real_axis_labels=REAL
        )


def test_parameter_file_reads_pseudo_and_marks_legacy_pseudo_contract_unresolved(tmp_path):
    par = tmp_path / "project.txt"
    par.write_text("dim = 3\npseudo = 1\n")
    state = ProjectState.from_parameter_file(par)
    assert state.dimension == 3
    assert state.pseudo_axis is True
    assert state.metadata["dimension_contract"] == "legacy_unresolved"


def test_loaded_axis_identity_is_preserved_in_canonical_topology():
    state = ProjectState(dimension=2, pseudo_axis=True)
    state.canonicalize_loaded_dimensions(
        3, ("15N", "time_T2", "1H"), real_axis_labels=REAL
    )
    topology = state.topology()
    assert topology.pseudo_axis.physical_index == 1
    assert topology.pseudo_axis.label == "time_T2"
    assert tuple(a.physical_index for a in topology.spectral_axes) == (0, 2)
    assert tuple(a.label for a in topology.spectral_axes) == ("15N", "1H")


def test_legacy_unresolved_pseudo2d_migrates_without_axis_label_guessing():
    state = ProjectState(dimension=2, pseudo_axis=True)
    state.metadata["dimension_contract"] = "legacy_unresolved"
    migrated = state.canonicalize_loaded_dimensions(
        2, ("unknown_real_axis", "1H"), real_axis_labels=REAL
    )
    assert migrated is True
    assert state.spectral_dimensions == 1
    assert state.physical_dimensions == 2
    assert state.topology().spectral_dim_count == 1
    assert state.topology().has_pseudo_axis is True


def test_pseudo_physical_index_compatibility_accessor():
    state = ProjectState(dimension=1, pseudo_axis=True)
    state.metadata['pseudo_physical_index'] = 1
    assert state.topology().pseudo_physical_index == 1


def test_usta_label_identifies_pseudo_axis_case_insensitively():
    state = ProjectState(dimension=1, pseudo_axis=True)
    state.metadata['dimension_contract'] = 'legacy_unresolved'
    migrated = state.canonicalize_loaded_dimensions(
        2, ('uSTA', '1H'), real_axis_labels=('usta',)
    )
    assert migrated is False
    assert state.spectral_dimensions == 1
    assert state.metadata['pseudo_physical_index'] == 0
    assert state.metadata['spectral_axis_labels'] == ('1H',)
