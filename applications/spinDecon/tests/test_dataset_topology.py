"""Contract tests for the canonical dataset dimensionality model."""
import pytest

from spinDecon.domain.topology import AxisKind, AxisSpec, DatasetTopology


@pytest.mark.parametrize(
    "spectral,pseudo,physical",
    [
        (1, False, 1),
        (1, True, 2),
        (2, False, 2),
        (2, True, 3),
        (3, False, 3),
        (3, True, 4),
        (4, False, 4),
        (4, True, 5),
    ],
)
def test_supported_topologies_have_one_unambiguous_count_contract(
    spectral, pseudo, physical
):
    topology = DatasetTopology.from_counts(spectral, pseudo)

    assert topology.spectral_dim_count == spectral
    assert topology.has_pseudo_axis is pseudo
    assert topology.pseudo_dim_count == int(pseudo)
    assert topology.physical_dim_count == physical
    assert len(topology.axes) == physical
    assert len(topology.spectral_axes) == spectral
    assert len(topology.pseudo_axes) == int(pseudo)


def test_physical_and_spectral_axis_indices_are_distinct_concepts():
    topology = DatasetTopology(
        spectral_dim_count=2,
        has_pseudo_axis=True,
        axes=(
            AxisSpec(0, AxisKind.PSEUDO_REAL, label="time"),
            AxisSpec(1, AxisKind.SPECTRAL, label="15N", spectral_index=0),
            AxisSpec(2, AxisKind.SPECTRAL, label="1H", spectral_index=1),
        ),
    )

    assert topology.pseudo_axis.physical_index == 0
    assert topology.spectral_axes[0].physical_index == 1
    assert topology.spectral_axes[0].spectral_index == 0
    assert [axis.label for axis in topology.physical_axes] == ["time", "15N", "1H"]


def test_pseudo_axis_can_appear_at_any_physical_position():
    topology = DatasetTopology.from_counts(
        2,
        True,
        pseudo_physical_index=1,
        spectral_labels=("15N", "1H"),
        pseudo_label="delay",
    )

    assert [axis.kind for axis in topology.physical_axes] == [
        AxisKind.SPECTRAL,
        AxisKind.PSEUDO_REAL,
        AxisKind.SPECTRAL,
    ]
    assert [axis.physical_index for axis in topology.spectral_axes] == [0, 2]
    assert topology.pseudo_axis.label == "delay"


def test_data_ndim_validation_uses_physical_dimension_count():
    topology = DatasetTopology.from_counts(2, True)
    topology.validate_data_ndim(3)

    with pytest.raises(ValueError, match="data.ndim=2"):
        topology.validate_data_ndim(2)


@pytest.mark.parametrize("spectral", [0, 5])
def test_spectral_dimension_count_is_limited_to_one_through_four(spectral):
    with pytest.raises(ValueError, match="spectral_dim_count"):
        DatasetTopology.from_counts(spectral)


def test_rejects_axis_count_that_disagrees_with_topology():
    with pytest.raises(ValueError, match="physical_dim_count"):
        DatasetTopology(
            spectral_dim_count=2,
            has_pseudo_axis=True,
            axes=(
                AxisSpec(0, AxisKind.SPECTRAL, spectral_index=0),
                AxisSpec(1, AxisKind.SPECTRAL, spectral_index=1),
            ),
        )


def test_rejects_duplicate_or_noncontiguous_physical_indices():
    with pytest.raises(ValueError, match="physical axis indices"):
        DatasetTopology(
            spectral_dim_count=2,
            has_pseudo_axis=False,
            axes=(
                AxisSpec(0, AxisKind.SPECTRAL, spectral_index=0),
                AxisSpec(0, AxisKind.SPECTRAL, spectral_index=1),
            ),
        )


def test_rejects_pseudo_axis_with_spectral_index():
    with pytest.raises(ValueError, match="cannot have a spectral_index"):
        AxisSpec(0, AxisKind.PSEUDO_REAL, spectral_index=0)
