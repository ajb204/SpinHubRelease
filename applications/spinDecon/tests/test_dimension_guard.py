import numpy as np
import pytest

from spinDecon.domain.topology import DatasetTopology
from spinDecon.domain.dimensions.guard import assert_full_dataset_contract


def test_full_dataset_guard_accepts_physical_dimension_count():
    topology = DatasetTopology.from_counts(2, True)
    assert assert_full_dataset_contract(topology, np.zeros((2, 3, 4))) is topology


def test_full_dataset_guard_rejects_spectral_count_used_as_physical_count():
    topology = DatasetTopology.from_counts(2, True)
    with pytest.raises(ValueError, match="physical_dim_count=3"):
        assert_full_dataset_contract(topology, np.zeros((3, 4)), where="test")
