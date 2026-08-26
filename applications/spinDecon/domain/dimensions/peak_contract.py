"""Peak dimensionality rules built on the canonical dataset topology.

Peak coordinates are spectral coordinates only.  A sampled real/pseudo axis
identifies a series observation/plane and never adds a coordinate column to a
peak list.
"""
from __future__ import annotations

from typing import Iterable, Tuple

from ..topology import DatasetTopology


def peak_coordinate_count(topology: DatasetTopology) -> int:
    return topology.spectral_dim_count


def peak_list_suffix(topology: DatasetTopology) -> str:
    return f'.{topology.spectral_dim_count}D.list'


def spectral_physical_indices(topology: DatasetTopology) -> Tuple[int, ...]:
    return tuple(axis.physical_index for axis in topology.spectral_axes)


def spectral_labels(topology: DatasetTopology) -> Tuple[str, ...]:
    return tuple(axis.label for axis in topology.spectral_axes)
