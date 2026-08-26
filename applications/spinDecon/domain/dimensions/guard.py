"""Runtime assertions for the canonical dataset dimensionality contract.

Use these at boundaries where a full physical spectrum is handed between
subsystems.  Derived projections/slices intentionally validate their own array
shape instead because their ndarray dimensionality no longer describes the
source experiment topology.
"""
from __future__ import annotations

from ..topology import DatasetTopology


def assert_full_dataset_contract(topology: DatasetTopology, data, *, where="dataset"):
    if data is None:
        raise ValueError(f"{where}: full dataset is missing")
    ndim = getattr(data, "ndim", None)
    if ndim is None:
        ndim = len(getattr(data, "shape", ()))
    try:
        topology.validate_data_ndim(ndim)
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    if len(topology.spectral_axes) != topology.spectral_dim_count:
        raise AssertionError(f"{where}: spectral axis count violates topology")
    if len(topology.pseudo_axes) != topology.pseudo_dim_count:
        raise AssertionError(f"{where}: pseudo axis count violates topology")
    return topology
