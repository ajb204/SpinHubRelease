"""Canonical dimension/axis access for spectrum viewers (Stage 6).

Viewers must not infer scientific dimensionality from ``data.ndim`` or from
pseudo-specific arithmetic.  This module gives them one route to the canonical
DatasetTopology owned by ProjectState and explicit physical/spectral axes.
"""
from __future__ import annotations

from ..topology import DatasetTopology


def topology_for(tab) -> DatasetTopology:
    state = getattr(tab, "state", None)
    if state is not None and getattr(state, "spectral_dimensions", 0):
        return state.topology()
    spectral = max(1, int(getattr(tab, "dim", 1)))
    pseudo = bool(getattr(tab, "pseudo", False))
    return DatasetTopology.from_counts(spectral, pseudo)


def spectral_dim_count(tab) -> int:
    return topology_for(tab).spectral_dim_count


def physical_dim_count(tab) -> int:
    return topology_for(tab).physical_dim_count


def spectral_physical_indices(tab):
    return tuple(a.physical_index for a in topology_for(tab).spectral_axes)


def pseudo_physical_index(tab):
    axis = topology_for(tab).pseudo_axis
    return None if axis is None else axis.physical_index
