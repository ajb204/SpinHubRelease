"""Dimension contracts shared by Process, Conversion and Processing windows.

GUI dimensionality is always spectral.  Conversion consumes physical axes;
Processing consumes spectral axes.  The legacy vpar encoding is isolated here
until that backend is migrated to DatasetTopology.
"""
from __future__ import annotations

from spinDecon.domain.topology import DatasetTopology


def conversion_axis_count(topology: DatasetTopology) -> int:
    """Number of conversion definitions: one per physical data axis."""
    return topology.physical_dim_count


def processing_axis_count(topology: DatasetTopology) -> int:
    """Number of processing definitions: one per spectral axis."""
    return topology.spectral_dim_count


def legacy_vpar_dimension(topology: DatasetTopology):
    """Temporary adapter for the un-migrated vendor conversion backend.

    This value must never be exposed as GUI dimensionality.
    """
    if topology.has_pseudo_axis and topology.spectral_dim_count in (1, 2):
        return f"{topology.physical_dim_count}p"
    return topology.spectral_dim_count
