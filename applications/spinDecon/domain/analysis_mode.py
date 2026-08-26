"""Compatibility-first interpretation of the legacy GUI dimensionality state.

This module deliberately does not change ProjectState, deconFrame, or any GUI
control.  It provides a read-only description that future workflow UI code can
use without duplicating the legacy ``dim``/``pseudo`` rules.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


class WorkflowKind(str, Enum):
    """Top-level scientific objective selected from dataset topology."""

    SPECTRAL_PEAK_LIST = "spectral_peak_list"
    PSEUDO_AXIS_SERIES = "pseudo_axis_series"


@dataclass(frozen=True)
class AnalysisMode:
    """Read-only, normalized description of an NMR dataset.

    ``legacy_dimension`` is the value currently used by ``deconFrame.dim`` and
    the dimension radio box.  For historical pseudo datasets this value counts
    spectral dimensions and the array has one additional physical pseudo axis.
    From Stage 3 onward the GUI dimension value has exactly one meaning:
    the number of spectral dimensions.  A pseudo axis, when present, is one
    additional physical dimension.
    """

    legacy_dimension: int
    physical_dimensions: int
    spectral_dimensions: int
    pseudo_dimensions: int
    has_pseudo_axis: bool
    workflow_kind: WorkflowKind

    @classmethod
    def from_legacy(
        cls,
        dimension: int,
        pseudo_axis: bool = False,
    ) -> "AnalysisMode":
        """Interpret the existing GUI state without modifying it.

        ``dimension`` is canonical spectral dimensionality. Legacy projects
        that persisted physical dimensionality are normalized at the project
        load boundary before they reach this class.
        """
        dim = int(dimension or 0)
        if dim < 1 or dim > 4:
            raise ValueError("dimension must be an integer from 1 to 4")

        pseudo = bool(pseudo_axis)
        if pseudo and dim == 4:
            # Canonical topology rule: pseudo4D is not supported.
            raise ValueError("pseudo-axis mode is not supported for legacy dimension 4")

        if not pseudo:
            return cls(
                legacy_dimension=dim,
                physical_dimensions=dim,
                spectral_dimensions=dim,
                pseudo_dimensions=0,
                has_pseudo_axis=False,
                workflow_kind=WorkflowKind.SPECTRAL_PEAK_LIST,
            )

        # Canonical rule: GUI/project dimension=N spectral dimensions while
        # a pseudo dataset contains one additional physical real axis.
        physical = dim + 1
        spectral = dim

        return cls(
            legacy_dimension=dim,
            physical_dimensions=physical,
            spectral_dimensions=spectral,
            pseudo_dimensions=1,
            has_pseudo_axis=True,
            workflow_kind=WorkflowKind.PSEUDO_AXIS_SERIES,
        )

    @classmethod
    def from_project_state(
        cls,
        state: Any,
    ) -> "AnalysisMode":
        """Build from ProjectState-like objects without creating a dependency."""
        return cls.from_legacy(
            getattr(state, "dimension", 0),
            getattr(state, "pseudo_axis", False),
        )
