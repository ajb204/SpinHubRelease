"""Canonical spectrum/peak-list policy for the supported analysis journeys.

This module turns dataset topology into architectural decisions.  It is pure
model code: no GUI, filesystem, or scientific calculation lives here.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .analysis_mode import AnalysisMode


class ReferencePolicy(str, Enum):
    NONE = "none"
    FULL_IS_REFERENCE = "full_is_reference"
    INDEPENDENT_PROJECTION = "independent_projection"


@dataclass(frozen=True)
class SpectrumPolicy:
    journey: str
    full_dimension: int
    reference_policy: ReferencePolicy
    reference_dimension: int | None
    uses_pseudo_axis: bool

    @property
    def has_independent_reference(self) -> bool:
        return self.reference_policy is ReferencePolicy.INDEPENDENT_PROJECTION

    @property
    def full_is_reference(self) -> bool:
        return self.reference_policy is ReferencePolicy.FULL_IS_REFERENCE

    @property
    def has_distinct_reference_peak_list(self) -> bool:
        """Whether the journey exposes a separately owned Reference peak list."""
        return self.reference_policy is ReferencePolicy.INDEPENDENT_PROJECTION

    @property
    def projection_peak_list_key(self) -> str | None:
        """Authoritative peak-list key used by the projection plane."""
        if self.reference_policy is ReferencePolicy.FULL_IS_REFERENCE:
            return "full"
        if self.reference_policy is ReferencePolicy.INDEPENDENT_PROJECTION:
            return "reference"
        return None

    @property
    def fitting_peak_list_key(self) -> str | None:
        """Authoritative peak-list key used by fitting/grouping operations."""
        return self.projection_peak_list_key


def spectrum_policy(mode: AnalysisMode) -> SpectrumPolicy:
    """Return peak-list ownership rules for a canonical journey.

    Projection peak lists used to draw markers on projected/protected views are
    derived display data and are deliberately not represented as an additional
    source of peak authority here.
    """
    n = mode.spectral_dimensions
    if mode.has_pseudo_axis:
        if n == 1:
            return SpectrumPolicy("pseudo2D", 1, ReferencePolicy.FULL_IS_REFERENCE, 1, True)
        if n == 2:
            return SpectrumPolicy("pseudo3D", 2, ReferencePolicy.INDEPENDENT_PROJECTION, 2, True)
        # Not currently a canonical journey, but keep the policy general for a
        # future pseudo-dimensional recovery.
        return SpectrumPolicy(f"pseudo{n + 1}D", n, ReferencePolicy.INDEPENDENT_PROJECTION, max(1, n - 1), True)

    if n == 1:
        return SpectrumPolicy("1D", 1, ReferencePolicy.NONE, None, False)
    if n == 2:
        # Physical 2D is the singular case: the projection/reference/full peak
        # concepts identify the same physical 2D peak collection.  There is no
        # independent lower-dimensional reference list to curate.
        return SpectrumPolicy("2D", 2, ReferencePolicy.FULL_IS_REFERENCE, 2, False)
    return SpectrumPolicy(f"{n}D", n, ReferencePolicy.INDEPENDENT_PROJECTION, n - 1, False)
