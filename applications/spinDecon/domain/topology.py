"""Canonical dimensionality vocabulary for NMR datasets.

This module is intentionally independent of wxPython, ProjectState, and the
legacy ``dim``/``pseudo`` representation.  It defines the dimensionality
contract that GUI and processing code can migrate to incrementally.

Terminology
-----------
``spectral_dim_count``
    Number of frequency / chemical-shift dimensions (1..4).
``has_pseudo_axis``
    Whether the dataset has one additional sampled real (non-spectral) axis.
``physical_dim_count``
    Number of axes in the physical data array.  This is derived as
    ``spectral_dim_count + int(has_pseudo_axis)``.

A bare integer axis number is deliberately not enough to describe an axis:
``AxisSpec`` records both its physical array index and, for spectral axes, its
spectral index.  This avoids conflating e.g. physical axis 1 with spectral
axis 1 when a pseudo axis is present.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Optional, Tuple


class AxisKind(str, Enum):
    """Scientific role of a physical data axis."""

    SPECTRAL = "spectral"
    PSEUDO_REAL = "pseudo_real"


@dataclass(frozen=True)
class AxisSpec:
    """Identity and scientific role of one physical data-array axis.

    ``physical_index`` is the index in the physical array.  ``spectral_index``
    is the zero-based index amongst spectral dimensions only and must be
    ``None`` for the real pseudo axis.
    """

    physical_index: int
    kind: AxisKind
    label: str = ""
    spectral_index: Optional[int] = None

    def __post_init__(self) -> None:
        if self.physical_index < 0:
            raise ValueError("physical_index must be >= 0")

        kind = AxisKind(self.kind)
        object.__setattr__(self, "kind", kind)

        if kind is AxisKind.SPECTRAL:
            if self.spectral_index is None or self.spectral_index < 0:
                raise ValueError("spectral axes require spectral_index >= 0")
        elif self.spectral_index is not None:
            raise ValueError("the pseudo-real axis cannot have a spectral_index")

    @property
    def is_spectral(self) -> bool:
        return self.kind is AxisKind.SPECTRAL

    @property
    def is_pseudo_real(self) -> bool:
        return self.kind is AxisKind.PSEUDO_REAL


@dataclass(frozen=True)
class DatasetTopology:
    """Authoritative description of spectral and physical dimensionality.

    The object validates its own invariants.  It does not inspect a NumPy
    array; callers can use :meth:`validate_data_ndim` at data-loading
    boundaries when they have an array available.
    """

    spectral_dim_count: int
    has_pseudo_axis: bool
    axes: Tuple[AxisSpec, ...]

    def __post_init__(self) -> None:
        spectral_count = int(self.spectral_dim_count)
        if spectral_count < 1 or spectral_count > 4:
            raise ValueError("spectral_dim_count must be an integer from 1 to 4")
        object.__setattr__(self, "spectral_dim_count", spectral_count)
        object.__setattr__(self, "has_pseudo_axis", bool(self.has_pseudo_axis))
        object.__setattr__(self, "axes", tuple(self.axes))

        if len(self.axes) != self.physical_dim_count:
            raise ValueError(
                "number of AxisSpec entries must equal physical_dim_count "
                f"({self.physical_dim_count})"
            )

        physical_indices = [axis.physical_index for axis in self.axes]
        expected_physical = list(range(self.physical_dim_count))
        if sorted(physical_indices) != expected_physical:
            raise ValueError(
                "physical axis indices must be unique and contiguous from 0 to "
                f"{self.physical_dim_count - 1}"
            )

        spectral_axes = self.spectral_axes
        if len(spectral_axes) != self.spectral_dim_count:
            raise ValueError("AxisSpec spectral-axis count disagrees with spectral_dim_count")

        spectral_indices = [axis.spectral_index for axis in spectral_axes]
        if sorted(spectral_indices) != list(range(self.spectral_dim_count)):
            raise ValueError(
                "spectral_index values must be unique and contiguous from 0 to "
                f"{self.spectral_dim_count - 1}"
            )

        if len(self.pseudo_axes) != self.pseudo_dim_count:
            raise ValueError("AxisSpec pseudo-axis count disagrees with has_pseudo_axis")

    @property
    def pseudo_dim_count(self) -> int:
        return int(self.has_pseudo_axis)

    @property
    def physical_dim_count(self) -> int:
        return self.spectral_dim_count + self.pseudo_dim_count

    @property
    def physical_axes(self) -> Tuple[AxisSpec, ...]:
        """Axes ordered by their physical array index."""
        return tuple(sorted(self.axes, key=lambda axis: axis.physical_index))

    @property
    def spectral_axes(self) -> Tuple[AxisSpec, ...]:
        """Spectral axes ordered by spectral index, independent of array order."""
        return tuple(
            sorted(
                (axis for axis in self.axes if axis.is_spectral),
                key=lambda axis: axis.spectral_index,
            )
        )

    @property
    def pseudo_axes(self) -> Tuple[AxisSpec, ...]:
        return tuple(axis for axis in self.physical_axes if axis.is_pseudo_real)

    @property
    def pseudo_axis(self) -> Optional[AxisSpec]:
        return self.pseudo_axes[0] if self.pseudo_axes else None

    @property
    def pseudo_physical_index(self) -> Optional[int]:
        """Physical array index of the pseudo-real axis, if present.

        Kept as a convenience/compatibility accessor while callers migrate to
        ``pseudo_axis.physical_index``.
        """
        axis = self.pseudo_axis
        return axis.physical_index if axis is not None else None

    def validate_data_ndim(self, ndim: int) -> None:
        """Raise ValueError unless a data array has the expected physical ndim."""
        actual = int(ndim)
        if actual != self.physical_dim_count:
            raise ValueError(
                f"data.ndim={actual} does not match topology physical_dim_count="
                f"{self.physical_dim_count}"
            )

    @classmethod
    def from_counts(
        cls,
        spectral_dim_count: int,
        has_pseudo_axis: bool = False,
        *,
        pseudo_physical_index: int = 0,
        spectral_labels: Optional[Iterable[str]] = None,
        pseudo_label: str = "",
    ) -> "DatasetTopology":
        """Construct a topology when only counts/order are known.

        This convenience constructor is primarily useful at UI/test boundaries.
        When a loader knows the real physical axis order it should construct
        explicit ``AxisSpec`` entries instead.  If a pseudo axis is requested,
        its physical position must therefore be supplied or explicitly accepted
        as the default (physical axis 0).
        """
        spectral_count = int(spectral_dim_count)
        pseudo = bool(has_pseudo_axis)
        physical_count = spectral_count + int(pseudo)

        if spectral_count < 1 or spectral_count > 4:
            raise ValueError("spectral_dim_count must be an integer from 1 to 4")
        if pseudo and not 0 <= int(pseudo_physical_index) < physical_count:
            raise ValueError("pseudo_physical_index is outside the physical axis range")

        labels = tuple(spectral_labels or ())
        if labels and len(labels) != spectral_count:
            raise ValueError("spectral_labels must contain one label per spectral axis")
        if not labels:
            labels = tuple("" for _ in range(spectral_count))

        axes = []
        spectral_index = 0
        for physical_index in range(physical_count):
            if pseudo and physical_index == int(pseudo_physical_index):
                axes.append(
                    AxisSpec(
                        physical_index=physical_index,
                        kind=AxisKind.PSEUDO_REAL,
                        label=pseudo_label,
                    )
                )
            else:
                axes.append(
                    AxisSpec(
                        physical_index=physical_index,
                        kind=AxisKind.SPECTRAL,
                        label=labels[spectral_index],
                        spectral_index=spectral_index,
                    )
                )
                spectral_index += 1

        return cls(
            spectral_dim_count=spectral_count,
            has_pseudo_axis=pseudo,
            axes=tuple(axes),
        )
