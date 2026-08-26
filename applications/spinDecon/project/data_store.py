"""Shared data model for spectra, projections, and peak lists."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Hashable
import traceback


@dataclass
class DataStore:
    """Single in-memory store for the GUI.

    Peak-list architecture: every current peak collection, irrespective of
    dimensionality, uses ``peak_lists`` and the same record schema.  Multiple
    semantic roles may coexist (for example reference/projection, full, and
    deconvolved lists).  Legacy connectivity/NOE ``conn_data`` is not part of
    the current data model.
    """

    spectra: dict[str, dict[str, Any]] = field(default_factory=dict)
    projections: dict[Hashable, dict[str, Any]] = field(default_factory=dict)
    views: dict[Hashable, dict[str, Any]] = field(default_factory=dict)
    projected_peak_lists: dict[Hashable, dict[str, Any]] = field(default_factory=dict)
    peak_lists: dict[str, dict[str, Any]] = field(default_factory=lambda: {
        "reference": {"peaks": [], "dimension": 2, "source_path": ""},
        "full": {"peaks": [], "rows": [], "dimension": 0, "source_path": ""},
    })
    metadata: dict[str, Any] = field(default_factory=dict)
    # Analysis-run metadata. Per-peak measurements live on the peak records.
    analysis: dict[str, Any] = field(default_factory=dict)

    # Compatibility fields used by the legacy widgets.
    dic: Any = None
    data: Any = None
    uc0: Any = None
    uc1: Any = None
    uc2: Any = None
    uc3: Any = None
    index0: Any = None
    index1: Any = None
    index2: Any = None
    index3: Any = None
    XX: Any = None
    YY: Any = None
    ZZ: Any = None
    XX2: Any = None
    YY2: Any = None
    uc0min: Any = None
    uc0max: Any = None
    uc1min: Any = None
    uc1max: Any = None
    uc2min: Any = None
    uc2max: Any = None
    uc3min: Any = None
    uc3max: Any = None
    dmax: Any = None
    noiseVal: Any = None
    projectedData: Any = None
    datadec: Any = None
    dicdec: Any = None
    peak: list[Any] = field(default_factory=list)
    spectrumfile: str = ""
    STD_raw_path: str = ""
    STD_std_path: str = ""
    STD: Any = None
    mixingTimes: Any = None
    labb: Any = None
    dim: int = 0
    pkIdx: list[Any] = field(default_factory=list)
    pkSlice1D: list[Any] = field(default_factory=list)
    pkSlice1Ddec: list[Any] = field(default_factory=list)
    Grps: dict[str, Any] = field(default_factory=dict)
    noeTags: list[Any] = field(default_factory=list)
    READ: int = 0
    DECON: int = 0
    PEAK: int = 0
    pseudo: bool = False
    # Retired uSTA-tab state retained solely for loading historical projects.
    uSTA: bool = False

    def _snapshot(self, **kwargs: Any) -> dict[str, Any]:
        return {k: v for k, v in kwargs.items() if v is not None}

    def save_spectrum(self, key: str = "raw", **kwargs: Any) -> None:
        payload = self._snapshot(**kwargs)
        self.spectra[key] = payload
        if key == "raw":
            # Plotting views are derived from the main spectrum.  Never allow a
            # Pseudo3D/Fitting pane to retain a view from a previously loaded
            # spectrum when the authoritative raw spectrum changes.
            for view_key in list(self.views):
                if isinstance(view_key, tuple) and view_key and view_key[0] == "pseudo3d":
                    self.views.pop(view_key, None)
            for name, value in payload.items():
                if hasattr(self, name):
                    setattr(self, name, value)
        elif key == "decon":
            # A Pseudo3D pane may have been opened before deconvolution.  Drop
            # any cached calculated view so its next lookup is built from this
            # newly published deconvolved spectrum rather than stale state.
            self.views.pop(("pseudo3d", "decon"), None)
            if "dic" in payload:
                self.dicdec = payload["dic"]
            if "data" in payload:
                self.datadec = payload["data"]

    def save_projection(self, key: Hashable, **kwargs: Any) -> None:
        self.projections[key] = self._snapshot(**kwargs)

    def save_view(self, key: Hashable, **kwargs: Any) -> None:
        """Store a plotting-ready view owned by the notebook data model."""
        self.views[key] = self._snapshot(**kwargs)

    def get_view(self, key: Hashable):
        return self.views.get(key)

    def save_peak_list(self, key: str = "full", **kwargs: Any) -> None:
        """Store a named peak-list payload without mutating legacy fields.

        Peak lists have distinct semantic roles.  In particular, writing the
        full nD list must never overwrite the 2D reference list used by the
        slice viewers and bore-mode deconvolution.
        """
        payload = self._snapshot(**kwargs)
        current = dict(self.peak_lists.get(key, {}))
        current.update(payload)
        self.peak_lists[key] = current

    def get_peak_list(self, key: str) -> dict[str, Any]:
        return self.peak_lists.setdefault(key, {})

    def reset(self) -> None:
        """Clear all cached spectra, projections, and peak lists in place."""
        pass
        pass
        fresh = DataStore()
        self.__dict__.clear()
        self.__dict__.update(fresh.__dict__)
        pass

    def set_metadata(self, **kwargs: Any) -> None:
        self.metadata.update(kwargs)

    # Workflow evidence API -------------------------------------------------
    #
    # These methods are the single write boundary for evidence consumed by
    # workflow_status.  Scientific widgets should call them only after an
    # operation has completed successfully; opening a tool is not evidence.

    def mark_peak_shape_determined(self, **details: Any) -> None:
        self.metadata["peak_shape_determined"] = True
        if details:
            self.metadata["peak_shape"] = self._snapshot(**details)

    def invalidate_pseudo_series_review(self) -> None:
        """Invalidate acceptance of the current pseudo-dimensional intensity series."""
        pass
        self.analysis.pop("pseudo_series_reviewed", None)
        self.analysis.pop("pseudo_series_review", None)
        pass

    def mark_pseudo_intensities_ready(self, *, invalidate_review: bool = True, **details: Any) -> None:
        # A newly generated extraction invalidates acceptance of the previous
        # intensity series.  Merely rediscovering already-existing fit files
        # during project restore must not do so, otherwise a persisted
        # pseudoSeriesInspected=1 is immediately lost on the next session.
        if invalidate_review:
            self.invalidate_pseudo_series_review()
        self.analysis["pseudo_intensities_ready"] = True
        if details:
            self.analysis["pseudo_intensities"] = self._snapshot(**details)

    def mark_pseudo_series_reviewed(self, **details: Any) -> None:
        pass
        self.analysis["pseudo_series_reviewed"] = True
        if details:
            self.analysis["pseudo_series_review"] = self._snapshot(**details)
        pass

    def invalidate_picked_peaks_review(self) -> None:
        """Invalidate explicit acceptance of the current full peak list."""
        self.analysis.pop("picked_peaks_reviewed", None)
        self.analysis.pop("picked_peaks_review", None)

    def mark_picked_peaks_reviewed(self, **details: Any) -> None:
        """Record explicit user acceptance of the current full peak list."""
        self.analysis["picked_peaks_reviewed"] = True
        if details:
            self.analysis["picked_peaks_review"] = self._snapshot(**details)

    def invalidate_fitting_review(self) -> None:
        """Invalidate acceptance of the current physical-2D fitting results."""
        self.analysis.pop("fitting_results_ready", None)
        self.analysis.pop("fitting_results", None)
        self.analysis.pop("fitting_results_reviewed", None)
        self.analysis.pop("fitting_results_review", None)

    def mark_fitting_results_ready(self, **details: Any) -> None:
        self.analysis["fitting_results_ready"] = True
        self.analysis.pop("fitting_results_reviewed", None)
        self.analysis.pop("fitting_results_review", None)
        if details:
            self.analysis["fitting_results"] = self._snapshot(**details)

    def invalidate_fitting_results_review(self) -> None:
        self.analysis.pop("fitting_results_reviewed", None)
        self.analysis.pop("fitting_results_review", None)

    def mark_fitting_results_reviewed(self, **details: Any) -> None:
        self.analysis["fitting_results_reviewed"] = True
        if details:
            self.analysis["fitting_results_review"] = self._snapshot(**details)

    def mark_pseudo_analysis_complete(self, **details: Any) -> None:
        self.analysis["pseudo_analysis_complete"] = True
        if details:
            self.analysis["pseudo_analysis"] = self._snapshot(**details)
