"""Application boundary for pseudo-axis workspaces.

This migration adapter centralises access to pseudo-axis data and project
persistence while the legacy wx workspaces are progressively decoupled from
``deconFrame``.
"""
import os
import numpy as np


class PseudoAxisService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def data(self):
        return np.asarray(self._legacy.data)

    @property
    def labels(self):
        return tuple(getattr(self._legacy, "labb", ()) or ())

    def spectral_axis(self, shape=None):
        shape = tuple(shape or self.data.shape)
        for name in ("index1", "index0"):
            arr = getattr(self._legacy, name, None)
            if arr is not None and len(arr) in shape:
                return np.asarray(arr, dtype=float)
        return np.arange(shape[-1], dtype=float)

    @property
    def downstream_analysis(self):
        return str(getattr(self._legacy, "downstream_analysis", "") or "")

    def set_downstream_analysis(self, selection):
        self._legacy.downstream_analysis = selection
        store = getattr(self._legacy, "store", None) or getattr(self._legacy, "data_store", None)
        if store is not None:
            try:
                store.metadata["downstream_analysis"] = selection
            except Exception:
                pass
        save = getattr(self._legacy, "OnButtonSave", None)
        if callable(save):
            try:
                save(None)
            except TypeError:
                save()
        self.notify_changed()

    def full_peak_payload(self):
        getter = getattr(self._legacy, "get_full_peak_payload", None)
        return getter() if callable(getter) else {}

    def notify_changed(self):
        notify = getattr(self._legacy, "_notify_analysis_changed", None)
        if callable(notify):
            notify()

    def view(self, kind="raw"):
        return self._legacy.get_pseudo3d_view(kind)

    def projection_view(self, *args, **kwargs):
        return self._legacy.get_projection_view(*args, **kwargs)

    def threshold_fraction(self):
        getter = getattr(self._legacy, "get_threshold_fraction", None)
        return float(getter()) if callable(getter) else float(self._legacy.threshBox.GetValue())

    def spectrum_path(self):
        getter = getattr(self._legacy, "get_spectrum_path", None)
        return getter() if callable(getter) else str(getattr(self._legacy, "spectrumfile", "") or "")

    def output_dir(self):
        getter = getattr(self._legacy, "_spec_output_dir", None)
        return getter() if callable(getter) else ""

    def fuda_dir(self): return self._legacy.get_fuda_dir()
    def fuda_peak_file(self): return self._legacy.get_fuda_peak_file()
    def fuda_parameter_file(self): return self._legacy.get_fuda_parameter_file()
    def parameter_value(self, name, default=""): return self._legacy.get_parameter_value(name, default=default)
    def parameter_float(self, name, default=0.0): return self._legacy.get_parameter_float(name, default=default)
    def update_parameters(self, values): return self._legacy.update_pseudo3d_parameters(values)
    def rebuild_projected_peaks(self): return self._legacy._rebuild_projected_peak_lists()
    def replace_groups(self, groups): return self._legacy.replace_pseudo3d_groups(groups)
    def add_group(self, name, peaks): return self._legacy.add_pseudo3d_group(name, peaks)
    def remove_group(self, name): return self._legacy.remove_pseudo3d_group(name)
    def add_peak_to_group(self, name, peak_name): return self._legacy.add_peak_to_pseudo3d_group(name, peak_name)
    def remove_peak_from_group(self, name, peak_name): return self._legacy.remove_peak_from_pseudo3d_group(name, peak_name)

    def fit_dir(self):
        """Return the canonical restrained-fit directory for pseudo-axis analyses."""
        path = self.fuda_dir()
        if path:
            return path
        root = self.output_dir() or os.path.dirname(self.spectrum_path())
        return os.path.join(root, "fit")

    def full_peak_file(self):
        """Resolve the authoritative Full Peak List file used by pseudo2D fitting."""
        box = getattr(self._legacy, "fullPeakBox", None)
        value = box.GetValue() if box is not None and hasattr(box, "GetValue") else ""
        resolver = getattr(self._legacy, "_resolve_spec_file", None)
        return resolver(value) if callable(resolver) else value

    def mark_series_reviewed(self, source="pseudo2d_fitting"):
        """Persist explicit review evidence and refresh workflow status."""
        store = getattr(self._legacy, "store", None) or getattr(self._legacy, "data_store", None)
        if store is None:
            return False
        store.mark_pseudo_series_reviewed(source=source)
        save = getattr(self._legacy, "OnButtonSave", None)
        if callable(save):
            try:
                save(True)
            except TypeError:
                save()
        self.notify_changed()
        return True

    def series_reviewed(self):
        store = getattr(self._legacy, "store", None) or getattr(self._legacy, "data_store", None)
        return bool(getattr(store, "analysis", {}).get("pseudo_series_reviewed")) if store is not None else False

    @property
    def unit_conversion_bounds(self):
        return (getattr(self._legacy, "uc0min", None), getattr(self._legacy, "uc0max", None))
