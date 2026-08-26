from .full_peak_service import FullPeakListService

"""Application boundary used by projection workspaces during GUI migration.

This adapter intentionally delegates to the legacy NMR workspace.  Projection
GUI code can therefore stop depending on the concrete ``deconFrame`` API while
scientific implementations are extracted independently.
"""

class ProjectionService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def data(self): return self._legacy.data

    @property
    def labels(self): return tuple(getattr(self._legacy, 'labb', ()) or ())

    @property
    def peaks(self): return self._legacy.peak

    def threshold_fraction(self): return float(self._legacy.threshBox.GetValue())
    def pseudo2d_data(self, ensure_file=False): return self._legacy.get_pseudo2d_projection_data(ensure_file=ensure_file)
    def full_peak_payload(self): return FullPeakListService(self._legacy).payload()
    def spectrum_view(self, *args, **kwargs): return self._legacy.get_spectrum_view(*args, **kwargs)
    def projection_view(self, *args, **kwargs): return self._legacy.get_projection_view(*args, **kwargs)
    def projected_peak_overlay(self, *args, **kwargs): return self._legacy.get_projected_peak_overlay(*args, **kwargs)
    def reference_peaks(self): return self._legacy.get_reference_peaks()
    def clear_peak_selection(self, *args, **kwargs): return self._legacy.clear_peak_selection(*args, **kwargs)
    def alias(self, *args, **kwargs): return self._legacy.alias(*args, **kwargs)

    def save_full_peak_list(self, *, peaks, records, rows, dimension, source_path=None):
        self._legacy.store.save_peak_list('full', peaks=peaks, records=records,
                                          rows=rows, dimension=dimension,
                                          source_path=source_path)
        rebuild = getattr(self._legacy, '_rebuild_projected_peak_lists', None)
        if callable(rebuild): rebuild()
        notify = getattr(self._legacy, '_notify_analysis_changed', None)
        if callable(notify): notify()
        refresh = getattr(self._legacy, 'refresh_full_peak_list_viewers', None)
        if callable(refresh): refresh()
        status = getattr(self._legacy, 'Status', None)
        if callable(status): status()

    def focus_full_peak_list_viewers(self, name=None):
        focus = getattr(self._legacy, 'focus_full_peak_list_viewers', None)
        if not callable(focus): return None
        return focus(name) if name is not None else focus()

    def pseudo2d_decon_projection(self):
        store = getattr(self._legacy, 'store', None)
        payload = store.spectra.get('pseudo2d_projection_decon', {}) if store is not None else {}
        data = payload.get('data')
        if data is None:
            cached = getattr(self._legacy, 'pseudo2d_projection_decon', None)
            data = None if cached is None else cached.get('data')
        return data

    def intensity_threshold(self):
        dmax = getattr(self._legacy, 'dmax', None)
        if dmax is None:
            import numpy as np
            loaded = getattr(self._legacy, 'data', None)
            if loaded is None:
                raise RuntimeError('Projection opened before spectrum data were loaded')
            dmax = float(np.max(np.fabs(np.asarray(loaded))))
            self._legacy.dmax = dmax
        return float(dmax) * self.threshold_fraction()

    def peak_shape_parameters(self, dimension=2):
        dim = int(dimension)
        return {
            'sigma': tuple(float(getattr(self._legacy, 'sig%dBox' % i).GetValue()) for i in range(1, dim + 1)),
            'lorentz': tuple(float(getattr(self._legacy, 'lorentz%dBox' % i).GetValue()) for i in range(1, dim + 1)),
            'voigt': tuple(float(getattr(self._legacy, 'voigt%dBox' % i).GetValue()) for i in range(1, dim + 1)),
        }
