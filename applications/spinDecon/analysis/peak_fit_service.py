"""Scientific-data boundary for the peak-shape fitting workspace."""
import numpy as np


class PeakFitService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def topology(self):
        return self._legacy._active_topology()

    @property
    def peaks(self):
        return self._legacy.peak

    @property
    def parameter_file(self):
        return getattr(self._legacy, "deconParFile", "")

    def threshold_fraction(self):
        getter = getattr(self._legacy, "get_threshold_fraction", None)
        return float(getter()) if callable(getter) else float(self._legacy.threshBox.GetValue())

    def fitting_payload(self):
        topology = self.topology
        dim = int(topology.spectral_dim_count)
        if topology.has_pseudo_axis:
            if dim == 2:
                view = self._legacy.get_pseudo3d_view("raw")
                if view is None:
                    raise ValueError("Cannot resolve spectral axes for pseudo-dimensional peak fitting.")
                data = np.asarray(view["data"][0])
                labels = [view["y_label"], view["x_label"]]
            elif dim == 1:
                projection = self._legacy.get_pseudo2d_projection_data(ensure_file=True)
                if projection is None:
                    raise ValueError("Cannot create the 1D spectral projection for pseudo2D peak fitting.")
                data = np.asarray(projection["data"]).squeeze()
                labels = [projection["label"]]
            else:
                axes = self._legacy._spectral_physical_axes()
                data = np.take(self._legacy.data, 0, axis=topology.pseudo_axis.physical_index)
                labels = [label for _, label in axes]
            indexes = [getattr(self._legacy, "index%d" % physical_index)
                       for physical_index, _ in self._legacy._spectral_physical_axes()]
        else:
            data = self._legacy.data
            labels = self._legacy.labb
            indexes = [getattr(self._legacy, "index%d" % i) for i in range(dim)]
        return {"topology": topology, "dimension": dim, "data": data,
                "labels": labels, "indexes": indexes, "peaks": self.peaks}

    def shape_parameters(self, dimension=None):
        """Return numeric peak-shape parameters without exposing wx controls."""
        dim = int(dimension if dimension is not None else self.topology.spectral_dim_count)
        sigmas = []
        voigt = []
        lorentz = []
        for which in range(1, dim + 1):
            sigmas.append(float(getattr(self._legacy, 'sig%dBox' % which).GetValue()))
            voigt.append(float(getattr(self._legacy, 'voigt%dBox' % which).GetValue()))
            lorentz.append(float(getattr(self._legacy, 'lorentz%dBox' % which).GetValue()))
        return {'sigmas': tuple(sigmas), 'voigt': tuple(voigt), 'lorentz': tuple(lorentz)}

    def fit_radius(self, which, dimension=None, default=0.1):
        """Resolve a persisted fitting radius behind the application boundary."""
        dim = int(dimension if dimension is not None else self.topology.spectral_dim_count)
        box = getattr(self._legacy, 'fitF%dBox' % int(which), None)
        if int(which) == 1 and dim == 1:
            legacy = getattr(self._legacy, 'fitRadBox', None)
            try:
                text = legacy.GetValue().strip() if legacy is not None else ''
                if text:
                    return abs(float(text))
            except Exception:
                pass
        if box is not None:
            try:
                text = box.GetValue().strip()
                if text:
                    return abs(float(text))
            except Exception:
                pass
        getter = getattr(self._legacy, 'get_parameter_float', None)
        if callable(getter):
            try:
                value = getter('3p_radF%d' % int(which))
                if value not in (None, 0, 0.0):
                    return abs(float(value))
            except Exception:
                pass
        return float(default)

    @property
    def data(self):
        return self._legacy.data


    @property
    def visible_axes(self):
        """Return the parent spectrum axes used only for optional view-region limiting."""
        return getattr(self._legacy, 'axes', None)

    def set_shape_parameters(self, sigmas, voigt, lorentz):
        for which, values in enumerate(zip(sigmas, voigt, lorentz), start=1):
            sigma, vfrac, lfrac = values
            getattr(self._legacy, 'sig%dBox' % which).SetValue(str('%.3g' % float(sigma)))
            getattr(self._legacy, 'voigt%dBox' % which).SetValue(str('%.3g' % float(vfrac)))
            getattr(self._legacy, 'lorentz%dBox' % which).SetValue(str('%.3g' % float(lfrac)))

    def set_fit_radius(self, which, value, dimension=None):
        dim = int(dimension if dimension is not None else self.topology.spectral_dim_count)
        box = getattr(self._legacy, 'fitF%dBox' % int(which), None)
        if box is not None:
            box.SetValue('%.6g' % float(value))
        if int(which) == 1 and dim == 1:
            legacy = getattr(self._legacy, 'fitRadBox', None)
            if legacy is not None:
                legacy.SetValue('%.6g' % float(value))

    def update_pseudo3d_parameters(self, values):
        """Persist pseudo-dimensional parameters using the workspace mapping contract."""
        updater = getattr(self._legacy, 'update_pseudo3d_parameters', None)
        return updater(dict(values)) if callable(updater) else None

    def save_project(self):
        save = getattr(self._legacy, 'OnButtonSave', None)
        return save(None) if callable(save) else None

    def save_fit_preferences(self, peak_count, link_widths):
        self._legacy.peak_fit_count = int(peak_count)
        self._legacy.peak_fit_link_widths = bool(link_widths)

    def mark_peak_shape_determined(self, dimension, sigmas, voigt, lorentz):
        store = getattr(self._legacy, 'store', None) or getattr(self._legacy, 'data_store', None)
        was_already_fitted = bool(getattr(self._legacy, 'peak_shape_fitted', False) or
                                  (store is not None and getattr(store, 'metadata', {}).get('peak_shape_determined')))
        if store is not None and hasattr(store, 'mark_peak_shape_determined'):
            store.mark_peak_shape_determined(source='peak_fit', dimension=int(dimension),
                                             sigma=list(sigmas), voigt=list(voigt), lorentz=list(lorentz))
        saved = getattr(self._legacy, 'peak_shape_saved', None)
        if callable(saved):
            saved(was_already_fitted=was_already_fitted)
        else:
            notebook = getattr(self._legacy, 'parent', None)
            notify = getattr(notebook, 'notify_analysis_changed', None)
            if callable(notify): notify()

    def sync_usta_shape(self, sigma, voigt, lorentz):
        if not bool(getattr(self._legacy, 'uSTA', False)):
            return
        tab = getattr(getattr(self._legacy, 'parent', None), 'tabSTD', None)
        if tab is None:
            return
        tab.sig1Box.SetValue(str('%.3g' % float(sigma)))
        tab.voigt1Box.SetValue(str('%.3g' % float(voigt)))
        tab.lorentz1Box.SetValue(str('%.3g' % float(lorentz)))
