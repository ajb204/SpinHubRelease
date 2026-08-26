"""Application boundary for peak editing and peak-analysis workspaces.

The adapter intentionally delegates to the legacy NMR workspace while peak
scientific operations are extracted.  GUI code should depend on this boundary
rather than on concrete ``deconFrame`` widgets.
"""
from pathlib import Path


class PeakService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def dimension(self):
        return int(getattr(self._legacy, "dim", 0) or 0)

    @property
    def labels(self):
        return tuple(getattr(self._legacy, "labb", ()) or ())

    @property
    def data(self):
        return self._legacy.data

    @property
    def peaks(self):
        return self._legacy.peak

    def threshold(self):
        return float(self._legacy.dmax) * float(self._legacy.threshBox.GetValue())

    def spectrum_path(self):
        path = getattr(self._legacy, "spectrumfile", "")
        if path:
            return str(path)
        try:
            return self._legacy._resolve_input_path(self._legacy.infileBox.GetValue())
        except Exception:
            return ""

    def is_pseudo3d(self):
        checker = getattr(self._legacy, "_is_pseudo3d_topology", None)
        return bool(checker()) if callable(checker) else False

    def pseudo3d_view(self, kind="raw"):
        return self._legacy.get_pseudo3d_view(kind)

    def spectrum_view(self, *args, **kwargs):
        return self._legacy.get_spectrum_view(*args, **kwargs)

    def projection_view(self, *args, **kwargs):
        return self._legacy.get_projection_view(*args, **kwargs)
    def commit_projection_peaks(self, peaks, *, full_2d=False, source_path=None,
                                projection_key=None, labels=None):
        """Persist edited projection peaks through the canonical controller APIs."""
        peaks = list(peaks or [])
        store = getattr(self._legacy, "store", None) or getattr(self._legacy, "data_store", None)
        state = getattr(self._legacy, "state", None)
        if store is None:
            return
        if full_2d:
            rows = [[str(pk.name), ("%.10g" % float(pk.y)), ("%.10g" % float(pk.x))]
                    for pk in peaks]
            records = self._legacy._full_peak_records(rows, dim=2)
            store.save_peak_list("full", peaks=records, records=records, rows=rows,
                                 dimension=2, source_path=source_path)
            self._legacy._rebuild_projected_peak_lists()
            if state is not None and source_path:
                rel_path = state._spec_relative(source_path)
                self._legacy.fullPeakBox.SetValue(rel_path)
                state.full_peak_file = rel_path
                state.dirty = True
            refresh = getattr(self._legacy, "refresh_full_peak_list_viewers", None)
            if callable(refresh):
                refresh()
            # Physical 2D has one authoritative peak collection.  Do not
            # manufacture a shadow Reference list: projection, fitting and
            # PeakFrame all consume the Full 2D list directly.
            self.notify_changed()
        else:
            setter = getattr(self._legacy, "set_reference_peaks", None)
            if callable(setter):
                setter(peaks)
                refresh = getattr(self._legacy, "refresh_reference_peak_views", None)
                if callable(refresh):
                    refresh()
            else:
                store.save_peak_list("reference", peaks=peaks, dimension=2,
                                     projection_key=projection_key, labels=labels)
        return peaks

    def notify_changed(self):
        notify = getattr(self._legacy, "_notify_analysis_changed", None)
        if callable(notify):
            notify()


    def view_labels(self, spectral_dim_count=None):
        """Return the canonical X/Y labels for the peak-editing plane."""
        dim = self.dimension if spectral_dim_count is None else int(spectral_dim_count)
        labels = list(self.labels)
        if not labels:
            return ('x', 'y')
        if self.is_pseudo3d():
            axes = self._legacy._spectral_physical_axes()
            if len(axes) == 2:
                return (axes[1][1], axes[0][1])
        if dim == 3 and len(labels) >= 3:
            pseudo_labels = {'time_T2', 'ID', 'ncyc', 'ncyc_cp', 'gzlvl5', 'gzlvl1'}
            spectral = [str(label) for label in labels[:3] if str(label) not in pseudo_labels]
            if len(spectral) == 2:
                return (spectral[-1], spectral[-2])
            return (labels[2], labels[1])
        if dim == 4 and len(labels) >= 4:
            return (labels[2], labels[3])
        if len(labels) >= 2:
            return (labels[1], labels[0])
        return (labels[0], labels[0])

    def bore_payload(self, spectral_dim_count=None):
        """Return (axis, cube, label) for the optional third-axis trace."""
        import numpy
        dim = self.dimension if spectral_dim_count is None else int(spectral_dim_count)
        if self.is_pseudo3d():
            view = self.pseudo3d_view('raw')
            if view is None:
                return None
            return (numpy.asarray(view['pseudo_axis']), numpy.asarray(view['data']),
                    str(view.get('pseudo_label') or 'pseudoaxis'))
        if dim == 3:
            labels = self.labels
            return (numpy.asarray(self._legacy.uc0.ppms_scale), numpy.asarray(self.data),
                    str(labels[0] if labels else 'z'))
        return None

    def cached_view(self, key):
        store = getattr(self._legacy, 'store', None) or getattr(self._legacy, 'data_store', None)
        return store.get_view(key) if store is not None else None

    def save_view(self, key, **payload):
        store = getattr(self._legacy, 'store', None) or getattr(self._legacy, 'data_store', None)
        if store is not None:
            store.save_view(key, **payload)

    def spectrum_payload_from_file(self, path, labels=None):
        """Load a spectrum file through the legacy payload builder during migration."""
        import nmrglue as ng
        dic, data = ng.pipe.read(path)
        return self._legacy._spectrum_view_payload(
            dic, data, source=path, labb=labels, transpose='n')

    def unit_converter(self, axis):
        return getattr(self._legacy, 'uc%d' % int(axis), None)

    def projection_payload(self, labels, *, decon=False, analysis_path=None):
        """Resolve the plotting-ready 2D plane for PeakFrame.

        This centralises DataStore/view selection while the legacy workspace still
        owns file decoding and projection generation.
        """
        import os
        labels = tuple(labels)
        store = getattr(self._legacy, 'store', None) or getattr(self._legacy, 'data_store', None)
        if store is None:
            raise RuntimeError('PeakService requires the application data store')
        if self.is_pseudo3d():
            if decon:
                key = ('peakframe_decon', labels[0], labels[1], 'n')
                view = store.get_view(key)
                decon_path = (analysis_path or '') + '.decon'
                if view is None and analysis_path and os.path.exists(decon_path):
                    try:
                        view = self.spectrum_payload_from_file(decon_path, labels=labels)
                        store.save_view(key, **view)
                    except Exception:
                        view = None
                if view is None:
                    raise RuntimeError('PeakFrame cannot resolve the pseudo3D projection deconvolution')
            else:
                pseudo_view = self.pseudo3d_view('raw')
                if pseudo_view is None:
                    raise RuntimeError('PeakFrame cannot resolve the pseudo-dimensional spectral plane')
                view = dict(pseudo_view)
                view['ZZ'] = pseudo_view['data'][0]
        elif self.dimension <= 2:
            view = self.spectrum_view(decon=decon, transpose='n')
        else:
            key = ('peakframe_decon', labels[0], labels[1], 'n')
            if decon:
                view = store.get_view(key)
                if view is None:
                    view = self.projection_view(labels[0], labels[1], decon=True, transpose='n')
            else:
                view = self.projection_view(labels[0], labels[1], decon=False, transpose='n')
        if view is None or view.get('ZZ') is None:
            kind = 'deconvolved ' if decon else 'raw '
            raise RuntimeError('PeakFrame requires the %s2D view in the application data store' % kind)
        return view

    def parameter(self, name, default=''):
        try:
            value = self._legacy.Parse(self._legacy.deconParFile, name)
            return default if value is None else value
        except Exception:
            return default

    def choose_file(self, event, target_box, *, full=True, save=None):
        """Delegate the legacy file chooser without exposing it to PeakFrame callers."""
        return self._legacy.onGetFile(event, target_box, full=full, save=save)

    def cache_external_2d_view(self, path, namespace='peakframe_overlay'):
        return self._legacy.cache_external_2d_view(path, namespace=namespace)

    def run_decon(self, **kwargs):
        return self._legacy.run_decon(**kwargs)

    def save_decon_parameters(self, values):
        """Persist deconvolution parameters through the project parameter file."""
        import os
        from ..project.parameter_store import update_parameter_file
        name = self._legacy.deconParFile
        target = os.path.join(self._legacy.dirBox.GetValue(), name)
        update_parameter_file(target, dict(values), source_path=name)
        state = getattr(self._legacy, 'state', None)
        if state is not None:
            state.set_parameter_file(name)
        return target

    def refresh_status(self):
        status = getattr(self._legacy, 'Status', None)
        if callable(status):
            return status()

    def alias_peak(self, peak, ppm, dimension):
        return self._legacy.alias(peak, ppm, dimension)


    @property
    def axis_limits(self):
        """Physical ppm limits keyed by zero-based spectrum axis."""
        result = []
        for axis in range(max(self.dimension, 2)):
            lo = getattr(self._legacy, 'uc%dmin' % axis, None)
            hi = getattr(self._legacy, 'uc%dmax' % axis, None)
            result.append((lo, hi))
        return tuple(result)

    def axis_scale(self, axis):
        import numpy
        uc = self.unit_converter(axis)
        return numpy.asarray(uc.ppms_scale) if uc is not None else numpy.asarray([])

    def resolve_spectrum_file(self, path):
        resolver = getattr(self._legacy, '_resolve_spec_file', None)
        return resolver(path) if callable(resolver) else path

    def open_reference_peak_list(self, event=None):
        opener = getattr(self._legacy, 'OnButtonReferencePeakList', None)
        return opener(event) if callable(opener) else None

    def set_reference_peaks(self, peaks):
        setter = getattr(self._legacy, 'set_reference_peaks', None)
        if callable(setter):
            result = setter(peaks)
            refresh = getattr(self._legacy, 'refresh_reference_peak_views', None)
            if callable(refresh):
                refresh()
            return result
        return None

    def physical_axis(self, dimension):
        return getattr(self._legacy, 'index%d' % int(dimension))
