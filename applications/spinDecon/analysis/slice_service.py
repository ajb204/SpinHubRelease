from .full_peak_service import FullPeakListService

"""Application boundary for slice viewers during legacy GUI migration."""


class SliceService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def labels(self):
        return tuple(getattr(self._legacy, "labb", ()) or ())

    @property
    def deconvolution_enabled(self):
        return bool(getattr(self._legacy, "DECON", 0))

    @property
    def peaks(self):
        return getattr(self._legacy, "peak", [])

    def reference_peaks(self):
        getter = getattr(self._legacy, "get_reference_peaks", None)
        return getter() if callable(getter) else self.peaks

    def threshold(self):
        return float(self._legacy.dmax) * float(self._legacy.threshBox.GetValue())

    def reference_1d_view(self, selection):
        return self._legacy.get_reference_1d_view(selection)

    @property
    def symmetry_enabled(self):
        control = getattr(self._legacy, "cb_grid", None)
        return bool(control.IsChecked()) if control is not None else False

    @property
    def spectrum_path(self):
        return str(getattr(self._legacy, "spectrumfile", "") or "")

    @property
    def dimension(self):
        return int(getattr(self._legacy, 'dim', len(self.labels)) or len(self.labels))

    def full_peak_payload(self):
        getter = getattr(self._legacy, 'get_full_peak_payload', None)
        return getter() if callable(getter) else {}

    def clear_peak_selection(self, redraw_full3d=True):
        clear = getattr(self._legacy, 'clear_peak_selection', None)
        return clear(redraw_full3d=redraw_full3d) if callable(clear) else None

    def select_full_peak(self, name, source_view=None, source_pane=None):
        select = getattr(self._legacy, 'select_full_peak', None)
        if callable(select):
            return select(name, source_view=source_view, source_pane=source_pane)
        return None

    def save_full_peak_records(self, records, rows, source_path=None):
        self._legacy.store.save_peak_list('full', peaks=records, records=records, rows=rows,
                                          dimension=self.dimension, source_path=source_path)
        rebuild = getattr(self._legacy, '_rebuild_projected_peak_lists', None)
        if callable(rebuild): rebuild()
        notify = getattr(self._legacy, '_notify_analysis_changed', None)
        if callable(notify): notify()
        refresh = getattr(self._legacy, 'refresh_full_peak_list_viewers', None)
        if callable(refresh): refresh()

    def redraw_full3d_if_open(self):
        viewer_getter = getattr(self._legacy, '_full3d_viewer', None)
        viewer = viewer_getter(ensure=False) if callable(viewer_getter) else None
        if viewer is not None:
            try: viewer.draw_figure(keepaxes=True)
            except Exception: pass

    def projection_view(self, *args, **kwargs):
        return self._legacy.get_projection_view(*args, **kwargs)

    @property
    def data(self):
        return getattr(self._legacy, 'data', None)

    @property
    def decon_data(self):
        return getattr(self._legacy, 'datadec', None)

    @property
    def datadec(self):
        """Compatibility alias for legacy slice viewers during migration."""
        return self.decon_data

    def slice_meshes(self, *, orth, pt_c=None, pt_h=None, pt_h_max=None, pt_h_min=None,
                     pt_c_max=None, pt_c_min=None):
        """Return calibrated X/Y meshes for a 2D strip from the active dataset."""
        if not orth:
            xs = getattr(self._legacy, 'XX')[:, int(pt_c), int(pt_h_max):int(pt_h_min)].transpose()
            ys = getattr(self._legacy, 'ZZ')[:, int(pt_c), int(pt_h_max):int(pt_h_min)].transpose()
        else:
            xs = getattr(self._legacy, 'XX')[:, int(pt_c_max):int(pt_c_min), int(pt_h)].transpose()
            ys = getattr(self._legacy, 'YY')[:, int(pt_c_max):int(pt_c_min), int(pt_h)].transpose()
        return xs, ys

    @property
    def peak_indices(self):
        return getattr(self._legacy, 'pkIdx', [])

    def axis(self, dimension):
        return getattr(self._legacy, 'index%d' % int(dimension))

    def mesh(self, name):
        return getattr(self._legacy, str(name))

    def peak(self, index):
        return self.peaks[int(index)]

    def peak_names(self):
        return [str(getattr(pk, 'name', '')) for pk in self.peaks]

    def peak_shape_width(self, dimension, multiplier=9.0):
        control = getattr(self._legacy, 'sig%dBox' % int(dimension), None)
        if control is None:
            return None
        return float(control.GetValue()) * float(multiplier)

    def sample(self, indices, decon=False):
        array = self.decon_data if decon else self.data
        return float(array[tuple(int(v) for v in indices)])

    def projection(self, x_label, y_label, decon=False, transpose='n'):
        return self._legacy.get_projection_view(x_label, y_label, decon=decon, transpose=transpose)

    def full_peak_slice_color(self, delta):
        getter = getattr(self._legacy, 'full_peak_slice_color', None)
        return getter(delta) if callable(getter) else None

    @property
    def decon_parameter_file(self):
        return getattr(self._legacy, 'deconParFile', None)

    def open_full_peak_list(self, event=None):
        return FullPeakListService(self._legacy).open(event)

    @property
    def spectral_dimension(self):
        topology = getattr(self._legacy, '_active_topology', None)
        if callable(topology):
            return int(topology().spectral_dim_count)
        return int(getattr(self._legacy, 'dim', self.dimension))

    def refresh_status(self):
        callback = getattr(self._legacy, 'Status', None)
        return callback() if callable(callback) else None

    @property
    def decon_enabled(self):
        return bool(getattr(self._legacy, 'DECON', False))

    def label(self, dimension):
        labels = self.labels
        index = int(dimension)
        return str(labels[index]) if 0 <= index < len(labels) else ''

    def axis_limits(self, dimension):
        index = int(dimension)
        return (float(getattr(self._legacy, 'uc%dmin' % index)),
                float(getattr(self._legacy, 'uc%dmax' % index)))

    @property
    def max_intensity(self):
        return float(getattr(self._legacy, 'dmax', 0.0))
