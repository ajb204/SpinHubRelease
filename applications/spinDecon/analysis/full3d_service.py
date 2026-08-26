"""Application boundary for Full-3D scientific view operations.

The first implementation delegates to the legacy NMR workspace.  Keeping that
dependency here prevents the Full3D GUI from depending on deconFrame's API and
allows the numerical implementation to move out independently later.
"""

class Full3DService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    def view_spec(self, bore_dim):
        return self._legacy.get_full3d_view_spec(bore_dim)

    def slice_view(self, bore_dim, slice_index):
        return self._legacy.get_full3d_slice_view(bore_dim, slice_index)

    def has_full_peaks(self):
        return bool((self._legacy.get_full_peak_payload().get('records') or []))

    def cross_sections(self, bore_dim, slice_index, x_ppm, y_ppm):
        return self._legacy.get_full3d_cross_sections(bore_dim, slice_index, x_ppm, y_ppm)

    def intensity_limits(self):
        return self._legacy.get_full3d_intensity_limits()

    def peak_overlay(self, bore_dim, slice_index):
        return self._legacy.get_full3d_peak_overlay(bore_dim, slice_index)

    def clear_peak_selection(self, redraw_full3d=False):
        return self._legacy.clear_peak_selection(redraw_full3d=redraw_full3d)
