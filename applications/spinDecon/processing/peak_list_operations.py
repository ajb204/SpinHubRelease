"""Reusable transformations for peak-list coordinates.

These functions deliberately know nothing about wx frames.  Spectrum-specific
metadata (for example alias widths) is supplied by the caller/controller.
"""


def transpose_2d_peaks(peaks):
    """Transpose the x/y coordinates of a 2D peak list in place."""
    for peak in peaks:
        peak.x, peak.y = peak.y, peak.x
    return peaks


def alias_peak_coordinate(peak, axis, direction, width_ppm):
    """Shift one displayed peak coordinate by one spectral width in place."""
    if axis not in ('x', 'y'):
        raise ValueError("axis must be 'x' or 'y'")
    if direction not in (-1, 1):
        raise ValueError('direction must be -1 or +1')
    width = abs(float(width_ppm))
    if width <= 0:
        raise ValueError('spectral width must be positive')
    before = float(getattr(peak, axis))
    after = before + direction * width
    setattr(peak, axis, after)
    return before, after
