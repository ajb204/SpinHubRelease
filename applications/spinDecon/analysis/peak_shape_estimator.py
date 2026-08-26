"""Robust deconvolution-filter estimation from representative NMR peaks.

This intentionally estimates a conservative *envelope* rather than a least-squares
peak fit.  Each side of each selected peak is treated independently so an overlapped
wing can be discarded while the clean wing remains useful.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class FilterEstimate:
    gaussian_fwhm: float
    lorentzian_fwhm: float
    voigt_fraction: float
    measured_fwhm: float
    recommended_fwhm: float
    clean_sides: int
    rejected_sides: int


def _crossing(distance, intensity, level):
    """Linear-interpolated first outward crossing of ``level``."""
    for i in range(1, len(intensity)):
        if intensity[i] <= level <= intensity[i - 1]:
            dy = intensity[i - 1] - intensity[i]
            f = 0.0 if dy == 0 else (intensity[i - 1] - level) / dy
            return float(distance[i - 1] + f * (distance[i] - distance[i - 1]))
    return None


def _analyse_side(axis, trace, centre, direction):
    """Return level radii for one wing, rejecting obvious shoulders/overlap."""
    sign = 1.0 if trace[centre] >= 0 else -1.0
    y = np.asarray(trace, float) * sign
    peak = float(y[centre])
    if not np.isfinite(peak) or peak <= 0:
        return None
    # A local baseline is intentionally conservative; only central crossings matter.
    edge_n = max(2, min(8, len(y) // 10))
    baseline = float(np.median(np.r_[y[:edge_n], y[-edge_n:]]))
    amp = peak - baseline
    if amp <= 0:
        return None
    yn = (y - baseline) / amp
    indices = np.arange(centre, len(y)) if direction > 0 else np.arange(centre, -1, -1)
    vals = yn[indices]
    dist = np.abs(np.asarray(axis, float)[indices] - float(axis[centre]))

    # Stop a wing when a meaningful rise/shoulder occurs after leaving the apex.
    kept = [0]
    rises = 0
    for j in range(1, len(vals)):
        if not np.isfinite(vals[j]):
            break
        if vals[j] > vals[j - 1] + 0.035 and vals[j - 1] < 0.85:
            rises += 1
        else:
            rises = 0
        if rises >= 2:
            break
        kept.append(j)
        if vals[j] <= 0.20:
            break
    vals = vals[kept]
    dist = dist[kept]
    if len(vals) < 3:
        return None
    crossings = {level: _crossing(dist, vals, level) for level in (0.75, 0.60, 0.50, 0.35)}
    # Half height is essential; 0.35 is needed to infer tail character.
    if crossings[0.50] is None or crossings[0.35] is None or crossings[0.50] <= 0:
        return None
    return crossings


def _pseudo_voigt_radius_ratio(level, fraction):
    # Equal Gaussian/Lorentzian FWHM, with W=1.  Binary search radius.
    lo, hi = 0.0, 5.0
    for _ in range(60):
        x = (lo + hi) * 0.5
        g = np.exp(-4.0 * np.log(2.0) * x * x)
        l = 1.0 / (1.0 + 4.0 * x * x)
        value = (1.0 - fraction) * g + fraction * l
        if value > level:
            lo = x
        else:
            hi = x
    return (lo + hi) * 0.5


def _estimate_fraction(sides):
    observed = []
    for c in sides:
        r50 = c[0.50]
        if c.get(0.35) is not None:
            observed.append((0.35, c[0.35] / r50))
        if c.get(0.75) is not None:
            observed.append((0.75, c[0.75] / r50))
        if c.get(0.60) is not None:
            observed.append((0.60, c[0.60] / r50))
    if not observed:
        return 0.5
    best_f, best_loss = 0.5, np.inf
    for f in np.linspace(0.0, 1.0, 101):
        r50 = _pseudo_voigt_radius_ratio(0.50, f)
        errors = []
        for level, ratio in observed:
            predicted = _pseudo_voigt_radius_ratio(level, f) / r50
            errors.append(abs(predicted - ratio))
        loss = float(np.median(errors))
        if loss < best_loss:
            best_f, best_loss = float(f), loss
    return best_f


def estimate_dimension(data, axis, peaks, dimension, widen_factor=1.10, width_quantile=0.65):
    """Estimate one dimension's safe UniDec filter from selected peak coordinates."""
    data = np.asarray(data)
    axis = np.asarray(axis, float)
    sides, rejected = [], 0
    for coord in np.asarray(peaks, dtype=int):
        selector = [int(v) for v in coord]
        centre = selector[dimension]
        selector[dimension] = slice(None)
        trace = np.asarray(data[tuple(selector)], float)
        for direction in (-1, 1):
            result = _analyse_side(axis, trace, centre, direction)
            if result is None:
                rejected += 1
            else:
                sides.append(result)
    if not sides:
        raise ValueError("No clean peak wings reached the half-height and 35% crossings")
    # Each clean wing supplies a half-width.  Use an upper-middle quantile and then
    # deliberately widen it: under-wide filters are more damaging in deconvolution.
    fwhms = np.asarray([2.0 * s[0.50] for s in sides], float)
    measured = float(np.median(fwhms))
    recommended = float(np.quantile(fwhms, width_quantile) * widen_factor)
    fraction = _estimate_fraction(sides)
    return FilterEstimate(recommended, recommended, fraction, measured, recommended,
                          len(sides), rejected)


def estimate_filter_shape(data, axes, peaks, link_widths=True):
    # The first implementation intentionally uses the stable equal-width model.  When
    # widths are unlinked we still provide the same robust starting estimate; users can
    # subsequently tune Gaussian and Lorentzian widths independently in the GUI.
    return [estimate_dimension(data, axes[d], peaks, d) for d in range(np.asarray(data).ndim)]


def estimate_level_radius(data, axes, peaks, level=0.10):
    """Estimate conservative extraction radii at ``level`` of peak height.

    Each peak wing is inspected independently. Wings that develop a sustained
    shoulder/rise before reaching the requested level are rejected, matching
    the overlap-avoidance philosophy used by the peak-shape Fit action.
    The returned radius for each dimension is the median of clean crossings.
    """
    data = np.asarray(data)
    peaks = np.asarray(peaks, dtype=int)
    if not len(peaks):
        raise ValueError("Find representative peaks before guessing extraction radii")
    radii=[]; diagnostics=[]
    for dimension, axis in enumerate(axes):
        axis=np.asarray(axis, float); clean=[]; rejected=0
        for coord in peaks:
            selector=[int(v) for v in coord]; centre=selector[dimension]; selector[dimension]=slice(None)
            trace=np.asarray(data[tuple(selector)], float)
            sign=1.0 if trace[centre] >= 0 else -1.0
            y=trace*sign; peak=float(y[centre])
            edge_n=max(2, min(8, len(y)//10)); baseline=float(np.median(np.r_[y[:edge_n], y[-edge_n:]]))
            amp=peak-baseline
            if not np.isfinite(amp) or amp <= 0:
                rejected += 2; continue
            yn=(y-baseline)/amp
            for direction in (-1,1):
                indices=np.arange(centre, len(y)) if direction>0 else np.arange(centre,-1,-1)
                vals=yn[indices]; dist=np.abs(axis[indices]-axis[centre]); kept=[0]; rises=0
                for j in range(1,len(vals)):
                    if not np.isfinite(vals[j]): break
                    if vals[j] > vals[j-1] + 0.035 and vals[j-1] < 0.85: rises += 1
                    else: rises=0
                    if rises >= 2: break
                    kept.append(j)
                    if vals[j] <= level: break
                crossing=_crossing(dist[kept], vals[kept], level) if len(kept)>=2 else None
                if crossing is None: rejected += 1
                else: clean.append(crossing)
        if not clean:
            raise ValueError("No clean F%d peak wings reached %.0f%% intensity; overlap may obscure the radius" % (dimension+1, level*100.0))
        radii.append(float(np.median(clean)))
        diagnostics.append((len(clean), rejected))
    return radii, diagnostics
