"""Controlled, cancellable peak picking used by the peak-shape fitting GUI."""
from dataclasses import dataclass
import time
import numpy as np


@dataclass
class PeakPickerSettings:
    threshold_fraction: float = 0.05
    max_peaks: int = 5
    min_separation: float = 2.0
    polarity: str = "both"
    max_candidates: int = 10000
    timeout_seconds: float = 10.0
    neighbourhood: int = 1
    region: tuple | None = None
    selection_mode: str = "representative"   # representative, intense
    isolation_radius: int = 5
    representative_low_percentile: float = 35.0
    representative_high_percentile: float = 90.0
    adaptive_threshold: bool = True


@dataclass
class PeakPickerResult:
    maxima: np.ndarray
    values: np.ndarray
    candidate_count: int
    examined_count: int
    status: str = "complete"
    message: str = ""
    representative_count: int = 0
    representative_widths: np.ndarray | None = None
    representative_indices: np.ndarray | None = None
    selected_widths: np.ndarray | None = None


class PeakPicker:
    def __init__(self, data, settings):
        self.data = np.asarray(data)
        self.settings = settings

    def _working_data(self):
        region = self.settings.region
        if region is None:
            return self.data, tuple(0 for _ in range(self.data.ndim))
        view = self.data[region]
        offsets = tuple((sl.start or 0) for sl in region)
        return view, offsets

    def _amplitude(self, array):
        if self.settings.polarity == "positive":
            return array
        if self.settings.polarity == "negative":
            return -array
        return np.abs(array)

    def _threshold(self, array):
        amp = self._amplitude(array)
        finite = amp[np.isfinite(amp)]
        if not finite.size:
            return np.inf
        base = float(np.max(np.abs(self.data))) * float(self.settings.threshold_fraction)
        if not self.settings.adaptive_threshold:
            return base
        # Robust noise estimate.  The median/MAD is insensitive to a modest peak population.
        raw = np.asarray(array, dtype=float)
        med = float(np.median(raw[np.isfinite(raw)]))
        mad = float(np.median(np.abs(raw[np.isfinite(raw)] - med)))
        noise = 1.4826 * mad
        return max(base, 5.0 * noise)

    def _candidate_mask(self, array):
        return self._amplitude(array) > self._threshold(array)

    def estimate_candidates(self):
        array, _ = self._working_data()
        return int(np.count_nonzero(self._candidate_mask(array)))

    def _local_maxima(self, array, mask, cancel_event=None, progress_callback=None):
        """Find maxima without enumerating every above-threshold point when crowded."""
        amp = self._amplitude(array)
        radius = max(1, int(self.settings.neighbourhood))
        # In automatic mode progressively trim only the *search workload*.  The scientific
        # threshold remains a lower bound; this replaces the old fatal max-candidates error.
        work_mask = mask.copy()
        count = int(np.count_nonzero(work_mask))
        limit = max(100, int(self.settings.max_candidates))
        if count > limit and self.settings.adaptive_threshold:
            vals = amp[work_mask]
            cutoff = float(np.partition(vals, max(0, vals.size - limit))[max(0, vals.size - limit)])
            work_mask &= amp >= cutoff
        coords = np.argwhere(work_mask)
        maxima, start, last = [], time.monotonic(), time.monotonic()
        total = len(coords)
        for examined, coord in enumerate(coords, 1):
            if cancel_event is not None and cancel_event.is_set():
                return maxima, examined, "cancelled"
            now = time.monotonic()
            if self.settings.timeout_seconds > 0 and now - start >= self.settings.timeout_seconds:
                return maxima, examined, "timeout"
            slices = tuple(slice(max(0, int(c)-radius), min(array.shape[d], int(c)+radius+1)) for d,c in enumerate(coord))
            window = amp[slices]
            local = np.unravel_index(np.argmax(window), window.shape)
            origin = tuple(sl.start for sl in slices)
            if all(local[d] == int(coord[d])-origin[d] for d in range(array.ndim)):
                maxima.append(tuple(int(c) for c in coord))
            if progress_callback is not None and (now-last > .08 or examined == total):
                progress_callback(examined, total, len(maxima)); last = now
        return maxima, total, "complete"

    def run(self, cancel_event=None, progress_callback=None):
        array, offsets = self._working_data()
        mask = self._candidate_mask(array)
        candidate_count = int(np.count_nonzero(mask))
        empty = (0, self.data.ndim)
        if candidate_count == 0:
            return PeakPickerResult(np.empty(empty,dtype=int), np.array([]), 0, 0, "no_candidates", "No usable peaks exceed the automatic detection threshold.")
        maxima, examined, status = self._local_maxima(array, mask, cancel_event, progress_callback)
        messages = {"complete":"Peak search complete.", "cancelled":"Peak search stopped by user.", "timeout":"Peak search reached the time limit."}
        return self._finish(maxima, array, offsets, candidate_count, examined, status, messages[status])

    @staticmethod
    def _half_widths(amp, coord):
        peak = float(amp[tuple(coord)])
        half = peak * .5
        widths = []
        for d in range(amp.ndim):
            lo = hi = int(coord[d])
            probe = list(coord)
            while lo > 0:
                probe[d] = lo-1
                if amp[tuple(probe)] < half: break
                lo -= 1
            probe[d] = int(coord[d])
            while hi < amp.shape[d]-1:
                probe[d] = hi+1
                if amp[tuple(probe)] < half: break
                hi += 1
            widths.append(max(1, hi-lo+1))
        return np.asarray(widths, dtype=float)

    def _representative_diagnostics(self, array, local, vals):
        """Return widths, representative mask and robust typicality ordering.

        A representative peak must be isolated from another detected maximum and lie in
        the requested intensity percentile band.  Among that clean population we rank
        widths by robust distance from the per-dimension median, so an exceptionally
        sharp peak cannot dominate merely because it is intense.
        """
        amp = self._amplitude(array)
        widths = np.asarray([self._half_widths(amp, c) for c in local], dtype=float)
        if len(local) > 1:
            delta = local[:, None, :] - local[None, :, :]
            distances = np.linalg.norm(delta, axis=2)
            distances[distances == 0] = np.inf
            nearest = np.min(distances, axis=1)
        else:
            nearest = np.full(len(local), np.inf)
        radius = max(2.0, float(self.settings.isolation_radius))
        isolated = nearest >= radius
        lo = np.percentile(vals, np.clip(self.settings.representative_low_percentile, 0, 100))
        hi = np.percentile(vals, np.clip(self.settings.representative_high_percentile, 0, 100))
        intensity_ok = (vals >= lo) & (vals <= hi)
        representative = isolated & intensity_ok
        # If a very sparse spectrum has no peak in the percentile band, isolation is the
        # scientifically more important criterion and provides a graceful fallback.
        pool = np.flatnonzero(representative)
        if not pool.size:
            pool = np.flatnonzero(isolated)
            representative = isolated.copy()
        if not pool.size:
            pool = np.arange(len(local))
            representative = np.ones(len(local), dtype=bool)
        centre = np.median(widths[pool], axis=0)
        mad = np.median(np.abs(widths[pool] - centre), axis=0)
        scale = np.where(mad > 0, 1.4826 * mad, np.maximum(centre * .20, 1.0))
        distance = np.sqrt(np.mean(((widths - centre) / scale) ** 2, axis=1))
        order = pool[np.argsort(distance[pool])]
        return widths, representative, order

    def _finish(self, maxima, array, offsets, candidate_count, examined, status, message):
        if not maxima:
            return PeakPickerResult(np.empty((0,self.data.ndim),dtype=int), np.array([]), candidate_count, examined, status, message)
        local = np.asarray(maxima,dtype=int).reshape((-1,self.data.ndim))
        vals = np.asarray([self._amplitude(array)[tuple(c)] for c in local],dtype=float)
        widths = np.asarray([self._half_widths(self._amplitude(array), c) for c in local], dtype=float)
        representative = np.ones(len(local), dtype=bool)
        if str(self.settings.selection_mode).lower() == "intense":
            order = np.argsort(vals)[::-1]
        else:
            widths, representative, order = self._representative_diagnostics(array, local, vals)
        selected=[]; min_sep=max(0.,float(self.settings.min_separation))
        for idx in order:
            coord=local[idx]
            if min_sep and any(np.linalg.norm(coord-local[j]) < min_sep for j in selected): continue
            selected.append(int(idx))
            if len(selected)>=max(1,int(self.settings.max_peaks)): break
        selected=np.asarray(selected,dtype=int)
        chosen_local=local[selected]; chosen_vals=vals[selected]; chosen_widths=widths[selected]
        global_coords=chosen_local+np.asarray(offsets,dtype=int)
        ascending=np.argsort(chosen_vals)
        kind = "representative isolated" if str(self.settings.selection_mode).lower() != "intense" else "most intense"
        rep_idx = np.flatnonzero(representative)
        return PeakPickerResult(global_coords[ascending], chosen_vals[ascending], candidate_count, examined, status,
                                message + " Selected %d %s peak%s." % (len(chosen_local), kind, "" if len(chosen_local)==1 else "s"),
                                representative_count=int(rep_idx.size), representative_widths=widths[rep_idx],
                                representative_indices=(local[rep_idx] + np.asarray(offsets,dtype=int)),
                                selected_widths=chosen_widths[ascending])

