"""Diagnostics for 3D deconvolution restricted by a 2D reference peak list.

Measurements are written onto the existing shared peak records.  Full 3D
records receive per-peak measurements in ``record['analysis']``; reference 2D
peaks receive source-level summaries in ``peak.analysis``.  ``DataStore.analysis``
contains run parameters/status only.
"""
from __future__ import annotations

from typing import Any, Callable
import numpy as np


def _nearest_index(scale: Any, value: float) -> int:
    arr = np.asarray(scale, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("Axis scale is empty or invalid")
    return int(np.nanargmin(np.abs(arr - float(value))))


def _axis_scale(store: Any, dim: int):
    uc = getattr(store, f"uc{dim}", None)
    scale = getattr(uc, "ppms_scale", None) if uc is not None else None
    if scale is None:
        scale = getattr(store, f"index{dim}", None)
    return np.asarray(scale, dtype=float)


def _record_array_index(store: Any, record: dict[str, Any]) -> tuple[int, int, int]:
    labels = list(getattr(store, "labb", []) or [])
    if len(labels) < 3:
        raise ValueError("3D axis labels are unavailable")
    axes = dict(record.get("axis_values") or {})
    indices = []
    for dim, label in enumerate(labels[:3]):
        key = str(label)
        if key not in axes:
            raise ValueError(f"Peak has no coordinate for axis {key}")
        indices.append(_nearest_index(_axis_scale(store, dim), axes[key]))
    return tuple(indices)


def _reference_indices(store: Any, peak: Any) -> tuple[int, int]:
    if hasattr(peak, "indexJ") and hasattr(peak, "indexK"):
        return int(peak.indexJ), int(peak.indexK)
    # Reference x/y correspond to array dimensions 2/1 in the established
    # data[:, indexJ, indexK] bore convention.
    return _nearest_index(_axis_scale(store, 1), float(peak.y)), _nearest_index(_axis_scale(store, 2), float(peak.x))


def _safe_noise(store: Any) -> float:
    try:
        value = abs(float(getattr(store, "noiseVal", 0.0)))
        return value if np.isfinite(value) else 0.0
    except (TypeError, ValueError):
        return 0.0


def _ppm_offset(scale: np.ndarray, a: int, b: int) -> float:
    try:
        return float(scale[b] - scale[a])
    except (IndexError, TypeError, ValueError):
        return float("nan")


def analyse_restricted_3d(
    store: Any,
    associate: Callable[[dict[str, Any]], tuple[int | None, Any | None]],
    *,
    overlap_fraction: float = 0.05,
    z_search_radius: int = 5,
    maximum_fraction_warning: float = 0.90,
    xy_search_radius_j: int = 3,
    xy_search_radius_k: int = 3,
    xy_displacement_warning: float = 1.5,
    residual_sigma_threshold: float = 3.0,
    residual_fraction_warning: float = 0.05,
    residual_xy_radius_j: int = 1,
    residual_xy_radius_k: int = 1,
    residual_min_support_traces: int = 3,
) -> dict[str, Any]:
    """Measure overlap, localisation and unexplained trace intensity."""
    data = np.asarray(getattr(store, "data", None))
    if data.ndim != 3:
        raise ValueError("Restricted 3D diagnostics require loaded 3D raw data")
    datadec = getattr(store, "datadec", None)
    datadec = np.asarray(datadec) if datadec is not None else None
    if datadec is not None and datadec.shape != data.shape:
        datadec = None

    references = list(store.get_peak_list("reference").get("peaks") or [])
    full_payload = store.get_peak_list("full")
    records = list(full_payload.get("records") or full_payload.get("peaks") or [])
    noise = _safe_noise(store)
    zscale, jscale, kscale = (_axis_scale(store, i) for i in range(3))

    for peak in references:
        if not hasattr(peak, "analysis") or not isinstance(peak.analysis, dict):
            peak.analysis = {}
        peak.analysis["restricted_3d"] = {
            "full_peak_indices": [], "count": 0,
            "overlap": {"overlapped_count": 0, "worst_fraction": None},
            "localisation": {"warning_count": 0, "clean_warning_count": 0,
                             "worst_maximum_fraction": None, "max_xy_displacement": None},
            "residual": {}, "classification": {},
        }

    measured = unmatched = invalid = 0
    for record_index, record in enumerate(records):
        diag = record.setdefault("analysis", {}).setdefault("restricted_3d", {})
        ref_index, ref_peak = associate(record)
        if ref_peak is None or ref_index is None:
            diag.update({"status": "unmatched", "reference_index": None})
            unmatched += 1
            continue
        try:
            i, _, _ = _record_array_index(store, record)
            j, k = _reference_indices(store, ref_peak)
            raw_intensity = float(data[i, j, k])
        except (ValueError, TypeError, IndexError):
            diag.update({"status": "invalid_coordinate", "reference_index": int(ref_index)})
            invalid += 1
            continue

        try:
            peak_intensity = float(record.get("intensity"))
        except (TypeError, ValueError):
            peak_intensity = None
        denominator = abs(raw_intensity)
        fraction = None
        if peak_intensity is not None and denominator > max(noise, np.finfo(float).eps):
            fraction = abs(peak_intensity) / denominator

        # Z localisation: how close the picked z position is to the strongest
        # point on the *restricted reference trace* in a local window.
        lo, hi = max(0, i-int(z_search_radius)), min(data.shape[0], i+int(z_search_radius)+1)
        local = np.abs(data[lo:hi, j, k])
        zmax_i = lo + int(np.nanargmax(local)) if local.size else i
        zmax = float(abs(data[zmax_i, j, k]))
        max_fraction = abs(raw_intensity) / zmax if zmax > max(noise, np.finfo(float).eps) else None

        # XY localisation must be measured on the Z slice selected by the
        # localisation test above, not blindly on the originally picked slice.
        # The latter can be off the local Z maximum (indeed that is one of the
        # conditions this diagnostic is designed to detect), and searching that
        # wrong plane can collapse otherwise distinct XY proposals onto a ridge.
        # Search around the original 2D source because the scientific question is
        # whether that restricted bore was placed at the correct XY location.
        xy_search_i = int(zmax_i)
        j0, j1 = max(0, j-int(xy_search_radius_j)), min(data.shape[1], j+int(xy_search_radius_j)+1)
        k0, k1 = max(0, k-int(xy_search_radius_k)), min(data.shape[2], k+int(xy_search_radius_k)+1)
        plane = np.abs(data[xy_search_i, j0:j1, k0:k1])
        if plane.size:
            # Keep the discrete maximum for display/debugging, but do not use it
            # as the proposed source position.  In the difficult case this tool
            # is meant to diagnose, two unresolved XY sources often form one
            # broad ridge: every Z slice can therefore have the *same* grid
            # maximum even though different 3D peaks lean to opposite sides.
            # A background-subtracted first moment preserves that sub-grid lean.
            relj, relk = np.unravel_index(int(np.nanargmax(plane)), plane.shape)
            best_j, best_k = j0 + int(relj), k0 + int(relk)
            finite = plane[np.isfinite(plane)]
            baseline = float(np.nanmin(finite)) if finite.size else 0.0
            weights = np.maximum(np.where(np.isfinite(plane), plane, baseline) - baseline, 0.0)
            # Squaring suppresses the broad low-level skirt without forcing the
            # estimate to snap to a (possibly shared) discrete maximum.
            weights = weights * weights
            total = float(np.sum(weights))
            if total > np.finfo(float).eps:
                jj, kk = np.indices(plane.shape, dtype=float)
                centroid_j = float(j0 + np.sum(weights * jj) / total)
                centroid_k = float(k0 + np.sum(weights * kk) / total)
            else:
                centroid_j, centroid_k = float(best_j), float(best_k)
        else:
            best_j, best_k = j, k
            centroid_j, centroid_k = float(j), float(k)
        # Localisation warnings and proposals use the continuous centroid.  This
        # is intentionally sensitive to shoulders/leaning peaks where no clean
        # secondary XY maximum exists.
        dj, dk = centroid_j-float(j), centroid_k-float(k)
        xy_disp = float(np.hypot(dj, dk))

        overlap_flag = fraction is not None and fraction <= float(overlap_fraction)
        max_flag = max_fraction is not None and max_fraction < float(maximum_fraction_warning)
        xy_flag = xy_disp > float(xy_displacement_warning)
        loc_flag = bool(max_flag or xy_flag)

        diag.update({
            "status": "measured", "reference_index": int(ref_index),
            "array_index": (int(i), int(j), int(k)),
            "overlap": {"raw_intensity": raw_intensity, "peak_intensity": peak_intensity,
                        "fraction": fraction, "is_overlapped": bool(overlap_flag)},
            "localisation": {
                "maximum_fraction": max_fraction, "z_max_index": int(zmax_i),
                "z_offset_points": int(zmax_i-i), "z_offset_ppm": _ppm_offset(zscale, i, zmax_i),
                "xy_search_z_index": int(xy_search_i),
                "xy_search_z_ppm": float(zscale[xy_search_i]),
                "xy_max_index": (int(best_j), int(best_k)),
                "xy_centroid_index": (float(centroid_j), float(centroid_k)),
                "xy_offset_points": (float(dj), float(dk)),
                "xy_offset_ppm": (float(np.interp(centroid_j, np.arange(jscale.size), jscale) - jscale[j]),
                                  float(np.interp(centroid_k, np.arange(kscale.size), kscale) - kscale[k])),
                "xy_displacement_points": xy_disp, "maximum_warning": bool(max_flag),
                "xy_warning": bool(xy_flag), "is_warning": loc_flag,
            },
        })

        src = ref_peak.analysis["restricted_3d"]
        src["full_peak_indices"].append(record_index); src["count"] += 1
        ov = src["overlap"]
        if overlap_flag: ov["overlapped_count"] += 1
        if fraction is not None:
            ov["worst_fraction"] = fraction if ov["worst_fraction"] is None else min(ov["worst_fraction"], fraction)
        loc = src["localisation"]
        if loc_flag: loc["warning_count"] += 1
        if loc_flag and not overlap_flag: loc["clean_warning_count"] += 1
        if max_fraction is not None:
            old = loc["worst_maximum_fraction"]
            loc["worst_maximum_fraction"] = max_fraction if old is None else min(old, max_fraction)
        old = loc["max_xy_displacement"]
        loc["max_xy_displacement"] = xy_disp if old is None else max(old, xy_disp)
        measured += 1

    # Source-level residual / unexplained signal.  A missed peak should be a
    # coherent feature in several neighbouring XY bores, whereas baseline/noise
    # excursions are much less likely to recur at the same Z position.
    for ref_index, peak in enumerate(references):
        src = peak.analysis["restricted_3d"]
        residual = src["residual"]
        try:
            j, k = _reference_indices(store, peak)
            raw_trace = np.asarray(data[:, j, k], dtype=float)
            if datadec is None:
                residual.update({"status": "no_decon_data", "fraction": None, "max_sigma": None})
            else:
                rj, rk = max(0, int(residual_xy_radius_j)), max(0, int(residual_xy_radius_k))
                j0, j1 = max(0, j-rj), min(data.shape[1], j+rj+1)
                k0, k1 = max(0, k-rk), min(data.shape[2], k+rk+1)
                raw_block = np.asarray(data[:, j0:j1, k0:k1], dtype=float)
                dec_block = np.asarray(datadec[:, j0:j1, k0:k1], dtype=float)
                res_block = raw_block - dec_block
                n_traces = int(res_block.shape[1] * res_block.shape[2])
                min_support = max(1, min(int(residual_min_support_traces), n_traces))
                threshold = float(residual_sigma_threshold) * noise

                if noise > 0:
                    excess = np.maximum(np.abs(res_block) - threshold, 0.0)
                    support = np.sum(excess > 0.0, axis=(1, 2))
                else:
                    excess = np.abs(res_block)
                    support = np.sum(excess > 0.0, axis=(1, 2))

                persistent_z = support >= min_support
                # Only intensity at Z positions supported by several nearby
                # traces contributes.  Use the median excess across supporting
                # traces so a single extreme bore cannot dominate the score.
                persistent_profile = np.zeros(data.shape[0], dtype=float)
                for zi in np.flatnonzero(persistent_z):
                    vals = excess[zi].ravel()
                    vals = vals[vals > 0.0]
                    if vals.size >= min_support:
                        persistent_profile[zi] = float(np.median(vals))

                raw_area = float(np.sum(np.abs(raw_trace)))
                res_fraction = float(np.sum(persistent_profile) / raw_area) if raw_area > 0 else 0.0
                central_res = np.asarray(res_block[:, j-j0, k-k0], dtype=float)
                max_sigma = float(np.max(np.abs(central_res)) / noise) if noise > 0 else None
                coherent_max_sigma = float(np.max(persistent_profile) / noise + float(residual_sigma_threshold)) if noise > 0 and np.any(persistent_profile > 0) else 0.0 if noise > 0 else None
                residual.update({
                    "status": "measured", "fraction": res_fraction,
                    "max_sigma": max_sigma, "coherent_max_sigma": coherent_max_sigma,
                    "max_abs": float(np.max(np.abs(central_res))),
                    "neighbourhood_shape": (int(j1-j0), int(k1-k0)),
                    "trace_count": n_traces, "min_support_traces": min_support,
                    "persistent_z_count": int(np.count_nonzero(persistent_z)),
                    "max_support_traces": int(np.max(support)) if support.size else 0,
                    "warning": bool(res_fraction > float(residual_fraction_warning)),
                })
        except (ValueError, TypeError, IndexError):
            residual.update({"status": "invalid_trace", "fraction": None, "max_sigma": None})

    result = classify_restricted_3d(
        store, overlap_fraction=overlap_fraction,
        maximum_fraction_warning=maximum_fraction_warning,
        xy_displacement_warning=xy_displacement_warning,
        residual_fraction_warning=residual_fraction_warning,
    )
    result.update({"measured_records": measured, "unmatched_records": unmatched, "invalid_records": invalid})
    params = {
        "overlap_fraction": float(overlap_fraction), "z_search_radius": int(z_search_radius),
        "maximum_fraction_warning": float(maximum_fraction_warning),
        "xy_search_radius_j": int(xy_search_radius_j), "xy_search_radius_k": int(xy_search_radius_k),
        "xy_displacement_warning": float(xy_displacement_warning),
        "residual_sigma_threshold": float(residual_sigma_threshold),
        "residual_fraction_warning": float(residual_fraction_warning),
        "residual_xy_radius_j": int(residual_xy_radius_j), "residual_xy_radius_k": int(residual_xy_radius_k),
        "residual_min_support_traces": int(residual_min_support_traces),
    }
    store.analysis["restricted_3d"] = {"status": "complete", "parameters": params, "summary": dict(result)}
    return result


def classify_restricted_3d(store: Any, *, overlap_fraction=0.05,
                           maximum_fraction_warning=0.90,
                           xy_displacement_warning=1.5,
                           residual_fraction_warning=0.05) -> dict[str, Any]:
    """Reclassify stored measurements without rescanning the 3D cube."""
    references = list(store.get_peak_list("reference").get("peaks") or [])
    records = list(store.get_peak_list("full").get("records") or store.get_peak_list("full").get("peaks") or [])
    overlapped_records = localisation_records = 0
    for record in records:
        diag = record.setdefault("analysis", {}).get("restricted_3d", {})
        ov, loc = diag.get("overlap", {}), diag.get("localisation", {})
        fraction = ov.get("fraction")
        ov["is_overlapped"] = bool(fraction is not None and float(fraction) <= float(overlap_fraction))
        mf, xd = loc.get("maximum_fraction"), loc.get("xy_displacement_points")
        loc["maximum_warning"] = bool(mf is not None and float(mf) < float(maximum_fraction_warning))
        loc["xy_warning"] = bool(xd is not None and float(xd) > float(xy_displacement_warning))
        loc["is_warning"] = bool(loc.get("maximum_warning") or loc.get("xy_warning"))
        overlapped_records += int(ov.get("is_overlapped", False))
        localisation_records += int(loc.get("is_warning", False))

    overlap_sources = localisation_sources = residual_sources = missing_xy_sources = attention_sources = 0
    for peak in references:
        src = getattr(peak, "analysis", {}).get("restricted_3d", {})
        inds = src.get("full_peak_indices", [])
        source_records = [records[i] for i in inds if isinstance(i, int) and 0 <= i < len(records)]
        overlap_count = sum(bool(r.get("analysis", {}).get("restricted_3d", {}).get("overlap", {}).get("is_overlapped")) for r in source_records)
        loc_count = sum(bool(r.get("analysis", {}).get("restricted_3d", {}).get("localisation", {}).get("is_warning")) for r in source_records)
        clean_loc_count = sum(bool(r.get("analysis", {}).get("restricted_3d", {}).get("localisation", {}).get("is_warning")) and not bool(r.get("analysis", {}).get("restricted_3d", {}).get("overlap", {}).get("is_overlapped")) for r in source_records)
        src.setdefault("overlap", {})["overlapped_count"] = overlap_count
        src.setdefault("localisation", {})["warning_count"] = loc_count
        src["localisation"]["clean_warning_count"] = clean_loc_count
        residual = src.setdefault("residual", {})
        rf = residual.get("fraction")
        residual["warning"] = bool(rf is not None and float(rf) > float(residual_fraction_warning))
        # Multiple clean, mislocalised 3D peaks are strong evidence that the
        # 2D source model deserves inspection; one is a weaker suggestion.
        possible_missing_xy = clean_loc_count >= 1
        strong_missing_xy = clean_loc_count >= 2
        flags = {
            "overlap": overlap_count > 0, "localisation": loc_count > 0,
            "unexplained": bool(residual.get("warning")),
            "possible_missing_xy": possible_missing_xy,
            "strong_missing_xy": strong_missing_xy,
        }
        flags["needs_attention"] = any((flags["overlap"], flags["localisation"], flags["unexplained"]))
        src["classification"] = flags
        overlap_sources += int(flags["overlap"]); localisation_sources += int(flags["localisation"])
        residual_sources += int(flags["unexplained"]); missing_xy_sources += int(flags["possible_missing_xy"])
        attention_sources += int(flags["needs_attention"])
    return {"reference_sources": len(references), "full_records": len(records),
            "overlapped_records": overlapped_records, "localisation_records": localisation_records,
            "overlapped_sources": overlap_sources, "localisation_sources": localisation_sources,
            "residual_sources": residual_sources, "possible_missing_xy_sources": missing_xy_sources,
            "attention_sources": attention_sources}


# Backward-compatible v1 entry points.
def analyse_overlap(store, associate, *, overlap_fraction=0.05):
    return analyse_restricted_3d(store, associate, overlap_fraction=overlap_fraction)


def classify_overlap(store, *, overlap_fraction=0.05):
    return classify_restricted_3d(store, overlap_fraction=overlap_fraction)
