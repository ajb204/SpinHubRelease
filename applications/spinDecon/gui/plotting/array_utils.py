"""Helpers for normalizing small coordinate arrays used by plot views."""

from __future__ import annotations

import numpy as np


def ensure_xy_points(points):
    """Return an (N, 2) float array for x/y coordinate pairs.

    Accepts:
    - empty input -> shape (0, 2)
    - a single pair like [x, y] -> shape (1, 2)
    - an existing (N, 2) array -> returned unchanged (as float)
    - a flat even-length vector -> reshaped to (-1, 2)

    The function is intentionally forgiving because many GUI paths build
    peak-coordinate lists incrementally and may end up with a single peak or
    a flat list after selection/filtering.
    """
    arr = np.asarray(points, dtype=float)
    if arr.size == 0:
        return np.empty((0, 2), dtype=float)
    if arr.ndim == 1:
        if arr.size % 2 == 0:
            arr = arr.reshape(-1, 2)
        else:
            arr = arr.reshape(1, -1)
    return arr


def scatter_xy_points(axis, points, **scatter_kwargs):
    """Scatter x/y points on an axis after normalizing the coordinate shape."""
    pts = ensure_xy_points(points)
    return axis.scatter(pts[:, 0], pts[:, 1], **scatter_kwargs)
