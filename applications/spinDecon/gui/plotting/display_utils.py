"""Shared helpers for high-frequency plot updates and cursor handling."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from matplotlib.widgets import Cursor, MultiCursor


@dataclass(frozen=True)
class DisplayMode:
    """Minimal display preferences used by interactive spectrum views."""

    use_blit: bool = True
    mouse_tracking: bool = True


def make_cursor(axis, *, horizOn=True, vertOn=True, useblit=True, **kwargs):
    """Create a matplotlib Cursor with a blit-friendly default."""
    return Cursor(axis, horizOn=horizOn, vertOn=vertOn, useblit=useblit, **kwargs)


def make_multi_cursor(canvas, axes, *, horizOn=True, vertOn=True, useblit=True, **kwargs):
    """Create a matplotlib MultiCursor with a blit-friendly default."""
    return MultiCursor(canvas, axes, horizOn=horizOn, vertOn=vertOn, useblit=useblit, **kwargs)


def blit_artists(canvas, axes, background, artists: Sequence):
    """Restore a cached background, draw artists, and blit the active axes.

    If no background is available, the function still draws the artists and
    falls back to a normal canvas refresh.
    """
    if background is not None:
        canvas.restore_region(background)
    for artist in artists:
        if artist is not None:
            axes.draw_artist(artist)
    if hasattr(canvas, 'blit'):
        canvas.blit(axes.bbox)
    else:
        canvas.draw_idle()
