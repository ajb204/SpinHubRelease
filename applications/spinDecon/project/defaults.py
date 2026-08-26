"""Canonical defaults for project settings that need not be persisted.

A parameter omitted from the system file means "use the application default".
Only user overrides should be written for these settings.
"""
from __future__ import annotations
import os

UNIDEC_DEFAULTS = {
    "thresh": 0.08,
    "fac": 1.4,
    "conv": 1e-7,
    "maxiter": 100,
}


def available_cpu_count() -> int:
    """Return the maximum number of CPUs currently available to this process."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def is_default_value(key: str, value) -> bool:
    """Compare a GUI/string value with its canonical default numerically."""
    if key == "ncpus":
        default = available_cpu_count()
    else:
        default = UNIDEC_DEFAULTS[key]
    try:
        return float(value) == float(default)
    except (TypeError, ValueError):
        return str(value).strip() == str(default)
