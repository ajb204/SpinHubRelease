"""Compatibility alias for :mod:`decon.domain.dimensions.labels`."""
import sys as _sys
from spinDecon.domain.dimensions import labels as _impl
_sys.modules[__name__] = _impl
