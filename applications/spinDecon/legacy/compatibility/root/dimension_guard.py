"""Compatibility alias for :mod:`decon.domain.dimensions.guard`."""
import sys as _sys
from spinDecon.domain.dimensions import guard as _impl
_sys.modules[__name__] = _impl
