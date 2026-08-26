"""Compatibility alias for :mod:`decon.domain.topology`."""
import sys as _sys
from spinDecon.domain import topology as _impl
_sys.modules[__name__] = _impl
