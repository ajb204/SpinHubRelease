"""Compatibility alias for :mod:`decon.domain.pseudo_axis`."""
import sys as _sys
from spinDecon.domain import pseudo_axis as _impl
_sys.modules[__name__] = _impl
