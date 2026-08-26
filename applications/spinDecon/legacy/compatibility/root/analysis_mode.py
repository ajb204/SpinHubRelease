"""Compatibility alias for :mod:`decon.domain.analysis_mode`."""
import sys as _sys
from spinDecon.domain import analysis_mode as _impl
_sys.modules[__name__] = _impl
