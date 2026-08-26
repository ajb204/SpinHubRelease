"""Compatibility alias for :mod:`decon.project.parameter_store`."""
import sys as _sys
from spinDecon.project import parameter_store as _impl
_sys.modules[__name__] = _impl
