"""Compatibility alias for :mod:`decon.project.state`."""
import sys as _sys
from spinDecon.project import state as _impl
_sys.modules[__name__] = _impl
