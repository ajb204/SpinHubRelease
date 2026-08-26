"""Compatibility alias for :mod:`decon.project.defaults`."""
import sys as _sys
from spinDecon.project import defaults as _impl
_sys.modules[__name__] = _impl
