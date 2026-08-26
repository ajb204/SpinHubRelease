"""Compatibility alias for :mod:`decon.project.service`."""
import sys as _sys
from spinDecon.project import service as _impl
_sys.modules[__name__] = _impl
