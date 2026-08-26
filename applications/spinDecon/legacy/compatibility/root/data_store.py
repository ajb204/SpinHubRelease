"""Compatibility alias for :mod:`decon.project.data_store`."""
import sys as _sys
from spinDecon.project import data_store as _impl
_sys.modules[__name__] = _impl
