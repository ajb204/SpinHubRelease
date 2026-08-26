"""Compatibility alias for :mod:`decon.project.summary`.

The module object is aliased rather than wildcard-imported so historical
private helpers remain available to callers and regression tests.
"""
import sys as _sys
from .project import summary as _canonical
_sys.modules[__name__] = _canonical
