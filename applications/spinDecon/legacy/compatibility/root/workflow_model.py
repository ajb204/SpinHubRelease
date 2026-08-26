"""Compatibility alias for :mod:`decon.workflow.model`."""
import sys as _sys
from spinDecon.workflow import model as _impl
_sys.modules[__name__] = _impl
