"""Compatibility alias for :mod:`decon.workflow.overview`."""
import sys as _sys
from spinDecon.workflow import overview as _impl
_sys.modules[__name__] = _impl
