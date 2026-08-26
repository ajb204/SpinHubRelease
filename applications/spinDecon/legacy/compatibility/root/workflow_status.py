"""Compatibility alias for :mod:`decon.workflow.status`."""
import sys as _sys
from spinDecon.workflow import status as _impl
_sys.modules[__name__] = _impl
