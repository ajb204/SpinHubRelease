"""Compatibility alias for :mod:`decon.workflow.legacy_registry`."""
import sys as _sys
from spinDecon.workflow import legacy_registry as _impl
_sys.modules[__name__] = _impl
