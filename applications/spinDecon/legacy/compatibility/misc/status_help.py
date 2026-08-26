"""Compatibility alias for :mod:`decon.gui.widgets.status_help`."""
import sys as _sys
from spinDecon.gui.widgets import status_help as _impl
_sys.modules[__name__] = _impl
