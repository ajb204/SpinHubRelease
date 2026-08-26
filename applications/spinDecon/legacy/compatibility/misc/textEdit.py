"""Compatibility alias for :mod:`decon.gui.dialogs.text_viewer`."""
import sys as _sys
from spinDecon.gui.dialogs import text_viewer as _impl
_sys.modules[__name__] = _impl
