"""Compatibility alias for :mod:`decon.gui.dialogs.shell_output`."""
import sys as _sys
from spinDecon.gui.dialogs import shell_output as _impl
_sys.modules[__name__] = _impl
