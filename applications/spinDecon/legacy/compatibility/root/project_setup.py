"""Compatibility alias for :mod:`decon.gui.dialogs.project_setup`."""
import sys as _sys
from spinDecon.gui.dialogs import project_setup as _impl
_sys.modules[__name__] = _impl
