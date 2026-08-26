"""Compatibility facade for plotting display helpers."""
import sys as _sys
from spinDecon.gui.plotting import display_utils as _impl
_sys.modules[__name__] = _impl
