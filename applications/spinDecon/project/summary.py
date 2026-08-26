"""Compatibility alias for :mod:`decon.gui.reporting.project_summary`.

Project summary generation renders live GUI state and therefore belongs to the
GUI reporting layer. The module object is aliased so historical private helper
imports continue to work.
"""
import sys as _sys
from spinDecon.gui.reporting import project_summary as _canonical
_sys.modules[__name__] = _canonical
