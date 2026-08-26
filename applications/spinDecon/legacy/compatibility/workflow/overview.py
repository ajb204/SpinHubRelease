"""Compatibility alias for :mod:`decon.gui.workspaces.workflow`.

The workflow model and status engine remain in :mod:`decon.workflow`; the wx
workflow overview is a GUI workspace.
"""
import sys as _sys
from spinDecon.gui.workspaces import workflow as _canonical
_sys.modules[__name__] = _canonical
