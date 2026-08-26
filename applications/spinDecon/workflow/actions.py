"""Compatibility alias for the application workflow controller.

Workflow models/status live in :mod:`decon.workflow`; GUI/application action
routing belongs to :mod:`decon.app.workflow_controller`.
"""
import sys as _sys
from spinDecon.app import workflow_controller as _canonical
_sys.modules[__name__] = _canonical
