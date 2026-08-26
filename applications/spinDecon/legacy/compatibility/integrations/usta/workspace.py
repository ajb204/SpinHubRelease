"""Compatibility import for the retired uSTA workspace.

The historical uSTA GUI is not connected to the active application workflow.
Its implementation is preserved under :mod:`decon.legacy.usta` so useful
scientific routines can be recovered later without treating the old tab as an
active integration.
"""
from spinDecon.legacy.usta.workspace_legacy import STDFrame

__all__ = ["STDFrame"]
