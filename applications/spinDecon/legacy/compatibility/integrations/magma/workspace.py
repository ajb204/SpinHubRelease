"""Compatibility facade for the retired MAGMA workspace.

MAGMA is not part of the active application architecture.  The implementation
is preserved under :mod:`decon.legacy.magma.workspace` for future recovery.
"""
from spinDecon.legacy.magma.workspace import *  # noqa: F401,F403
