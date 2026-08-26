"""Legacy compatibility import for the retired NOE workspace.

Connectivity functionality may return later behind a dedicated connection model;
this module must not become a second peak authority.
"""
from spinDecon.legacy.noe.workspace import *  # noqa: F401,F403
