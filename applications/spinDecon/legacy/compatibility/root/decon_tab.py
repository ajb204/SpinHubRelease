"""Compatibility entry point for :mod:`decon.app.notebook`.

The application shell is canonical under ``decon.app``; this historical module
remains for launch scripts importing ``decon.decon_tab.MyApp``.
"""
from .app.notebook import *  # noqa: F401,F403
