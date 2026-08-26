"""Compatibility import for the migrated pseudo-2D diffusion workspace.

New code should import :mod:`decon.gui.workspaces.pseudo2d_diffusion`.
"""
from spinDecon.gui.workspaces.pseudo2d_diffusion import DiffusionROIFrame, Pseudo2DDiffusion

__all__ = ["DiffusionROIFrame", "Pseudo2DDiffusion"]
