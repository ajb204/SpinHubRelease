"""Shared application dependencies passed to GUI workspaces.

`legacy_nmr_workspace` is an explicit migration bridge. New code should use
project/data/parameters/services instead of reaching through the GUI tree.
"""
from dataclasses import dataclass
from typing import Any

@dataclass
class ApplicationContext:
    project: Any = None
    data: Any = None
    parameters: Any = None
    decon: Any = None
    workflow: Any = None
    full3d: Any = None
    one_d: Any = None
    projection: Any = None
    peaks: Any = None
    full_peaks: Any = None
    peak_fit: Any = None
    phasing: Any = None
    slices: Any = None
    pseudo: Any = None
    diffusion: Any = None
    nmr_workspace: Any = None
    legacy_nmr_workspace: Any = None  # deprecated compatibility alias
