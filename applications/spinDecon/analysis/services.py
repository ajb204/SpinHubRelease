"""Composition root for GUI-facing analysis services.

The notebook should construct the NMR workspace, not know every scientific
service implementation.  This module centralises that temporary legacy bridge
until services can depend only on project/data/domain objects.
"""
from __future__ import annotations

from .diffusion_service import DiffusionService
from .full3d_service import Full3DService
from .full_peak_service import FullPeakListService
from .oned_service import OneDService
from .peak_fit_service import PeakFitService
from .peak_service import PeakService
from .phasing_service import PhasingService
from .projection_service import ProjectionService
from .pseudo_service import PseudoAxisService
from .slice_service import SliceService


def attach_analysis_services(context, legacy_workspace):
    """Attach the analysis service set to an ``ApplicationContext``.

    ``legacy_workspace`` is intentionally explicit: it is the migration input,
    not an ownership API.  Individual services can be modernised independently
    without expanding the notebook's dependency surface.
    """
    context.nmr_workspace = legacy_workspace
    context.legacy_nmr_workspace = legacy_workspace  # compatibility for external/older hosts
    context.full3d = Full3DService(legacy_workspace)
    context.one_d = OneDService(legacy_workspace)
    context.projection = ProjectionService(legacy_workspace)
    context.peaks = PeakService(legacy_workspace)
    context.full_peaks = FullPeakListService(legacy_workspace)
    context.peak_fit = PeakFitService(legacy_workspace)
    context.phasing = PhasingService(legacy_workspace)
    context.slices = SliceService(legacy_workspace)
    context.pseudo = PseudoAxisService(legacy_workspace)
    context.diffusion = DiffusionService(legacy_workspace)
    return context
