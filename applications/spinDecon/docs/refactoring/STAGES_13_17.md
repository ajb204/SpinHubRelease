# Refactor stages 13-17

All stages were regression-gated against the recorded baseline. The baseline
remains **297 passed, 7 failed, 1 skipped**; no new failing test was introduced.

## Stage 13 - application-context bridge for legacy Frames

Added `gui/context.py` and attached `ApplicationContext` to the legacy NMR
workspace. Migrated initial viewer state/data ownership lookups away from
parent-chain inference in Full3D, OneDView, slicePlot, Pseudo2D, Pseudo3D and
PhasingSpectra. Scientific callbacks remain on the legacy workspace until a
service boundary exists for each family.

## Stage 14 - Full3D application service

Added `analysis/full3d_service.py`. Full3D presentation code now calls a
context-owned service for view specifications, slices, cross-sections,
intensity limits, overlays and peak-selection clearing. The service currently
delegates to the legacy workspace, providing a controlled seam for later
extraction of the numerical implementation.

## Stage 15 - first workspace migration

Moved the active Full3D implementation to `gui/workspaces/full3d.py`.
`Frames/Full3D.py` is now a compatibility import. The application shell imports
the new location directly.

## Stage 16 - phasing workspace migration

Moved the active phasing implementation to `gui/workspaces/phasing.py` with a
compatibility import at `Frames/PhasingSpectra.py`. The application shell uses
the new location.

## Stage 17 - legacy quarantine begins

Moved the superseded Pseudo2D implementation to `legacy/pseudo2d_old.py`, with
an old-path compatibility import. Other probable legacy shells are documented
but deliberately not moved until standalone/external use can be ruled out.

## Next high-value work

1. Introduce equivalent service boundaries for peak lists and projection/slice
   view data before moving those high-coupling Frames.
2. Migrate `Projection`, `OneDView`, then slice viewers while preserving source-
   inspection tests at compatibility paths where required.
3. Move peak GUI only after peak repository/service ownership is explicit.
4. Leave Pseudo3D and final `deconFrame` -> `NMRWorkspace` conversion until the
   above dependencies have been removed.
