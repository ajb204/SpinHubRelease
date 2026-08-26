# Refactor stages 58-61

## Stage 58 - diffusion application boundary

Added `analysis.diffusion_service.DiffusionService` and wired it into
`ApplicationContext`.  The pseudo-2D diffusion workspace now obtains spectrum
arrays, calibrated axes, labels, thresholds, noise and parameter-file access
through this boundary.  It no longer uses `self.tabOne` as its scientific data
API.

## Stage 59 - regression suite repair

The four remaining failures were source-contract tests rather than behavioural
regressions.  One test used a working-directory-relative path; three diffusion
tests asserted implementation details that had already been superseded by the
normalised ROI implementation/service migration.  The tests were updated to
assert the current contracts.  This establishes a clean regression gate:
`323 passed, 1 skipped` at this stage.

## Stage 60 - PeakService operational boundary

Expanded `PeakService` with parameter access, deconvolution dispatch, external
2D view caching, project parameter persistence, status refresh, aliasing and
physical-axis access.  PeakFrame now uses those operations instead of reaching
into deconFrame.  Direct textual `tabOne` occurrences in `Frames/peakFrame.py`
fall from 141 at Stage 57 to 54; remaining occurrences are concentrated in
legacy fallback paths and specialist operations rather than basic data access.

## Stage 61 - canonical diffusion workspace

Moved the active pseudo-2D diffusion implementation to
`gui/workspaces/pseudo2d_diffusion.py`.  `Frames/pseudo2Ddiffusion.py` is now a
small compatibility import.  Application and Pseudo2D callers use the canonical
workspace path, and source-inspection tests follow the canonical implementation.

## Architectural status

- Full Peak List remains the authoritative complete peak collection.
- `conn_data` remains legacy compatibility state and is not exposed by the
  Full Peak List service.
- New scientific/service modules remain wx-independent.
- New GUI workspace implementations are progressively leaving `Frames/` only
  after their state/data dependencies have been extracted.

## Next safe targets

1. Finish the remaining PeakFrame specialist-operation boundary.
2. Extract Pseudo2D's remaining direct project/data operations before moving
   its active implementation out of `Frames/`.
3. Revisit PeakFitFrame's residual compatibility/control accesses.
4. Move Slice2D/4D implementations only after source-sensitive contracts have
   been redirected to canonical paths.
