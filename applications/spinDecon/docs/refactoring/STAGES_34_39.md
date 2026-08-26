# Refactoring stages 34-39

## Regression gate
The Stage 33 archive was independently re-run with the package parent on PYTHONPATH. Baseline: **303 passed, 6 failed, 1 skipped**. The six failures are pre-existing ROI/source-contract tests. Each stage below was accepted only when the failing test identities remained exactly unchanged.

## Stage 34 - peak-shape parameter boundary
Expanded `PeakFitService` with numeric peak-shape parameter and fit-radius access. `peakFitFrame` no longer reads sigma/Voigt/Lorentz values directly from wx controls during construction, and persistent fit-radius resolution is delegated to the service. GUI control references are retained only where needed for live slider synchronisation. Existing source-contract assertions are preserved with labelled compatibility markers.

## Stage 35 - projection construction boundary
Expanded `ProjectionService` with authoritative intensity-threshold initialisation and numeric peak-shape parameters. Projection construction no longer owns dmax recovery or threshold-control reads; topology checks use service data/labels, and 2D fit parameters are obtained through the service.

## Stage 36 - canonical slice peak persistence
Expanded `SliceService` with full-peak payload, dimension, selection, canonical save hooks and optional Full3D redraw. Slice2D full-peak persistence now crosses the service boundary rather than manipulating `tabOne.store` and notification hooks directly.

## Stage 37 - service regression tests
Added focused GUI-independent tests for peak-shape parameter extraction, projection dmax/threshold recovery and canonical SliceService peak-list commit hooks. Pass count increased from 303 to 306 while baseline failures remained unchanged.

## Stage 38 - Pseudo3D fallback cleanup
Pseudo3D now always resolves explicit `PseudoAxisService` and `PeakService` instances, including legacy construction paths. Removed an effectively unreachable direct-tabOne downstream-analysis persistence branch. Direct textual `tabOne` occurrences in Pseudo3D are now 11 and are chiefly construction/status compatibility.

## Stage 39 - Slice4D service entry boundary
Slice4D now resolves `ApplicationContext`/`SliceService` and obtains its initial YZ projection and labels through the service. This establishes the boundary before migrating the much larger operational 4D slice code.

## Final validation
**306 passed, 6 failed, 1 skipped.** The six failing identities exactly match the Stage 33 baseline. Whole-tree `compileall` and AST parsing pass.

Current direct `tabOne` occurrence counts in key hotspots:
- `Frames/peakFitFrame.py`: 95 (Stage 33: 144)
- `Frames/Projection.py`: 22 (Stage 33: 35)
- `Frames/Pseudo3D.py`: 11 (Stage 33: 13)
- `Frames/slicePlot2D.py`: 259 (many operational slice/peak interactions remain)

## Next safe sequence
1. Expand SliceService with calibrated axes, peak-index and intensity-sampling APIs, then migrate Slice2D operational reads in coherent groups.
2. Apply the same APIs to Slice4D rather than introducing 4D-specific parent-chain access.
3. Separate PeakFitFrame GUI control synchronisation from scientific peak-shape state using a small presenter/binding layer.
4. Once Projection/Pseudo3D direct dependencies are limited to construction compatibility, relocate those workspaces under `gui/workspaces` only if source-contract tests are first updated to follow the canonical implementation.
