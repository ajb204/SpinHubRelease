# Refactoring Stage 51

## Regression gate
Stage 50 baseline, invoked from the package directory with the project parent on `PYTHONPATH`: **309 passed, 6 failed, 1 skipped**. The six failures are the established pseudo2D diffusion/source-cleanup failures. Stage 51 preserves exactly those six failing identities.

## Stage 51 - Slice2D calibrated-mesh boundary
Expanded `SliceService` so calibrated 2D strip meshes are obtained through the application service rather than directly from `parent.tabOne.XX/YY/ZZ`. `Slice2D.ReSlice2d` now asks `SliceService.slice_meshes(...)` for both orientations. The service also exposes a temporary `datadec` compatibility alias because the migrating viewer still uses that historic spelling in a few locations.

Moved the remaining active 4D spectral-dimension checks in the Slice2D connection-list UI to `SliceService.spectral_dimension`, and routed the Full Peak List toolbar action through `SliceService.open_full_peak_list` rather than directly discovering `tabOne.OnButtonFullPeakList`.

Added a regression test covering both strip-mesh orientations and the deconvolution-data compatibility alias.

## Validation
**310 passed, 6 failed, 1 skipped.** No new failing identity relative to Stage 50. Whole-tree `compileall` and AST parsing pass.

## Next safe boundary
The remaining Slice2D `tabOne` occurrences are now predominantly constructor compatibility/fallbacks, comments, and the deprecated local NOE/connection manager. The next architectural step should reconcile the source-sensitive Slice2D cleanup tests and then remove/quarantine that deprecated local manager, rather than continuing mechanical substitutions.
