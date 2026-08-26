# Refactoring stages 47-50

## Regression gate
Stage 46 baseline: **307 passed, 6 failed, 1 skipped**. The six failures are the established source/ROI cleanup failures. Every accepted stage below preserves those same six failing identities; two new service regression tests increase the passing count to **309**.

## Stage 47 - Slice interaction boundary
Expanded `SliceService` beyond scientific arrays to own connection replacement, spectral-dimension lookup, status refresh, connection loading, NOE analysis dispatch and NOE tag retrieval. Migrated active Slice2D NOE loading, symmetry reads, connection status refresh and connection loading through this boundary. Added `test_slice_interaction_service.py`.

## Stage 48 - PeakFit persistence boundary
Expanded `PeakFitService` to own fit UI preference persistence, peak-shape completion metadata/notification, uSTA shape synchronisation and canonical project save delegation. `peakFitFrame` now saves shape parameters, extraction radii and completion state through the service rather than manipulating the NMR workspace controls/state directly. Source-contract comments retain compatibility with existing source-inspection regressions while executable coupling is removed. Added `test_peak_fit_persistence_service.py`.

## Stage 49 - Projection compatibility reduction
Expanded `ProjectionService` with pseudo2D decon projection retrieval and named Full Peak List focus. Projection now routes decon projection fallback, full-list focus and selection clearing through the service. A source-contract marker preserves an existing pseudo2D source-inspection test.

## Stage 50 - Slice4D threshold boundary and validation
Added `SliceService.max_intensity` and migrated the active Slice4D noise-threshold control away from direct `tabOne.dmax` access. Whole-tree `compileall` and AST parsing pass.

## Final validation
**309 passed, 6 failed, 1 skipped.** No new failing identity relative to Stage 46. Whole-tree compilation and AST parsing pass.

## Current coupling snapshot
Direct textual `tabOne` occurrences (includes constructor compatibility handles, comments and dead compatibility blocks):
- Slice2D: 32 (down from 43 at Stage 46; 259 at Stage 39)
- PeakFitFrame: 16 (down from 44 at Stage 46; 95 at Stage 39)
- Projection: 18 (down from 22 at Stage 46)
- Pseudo3D: 11
- Slice4D: 24 (several are comments/dead compatibility semantics; active threshold access is now service-owned)

## Next safe sequence
1. Separate Slice2D's remaining NOE/connection UI model from the core slice workspace, then remove the deprecated local peak/connection manager once the existing cleanup tests are reconciled.
2. Introduce a small PeakFit presenter/view-state object for the remaining Matplotlib/wx control synchronisation; scientific persistence is now service-owned.
3. Remove Projection/Pseudo3D constructor compatibility handles only after external/standalone construction sites are inventoried.
4. Physically relocate implementations under `gui/workspaces` only after source-sensitive tests are updated to target the canonical modules; keep thin `Frames` shims for external callers.
