# Refactoring stages 40-46

## Regression gate
Stage 39 was independently re-run before modification: **306 passed, 6 failed, 1 skipped**. The six failing identities are the same pre-existing ROI/source-contract failures recorded by earlier stages. Every stage below was accepted only with exactly those six failures and no new failures.

## Stage 40 - SliceService scientific-data API
Expanded `SliceService` with calibrated axis access, raw/deconvolved arrays, peak indices, connections, meshes, peak names, peak-shape widths, point sampling, projections, Full3D slice colour and full-peak-list hooks. Added a GUI-independent service regression test. Result: **307 passed, 6 failed, 1 skipped**.

## Stage 41 - Slice2D calibrated geometry boundary
Migrated coherent Slice2D geometry groups to `SliceService`: construction widths/names, 3D validation, calibrated axes, peak indices, intensity sampling, full-peak overlay colour and projection retrieval. Direct textual `tabOne` occurrences fell from 259 at Stage 39 to 194 at this checkpoint.

## Stage 42 - Slice4D construction/slicing boundary
Migrated peak-list choices, contour threshold initialisation and core raw 4D slicing data/mesh reads through `SliceService`.

## Stage 43 - PeakFit scientific array boundary
Expanded `PeakFitService` with raw fitting-array ownership plus shape/radius update hooks and project-save delegation. Migrated all operational `self.tabOne.data` reads in PeakFitFrame to the service. Direct textual `tabOne` occurrences fell from 95 at Stage 39 to 44.

## Stage 44 - Slice peak/data/connection ownership
Migrated direct Slice2D and Slice4D reads of raw/deconvolved arrays, peak lists, connection lists, labels and peak-index arrays to `SliceService`. This was intentionally done as property-equivalent substitutions only; behavioural callbacks were left untouched.

## Stage 45 - Slice calibrated-axis and decon state ownership
Migrated calibrated axis reads, deconvolution state, reference 1D retrieval, decon parameter path and full-peak-list opening through `SliceService`. Current direct textual `tabOne` occurrences: Slice2D **43**, Slice4D **20**. Remaining accesses are primarily legacy construction, NOE/MAGMA callbacks, status/connection commands and source-contract comments.

## Stage 46 - canonical GUI workspace imports
Added canonical `gui.workspaces.projection`, `gui.workspaces.pseudo3d`, and `gui.workspaces.slices` import boundaries and changed the application notebook to import through them. Implementations deliberately remain under `Frames/` until their legacy parent-chain callbacks are removed, preserving source-sensitive regression tests and standalone compatibility.

## Final validation
**307 passed, 6 failed, 1 skipped.** No new failing identity relative to Stage 39. Whole-tree `compileall` and AST parsing pass.

## Next safe sequence
1. Extract Slice NOE/connection/status operations into a dedicated interaction/service boundary rather than moving those callbacks mechanically.
2. Finish PeakFit control synchronisation through `PeakFitService` or a small presenter; remaining direct references are now mostly GUI state/persistence rather than scientific array reads.
3. Reduce Projection and Pseudo3D construction/status compatibility references, then move their implementations behind the already-created canonical workspace modules.
4. Once Slice2D/4D legacy callbacks are service-owned, physically relocate their implementations under `gui/workspaces` while leaving thin `Frames` compatibility modules for external callers.
