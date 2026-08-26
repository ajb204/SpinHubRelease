# Refactor stages 25-29

## Regression gate
Baseline for this supplied Stage 24 archive in the current environment: **301 passed, 6 failed, 1 skipped**. The six failures are pre-existing source-contract regressions in `pseudo2Ddiffusion.py` and `slicePlot2D.py`. Every stage below was required to preserve that exact failure set.

## Stage 25 - pseudo-axis application boundary
- Added `analysis/pseudo_service.py` (`PseudoAxisService`).
- Added `pseudo` to `ApplicationContext` and initialise it beside the other migration services.
- Pseudo2D now obtains its processed array, spectral axis and display labels through this boundary when an application context is available.

## Stage 26 - peak persistence boundary
- Expanded `PeakService` with canonical projection-peak persistence.
- `PeakFrame._commit_projection_peaks()` now delegates persistence/controller synchronisation to `PeakService` rather than directly orchestrating `deconFrame` internals.
- Legacy fallback remains for standalone/older construction paths.

## Stage 27 - Slice2D construction boundary
- Slice2D now resolves `ApplicationContext` explicitly.
- Initial reference peaks, threshold, symmetry setting and spectrum path are obtained through `SliceService` where available.
- Extended `SliceService` with symmetry and spectrum-path accessors.

## Stage 28 - Pseudo3D shared-state boundary
- Pseudo3D now resolves both `PseudoAxisService` and `PeakService`.
- Initial threshold uses `PeakService`.
- Analysis-change notification and downstream-analysis persistence use `PseudoAxisService`.

## Stage 29 - Pseudo2D persistence boundary and validation
- Pseudo2D downstream-analysis notification/persistence now uses `PseudoAxisService`.
- Architecture regression now requires `ApplicationContext.pseudo`.
- Full tree passes `compileall` and AST parsing.

## Validation
- pytest: **301 passed, 6 failed, 1 skipped** after each stage; no new failures.
- `python -m compileall -q .`: pass.
- AST parse of all Python source: pass.
- A live `decon.decon_tab` import cannot be executed in the refactoring container because wxPython is not installed. The source/compile checks are therefore retained in addition to pytest. The previously reported malformed Full3D import remains fixed.

## Next safe sequence
1. Expand `SliceService` around full-peak selection/persistence before editing the dense Slice2D/4D callback paths.
2. Extract pseudo group/FUDA path operations behind `PseudoAxisService` before further Pseudo3D migration.
3. Extract peak projection/view payload fallbacks still remaining in PeakFrame.
4. Only after those boundaries stabilise, migrate active GUI modules from `Frames/` to `gui/workspaces/` while preserving source-sensitive compatibility files until their regression tests are modernised.
