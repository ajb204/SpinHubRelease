# Refactor stages 68-73

## Validation gate

Stage 67 baseline: **332 passed, 1 skipped, 0 failed**.
Final Stage 73 result: **335 passed, 1 skipped, 0 failed**. The three additional passes are new architecture regressions added in Stage 73. Whole-tree compile and AST parsing also pass.

## Stage 68 - project summary compatibility boundary

`decon/project/summary.py` is now the canonical implementation. The historical `decon/project_summary.py` is a module alias so private compatibility helpers remain available. Relative imports in the moved implementation were corrected to package-stable imports. This removes a duplicated 100k source file without changing the public API.

## Stage 69 - Projection physical migration

The active Projection implementation moved from `Frames/Projection.py` to `gui/workspaces/projection.py`. The historical path is now a small compatibility import. Source-inspection regressions were redirected to the canonical implementation and remained green.

## Stage 70 - Pseudo3D physical migration

The active Pseudo3D implementation moved to `gui/workspaces/pseudo3d.py`; `Frames/Pseudo3D.py` is now compatibility-only. Existing decon toggle and dimension regressions now inspect the canonical source.

## Stage 71 - slice workspace physical migration

The three slice implementations moved to `gui/workspaces/slice1d.py`, `slice2d.py`, and `slice4d.py`. `gui/workspaces/slices.py` is the application aggregator. Historical `Frames/slicePlot*.py` paths remain compatibility imports. The 4D NOE `AssMan` compatibility bridge remains available from the canonical Slice2D implementation; it has not been promoted to authoritative peak ownership.

## Stage 72 - peak workspace physical migration

Peak review, peak fitting, and the authoritative Full Peak List moved to `gui/workspaces/peak_review.py`, `peak_fit.py`, and `full_peak_list.py`. `gui/workspaces/peaks.py` is now the canonical aggregator. Historical Frames modules are compatibility-only. The Full Peak List remains authoritative; no `conn_data` API was introduced.

## Stage 73 - generic GUI infrastructure migration

The shared Matplotlib toolbar moved to `gui/plotting/toolbar.py` and reusable controls to `gui/widgets/common.py`. Canonical GUI modules no longer import these facilities through `Frames`. Compatibility wrappers remain for older imports.

## Next safe sequence

1. Move processing/conversion dialogs from `Frames` into `gui/dialogs/processing` while preserving wrappers.
2. Isolate MAGMA, uSTA and UniDec under `integrations/`; uSTA has the largest remaining direct `tabOne` hotspot and should receive a service boundary before a deeper rewrite.
3. Continue reducing residual `tabOne` fallbacks in canonical peak/slice/projection workspaces. Keep fallbacks only for intentionally supported standalone construction.
4. Delay `deconFrame` physical migration until optional integrations and processing dialogs no longer import it for helpers such as `peakEntry`/`ParseFlt`; extract those helpers first.
