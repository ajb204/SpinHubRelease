# Refactor stages 22-24

## Stage 22 - launch regression repair and import hygiene

Repaired a migration error in `gui/workspaces/full3d.py`.  During the Full3D
move, the original standard-library/scientific imports
`string, copy, math, numpy, os` were accidentally appended to the
`decon.gui.context` import.  This caused `spinDecon` to fail during import.
The imports are now restored to their original ownership.

Added `test_gui_context_import_hygiene.py`, which parses every Python module
and ensures that imports from `decon.gui.context` are restricted to the three
migration helpers (`context_for`, `project_for`, `data_for`).  This directly
regresses the launch failure without requiring wxPython in the headless test
environment.

A complete `compileall`/AST parse of the package succeeds.

## Stage 23 - peak application boundary

Added `analysis/peak_service.py` and exposed it as `ApplicationContext.peaks`.
PeakFrame now obtains its initial/synchronised threshold and spectrum path via
the service where an application context is available.  Legacy fallbacks are
retained deliberately.

## Stage 24 - slice application boundary

Added `analysis/slice_service.py` and exposed it as `ApplicationContext.slices`.
The 1D slice viewer now obtains reference peaks, threshold and reference 1D
view payloads through this boundary, while preserving legacy fallback paths.

Added architecture regression tests requiring analysis service modules to
remain independent of wx and requiring the expected migration services on
`ApplicationContext`.

## Regression status

The suite remains at the pre-existing failure baseline.  New architecture
checks pass.  No new functional regression was accepted during these stages.

## Next migration order

1. Expand PeakService around canonical peak-list persistence and selection.
2. Expand SliceService for 2D/4D payloads before modifying the source-sensitive
   `slicePlot2D.py` tests.
3. Introduce a pseudo-axis service for Pseudo2D/Pseudo3D downstream-analysis
   state and data payloads.
4. Only then migrate the corresponding Frames into `gui/workspaces`.
