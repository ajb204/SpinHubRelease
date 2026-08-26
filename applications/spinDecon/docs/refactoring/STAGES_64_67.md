# Refactor stages 64-67

## Stage 64 - Canonical peak workspace/report boundary

Added `gui/workspaces/peaks.py` as the application-facing import surface for
peak review, peak-shape fitting and the authoritative Full Peak List UI.

Corrected `project/summary.py`, which had retained invalid package-relative
`.Frames` imports after the earlier project-package migration. Report generation
now imports peak workspaces through the canonical GUI boundary.

Regression gate: **328 passed, 1 skipped**.

## Stage 65 - Quarantine duplicate top-level STD implementation

Audited `Frames/STD_frame.py` against the maintained
`Frames/uSTA/STD_frame.py`. The notebook imports the uSTA implementation and no
active internal code imports the old top-level implementation. The old source
has therefore been retained under `legacy/usta/std_frame_legacy.py`, while
`Frames/STD_frame.py` is now a compatibility import to the maintained class.

This is quarantine, not deletion: historical functionality remains available
for deliberate future recovery.

Regression gate: **330 passed, 1 skipped**.

## Stage 66 - Physical migration of the 1D workspace

Moved the active `OneDFrame` implementation from `Frames/OneDView.py` to
`gui/workspaces/oned.py`. The historical path is now a compatibility import.
Updated the source-contract regression to inspect the canonical implementation.

Regression gate: **331 passed, 1 skipped**.

## Stage 67 - Physical migration of the Pseudo2D workspace

Moved the active `Pseudo2D` and `Pseudo2DFittingFrame` implementation from
`Frames/Pseudo2D.py` to `gui/workspaces/pseudo2d.py`. The old path is now a
compatibility import. Source-sensitive Pseudo2D regressions were redirected to
the canonical implementation rather than weakened or removed.

Regression gate: **332 passed, 1 skipped, 0 failed**.

## Architectural state

The normal application now has physical implementations under
`gui/workspaces` for 1D, Pseudo2D, Pseudo2D diffusion, Full3D and phasing, with
canonical application import surfaces for the remaining migrated workspaces.

The Full Peak List remains authoritative. `conn_data` remains legacy
compatibility state only and was not added to the canonical peak workspace API.

## Next safe sequence

1. Continue reducing PeakFrame's remaining legacy construction/data fallbacks,
   then physically migrate peak/peak-fit implementations when their
   source-contract tests can follow the canonical paths.
2. Establish an explicit uSTA service/context boundary before moving the active
   `Frames/uSTA` implementation; it remains one of the largest `tabOne`
   coupling hotspots.
3. Continue reducing Slice2D/Slice4D operational callbacks behind SliceService.
4. Audit CATIA and `multi_tab.py` as parallel/standalone application shells and
   quarantine them only after confirming they are not supported entry points.
5. Preserve legacy NOE/connectivity code for future restoration behind a
   dedicated connection model; do not restore `conn_data` as peak authority.
