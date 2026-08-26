# Stages 82-85

## Stage 82 - Post-migration import hardening

Audited canonical modules for stale sibling-relative imports left behind by earlier
physical moves. Fixed three latent runtime paths that the source-only regression suite
could not exercise without wxPython: 3D peak-shape optimisation now imports its helper
from the decon analysis package, MAGMA PDF viewing imports the root compatibility viewer,
and uSTA peak fitting opens the canonical peak-fit workspace. Added regression guards
against reintroducing those stale relative imports.

## Stage 83 - Scientific helper ownership

Moved `PeakShapeOptimizer.py` to `analysis/peak_shape_optimizer.py` and
`shiftXPostFilter.py` to `analysis/shiftx_post_filter.py`. Historical root paths are now
compatibility imports. Canonical Peak Fit, Projection, and MAGMA code imports the analysis
owners directly.

## Stage 84 - Retired uSTA simulation quarantine

Moved the old `Frames/uSTA/uSTA_sims_frame.py` source to `legacy/usta/`. The preserved
file contains no executable definitions (the prototype had already been fully commented
out) and its application-shell entry points were already disabled. The historical Frames
path remains as a compatibility placeholder.

## Stage 85 - CATIA duplicate removal

Removed the duplicate active copies of the old CATIA application. The preserved
implementation now lives only under `legacy/catia/`; the root `catia_tab.py` and
`Frames/catiaApp.py` are compatibility entry points. Fixed the legacy application's own
relative import so the quarantine copy is internally coherent if it is deliberately run.

## Validation

- Full regression suite: 343 passed, 1 skipped, 0 failed.
- Whole-tree byte compilation passes.
- Whole-tree AST parsing passes.
- Full Peak List remains the authoritative peak collection.
- No new `conn_data` ownership was introduced.

## Next safe targets

1. Audit the remaining root-level standalone/prototype modules (`multi_tab.py`, `multi.py`,
   `notepad.py`, `LoadFilePopUp.py`) and quarantine only those confirmed disconnected.
2. Continue reducing root compatibility duplicates for domain/project/workflow modules,
   but only after migrating source-contract tests to canonical owners.
3. Audit optional integration launch paths under real wxPython because container tests do
   not import GUI modules.
