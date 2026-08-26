# Stages 52-53: Full Peak List authority and conn_data quarantine

## User-confirmed architectural rule

The Full Peak List is the authoritative complete peak list. `conn_data` is legacy.
Some connectivity/NOE functions may be restored later, but restoration must use a
dedicated connection model/service linked to the Full Peak List; it must not make
`conn_data` or Slice2D-local peak state authoritative again.

## Stage 52 - remove Slice2D-local peak manager

- Removed the deprecated modeless `peaks_box` and its Select/Deselect/Delete/Add controls.
- The Slice2D Peaks toolbar button now opens the canonical Full Peak List through `SliceService`.
- Removed the obsolete Slice2D `on_NOE_button` route into the old conn_data viewer.
- Retained `AssMan/AssManFrame` only as an explicitly labelled 4D NOE compatibility bridge.
- Preserved deferred NOE semantics as a future-migration marker beside `peakListFrame`, rather than as an active Slice2D API.

## Stage 53 - remove orphaned conn_data UI callbacks

- Removed the orphaned local connectivity search/save/load callbacks left behind after the peak window was removed.
- Removed the no-op `update_conn_data` shim and stale NOE-button enablement.
- Documented `SliceService.connections` as compatibility-only legacy connectivity data, never peak storage.
- Added `test_full_peak_authority_boundary.py` to protect the new ownership rule.
- Updated the old source-contract test so it expects the service boundary rather than a direct `tabOne` callback.

## Regression gate

Final suite: **314 passed, 4 failed, 1 skipped**.

The four remaining failures are pre-existing/source-contract issues unrelated to this change:
1. fit-radius test assumes a historical working-directory path;
2-4. pseudo2D diffusion source-text expectations predate its service/data migration.

No new failing test was introduced. Two previously failing Slice2D cleanup tests now pass because the intended cleanup has actually been completed.

Whole-tree `compileall` and AST parsing pass.
