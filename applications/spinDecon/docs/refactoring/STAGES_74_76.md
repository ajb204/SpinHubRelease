# Stages 74-76

## Stage 74 - Processing GUI canonical boundary

Created `gui/dialogs/processing/` with canonical copies of the processing, conversion,
projection-processing, and Bruker-combiner GUI implementations. Internal application
imports now target the canonical package. The historical `Frames/` sources are retained
in full temporarily because a large set of regression tests deliberately inspect their
source text; replacing them with wrappers caused 43 source-contract failures and was
reverted. Runtime behavior remained green after the import migration.

## Stage 75 - Optional integration isolation

Created `integrations/magma/`, `integrations/unidec/`, and `integrations/usta/`.
MAGMA and uSTA are now imported by the application from these canonical locations.
Historical `Frames` paths are compatibility wrappers for MAGMA, UniDec, and the maintained
uSTA implementation. MAGMA's results import was made package-explicit.

## Stage 76 - CATIA quarantine checkpoint

Preserved the old parallel CATIA application sources under `legacy/catia/`. The original
paths remain untouched for compatibility because CATIA is not part of the active notebook
workflow and should not be deleted until external launch/use is ruled out.

## Regression gate

- pytest: 335 passed, 1 skipped, 0 failed
- Full Peak List remains the authoritative peak list.
- No new `conn_data` ownership was introduced.
- Processing source-sensitive tests remain intact rather than being weakened.

## Next safe sequence

1. Move processing source-contract tests to canonical paths in a dedicated test migration,
   then replace `Frames` processing files with compatibility wrappers.
2. Audit remaining active imports of `Frames/deconFrame.py`; extract helper functions that
   are still imported by canonical modules.
3. Introduce the final `NMRWorkspace` canonical module only after those consumers are gone.
4. Audit dormant `NOEframe`, `CPMGframe`, `DecayFrame`, `SetDataStoreFrame`, and uSTA sims
   for active callbacks before quarantine/deletion decisions.
