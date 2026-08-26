# Refactor stages 54-57

## Regression gate

Stage 53 baseline in this environment: **314 passed, 4 failed, 1 skipped**.
The four failing identities are pre-existing source-contract/cwd-sensitive tests.

Final Stage 57 result: **317 passed, 4 failed, 1 skipped**. No new failing identity was accepted.
Whole-tree `compileall` and AST parsing also pass.

## Stage 54 - PeakFrame view boundary

Expanded `PeakService` with canonical peak-plane labels, bore payloads, cached views,
file-to-view adaptation and unit-converter access. PeakFrame now prefers this boundary
for construction-time labels, pseudo3D bore data and unit converters.

## Stage 55 - PeakFrame projection payload boundary

Moved 2D/pseudo3D/projection view selection into `PeakService.projection_payload()`.
PeakFrame retains its old implementation only as a standalone-construction fallback.
Added a service regression test for 2D raw/deconvolved plane selection.

## Stage 56 - Full Peak List authority

Added `FullPeakListService` and `ApplicationContext.full_peaks`. This service is the
explicit canonical boundary for the complete-spectrum peak collection. It deliberately
has no `conn_data`/connections API. Saving a full list rebuilds projected lists, notifies
analysis state, refreshes open Full Peak List viewers and refreshes status.

Added regression tests proving that Full Peak List persistence does not read or mutate
legacy `conn_data`.

## Stage 57 - Route consumers through Full Peak List authority

ProjectionService full-list payload/save operations and SliceService Full Peak List
opening now delegate to `FullPeakListService`. This prevents those services from growing
independent notions of full-list ownership.

## Architectural decision

`conn_data` remains a compatibility payload for legacy connectivity/NOE behaviour only.
It is not peak storage and must not become an alternative to the Full Peak List. Future
restored connection functions should be implemented behind a dedicated connection model
or service that references canonical peak identities from the Full Peak List.

## Next safe targets

1. Continue PeakFrame operational extraction (decon run, parameter-file/file-picker and alias/index operations).
2. Introduce a dedicated ConnectionService before restoring or migrating NOE/connectivity functionality.
3. Reduce remaining Slice2D/4D compatibility accesses without reintroducing local peak ownership.
4. Reconcile the remaining source-contract tests before physical relocation of source-sensitive modules.
