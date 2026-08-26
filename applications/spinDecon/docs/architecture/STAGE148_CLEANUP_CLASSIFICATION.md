# Stage 148 - Cleanup and quarantine classification

## Decision boundary

The active project should contain only current implementation and current runtime entry points. Historical implementation is retained under `legacy/`, but active code must not import from `legacy/`. Historical compatibility paths are not considered part of the current API unless a specific release requirement says otherwise.

## Classification rules

- **ACTIVE** - used by the current application/workflow or a current subsystem.
- **ACTIVE SUPPORT** - current helper/dialog/reporting/processing code reached by ACTIVE code.
- **PACKAGE MARKER** - `__init__.py` belonging to an active package; retain.
- **CANONICAL COMPAT** - old import path forwarding to current implementation. Candidate to move to `legacy/compatibility/` and remove from active namespace.
- **LEGACY COMPAT** - old import path forwarding to `legacy/`. Remove from active namespace; the underlying legacy source is already retained.
- **QUARANTINE BRIDGE** - apparently current path whose only implementation is in `legacy/`. Move/remove from active namespace so legacy is genuinely quarantined.
- **QUARANTINED** - source under `legacy/`; retain.
- **REDUNDANT** - empty package/directory with no current role; remove.
- **GENERATED** - caches/OS metadata; remove.

## Workflow-driven active workspace set

The current notebook and workflow controller directly construct or route to:

- `gui/workspaces/nmr.py` - main NMR workspace.
- `gui/workspaces/projection.py` - projection inspection.
- `gui/workspaces/slice1d.py` and `slice2d.py` - slice review.
- `gui/workspaces/oned.py` - 1D view.
- `gui/workspaces/full3d.py` - true 3D workspace.
- `gui/workspaces/pseudo2d.py` - pseudo2D extraction/review.
- `gui/workspaces/pseudo2d_diffusion.py` - pseudo2D diffusion analysis.
- `gui/workspaces/pseudo3d.py` - fitting/pseudo3D path.
- `gui/workspaces/phasing.py` - phasing.
- `gui/workspaces/workflow.py` - workflow overview.

Current workspace support also includes:

- `full_peak_list.py`, `peak_review.py`, `peak_fit.py`, `peaks.py` - current peak-list/review/fitting infrastructure.
- `cpmg.py` and `decay.py` - launched from the active pseudo3D analysis selector; Decay is also used by project reporting. These are **ACTIVE**, despite their historical appearance.
- `slices.py` - current notebook import/support module.

### Workspace quarantine candidates

- `gui/workspaces/slice4d.py` - **QUARANTINE BRIDGE** only; it imports `legacy.slice4d.workspace`. The active notebook explicitly retires `AddTabFour4D`. Remove this active-path facade and retain `legacy/slice4d/`.

## Integrations

- `integrations/unidec/` - **ACTIVE**. `integrations/unidec/workspace.py` is used by current NMR/notebook integration.
- `integrations/magma/` - **QUARANTINE BRIDGE**. Both current-looking modules only forward to `legacy/magma`. Remove the active `integrations/magma/` package; retain `legacy/magma/`.
- `integrations/usta/` - **QUARANTINE BRIDGE**. `workspace.py` forwards to `legacy/usta/workspace_legacy.py`. Remove the active `integrations/usta/` package unless uSTA is deliberately restored to the workflow; retain `legacy/usta/`.

## Root modules

### Keep as current package/runtime

- `__init__.py` - package/version marker.

### CANONICAL COMPAT - remove from active root; preserve under `legacy/compatibility/root/` if historical source retention is desired

- `PeakShapeOptimizer.py` -> `analysis/peak_shape_optimizer.py`
- `analysis_mode.py` -> `domain/analysis_mode.py`
- `data_store.py` -> `project/data_store.py`
- `dataset_topology.py` -> `domain/topology.py`
- `decon_service.py` -> `project/decon_service.py`
- `decon_tab.py` -> `app/notebook.py`
- `dimension_guard.py` -> `domain/dimensions/guard.py`
- `dimension_labels.py` -> `domain/dimensions/labels.py`
- `parameter_store.py` -> `project/parameter_store.py`
- `pdfViewer.py` -> `gui/dialogs/pdf_viewer.py`
- `peak_dimension_contract.py` -> `domain/dimensions/peak_contract.py`
- `peak_picker.py` -> `analysis/peak_picker.py`
- `peak_shape_estimator.py` -> `analysis/peak_shape_estimator.py`
- `project_defaults.py` -> `project/defaults.py`
- `project_service.py` -> `project/service.py`
- `project_setup.py` -> `gui/dialogs/project_setup.py`
- `project_state.py` -> `project/state.py`
- `project_summary.py` -> `project/summary.py`
- `pseudo_axis_table.py` -> `domain/pseudo_axis.py`
- `shiftXPostFilter.py` -> `analysis/shiftx_post_filter.py`
- `viewer_dimension_contract.py` -> `domain/dimensions/viewer_contract.py`
- `workflow_model.py` -> `workflow/model.py`
- `workflow_overview.py` -> `workflow/overview.py`
- `workflow_registry.py` -> `workflow/legacy_registry.py`
- `workflow_status.py` -> `workflow/status.py`

No active non-test source imports these historical root paths. Their remaining dependency is predominantly regression tests that deliberately preserve compatibility.

### LEGACY COMPAT - remove from active root; implementation is already quarantined

- `LoadFilePopUp.py` -> `legacy/tools/LoadFilePopUp.py`
- `catia_tab.py` -> `legacy/catia/catia_tab.py`
- `multi.py`, `multi_tab.py`, `multiPlot.py`, `multiPlot2D.py`, `multiPlot2Dnorm.py`, `slicePlot2Dnorm.py` -> `legacy/manco/`
- `notepad.py` -> `legacy/tools/notepad.py`

These wrappers add no preservation value once old import compatibility is intentionally dropped.

## `Frames/`

`Frames/` contains no current implementation ownership. Active non-test source does not import it. Therefore the **whole namespace is removable from the active tree** once compatibility tests are updated.

### CANONICAL COMPAT

Move to `legacy/compatibility/Frames/` only if retaining historical import-adapter source is useful:

- `CPMGframe.py` -> `gui/workspaces/cpmg.py` (**underlying CPMG remains ACTIVE**)
- `DecayFrame.py` -> `gui/workspaces/decay.py` (**underlying Decay remains ACTIVE**)
- `Full3D.py` -> `gui/workspaces/full3d.py`
- `OneDView.py` -> `gui/workspaces/oned.py`
- `PhasingSpectra.py` -> `gui/workspaces/phasing.py`
- `Projection.py` -> `gui/workspaces/projection.py`
- `Pseudo2D.py` -> `gui/workspaces/pseudo2d.py`
- `Pseudo3D.py` -> `gui/workspaces/pseudo3d.py`
- `combineBrukerFrame.py`, `conversionFrame.py`, `processFrame.py`, `processProjectionsFrame.py`, `processingFrame.py` -> current `gui/dialogs/processing/`
- `deconFrame.py` -> `gui/workspaces/nmr.py`
- `matplotlib_toolbar.py` -> `gui/plotting/toolbar.py`
- `peakFitFrame.py`, `peakFrame.py`, `peakListFrame.py` -> current peak workspaces
- `pseudo2Ddiffusion.py` -> `gui/workspaces/pseudo2d_diffusion.py`
- `slicePlot.py`, `slicePlot2D.py` -> current slice workspaces
- `uindecNMRFrame.py` -> `integrations/unidec/workspace.py`
- `widgets.py` -> `gui/widgets/common.py`

### LEGACY COMPAT / QUARANTINE BRIDGES

Remove from active tree; legacy implementation remains:

- `NOEframe.py`
- `Pseudo2Dold.py`
- `SetDataStoreFrame.py`
- `catiaApp.py`
- `frameFeatures.py`
- `magmaFrame.py`, `magmaResults.py`
- `slicePlot4D.py`
- `STD_frame.py`, `uSTA/STD_frame.py`, `uSTA/uSTA_sims_frame.py`

`uSTA/uSTA_sims_frame.py` is especially redundant: it contains documentation only and no executable classes.

## `misc/`

All executable files in `misc/` are compatibility facades. No active non-test source imports `decon.misc`. Remove `misc/` from the active tree after moving the facades to `legacy/compatibility/misc/` if historical adapter source is worth retaining.

Mappings:

- `array_utils.py` -> `gui/plotting/array_utils.py`
- `display_utils.py` -> `gui/plotting/display_utils.py`
- `errors.py` -> `gui/dialogs/errors.py`
- `peak_io.py` -> `processing/peak_io.py`
- `peaks.py` -> `domain/peaks.py`
- `shell_output.py` -> `gui/dialogs/shell_output.py`
- `status_help.py` -> `gui/widgets/status_help.py`
- `textEdit.py` -> `gui/dialogs/text_viewer.py`

## `INDIANA/`

`INDIANA/cellDiff.py` is a compatibility facade for `legacy/indiana/cell_diff.py`. Remove the active `INDIANA/` package and retain `legacy/indiana/`.

## Empty/redundant structure

- `archive/` - empty: **REMOVE**.
- `tools/` - contains only zero-byte `__init__.py`: **REMOVE**.
- Empty `__init__.py` files inside active packages are **PACKAGE MARKERS**, not cleanup debris.

## Generated/non-source debris

Remove from the project/archive:

- `__MACOSX/`
- every `__pycache__/`
- every `*.pyc`
- `.pytest_cache/`

## Legacy

Everything under `legacy/` is **QUARANTINED / KEEP**. Do not delete it in Stage 148. The desired invariant is stronger than the current tree: **active packages should not import from `legacy/`**. The current violations are the Slice4D, MAGMA and uSTA bridge modules identified above.

## Tests that must change with compatibility removal

A number of tests intentionally enforce the old compatibility surfaces. These are architecture-history tests, not evidence that the current runtime needs those paths. In particular, Stage 68-147 architecture tests, canonical-module-identity tests, and several domain/project tests import historical root aliases.

The cleanup should update tests to import canonical modules and replace "facade must exist" assertions with the new invariant:

1. active source must not import `decon.legacy`;
2. active source must not import `decon.Frames`, `decon.misc`, or historical root aliases;
3. `legacy/` remains import-isolated from application startup;
4. workflow/notebook routes only to the current workspace set;
5. removed compatibility namespaces are absent from the active tree.

## Recommended execution order

### Pass 148A - no scientific behavior change

1. Remove generated debris, `archive/`, and empty `tools/`.
2. Move canonical compatibility facade source into `legacy/compatibility/{root,Frames,misc}/` for historical reference.
3. Remove legacy-facing facades from root/`Frames`/`INDIANA` because the actual implementations already exist in `legacy/`.
4. Remove `integrations/magma/`, `integrations/usta/`, and `gui/workspaces/slice4d.py`; they are bridges into quarantined code, not current implementations.
5. Keep `integrations/unidec/`.
6. Update active imports/tests/docs to canonical paths.
7. Run the full headless suite.

### Pass 148B - dead-current-code audit

After 148A is green, use import reachability plus workflow/runtime entry points to examine current-looking modules not reached by the workflow. Do **not** classify by filename or age alone: CPMG and Decay demonstrate why. For each candidate, inspect dynamic imports/callbacks and reporting/export use before moving it to `legacy/`.

### Pass 148C - native GUI golden journeys

Run the documented 1D, 2D, 3D, pseudo2D, pseudo3D, 4D and phasing journeys on wx/macOS. This is the point at which Stage 148 can be considered behavior-preserving.

## Proposed final top-level source layout

```text
decon/
    __init__.py
    analysis/
    app/
    domain/
    gui/
    integrations/
        unidec/
    legacy/
        compatibility/   # optional historical facade source, not import API
        ...preserved historical implementations...
    line_fitting/
    processing/
    project/
    workflow/
    tests/
    docs/
```

The key cleanup principle is: **being unused by the workflow is a strong signal, but an active dynamic caller or report/export path overrides that signal.** Conversely, a current-looking path whose only job is to import `legacy/` is not current code and should be moved behind the quarantine boundary.

## Pass 148A execution result (2026-08-23)

Pass 148A has now been applied.

- Root compatibility modules were moved to `legacy/compatibility/root/`.
- The complete historical `Frames/` compatibility namespace was moved to `legacy/compatibility/Frames/`.
- The complete `misc/` compatibility namespace was moved to `legacy/compatibility/misc/`.
- The top-level `INDIANA/` compatibility surface was moved to `legacy/compatibility/INDIANA/`; `legacy/indiana/` remains preserved.
- The current-looking MAGMA and uSTA bridges were moved to `legacy/compatibility/integrations/`; their implementations remain in `legacy/magma/` and `legacy/usta/`.
- `gui/workspaces/slice4d.py` was moved to `legacy/compatibility/gui/workspaces/`; the implementation remains in `legacy/slice4d/`.
- Empty `archive/` and `tools/` were removed.
- Generated Python/test caches were removed.
- Tests were changed to use canonical imports rather than historical root/misc compatibility paths.
- Compatibility-preservation assertions superseded by the final quarantine policy were retired and replaced by `tests/test_stage148_quarantine_boundary.py`.

The canonical workflow journeys retained in the active tree are 1D, pseudo2D, 2D, pseudo3D and 3D, together with their active supporting workspaces/services. uSTA, INDIANA, MAGMA and Slice4D are explicitly quarantined.

Validation after the move: **405 passed, 1 skipped** in the headless suite (`PYTHONPATH=.. python -m pytest -q tests`).

The next cleanup should be Pass 148B: audit current-looking modules for workflow/runtime reachability. No substantial module should move merely because its filename appears historical; dynamic workspace launchers, reporting/export paths and callback references must be checked first.

## Stage 148B execution

Stage 148B tightened the quarantine boundary while preserving the tested 1D,
pseudo2D, 2D, pseudo3D and 3D journeys.

Changes made:

- Renamed active `workflow/legacy_registry.py` to `workflow/registry.py`. The
  registry is actively imported by the notebook and project service and is
  therefore canonical despite its historical filename. A historical copy is
  retained under `legacy/compatibility/workflow/legacy_registry.py`.
- Removed active `workflow/overview.py`, which was only a compatibility alias
  for `gui/workspaces/workflow.py`; retained a historical copy under
  `legacy/compatibility/workflow/overview.py`.
- Removed the retired `AddMagmaTab` and `AddTabFour4D` compatibility hooks from
  the active notebook. MAGMA and Slice4D implementations remain quarantined.
- Extracted the unreachable `AssMan`, `AssManFrame`, and their private
  `SortedListCtrl` helper from active `gui/workspaces/slice2d.py`. They were an
  explicit legacy `conn_data`/4D-NOE viewer and are now retained at
  `legacy/compatibility/gui/slice2d_assignment_viewer.py`.
- Updated architecture tests to assert quarantine rather than preservation of
  these obsolete active surfaces.

### Current-workflow modules retained after audit

The workflow/notebook directly requires NMR, Projection, Slice1D, Slice2D,
OneD, Full3D, Pseudo2D, Pseudo2D Diffusion, Pseudo3D/Fitting, Phasing and the
Workflow overview. Full Peak List, Peak Review, Peak Fit, CPMG and Decay are
also retained because current journeys call them dynamically or reporting uses
them. `gui/workspaces/slices.py` and `gui/workspaces/peaks.py` are small but
active canonical aggregation boundaries, not compatibility residue.

### Deferred internal cleanup (do not quarantine yet)

Some canonical modules still contain *internal* migration fields and historical
formats (`legacy_nmr_workspace`, `tabOne`, `conn_data`, legacy parameter-file
parsing, and legacy vpar dimensional encoding). These are used by current
journeys, including CPMG/Decay and project loading, so they are not safe
file-level quarantine candidates. They should be reduced incrementally behind
current services rather than moved wholesale into `legacy/`.

The active Python source has no imports of `decon.legacy`; occurrences of
`legacy_vpar_dimension` refer to a current load/format migration helper in the
processing package, not a dependency on the quarantined package.

Validation after Stage 148B: **403 passed, 1 skipped**.

## Stage 148C - Slice2D conn_data quarantine

Stage 148C applies the final authority decision for the modern Slice2D journey:
the Full Peak List is authoritative and `conn_data` is not part of the current
Slice2D model or service boundary.

Changes:

- Preserved the pre-148C Slice2D implementation at
  `legacy/slice2d/conn_data_workspace_snapshot.py` so the historical NOE,
  connectivity, selection, and MAGMA behaviours remain available for later
  recovery behind a dedicated modern connection/NOE model.
- Removed the obsolete `RunFrame(..., conn_data, ...)` compatibility launcher
  from the active Slice2D workspace.
- Removed Slice2D conn_data overlays, conn_data click-selection/creation,
  NOE-test/fishing/addition routines, and the embedded MAGMA interaction path.
- Removed `connections`, `load_connections`, `analyse_connections`, and NOE-tag
  access from `SliceService`.  The current Slice2D service no longer exposes a
  connectivity API.
- Removed reconstruction of Slice2D 1D markers from `DataStore.conn_data` in
  `NMRWorkspace.get_reference_1d_view`.  Slice2D resolves peak markers from the
  authoritative Full Peak List instead.
- Updated boundary tests to require the absence of the retired SliceService
  connectivity API.

This pass deliberately does **not** delete the wider historical `conn_data`
payload yet.  It is still referenced by quarantined functionality and by some
older current-module compatibility/loading paths that require a separate audit.
It is no longer reachable through the modern Slice2D journey.

Validation after Stage 148C: **403 passed, 1 skipped**.

## Stage 148D - conn_data quarantine expansion

Stage 148D follows the architectural decision that the Full Peak List is the
authoritative peak collection for the modern 1D, pseudo2D, 2D, pseudo3D and 3D
journeys. `conn_data` is a deprecated connectivity/NOE representation and must
not be exposed as a parallel peak source by current workspaces.

Changes:

- Removed the `connections` / `conn_data` API from `OneDService` and
  `ProjectionService`.
- Removed the historical conn_data stick overlay from the current 1D viewer.
- Converted current Projection 2D peak-fitting input to the authoritative Full
  Peak List overlay.
- Removed remaining conn_data-based 4D projection marker hooks. Slice4D is
  already quarantined and these hooks had no role in the canonical journeys.
- Removed the dormant NOE manager source blocks from current CPMG and Decay
  modules. Complete pre-148D snapshots are retained under `legacy/conn_data/`.
- Added architecture regression tests preventing these GUI/service surfaces
  from reacquiring conn_data dependencies.

Remaining migration debt is intentionally confined primarily to `nmr.py`,
`project/data_store.py`, and `domain/peaks.py`. The NMR controller still parses
some deconvolution result files through `connEntry`/`conn_data` for current
2D/3D result-loading paths. This representation should not be deleted until a
canonical deconvolution-result record replaces it. The quarantine pass therefore
stops at that behavioural boundary rather than risking the tested journeys.

Validation after Stage 148D: **410 passed, 1 skipped**.

## Stage 148E - canonical nD peak-list model

The final active `conn_data` dependency has been removed from NMR result loading.
Calculated/deconvolved nD peak files are now loaded into the same dimension-
independent record schema used by the Full Peak List and stored as the distinct
`decon` peak-list role. Reference/projection, full, and deconvolved lists may
coexist; their role is semantic and does not imply a different storage type.

`DataStore.conn_data` and active `connEntry` have been removed. Pre-migration NMR,
domain-record, and DataStore implementations are preserved under
`legacy/conn_data/`. The governing architectural rule is documented in
`docs/architecture/PEAK_LIST_MODEL.md`.
