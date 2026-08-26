# Refactor Stages 86-87

## Stage 86 - quarantine disconnected MANCO-era and standalone GUI prototypes

Static dependency inspection confirmed that `multi_tab.py`, `multi.py`, `multiPlot.py`, `multiPlot2D.py`, `multiPlot2Dnorm.py`, and `slicePlot2Dnorm.py` form a self-contained historical MANCO-era application and are not imported by the active spinDecon application. Their implementations now live under `legacy/manco/`; historical root paths are compatibility imports.

`LoadFilePopUp.py` and `notepad.py` likewise had no active imports. Their implementations are retained under `legacy/tools/` with compatibility imports at the historical paths.

No functionality was deleted. The MANCO notebook was adjusted only so its internal imports resolve within the quarantined package and to canonical current workspaces where appropriate.

Regression gate after Stage 86: 343 passed, 1 skipped, 0 failed.

## Stage 87 - give peak algorithms canonical analysis ownership

`peak_picker.py` and `peak_shape_estimator.py` are scientific algorithms rather than application-shell modules. Their canonical implementations now live at:

- `analysis/peak_picker.py`
- `analysis/peak_shape_estimator.py`

`gui/workspaces/peak_fit.py` imports these canonical modules. Historical root paths remain compatibility imports so downstream code and existing tests remain valid.

Architecture tests were added for the legacy quarantine and analysis ownership boundaries.

Final regression gate: 346 passed, 1 skipped, 0 failed.

Whole-tree `compileall` passes. AST parsing passes across 268 Python files.

## Next safe sequence

1. Audit the remaining root-level GUI/helper modules (`pdfViewer.py` in particular) and move active GUI helpers to `gui/` while retaining compatibility paths.
2. Continue grouping root compatibility facades and document the public compatibility surface.
3. Audit remaining `Frames/` wrappers for any implementation that has not yet been physically migrated or deliberately quarantined.
4. Only remove compatibility modules once external launch/import expectations are known.
