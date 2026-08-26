# Refactor stages 88-90

## Stage 88 - PDF viewer ownership

Moved the active embedded PDF viewer implementation from the package root to
`gui/dialogs/pdf_viewer.py`. The historical `pdfViewer.py` path is now a small
compatibility import. NMR project-summary and MAGMA consumers import the
canonical dialog directly.

Regression protection was added for canonical ownership and active imports.

## Stage 89 - Workflow package activation

The canonical `workflow/` package had been created earlier but some of its
internal imports still referenced filenames from the old flat layout. Those
imports were repaired and the running notebook now imports the workflow overview
and legacy registry through `workflow/`.

The root `AnalysisMode` compatibility module remains the temporary identity
boundary because existing callers/tests monkeypatch that class directly. This is
intentional; changing class identity during a structural migration caused a
workflow-summary regression and was reverted before accepting the stage.

## Stage 90 - Frames isolation hardening

Removed the last active-tree references to `decon.Frames`. The historical STD
compatibility entry point now targets the canonical uSTA integration directly.
Added an architecture test requiring non-legacy/non-test code to remain free of
Frames imports.

## Validation

Final regression gate: **350 passed, 1 skipped, 0 failed**.

Additional validation:

- whole-tree `compileall`: pass
- AST parse: 270 Python files pass
- Full Peak List remains authoritative
- no new `conn_data` ownership introduced

## Next safe sequence

1. Convert duplicated root workflow modules to compatibility aliases once
   source-contract tests are redirected to `workflow/`.
2. Do the same for duplicated domain/project root modules in small groups,
   preserving class identity where compatibility requires it.
3. Audit the remaining root modules and `Frames/` wrappers for external API
   compatibility, then document a deprecation/removal policy rather than
   deleting wrappers abruptly.
4. Add an optional wx-enabled launch/import smoke test for environments where
   wxPython is available.
