# Refactoring stages 95-98

## Stage 95 - Canonical DataStore ownership

`project/data_store.py` is now the single DataStore implementation. The historical
`decon.data_store` module is a module alias retained for compatibility and class
identity. Active application code imports DataStore from `project.data_store`.
The source-sensitive deconvolution regression test now inspects the canonical
implementation.

## Stage 96 - Canonical project imports

Active app, GUI, integration, domain and workflow code was redirected away from
root project compatibility modules (`parameter_store`, `project_defaults`,
`project_summary`, `decon_service`, etc.) and toward the `project/` package.
Compatibility modules remain available for external callers and older scripts.

## Stage 97 - Canonical application shell

The live notebook/application shell moved from `decon_tab.py` to
`app/notebook.py`. `decon_tab.py` is now a compatibility entry point, preserving
`from decon.decon_tab import MyApp` used by the spinDecon launcher. Source-level
regression tests now inspect the canonical application shell.

## Stage 98 - Root compatibility boundary hardening

The remaining active project-summary workflow imports were redirected to the
canonical `workflow/` package. Architecture tests now prevent active packages
from importing the root compatibility layer.

## Validation

- 358 passed, 1 skipped, 0 failed.
- Whole-tree `compileall` passes.
- AST parsing passes across 273 Python files (including the new architecture test).

## Next safe sequence

1. Audit the remaining root compatibility modules and document the supported
   public compatibility surface.
2. Audit `Frames/` wrappers and add a test ensuring no active implementation can
   accidentally grow there again.
3. Split startup/environment preparation from `app/notebook.py` into a small
   application bootstrap module, while retaining `MyApp` compatibility.
4. Continue isolating optional legacy connectivity/NOE code without promoting
   `conn_data` back into canonical peak ownership.
