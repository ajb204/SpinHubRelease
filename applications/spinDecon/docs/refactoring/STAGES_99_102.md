# Stages 99-102

## Stage 99 - ProjectService canonical boundaries

`project/service.py` no longer reaches back through root compatibility modules. It now uses
`workflow/legacy_registry.py` directly and opens the GUI through `app/notebook.py`.
This also fixes a latent relative import (`.decon_tab`) that could fail when `ProjectService.open()`
was exercised outside the normal launcher path.

## Stage 100 - Explicit workflow legacy bridge

Workflow status evaluation now prefers `ApplicationContext.legacy_nmr_workspace`, while still
accepting `context.tabOne` for old hosts/tests. This aligns workflow evaluation with the explicit
migration boundary introduced earlier instead of inventing a second implicit context API.

## Stage 101 - Analysis service composition root

The notebook no longer imports and constructs every analysis service individually. New
`analysis/services.py` owns service composition and attaches the current service set to
`ApplicationContext`. The NMR workspace remains an explicit temporary migration input.

This reduces the application shell's dependency surface and provides one place to replace
legacy-workspace-backed services with domain/data-backed implementations over time.

## Stage 102 - Dependency-direction guard

Added an AST-based architecture regression test ensuring active packages do not import from
`decon.Frames` or from the root compatibility modules. Root modules and `Frames/` can therefore
remain stable backwards-compatible surfaces without becoming implementation dependencies again.

## Validation

The complete regression suite is green after every accepted boundary. Final gate for this block:
360 passed, 1 skipped, 0 failed before the dependency guard; the guard adds one further passing test.
Full-tree byte compilation and AST parsing are also required before packaging.
