# Compatibility surfaces

The refactored application has canonical implementation packages under
`app/`, `analysis/`, `domain/`, `gui/`, `processing/`, `project/`, `workflow/`,
and `integrations/`.

Historical root modules, `Frames/`, and `misc/` are retained as thin import
facades where external callers may still depend on old module paths. They are
not architectural ownership locations and active code must not import them.

`legacy/` contains quarantined implementations retained for historical or
possible future recovery. In particular, CATIA, MANCO, old NOE code, old
Pseudo2D code, uSTA prototypes, and INDIANA are not active application
subsystems.

Compatibility facades should only be removed as an explicit release/API
compatibility decision. New functionality must never be added to them.
