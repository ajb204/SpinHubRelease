# Stage 3 - Main NMR and Workflow canonical dimensionality

Stage 3 makes the main NMR dataset selector and Workflow selector use the same
contract established by Stages 1-2.

## Canonical UI contract

- The NMR dimension radio box is the number of **spectral dimensions**.
- Workflow "Spectral dimensions" is the same value.
- `pseudoAxis` / "Contains a real pseudo-axis" adds exactly one **physical**
  real axis.
- `physical_dimensions = spectral_dimensions + int(pseudo_axis)`.

Workflow now displays the derived physical dimension count read-only so the
three concepts are visible without changing the user's source-of-truth inputs.

## Removed compatibility interpretation

`AnalysisMode` no longer has `canonical_physical_pseudo`. It never interprets
GUI/project dimension as a physical count. Legacy projects are canonicalized
at the Stage 2 load boundary.

Workflow pseudo-3D routing now requires exactly 2 spectral dimensions plus a
pseudo axis. The former `(2, 3)` compatibility route has been removed.

## Deferred restriction

The existing GUI still rejects 4 spectral dimensions + pseudo in
`deconFrame.pseudoBoxCheck`. Stage 1 can represent this topology, but enabling
it in the live GUI is intentionally deferred until Process/Conversion can
represent five physical axes safely. Main NMR and Workflow reject it
consistently in the meantime.
