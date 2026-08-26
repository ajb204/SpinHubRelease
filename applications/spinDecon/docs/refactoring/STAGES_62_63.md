# Refactor stages 62-63

## Stage 62 - Pseudo2D fitting/review boundary

Expanded `PseudoAxisService` so the pseudo2D fitting inspector no longer needs
to own project-level operations when running under `ApplicationContext`.

New service responsibilities:

- canonical restrained-fit directory resolution;
- authoritative Full Peak List file resolution;
- explicit pseudo-series review persistence;
- workflow/status notification after review;
- readback of persisted review state.

`Pseudo2DFittingFrame` retains a legacy `tabOne` construction fallback only for
standalone compatibility. The normal application path uses `PseudoAxisService`.

Regression gate after this stage: **324 passed, 1 skipped**.

## Stage 63 - Canonical workspace imports

Added canonical application import modules:

- `gui/workspaces/oned.py`
- `gui/workspaces/pseudo2d.py`

`decon_tab.py` now imports `OneDFrame` and `Pseudo2D` through these canonical
workspace paths. Their implementations deliberately remain in `Frames/` while
legacy construction and source-contract tests are retired incrementally.

Added regression coverage for the pseudo fitting/review service boundary and
for the canonical notebook imports.

Final regression gate: **326 passed, 1 skipped, 0 failed**.

## Architectural state

The normal application import surface now treats 1D, pseudo2D, pseudo2D
diffusion, pseudo3D, projection, slices, Full3D and phasing as workspaces under
`gui/workspaces`, even where a compatibility implementation still resides in
`Frames`.

The Full Peak List remains authoritative. No connection/`conn_data` ownership
was added to the pseudo service.

## Next safe sequence

1. Continue reducing constructor-only `tabOne` access in OneD and Pseudo2D.
2. Introduce canonical peak/peak-fit workspace import surfaces before physical
   relocation of their source-sensitive implementations.
3. Audit the duplicate top-level `Frames/STD_frame.py` against the active
   `Frames/uSTA/STD_frame.py`; quarantine only after behavioural equivalence or
   obsolescence is established.
4. Keep `conn_data` as legacy compatibility state pending a future dedicated
   connection/NOE model.
