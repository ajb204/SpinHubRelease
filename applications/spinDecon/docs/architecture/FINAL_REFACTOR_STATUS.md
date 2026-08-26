# Refactor status after Stage 147

## Status

The active code-level architectural refactor is complete for the current scope.
The modern application is organised around `app/`, `project/`, `domain/`,
`workflow/`, `analysis/`, `processing/`, and `gui/`. Historical root and `Frames/`
paths are compatibility surfaces rather than implementation owners.

The remaining textual `tabOne` occurrences in the workspaces addressed by the
final pass are constructor argument names, source-contract comments, or deliberate
compatibility adapters; Peak Review, PeakFit and Slice2D no longer store `tabOne`
as their scientific-state API.

## Explicitly deferred legacy scope

The following are intentionally not part of the active architecture and should be
revisited as separate feature-recovery projects:

- MAGMA
- Slice4D
- old uSTA tab
- old NOE/connectivity manager
- `conn_data`-centric behavior

`conn_data` remains compatibility state only. The Full Peak List is authoritative.
A future connection/NOE model should reference canonical Full Peak List identities.

## Required runtime validation before release

The headless suite cannot validate native wxPython lifecycle and interactive
Matplotlib behavior. Run these golden journeys on the supported macOS/wx setup:

1. Launch `spinDecon` and create/open a project.
2. 1D: load/process, inspect 1D view, peak-shape controls, save/reopen.
3. 2D: load/process, Peak Review, Full Peak List edits, Slice1D/Slice2D, save/reopen.
4. 3D: projections, reference peaks, Slice1D/Slice2D, fitting, save/reopen.
5. Pseudo2D: axis selection, reference peaks, intensity review, Diffusion/Decay.
6. Pseudo3D: spectral plane, reference peaks, fitting/review.
7. Phasing: open workspace, alter phase, verify data/project state survives save/reload.
8. Compatibility smoke: historical launcher `from decon.decon_tab import MyApp`.

Any runtime defect found here should be fixed without reopening ownership of
quarantined legacy modules unless that feature is intentionally being restored.

## Canonical peak-list rule (Stage 148E)

All current peak collections, at every supported dimensionality, use the common
peak-list model documented in `PEAK_LIST_MODEL.md`. A journey may own multiple
lists (for example reference/projection, full, and deconvolved), but these are
roles of the same representation. `conn_data`/`connEntry` are legacy-only and
must not be reintroduced into active peak ownership.

## Final canonical-journey architecture (Stage 148F+)

The active architecture is now defined around five supported journeys: 1D,
pseudo2D, physical 2D, pseudo3D and physical 3D.  Their spectrum/view and
peak-list ownership rules are specified in `CANONICAL_JOURNEYS.md` and encoded
in `domain.spectrum_policy`.

The decisive rule is authority rather than display: Full and Deconvolved lists
use the common peak schema; independent Reference lists exist only for journeys
that scientifically require a lower-dimensional projection; projection marker
lists used by protected/projected views are derived and non-authoritative.
Physical 2D is the singular case where Full/reference/projection concepts refer
to the same 2D peak collection, with Full as the sole authority.

## Projection fitting boundary

Peak fitting is a critical canonical operation and is owned by the PeakFit
workspace launched from the UniDec/NMR workspace. Spectrum deconvolution itself
uses the supported FUDA or decon reconstruction modes as appropriate.

The former Projection-workspace line-fitting UI is retired. Projection spectra
remain views/reference spectra according to the spectrum policy; they do not
own an independent fitting workflow. The historical Projection implementation
is preserved under `legacy/projection_line_fitting/` for reference.

The scientific `line_fitting/` package remains current support code where it is
used by active fitting workflows. Quarantining the old Projection UI does not
quarantine the scientific fitting routines themselves.
