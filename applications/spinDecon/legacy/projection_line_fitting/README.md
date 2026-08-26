# Historical projection-spectrum line fitting

This snapshot preserves the former Projection workspace line-fitting mode.
That mode fitted lines directly in projection spectra using `Unidec_line_fitting`.
It has been superseded by the canonical PeakFit workflow launched from the
UniDec/NMR workspace and is not part of the current Projection workspace.

`line_fitting/` itself remains active because current fitting workflows still
use its scientific routines; only the historical Projection-owned UI/path is
quarantined here.
