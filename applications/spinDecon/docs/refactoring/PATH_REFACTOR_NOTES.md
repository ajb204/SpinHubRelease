# Project path refactor test build

This build centralises project path ownership in `ProjectState` and exposes the three project directories on the NMR tab:

- Working dir: project/CWD base
- OutPath: raw/FID data directory, relative to Working dir unless absolute
- SpecPath: processed spectrum directory, relative to Working dir unless absolute

Spectrum-associated GUI values (`nmrPipe file`, Reference 2D peak list, Full nD peak list) are stored as filenames and resolved under SpecPath for I/O.

Process no longer creates/owns OutPath or SpecPath text controls. Process, conversion, processing scripts and projection processing obtain resolved paths from shared state.

Decon receives resolved SpecPath-qualified spectrum/reference paths. Its `.decon` spectrum is loaded by the existing decon-output loader and the dimensionality-matched full peak list is loaded from the same SpecPath namespace. The visible Full nD control remains a filename rather than being replaced by an absolute path.

Legacy parameter keys remain compatible: `indir`, `fiddir`, `specPath`, `infile`, `peakfile`.


## Full nD path semantics

`Full nD` is now persisted as `fullPeakFile` and, like the nmrPipe and Reference 2D controls, is a path relative to SpecPath. Subdirectories are preserved. On an ordinary spectrum read an explicitly configured Full nD value is retained; if it is empty, the conventional `<spectrum>.<n>D.list` value is derived. After a successful decon run the control/state are intentionally updated to the newly generated full-dimensional peak list and that list is loaded using its resolved SpecPath path.
