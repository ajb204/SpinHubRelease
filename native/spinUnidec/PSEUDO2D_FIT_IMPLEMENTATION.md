# Pseudo2D restrained fitting

`pseudo2DFit=1`, `recon=1`, and `FIT=1` dispatch physical 2D data to
`Protocol2PFit`. The slow NMRPipe axis is treated as the real/pseudo series and
the fast axis as the single spectral F1 dimension. The supplied 1D peak list
provides fixed peak identifiers and ppm positions. `FitF1` is the extraction
radius. Overlapping peak windows are fitted simultaneously by linear least
squares for each pseudo slice, producing independent slice amplitudes and ESDs.
Results are written to `fit/<peak>.out` and `fit/<peak>.dat` without changing
Protocol1D/2D/3P behaviour.
