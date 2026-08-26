Pseudo2D independent absorption + dispersion test

Enabled by FitPhase 1 (reusing the existing experimental switch).

Stage 1: fit the absorptive pseudo-Voigt exactly as in the successful staged fitter.
Stage 2: freeze absorptive centre/width/amplitudes; add a dispersive component with signed per-slice amplitude, independent shared linewidth scale, and independent pseudo-Voigt mixing per resonance.
Stage 3: release absorptive amplitudes/widths/centres and dispersive amplitudes/widths together; solve A/D amplitudes simultaneously per pseudo slice and coordinate-refine shape parameters.

Output .out files report absorptive width/mixing, dispersive width/mixing, and separate Abs.Intensity / Disp.Intensity columns.
Console diagnostics report SSE input, zero-A, fixed-A+fit-D, and free-A+D.
