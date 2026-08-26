# Pseudo2D phased peak fitting

`Protocol2PFit` can now relax one spectral phase per restrained resonance after
first solving the existing absorptive (`phase = 0`) fit.

Enable in `decon.init` with:

    FitPhase 1

The phase is shared across all pseudo slices for a resonance.  At every phase
trial the per-slice amplitudes are re-solved exactly by linear least squares.
The phase is restricted to +/-90 degrees to remove the 180-degree degeneracy
with signed amplitudes.

The dispersive pseudo-Voigt is evaluated analytically as the same linear
Gaussian/Lorentzian mixture as the absorptive peak.  The Lorentzian quadrature
is elementary; the Gaussian quadrature uses Dawson's integral, evaluated with
a compact high-accuracy approximation (no FFT Hilbert transform in the fit).

Fitted phase is written as `phase(deg)` in each `fit/<peak>.out` file.  With
`FitPhase 0` (the default), pseudo2D fitting retains the previous zero-phase
behaviour.
