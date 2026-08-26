# Pseudo2D shared Gaussian/Lorentzian test

This numerical test keeps the shared signed D/A ratio across pseudo slices and the independent absorptive/dispersive linewidths, but replaces separate absorptive/dispersive pseudo-Voigt mixing with one fitted `g1(shared)` per resonance.

Stages:
1. Zero-dispersion absorptive centre/linewidth seed (existing behaviour).
2. Frozen absorption + dispersive linewidth and shared D/A fit, initially using the configured `voigt1`.
3. Released A+D fit where absorptive width, dispersive width, D/A, centre, and one common Gaussian/Lorentzian mixing parameter are refined. The same `g1(shared)` is used in both A and D functions.

The console and `.out` files report `g(shared)` / `g1(shared)` rather than separate `g1(abs)` and `g1(disp)`.
