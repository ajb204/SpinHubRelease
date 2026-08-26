# Pseudo2D staged phase/shape fitting

`Protocol2PFit` now separates the F1 extraction/overlap radius from the spectral linewidth.

For each overlap group it performs:

1. zero-phase shape seeding from the physical `sig1`/`lor1` widths, including narrower multi-starts;
2. restrained coordinate refinement of per-peak width scale and centre;
3. one-time caching of absorptive and analytic-dispersive basis arrays;
4. phase refinement using only cached arrays;
5. a short width refinement with the fitted phases present;
6. one cache rebuild and final cached phase refinement.

The width scale multiplies `sig1` and `lor1` together, preserving their configured ratio. Widths are bounded below by the digital spacing and searched up to 2x the configured physical width. Peak centres are restrained to the smaller of 0.25 x FitF1 and 0.75 x the configured physical linewidth.

Diagnostic output reports input SSE, zero-phase shape-refined SSE, first phase SSE, final SSE, fitted widths in Hz and phases in degrees.

`FitPhase 0` retains zero-phase fitting but still uses the new shape seeding/refinement. `FitPhase 1` enables the staged cached phase passes.
