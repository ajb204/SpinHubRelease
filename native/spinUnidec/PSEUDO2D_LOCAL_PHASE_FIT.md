# Pseudo2D linewidth-local phase fit diagnostic

When `FitPhase 1` is enabled, `Protocol2PFit` first performs the existing zero-phase centre/linewidth fit. Centres and widths are then frozen and absorptive/dispersive peak shapes are cached.

The broad F1 fit radius continues to define overlap groups and the global amplitude least-squares solve. Phase is optimized with a separate local objective containing points within +/- 5 fitted linewidths of any resonance in the overlap group. Amplitudes are still re-solved using the full overlap region for every phase trial. This gives negative/dispersive structure near resonances substantially more leverage without changing the peak model or the successful linewidth fit.

The console reports the number of local phase points plus local SSE before/after phase optimization and the final global SSE.
