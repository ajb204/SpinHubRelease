# Pseudo2D shared dispersive contribution test

With `FitPhase 1`, the independent absorptive/dispersive pseudo2D model now constrains each resonance to one signed dispersive-to-absorptive ratio across the complete pseudo dimension.

The model is `I[p,z] * (A[p,x] + r[p] * D[p,x])`, where `r[p]` is shared by all pseudo slices. Absorptive and dispersive linewidth/shape parameters remain independent. Stage 2 fits the shared ratio with the absorptive solution frozen; stage 3 releases absorptive amplitudes and shape parameters while retaining the shared ratio constraint.

Output reports `D/A(shared)` and the per-slice `Disp.Intensity` column is derived as `D/A(shared) * Abs.Intensity`.
