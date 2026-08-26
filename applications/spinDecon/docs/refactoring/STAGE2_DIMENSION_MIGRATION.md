# Stage 2: canonical project/load dimensionality

Implemented boundary contract:

- `ProjectState.dimension` remains for compatibility but now **always means spectral dimensions**.
- New code can use `ProjectState.spectral_dimensions` explicitly.
- `ProjectState.physical_dimensions` is derived as spectral + pseudo.
- `ProjectState.topology()` exposes the Stage 1 `DatasetTopology`.
- Parameter-file loading now reads the persisted `pseudo` flag.
- Legacy pseudo projects are marked `legacy_unresolved` until spectrum metadata is available.
- `canonicalize_loaded_dimensions()` resolves the historical physical `dim=3,pseudo=1`
  representation only when loaded physical dimensionality and a recognised real-axis label
  provide evidence. It normalises the state to `spectral=2,pseudo=True,physical=3`.
- `deconFrame.makeinp()` uses that load boundary and no longer promotes a 2D+pseudo dataset
  to GUI/project `dim=3`.

This stage intentionally leaves downstream GUI compatibility code in place. Stage 3 can now
make Main NMR and Workflow consume the canonical state without having to reinterpret loaded data.
