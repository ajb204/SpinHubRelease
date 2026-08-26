# Stage 7 dimensionality contract

The project now has one scientific definition system:

- `spectral_dim_count`: number of frequency/chemical-shift axes.
- `has_pseudo_axis`: whether one additional sampled real axis exists.
- `physical_dim_count`: `spectral_dim_count + int(has_pseudo_axis)`.
- `AxisSpec.physical_index`: ndarray/NMRPipe axis identity.
- `AxisSpec.spectral_index`: position amongst spectral axes only; absent for the pseudo axis.

## GUI ownership

Main NMR and Workflow edit spectral count and pseudo presence. Conversion consumes physical axes. Processing, peaks, fitting and phasing consume spectral axes. Viewers/projections consume explicit `AxisSpec` identities. `data.ndim` describes an actual array and is not a synonym for spectral dimensionality.

## Compatibility boundary

Old project dimensionality is decoded once by `ProjectState.canonicalize_loaded_dimensions`. The old vendor `vpar_decon` backend still has a private `2p`/`3p` wire encoding; `processing.dimension_contract.legacy_vpar_dimension` is the only adapter from canonical topology into that backend. It must not be used as GUI state.

## Enforcement

`dimension_guard.assert_full_dataset_contract` validates full loaded spectra against `DatasetTopology`. Derived projections/slices validate their own local ndarray shapes. New GUI code must not add pseudo arithmetic to a spectral count; ask the topology for physical counts/axes instead.
