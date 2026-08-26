# Stage 6 viewer dimensionality migration

Viewer code now obtains dimensionality through `viewer_dimension_contract` and `DatasetTopology`.

Rules:
- `spectral_dim_count` controls spectral/view dimensionality.
- `physical_dim_count` describes the backing array.
- `spectral_axes` map spectral coordinates to physical array indices.
- `pseudo_axis` identifies the real sampled axis explicitly.
- `data.ndim` is used only to validate/describe an actual array, not to infer spectral dimensionality.

Migrated entry points: Projection, Pseudo2D, Pseudo3D, OneDView, slicePlot, slicePlot2D and PhasingSpectra.
Projection's physical pseudo-3D case is now identified as 2 spectral + pseudo = 3 physical, rather than legacy `dim == 3`, and spectral labels are selected through topology axis identities.
