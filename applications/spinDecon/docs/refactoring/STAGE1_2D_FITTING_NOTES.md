# Stage 1: native 2D fitting through the Pseudo3D workspace

Implemented for testing:

- True 2D spectra now receive a canonical logical fitting view with shape `[1, y, x]`.
- The physical spectrum remains 2D; no temporary/fake 3D NMRPipe file is created.
- Normal 2D loading opens both the existing Projections tab and the Fitting tab.
- Pseudo-axis navigation and Analysis controls are disabled for true 2D data.
- FUDA parameter generation writes `ZCOOR=2D` for physical 2D spectra while retaining the existing pseudo-axis label for physical pseudo3D data.
- The original spectrum path is passed unchanged as `SPECFILE`, so nmrPipeFit.py uses its native 2D path.

The existing pseudo3D behaviour is retained.

## Fix 2: shared contour threshold and 2D peak overlay
- Committing Threshold on the main NMR page now propagates the absolute threshold to the Pseudo3D/Fitting Contours `Min` field and redraws it. This is dimension-independent and applies to physical pseudo3D as well as the 2D adapter.
- The initial Pseudo3D contour minimum now uses the main threshold fraction rather than the raw spectrum maximum.
- Projected reference peak overlays are now published for physical 2D data using the canonical fitting-view X/Y labels. This fixes the Pseudo3D `Peaks` checkbox showing labels/markers for pseudo3D but no scatter markers for 2D.
