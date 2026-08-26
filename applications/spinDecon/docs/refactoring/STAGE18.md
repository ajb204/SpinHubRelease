# Stage 18 - 1D workspace service boundary

Introduced `analysis/oned_service.py` and exposed it through `ApplicationContext`.
`Frames/OneDView.py` now obtains spectrum arrays, axes, labels, threshold,
deconvolution data, peaks and connection overlays through the service instead
of treating `tabOne` as its scientific data API. The constructor still accepts
`tabOne` for backwards-compatible construction and context discovery.

Regression gate: **297 passed, 7 baseline failures, 1 skipped**.
