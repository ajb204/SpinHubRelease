# Stage 21 - peak workspace state ownership

`Frames/peakFrame.py` now resolves `ApplicationContext` on construction and
uses the context-owned `ProjectState` and `DataStore` for project/data ownership.
Topology and pseudo-axis decisions that previously reached through
`tabOne.state` now use the explicit project state. Legacy scientific callbacks
remain on `tabOne` until the peak-analysis service boundary is introduced.

Regression gate: **297 passed, 7 baseline failures, 1 skipped**.
