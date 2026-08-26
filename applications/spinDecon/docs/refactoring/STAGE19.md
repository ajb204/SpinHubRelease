# Stage 19 - projection service boundary

Added `analysis/projection_service.py` and exposed it through `ApplicationContext`.
The Projection workspace now obtains central spectrum/projection views, pseudo-2D
projection payloads, peak payloads, labels, spectrum data, overlays, reference
peaks and common peak/connection state through that boundary. Direct `tabOne`
references in the module were reduced substantially while preserving legacy
callbacks that do not yet have an application service owner.

Regression gate: **297 passed, 7 baseline failures, 1 skipped**.
