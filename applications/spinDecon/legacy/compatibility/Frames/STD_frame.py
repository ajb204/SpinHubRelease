"""Compatibility import for the historical top-level STD workspace.

The maintained uSTA implementation lives in ``Frames/uSTA/STD_frame.py``.
The superseded source is retained under ``legacy/usta/std_frame_legacy.py``
for reference while useful functionality is migrated deliberately.
"""
from spinDecon.integrations.usta.workspace import STDFrame

__all__ = ["STDFrame"]
