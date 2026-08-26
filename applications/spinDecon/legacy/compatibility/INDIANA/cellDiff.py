"""Compatibility facade for the quarantined historical INDIANA backend."""
import sys as _sys
from spinDecon.legacy.indiana import cell_diff as _impl
_sys.modules[__name__] = _impl
