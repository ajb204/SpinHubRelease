"""Compatibility facade for the retired 4D slice workspace.

The historical 4D viewer is tightly coupled to MAGMA and legacy ``conn_data``.
It is preserved under :mod:`decon.legacy.slice4d.workspace` until a future 4D
viewer can be rebuilt around the authoritative Full Peak List and a dedicated
connection/NOE model.
"""
from spinDecon.legacy.slice4d.workspace import SliceFrame4D
__all__ = ["SliceFrame4D"]
