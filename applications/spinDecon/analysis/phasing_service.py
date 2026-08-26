"""GUI-independent data boundary for the phasing workspace."""
from __future__ import annotations

from spinDecon.domain.dimensions.viewer_contract import topology_for


class PhasingService:
    def __init__(self, workspace):
        self.workspace = workspace

    @property
    def topology(self):
        return topology_for(self.workspace)

    @property
    def labels(self):
        return list(getattr(self.workspace, "labb", []) or [])

    @property
    def working_directory(self):
        control = getattr(self.workspace, "dirBox", None)
        return str(control.GetValue()) if control is not None else ""

    @property
    def peaks(self):
        return getattr(self.workspace, "peak", [])

    @property
    def pseudo_spectrum(self):
        return getattr(self.workspace, "pseudo_spectrum", None)

    def axis_limits(self, dimension):
        return (
            float(getattr(self.workspace, f"uc{dimension}min")),
            float(getattr(self.workspace, f"uc{dimension}max")),
        )
