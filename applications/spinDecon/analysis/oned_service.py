"""Application boundary for the 1D spectrum workspace.

The current implementation deliberately delegates to the legacy NMR workspace.
This gives the 1D GUI a stable application API while scientific state is moved
out of ``deconFrame`` incrementally.
"""

class OneDService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def data(self):
        return self._legacy.data

    @property
    def index(self):
        return self._legacy.index0

    @property
    def peaks(self):
        return self._legacy.peak

    @property
    def labels(self):
        return tuple(self._legacy.labb)

    @property
    def deconvolved_data(self):
        return self._legacy.datadec

    @property
    def deconvolution_enabled(self):
        return self._legacy.DECON == 1

    @property
    def axis_limits(self):
        return self._legacy.uc0min, self._legacy.uc0max

    def threshold(self):
        return float(self._legacy.dmax) * float(self._legacy.threshBox.GetValue())
