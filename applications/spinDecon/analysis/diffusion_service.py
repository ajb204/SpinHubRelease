"""Application boundary for pseudo-2D diffusion analysis.

This adapter keeps the diffusion workspace independent of the concrete NMR GUI
while the underlying project/data APIs are incrementally extracted.
"""
import os
import numpy


class DiffusionService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    @property
    def data(self):
        return numpy.asarray(self._legacy.data)

    @property
    def spectral_axis(self):
        return numpy.asarray(self._legacy.index1)

    @property
    def pseudo_axis(self):
        return numpy.asarray(self._legacy.index0)

    @property
    def peaks(self):
        return self._legacy.peak

    @property
    def labels(self):
        return tuple(getattr(self._legacy, 'labb', ()) or ())

    @property
    def spectral_label(self):
        labels = self.labels
        return str(labels[1]) if len(labels) > 1 else 'ppm'

    @property
    def spectral_bounds(self):
        return float(self._legacy.uc0min), float(self._legacy.uc0max)

    def threshold(self):
        return float(self._legacy.dmax) * float(self._legacy.threshBox.GetValue())

    def data_maximum(self):
        return float(self._legacy.dmax)

    def noise_value(self):
        return float(getattr(self._legacy, 'noiseVal', 0.0) or 0.0)

    def parameter_file(self):
        directory = self._legacy.dirBox.GetValue()
        name = self._legacy.deconParFile
        return os.path.join(directory, name), name

    def parameter(self, name, *, numeric=False, default=0):
        try:
            if numeric:
                return self._legacy.ParseFlt(self._legacy.deconParFile, name)
            return self._legacy.Parse(self._legacy.deconParFile, name)
        except Exception:
            return default
