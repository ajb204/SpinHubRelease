from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from spinDecon.processing.vpar_decon import AcquisitionInfo
from spinDecon.project.state import ProjectState


class DatasetStatus(str, Enum):
    ACQUISITION_ONLY = 'acquisition_only'
    READY_TO_PROCESS = 'ready_to_process'
    PROCESSED = 'processed'
    PEAKS_AVAILABLE = 'peaks_available'
    SPECTRUM_ONLY = 'spectrum_only'
    RESOURCES_UNAVAILABLE = 'resources_unavailable'
    INVALID_PROJECT = 'invalid_project'


class ResourceState(str, Enum):
    NOT_CONFIGURED = 'not_configured'
    AVAILABLE = 'available'
    MISSING = 'missing'


@dataclass(frozen=True)
class AcquisitionRecord:
    path: Path
    info: AcquisitionInfo


@dataclass(frozen=True)
class ProjectResources:
    raw_path: Optional[Path]
    raw_available: bool
    spectrum_path: Optional[Path]
    spectrum_available: bool
    reference_peak_path: Optional[Path]
    reference_peaks_available: bool
    full_peak_path: Optional[Path]
    full_peaks_available: bool

    @staticmethod
    def _state(path: Optional[Path], available: bool) -> ResourceState:
        if path is None:
            return ResourceState.NOT_CONFIGURED
        return ResourceState.AVAILABLE if available else ResourceState.MISSING

    @property
    def raw_state(self) -> ResourceState:
        return self._state(self.raw_path, self.raw_available)

    @property
    def spectrum_state(self) -> ResourceState:
        return self._state(self.spectrum_path, self.spectrum_available)

    @property
    def reference_peaks_state(self) -> ResourceState:
        return self._state(self.reference_peak_path, self.reference_peaks_available)

    @property
    def full_peaks_state(self) -> ResourceState:
        return self._state(self.full_peak_path, self.full_peaks_available)

    @property
    def any_peaks_available(self) -> bool:
        return self.reference_peaks_available or self.full_peaks_available

    @classmethod
    def empty(cls) -> 'ProjectResources':
        return cls(None, False, None, False, None, False, None, False)


@dataclass
class ProjectRecord:
    parameter_file: Path
    state: Optional[ProjectState]
    resources: ProjectResources
    error: Optional[str] = None

    @property
    def valid(self) -> bool:
        return self.state is not None and self.error is None

    # Compatibility properties for existing callers.
    @property
    def raw_path(self): return self.resources.raw_path
    @property
    def spectrum_path(self): return self.resources.spectrum_path
    @property
    def raw_available(self): return self.resources.raw_available
    @property
    def spectrum_available(self): return self.resources.spectrum_available


@dataclass(frozen=True)
class DatasetCapabilities:
    can_create_project: bool
    can_open_project: bool
    can_process_raw: bool
    can_reprocess: bool
    can_view_spectrum: bool
    can_deconvolve: bool
    can_view_peaks: bool


@dataclass
class NMRDataset:
    acquisition: Optional[AcquisitionRecord] = None
    project: Optional[ProjectRecord] = None
    capabilities: Optional[DatasetCapabilities] = None
    status: Optional[DatasetStatus] = None
    status_text: str = ''
    # Relationship metadata is populated by the resolver.  A dataset remains
    # one actionable project, while related_projects exposes alternative
    # Decon projects that reference the same raw acquisition.
    related_projects: tuple[ProjectRecord, ...] = ()
    project_index: Optional[int] = None

    @property
    def project_count(self) -> int:
        return len(self.related_projects)

    @property
    def has_alternative_projects(self) -> bool:
        return self.project_count > 1
