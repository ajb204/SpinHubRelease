from __future__ import annotations
from pathlib import Path

from spinDecon.processing.vpar_decon import inspect_acquisition
from spinDecon.project.state import ProjectState

from .models import AcquisitionRecord, ProjectRecord, ProjectResources

_RAW_NAMES = ('fid', 'ser', 'fid.gz', 'ser.gz')


def discover_acquisitions(root) -> list[AcquisitionRecord]:
    """Discover acquisitions independently of Decon projects.

    The filesystem walk only generates candidates. The authoritative decision
    about spectrometer type belongs to vpar_decon.inspect_acquisition().
    """
    root = Path(root)
    candidate_dirs = {p.parent for name in _RAW_NAMES for p in root.rglob(name) if p.is_file()}
    records = []
    for path in sorted(candidate_dirs, key=lambda p: str(p)):
        info = inspect_acquisition(path)
        if info is not None:
            records.append(AcquisitionRecord(path=path, info=info))
    return records


def inspect_project(parameter_file) -> ProjectRecord:
    """Load a deconParFile without requiring referenced resources to exist."""
    parameter_file = Path(parameter_file)
    state = ProjectState.from_parameter_file(parameter_file)
    raw_path = Path(state.raw_dir()) if state.raw_path else None
    spectrum_text = state.spectrum_path() if state.input_file else ''
    spectrum_path = Path(spectrum_text) if spectrum_text else None
    reference_text = state.reference_peak_path() if state.reference_peak_file else ''
    reference_peak_path = Path(reference_text) if reference_text else None
    full_text = state.full_peak_path() if state.full_peak_file else ''
    full_peak_path = Path(full_text) if full_text else None
    resources = ProjectResources(
        raw_path=raw_path,
        raw_available=bool(raw_path and raw_path.exists()),
        spectrum_path=spectrum_path,
        spectrum_available=bool(spectrum_path and spectrum_path.is_file()),
        reference_peak_path=reference_peak_path,
        reference_peaks_available=bool(reference_peak_path and reference_peak_path.is_file()),
        full_peak_path=full_peak_path,
        full_peaks_available=bool(full_peak_path and full_peak_path.is_file()),
    )
    return ProjectRecord(parameter_file=parameter_file, state=state, resources=resources)


def invalid_project(parameter_file, error) -> ProjectRecord:
    """Represent an unreadable/malformed deconParFile instead of hiding it."""
    return ProjectRecord(
        parameter_file=Path(parameter_file), state=None,
        resources=ProjectResources.empty(), error=str(error),
    )


def discover_projects(root) -> list[ProjectRecord]:
    """Discover every deconParFile, including malformed projects."""
    root = Path(root)
    records = []
    for p in sorted(root.rglob('deconParFile'), key=lambda p: str(p)):
        if not p.is_file():
            continue
        try:
            records.append(inspect_project(p))
        except (OSError, ValueError) as exc:
            records.append(invalid_project(p, exc))
    return records
