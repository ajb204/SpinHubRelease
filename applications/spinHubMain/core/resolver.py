from __future__ import annotations
from pathlib import Path

from .models import AcquisitionRecord, ProjectRecord, NMRDataset
from .capabilities import annotate_dataset
from spinDecon.processing.vpar_decon import inspect_acquisition


def _normal(path: Path) -> Path:
    """Normalise for matching without requiring the path to exist."""
    return Path(path).expanduser().resolve(strict=False)


def resolve_datasets(acquisitions: list[AcquisitionRecord], projects: list[ProjectRecord]) -> list[NMRDataset]:
    """Link projects to acquisitions using the raw path stored in deconParFile.

    Filesystem proximity is intentionally ignored. Unmatched projects and
    acquisitions remain first-class datasets.  When several projects point to
    one acquisition, each project remains independently actionable but carries
    relationship metadata describing its sibling projects.
    """
    by_path = {_normal(a.path): a for a in acquisitions}
    matched = set()
    datasets = []

    for project in projects:
        acquisition = None
        if project.raw_path is not None:
            raw_key = _normal(project.raw_path)
            acquisition = by_path.get(raw_key)
            # A project may reference raw data outside the selected scan root.
            # Inspect that exact acquisition, but never recursively traverse its
            # parent/archive: project references are direct links, not scan roots.
            if acquisition is None and project.raw_available:
                try:
                    info = inspect_acquisition(project.raw_path)
                except (OSError, ValueError):
                    info = None
                if info is not None:
                    acquisition = AcquisitionRecord(path=Path(project.raw_path), info=info)
                    by_path[raw_key] = acquisition
        if acquisition is not None:
            matched.add(_normal(acquisition.path))
        datasets.append(annotate_dataset(NMRDataset(acquisition=acquisition, project=project)))

    # Build relationships only after all project references have been resolved.
    # Group by canonical raw path rather than object identity so discovered and
    # directly-inspected acquisitions behave identically.
    projects_by_raw = {}
    for dataset in datasets:
        if dataset.acquisition is None or dataset.project is None:
            continue
        key = _normal(dataset.acquisition.path)
        projects_by_raw.setdefault(key, []).append(dataset.project)

    for dataset in datasets:
        if dataset.acquisition is None or dataset.project is None:
            continue
        siblings = tuple(projects_by_raw[_normal(dataset.acquisition.path)])
        dataset.related_projects = siblings
        dataset.project_index = siblings.index(dataset.project) + 1

    for acquisition in acquisitions:
        if _normal(acquisition.path) not in matched:
            datasets.append(annotate_dataset(NMRDataset(acquisition=acquisition, project=None)))

    return datasets
