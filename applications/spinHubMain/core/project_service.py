"""SpinHub adapter for the public Decon project integration API."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from spinDecon.project.service import ProjectService as DeconProjectService


class ProjectService:
    def __init__(self, backend: Optional[DeconProjectService] = None):
        self.backend = backend or DeconProjectService()

    def inspect(self, parameter_file):
        return self.backend.inspect(parameter_file)

    def create_for_acquisition(self, acquisition_path, project_dir=None, *, dimension=None):
        raw = Path(acquisition_path).expanduser().resolve(strict=False)
        destination = (Path(project_dir).expanduser() if project_dir is not None else raw / 'decon')
        return self.backend.create(raw, destination, dimension=dimension)

    def relink_raw(self, parameter_file, acquisition_path, *, expected_old_path=None):
        """Validate an NMR acquisition, then update the project's raw reference."""
        from spinDecon.processing.vpar_decon import inspect_acquisition
        raw = Path(acquisition_path).expanduser().resolve(strict=False)
        if inspect_acquisition(raw) is None:
            raise ValueError(f'Not a recognised NMR acquisition: {raw}')
        return self.backend.relink_raw(
            parameter_file, raw, expected_old_path=expected_old_path)

    def relink_spectrum(self, parameter_file, spectrum_path, *, expected_old_path=None):
        """Update the project's configured main spectrum through Decon's API."""
        spectrum = Path(spectrum_path).expanduser().resolve(strict=False)
        return self.backend.relink_spectrum(
            parameter_file, spectrum, expected_old_path=expected_old_path)

    def open_project(self, parameter_file, *, workflow=None, show=True):
        return self.backend.open(parameter_file, workflow=workflow, show=show)

    def open_dataset(self, dataset, *, workflow=None, show=True):
        """Open a dataset at the centrally recommended workflow."""
        if dataset.project is None:
            raise ValueError('Dataset has no Decon project')
        from .workflows import recommended_workflow
        selected = workflow if workflow is not None else recommended_workflow(dataset)
        return self.open_project(dataset.project.parameter_file, workflow=selected, show=show)

    @staticmethod
    def parameter_file_for_legacy_row(path) -> Path:
        """Resolve a parameter file from the path stored by the legacy table."""
        p = Path(path).expanduser().resolve(strict=False)
        candidates = [p / 'deconParFile', p.parent / 'deconParFile']
        if p.name in ('fid', 'fid.gz', 'ser', 'ser.gz'):
            candidates.insert(0, p.parent / 'deconParFile')
        # Historical SpinHub rows commonly point at <project>/raw.
        if p.name == 'raw':
            candidates.insert(0, p.parent / 'deconParFile')
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return candidates[0]
