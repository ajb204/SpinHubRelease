"""Public, path-safe integration boundary for opening and creating Decon projects.

This module deliberately has no wx import at module import time.  File managers
such as SpinHub can inspect/create projects without changing process cwd or
starting a shell command.  GUI imports happen only when ``open`` is requested.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
from typing import Optional

from .defaults import available_cpu_count
from .parameter_store import update_parameter_file
from .state import ProjectState
from spinDecon.workflow.registry import WORKFLOW_BY_KEY, WORKFLOW_REGISTRY


@dataclass(frozen=True)
class ProjectOpenResult:
    parameter_file: Path
    workflow: Optional[str]
    frame: object


DEFAULT_PARAMETER_NAME = "spinHub.par"
LEGACY_PARAMETER_NAMES = ("deconParFile",)


class ProjectService:
    """Stable API used by external project browsers such as SpinHub."""

    def discover_parameter_file(self, directory: str | Path) -> Path | None:
        """Return the preferred project parameter file in *directory*.

        Legacy ``deconParFile`` takes precedence when both names exist so an
        established project is never silently switched to the new filename.
        """
        root = Path(directory).expanduser().resolve(strict=False)
        for name in (*LEGACY_PARAMETER_NAMES, DEFAULT_PARAMETER_NAME):
            candidate = root / name
            if candidate.is_file():
                return candidate
        return None

    def create_initial_parameter_file(
        self, directory: str | Path, acquisition_path: str | Path, *,
        dimension: int, pseudo_axis: bool = False, spec_path: str = "./spec",
        parameter_name: str = DEFAULT_PARAMETER_NAME,
    ) -> ProjectState:
        """Create the minimal, valid system file required before GUI startup."""
        dim = int(dimension)
        if dim < 1:
            raise ValueError("spectral dimensionality must be at least 1")
        root = Path(directory).expanduser().resolve(strict=False)
        root.mkdir(parents=True, exist_ok=True)
        raw = Path(acquisition_path).expanduser().resolve(strict=False)
        parameter_file = root / parameter_name
        update_parameter_file(parameter_file, {
            "indir": str(root), "fiddir": str(raw), "specPath": spec_path,
            "dim": str(dim), "pseudo": int(bool(pseudo_axis)),
        })
        spec = Path(spec_path).expanduser()
        if not spec.is_absolute():
            spec = root / spec
        spec.resolve(strict=False).mkdir(parents=True, exist_ok=True)
        return self.inspect(parameter_file)

    def prepare_open(self, parameter_file: str | Path, *, change_cwd: bool = True) -> ProjectState:
        """Load the canonical project state and establish its runtime context.

        ``ProjectState`` is the authority for project paths.  Changing cwd is a
        deliberate compatibility step for legacy UniDecNMR code; new code must
        resolve resources through the returned state instead of cwd.
        """
        path = Path(parameter_file).expanduser().resolve(strict=False)
        state = self.inspect(path)
        if change_cwd:
            os.chdir(state.working_dir)
        return state

    def inspect(self, parameter_file: str | Path) -> ProjectState:
        return ProjectState.from_parameter_file(Path(parameter_file).expanduser())

    def create(
        self,
        acquisition_path: str | Path,
        destination: str | Path,
        *,
        dimension: int | None = None,
        parameter_name: str = DEFAULT_PARAMETER_NAME,
    ) -> ProjectState:
        """Create a project that *references* an acquisition without moving it.

        ``destination`` is the project directory.  Existing files and the raw
        acquisition are never reorganised.  If the parameter file already
        exists only the supplied core path settings are updated.
        """
        raw = Path(acquisition_path).expanduser().resolve(strict=False)
        project_dir = Path(destination).expanduser().resolve(strict=False)
        project_dir.mkdir(parents=True, exist_ok=True)
        parameter_file = project_dir / parameter_name

        updates = {
            "indir": str(project_dir),
            "fiddir": str(raw),
            "specPath": "./spec",
        }
        if dimension is not None and int(dimension) > 0:
            updates["dim"] = str(int(dimension))
        update_parameter_file(parameter_file, updates)
        return self.inspect(parameter_file)

    def relink_raw(self, parameter_file: str | Path, acquisition_path: str | Path, *, expected_old_path: str | Path | None = None) -> ProjectState:
        """Update only ``fiddir`` after an explicit caller confirmation.

        ``expected_old_path`` provides optimistic concurrency protection: if
        the project was changed since it was inspected, the relink is refused.
        """
        parameter_file = Path(parameter_file).expanduser().resolve(strict=False)
        current = self.inspect(parameter_file)
        if expected_old_path is not None:
            expected = Path(expected_old_path).expanduser().resolve(strict=False)
            actual = Path(current.raw_path).expanduser().resolve(strict=False) if current.raw_path else None
            if actual != expected:
                raise RuntimeError(f'Configured raw path changed since inspection: {actual}')
        raw = Path(acquisition_path).expanduser().resolve(strict=False)
        if not raw.exists():
            raise FileNotFoundError(raw)
        update_parameter_file(parameter_file, {"fiddir": str(raw)})
        return self.inspect(parameter_file)

    def relink_spectrum(self, parameter_file: str | Path, spectrum_path: str | Path, *, expected_old_path: str | Path | None = None) -> ProjectState:
        """Point a project at an existing main spectrum file.

        Both ``specPath`` and ``infile`` are updated together.  The optional
        expected path provides the same optimistic concurrency protection as
        :meth:`relink_raw`.
        """
        parameter_file = Path(parameter_file).expanduser().resolve(strict=False)
        current = self.inspect(parameter_file)
        if expected_old_path is not None:
            expected = Path(expected_old_path).expanduser().resolve(strict=False)
            actual_path = current.spectrum_path()
            actual = Path(actual_path).expanduser().resolve(strict=False) if actual_path else None
            if actual != expected:
                raise RuntimeError(f'Configured spectrum changed since inspection: {actual}')
        spectrum = Path(spectrum_path).expanduser().resolve(strict=False)
        if not spectrum.is_file():
            raise FileNotFoundError(spectrum)
        update_parameter_file(parameter_file, {
            "specPath": str(spectrum.parent),
            "infile": spectrum.name,
        })
        return self.inspect(parameter_file)

    def workflows(self) -> tuple[dict, ...]:
        return tuple(dict(item) for item in WORKFLOW_REGISTRY)

    def open(self, parameter_file: str | Path, *, workflow: str | None = None, show: bool = True) -> ProjectOpenResult:
        """Open Decon in-process using one canonical project-opening path."""
        path = Path(parameter_file).expanduser().resolve(strict=False)
        if workflow is not None and workflow not in WORKFLOW_BY_KEY:
            raise ValueError(f"Unknown Decon workflow: {workflow}")

        state = self.prepare_open(path, change_cwd=True)

        # GUI application ownership lives in the app layer; this project service
        # remains importable and inspectable without wxPython.
        from spinDecon.app.launcher import open_project_gui
        frame = open_project_gui(path, state=state, workflow=workflow, show=show)
        return ProjectOpenResult(path, workflow, frame)
