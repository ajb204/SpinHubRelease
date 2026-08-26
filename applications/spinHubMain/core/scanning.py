"""Incremental, cancellable recursive scanning for SpinHub.

The scanner has no wx dependency. GUI callers may run ``scan_tree`` in a worker
thread and marshal callbacks onto the UI thread.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from threading import Event
from typing import Callable

from .discovery import inspect_project, invalid_project
from .models import AcquisitionRecord, ProjectRecord
from spinDecon.processing.vpar_decon import inspect_acquisition

_RAW_NAMES = frozenset(('fid', 'ser', 'fid.gz', 'ser.gz'))


@dataclass(frozen=True)
class ScanProgress:
    visited: int = 0
    acquisition_candidates: int = 0
    acquisitions: int = 0
    projects: int = 0
    current_path: str = ''


@dataclass(frozen=True)
class ScanResult:
    acquisitions: tuple[AcquisitionRecord, ...]
    projects: tuple[ProjectRecord, ...]
    progress: ScanProgress
    cancelled: bool = False


def scan_tree(root, cancel_event: Event | None = None,
              progress_callback: Callable[[ScanProgress], None] | None = None) -> ScanResult:
    """Scan *root* once, reporting progress and honouring cancellation.

    Candidate acquisition directories are deduplicated before inspection. A
    deconParFile is inspected independently, so projects remain discoverable
    even when their configured raw data is offline.
    """
    root = Path(root)
    cancel_event = cancel_event or Event()
    acquisitions: list[AcquisitionRecord] = []
    projects: list[ProjectRecord] = []
    candidate_dirs: set[Path] = set()
    progress = ScanProgress()

    def report(path=''):
        nonlocal progress
        progress = ScanProgress(
            visited=progress.visited,
            acquisition_candidates=len(candidate_dirs),
            acquisitions=len(acquisitions), projects=len(projects),
            current_path=str(path),
        )
        if progress_callback:
            progress_callback(progress)

    for path in root.rglob('*'):
        if cancel_event.is_set():
            report(path)
            return ScanResult(tuple(acquisitions), tuple(projects), progress, True)
        progress = ScanProgress(progress.visited + 1, len(candidate_dirs),
                                len(acquisitions), len(projects), str(path))
        if not path.is_file():
            if progress_callback and progress.visited % 100 == 0:
                progress_callback(progress)
            continue
        if path.name in _RAW_NAMES:
            candidate_dirs.add(path.parent)
        if path.name == 'deconParFile':
            try:
                projects.append(inspect_project(path))
            except (OSError, ValueError) as exc:
                # Broken projects remain visible and actionable in the browser.
                projects.append(invalid_project(path, exc))
        if progress_callback and progress.visited % 50 == 0:
            report(path)

    for candidate in sorted(candidate_dirs, key=lambda p: str(p)):
        if cancel_event.is_set():
            report(candidate)
            return ScanResult(tuple(acquisitions), tuple(projects), progress, True)
        try:
            info = inspect_acquisition(candidate)
        except (OSError, ValueError):
            info = None
        if info is not None:
            acquisitions.append(AcquisitionRecord(path=candidate, info=info))
        report(candidate)

    report(root)
    return ScanResult(tuple(acquisitions), tuple(projects), progress, False)
