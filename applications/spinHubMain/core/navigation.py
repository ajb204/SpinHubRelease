"""Logical navigation model for the SpinHub browser.

The tree is intentionally independent of wx.  Decon projects are the primary
navigation objects: every discovered ``deconParFile`` appears under
``Analysed`` and, when it references raw NMR data, that acquisition is shown as
a child.  Raw acquisitions that are not referenced by any discovered project
appear separately under ``NMR data``.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .models import NMRDataset, ResourceState


@dataclass
class NavigationNode:
    label: str
    kind: str
    dataset: Optional[NMRDataset] = None
    children: list['NavigationNode'] = field(default_factory=list)
    status: str = ''


def _normal(path: Path) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def _display_path(path: Path, root: Path | None) -> str:
    """Prefer a path relative to the browser root, falling back to full path."""
    try:
        return str(path.relative_to(root)) if root else str(path)
    except ValueError:
        return str(path)


def _project_label(dataset: NMRDataset) -> str:
    project = dataset.project
    if project is None:
        return ''
    parent = project.parameter_file.parent
    label = parent.name or str(parent)
    if not project.valid:
        return f'{label} [invalid]'
    if project.resources.raw_state == ResourceState.MISSING or project.resources.spectrum_state == ResourceState.MISSING:
        return f'{label} [needs attention]'
    return label


def _project_raw_node(dataset: NMRDataset, root: Path | None) -> NavigationNode | None:
    """Represent the raw acquisition referenced by a Decon project.

    Use the resolved acquisition when available.  If the referenced raw data
    are missing, retain the configured path in the tree so the relationship is
    still visible and can be repaired from the selected project.
    """
    project = dataset.project
    if project is None:
        return None

    if dataset.acquisition is not None:
        path = dataset.acquisition.path
        label = _display_path(path, root)
    elif project.resources.raw_path is not None:
        path = project.resources.raw_path
        label = f'{_display_path(path, root)} [missing]'
    else:
        return None

    return NavigationNode(label, 'project_raw', dataset, status=dataset.status_text)


def build_navigation(datasets: list[NMRDataset], root: Path | None = None) -> NavigationNode:
    """Build the project-first browser hierarchy.

    ``Analysed`` contains every discovered Decon project, with its referenced
    raw acquisition as a child where configured.  ``NMR data`` contains only
    acquisition-only datasets, i.e. raw NMR data not referenced by any
    discovered project.  This changes navigation presentation only; the
    underlying datasets and project/acquisition relationships are untouched.
    """
    root_node = NavigationNode('Project navigator', 'root')
    analysed = NavigationNode('Analysed', 'group')
    nmr_data = NavigationNode('NMR data', 'group')

    for dataset in datasets:
        if dataset.project is not None:
            project_node = NavigationNode(
                _project_label(dataset), 'project', dataset, status=dataset.status_text
            )
            raw_node = _project_raw_node(dataset, root)
            if raw_node is not None:
                project_node.children.append(raw_node)
            analysed.children.append(project_node)
        elif dataset.acquisition is not None:
            nmr_data.children.append(NavigationNode(
                _display_path(dataset.acquisition.path, root),
                'acquisition', dataset, status=dataset.status_text
            ))

    # Deterministic presentation irrespective of scan order.  Keep Analysed
    # first so processed Decon projects have explicit navigation priority.
    analysed.children.sort(key=lambda n: n.label.casefold())
    nmr_data.children.sort(key=lambda n: n.label.casefold())
    root_node.children.extend((analysed, nmr_data))
    return root_node


def browser_summary(datasets: list[NMRDataset]) -> dict[str, int]:
    """Small dashboard counts for the browser header."""
    acquisitions = len({_normal(d.acquisition.path) for d in datasets if d.acquisition is not None})
    projects = sum(1 for d in datasets if d.project is not None)
    attention = sum(1 for d in datasets if d.project is not None and (
        not d.project.valid or
        d.project.resources.raw_state == ResourceState.MISSING or
        d.project.resources.spectrum_state == ResourceState.MISSING
    ))
    ready = sum(1 for d in datasets if d.project is not None and d.project.valid and
                d.project.resources.raw_state == ResourceState.AVAILABLE and
                d.project.resources.spectrum_state == ResourceState.AVAILABLE)
    return {'acquisitions': acquisitions, 'projects': projects, 'ready': ready, 'attention': attention}
