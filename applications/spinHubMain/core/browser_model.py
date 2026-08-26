"""Presentation helpers for the SpinHub dataset browser.

This module deliberately has no wx dependency so browser behaviour can be
unit-tested independently of the GUI.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from .models import NMRDataset, ResourceState
from .workflows import recommended_workflow, workflow_actions, workflow_title


@dataclass(frozen=True)
class BrowserRow:
    name: str
    path: str
    source: str
    vendor: str
    dimension: str
    sequence: str
    raw: str
    spectrum: str
    projects: str
    status: str


def _resource_label(state: ResourceState) -> str:
    return {
        ResourceState.NOT_CONFIGURED: 'Not configured',
        ResourceState.AVAILABLE: 'Available',
        ResourceState.MISSING: 'Missing',
    }[state]


def row_for(dataset: NMRDataset) -> BrowserRow:
    acquisition = dataset.acquisition
    project = dataset.project
    if project is not None:
        display_path = project.parameter_file.parent
        name = display_path.name or str(display_path)
        source = 'Decon project'
    elif acquisition is not None:
        display_path = acquisition.path
        name = display_path.name or str(display_path)
        source = 'Acquisition'
    else:  # defensive; resolver should never produce an empty dataset
        display_path = Path('.')
        name = '(unknown)'
        source = 'Unknown'

    info = acquisition.info if acquisition is not None else None
    raw_state = project.resources.raw_state if project else (ResourceState.AVAILABLE if acquisition is not None else ResourceState.NOT_CONFIGURED)
    spectrum_state = project.resources.spectrum_state if project else ResourceState.NOT_CONFIGURED
    return BrowserRow(
        name=name,
        path=str(display_path),
        source=source,
        vendor=(info.vendor if info else ''),
        dimension=(str(info.dimension) if info and info.dimension is not None else ''),
        sequence=(info.sequence or '' if info else ''),
        raw=_resource_label(raw_state),
        spectrum=_resource_label(spectrum_state),
        projects=(f'{dataset.project_index} of {dataset.project_count}' if dataset.project and dataset.project_count > 1 else
                  ('1' if dataset.project else '0')),
        status=dataset.status_text,
    )



@dataclass(frozen=True)
class ResourceCard:
    key: str
    title: str
    state: str
    path: str
    action: str
    action_kind: str
    enabled: bool


def resource_cards(dataset: NMRDataset) -> tuple[ResourceCard, ...]:
    """Return the three workflow stages shown in the project dashboard.

    Keeping this wx-free makes the UX policy independently testable.
    """
    p = dataset.project
    c = dataset.capabilities
    if p is None:
        raw_path = str(dataset.acquisition.path) if dataset.acquisition else ''
        return (
            ResourceCard('raw', '1. Raw data', 'Available' if dataset.acquisition else 'Not configured', raw_path,
                         'Create Decon project', 'create', bool(c and c.can_create_project)),
            ResourceCard('spectrum', '2. Spectrum', 'Not configured', '', 'Process raw data', 'prepare', False),
            ResourceCard('peaks', '3. Peaks / analysis', 'Not configured', '', 'Analyse spectrum', 'decon', False),
        )
    r = p.resources
    raw_state = _resource_label(r.raw_state)
    spec_state = _resource_label(r.spectrum_state)
    peak_state = 'Available' if r.any_peaks_available else ('Not configured' if not (r.reference_peak_path or r.full_peak_path) else 'Missing')
    if not p.valid:
        return (
            ResourceCard('raw', '1. Raw data', 'Unavailable', '', 'Repair project first', 'none', False),
            ResourceCard('spectrum', '2. Spectrum', 'Unavailable', '', 'Repair project first', 'none', False),
            ResourceCard('peaks', '3. Peaks / analysis', 'Unavailable', '', 'Repair project first', 'none', False),
        )
    raw_action = ('Locate raw data…', 'locate_raw', True) if r.raw_state is ResourceState.MISSING else ('Process raw data', 'prepare', bool(c and c.can_process_raw))
    spec_action = ('Locate spectrum…', 'locate_spectrum', True) if r.spectrum_state is ResourceState.MISSING else ('Open / deconvolve', 'decon', bool(c and c.can_deconvolve))
    peak_action = ('Inspect peaks', 'inspect', bool(c and c.can_view_peaks))
    peak_path = r.full_peak_path or r.reference_peak_path
    return (
        ResourceCard('raw', '1. Raw data', raw_state, str(r.raw_path or ''), *raw_action),
        ResourceCard('spectrum', '2. Spectrum', spec_state, str(r.spectrum_path or ''), *spec_action),
        ResourceCard('peaks', '3. Peaks / analysis', peak_state, str(peak_path or ''), *peak_action),
    )

def detail_lines(dataset: NMRDataset) -> list[str]:
    lines = []
    if dataset.project:
        p = dataset.project
        lines.append(f'Project: {p.parameter_file}')
        if not p.valid:
            lines.append('Project state: invalid')
            lines.append(f'Error: {p.error or "could not read deconParFile"}')
        if p.resources.raw_path:
            lines.append(f'Raw acquisition: {p.resources.raw_path} ({_resource_label(p.resources.raw_state).lower()})')
        else:
            lines.append('Raw acquisition: not configured')
        if p.resources.spectrum_path:
            lines.append(f'Spectrum: {p.resources.spectrum_path} ({_resource_label(p.resources.spectrum_state).lower()})')
        else:
            lines.append('Spectrum: not configured')
        if p.resources.any_peaks_available:
            lines.append('Peak data: available')
        if dataset.has_alternative_projects:
            lines.extend(['', f'Decon projects for this acquisition: {dataset.project_count}'])
            for index, sibling in enumerate(dataset.related_projects, 1):
                marker = ' (selected)' if sibling is p else ''
                lines.append(f'  {index}. {sibling.parameter_file}{marker}')
    elif dataset.acquisition:
        lines.append(f'Acquisition: {dataset.acquisition.path}')

    if dataset.acquisition:
        i = dataset.acquisition.info
        bits = [i.vendor]
        if i.dimension is not None:
            bits.append(f'{i.dimension}D')
        if i.sequence:
            bits.append(i.sequence)
        lines.append('Spectrometer: ' + ' | '.join(x for x in bits if x))
        if i.nuclei:
            lines.append('Nuclei: ' + ', '.join(i.nuclei))
        if i.observation_frequency_mhz is not None:
            lines.append(f'Observation frequency: {i.observation_frequency_mhz:.2f} MHz')
        if i.temperature_k is not None:
            lines.append(f'Temperature: {i.temperature_k:.2f} K')
        if i.acquired_at:
            lines.append(f'Acquired: {i.acquired_at}')
    lines.append('Status: ' + dataset.status_text)

    actions = workflow_actions(dataset)
    if actions:
        lines.extend(['', 'Available workflows:'])
        for action in actions:
            marker = 'Available' if action.available else 'Unavailable'
            lines.append(f'  {action.title}: {marker} - {action.reason}')
    return lines


def primary_action_label(dataset: NMRDataset) -> str:
    c = dataset.capabilities
    if c and c.can_create_project:
        return 'Create Decon project'
    if c and c.can_open_project:
        key = recommended_workflow(dataset)
        if key is None:
            return 'Open project summary'
        return 'Continue: ' + workflow_title(key)
    return 'Open analysis'


def filter_datasets(datasets, query: str = '', source: str = 'All', status: str = 'All'):
    """Filter datasets using the same display information shown in the browser."""
    query = (query or '').strip().casefold()
    source = (source or 'All').strip()
    status = (status or 'All').strip()
    result = []
    for dataset in datasets:
        row = row_for(dataset)
        if source != 'All' and row.source != source:
            continue
        if status != 'All' and row.status != status:
            continue
        haystack = ' '.join((row.name, row.path, row.source, row.vendor,
                             row.dimension, row.sequence, row.raw,
                             row.spectrum, row.projects, row.status)).casefold()
        if query and query not in haystack:
            continue
        result.append(dataset)
    return result


def sort_datasets(datasets, column: str = 'Dataset', descending: bool = False):
    """Stable, case-insensitive sorting by a browser column."""
    fields = {
        'Dataset': 'name', 'Path': 'path', 'Source': 'source',
        'Spectrometer': 'vendor', 'Dim': 'dimension',
        'Pulse sequence': 'sequence', 'Raw': 'raw', 'Spectrum': 'spectrum',
        'Projects': 'projects', 'Status': 'status',
    }
    field = fields.get(column, 'name')
    def key(dataset):
        value = getattr(row_for(dataset), field)
        if field == 'dimension':
            try:
                return (0, int(value))
            except (TypeError, ValueError):
                return (1, str(value).casefold())
        return str(value).casefold()
    return sorted(datasets, key=key, reverse=descending)
