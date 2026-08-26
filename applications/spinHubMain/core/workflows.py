"""Contextual Decon workflow policy for SpinHub.

This module is deliberately GUI-free.  It converts dataset capabilities into
user-facing actions and a recommended continuation workflow.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from .models import NMRDataset


@dataclass(frozen=True)
class WorkflowAction:
    key: str
    title: str
    available: bool
    reason: str


_WORKFLOWS = (
    ('prepare', 'Prepare raw data'),
    ('inspect', 'Inspect projections'),
    ('decon', 'Detect peaks / deconvolve'),
    ('slices', 'Explore slices'),
    ('special', 'Specialized analyses'),
)


def recommended_workflow(dataset: NMRDataset) -> Optional[str]:
    """Return the most useful next Decon workflow, or None for project summary."""
    c = dataset.capabilities
    if dataset.project is None or c is None:
        return None
    if c.can_view_peaks:
        return 'inspect'
    if c.can_deconvolve:
        return 'decon'
    if c.can_process_raw:
        return 'prepare'
    return None


def workflow_actions(dataset: NMRDataset) -> tuple[WorkflowAction, ...]:
    """Describe every Decon workflow, including why unavailable ones are disabled."""
    c = dataset.capabilities
    if dataset.project is None or c is None:
        return ()

    spectrum_reason = 'Requires an available processed spectrum.'
    raw_reason = 'Requires the configured raw acquisition to be available.'
    actions = []
    for key, title in _WORKFLOWS:
        if key == 'prepare':
            ok, reason = c.can_process_raw, ('Raw acquisition is available.' if c.can_process_raw else raw_reason)
        elif key == 'decon':
            ok, reason = c.can_deconvolve, ('Processed spectrum is available.' if c.can_deconvolve else spectrum_reason)
        elif key in ('inspect', 'slices', 'special'):
            ok, reason = c.can_view_spectrum, ('Processed spectrum is available.' if c.can_view_spectrum else spectrum_reason)
        else:
            ok, reason = False, 'Unavailable.'
        actions.append(WorkflowAction(key, title, ok, reason))
    return tuple(actions)


def workflow_title(key: Optional[str]) -> str:
    if key is None:
        return 'Project summary'
    return dict(_WORKFLOWS).get(key, key)
