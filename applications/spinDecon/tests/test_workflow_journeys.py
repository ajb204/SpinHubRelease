"""Behavioural specifications for the five supported guided-workflow journeys.

These tests deliberately exercise only AnalysisMode + WorkflowPlan + DataStore
+ workflow_status.  They are intended to remain independent of wxPython so the
scientific workflow can be changed safely before the GUI is simplified.
"""
from __future__ import annotations

import pytest

from spinDecon.domain.analysis_mode import AnalysisMode, WorkflowKind
from spinDecon.project.data_store import DataStore
from spinDecon.project.state import ProjectState
from spinDecon.workflow.model import StageRequirement, build_workflow_plan
from spinDecon.workflow.status import (
    StageStatus,
    available_actions,
    evaluate_workflow,
    recommended_action,
)


def _workflow(dimension: int, *, pseudo: bool = False):
    state = ProjectState(dimension=dimension, pseudo_axis=pseudo)
    store = DataStore()
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    return state, store, plan


def _states(plan, state, store):
    states = evaluate_workflow(plan, state, store)
    return states, {item.key: item.status for item in states}


def _recommendation(plan, state, store):
    states, by_key = _states(plan, state, store)
    return recommended_action(plan, states), by_key, set(available_actions(states))


def _prepare_spectrum_and_shape(store):
    store.save_spectrum("raw", data=[1])
    store.mark_peak_shape_determined(source="journey-test")


def test_1d_spectral_journey_reference_is_optional():
    state, store, plan = _workflow(1)
    assert plan.mode.workflow_kind is WorkflowKind.SPECTRAL_PEAK_LIST
    assert plan.stage("reference_peaks").requirement is StageRequirement.OPTIONAL

    action, status, _ = _recommendation(plan, state, store)
    assert action == "spectrum"
    assert status["peak_shape"] is StageStatus.BLOCKED

    store.save_spectrum("raw", data=[1])
    action, status, _ = _recommendation(plan, state, store)
    assert action == "peak_shape"
    assert status["reference_peaks"] is StageStatus.OPTIONAL
    assert status["peak_pick"] is StageStatus.READY

    store.mark_peak_shape_determined(source="journey-test")
    action, status, available = _recommendation(plan, state, store)
    assert action == "peak_pick"
    assert "reference_peaks" not in available

    store.save_peak_list("full", peaks=[object()], dimension=1)
    action, status, _ = _recommendation(plan, state, store)
    assert status["peak_pick"] is StageStatus.COMPLETE
    assert status["review_peaks"] is StageStatus.READY
    assert status["fit_spectrum"] is StageStatus.OPTIONAL
    assert action == "review_peaks"


def test_2d_spectral_journey_full_list_is_the_reference_authority():
    state, store, plan = _workflow(2)
    assert plan.stage("reference_peaks").requirement is StageRequirement.OPTIONAL

    _prepare_spectrum_and_shape(store)
    action, status, available = _recommendation(plan, state, store)
    assert status["reference_peaks"] is StageStatus.OPTIONAL
    assert status["peak_pick"] is StageStatus.READY
    assert action == "peak_pick"

    # Physical 2D is the singularity: Full is also the reference authority.
    store.save_peak_list("full", peaks=[object()], dimension=2)
    action, status, _ = _recommendation(plan, state, store)
    assert status["peak_pick"] is StageStatus.COMPLETE
    assert status["reference_peaks"] is StageStatus.COMPLETE
    assert status["review_peaks"] is StageStatus.READY
    assert action == "review_peaks"


def test_3d_spectral_journey_reference_is_a_hard_gate():
    state, store, plan = _workflow(3)
    assert plan.stage("reference_peaks").requirement is StageRequirement.REQUIRED

    _prepare_spectrum_and_shape(store)
    action, status, available = _recommendation(plan, state, store)
    assert status["reference_peaks"] is StageStatus.READY
    assert status["peak_pick"] is StageStatus.BLOCKED
    assert "peak_pick" not in available
    assert action == "reference_peaks"

    store.save_peak_list("reference", peaks=[object()], dimension=2)
    action, status, _ = _recommendation(plan, state, store)
    assert status["peak_pick"] is StageStatus.READY
    assert action == "peak_pick"

    store.save_peak_list("full", peaks=[object()], dimension=3)
    action, status, _ = _recommendation(plan, state, store)
    assert status["peak_pick"] is StageStatus.COMPLETE
    assert status["review_peaks"] is StageStatus.READY
    assert "fit_spectrum" not in status
    assert action == "review_peaks"


@pytest.mark.parametrize(
    "dimension, expected_spectral_dimensions",
    [(1, 1), (2, 2)],
    ids=["1d+pseudo", "2d+pseudo"],
)
def test_pseudo_series_journey(dimension, expected_spectral_dimensions):
    state, store, plan = _workflow(dimension, pseudo=True)
    assert plan.mode.workflow_kind is WorkflowKind.PSEUDO_AXIS_SERIES
    assert plan.mode.spectral_dimensions == expected_spectral_dimensions
    assert plan.stage("reference_peaks").requirement is StageRequirement.REQUIRED

    action, status, _ = _recommendation(plan, state, store)
    assert action == "spectrum"
    assert status["extract_intensities"] is StageStatus.BLOCKED

    store.save_spectrum("raw", data=[1])
    action, status, _ = _recommendation(plan, state, store)
    assert action == "peak_shape"
    assert status["reference_peaks"] is StageStatus.BLOCKED
    assert status["extract_intensities"] is StageStatus.BLOCKED

    store.mark_peak_shape_determined(source="journey-test")
    action, status, _ = _recommendation(plan, state, store)
    assert action == "reference_peaks"

    # Pseudo2D (1 spectral + pseudo) uses Full 1D as its authoritative
    # reference-frequency set; higher-dimensional pseudo workflows use the
    # separate reference list.
    peak_key = "full" if dimension == 1 else "reference"
    store.save_peak_list(peak_key, peaks=[object()], dimension=expected_spectral_dimensions)
    action, status, _ = _recommendation(plan, state, store)
    assert status["extract_intensities"] is StageStatus.READY
    assert action == "extract_intensities"

    store.mark_pseudo_intensities_ready(source="journey-test")
    action, status, available = _recommendation(plan, state, store)
    assert status["extract_intensities"] is StageStatus.COMPLETE
    assert "review_series" in available
    if dimension == 2:
        assert "analyse_series" in available
    else:
        assert "analyse_series" not in status
    assert action == "review_series"

    store.mark_pseudo_series_reviewed(source="journey-test")
    action, status, _ = _recommendation(plan, state, store)
    assert status["review_series"] is StageStatus.COMPLETE
    if dimension == 2:
        assert action == "analyse_series"
        store.mark_pseudo_analysis_complete(source="journey-test", model="decay")
        action, status, available = _recommendation(plan, state, store)
        assert status["analyse_series"] is StageStatus.COMPLETE
        assert action is None
        assert not available
    else:
        # Pseudo2D ends after review; downstream analysis is handled by its
        # specialist workspace rather than a separate guided-workflow stage.
        assert action is None
