"""Pure workflow definitions for the guided NMR analysis interface.

There are no wx imports and no scientific operations here.  The model only
states which scientific stages apply to an AnalysisMode.  This keeps milestone
1 inert with respect to the existing GUI and analysis engine.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from ..domain.analysis_mode import AnalysisMode, WorkflowKind
from ..domain.spectrum_policy import ReferencePolicy, spectrum_policy


class StageRequirement(str, Enum):
    REQUIRED = "required"
    CONDITIONAL = "conditional"
    OPTIONAL = "optional"


@dataclass(frozen=True)
class WorkflowStage:
    key: str
    title: str
    description: str
    requirement: StageRequirement = StageRequirement.REQUIRED
    requires: tuple[str, ...] = ()
    recommendation_rank: int = 100


@dataclass(frozen=True)
class WorkflowPlan:
    mode: AnalysisMode
    stages: tuple[WorkflowStage, ...]

    @property
    def objective(self) -> str:
        if self.mode.workflow_kind is WorkflowKind.PSEUDO_AXIS_SERIES:
            return "Extract peak intensities across the pseudo-axis for downstream analysis."
        return "Produce a peak list for the frequency-domain spectrum."

    def stage(self, key: str) -> WorkflowStage:
        for item in self.stages:
            if item.key == key:
                return item
        raise KeyError(key)


_SPECTRUM = WorkflowStage(
    "spectrum", "Prepare spectrum",
    "Convert vendor time-domain data to a frequency-domain spectrum, or load an existing processed spectrum.",
    recommendation_rank=10,
)
_PEAK_SHAPE = WorkflowStage(
    "peak_shape", "Determine peak shape",
    "Determine the peak-shape parameters used by subsequent peak fitting.",
    requires=("spectrum",), recommendation_rank=20,
)
_REFERENCE = WorkflowStage(
    "reference_peaks", "Establish reference peaks",
    "Load or create the reference frequencies/peak list used to constrain the analysis.",
    StageRequirement.CONDITIONAL,
    requires=("spectrum",), recommendation_rank=30,
)
_PICK = WorkflowStage(
    "peak_pick", "Pick peaks",
    "Run the main peak-picking calculation for the spectral dataset.",
    requires=("spectrum",), recommendation_rank=40,
)
_FIT = WorkflowStage(
    "fit_spectrum", "Fit spectrum",
    "Fit selected peaks to the peak-shape model when a fitted spectrum is required.",
    StageRequirement.OPTIONAL,
    requires=("peak_pick",), recommendation_rank=60,
)
_REVIEW_FIT = WorkflowStage(
    "review_fitting", "Review fitting results",
    "Inspect the fitted 2D peaks in the Fitting workspace and fitting-results window.",
    requires=("fit_spectrum",), recommendation_rank=70,
)
_REVIEW_PEAKS = WorkflowStage(
    "review_peaks", "Review picked peaks",
    "Open the full-dimensional peak list and the 2D slice viewer together to assess the picked peaks against the spectrum.",
    requires=("peak_pick",), recommendation_rank=50,
)
_EXTRACT = WorkflowStage(
    "extract_intensities", "Extract pseudo-axis intensities",
    "Fit or extract each selected peak through the pseudo-axis stack to obtain an intensity series.",
    requires=("peak_shape", "reference_peaks"), recommendation_rank=40,
)
_REVIEW_SERIES = WorkflowStage(
    "review_series", "Review intensity series",
    "Inspect peak intensity versus pseudo-axis slice before downstream analysis.",
    requires=("extract_intensities",), recommendation_rank=50,
)
_ANALYSE_SERIES = WorkflowStage(
    "analyse_series", "Analyse intensity series",
    "Pass intensity series to experiment-specific analysis such as diffusion or relaxation fitting.",
    requires=("extract_intensities",), recommendation_rank=60,
)


def build_workflow_plan(mode: AnalysisMode) -> WorkflowPlan:
    """Return the scientific workflow implied by ``mode``.

    The plan is descriptive only.  It intentionally contains no GUI callbacks,
    file tests, or status flags; those belong to later implementation phases.
    """
    if mode.workflow_kind is WorkflowKind.PSEUDO_AXIS_SERIES:
        # A pseudo-axis analysis is anchored to selected spectral frequencies.
        # Make that dependency explicit even though the legacy specialist tabs
        # may also allow manual region selection.
        reference = WorkflowStage(_REFERENCE.key, _REFERENCE.title, _REFERENCE.description, StageRequirement.REQUIRED, ("peak_shape",), _REFERENCE.recommendation_rank)
        stages = (
            _SPECTRUM, _PEAK_SHAPE, reference,
            _EXTRACT, _REVIEW_SERIES,
        ) if mode.spectral_dimensions == 1 else (
            _SPECTRUM, _PEAK_SHAPE, reference,
            _EXTRACT, _REVIEW_SERIES, _ANALYSE_SERIES,
        )
    else:
        # Peak-list ownership comes from the canonical spectrum policy.  In
        # physical 2D, Full == reference conceptually, so there is no separate
        # reference-list gate.  Higher-dimensional spectra use an independent
        # lower-dimensional reference projection.
        policy = spectrum_policy(mode)
        if policy.reference_policy is ReferencePolicy.INDEPENDENT_PROJECTION:
            reference = WorkflowStage(_REFERENCE.key, _REFERENCE.title, _REFERENCE.description, StageRequirement.REQUIRED, _REFERENCE.requires, _REFERENCE.recommendation_rank)
        else:
            reference = WorkflowStage(_REFERENCE.key, _REFERENCE.title,
                "No independent reference list is required; the Full Peak List supplies the peak identities for this spectrum.",
                StageRequirement.OPTIONAL, _REFERENCE.requires, _REFERENCE.recommendation_rank)
        pick = WorkflowStage(_PICK.key, _PICK.title, _PICK.description, reference.requirement, _PICK.requires + (("reference_peaks",) if reference.requirement is StageRequirement.REQUIRED else ()), _PICK.recommendation_rank)
        core = (_SPECTRUM, _PEAK_SHAPE, reference, _PICK if reference.requirement is not StageRequirement.REQUIRED else pick)
        # True 3D spectral analysis ends with review of the generated peak list
        # and inspection in 2D slices; there is no separate whole-spectrum fit.
        if mode.spectral_dimensions == 3:
            stages = core + (_REVIEW_PEAKS,)
        elif mode.spectral_dimensions == 2:
            # Mirror the pseudo3D structure: pick -> explicitly check peaks ->
            # restrained fit -> explicitly inspect fitting results.
            fit_2d = WorkflowStage(
                _FIT.key, _FIT.title,
                "Fit the checked 2D peak list to the peak-shape model using restrained reconstruction.",
                StageRequirement.REQUIRED, ("review_peaks",), _FIT.recommendation_rank)
            stages = core + (_REVIEW_PEAKS, fit_2d, _REVIEW_FIT)
        else:
            stages = core + (_FIT, _REVIEW_PEAKS)
    return WorkflowPlan(mode=mode, stages=stages)
