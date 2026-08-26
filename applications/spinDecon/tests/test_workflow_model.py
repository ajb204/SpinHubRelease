import unittest

from spinDecon.domain.analysis_mode import AnalysisMode
from spinDecon.workflow.model import StageRequirement, build_workflow_plan


class WorkflowModelTests(unittest.TestCase):
    def test_spectral_plan_targets_peak_list(self):
        plan = build_workflow_plan(AnalysisMode.from_legacy(3, False))
        keys = [stage.key for stage in plan.stages]
        self.assertEqual(keys, [
            "spectrum", "peak_shape", "reference_peaks",
            "peak_pick", "review_peaks",
        ])
        self.assertEqual(plan.stage("reference_peaks").requirement, StageRequirement.REQUIRED)
        self.assertNotIn("fit_spectrum", keys)

    def test_pseudo_plan_targets_intensity_series(self):
        plan = build_workflow_plan(AnalysisMode.from_legacy(2, True))
        keys = [stage.key for stage in plan.stages]
        self.assertEqual(keys, [
            "spectrum", "peak_shape", "reference_peaks",
            "extract_intensities", "review_series", "analyse_series",
        ])
        self.assertNotIn("peak_pick", keys)
        self.assertIn("pseudo-axis", plan.objective)


if __name__ == "__main__":
    unittest.main()


def test_stage_dependencies_are_explicit():
    spectral = build_workflow_plan(AnalysisMode.from_legacy(3, False))
    assert spectral.stage("peak_shape").requires == ("spectrum",)
    assert spectral.stage("peak_pick").requires == ("spectrum", "reference_peaks")
    assert spectral.stage("review_peaks").requires == ("peak_pick",)

    pseudo = build_workflow_plan(AnalysisMode.from_legacy(2, True))
    assert pseudo.stage("reference_peaks").requires == ("peak_shape",)
    assert pseudo.stage("extract_intensities").requires == ("peak_shape", "reference_peaks")
    assert pseudo.stage("review_series").requires == ("extract_intensities",)
    # Analysis is scientifically available after extraction; review remains the
    # preferred recommendation rather than a hard gate.
    assert pseudo.stage("analyse_series").requires == ("extract_intensities",)
