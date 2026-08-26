import unittest

from spinDecon.domain.analysis_mode import AnalysisMode, WorkflowKind


class AnalysisModeTests(unittest.TestCase):
    def test_normal_dimensions_are_all_spectral(self):
        for dim in range(1, 5):
            mode = AnalysisMode.from_legacy(dim, False)
            self.assertEqual(mode.physical_dimensions, dim)
            self.assertEqual(mode.spectral_dimensions, dim)
            self.assertEqual(mode.pseudo_dimensions, 0)
            self.assertEqual(mode.workflow_kind, WorkflowKind.SPECTRAL_PEAK_LIST)

    def test_legacy_pseudo_adds_one_physical_axis(self):
        mode = AnalysisMode.from_legacy(2, True)
        self.assertEqual(mode.legacy_dimension, 2)
        self.assertEqual(mode.spectral_dimensions, 2)
        self.assertEqual(mode.physical_dimensions, 3)
        self.assertEqual(mode.pseudo_dimensions, 1)
        self.assertEqual(mode.workflow_kind, WorkflowKind.PSEUDO_AXIS_SERIES)

    def test_dimension_always_means_spectral_dimensions(self):
        mode = AnalysisMode.from_legacy(3, True)
        self.assertEqual(mode.spectral_dimensions, 3)
        self.assertEqual(mode.physical_dimensions, 4)
        self.assertEqual(mode.pseudo_dimensions, 1)

    def test_pseudo4d_matches_existing_gui_restriction(self):
        with self.assertRaises(ValueError):
            AnalysisMode.from_legacy(4, True)

    def test_invalid_dimension_is_rejected(self):
        for dim in (0, 5, -1):
            with self.assertRaises(ValueError):
                AnalysisMode.from_legacy(dim)


if __name__ == "__main__":
    unittest.main()
