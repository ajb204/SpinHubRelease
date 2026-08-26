from spinDecon.domain.analysis_mode import AnalysisMode
from spinDecon.domain.spectrum_policy import ReferencePolicy, spectrum_policy


def test_canonical_journey_peak_ownership():
    cases = [
        (AnalysisMode.from_legacy(1, False), "1D", ReferencePolicy.NONE, 1, None),
        (AnalysisMode.from_legacy(1, True), "pseudo2D", ReferencePolicy.FULL_IS_REFERENCE, 1, 1),
        (AnalysisMode.from_legacy(2, False), "2D", ReferencePolicy.FULL_IS_REFERENCE, 2, 2),
        (AnalysisMode.from_legacy(2, True), "pseudo3D", ReferencePolicy.INDEPENDENT_PROJECTION, 2, 2),
        (AnalysisMode.from_legacy(3, False), "3D", ReferencePolicy.INDEPENDENT_PROJECTION, 3, 2),
    ]
    for mode, journey, reference_policy, full_dim, reference_dim in cases:
        policy = spectrum_policy(mode)
        assert policy.journey == journey
        assert policy.reference_policy is reference_policy
        assert policy.full_dimension == full_dim
        assert policy.reference_dimension == reference_dim


def test_physical_2d_has_no_independent_reference_authority():
    policy = spectrum_policy(AnalysisMode.from_legacy(2, False))
    assert policy.full_is_reference
    assert not policy.has_independent_reference


def test_physical_2d_exposes_only_full_peak_authority():
    policy = spectrum_policy(AnalysisMode.from_legacy(2, False))
    assert policy.full_is_reference
    assert not policy.has_distinct_reference_peak_list
    assert policy.projection_peak_list_key == 'full'
    assert policy.fitting_peak_list_key == 'full'


def test_pseudo3d_keeps_independent_reference_authority():
    policy = spectrum_policy(AnalysisMode.from_legacy(2, True))
    assert policy.has_distinct_reference_peak_list
    assert policy.projection_peak_list_key == 'reference'
    assert policy.fitting_peak_list_key == 'reference'
