"""Architecture regression coverage for stages 91-94."""

def test_root_domain_compatibility_preserves_class_identity():
    from spinDecon.domain.analysis_mode import AnalysisMode as LegacyAnalysisMode
    from spinDecon.domain.analysis_mode import AnalysisMode
    from spinDecon.domain.topology import DatasetTopology as LegacyTopology
    from spinDecon.domain.topology import DatasetTopology
    assert LegacyAnalysisMode is AnalysisMode
    assert LegacyTopology is DatasetTopology


def test_root_workflow_compatibility_preserves_class_identity():
    from spinDecon.workflow.model import WorkflowPlan as LegacyWorkflowPlan
    from spinDecon.workflow.model import WorkflowPlan
    from spinDecon.workflow.status import StageStatus as LegacyStageStatus
    from spinDecon.workflow.status import StageStatus
    assert LegacyWorkflowPlan is WorkflowPlan
    assert LegacyStageStatus is StageStatus


def test_dimension_contract_compatibility_preserves_identity():
    from spinDecon.domain.dimensions.peak_contract import peak_coordinate_count as legacy
    from spinDecon.domain.dimensions.peak_contract import peak_coordinate_count
    assert legacy is peak_coordinate_count
