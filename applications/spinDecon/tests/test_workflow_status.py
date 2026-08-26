import tempfile, unittest
from pathlib import Path
from spinDecon.domain.analysis_mode import AnalysisMode
from spinDecon.project.data_store import DataStore
from spinDecon.project.state import ProjectState
from spinDecon.workflow.model import StageRequirement, build_workflow_plan
from spinDecon.workflow.status import StageStatus, evaluate_workflow, next_action
class WorkflowStatusTests(unittest.TestCase):
 def test_requirements(self):
  req=lambda d: build_workflow_plan(AnalysisMode.from_legacy(d)).stage('reference_peaks').requirement
  self.assertEqual([req(i) for i in range(1,5)],[StageRequirement.OPTIONAL,StageRequirement.OPTIONAL,StageRequirement.REQUIRED,StageRequirement.REQUIRED])
 def test_empty(self):
  s=ProjectState(dimension=2); st=evaluate_workflow(build_workflow_plan(AnalysisMode.from_project_state(s)),s,DataStore()); d={x.key:x for x in st}
  self.assertEqual(d['spectrum'].status,StageStatus.READY); self.assertEqual(d['peak_shape'].status,StageStatus.BLOCKED); self.assertEqual(next_action(st),'spectrum')
 def test_spectrum_and_shape(self):
  s=ProjectState(dimension=2); store=DataStore(); store.save_spectrum('raw',data=[1]); plan=build_workflow_plan(AnalysisMode.from_project_state(s)); d={x.key:x for x in evaluate_workflow(plan,s,store)}
  self.assertEqual(d['spectrum'].status,StageStatus.COMPLETE); self.assertEqual(d['peak_shape'].status,StageStatus.READY)
  store.metadata['peak_shape_determined']=True; d={x.key:x for x in evaluate_workflow(plan,s,store)}; self.assertEqual(d['peak_shape'].status,StageStatus.COMPLETE)
 def test_3d_reference_gate(self):
  s=ProjectState(dimension=3); store=DataStore(); store.data=[1]; plan=build_workflow_plan(AnalysisMode.from_project_state(s)); d={x.key:x for x in evaluate_workflow(plan,s,store)}; self.assertEqual(d['peak_pick'].status,StageStatus.BLOCKED)
  store.save_peak_list('reference',peaks=[object()]); d={x.key:x for x in evaluate_workflow(plan,s,store)}; self.assertEqual(d['peak_pick'].status,StageStatus.READY)
 def test_full_list(self):
  s=ProjectState(dimension=3); store=DataStore(); store.data=[1]; store.save_peak_list('full',peaks=[object()]); plan=build_workflow_plan(AnalysisMode.from_project_state(s)); d={x.key:x for x in evaluate_workflow(plan,s,store)}
  self.assertEqual(d['peak_pick'].status,StageStatus.COMPLETE); self.assertEqual(d['review_peaks'].status,StageStatus.READY); self.assertNotIn('fit_spectrum',d)
 def test_files(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'spec'; p.mkdir(); (p/'x.ft2').write_text('x'); (p/'ref.list').write_text('x'); s=ProjectState(working_dir=td,spec_path='spec',input_file='x.ft2',reference_peak_file='ref.list',dimension=2); d={x.key:x for x in evaluate_workflow(build_workflow_plan(AnalysisMode.from_project_state(s)),s,DataStore())}; self.assertEqual(d['spectrum'].status,StageStatus.COMPLETE); self.assertEqual(d['reference_peaks'].status,StageStatus.OPTIONAL)
if __name__=='__main__': unittest.main()


def test_pseudo_workflow_blocks_extraction_until_reference(tmp_path):
    from spinDecon.domain.analysis_mode import AnalysisMode
    from spinDecon.workflow.model import build_workflow_plan
    from spinDecon.workflow.status import evaluate_workflow, StageStatus
    from spinDecon.project.state import ProjectState
    from spinDecon.project.data_store import DataStore
    state = ProjectState(working_dir=str(tmp_path), dimension=3, pseudo_axis=True)
    store = DataStore()
    store.save_spectrum('raw', data=[1])
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert states['extract_intensities'].status is StageStatus.BLOCKED
    store.save_peak_list('reference', peaks=[object()])
    store.mark_peak_shape_determined(source='test')
    states = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert states['extract_intensities'].status is StageStatus.READY


def test_pseudo_workflow_recognises_intensity_evidence(tmp_path):
    from spinDecon.domain.analysis_mode import AnalysisMode
    from spinDecon.workflow.model import build_workflow_plan
    from spinDecon.workflow.status import evaluate_workflow, StageStatus
    from spinDecon.project.state import ProjectState
    from spinDecon.project.data_store import DataStore
    state = ProjectState(working_dir=str(tmp_path), dimension=3, pseudo_axis=True)
    store = DataStore()
    store.save_spectrum('raw', data=[1])
    store.save_peak_list('reference', peaks=[object()])
    store.analysis['pseudo_intensities_ready'] = True
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert states['extract_intensities'].status is StageStatus.COMPLETE
    assert states['review_series'].status is StageStatus.READY
    assert states['analyse_series'].status is StageStatus.READY


def test_datastore_workflow_evidence_api():
    store = DataStore()
    store.mark_peak_shape_determined(model="voigt", source="fit")
    store.mark_pseudo_intensities_ready(count=12)
    store.mark_pseudo_series_reviewed(source="review")
    store.mark_pseudo_analysis_complete(model="relaxation")

    assert store.metadata["peak_shape_determined"] is True
    assert store.metadata["peak_shape"] == {"model": "voigt", "source": "fit"}
    assert store.analysis["pseudo_intensities_ready"] is True
    assert store.analysis["pseudo_intensities"] == {"count": 12}
    assert store.analysis["pseudo_series_reviewed"] is True
    assert store.analysis["pseudo_series_review"] == {"source": "review"}
    assert store.analysis["pseudo_analysis_complete"] is True
    assert store.analysis["pseudo_analysis"] == {"model": "relaxation"}


def test_workflow_evidence_api_drives_status(tmp_path):
    state = ProjectState(working_dir=str(tmp_path), dimension=3, pseudo_axis=True)
    store = DataStore()
    store.save_spectrum("raw", data=[1])
    store.save_peak_list("reference", peaks=[object()])
    store.mark_peak_shape_determined(source="test")
    store.mark_pseudo_intensities_ready(source="test")
    store.mark_pseudo_series_reviewed(source="test")
    store.mark_pseudo_analysis_complete(source="test")

    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert states["peak_shape"].status is StageStatus.COMPLETE
    assert states["extract_intensities"].status is StageStatus.COMPLETE
    assert states["review_series"].status is StageStatus.COMPLETE
    assert states["analyse_series"].status is StageStatus.COMPLETE


def test_3d_spectral_journey_reference_then_full_list(tmp_path):
    state = ProjectState(working_dir=str(tmp_path), dimension=3)
    store = DataStore()
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))

    store.save_spectrum("raw", data=[1])
    store.mark_peak_shape_determined(source="journey")
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key["reference_peaks"].status is StageStatus.READY
    assert by_key["peak_pick"].status is StageStatus.BLOCKED
    assert next_action(states) == "reference_peaks"

    store.save_peak_list("reference", peaks=[object()], dimension=2)
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key["reference_peaks"].status is StageStatus.COMPLETE
    assert by_key["peak_pick"].status is StageStatus.READY
    assert next_action(states) == "peak_pick"

    store.save_peak_list("full", peaks=[object()], dimension=3)
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key["peak_pick"].status is StageStatus.COMPLETE
    assert by_key["review_peaks"].status is StageStatus.READY
    assert next_action(states) == "review_peaks"


def test_pseudo_journey_extraction_review_analysis(tmp_path):
    state = ProjectState(working_dir=str(tmp_path), dimension=3, pseudo_axis=True)
    store = DataStore()
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))

    store.save_spectrum('raw', data=[1])
    store.mark_peak_shape_determined(source='journey')
    store.save_peak_list('reference', peaks=[object()], dimension=2)
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key['extract_intensities'].status is StageStatus.READY
    assert next_action(states) == 'extract_intensities'

    store.mark_pseudo_intensities_ready(source='fuda')
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key['extract_intensities'].status is StageStatus.COMPLETE
    assert by_key['review_series'].status is StageStatus.READY
    assert next_action(states) == 'review_series'

    store.mark_pseudo_series_reviewed(source='fuda_overview')
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key['review_series'].status is StageStatus.COMPLETE
    assert by_key['analyse_series'].status is StageStatus.READY
    assert next_action(states) == 'analyse_series'

    store.mark_pseudo_analysis_complete(model='decay', source='catia')
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key['analyse_series'].status is StageStatus.COMPLETE
    assert next_action(states) is None


def test_recommendation_is_independent_of_state_tuple_order(tmp_path):
    from spinDecon.workflow.status import available_actions, recommended_action
    state = ProjectState(working_dir=str(tmp_path), dimension=3, pseudo_axis=True)
    store = DataStore()
    store.save_spectrum("raw", data=[1])
    store.save_peak_list("reference", peaks=[object()])
    store.mark_pseudo_intensities_ready(source="test")
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = evaluate_workflow(plan, state, store)

    assert set(available_actions(states)) == {"peak_shape", "review_series", "analyse_series"}
    # Peak shape is unfinished and therefore remains the preferred scientific
    # next step even though later pseudo-series actions are technically usable.
    assert recommended_action(plan, tuple(reversed(states))) == "peak_shape"

    store.mark_peak_shape_determined(source="test")
    states = evaluate_workflow(plan, state, store)
    assert set(available_actions(states)) == {"review_series", "analyse_series"}
    assert recommended_action(plan, tuple(reversed(states))) == "review_series"


def test_existing_pseudo_fit_evidence_preserves_persisted_review():
    store = DataStore()
    store.mark_pseudo_series_reviewed(source="system_file")
    store.mark_pseudo_intensities_ready(
        source="existing_protocol3p_fit_cold_start",
        invalidate_review=False,
    )
    assert store.analysis["pseudo_series_reviewed"] is True


def test_new_pseudo_fit_invalidates_previous_review():
    store = DataStore()
    store.mark_pseudo_series_reviewed(source="previous_session")
    store.mark_pseudo_intensities_ready(source="protocol3p_recon")
    assert "pseudo_series_reviewed" not in store.analysis


def test_2d_pseudo_analysis_selection_completes_terminal_workflow(tmp_path):
    state = ProjectState(working_dir=str(tmp_path), dimension=2, pseudo_axis=True)
    store = DataStore()
    store.save_spectrum("raw", data=[1])
    store.save_peak_list("reference", peaks=[object()], dimension=2)
    store.mark_peak_shape_determined(source="test")
    store.mark_pseudo_intensities_ready(source="test")
    store.mark_pseudo_series_reviewed(source="test")
    store.analysis["downstream_analysis"] = "Decay"

    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = evaluate_workflow(plan, state, store)
    by_key = {item.key: item for item in states}
    assert by_key["analyse_series"].status is StageStatus.COMPLETE
    assert "Decay" in by_key["analyse_series"].detail
    assert next_action(states) is None


def test_pseudo2d_extract_requires_fit_pair_for_every_full_1d_peak(tmp_path):
    from spinDecon.workflow.status import _pseudo2d_fit_files_complete
    state = ProjectState(working_dir=str(tmp_path), dimension=2, pseudo_axis=True)
    store = DataStore()
    store.save_peak_list('full', rows=[['10', '8.10', '1.0'], ['20', '7.20', '1.0']], dimension=1)
    fit = tmp_path / 'fit'
    fit.mkdir()

    class Tab:
        def get_fuda_dir(self): return str(fit)
    class Context:
        tabOne = Tab()

    (fit / '10.dat').write_text('x')
    (fit / '10.out').write_text('x')
    (fit / '20.dat').write_text('x')
    assert not _pseudo2d_fit_files_complete(state, store, Context())

    (fit / '20.out').write_text('x')
    assert _pseudo2d_fit_files_complete(state, store, Context())

def test_loaded_spectrum_without_backing_file_does_not_complete_prepare_spectrum(tmp_path):
    state = ProjectState(working_dir=str(tmp_path), spec_path='spec', input_file='missing.ft2', dimension=2)
    store = DataStore()
    store.save_spectrum('raw', data=[1])
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = {item.key: item for item in evaluate_workflow(plan, state, store)}
    assert states['spectrum'].status is StageStatus.READY


def test_backing_spectrum_file_completes_prepare_spectrum_without_loading(tmp_path):
    spec = tmp_path / 'spec'
    spec.mkdir()
    (spec / 'made.ft2').write_text('spectrum')
    state = ProjectState(working_dir=str(tmp_path), spec_path='spec', input_file='made.ft2', dimension=2)
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))
    states = {item.key: item for item in evaluate_workflow(plan, state, DataStore())}
    assert states['spectrum'].status is StageStatus.COMPLETE


def test_peak_pick_stale_and_review_checked_are_workflow_evidence(tmp_path):
    spectrum = tmp_path / "spectrum.ft3"
    spectrum.write_bytes(b"test")
    state = ProjectState(working_dir=str(tmp_path), spec_path=str(tmp_path), input_file="spectrum.ft3", dimension=3)
    store = DataStore()
    store.save_spectrum('raw', data=[1])
    store.mark_peak_shape_determined(source='test')
    store.save_peak_list('reference', peaks=[object()], dimension=2)
    store.save_peak_list('full', peaks=[object()], dimension=3)
    plan = build_workflow_plan(AnalysisMode.from_project_state(state))

    by_key = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert by_key['peak_pick'].status is StageStatus.COMPLETE
    assert by_key['review_peaks'].status is StageStatus.READY

    store.mark_picked_peaks_reviewed(source='test')
    by_key = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert by_key['review_peaks'].status is StageStatus.COMPLETE

    store.analysis['peak_pick_stale'] = True
    store.invalidate_picked_peaks_review()
    by_key = {x.key: x for x in evaluate_workflow(plan, state, store)}
    assert by_key['peak_pick'].status is StageStatus.READY
    assert by_key['review_peaks'].status is StageStatus.BLOCKED
