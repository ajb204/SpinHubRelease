from types import SimpleNamespace

from spinDecon.domain.analysis_mode import AnalysisMode
from spinDecon.project.data_store import DataStore
from spinDecon.project.summary import _workflow_report_data


def _frame_for_2d_pseudo(store):
    state = SimpleNamespace(dimension=2, pseudo_axis=True)
    # AnalysisMode accepts the normal ProjectState attributes; supply aliases
    # used by older versions defensively.
    state.spectral_dimensions = 2
    parent = SimpleNamespace(data_store=store, tabOne=SimpleNamespace(
        pseudo_intensities_stale=False,
        downstream_analysis='Decay',
    ))
    parent.get_page_by_title = lambda title: None
    return SimpleNamespace(state=state, store=store, parent=parent)


def test_summary_marks_review_complete_when_terminal_2d_pseudo_workflow_complete(monkeypatch):
    store = DataStore()
    store.analysis['pseudo_intensities_ready'] = True
    store.analysis['downstream_analysis'] = 'Decay'
    frame = _frame_for_2d_pseudo(store)

    # Avoid coupling this regression test to ProjectState construction details.
    monkeypatch.setattr(AnalysisMode, 'from_project_state', classmethod(
        lambda cls, state: cls.from_legacy(2, True)
    ))

    report = _workflow_report_data(frame, [])
    rows = {title: (complete, status, detail) for title, complete, status, detail in report['rows']}
    review = next(value for title, value in rows.items() if 'review' in title.lower() and 'intensity' in title.lower())
    analyse = next(value for title, value in rows.items() if 'analyse' in title.lower() and 'intensity' in title.lower())
    assert analyse[0] is True
    assert review[0] is True
    assert review[1] == 'complete'
