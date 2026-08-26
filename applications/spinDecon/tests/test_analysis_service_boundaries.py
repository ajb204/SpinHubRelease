"""Architecture regression tests for application-service extraction."""
import ast
from pathlib import Path


def test_analysis_services_do_not_import_wx():
    root = Path(__file__).resolve().parents[1] / "analysis"
    service_files = list(root.glob("*_service.py"))
    assert service_files
    offenders = []
    for path in service_files:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import) and any(a.name == "wx" or a.name.startswith("wx.") for a in node.names):
                offenders.append(path.name)
            if isinstance(node, ast.ImportFrom) and node.module and (node.module == "wx" or node.module.startswith("wx.")):
                offenders.append(path.name)
    assert not offenders, sorted(set(offenders))


def test_application_context_exposes_migration_services():
    from spinDecon.app.context import ApplicationContext
    fields = ApplicationContext.__dataclass_fields__
    for name in ("full3d", "one_d", "projection", "peaks", "slices", "pseudo"):
        assert name in fields


def test_peak_fit_service_is_gui_independent():
    from pathlib import Path
    source = (Path(__file__).parents[1] / 'analysis' / 'peak_fit_service.py').read_text()
    assert 'import wx' not in source
    assert 'from wx' not in source


def test_application_context_exposes_peak_fit_boundary():
    from spinDecon.app.context import ApplicationContext
    assert 'peak_fit' in ApplicationContext.__dataclass_fields__

class _ValueControl:
    def __init__(self, value): self.value = value
    def GetValue(self): return str(self.value)


def test_peak_fit_shape_parameters_are_numeric_and_dimension_bounded():
    from spinDecon.analysis.peak_fit_service import PeakFitService
    class Legacy:
        sig1Box=_ValueControl(1.1); sig2Box=_ValueControl(2.2)
        voigt1Box=_ValueControl(.1); voigt2Box=_ValueControl(.2)
        lorentz1Box=_ValueControl(3.3); lorentz2Box=_ValueControl(4.4)
    values = PeakFitService(Legacy()).shape_parameters(2)
    assert values == {'sigmas': (1.1, 2.2), 'voigt': (0.1, 0.2), 'lorentz': (3.3, 4.4)}


def test_projection_threshold_initialises_missing_dmax_from_data():
    import numpy as np
    from spinDecon.analysis.projection_service import ProjectionService
    class Legacy:
        dmax=None
        data=np.asarray([-2.0, 5.0, 3.0])
        threshBox=_ValueControl(.2)
    legacy=Legacy()
    assert ProjectionService(legacy).intensity_threshold() == 1.0
    assert legacy.dmax == 5.0


def test_slice_service_full_peak_commit_runs_canonical_hooks():
    from spinDecon.analysis.slice_service import SliceService
    calls=[]
    class Store:
        def save_peak_list(self, name, **kwargs): calls.append(('save', name, kwargs['dimension']))
    class Legacy:
        dim=3; labb=('A','B','C'); store=Store()
        def _rebuild_projected_peak_lists(self): calls.append(('rebuild',))
        def _notify_analysis_changed(self): calls.append(('notify',))
        def refresh_full_peak_list_viewers(self): calls.append(('refresh',))
    SliceService(Legacy()).save_full_peak_records([], [], source_path='x')
    assert calls == [('save','full',3), ('rebuild',), ('notify',), ('refresh',)]


def test_slice_service_exposes_axes_peaks_and_sampling_without_gui_imports():
    import numpy as np
    from spinDecon.analysis.slice_service import SliceService
    class Control:
        def GetValue(self): return '2.0'
    class Peak:
        def __init__(self, name): self.name = name
    class Legacy:
        labb = ['H', 'N', 'C']
        peak = [Peak('A'), Peak('B')]
        pkIdx = [[1, 2, 3], [0, 1, 2]]
        index0 = np.array([1., 2.])
        data = np.arange(8.).reshape(2,2,2)
        datadec = data + 10
        sig2Box = Control()
    svc = SliceService(Legacy())
    assert svc.peak_names() == ['A', 'B']
    assert svc.axis(0).tolist() == [1., 2.]
    assert svc.sample((1,1,1)) == 7.0
    assert svc.sample((1,1,1), decon=True) == 17.0
    assert svc.peak_shape_width(2) == 18.0
    assert not hasattr(svc, 'connections')
