import numpy as np

from spinDecon.analysis.cpmg_service import _safe_residual, baldwin_r2eff, fit_local


def test_safe_residual_replaces_nonfinite_values():
    r = _safe_residual(np.array([1.0, np.nan, np.inf]), np.zeros(3))
    assert np.all(np.isfinite(r))
    assert r[1] == 1e6
    assert r[2] == 1e6


def test_local_fit_returns_finite_parameters_for_synthetic_curve():
    x = np.array([50., 100., 150., 250., 400., 600., 800.])
    y = baldwin_r2eff(x, 0.04, 0.03, 900., 12., 0.8 * 60.8)
    result = fit_local(x, y, 0.04, 60.8)
    assert result['valid'] is True
    assert result['success'] is True
    assert np.all(np.isfinite([result['pb'], result['kex'], result['R0'], result['dw'], result['Rex']]))


def test_fit_reports_r2_infinity():
    import numpy as np
    from spinDecon.analysis.cpmg_service import baldwin_r2eff, fit_local
    x=np.array([50.,100.,200.,400.,800.,1200.])
    y=baldwin_r2eff(x,0.04,0.03,800.,12.,90.)
    result=fit_local(x,y,0.04,75.)
    assert result['valid']
    assert np.isfinite(result['R2inf'])
    assert abs(result['R2inf'] - 12.0) < 0.1


def test_global_fit_exposes_report_metrics_for_each_peak():
    from spinDecon.analysis.cpmg_service import baldwin_r2eff, fit_global
    x = np.asarray([50., 100., 200., 400., 800., 1000.])
    curves = {}
    for name, r0, dw in [('P1', 12.0, 1.2), ('P2', 16.0, 0.8)]:
        y = baldwin_r2eff(x, 0.04, 0.04, 900.0, r0, dw * 60.0)
        curves[name] = {'x': x, 'y': y, 'e': np.ones_like(x) * 0.3}
    result = fit_global(curves, 0.04, 60.0)
    assert result['success']
    for peak in result['peaks'].values():
        assert peak['valid']
        for key in ('Rex', 'R2inf', 'R0line', 'chi2Line', 'chi2local', 'improvement'):
            assert key in peak
            assert np.isfinite(peak[key])
