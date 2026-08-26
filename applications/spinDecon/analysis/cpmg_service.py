"""Numerical helpers for pseudo-3D CPMG relaxation-dispersion analysis.

This module deliberately has no wx dependency.  The GUI owns presentation and
project I/O; curve construction and the Baldwin two-state model live here so
they can be tested and reused independently.
"""
from __future__ import annotations

import math
import numpy as np

GAMMA = {'15N': 27.116, '1H': 267.522, '19F': 251.815, '13C': 67.2828,
         '2H': 41.065, '31P': 108.291}


def observe_frequency_mhz(proton_frequency_mhz, nucleus):
    nucleus = str(nucleus).strip()
    if nucleus not in GAMMA:
        raise ValueError("Unsupported nucleus %r (known: %s)" % (nucleus, ', '.join(sorted(GAMMA))))
    return float(proton_frequency_mhz) / GAMMA['1H'] * GAMMA[nucleus]


def build_r2eff(axis, intensities, time_t2, errors=None, reference_tolerance=1e-8):
    axis = np.asarray(axis, dtype=float)
    intensities = np.asarray(intensities, dtype=float)
    if axis.size != intensities.size:
        raise ValueError('Pseudo-axis and intensity arrays have different lengths')
    time_t2 = float(time_t2)
    if not np.isfinite(time_t2) or time_t2 <= 0:
        raise ValueError('CPMG relaxation time must be greater than zero')
    refs = np.flatnonzero(np.abs(axis) <= reference_tolerance)
    if refs.size == 0:
        raise ValueError('The selected pseudo-axis has no zero-frequency reference plane')
    ref = float(np.mean(intensities[refs]))
    if not np.isfinite(ref) or abs(ref) <= np.finfo(float).eps:
        raise ValueError('Reference intensity is zero or invalid')
    mask = np.abs(axis) > reference_tolerance
    x, vals = axis[mask], intensities[mask]
    if np.any(np.abs(vals) <= np.finfo(float).eps):
        raise ValueError('A non-reference CPMG intensity is zero')
    y = np.log(np.abs(ref / vals)) / time_t2
    if errors is None:
        # Estimate repeat noise from duplicate frequencies; retain the historical
        # 0.3 fallback when no repeats are available.
        diffs = []
        for freq in np.unique(x):
            group = vals[np.isclose(x, freq)]
            if group.size > 1:
                diffs.extend((group[1:] - group[0]) ** 2)
        sigma = math.sqrt(float(np.sum(diffs))) / len(diffs) if diffs else 0.3
        e = np.sqrt((sigma / ref) ** 2 + (sigma / vals) ** 2) / time_t2
    else:
        err = np.asarray(errors, dtype=float)
        if err.size != axis.size:
            raise ValueError('Intensity-error and pseudo-axis arrays have different lengths')
        ref_err = float(np.sqrt(np.mean(err[refs] ** 2)))
        e = np.sqrt((ref_err / ref) ** 2 + (err[mask] / vals) ** 2) / time_t2
    order = np.argsort(x)
    return x[order], y[order], np.abs(e[order])


def baldwin_r2eff(nu_cpmg, time_t2, pb, kex, r0, dw):
    x = np.asarray(nu_cpmg, dtype=float)
    if np.any(x <= 0):
        raise ValueError('Baldwin CPMG model requires strictly positive nu_CPMG values')
    pb, kex, r0, dw, trel = map(float, (pb, kex, r0, dw, time_t2))
    delta_r2 = 0.0
    keg, kge = kex * (1-pb), kex * pb
    ncyc = trel * x
    tcp = trel / (4.0*ncyc)
    g1 = 2*dw*(delta_r2+keg-kge)
    g2 = (delta_r2+keg-kge)**2 + 4*keg*kge-dw**2
    root = (g1**2+g2**2)**0.25
    g3 = np.cos(0.5*np.arctan2(g1,g2))*root
    g4 = np.sin(0.5*np.arctan2(g1,g2))*root
    N = complex(g3, g4)
    NNc = g3**2+g4**2
    f0=(dw**2+g3**2)/NNc; f2=(dw**2-g4**2)/NNc
    t2=(dw+g4)*complex(dw,-g3)/NNc
    t1pt2=complex(2*dw**2,-g1)/NNc
    oGt2=complex((delta_r2+keg-kge-g3),(dw-g4))*t2
    rpre=(2*r0+kex)/2.0
    E0=2*tcp*g3; E2=2*tcp*g4; E1=complex(g3,-g4)*tcp
    ex0b=f0*np.cosh(E0)-f2*np.cos(E2)
    ex0c=f0*np.sinh(E0)-f2*np.sin(E2)*1j
    ex1c=np.sinh(E1)
    v3=np.sqrt(ex0b**2-1)
    y=np.power((ex0b-v3)/(ex0b+v3),ncyc)
    v2pPdN=(complex(delta_r2+kex,dw)*ex0c+(-oGt2-kge*t1pt2)*2*ex1c)
    Tog=((1+y)/2+(1-y)/(2*v3)*v2pPdN/N)
    return rpre-ncyc/trel*np.arccosh(ex0b.real)-1/trel*np.log(Tog.real)


def _safe_residual(model, y, err=None, invalid_penalty=1e6):
    """Return a finite residual vector, even for invalid model trials.

    The Baldwin expression can overflow for extreme trial parameters explored by
    least_squares.  Optimisers require finite residuals, so non-finite model
    points are converted to a large finite penalty rather than being allowed to
    poison the fit with NaN/Inf values.
    """
    yy = np.asarray(y, float)
    mm = np.asarray(model, float)
    r = mm - yy
    bad = ~np.isfinite(r)
    if err is not None:
        e = np.asarray(err, float)
        good = np.isfinite(e) & (e > 0)
        if np.any(good):
            fallback = float(np.nanmedian(e[good]))
            scale = np.where(good, e, fallback)
            r = r / scale
            bad |= ~np.isfinite(r)
    if np.any(bad):
        r = np.asarray(r, float).copy()
        r[bad] = float(invalid_penalty)
    return r



def r2_infinity(time_t2, pb, kex, r0, dw_rad_s, start_hz=1.0e4, rtol=1.0e-7, max_steps=12):
    """Return the converged high-nu_CPMG limit of the fitted R2eff model.

    The limit is evaluated numerically because it is a property of the fitted
    Baldwin model, not simply the fitted R0 parameter.
    """
    nu = float(start_hz)
    previous = None
    for _ in range(int(max_steps)):
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            value = np.asarray(baldwin_r2eff(np.asarray([nu], float), time_t2, pb, kex, r0, dw_rad_s), float).ravel()[0]
        if np.isfinite(value):
            value = float(value)
            if previous is not None and abs(value - previous) <= rtol * max(1.0, abs(value)):
                return value
            previous = value
        nu *= 2.0
    return float(previous) if previous is not None else float('nan')


def _parameter_errors(opt):
    """Approximate 1-sigma parameter errors from the least-squares Jacobian.

    Uses s^2 (J^T J)^-1 with s^2 = RSS/(m-n).  A pseudo-inverse keeps
    rank-deficient fits reportable; non-identifiable parameters return NaN.
    """
    try:
        jac=np.asarray(opt.jac,float)
        m,n=jac.shape
        if m <= n or not np.all(np.isfinite(jac)):
            return np.full(n,np.nan)
        rss=float(2.0*opt.cost)
        cov=np.linalg.pinv(jac.T.dot(jac))*rss/float(m-n)
        diag=np.diag(cov)
        return np.sqrt(np.where(diag >= 0.0,diag,np.nan))
    except Exception:
        try: return np.full(len(opt.x),np.nan)
        except Exception: return np.asarray([],float)

def fit_local(nu_cpmg, r2eff, time_t2, observe_mhz, errors=None, initial=None):
    """Fit one CPMG dispersion. Returns a GUI-independent result dictionary."""
    from scipy.optimize import least_squares
    x = np.asarray(nu_cpmg, float); y = np.asarray(r2eff, float)
    if x.size < 4:
        raise ValueError('At least four non-reference CPMG points are required for a local fit')
    flat = float(np.mean(y)); flat_res = _safe_residual(np.full(y.shape, flat), y, errors)
    chi2_line = float(np.mean(flat_res ** 2))
    guess = dict(pb=0.02, kex=1000.0, r0=float(np.min(y)), dw_ppm=1.0)
    if initial: guess.update({k: float(v) for k, v in initial.items() if k in guess})
    p0 = [guess['pb'], guess['kex'], guess['r0'], guess['dw_ppm']]
    # pb is a minor-state population; positive kex/R0/dw keep the model physical.
    lower = [1e-6, 1e-3, 0.0, 1e-8]
    upper = [0.499999, 1e7, max(1e4, float(np.max(y))*20.0), 1e4]
    # Keep an imported/remembered initial guess inside the admissible region.
    p0 = np.minimum(np.maximum(np.asarray(p0, float), np.asarray(lower) + 1e-12),
                    np.asarray(upper) - 1e-12)
    def residual(p):
        pb, kex, r0, dw_ppm = p
        try:
            with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
                model = baldwin_r2eff(x, time_t2, pb, kex, r0,
                                      dw_ppm * float(observe_mhz))
        except (ArithmeticError, FloatingPointError, ValueError, OverflowError):
            return np.full(y.shape, 1e6, dtype=float)
        return _safe_residual(model, y, errors)
    opt = least_squares(residual, p0, bounds=(lower, upper), max_nfev=10000,
                        x_scale='jac')
    pb, kex, r0, dw_ppm = opt.x
    perr = _parameter_errors(opt)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        model = baldwin_r2eff(x, time_t2, pb, kex, r0, dw_ppm * float(observe_mhz))
    finite_model = bool(np.all(np.isfinite(model)))
    finite_parameters = bool(np.all(np.isfinite([pb, kex, r0, dw_ppm])))
    valid = bool(opt.success and finite_model and finite_parameters)
    chi2 = float(np.mean(residual(opt.x) ** 2))
    rex = float(np.max(model) - np.min(model)) if finite_model else float('nan')
    r2inf = r2_infinity(time_t2, pb, kex, r0, dw_ppm * float(observe_mhz)) if valid else float('nan')
    improvement = float(1.0 - chi2/chi2_line) if valid and chi2_line > 0 else float('nan')
    return {'pb': float(pb), 'kex': float(kex), 'R0': float(r0), 'R2inf': r2inf, 'dw': float(dw_ppm),
            'R0line': flat, 'chi2Line': chi2_line, 'chi2local': chi2,
            'Rex': rex, 'improvement': improvement, 'model': np.asarray(model),
            'pb_error': float(perr[0]), 'kex_error': float(perr[1]), 'R0_error': float(perr[2]), 'dw_error': float(perr[3]),
            'success': valid, 'valid': valid,
            'message': str(opt.message) if valid else 'Non-finite or unsuccessful CPMG fit: %s' % opt.message}


def fit_global(curves, time_t2, observe_mhz, initial=None):
    """Global fit with shared kex/pb and local R0/dw for each supplied peak."""
    from scipy.optimize import least_squares
    names = list(curves)
    if not names:
        raise ValueError('No CPMG curves selected for global fitting')
    guess = dict(pb=0.02, kex=1000.0)
    if initial: guess.update({k: float(v) for k, v in initial.items() if k in guess})
    p0 = [guess['pb'], guess['kex']]
    lower = [1e-6, 1e-3]; upper = [0.499999, 1e7]
    for name in names:
        y = np.asarray(curves[name]['y'], float)
        finite_y = y[np.isfinite(y)]
        if finite_y.size == 0:
            raise ValueError('CPMG curve %s contains no finite R2,eff values' % name)
        # R0 is physically non-negative, but noisy/processed data can have a
        # slightly negative minimum.  Never let that make p0 violate bounds.
        r0_guess = max(0.0, float(np.min(finite_y)))
        p0 += [r0_guess, 1.0]
        lower += [0.0, 1e-8]
        upper += [max(1e4, float(np.max(finite_y))*20.0), 1e4]

    # scipy.least_squares rejects an initial point on the wrong side of any
    # bound before optimisation begins.  Imported/user starting values and
    # data-derived local guesses are therefore normalised into the admissible
    # interval.  Keep a tiny margin for strict-bound numerical edge cases.
    p0=np.asarray(p0,float); lower=np.asarray(lower,float); upper=np.asarray(upper,float)
    if not np.all(np.isfinite(p0)):
        raise ValueError('Global CPMG initial guesses must be finite numeric values')
    eps=np.maximum(1e-12,1e-10*np.maximum(1.0,np.abs(upper-lower)))
    p0=np.minimum(np.maximum(p0,lower+eps),upper-eps)

    def residual(p):
        pb, kex = p[:2]; out=[]; j=2
        for name in names:
            c=curves[name]; r0, dw_ppm = p[j:j+2]; j += 2
            model=baldwin_r2eff(c['x'], time_t2, pb, kex, r0, dw_ppm*float(observe_mhz))
            out.extend(_safe_residual(model, c['y'], c.get('e')))
        return np.asarray(out)
    opt=least_squares(residual, p0, bounds=(lower, upper), max_nfev=30000)
    pb, kex=opt.x[:2]; perr=_parameter_errors(opt); j=2; peak_results={}
    for name in names:
        c=curves[name]; r0,dw_ppm=opt.x[j:j+2]; j+=2
        model=baldwin_r2eff(c['x'],time_t2,pb,kex,r0,dw_ppm*float(observe_mhz))
        rr=_safe_residual(model,c['y'],c.get('e'))
        chi2 = float(np.mean(rr**2))
        flat = float(np.mean(np.asarray(c['y'], float)))
        flat_res = _safe_residual(np.full(np.asarray(c['y']).shape, flat), c['y'], c.get('e'))
        chi2_line = float(np.mean(flat_res ** 2))
        finite_model = bool(np.all(np.isfinite(model)))
        valid = bool(opt.success and finite_model and np.all(np.isfinite([pb, kex, r0, dw_ppm])))
        rex = float(np.max(model)-np.min(model)) if finite_model else float('nan')
        r2inf = r2_infinity(time_t2, pb, kex, r0, dw_ppm*float(observe_mhz)) if valid else float('nan')
        improvement = float(1.0-chi2/chi2_line) if valid and chi2_line > 0 else float('nan')
        peak_results[name]={'pb':float(pb),'kex':float(kex),'R0':float(r0),'R2inf':r2inf,
                            'dw':float(dw_ppm),'R0line':flat,'chi2Line':chi2_line,
                            'chi2local':chi2,'Rex':rex,'improvement':improvement,
                            'model':np.asarray(model),'pb_error':float(perr[0]),'kex_error':float(perr[1]),
                            'R0_error':float(perr[j-2]),'dw_error':float(perr[j-1]),'success':valid,'valid':valid}
    return {'pb':float(pb),'kex':float(kex),'pb_error':float(perr[0]),'kex_error':float(perr[1]),'peaks':peak_results,
            'success':bool(opt.success),'message':str(opt.message),
            'chi2':float(np.mean(residual(opt.x)**2))}
