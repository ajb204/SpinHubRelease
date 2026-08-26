#!/usr/bin/python

"""
Academic Use Licence

These licence terms apply to all licences granted by THE CHANCELLOR, MASTERS AND SCHOLARS OF THE UNIVERSITY OF OXFORD whose administrative offices are at University Offices, Wellington Square, Oxford OX1 2JD, United Kingdom (the "University") for use of UniDecNMR ("the Software") downloaded from the following website: https://github.com/charliebuchanan/UniDecNMR ("the Website")
By downloading the Software through the Source, you (the "Licensee") are confirming that you agree that your use of the Software is subject to these licence terms.

PLEASE READ THESE LICENCE TERMS CAREFULLY BEFORE DOWNLOADING THE SOFTWARE THROUGH THIS WEBSITE.  IF YOU DO NOT AGREE TO THESE LICENCE TERMS YOU SHOULD NOT DOWNLOAD THE SOFTWARE.

THE SOFTWARE IS INTENDED FOR USE BY ACADEMICS CARRYING OUT RESEARCH AND NOT FOR USE BY CONSUMERS OR COMMERCIAL BUSINESSES.

1.	Academic Use Licence
1.1	The Licensee is granted a limited non-exclusive and non-transferable royalty free licence to download and use the Software provided that the Licensee will:
(a)	limit their use of the Software to their own internal academic non-commercial research which is undertaken for the purposes of education or other scholarly use; 
(b)	not use the Software for or on behalf of any third party or to provide a service or integrate all or part of the Software into a product for sale or license to third parties;
(c)	use the Software in accordance with the prevailing instructions and guidance for use given on the Website and comply with procedures on the Website for user identification, authentication and access;
(d)	comply with all applicable laws and regulations with respect to their use of the Software; and 
(e)	ensure that the Copyright Notice "Copyright (c) 2022, University of Oxford" appears prominently wherever the Software is reproduced and on any documents or other material created using the Software.
1.2	The Licensee may only reproduce, modify, transmit or transfer the Software where:
(a)	such reproduction, modification, transmission or transfer is for academic, research or other scholarly use;
(b)	the conditions of this Licence are imposed upon the receiver of the Software or any modified Software;
(c)	all original and modified Source Code is included in any transmitted software program; and
(d)	the Licensee grants the University an irrevocable, indefinite, royalty free, non-exclusive unlimited licence to use and sub-licence any modified Source Code as part of the Software.

1.3	The University reserves the right at any time and without liability or prior notice to the Licensee to revise, modify and replace the functionality and performance of the access to and operation of the Software.
1.4	The Licensee acknowledges and agrees that the University owns all intellectual property rights in the Software.  The Licensee shall not have any right, title or interest in the Software.
1.5	This Licence will terminate immediately and the Licensee will no longer have any right to use the Software or exercise any of the rights granted to the Licensee upon any breach of the conditions in Section 1 of this Licence.

2.	Indemnity and Liability 
2.1	The Licensee shall defend, indemnify and hold harmless the University against any claims, actions, proceedings, losses, damages, expenses and costs (including without limitation court costs and reasonable legal fees) arising out of or in connection with the Licensee's possession or use of the Software, or any breach of these terms by the Licensee. 
2.2	The Software is provided on an 'as is' basis and the Licensee uses the Software at their own risk. No representations, conditions, warranties or other terms of any kind are given in respect of the the Software and all statutory warranties and conditions are excluded to the fullest extent permitted by law. Without affecting the generality of the previous sentences, the University gives no implied or express warranty and makes no representation that the Software or any part of the Software: (a) will enable specific results to be obtained; or (b) meets a particular specification or is comprehensive within its field or that it is error free or will operate without interruption; or (c) is suitable for any particular, or the Licensee's specific purposes. 
2.3	Except in relation to fraud, death or personal injury, the University's liability to the Licensee for any use of the Software, in negligence or arising in any other way out of the subject matter of these licence terms, will not extend to any incidental or consequential damages or losses, or any loss of profits, loss of revenue, loss of data, loss of contracts or opportunity, whether direct or indirect.
2.4	The Licensee hereby irrevocably undertakes to the University not to make any claim against any employee, student, researcher or other individual engaged by the University, being a claim which seeks to enforce against any of them any liability whatsoever in connection with these licence terms or their subject-matter. 

3.	General 
3.1	Severability - If any provision (or part of a provision) of these licence terms is found by any court or administrative body of competent jurisdiction to be invalid, unenforceable or illegal, the other provisions shall remain in force.
3.2	Entire Agreement - These licence terms constitute the whole agreement between the parties and supersede any previous arrangement, understanding or agreement between them relating to the Software. 
3.3	Law and Jurisdiction - These licence terms and any disputes or claims arising out of or in connection with them shall be governed by, and construed in accordance with, the law of England. The Licensee irrevocably submits to the exclusive jurisdiction of the English courts for any dispute or claim that arises out of or in connection with these licence terms.

If you are interested in using the Software commercially, please contact Oxford University Innovation Limited to negotiate a licence. Contact details are enquiries@innovation.ox.ac.uk 

"""

import numpy,sys,os,string,math,wx
import subprocess
import tempfile
import nmrglue as ng
from spinDecon.gui.dialogs.errors import errorMessage
from spinDecon.project.parameter_store import update_parameter_file
from spinDecon.gui.dialogs.shell_output import ShellOutputFrame, run_command_with_output
from spinDecon.processing.dimension_contract import legacy_vpar_dimension
from spinDecon.processing.nmrpipe_scripts import (
    nmrPipe as _nmrPipe,
    MakeProj2D,
    MakeProj3D,
    MakeProj3P,
    MakeProj4D,
)
from spinDecon.processing.script_context import ProcessingScriptState
from spinDecon.domain.dimensions.labels import clean_dimension_label, canonical_spectral_labels, discover_bruker_labels



import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar
import matplotlib.cm as cm
import nmrglue as ng
from matplotlib.figure import Figure
import posixpath
from wx.lib.mixins.listctrl import ColumnSorterMixin
from matplotlib.widgets import Cursor
import re

# Choose the direct-spectrum phasing backend in code, not in the GUI.
# Set to 'pipe' to use the nmrPipe phase script, or 'glue' to use nmrglue.
DIRECT_PHASE_BACKEND = 'glue'

def path_escape(indir):
    fields = indir.split(' ')
    pass
    final_string=''
    for x in range(len(fields)-1):
        if fields[x][-1] != '\\':
           final_string+=fields[x]+'\\ '
        else:
           final_string+=fields[x]

    final_string+=fields[-1]
    return final_string



def _autophase_score_trace(data, p0, p1=0.0):
    """Score a complex 1D spectrum for absorptive phasing.

    The function favors spectra whose displayed real part is mostly positive
    and whose imaginary part is small. It is intentionally simple and robust
    so that the GUI can use it repeatedly during a grid search.
    """
    arr = numpy.asarray(data)
    if arr.size == 0:
        return float('-inf')

    try:
        phased = ng.proc_base.ps(arr, p0=float(p0), p1=float(p1))
    except Exception:
        return float('-inf')

    real = -numpy.real(numpy.asarray(phased, dtype=numpy.complex128))
    imag = numpy.imag(numpy.asarray(phased, dtype=numpy.complex128))

    if not numpy.isfinite(real).all() or not numpy.isfinite(imag).all():
        return float('-inf')

    # Remove any small DC offset so the score is not dominated by the baseline.
    real = real - numpy.median(real)
    imag = imag - numpy.median(imag)

    pos = numpy.sum(numpy.clip(real, 0.0, None))
    neg = numpy.sum(numpy.clip(-real, 0.0, None))
    imag_penalty = numpy.sum(numpy.abs(imag))

    amplitude = numpy.sum(numpy.abs(real)) + 0.5 * numpy.sum(numpy.abs(imag)) + 1e-12
    # Reward a large positive absorptive signal while penalizing negative lobes
    # and residual dispersion. The constants are intentionally conservative so
    # that the routine prefers well-formed peaks without overfitting to noise.
    score = (pos - 1.75 * neg - 0.35 * imag_penalty) / amplitude
    return float(score)


def _autophase_grid_candidates(center, span, step):
    if step <= 0:
        return [float(center)]
    lo = float(center) - float(span)
    hi = float(center) + float(span)
    values = numpy.arange(lo, hi + (step * 0.5), float(step), dtype=float)
    if values.size == 0:
        return [float(center)]
    return [float(v) for v in values]


def _autophase_grid_search(data, p1=0.0, initial=0.0):
    """Return the best P0 candidate for a phased complex trace."""
    best_p0 = float(initial)
    best_score = float('-inf')

    # Three-stage search: broad, local, then fine. This keeps the logic robust
    # while avoiding expensive per-point fitting.
    stages = [
        (0.0, 180.0, 15.0),
        (None, 15.0, 3.0),
        (None, 3.0, 0.3),
    ]

    center = float(initial)
    for stage_idx, (stage_center, span, step) in enumerate(stages):
        if stage_center is not None:
            center = float(stage_center)
        candidates = _autophase_grid_candidates(center, span, step)
        stage_best_p0 = center
        stage_best_score = float('-inf')
        for cand in candidates:
            score = _autophase_score_trace(data, cand, p1=p1)
            if (score > stage_best_score + 1e-12) or (
                abs(score - stage_best_score) <= 1e-12 and abs(cand - center) < abs(stage_best_p0 - center)
            ):
                stage_best_score = score
                stage_best_p0 = cand
        center = stage_best_p0
        if stage_best_score > best_score or stage_idx == 0:
            best_score = stage_best_score
            best_p0 = stage_best_p0
    return float(best_p0), float(best_score)



def _pseudo_voigt_model(x, amp, center, fwhm, eta, offset, slope):
    x = numpy.asarray(x, dtype=float)
    fwhm = max(float(fwhm), 1e-9)
    eta = min(max(float(eta), 0.0), 1.0)
    dx = (x - float(center)) / fwhm
    gauss = numpy.exp(-4.0 * numpy.log(2.0) * dx * dx)
    lorentz = 1.0 / (1.0 + 4.0 * dx * dx)
    return float(offset) + float(slope) * (x - float(center)) + float(amp) * ((1.0 - eta) * gauss + eta * lorentz)


def _fit_pseudo_voigt_window(x, y):
    """Fit a simple pseudo-Voigt model to a 1D peak window.

    Returns None when the fit is not stable enough to trust.
    """
    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    if x.size < 7 or y.size < 7 or x.size != y.size:
        return None

    finite = numpy.isfinite(x) & numpy.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size < 7:
        return None

    edge_n = max(2, int(round(x.size * 0.12)))
    edge_n = min(edge_n, max(2, x.size // 2))
    edge = numpy.concatenate([y[:edge_n], y[-edge_n:]]) if edge_n > 0 else y
    baseline = float(numpy.median(edge)) if edge.size else 0.0
    y0 = y - baseline

    peak_idx = int(numpy.nanargmax(y0))
    peak_val = float(y0[peak_idx])
    if not numpy.isfinite(peak_val) or peak_val <= 0:
        peak_idx = int(numpy.nanargmax(numpy.abs(y0)))
        peak_val = float(abs(y0[peak_idx]))
    if not numpy.isfinite(peak_val) or peak_val <= 0:
        return None

    center = float(x[peak_idx])
    span = float(max(numpy.max(x) - numpy.min(x), 1e-6))
    if x.size > 1:
        dx = float(numpy.median(numpy.abs(numpy.diff(x))))
    else:
        dx = span / 10.0
    fwhm0 = max(dx * 8.0, span / 14.0, 1e-4)
    p0 = [peak_val, center, fwhm0, 0.5, baseline, 0.0]
    lower = [0.0, float(numpy.min(x)), max(dx * 0.5, 1e-9), 0.0, -numpy.inf, -numpy.inf]
    upper = [numpy.inf, float(numpy.max(x)), span, 1.0, numpy.inf, numpy.inf]

    try:
        from scipy.optimize import curve_fit
        popt, _pcov = curve_fit(_pseudo_voigt_model, x, y, p0=p0, bounds=(lower, upper), maxfev=20000)
        fit = _pseudo_voigt_model(x, *popt)
        resid = float(numpy.sqrt(numpy.mean((y - fit) ** 2)))
        amp = float(max(abs(popt[0]), 1e-12))
        span_norm = float(max(numpy.max(y) - numpy.min(y), amp, 1e-12))
        return {
            'fit': numpy.asarray(fit, dtype=float),
            'popt': popt,
            'residual': resid / span_norm,
            'fwhm': float(abs(popt[2])),
            'eta': float(min(max(popt[3], 0.0), 1.0)),
            'amp': amp,
            'baseline': float(popt[4]),
            'slope': float(popt[5]),
        }
    except Exception:
        # Moment-based fallback: enough to keep the grid-search robust when the
        # pseudo-Voigt fit is unstable because the peak is too flat or noisy.
        y_pos = numpy.clip(y0, 0.0, None)
        total = float(numpy.sum(y_pos))
        if total <= 0:
            return None
        mu = float(numpy.sum(x * y_pos) / total)
        var = float(numpy.sum(((x - mu) ** 2) * y_pos) / total)
        sigma = max(math.sqrt(max(var, 1e-12)), dx)
        fwhm = 2.355 * sigma
        fit = _pseudo_voigt_model(x, peak_val, mu, fwhm, 0.5, baseline, 0.0)
        resid = float(numpy.sqrt(numpy.mean((y - fit) ** 2)))
        span_norm = float(max(numpy.max(y) - numpy.min(y), peak_val, 1e-12))
        return {
            'fit': numpy.asarray(fit, dtype=float),
            'popt': numpy.array([peak_val, mu, fwhm, 0.5, baseline, 0.0], dtype=float),
            'residual': resid / span_norm,
            'fwhm': float(fwhm),
            'eta': 0.5,
            'amp': float(peak_val),
            'baseline': float(baseline),
            'slope': 0.0,
        }


def _apodization_candidate_score(x, y_real, y_imag, natural_fwhm=None):
    """Score an apodized 1D spectrum.

    Lower is better. The score prefers a compact pseudo-Voigt-like peak with
    a small imaginary component and minimal truncation ripple.
    """
    x = numpy.asarray(x, dtype=float)
    y_real = numpy.asarray(y_real, dtype=float)
    y_imag = numpy.asarray(y_imag, dtype=float)
    if x.size < 7 or y_real.size != x.size or y_imag.size != x.size:
        return {
            'score': float('inf'),
            'peak_index': 0,
            'sign': 1.0,
            'fwhm': float('nan'),
            'residual': float('inf'),
            'imag_penalty': float('inf'),
            'neg_penalty': float('inf'),
            'tail_penalty': float('inf'),
            'fit': None,
            'popt': None,
        }

    peak_index = int(numpy.nanargmax(numpy.abs(y_real)))
    sign = 1.0 if float(y_real[peak_index]) >= 0.0 else -1.0
    y = sign * y_real
    imag = sign * y_imag

    if not numpy.isfinite(y).all() or not numpy.isfinite(imag).all():
        return {
            'score': float('inf'),
            'peak_index': peak_index,
            'sign': sign,
            'fwhm': float('nan'),
            'residual': float('inf'),
            'imag_penalty': float('inf'),
            'neg_penalty': float('inf'),
            'tail_penalty': float('inf'),
            'fit': None,
            'popt': None,
        }

    local_peak = int(numpy.nanargmax(y))
    peak_height = float(max(y[local_peak], numpy.max(numpy.abs(y)), 1e-12))
    threshold = peak_height * 0.20
    left = local_peak
    while left > 1 and y[left] > threshold:
        left -= 1
    right = local_peak
    last = y.size - 1
    while right < last - 1 and y[right] > threshold:
        right += 1
    pad = max(4, int(round(y.size * 0.06)))
    lo = max(0, left - pad)
    hi = min(y.size, right + pad + 1)
    if hi - lo < 9:
        half = max(6, y.size // 6)
        lo = max(0, local_peak - half)
        hi = min(y.size, local_peak + half + 1)

    xw = x[lo:hi]
    yw = y[lo:hi]
    iw = imag[lo:hi]
    fit = _fit_pseudo_voigt_window(xw, yw)

    if fit is None:
        # Fall back to a simple robustness score when the peak is too messy to
        # fit cleanly; this still allows the grid search to continue.
        peak = float(max(numpy.max(yw), 1e-12))
        residual = float(numpy.sqrt(numpy.mean((yw - numpy.median(yw)) ** 2)) / peak)
        imag_penalty = float(numpy.mean(numpy.abs(iw)) / peak)
        neg_penalty = float(numpy.sum(numpy.clip(-yw, 0.0, None)) / (numpy.sum(numpy.abs(yw)) + 1e-12))
        tail_penalty = float(numpy.mean(numpy.abs(yw[:max(2, yw.size // 6)]) ) / peak)
        fwhm = float('nan')
        score = residual + 0.45 * imag_penalty + 0.75 * neg_penalty + 0.20 * tail_penalty
        return {
            'score': float(score),
            'peak_index': peak_index,
            'sign': sign,
            'fwhm': fwhm,
            'residual': residual,
            'imag_penalty': imag_penalty,
            'neg_penalty': neg_penalty,
            'tail_penalty': tail_penalty,
            'fit': None,
            'popt': None,
            'window': (lo, hi),
        }

    fit_y = numpy.asarray(fit['fit'], dtype=float)
    peak = float(max(fit['amp'], numpy.max(numpy.abs(yw)), 1e-12))
    residual = float(fit['residual'])
    imag_penalty = float(numpy.mean(numpy.abs(iw)) / peak)
    neg_penalty = float(numpy.sum(numpy.clip(-yw, 0.0, None)) / (numpy.sum(numpy.abs(yw)) + 1e-12))
    edge_n = max(2, int(round(yw.size * 0.18)))
    tail_idx = numpy.r_[0:edge_n, max(edge_n, yw.size - edge_n):yw.size]
    tail_penalty = float(numpy.mean(numpy.abs(yw[tail_idx] - fit_y[tail_idx])) / peak)

    fwhm = float(abs(fit['fwhm'])) if numpy.isfinite(fit['fwhm']) else float('nan')
    broad_penalty = 0.0
    if natural_fwhm is not None and numpy.isfinite(natural_fwhm) and natural_fwhm > 0 and numpy.isfinite(fwhm):
        broad_penalty = max(0.0, (fwhm / natural_fwhm) - 1.0)

    score = residual + 0.40 * imag_penalty + 0.65 * neg_penalty + 0.25 * tail_penalty + 0.90 * broad_penalty
    return {
        'score': float(score),
        'peak_index': peak_index,
        'sign': sign,
        'fwhm': fwhm,
        'residual': residual,
        'imag_penalty': imag_penalty,
        'neg_penalty': neg_penalty,
        'tail_penalty': tail_penalty,
        'fit': fit_y,
        'popt': fit['popt'],
        'window': (lo, hi),
    }

class ProcMan(wx.App):
    def __init__(self,inherit,showFlg=True):
        self.frame_ProcessFrame=ProcessFrame(None,20,'Process',inherit,showFlg=showFlg)
        if(showFlg):
            self.frame_ProcessFrame.Show(True)
#        return Frame1(parent)

        

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]

class ProcessFrame(wx.Frame):

    def __init__(self,parent,id,title,inherit,showFlg=True):
        #wx.Panel.__init__(self, parent=parent)
        #self.parent=parent
        #self.tabOne=parent.tabOne
        #self.create_main_panel()
        #self.draw_figure()
        #self.canvas.draw()

        self.parent=inherit
        self.state = getattr(inherit, "state", None)
        self.showFlg=showFlg
        self.conv_frame = None
        self.processing_frame = None
        self.projections_frame = None
        #try:
        self.ncpus=float(self.parent.coreBox.GetValue())
        self.p0_orig = 0
        self.p1_orig = 0
        if(self.ncpus==0):
            self.ncpus=1
        #except:
        #    self.ncpus=1
        self.ind_phase = False
        self._direct_phase_p0 = 0.0
        self._direct_phase_p1 = 0.0
        self._syncing_direct_phase_controls = False
        self.FIDsel = 1
        #self.spectral_dim_count=int(decon.Parse(self.parent.deconParFile,'dim'))
        # Canonical Stage-4 dimensionality contract.  ProcessFrame.dim now
        # always means *spectral* dimensions.  Pseudo/physical dimensionality
        # is represented explicitly and never encoded as '2p'/'3p'.
        # The main-window controls are the UI boundary at which dimensionality
        # may change.  Synchronise them into ProjectState once, then derive every
        # ProcessFrame dimensionality value from the canonical topology.  This
        # prevents stale state and GUI values from becoming competing sources of
        # truth (notably for 2D + pseudo-axis datasets).
        gui_spectral_dim_count = max(1, int(self.parent.dim))
        gui_has_pseudo_axis = bool(self.state.pseudo_axis) if self.state is not None else bool(getattr(self.parent, 'pseudo', False))

        if self.state is not None:
            self.state.sync_from_values(
                spectral_dimensions=gui_spectral_dim_count,
                pseudo_axis=gui_has_pseudo_axis,
            )
            self.topology = self.state.topology()
        else:
            # Compatibility for isolated/tests callers without ProjectState.
            from spinDecon.domain.topology import DatasetTopology
            self.topology = DatasetTopology.from_counts(
                gui_spectral_dim_count, gui_has_pseudo_axis
            )

        # Compatibility aliases are read-only snapshots of topology.  New code
        # should consume self.topology directly.
        self.spectral_dim_count = self.topology.spectral_dim_count
        self.has_pseudo_axis = self.topology.has_pseudo_axis
        self.physical_dim_count = self.topology.physical_dim_count
        self.dim = self.topology.spectral_dim_count  # spectral only
        # FT selection choices are required by the nmrPipe script builders even
        # when the indirect-dimensions box is not shown in this frame. Keep the
        # canonical list here so script generation does not depend on GUI state.
        self.ftlisty=['Auto', 'Neg', 'Alt', 'AltNeg', 'Real']
        self.ftdic={ft:i for i, ft in enumerate(self.ftlisty)}
        self.ftdic['y']=1
        self.ftdic['n']=0
        self.statusPanel = None
        # The raw-data path is fixed for the lifetime of this Process window.
        # Detect its spectrometer format once and let conversion/processing
        # consumers reuse self.tp rather than repeatedly probing the filesystem.
        self.tp = None

        # Do not mutate dim for pseudo datasets.  All children consume the
        # explicit spectral/physical counts above.
        
        pass

        self.READ1D=0   #have we read in the 1Dslice?
        self.nmrPipe = getattr(self.parent, "nmrPipe", _nmrPipe(self.parent))
        if getattr(self.parent, "nmrPipe", None) is None:
            self.parent.nmrPipe = self.nmrPipe

        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title='Process NMR spectrum ...')
        self.SetClientSize(wx.Size(900, 280))
        self.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.GetSpectrometerType()
        try:
            self.CreateStatusBar(1)
            self.SetStatusText('Ready')
        except Exception:
            pass

        panel=wx.Panel(self,-1)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    

        self.create_main_panel()

        if(showFlg):
            self.Show(True)
        self.Fit()
        # A Process session owns the canonical live dimension labels.  Hydrate
        # them from the saved system file before any child window can consume
        # defaults.  Conversion subsequently edits these shared values; it is
        # no longer responsible for discovering them.
        self._hydrate_shared_dimension_labels()
        self.GetLabs()
        self.SetLab()
        # Register the live Process family with the owning project and capture
        # the persistence baseline only after the Process controls are loaded.
        try:
            self.parent.process_frame = self
        except Exception:
            pass
        self._capture_parameter_baseline()

    def SetWin(self,win,win2val,win3val, firstPoint): #reset sin bell parameters if on defaults
        if(win.GetValue()=='SP'):
            if(win2val.GetValue()==str(0)):
                win2val.SetValue(str(0.5))
            if(win3val.GetValue()==str(0)):
                win3val.SetValue(str(2.))
            if(str(firstPoint.GetValue())==str(0)):
                firstPoint.SetValue(str(0.5))
        if(win.GetValue()=='GM'):
            if(win2val.GetValue()==str(0)):
                win2val.SetValue(str(2.))
            if(win3val.GetValue()==str(0)):
                win3val.SetValue(str(20))
            if(firstPoint.GetValue()==str(0)):
                firstPoint.SetValue(str(0.5))

        if(win.GetValue()=='EM'):
            if(win2val.GetValue()==str(0)):
                win2val.SetValue(str(1.))
            if(win3val.GetValue()==str(0)):
                win3val.SetValue(str(0.))
            if(str(firstPoint.GetValue())==str(0)):
                firstPoint.SetValue(str(0.5))

    #FGA added
    def onGetFile(self, e, textBox,store=False):

        if(store):
            cwd = str(self._raw_output_dir()).strip()
            if not os.path.exists(cwd):
                pass
                
        else:
            #get dialog box here
            cwd = os.getcwd()
            cwd = os.path.join(cwd, self._raw_output_dir().replace('./', ''))
            cwd = os.path.join(cwd, '')

        # cwd=self._raw_output_dir()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=cwd, defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            pass
            textBox.SetValue(splitPath[1])
            pass
            pass
        dlg.Destroy()
    """
def onGetFile(self, e, textBox,full=False,save=False):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.list)|*.list|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if(full==False):
                splitPath = path.split(cwd)
                    textBox.SetValue('.' + splitPath[1])
                print("You chose the following file(s):")
                print(path)
            else:
                print(path)
                textBox.SetValue(str(path))
                print("You chose the following file(s):")
                print(path)
                save(True) #execute save function (peaklist file)

        dlg.Destroy()
    """
    def onGetDir(self, e, textBox,full=False,default=''):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.DirDialog(self, message="Choose a folder",         style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,defaultPath=default)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            textBox.SetValue('.' + splitPath[1])
            #try:
            #    textBox.SetValue('.' + splitPath[1])
            #except:
            #    textBox.SetValue(path)
            pass
            pass
            #os.chdir(path)
            #self.dirBox.SetValue(path)
            #print "CWD: ",os.getcwd()

        dlg.Destroy()

    def OnSliderScroll(self,event):
        if getattr(self, '_autophase_in_progress', False):
            return
        if self.ind_phase == False  and type(self.ind_phase)==bool:

            self._set_direct_phase_values(self._phase_slider_value(self.sld_0, self.sld_0_mode_btn), self._phase_slider_value(self.sld_1, self.sld_1_mode_btn))
            show_fid = False
            try:
                show_fid = bool(self.cb_show_fid.IsChecked())
            except Exception:
                show_fid = False

            pass

            if show_fid:
                try:
                    self.draw_figure()
                    return
                except Exception:
                    pass

            try:
                p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
                p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
                if not hasattr(self, 'direct_ft_data'):
                    self._build_direct_frequency_data()
                x, y = self._render_direct_frequency(p0, p1)
                self.phasing.set_data(x, y)
                self._debug_plot_state('after slider spectrum update')
                # Phase changes are data-only updates.  In particular, do not
                # touch either axis limit here: the toolbar Redraw action is the
                # only operation allowed to reset the Y zoom.
                self._update_phase_readouts(p0, p1, show_fid=False)
                self._blit_phase_preview()
            except Exception:
                pass
        else:
            p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
            p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
            self.p0s[self.ind_phase].SetValue(str(p0))
            self.p1s[self.ind_phase].SetValue(str(p1))
            self.change_ind_phasing(p0, p1)

    def on_xrange_enter(self, event):
        """Apply Min/Max ppm edits immediately while preserving Y zoom."""
        self.draw_figure(reset_y=False)

    @staticmethod
    def _estimate_trace_noise(values):
        """Return robust (centre, sigma) using the deconMain noise model."""
        sample = numpy.asarray(values, dtype=float).ravel()
        sample = sample[numpy.isfinite(sample)]
        if sample.size < 20:
            raise ValueError('At least 20 finite points are required for automatic ranging.')

        mu = float(numpy.median(sample))
        mad = float(numpy.median(numpy.abs(sample - mu)))
        sigma = 1.482602218505602 * mad
        if not numpy.isfinite(sigma) or sigma <= 0:
            sigma = float(numpy.std(sample))
        if not numpy.isfinite(sigma) or sigma <= 0:
            raise ValueError('Spectrum has zero/undefined noise.')

        for _ in range(8):
            core = sample[numpy.abs(sample - mu) <= 2.5 * sigma]
            if core.size < 20:
                break
            new_mu = float(numpy.mean(core))
            new_sigma = float(numpy.std(core)) / numpy.sqrt(0.9112563609)
            if not numpy.isfinite(new_sigma) or new_sigma <= 0:
                break
            converged = abs(new_sigma - sigma) <= 1e-5 * sigma
            mu, sigma = new_mu, new_sigma
            if converged:
                break
        return mu, sigma

    def _current_full_phased_trace(self):
        """Return ppm axis and full direct trace with the current GUI phase."""
        dic, data = self._build_direct_frequency_data(include_phased=False)
        pipe_proc = self._pipe_proc_module()
        if pipe_proc is None:
            raise RuntimeError('nmrglue.pipe_proc is not available')
        p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
        p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
        dic, data = pipe_proc.ps(dic, data, p0=p0, p1=p1)
        uc = ng.pipe.make_uc(dic, data)
        return numpy.asarray(uc.ppm_scale(), dtype=float), -numpy.real(numpy.asarray(data, dtype=complex))

    def _automatic_xrange(self, amide_only=False):
        """Estimate a useful 1H display range from the currently phased trace."""
        ppm, trace = self._current_full_phased_trace()
        valid = numpy.isfinite(ppm) & numpy.isfinite(trace)
        if amide_only:
            # Keep AutoAmide clear of the water region: only assess signal at
            # 5.0 ppm and above.
            valid &= ppm >= 5.0
        if numpy.count_nonzero(valid) < 20:
            raise ValueError('Not enough spectral points in the requested region.')

        region_values = trace[valid]
        mu, sigma = self._estimate_trace_noise(region_values)
        # A phased spectrum should have a noise baseline.  Treat a point as
        # signal only at SNR >= 5, using absolute departure from the robust
        # noise centre so negative residual peaks are not lost.  The higher
        # threshold prevents broad/low-level baseline distortion from setting
        # the automatic display limits.
        signal = valid & (numpy.abs(trace - mu) >= 5.0 * sigma)

        # Reject isolated Gaussian-noise excursions.  A real resonance must
        # occupy at least two adjacent sampled points; one-point gaps are
        # bridged so a noisy peak top does not split a resonance.
        sig = signal.astype(bool)
        if sig.size >= 3:
            sig[1:-1] |= sig[:-2] & sig[2:]
        neighbours = numpy.zeros_like(sig, dtype=int)
        neighbours[1:] += sig[:-1]
        neighbours[:-1] += sig[1:]
        sig &= neighbours > 0
        idx = numpy.flatnonzero(sig)
        if idx.size == 0:
            raise ValueError('No signal could be distinguished from the noise.')

        signal_low = float(numpy.min(ppm[idx]))
        signal_high = float(numpy.max(ppm[idx]))

        if amide_only:
            # AutoAmide searches only at >= 5 ppm, but 5 ppm is a search
            # boundary rather than a forced display boundary.  Walking from
            # 5 ppm towards higher ppm, the first detected SNR>=5 signal sets
            # the low edge of the region of interest.  Padding is then based
            # on the detected signal width, so a spectrum whose first real
            # resonance starts at (for example) 6.5 ppm will zoom in towards
            # that resonance instead of remaining stuck at 5.0 ppm.
            roi_low = signal_low
            roi_high = signal_high
        else:
            # AutoCentre retains the expected 4.6--4.8 ppm water region.
            roi_low = min(signal_low, 4.6)
            roi_high = max(signal_high, 4.8)
        width = roi_high - roi_low
        if not numpy.isfinite(width) or width <= 0:
            raise ValueError('Could not determine a finite signal width.')
        pad = 0.10 * width
        acquired_low = float(numpy.nanmin(ppm))
        acquired_high = float(numpy.nanmax(ppm))
        xmin = max(acquired_low, roi_low - pad)
        xmax = min(acquired_high, roi_high + pad)
        if amide_only:
            # 5 ppm remains the hard lower bound of AutoAmide: padding may
            # extend towards it, but can never re-enter the water region.
            xmin = max(xmin, 5.0)
        return xmin, xmax

    def _apply_automatic_xrange(self, amide_only=False):
        try:
            xmin, xmax = self._automatic_xrange(amide_only=amide_only)
            self.xminBox.SetValue(f'{xmin:.3f}')
            self.xmaxBox.SetValue(f'{xmax:.3f}')
            # Apply the new limits immediately.  Preserve the user's Y zoom.
            self.draw_figure(reset_y=False)
        except Exception:
            errorMessage('Could not determine an automatic spectral range.')

    def OnAutoCentre(self, event):
        self._apply_automatic_xrange(amide_only=False)

    def OnAutoAmide(self, event):
        self._apply_automatic_xrange(amide_only=True)

    def on_show_fid(self, event):
        # Switching between the frequency-domain spectrum and the raw FID
        # changes the plotted data (and often its amplitude by orders of
        # magnitude).  Treat either checkbox transition as a new view and
        # recompute the Y limits from the data that are about to be shown.
        # Other redraw paths continue to preserve the user's Y zoom.
        self.draw_figure(reset_y=True)

    def on_fid_select(self, event):
        self._set_fid_selection(self._current_fid_selection())
        self.draw_figure()

    def on_phase_slice_mode(self, event):
        """Remember which pseudo-2D time-domain trace is used for phasing."""
        value = self.phaseSliceMode.GetStringSelection() or 'First'
        self.phaseSliceModeValue = value
        if getattr(self, 'state', None) is not None:
            self.state.metadata['phaseSliceMode'] = value

    def GetReferenceMode(self, event=None): ###NOT IN CHARLIE'S VERSION###
        ref_map = {0: 'Water', 1: 'Auto', 2: 'Manual'}
        try:
            ref_idx = int(self.parent.Parse(self.parent.deconParFile, 'refBox'))
        except Exception:
            ref_idx = 0
        self.reference_mode = ref_map.get(ref_idx, 'Water')

    def _make_status_lamp(self, parent, colour):
        # Match the compact bordered status lamps used on deconMain.  The
        # border gives the indicator a distinct lamp-like face rather than a
        # flat block of colour.
        lamp = wx.Panel(parent, size=(14, 14), style=wx.BORDER_SIMPLE)
        lamp.SetMinSize((14, 14))
        lamp.SetBackgroundColour(wx.Colour(colour))
        return lamp

    def _dimension_status_text(self):
        return (f'{self.topology.spectral_dim_count} spectral / '
                f'{self.topology.physical_dim_count} physical '
                f'(pseudo-axis: {"yes" if self.topology.has_pseudo_axis else "no"})')

    def _format_status_text(self):
        tp = getattr(self, 'tp', '')
        if tp == 'bruk':
            return 'bruker'
        if tp == 'var':
            return 'varian'
        return 'unknown'

    def _raw_output_dir(self):
        if getattr(self, 'state', None) is not None:
            return self.state.raw_dir()
        return './raw'

    def _spec_output_dir(self):
        spec = self.state.spec_dir() if getattr(self, 'state', None) is not None else './spec'
        try:
            os.makedirs(spec, exist_ok=True)
        except Exception:
            pass
        return spec


    def _output_path(self, *parts):
        base = self._spec_output_dir()
        return os.path.join(base, *parts) if parts else base

    def _processed_spectrum_name(self, pipefile=None):
        base = self._spec_output_dir()
        if not base:
            return pipefile or ''

        candidates = []
        for name in (pipefile, 'slice.phased.ft1', 'slice.ft1', 'test.ft1', 'test.ft2', 'test.ft3', 'test.ft4', 'test.ft'):
            if name:
                candidates.extend([name, name + '.gz'])

        seen = set()
        for name in candidates:
            if not name or name in seen:
                continue
            seen.add(name)
            path = os.path.join(base, name)
            if os.path.exists(path):
                return os.path.basename(name[:-3]) if name.endswith('.gz') else os.path.basename(name)

        import glob
        hits = [path for path in glob.glob(os.path.join(base, '*.ft*')) if os.path.isfile(path)]
        if hits:
            hits.sort(key=lambda path: os.path.getmtime(path), reverse=True)
            return os.path.basename(hits[0])
        return pipefile or ''

    def _update_nmrpipe_file_box(self, pipefile=None):
        pipefile = pipefile or self._process_pipefile()
        spectrum_name = self._processed_spectrum_name(pipefile)
        if not spectrum_name:
            spectrum_name = pipefile or ''
        parent = getattr(self, 'parent', None)
        if parent is not None:
            try:
                if hasattr(parent, 'infileBox'):
                    parent.infileBox.SetValue(os.path.basename(spectrum_name))
            except Exception:
                pass
            try:
                if hasattr(parent, 'state') and parent.state is not None:
                    parent.state.sync_from_values(
                        input_file=os.path.basename(spectrum_name),
                        spec_path=(parent.specPathBox.GetValue() if hasattr(parent, 'specPathBox') else self.state.spec_path),
                    )
            except Exception:
                pass
            # SetValue normally emits EVT_TEXT, but processing completion is a
            # state boundary: refresh the lamp and workflow explicitly so the
            # GUI advances even if a platform suppresses/programmatically
            # coalesces that event.
            try:
                refresh_lamps = getattr(parent, 'update_project_lamps', None)
                if callable(refresh_lamps):
                    refresh_lamps()
            except Exception:
                pass
            try:
                notebook = getattr(parent, 'parent', None)
                notify = getattr(notebook, 'notify_analysis_changed', None)
                if callable(notify):
                    notify()
            except Exception:
                pass
        return spectrum_name

    def _conversion_outputs_exist(self):
        # Conversion output layout depends on *physical* topology.  In
        # particular, 2 spectral + 1 pseudo is converted into a series under
        # spec/fids rather than the ordinary 2D spec/test.fid file.  Re-use
        # the same candidate discovery as the direct-FID reader so the status
        # lamp and the code which actually opens converted data cannot drift.
        try:
            candidates = self._direct_fid_candidates(show_fid=True)
        except Exception:
            candidates = []
        return any(os.path.isfile(path) for path in candidates)

    def _process_pipefile(self):
        return self.nmrPipe.pipefile_for(self)

    def _processed_output_exists(self):
        base = self._spec_output_dir()
        if not base:
            return False

        # The canonical processed filename is spectral-dimensional: pseudo-3D
        # is therefore test.ft2, not test.ft3.  Also accept the same fallback
        # products used by _processed_spectrum_name() so the lamp reflects a
        # successfully created spectrum even for legacy/project-specific names.
        pipefile = self._process_pipefile()
        candidates = [pipefile, 'slice.phased.ft1', 'slice.ft1']
        for name in candidates:
            if not name:
                continue
            path = os.path.join(base, name)
            if os.path.isfile(path) or os.path.isfile(path + '.gz'):
                return True

        try:
            import glob
            return any(os.path.isfile(path) for path in glob.glob(os.path.join(base, '*.ft*')))
        except Exception:
            return False

    def _set_lamp(self, ctrl, on):
        if ctrl is None:
            return
        try:
            ctrl.SetBackgroundColour(wx.Colour(46, 160, 67) if on else wx.Colour(210, 55, 55))
            ctrl.Refresh()
        except Exception:
            pass

    def UpdateLampLights(self):
        if not hasattr(self, 'dimensionStatusValue'):
            return
        try:
            self.dimensionStatusValue.SetLabel(self._dimension_status_text())
        except Exception:
            pass
        try:
            self.formatStatusValue.SetLabel(self._format_status_text())
        except Exception:
            pass
        self._set_lamp(self.convertedLamp, self._conversion_outputs_exist())
        self._set_lamp(self.processedLamp, self._processed_output_exists())
        try:
            if getattr(self, 'statusPanel', None) is not None:
                self.statusPanel.Layout()
            self.Layout()
            self.Refresh()
        except Exception:
            pass


    def _set_status_text(self, text):
        try:
            if getattr(self, 'status_bar', None) is not None:
                self.status_bar.SetStatusText(str(text or ''))
        except Exception:
            pass

    def _hover_lookup_text(self, widget):
        hover_map = getattr(self, '_hover_help_map', {})
        current = widget
        while current is not None:
            try:
                text = hover_map.get(current)
            except Exception:
                text = None
            if text:
                return text
            try:
                current = current.GetParent()
            except Exception:
                current = None
        return self._hover_default_status

    def _on_hover_event(self, event):
        try:
            widget = event.GetEventObject()
        except Exception:
            widget = None
        self._set_status_text(self._hover_lookup_text(widget))
        try:
            event.Skip()
        except Exception:
            pass

    def _on_hover_leave_frame(self, event):
        self._set_status_text(self._hover_default_status)
        try:
            event.Skip()
        except Exception:
            pass

    def _install_hover_help(self, widget, text):
        if widget is None:
            return
        try:
            if not hasattr(self, '_hover_help_map'):
                self._hover_help_map = {}
            self._hover_help_map[widget] = text
        except Exception:
            pass
        try:
            widget.SetToolTip(text)
        except Exception:
            pass
        try:
            widget.Bind(wx.EVT_ENTER_WINDOW, self._on_hover_event)
            widget.Bind(wx.EVT_MOTION, self._on_hover_event)
        except Exception:
            pass

    def _install_hover_map(self, mapping):
        for widget, text in mapping:
            self._install_hover_help(widget, text)


    def _hover_widget_name(self, widget):
        for name, value in self.__dict__.items():
            if value is widget:
                return name
        return ''

    def _guess_hover_text(self, widget, name=''):
        name_l = str(name or '').lower()
        label = ''
        try:
            label = str(widget.GetLabel()).strip()
        except Exception:
            label = ''
        label_l = label.lower()

        label_map = {
            'p0': 'Zero-order phase correction.',
            'p1': 'First-order phase correction.',
            'lp': 'Linear prediction.',
            'poly': 'Polynomial baseline correction.',
            'f1180': 'Apply first-point correction.',
            'flip': 'Signal polarity / flip handling.',
            'window': 'Select the apodization window function.',
            'op1': 'Primary parameter for the selected window function.',
            'op2': 'Secondary parameter for the selected window function.',
            'first pt': 'First-point scaling factor used by the window.',
            'direct dimension:': 'Settings for the direct dimension.',
            'dimensions:': 'Configure nucleus labels for each dimension.',
            'referencing:': 'Referencing settings.',
            'sparse sampling:': 'Sparse-sampling schedule settings.',
            'status': 'Status summary.',
            'dimension:': 'Current data dimensionality.',
            'format:': 'Detected input data format.',
            'converted': 'Converted output status indicator.',
            'processed': 'Processed output status indicator.',
            'script target': 'Choose the target script mode.',
            'smile maxiter:': 'Maximum number of SMILE iterations.',
            'outpath:': 'Working output directory.',
            'specpath:': 'Processed spectrum output directory.',
            'fid select:': 'Choose which FID trace or slice is used for the preview.',
            'x ranges (ppm):': 'Set the displayed ppm window.',
            'min:': 'Lower ppm limit of the display window.',
            'max:': 'Upper ppm limit of the display window.',
        }
        if label_l in label_map:
            return label_map[label_l]

        class_name = type(widget).__name__.lower()

        if 'showscript' in name_l:
            return 'Preview the generated script.'
        if 'guessbtn' in name_l:
            return 'Guess the processing settings from the current file.'
        if 'runbtn' in name_l:
            return 'Run the current processing or conversion step.'
        if 'savebtn' in name_l:
            return 'Save the current settings.'
        if 'closebtn' in name_l:
            return 'Close this window.'
        if 'buttonconv' in name_l:
            return 'Open the conversion settings.'
        if 'buttonproc' in name_l:
            return 'Open the processing settings.'
        if 'opendirfilebtn' in name_l:
            return 'Choose the output directory.'
        if 'openspecdirbtn' in name_l:
            return 'Choose the spectrum output folder.'
        if 'nusbrowse' in name_l:
            return 'Browse for a sparse-sampling schedule.'
        if 'sld_0_auto_btn' in name_l or 'sld_1_auto_btn' in name_l:
            return 'Estimate the phase automatically from the current signal.'
        if 'sld_0_mode_btn' in name_l:
            return 'Switch between coarse and fine P0 control.'
        if 'sld_1_mode_btn' in name_l:
            return 'Switch between coarse and fine P1 control.'
        if 'scriptmodebox' in name_l:
            return 'Choose whether the generated script targets standard processing or SMILE.'
        if 'maxiterbox' in name_l:
            return 'Maximum number of SMILE iterations.'
        if 'refbox' in name_l:
            return 'Choose the referencing mode.'
        if 'xcenbox' in name_l:
            return 'Enter the manual reference position in ppm.'
        if 'realname' in name_l:
            return 'Enter the real-axis label used for pseudo-dimension data.'
        if 'nusfil' in name_l:
            return 'Path to the sparse-sampling schedule file.'
        if 'path_label' in name_l:
            return 'Current processing script path.'
        if 'scriptbox' in name_l:
            return 'Editable processing script text.'
        if 'convertedautobtn' in name_l:
            return 'Automatically refresh the converted status.'
        if 'processedautobtn' in name_l:
            return 'Automatically refresh the processed status.'
        if 'convertedlamp' in name_l:
            return 'Converted output status light.'
        if 'processedlamp' in name_l:
            return 'Processed output status light.'
        if 'cb_show_fid' in name_l:
            return 'Display the raw FID instead of the phased spectrum.'
        if 'fidselect' in name_l:
            return 'Choose which FID trace or slice is used for the preview.'
        if 'xminbox' in name_l:
            return 'Lower ppm limit of the display window.'
        if 'xmaxbox' in name_l:
            return 'Upper ppm limit of the display window.'

        if 'check' in class_name:
            if 'base' in name_l and 'lin' in name_l:
                return 'Enable linear baseline correction.'
            if 'base' in name_l and 'sol' in name_l:
                return 'Enable digital solvent suppression.'
            if 'base' in name_l and 'pol' in name_l:
                return 'Apply polynomial baseline correction.'
            if 'f1180' in name_l:
                return 'Apply first-point correction.'
            if 'lp' in name_l:
                return 'Enable linear prediction.'
            if 'rk' in name_l:
                return 'Enable the Rance-Kay option for this indirect dimension.'
            if 'show_fid' in name_l:
                return 'Display the raw FID instead of the phased spectrum.'
            return 'Toggle this option.'

        if 'togglebutton' in class_name:
            if 'auto' in label_l:
                return 'Automatically refresh the related status value.'
            if label_l == 'c':
                return 'Switch between coarse and fine control.'
            return 'Toggle this option.'

        if 'radiobox' in class_name:
            return 'Choose one of the available modes.'

        if 'slider' in class_name:
            if '0' in name_l:
                return 'Adjust the zero-order phase.'
            if '1' in name_l:
                return 'Adjust the first-order phase.'
            return 'Adjust the current value.'

        if 'combobox' in class_name:
            if 'windowbox' in name_l:
                return 'Choose the apodization window function.'
            if 'cb_ft' in name_l:
                return 'Choose the flip mode for this dimension.'
            if 'refbox' in name_l:
                return 'Choose the referencing mode.'
            return 'Choose a value from the list.'

        if 'textctrl' in class_name:
            if 'p0' in name_l:
                return 'Enter the zero-order phase correction.'
            if 'p1' in name_l:
                return 'Enter the first-order phase correction.'
            if 'firstpoint' in name_l:
                return 'Enter the first-point scaling factor.'
            if 'win2val' in name_l:
                return 'Enter the primary window parameter.'
            if 'win3val' in name_l:
                return 'Enter the secondary window parameter.'
            if 'label' in name_l:
                return 'Enter the nucleus label for this dimension.'
            if 'maxiterbox' in name_l:
                return 'Maximum number of SMILE iterations.'
            if 'xcenbox' in name_l:
                return 'Enter the manual reference position in ppm.'
            if 'nusfil' in name_l:
                return 'Path to the sparse-sampling schedule file.'
            if 'scriptbox' in name_l:
                return 'Editable processing script text.'
            if 'path' in name_l:
                return 'Current file or folder path.'
            return 'Enter a value.'

        if 'button' in class_name:
            if label_l == '...':
                if 'spec' in name_l:
                    return 'Choose the spectrum output folder.'
                if 'dir' in name_l:
                    return 'Choose the output directory.'
                return 'Browse for a file or folder.'
            if label_l == 'autophase':
                return 'Estimate the phase automatically from the current signal.'
            if label_l == 'load':
                return 'Load the sparse-sampling schedule file.'
            if label_l == 'advanced':
                if 'conv' in name_l:
                    return 'Open advanced conversion settings.'
                if 'proc' in name_l:
                    return 'Open advanced processing settings.'
                return 'Open advanced settings.'
            if label_l == 'guess':
                return 'Guess the processing settings from the current file.'
            return f'{label or "Button"} action.'

        if 'statictext' in class_name:
            if label_l == 'p0:':
                return 'Zero-order phase correction.'
            if label_l == 'p1:':
                return 'First-order phase correction.'
            if label_l == 'f1180':
                return 'Apply first-point correction.'
            if label_l == 'lp':
                return 'Linear prediction.'
            if label_l == 'poly':
                return 'Polynomial baseline correction.'
            if label_l == 'flip':
                return 'Signal polarity / flip handling.'
            if label_l == 'window':
                return 'Apodization window function.'
            if label_l == 'op1':
                return 'Primary parameter for the selected window function.'
            if label_l == 'op2':
                return 'Secondary parameter for the selected window function.'
            if label_l == 'first pt':
                return 'First-point scaling factor.'
            if label_l.startswith('dim '):
                return f'Processing settings for {label_l}.'
            if label_l.startswith('label '):
                return 'Enter the nucleus label for this dimension.'
            if label_l == 'direct dimension:':
                return 'Settings for the direct dimension.'
            if label_l == 'referencing:':
                return 'Referencing settings.'
            if label_l == 'sparse sampling:':
                return 'Sparse-sampling schedule settings.'
            if label_l == 'status':
                return 'Status summary.'
        return ''

    def _install_default_hover_help(self):
        try:
            items = list(self.__dict__.items())
        except Exception:
            items = []
        for name, widget in items:
            if widget is None:
                continue
            if name.startswith('_'):
                continue
            if name in {'panel', 'status_bar', 'statusBox', 'statusSizer', 'statusPanel', 'canvas', 'toolbar'}:
                continue
            try:
                if not isinstance(widget, wx.Window):
                    continue
            except Exception:
                continue
            try:
                if hasattr(self, '_hover_help_map') and widget in self._hover_help_map:
                    continue
            except Exception:
                pass
            text = self._guess_hover_text(widget, name)
            if not text:
                continue
            self._install_hover_help(widget, text)

    def _on_canvas_motion(self, event):
        try:
            if getattr(event, 'inaxes', None) is not getattr(self, 'axes', None):
                self._set_status_text(self._hover_default_status)
                return
        except Exception:
            self._set_status_text(self._hover_default_status)
            return
        try:
            if bool(self.cb_show_fid.IsChecked()):
                self._set_status_text('Showing FID')
                return
        except Exception:
            pass
        try:
            p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
            p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
            self._set_status_text(f'P0={p0:.2f}  P1={p1:.2f}')
        except Exception:
            self._set_status_text('Spectrum')

    def _on_canvas_leave(self, event):
        self._set_status_text(self._hover_default_status)

    def _build_status_panel(self):
        self.statusBox = wx.StaticBox(self, -1, 'Status')
        self.statusSizer = wx.StaticBoxSizer(self.statusBox, wx.VERTICAL)
        grid = wx.FlexGridSizer(rows=0, cols=4, vgap=6, hgap=10)
        grid.AddGrowableCol(1, 1)

        self.spectral_dim_countensionStatusLab = wx.StaticText(self.statusBox, label='Dimension:')
        self.dimensionStatusValue = wx.StaticText(self.statusBox, label='')
        self.formatStatusLab = wx.StaticText(self.statusBox, label='Format:')
        self.formatStatusValue = wx.StaticText(self.statusBox, label='')
        self.convertedLamp = self._make_status_lamp(self.statusBox, wx.Colour(200, 0, 0))
        self.convertedStatusLab = wx.StaticText(self.statusBox, label='converted')
        self.convertedAutoBtn = wx.ToggleButton(self.statusBox, label='auto', size=(55, 22))
        self.buttonConv = wx.Button(self.statusBox, label='advanced', size=(-1, 22))
        self.processedLamp = self._make_status_lamp(self.statusBox, wx.Colour(200, 0, 0))
        self.processedStatusLab = wx.StaticText(self.statusBox, label='processed')
        self.processedAutoBtn = wx.ToggleButton(self.statusBox, label='auto', size=(55, 22))
        self.buttonProc = wx.Button(self.statusBox, label='advanced', size=(-1, 22))
        self.combineBrukerBtn = wx.Button(self.statusBox, label='combine', size=(-1, 22))
        self.combineBrukerBtn.Show(self.has_pseudo_axis and self.spectral_dim_count == 1)

        grid.Add(self.spectral_dim_countensionStatusLab, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.dimensionStatusValue, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add((1, 1))
        grid.Add((1, 1))
        grid.Add(self.formatStatusLab, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.formatStatusValue, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add((1, 1))
        grid.Add(self.combineBrukerBtn, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.convertedLamp, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.convertedStatusLab, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.convertedAutoBtn, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.buttonConv, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.processedLamp, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.processedStatusLab, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.processedAutoBtn, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.buttonProc, 0, wx.ALIGN_CENTER_VERTICAL)

        self.buttonConv.Bind(wx.EVT_BUTTON, self.OnButtonConversion)
        self.buttonProc.Bind(wx.EVT_BUTTON, self.OnButtonProcessing)
        self.combineBrukerBtn.Bind(wx.EVT_BUTTON, self.OnCombineBruker)
        self.processedAutoBtn.Bind(wx.EVT_TOGGLEBUTTON, self.OnButtonProcessingAuto)
        self.convertedAutoBtn.Bind(wx.EVT_TOGGLEBUTTON, self.OnButtonConversionAuto)

        self._install_hover_map([
            (self.spectral_dim_countensionStatusLab, 'Current data dimensionality, including whether a pseudo axis is present.'),
            (self.dimensionStatusValue, 'Current data dimensionality, including whether a pseudo axis is present.'),
            (self.formatStatusLab, 'Detected input data format for this dataset.'),
            (self.formatStatusValue, 'Detected input data format for this dataset.'),
            (self.convertedLamp, 'Green means the converted output files are present in the output folder.'),
            (self.convertedStatusLab, 'Converted output status indicator.'),
            (self.convertedAutoBtn, 'Automatically refresh the converted status.'),
            (self.buttonConv, 'Open the conversion settings.'),
            (self.processedLamp, 'Green means the processed output files are present in the spectrum folder.'),
            (self.processedStatusLab, 'Processed output status indicator.'),
            (self.processedAutoBtn, 'Automatically refresh the processed status.'),
            (self.buttonProc, 'Open the processing settings.'),
            (self.combineBrukerBtn, 'Combine numbered Bruker experiment folders in the fid path into one pseudo-2D SER.'),
        ])

        self.statusPanel = self
        self.statusSizer.Add(grid, 0, wx.ALL | wx.EXPAND, 6)


    def direct_dimension_box(self):
        self.refList = ['Water', 'Auto', 'Manual']
        self.reference_mode = 'Water'

        self.dataLbl = wx.StaticBox(self, -1, 'Direct dimension:')
        self.dataSizer = wx.StaticBoxSizer(self.dataLbl, wx.VERTICAL)

        self.sld_0_lbl = wx.StaticText(self.dataLbl, label='P0:')
        self.sld_1_lbl = wx.StaticText(self.dataLbl, label='P1:')
        self.sld_0 = wx.Slider(self.dataLbl, value=1000, minValue=-36000, maxValue=36000, size=(250, -1),
                               style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        self.sld_1 = wx.Slider(self.dataLbl, value=1000, minValue=-36000, maxValue=36000, size=(250, -1),
                               style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        self.sld_0_mode_btn = wx.ToggleButton(self.dataLbl, label='C', size=(26, 22))
        self.sld_1_mode_btn = wx.ToggleButton(self.dataLbl, label='C', size=(26, 22))
        self.sld_0_auto_btn = wx.Button(self.dataLbl, label='Autophase', size=(86, 22))
        self.sld_1_auto_btn = wx.Button(self.dataLbl, label='Autophase', size=(86, 22))
        self.sld_0_mode_btn.SetValue(False)
        self.sld_1_mode_btn.SetValue(False)
        self.sld_0_mode_btn.SetLabel('C')
        self.sld_1_mode_btn.SetLabel('C')

        self.cb_show_fid = wx.CheckBox(self.dataLbl, -1, '')
        self.cb_show_fid.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.on_show_fid, self.cb_show_fid)

        self.xLab = wx.StaticText(self.dataLbl, label='X ranges (ppm):')
        self.xminLab = wx.StaticText(self.dataLbl, label='Min:')
        self.xmaxLab = wx.StaticText(self.dataLbl, label='Max:')
        # TE_PROCESS_ENTER lets the ppm range fields apply immediately when
        # Enter is pressed, without stealing the Y zoom selected by the user.
        self.xmaxBox = wx.TextCtrl(self.dataLbl, size=(40, 22), style=wx.TE_PROCESS_ENTER)
        self.xminBox = wx.TextCtrl(self.dataLbl, size=(40, 22), style=wx.TE_PROCESS_ENTER)
        self.xminBox.Bind(wx.EVT_TEXT_ENTER, self.on_xrange_enter)
        self.xmaxBox.Bind(wx.EVT_TEXT_ENTER, self.on_xrange_enter)
        self.autoCentreBtn = wx.Button(self.dataLbl, label='AutoCentre', size=(86, 22))
        self.autoAmideBtn = wx.Button(self.dataLbl, label='AutoAmide', size=(86, 22))
        self.autoCentreBtn.Bind(wx.EVT_BUTTON, self.OnAutoCentre)
        self.autoAmideBtn.Bind(wx.EVT_BUTTON, self.OnAutoAmide)

        self.sld_0.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
        self.sld_1.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
        self.sld_0_mode_btn.Bind(wx.EVT_TOGGLEBUTTON, self.OnP0RangeToggle)
        self.sld_1_mode_btn.Bind(wx.EVT_TOGGLEBUTTON, self.OnP1RangeToggle)
        self.sld_0_auto_btn.Bind(wx.EVT_BUTTON, self.OnP0AutoPhase)

        self.sizer = wx.GridBagSizer(10, 8)
        self.phaseSizer = wx.BoxSizer(wx.VERTICAL)

        self.p0Row = wx.BoxSizer(wx.HORIZONTAL)
        self.p0Row.Add(self.sld_0_lbl, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.p0Row.Add(self.sld_0, 1, flag=wx.EXPAND)
        self.p0Row.Add(self.sld_0_mode_btn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)
        self.p0Row.Add(self.sld_0_auto_btn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)

        self.p1Row = wx.BoxSizer(wx.HORIZONTAL)
        self.p1Row.Add(self.sld_1_lbl, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.p1Row.Add(self.sld_1, 1, flag=wx.EXPAND)
        self.p1Row.Add(self.sld_1_mode_btn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)
        self.p1Row.Add(self.sld_1_auto_btn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)

        self.fidSelectLab = wx.StaticText(self.dataLbl, label='FID select:')
        self.fidSelectOptions = [str(i) for i in range(1, self._fid_selection_count() + 1)]
        self.fidSelect = wx.Choice(self.dataLbl, choices=self.fidSelectOptions)
        self.Bind(wx.EVT_CHOICE, self.on_fid_select, self.fidSelect)
        self._set_fid_selection(self._normalize_fid_selection(getattr(self, 'FIDsel', 1)))

        self.fidSelectSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fidSelectSizer.Add(self.fidSelectLab, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=6)
        self.fidSelectSizer.Add(self.fidSelect, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.fidSelectSizer.Add(self.cb_show_fid, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        self.phaseSliceModeLab = wx.StaticText(self.dataLbl, label='Pseudo-2D phasing:')
        self.phaseSliceMode = wx.Choice(self.dataLbl, choices=['First', 'Summed'])
        self.phaseSliceMode.SetSelection(0)
        self.phaseSliceMode.Bind(wx.EVT_CHOICE, self.on_phase_slice_mode)
        self.phaseSliceModeSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phaseSliceModeSizer.Add(self.phaseSliceModeLab, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=6)
        self.phaseSliceModeSizer.Add(self.phaseSliceMode, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        pseudo2d = bool(getattr(self, 'has_pseudo_axis', False) and self._spectral_dimension_count() == 1)
        self.phaseSliceModeLab.Show(pseudo2d)
        self.phaseSliceMode.Show(pseudo2d)

        self.xRangeRow = wx.BoxSizer(wx.HORIZONTAL)
        self.xRangeRow.Add(self.xLab, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=8)
        self.xRangeRow.Add(self.xminLab, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=4)
        self.xRangeRow.Add(self.xminBox, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=10)
        self.xRangeRow.Add(self.xmaxLab, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, border=4)
        self.xRangeRow.Add(self.xmaxBox, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        # Keep the automatic range controls at the right of the row, aligned
        # with the automatic phase controls above as the frame expands.
        self.xRangeRow.AddStretchSpacer(1)
        self.xRangeRow.Add(self.autoCentreBtn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)
        self.xRangeRow.Add(self.autoAmideBtn, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=8)

        self.phaseSizer.Add(self.p0Row, 0, wx.EXPAND)
        self.phaseSizer.Add(self.p1Row, 0, wx.TOP | wx.EXPAND, 6)
        self.phaseSizer.Add(self.xRangeRow, 0, wx.TOP, 6)
        self.phaseSizer.Add(self.fidSelectSizer, 0, wx.TOP, 6)
        self.phaseSizer.Add(self.phaseSliceModeSizer, 0, wx.TOP, 6)
        self.sizer.Add(self.phaseSizer, (0, 0), flag=wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)

        self._install_hover_map([
            (self.sld_0_lbl, 'Adjust zeroth-order phase for the direct dimension.'),
            (self.sld_1_lbl, 'Adjust first-order phase for the direct dimension.'),
            (self.sld_0, 'Adjust zeroth-order phase for the direct dimension.'),
            (self.sld_1, 'Adjust first-order phase for the direct dimension.'),
            (self.sld_0_mode_btn, 'Switch between coarse and fine P0 control.'),
            (self.sld_1_mode_btn, 'Switch between coarse and fine P1 control.'),
            (self.sld_0_auto_btn, 'Estimate P0 automatically from the direct-dimension signal.'),
            (self.sld_1_auto_btn, 'Estimate P1 automatically from the direct-dimension signal.'),
            (self.cb_show_fid, 'Display the raw FID instead of the phased spectrum.'),
            (self.fidSelectLab, 'Choose which FID trace or slice is used for the preview.'),
            (self.fidSelect, 'Choose which FID trace or slice is used for the preview.'),
            (self.xLab, 'Set the displayed ppm window.'),
            (self.xminLab, 'Lower ppm limit of the display window.'),
            (self.xmaxLab, 'Upper ppm limit of the display window.'),
            (self.xminBox, 'Lower ppm limit of the display window.'),
            (self.xmaxBox, 'Upper ppm limit of the display window.'),
            (self.autoCentreBtn, 'Automatically zoom to the full 1H signal region, including water, with baseline padding.'),
            (self.autoAmideBtn, 'Automatically zoom to the amide 1H region (5.0 ppm and above), excluding water, with baseline padding.'),
        ])

        self._set_phase_slider_mode(self.sld_0, self.sld_0_mode_btn, False, current_value=self._phase_slider_value_for_mode(self.sld_0, fine=False))
        self._set_phase_slider_mode(self.sld_1, self.sld_1_mode_btn, False, current_value=self._phase_slider_value_for_mode(self.sld_1, fine=False))

        self.border = wx.BoxSizer(wx.VERTICAL)
        self.border.Add(self.sizer, 1, wx.ALL | wx.EXPAND, 7)
        self.border.AddSpacer(5)
        self.dataSizer.Add(self.border, flag=wx.ALIGN_CENTER_HORIZONTAL)

    def FindData(self):

        from pathlib import Path
        for test in ('fid', 'ser', 'fid.gz', 'ser.gz'):
            for path in Path('./').rglob(test):
                return str(path).split(test)[0]
        return ''

    def _phase_slider_scale(self, button=None):
        return 100

    def _phase_slider_value_for_mode(self, slider, fine=False):
        try:
            raw = float(slider.GetValue())
        except Exception:
            raw = 0.0
        scale = 100.0
        return raw / scale

    def _phase_slider_value(self, slider, button=None):
        try:
            raw = float(slider.GetValue())
        except Exception:
            raw = 0.0
        scale = self._phase_slider_scale(button)
        return raw / float(scale)

    def _set_phase_slider_value(self, slider, button, value):
        try:
            scale = self._phase_slider_scale(button)
            slider.SetValue(int(round(float(value) * float(scale))))
        except Exception:
            pass

    def _set_phase_slider_mode(self, slider, button, fine=False, current_value=None):
        try:
            if current_value is None:
                current_value = self._phase_slider_value_for_mode(slider, fine=not fine)
            cur = float(current_value)
            span = 10.0 if fine else 180.0
            scale = 100
            slider.SetRange(int(round((cur - span) * scale)), int(round((cur + span) * scale)))
            try:
                slider.SetLineSize(1)
                slider.SetPageSize(1 if fine else 10)
            except Exception:
                pass
            button.SetLabel('F' if fine else 'C')
            self._set_phase_slider_value(slider, button, cur)
            pass
        except Exception:
            pass

    def OnP0RangeToggle(self, event):
        new_fine = bool(self.sld_0_mode_btn.GetValue())
        current = self._phase_slider_value_for_mode(self.sld_0, fine=not new_fine)
        self._set_phase_slider_mode(self.sld_0, self.sld_0_mode_btn, new_fine, current_value=current)
        self._set_direct_phase_values(self._phase_slider_value(self.sld_0, self.sld_0_mode_btn), self._phase_slider_value(self.sld_1, self.sld_1_mode_btn))

    def OnP1RangeToggle(self, event):
        new_fine = bool(self.sld_1_mode_btn.GetValue())
        current = self._phase_slider_value_for_mode(self.sld_1, fine=not new_fine)
        self._set_phase_slider_mode(self.sld_1, self.sld_1_mode_btn, new_fine, current_value=current)
        self._set_direct_phase_values(self._phase_slider_value(self.sld_0, self.sld_0_mode_btn), self._phase_slider_value(self.sld_1, self.sld_1_mode_btn))

    def _autophase_direct_region(self):
        """Return the unphased direct-dimension trace used for autophase."""
        dic, data = self._build_direct_frequency_data(include_phased=False)
        dic, data, uc = self._extract_direct_region(dic, data)
        return dic, data, uc

    def OnP0AutoPhase(self, event):
        if getattr(self, '_autophase_in_progress', False):
            return

        show_fid = False
        try:
            show_fid = bool(self.cb_show_fid.IsChecked())
        except Exception:
            show_fid = False
        if show_fid:
            # Autophase only applies to the frequency-domain spectrum.
            return

        self._autophase_in_progress = True
        try:
            _, data, _ = self._autophase_direct_region()
            current_p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
            current_p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
            best_p0, best_score = _autophase_grid_search(data, p1=current_p1, initial=current_p0)
            pass

            # Update the slider and text controls only after the search completes.
            self._set_direct_phase_values(best_p0, current_p1, sync_sliders=True)
            self.draw_figure()
        except Exception:
            pass
            errorMessage('Could not determine an automatic phase.')
        finally:
            self._autophase_in_progress = False

    def _direct_apodization_preview(self, g1, g2):
        """Return a processed 1D spectrum for candidate GM apodization values.

        This uses the same direct-dimension data path as the main processing run
        but keeps the GUI untouched so the grid search can evaluate candidates
        without redrawing.
        """
        pipe_proc = self._pipe_proc_module()
        if pipe_proc is None:
            raise RuntimeError('nmrglue.pipe_proc is not available')

        if not hasattr(self, 'arr') or getattr(self, 'arr', None) is None:
            if not self._load_direct_fid():
                raise RuntimeError('Could not load the direct-dimension FID')

        dic = dict(getattr(self, 'arrdic', {}) or {})
        data = numpy.array(self._first_fid_slice(self.arr), copy=True)

        if self._ctrl_checked('cb_baseSol', False):
            dic, data = pipe_proc.sol(dic, data)

        _, _, _, c = self._direct_window_args()
        dic, data = pipe_proc.apod(dic, data, qName='GM', q1=float(g1), q2=float(g2), c=float(c))
        dic, data = pipe_proc.zf(dic, data, zf=2)
        dic, data = pipe_proc.ft(dic, data, **self._direct_ft_flags())

        try:
            p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
            p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
        except Exception:
            p0 = 0.0
            p1 = 0.0
        try:
            data = ng.proc_base.ps(data, p0=float(p0), p1=float(p1))
        except Exception:
            pass

        dic, data, uc = self._extract_direct_region(dic, data)
        x = numpy.asarray(uc.ppm_scale())
        y_real = -numpy.real(numpy.asarray(data))
        y_imag = numpy.imag(numpy.asarray(data))
        return x, y_real, y_imag

    def _apodization_grid_candidates(self, current_g1, current_g2):
        current_g1 = max(float(current_g1), 0.0)
        current_g2 = max(float(current_g2), 0.0)
        g1_max = max(6.0, current_g1 * 2.0, 4.0)
        g2_max = max(30.0, current_g2 * 2.0, 20.0)
        g1_vals = numpy.linspace(0.0, g1_max, 7)
        g2_vals = numpy.linspace(0.0, g2_max, 7)
        return [(float(g1), float(g2)) for g1 in g1_vals for g2 in g2_vals if g2 + 1e-12 >= g1]

    def _refine_apodization_candidates(self, best_g1, best_g2, g1_step, g2_step):
        g1_vals = sorted({
            max(0.0, float(best_g1) - float(g1_step)),
            max(0.0, float(best_g1) - float(g1_step) / 2.0),
            float(best_g1),
            float(best_g1) + float(g1_step) / 2.0,
            float(best_g1) + float(g1_step),
        })
        g2_vals = sorted({
            max(0.0, float(best_g2) - float(g2_step)),
            max(0.0, float(best_g2) - float(g2_step) / 2.0),
            float(best_g2),
            float(best_g2) + float(g2_step) / 2.0,
            float(best_g2) + float(g2_step),
        })
        return [(float(g1), float(g2)) for g1 in g1_vals for g2 in g2_vals if g2 + 1e-12 >= g1]

    def _autoapodise_report(self, progress_frame, message, status=None):
        try:
            pass
        except Exception:
            pass
        if progress_frame is not None:
            try:
                progress_frame.append_text(message.rstrip() + '\n')
            except Exception:
                pass
            if status is not None:
                try:
                    progress_frame.set_status(status)
                except Exception:
                    pass
            try:
                wx.YieldIfNeeded()
            except Exception:
                try:
                    wx.Yield()
                except Exception:
                    pass

    def _evaluate_apodization_candidates(self, candidates, natural_fwhm=None, progress_frame=None, stage='coarse'):
        evaluations = []
        total = len(candidates)
        for idx, (g1, g2) in enumerate(candidates, start=1):
            try:
                x, y_real, y_imag = self._direct_apodization_preview(g1, g2)
                metrics = _apodization_candidate_score(x, y_real, y_imag, natural_fwhm=natural_fwhm)
                metrics.update({'g1': float(g1), 'g2': float(g2), 'x': x, 'y_real': y_real, 'y_imag': y_imag})
                evaluations.append(metrics)
                pass
                if progress_frame is not None and (idx == 1 or idx == total or idx % max(1, total // 10) == 0):
                    self._autoapodise_report(
                        progress_frame,
                        f'[{stage}] {idx}/{total} g1={g1:.3f} g2={g2:.3f} score={metrics["score"]:.6f}',
                        status=f'{stage.capitalize()} search {idx}/{total}'
                    )
            except Exception:
                pass
                if progress_frame is not None:
                    self._autoapodise_report(progress_frame, f'[{stage}] failed g1={g1:.3f} g2={g2:.3f}')
        return evaluations

    def OnAutoApodise(self, event=None, progress_frame=None):
        if getattr(self, '_autoapodise_in_progress', False):
            return
        if getattr(self, 'ind_phase', False):
            errorMessage('Automatic apodization is only available for the direct dimension.')
            return

        if progress_frame is None:
            try:
                progress_frame = ShellOutputFrame(self, title='AutoApodise Progress')
                progress_frame.set_status('Starting automatic apodization...')
                progress_frame.append_text('Starting automatic apodization...\n')
                progress_frame.Show()
            except Exception:
                import traceback
                pass
                progress_frame = None

        self._autoapodise_in_progress = True
        final_run_started = False
        try:
            self._autoapodise_report(progress_frame, 'Preparing automatic apodization search...', 'Preparing spectrum')
            current_window = self._ctrl_value('windowBox0', default='GM') or 'GM'
            if current_window != 'GM':
                self._set_processing_widget_value('windowBox0', 'GM')
                try:
                    self.SetWin(
                        self._processing_widget('windowBox0'), self._processing_widget('win3Val0'),
                        self._processing_widget('win2Val0'), self._processing_widget('firstPoint0')
                    )
                except Exception:
                    pass

            # Use the currently processed spectrum to identify the strongest peak
            # and its width, then search around the current GM settings for the
            # best compromise between sharpening and artefact suppression.
            self._autoapodise_report(progress_frame, 'Extracting current spectrum and ROI...', 'Extracting current spectrum')
            dic, data = self._build_direct_frequency_data(include_phased=True)
            dic, data, uc = self._extract_direct_region(dic, data)
            x_cur = numpy.asarray(uc.ppm_scale())
            y_cur = -numpy.real(numpy.asarray(data))
            if x_cur.size < 8:
                raise RuntimeError('Not enough points in the current spectrum to optimise apodization.')

            peak_idx = int(numpy.nanargmax(numpy.abs(y_cur)))
            half_width = max(8, int(round(x_cur.size * 0.12)))
            lo = max(0, peak_idx - half_width)
            hi = min(x_cur.size, peak_idx + half_width + 1)
            x_roi = x_cur[lo:hi]
            y_roi = y_cur[lo:hi]
            if x_roi.size < 8:
                x_roi = x_cur
                y_roi = y_cur
            self._autoapodise_report(progress_frame, f'ROI length={x_roi.size} points; peak index={peak_idx}', 'Estimating linewidth')

            # A quick estimate of the natural linewidth is the narrowest well-fit
            # peak from a broad coarse scan; later we penalise candidates that are
            # broader than this by a comfortable margin.
            cur_g1 = float(self._ctrl_value('win3Val0', default='2.0') or 2.0)
            cur_g2 = float(self._ctrl_value('win2Val0', default='20.0') or 20.0)
            coarse_candidates = self._apodization_grid_candidates(cur_g1, cur_g2)
            self._autoapodise_report(progress_frame, f'Running coarse grid search with {len(coarse_candidates)} candidates...', 'Coarse search')
            coarse_eval = self._evaluate_apodization_candidates(coarse_candidates, natural_fwhm=None, progress_frame=progress_frame, stage='coarse')
            valid_widths = [e['fwhm'] for e in coarse_eval if numpy.isfinite(e.get('fwhm', numpy.nan)) and e.get('fwhm', 0.0) > 0.0]
            natural_fwhm = float(numpy.percentile(valid_widths, 20)) if valid_widths else None
            if natural_fwhm is not None:
                self._autoapodise_report(progress_frame, f'Estimated natural linewidth FWHM = {natural_fwhm:.6f}', 'Estimating natural linewidth')
            else:
                self._autoapodise_report(progress_frame, 'Could not estimate a stable natural linewidth from the coarse scan.', 'Estimating natural linewidth')

            # Refine around the best coarse candidate with a smaller grid.
            best_coarse = min(coarse_eval, key=lambda e: e['score']) if coarse_eval else None
            if best_coarse is None:
                raise RuntimeError('Could not evaluate any apodization candidates.')

            g1_step = max(0.25, (max(6.0, cur_g1 * 2.0, 4.0) / 6.0))
            g2_step = max(1.0, (max(30.0, cur_g2 * 2.0, 20.0) / 6.0))
            refined_candidates = self._refine_apodization_candidates(best_coarse['g1'], best_coarse['g2'], g1_step, g2_step)
            self._autoapodise_report(progress_frame, f'Running fine grid search with {len(refined_candidates)} candidates...', 'Fine search')
            refined_eval = self._evaluate_apodization_candidates(refined_candidates, natural_fwhm=natural_fwhm, progress_frame=progress_frame, stage='fine')
            all_eval = coarse_eval + refined_eval
            if not all_eval:
                raise RuntimeError('Could not optimise apodization: no valid candidate spectra were produced.')

            best = min(all_eval, key=lambda e: e['score'])
            self._autoapodise_report(
                progress_frame,
                'Best candidate: g1={:.3f}, g2={:.3f}, score={:.6f}, fwhm={}'.format(
                    best['g1'], best['g2'], best['score'],
                    '{:.6f}'.format(best['fwhm']) if numpy.isfinite(best.get('fwhm', numpy.nan)) else 'nan'
                ),
                'Applying best apodization'
            )

            # Update the apodization controls only once the search is complete.
            self._set_processing_widget_value('windowBox0', 'GM')
            self._set_processing_widget_value('win3Val0', '{:.6g}'.format(best['g1']))
            self._set_processing_widget_value('win2Val0', '{:.6g}'.format(best['g2']))
            try:
                self.SetWin(
                    self._processing_widget('windowBox0'), self._processing_widget('win3Val0'),
                    self._processing_widget('win2Val0'), self._processing_widget('firstPoint0')
                )
            except Exception:
                pass

            self._autoapodise_report(progress_frame, 'Re-running processing with selected apodization...', 'Reprocessing spectrum')

            # Re-run the processing script once with the selected apodization so
            # the plot and downstream outputs update together.
            def _finish_autoapodise():
                try:
                    self.draw_figure()
                finally:
                    self._autoapodise_in_progress = False
                    if progress_frame is not None:
                        try:
                            progress_frame.set_status('AutoApodise complete')
                            progress_frame.append_text('AutoApodise complete\n')
                        except Exception:
                            pass

            self._run_processing_auto(on_finish=_finish_autoapodise, output_frame=progress_frame)
            final_run_started = True
        except Exception:
            import traceback
            pass
            pass
            errorMessage('Could not determine an automatic apodization setting.')
            if progress_frame is not None:
                try:
                    progress_frame.append_text('Automatic apodization failed.\n')
                    progress_frame.set_status('AutoApodise failed')
                except Exception:
                    import traceback
                    pass
            if not final_run_started:
                self._autoapodise_in_progress = False


    AutoApodise = OnAutoApodise

    def button_box(self):
        self.controls_lbl = wx.StaticBox(self, -1, 'Controls:', size=(240, 140))
        self.controls_sizer = wx.StaticBoxSizer(self.controls_lbl, wx.VERTICAL)

        self.buttonProj = wx.Button(self.controls_lbl, label='Projections', size=(-1, 22))
        self.buttonClean = wx.Button(self.controls_lbl, label='Clean', size=(-1, 22))
        self.projectionsLabel = wx.StaticText(self.controls_lbl, label='Projections:')
        self.projections = wx.ComboBox(
            self.controls_lbl, choices=['sum', 'skyline'], value='skyline',
            style=wx.CB_READONLY, size=(100, 22)
        )
        self.sizer2c = wx.GridBagSizer(3, 7)
        self.sizer2c.Add((1, 1), (0, 0))
        self.sizer2c.Add(self.buttonProj, (1, 1), border=20, flag=wx.RIGHT | wx.LEFT)
        self.sizer2c.Add(self.buttonClean, (1, 2), border=20, flag=wx.RIGHT | wx.LEFT)
        self.sizer2c.Add(self.projectionsLabel, (1, 4), border=6, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL)
        self.sizer2c.Add(self.projections, (1, 5), border=20, flag=wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.controls_sizer.Add(self.sizer2c, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.buttonProj.Bind(wx.EVT_BUTTON, self.OnButtonProjections)
        self.buttonClean.Bind(wx.EVT_BUTTON, self.OnButtonClean)

        self._install_hover_map([
            (self.buttonProj, 'Generate or refresh projection plots.'),
            (self.buttonClean, 'Remove or reset temporary processing outputs.'),
        ])

    def OnButtonConversion(self, event):
        self.sync_current_gui_state()
        if getattr(self, 'conv_frame', None) is not None:
            try:
                self.conv_frame.Raise()
                self.conv_frame.SetFocus()
                return
            except Exception:
                self.conv_frame = None
        from spinDecon.gui.dialogs.processing.conversion import ConversionFrame
        self.conv_frame = ConversionFrame(self)
        self.conv_frame.Show(True)

    def OnButtonConversionAuto(self, event):
        btn = getattr(self, 'convertedAutoBtn', None)
        if btn is not None:
            try:
                btn.SetValue(True)
                btn.Refresh()
                btn.Update()
                try:
                    wx.Yield()
                except Exception:
                    pass
            except Exception:
                pass
        def _auto_done(*_args, **_kwargs):
            if btn is not None:
                try:
                    btn.SetValue(False)
                    btn.Refresh()
                    btn.Update()
                except Exception:
                    pass

        try:
            self.MakeConvScript(on_finish=_auto_done)
        except Exception as exc:
            _auto_done()
            # Do not silently swallow Auto-conversion failures.  At this point
            # no subprocess output window may exist yet, so use the normal GUI
            # error channel and leave the traceback available to developers.
            import traceback
            traceback.print_exc()
            errorMessage('Automatic conversion could not be started: %s' % exc)
        if event is not None:
            event.Skip()

    def _ensure_processing_frame(self, show=False, reload_from_file=False):
        frame = getattr(self, 'processing_frame', None)
        if frame is not None:
            try:
                if reload_from_file:
                    loader = getattr(frame, 'reload_from_file', None)
                    if callable(loader):
                        loader()
                    elif hasattr(frame, '_load_from_file'):
                        frame._load_from_file()
                return frame
            except Exception:
                self.processing_frame = None
        from spinDecon.gui.dialogs.processing.settings import ProcessingFrame
        frame = ProcessingFrame(self)
        self.processing_frame = frame
        if reload_from_file:
            try:
                loader = getattr(frame, 'reload_from_file', None)
                if callable(loader):
                    loader()
                elif hasattr(frame, '_load_from_file'):
                    frame._load_from_file()
            except Exception:
                pass
        if not show:
            try:
                frame.Hide()
            except Exception:
                pass
        return frame

    def _run_process_script_silently(self, script_path, on_finish=None, output_frame=None):
        import subprocess
        import threading

        pipefile = self._process_pipefile()
        self.ResetReads()

        self._spec_output_dir()

        def _finish(*_args, **_kwargs):
            try:
                self.DoProjections(pipefile)
            except Exception:
                pass
            try:
                self.SetLab(refresh=False)
            except Exception:
                pass
            try:
                self._update_nmrpipe_file_box(pipefile)
            except Exception:
                pass
            try:
                self._maybe_invert_test_ft_after_processing(script_path=script_path, pipefile=pipefile)
            except Exception:
                pass
            try:
                self.UpdateLampLights()
            except Exception:
                pass
            if on_finish is not None:
                try:
                    on_finish()
                except Exception:
                    pass

        def worker():
            try:
                proc = subprocess.Popen(['csh', script_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                if proc.stdout is not None:
                    for _line in proc.stdout:
                        if output_frame is not None:
                            try:
                                wx.CallAfter(output_frame.append_text, _line)
                            except Exception:
                                pass
                proc.wait()
                try:
                    self._maybe_invert_test_ft_after_processing(script_path=script_path, pipefile=pipefile)
                except Exception:
                    pass
                if output_frame is not None:
                    try:
                        wx.CallAfter(output_frame.append_text, f'\n[process exited with code {proc.returncode}]\n')
                        wx.CallAfter(output_frame.set_status, f'Complete (exit code {proc.returncode})')
                    except Exception:
                        pass
                wx.CallAfter(_finish)
            except Exception:
                pass
                if output_frame is not None:
                    try:
                        wx.CallAfter(output_frame.append_text, '\n[error] processing script failed\n')
                        wx.CallAfter(output_frame.set_status, 'Failed')
                    except Exception:
                        pass

        threading.Thread(target=worker, daemon=True).start()
        return script_path

    def _run_processing_auto(self, on_finish=None, output_frame=None):
        # External processing must see the current widgets, not stale disk state.
        self.save_current_gui_state(reason='processing-auto')
        frame = None
        temp_frame = False
        try:
            self._spec_output_dir()
            frame = self._ensure_processing_frame(show=False, reload_from_file=False)
            try:
                temp_frame = not bool(frame.IsShown())
            except Exception:
                temp_frame = False
            try:
                self._bind_processing_controls(frame)
            except Exception:
                pass
            try:
                # Regenerate the script from the freshly reloaded controls so
                # the run uses the updated phase values from the projection window.
                guess = getattr(frame, 'on_guess', None)
                if callable(guess):
                    guess(None)
                else:
                    refresher = getattr(frame, 'refresh_from_target', None)
                    if callable(refresher):
                        refresher()
            except Exception:
                pass
            try:
                lp = frame._current_lp_flag()
            except Exception:
                lp = self._processing_target_value()
                lp = 'm' if str(lp).strip().lower() in {'mddnmr', 'mdd', 'm', '2'} else ('y' if str(lp).strip().lower() in {'smile', 'y', 'yes', '1', 'lp'} else 'n')
            try:
                frame._sync_ncpus_to_process_frame()
            except Exception:
                pass
            try:
                frame._refresh_status_text()
            except Exception:
                pass
            script_path = frame._script_target_path(lp)
            try:
                # Preserve the asynchronous completion callback supplied by the
                # caller (e.g. the ProcessProjections window).  The processing
                # window's public event handler has no callback parameter, so
                # pass it explicitly into the script runner here.
                frame.on_run(None, on_finish=on_finish)
            except Exception:
                pass
                raise
            return script_path
        finally:
            if temp_frame and frame is not None:
                try:
                    if getattr(self, 'processing_frame', None) is frame:
                        self._unbind_processing_controls()
                except Exception:
                    import traceback
                    pass
                try:
                    frame.Destroy()
                except Exception:
                    import traceback
                    pass

    def OnButtonProcessingAuto(self, event):
        btn = getattr(self, 'processedAutoBtn', None)
        if btn is not None:
            try:
                btn.SetValue(True)
                btn.Refresh()
                btn.Update()
                try:
                    wx.Yield()
                except Exception:
                    pass
            except Exception:
                pass
        try:
            self._run_processing_auto(on_finish=lambda: None)
        except Exception:
            pass
        finally:
            if btn is not None:
                try:
                    btn.SetValue(False)
                    btn.Refresh()
                    btn.Update()
                except Exception:
                    pass
        if event is not None:
            event.Skip()

    def _pipe_proc_module(self):
        pipe_proc = getattr(ng, 'pipe_proc', None)
        if pipe_proc is not None:
            return pipe_proc
        process_mod = getattr(ng, 'process', None)
        return getattr(process_mod, 'pipe_proc', None)

    def _direct_phase_backend(self):
        backend = str(self._parse_param('directPhaseBackend', default=globals().get('DIRECT_PHASE_BACKEND', 'glue'))).strip().lower()
        if backend not in ('pipe', 'glue'):
            pass
            return 'glue'
        return backend

    def _spectral_dimension_count(self):
        dim = getattr(self.parent, 'dim', None)
        try:
            return max(1, int(dim))
        except Exception:
            dim = getattr(self, 'dim', 1)
            if isinstance(dim, str) and dim.endswith('p') and dim[:-1].isdigit():
                return max(1, int(dim[:-1]) - 1)
            try:
                return max(1, int(dim))
            except Exception:
                return 1

    def _fid_selection_count(self):
        spectral_dim_count = self._spectral_dimension_count()
        if spectral_dim_count < 2:
            return 1
        return 2 ** (spectral_dim_count - 1)

    def _normalize_fid_selection(self, selection=None):
        max_sel = self._fid_selection_count()
        if max_sel <= 1:
            return 1
        if selection is None:
            selection = getattr(self, 'FIDsel', 1)
        try:
            sel = int(float(selection))
        except Exception:
            sel = 1
        if sel < 1:
            sel = 1
        if sel > max_sel:
            sel = max_sel
        return sel

    def _fid_selection_indices(self, selection=None):
        spectral_dim_count = self._spectral_dimension_count()
        lead_axes = max(0, spectral_dim_count - 1)
        if lead_axes == 0:
            return ()
        sel = self._normalize_fid_selection(selection) - 1
        return tuple((sel >> axis) & 1 for axis in range(lead_axes))

    def _set_fid_selection(self, selection=None, redraw=False):
        sel = self._normalize_fid_selection(selection)
        self.FIDsel = sel
        if hasattr(self, 'fidSelect'):
            try:
                self.fidSelect.SetSelection(sel - 1)
            except Exception:
                pass
        if getattr(self, 'state', None) is not None:
            try:
                self.state.fid_selection = sel
                self.state.metadata['FIDsel'] = sel
            except Exception:
                pass
        if redraw:
            try:
                self.draw_figure()
            except Exception:
                pass

    def _current_fid_selection(self):
        if hasattr(self, 'fidSelect'):
            try:
                return int(self.fidSelect.GetSelection()) + 1
            except Exception:
                pass
        return self._normalize_fid_selection(getattr(self, 'FIDsel', 1))

    def _processing_target_value(self):
        frame = getattr(self, 'processing_frame', None)
        if frame is not None:
            try:
                return str(frame._script_target_value())
            except Exception:
                pass
        try:
            if frame is not None:
                return {'y': 'SMILE', 'm': 'MDDNMR'}.get(frame._current_lp_flag(), 'Process')
        except Exception:
            pass
        return 'Process'

    def _first_fid_slice(self, data):
        arr = numpy.asarray(data)
        if arr.ndim <= 1:
            return arr
        lead_axes = arr.ndim - 1
        indices = list(self._fid_selection_indices())
        if len(indices) < lead_axes:
            indices.extend([0] * (lead_axes - len(indices)))
        elif len(indices) > lead_axes:
            indices = indices[:lead_axes]
        return numpy.asarray(arr[tuple(indices) + (slice(None),)])

    def _debug_print_fid_stack(self, infile, data):
        try:
            arr = numpy.asarray(data)
        except Exception:
            pass
            return

        pass
        pass
        if arr.ndim <= 1:
            absvals = numpy.abs(arr)
            pass
            return

        for idx in range(arr.shape[0]):
            slice_arr = numpy.asarray(arr[idx])
            while slice_arr.ndim > 1:
                slice_arr = slice_arr[0]
            absvals = numpy.abs(slice_arr)
            absmax = float(numpy.max(absvals)) if absvals.size else 0.0
            absmean = float(numpy.mean(absvals)) if absvals.size else 0.0
            abssum = float(numpy.sum(absvals)) if absvals.size else 0.0
            nz = int(numpy.count_nonzero(absvals))
            preview = numpy.array2string(absvals[:8], precision=6, separator=', ')
            pass
        try:
            sys.stdout.flush()
        except Exception:
            pass

    def _direct_fid_candidates(self, show_fid=False):
        import glob
        base = self._spec_output_dir()
        if not base:
            return []

        raw_root = os.path.join(base, 'raw')
        roots = [base, raw_root]

        def add_exact(cands, *parts):
            for root in roots:
                cands.append(os.path.join(root, *parts))

        def add_glob(cands, *parts):
            for root in roots:
                cands.extend(sorted(glob.glob(os.path.join(root, *parts))))

        candidates = []
        add_exact(candidates, 'slice.fid')
        add_exact(candidates, 'slice.fid.gz')

        # Retain legacy fallbacks for older layouts and debugging.
        if show_fid or self.physical_dim_count >= 3:
            if self.spectral_dim_count == 4:
                add_exact(candidates, 'fids', 'test001001.fid')
                add_exact(candidates, 'fids', 'test001001.fid.gz')
                add_glob(candidates, 'fids', 'test*.fid')
                add_glob(candidates, 'fids', 'test*.fid.gz')
            else:
                add_exact(candidates, 'fids', 'test001.fid')
                add_exact(candidates, 'fids', 'test001.fid.gz')
                add_glob(candidates, 'fids', 'test*.fid')
                add_glob(candidates, 'fids', 'test*.fid.gz')
        else:
            add_exact(candidates, 'test.fid')
            add_exact(candidates, 'test.fid.gz')
            add_exact(candidates, 'fid')
            add_exact(candidates, 'fid.gz')
            add_exact(candidates, 'ser')
            add_exact(candidates, 'ser.gz')
            add_glob(candidates, 'fids', 'test*.fid')
            add_glob(candidates, 'fids', 'test*.fid.gz')

        add_exact(candidates, 'test.fid')
        add_exact(candidates, 'test.fid.gz')
        add_exact(candidates, 'fid')
        add_exact(candidates, 'fid.gz')
        add_exact(candidates, 'ser')
        add_exact(candidates, 'ser.gz')

        seen = set()
        ordered = []
        for item in candidates:
            if item and item not in seen:
                seen.add(item)
                ordered.append(item)
        return ordered

    def _load_direct_fid(self, show_fid=False):
        self.READ1D = 0
        self.arrdic = {}
        self.arr = numpy.zeros(256, dtype=complex)
        self.index = numpy.arange(self.arr.shape[0])
        self.direct_fid_path = None
        for infile in self._direct_fid_candidates(show_fid=show_fid):
            if not infile or not os.path.exists(infile):
                continue
            try:
                dic, data = ng.pipe.read(infile)
                self._debug_print_fid_stack(infile, data)
                data = self._first_fid_slice(data)
                self.arrdic = dic
                self.arr = numpy.asarray(data)
                self.index = numpy.arange(self.arr.shape[-1])
                self.direct_fid_path = infile
                self.READ1D = 1
                return True
            except Exception:
                pass
        return False

    def _direct_dimension_is_proton(self):
        """Return True when the displayed/direct axis is labelled as 1H."""
        labels = []
        try:
            labels.extend(getattr(self, 'labb', []) or [])
        except Exception:
            pass
        try:
            labels.append(self._conversion_value('label1', key='label1', default=''))
        except Exception:
            pass
        for label in labels[:1] or labels:
            text = re.sub(r'[^a-z0-9]', '', str(label).lower())
            if text in ('h', 'h1', '1h', 'proton') or text.startswith('h1'):
                return True
        return False

    def auto_select_strongest_fid(self):
        """Select the preview trace with the largest absolute spectral signal.

        For a proton direct dimension the water region (4.6--5.0 ppm) is
        excluded from the score, so a dominant solvent peak cannot choose the
        phasing trace.  The selection is made from the processed preview stack
        produced by conversion and is then persisted through _set_fid_selection.
        """
        if self._fid_selection_count() <= 1:
            self._set_fid_selection(1)
            return 1

        source = None
        for infile in self._direct_spectrum_candidates(include_phased=False):
            if infile and os.path.exists(infile):
                source = infile
                break
        if source is None:
            return self._current_fid_selection()

        try:
            dic, data = ng.pipe.read(source)
            data = numpy.asarray(data)
            if data.ndim <= 1:
                return self._current_fid_selection()
            uc = ng.pipe.make_uc(dic, data)
            ppm = numpy.asarray(uc.ppm_scale())
            proton = self._direct_dimension_is_proton()
            mask = numpy.ones(ppm.shape, dtype=bool)
            if proton and ppm.size == data.shape[-1]:
                mask &= ~((ppm >= 4.6) & (ppm <= 5.0))
            if not numpy.any(mask):
                mask = numpy.ones(ppm.shape, dtype=bool)

            best_sel = 1
            best_score = -numpy.inf
            for sel in range(1, self._fid_selection_count() + 1):
                indices = list(self._fid_selection_indices(sel))
                lead_axes = data.ndim - 1
                if len(indices) < lead_axes:
                    indices.extend([0] * (lead_axes - len(indices)))
                indices = indices[:lead_axes]
                trace = numpy.asarray(data[tuple(indices) + (slice(None),)])
                values = numpy.abs(trace[mask]) if trace.shape[-1] == mask.size else numpy.abs(trace)
                finite = values[numpy.isfinite(values)]
                score = float(numpy.max(finite)) if finite.size else -numpy.inf
                if score > best_score:
                    best_score = score
                    best_sel = sel
            self._set_fid_selection(best_sel)
            self._apply_fid_dependent_indirect_phase(best_sel)
            return best_sel
        except Exception:
            return self._current_fid_selection()

    def _apply_fid_dependent_indirect_phase(self, selection=None):
        """Apply the FID quadrature convention to the first indirect phase.

        For 2D data, preview trace 2 requires the same +90 degree zero-order
        phase correction on the first indirect dimension that 3D preview
        traces 3 and 4 require.  Keep both the live Processing control (when
        that window exists) and shared GUI state in sync.  The normal Process
        save transaction remains responsible for writing the parameter file.
        """
        ndim = self._spectral_dimension_count()
        sel = self._normalize_fid_selection(selection)
        needs_phase = (ndim == 2 and sel == 2) or (ndim == 3 and sel in (3, 4))
        if not needs_phase:
            return False

        value = '90'
        frame = getattr(self, 'processing_frame', None)
        ctrl = getattr(frame, 'p0_1', None) if frame is not None else None
        if ctrl is not None:
            try:
                ctrl.SetValue(value)
            except Exception:
                pass

        state = getattr(self, 'state', None)
        if state is not None:
            try:
                state.update_gui_settings({'p0_1': value})
            except Exception:
                try:
                    state.gui_settings['p0_1'] = value
                except Exception:
                    pass
        return True

    def _debug_plot_state(self, context):
        try:
            line_count = len(getattr(self.axes, 'lines', []))
        except Exception:
            line_count = 'unknown'
        pass
        try:
            for i, line in enumerate(getattr(self.axes, 'lines', [])[:3]):
                try:
                    xdata = numpy.asarray(line.get_xdata())
                    ydata = numpy.asarray(line.get_ydata())
                    pass
                    pass
                    pass
                except Exception:
                    pass
        except Exception:
            pass

    def _direct_window_args(self):
        window = self._ctrl_value('windowBox0', default='GM') or 'GM'
        try:
            q1 = float(self._ctrl_value('win3Val0', default='2.0'))
        except Exception:
            q1 = 2.0
        try:
            q2 = float(self._ctrl_value('win2Val0', default='20.0'))
        except Exception:
            q2 = 20.0
        try:
            c = float(self._ctrl_value('firstPoint0', default='0.5'))
        except Exception:
            c = 0.5
        return window, q1, q2, c

    def _direct_apodization_settings(self):
        """Return the direct-dimension window settings used for slice.ft1 generation."""
        window_map = {0: 'GM', 1: 'SP', 2: 'EM'}
        window = str(self._ctrl_value('windowBox0', default='')).strip().upper()
        if window not in window_map.values():
            try:
                window = window_map.get(int(self._parse_param('window0', default=0)), 'GM')
            except Exception:
                window = 'GM'

        def _float_from_control(name, param_key, default):
            raw = self._ctrl_value(name, default='')
            try:
                return float(raw)
            except Exception:
                try:
                    return float(self._parse_param(param_key, default=default))
                except Exception:
                    return float(default)

        q1 = _float_from_control('win3Val0', 'win3Val0', 2.0)
        q2 = _float_from_control('win2Val0', 'win2Val0', 20.0)
        c = _float_from_control('firstPoint0', 'firstPoint0', 0.5)
        return window, q1, q2, c

    def _direct_phase_script_path(self):
        base = self._spec_output_dir()
        return os.path.join(base, 'nmrproc.1D.phase.com') if base else 'nmrproc.1D.phase.com'

    def _direct_phased_spectrum_path(self):
        base = self._spec_output_dir()
        return os.path.join(base, 'slice.phased.ft1') if base else 'slice.phased.ft1'

    def _write_direct_phase_script(self, p0, p1):
        return self.nmrPipe.write_direct_phase_script(self, p0, p1)

    def _run_direct_phase_script(self, p0, p1):
        return self.nmrPipe.run_direct_phase_script(self, p0, p1)

    def _direct_spectrum_candidates(self, include_phased=True):
        import glob
        base = self._spec_output_dir()
        if not base:
            return []

        candidates = []
        if include_phased:
            for name in ('slice.phased.ft1', 'slice.phased.ft1.gz'):
                candidates.append(os.path.join(base, name))
        for name in ('slice.ft1', 'slice.ft1.gz', 'test.ft', 'test.ft.gz', 'test.ft2', 'test.ft2.gz', 'test.ft3', 'test.ft3.gz', 'test.ft4', 'test.ft4.gz'):
            candidates.append(os.path.join(base, name))
        for pat in ('test*.ft', 'test*.ft.gz'):
            candidates.extend(sorted(glob.glob(os.path.join(base, pat))))
        seen = set()
        ordered = []
        for item in candidates:
            if item and item not in seen:
                seen.add(item)
                ordered.append(item)
        return ordered

    def _load_direct_spectrum_file(self, infile, label='spectrum'):
        dic, data = ng.pipe.read(infile)
        data = numpy.asarray(data)
        pass
        data = self._first_fid_slice(data)
        pass
        pass
        pass
        return dic, data

    def _legacy_build_direct_frequency_data(self):
        pipe_proc = self._pipe_proc_module()
        if pipe_proc is None:
            raise RuntimeError('nmrglue.pipe_proc is not available')
        if not hasattr(self, 'arr'):
            self._load_direct_fid()
        dic = dict(getattr(self, 'arrdic', {}) or {})
        data = numpy.array(self._first_fid_slice(self.arr), copy=True)

        if self._ctrl_checked('cb_baseSol', False):
            dic, data = pipe_proc.sol(dic, data)

        window, q1, q2, c = self._direct_window_args()
        if window:
            dic, data = pipe_proc.apod(dic, data, qName=window, q1=q1, q2=q2, c=c)

        dic, data = pipe_proc.zf(dic, data, zf=2)
        ft_flags = self._direct_ft_flags()
        dic, data = pipe_proc.ft(dic, data, **ft_flags)
        return dic, data

    def _build_direct_frequency_data(self, include_phased=True):
        # Varian 3p data have two spectral dimensions plus a real pseudo axis.
        # The direct phasing preview for this acquisition is deliberately made
        # by fid.test.slice.com -> slice.fid -> nmrproc.1D.com -> slice.ft1.
        # Keep that processed NMRPipe spectrum authoritative: loading the raw
        # slice FID first can silently put this case back onto the legacy
        # nmrglue processing path instead of displaying the spectrum that the
        # normal preview script just generated.
        strict_slice = (getattr(self, 'tp', None) == 'var' and
                        getattr(self, 'has_pseudo_axis', False) and getattr(self, 'spectral_dim_count', 1) == 2)
        if strict_slice:
            base = self._spec_output_dir()
            slice_paths = [os.path.join(base, 'slice.ft1'),
                           os.path.join(base, 'slice.ft1.gz')]
            for infile in slice_paths:
                if not os.path.exists(infile):
                    continue
                # Do not suppress an error here.  If nmrproc.1D.com produced
                # slice.ft1, a read failure is a real preview-path error and
                # should not be hidden by rebuilding a different spectrum from
                # slice.fid.
                dic, data = self._load_direct_spectrum_file(infile, label='spectrum')
                self.direct_ft_dic = dic
                self.direct_ft_data = numpy.array(data, copy=True)
                self.direct_ft_uc = ng.pipe.make_uc(dic, data)
                self.direct_ft_path = infile
                return dic, data

        # Preserve the established behaviour for all other acquisition types.
        if not hasattr(self, 'arr') or self.arr is None:
            self._load_direct_fid()
        for infile in self._direct_spectrum_candidates(include_phased=include_phased):
            if not infile or not os.path.exists(infile):
                continue
            try:
                dic, data = self._load_direct_spectrum_file(infile, label='spectrum')
                self.direct_ft_dic = dic
                self.direct_ft_data = numpy.array(data, copy=True)
                self.direct_ft_uc = ng.pipe.make_uc(dic, data)
                self.direct_ft_path = infile
                return dic, data
            except Exception:
                pass
        dic, data = self._legacy_build_direct_frequency_data()
        self.direct_ft_dic = dic
        self.direct_ft_data = numpy.array(data, copy=True)
        self.direct_ft_uc = ng.pipe.make_uc(dic, data)
        self.direct_ft_path = None
        return dic, data

    def _extract_direct_region(self, dic, data):
        xmin, xmax = self.nmrPipe.get_xmin_xmax(self)
        pass
        if xmin == '*' and xmax == '*':
            pass
            return dic, data, ng.pipe.make_uc(dic, data)

        pipe_proc = self._pipe_proc_module()
        if pipe_proc is None:
            pass
            return dic, data, ng.pipe.make_uc(dic, data)

        uc = ng.pipe.make_uc(dic, data)
        x1 = self._ppm_to_point(uc, xmin)
        xn = self._ppm_to_point(uc, xmax)
        pass
        if x1 != 'default' and xn != 'default' and x1 > xn:
            x1, xn = xn, x1
            pass

        try:
            dic, data = pipe_proc.ext(dic, data, x1=x1, xn=xn, sw=True)
            pass
        except Exception as exc:
            pass
            pass
        return dic, data, ng.pipe.make_uc(dic, data)

    def _get_direct_phase_values(self):
        p0 = self._direct_phase_p0
        p1 = self._direct_phase_p1
        try:
            p0 = float(self._processing_live_value('p0', p0))
        except Exception:
            pass
        try:
            p1 = float(self._processing_live_value('p1', p1))
        except Exception:
            pass
        return p0, p1

    def _saved_direct_phase_values(self):
        if hasattr(self, 'sld_0') and hasattr(self, 'sld_1'):
            try:
                return (
                    self._phase_slider_value(self.sld_0, self.sld_0_mode_btn),
                    self._phase_slider_value(self.sld_1, self.sld_1_mode_btn),
                )
            except Exception:
                pass
        try:
            p0, p1 = self._get_direct_phase_values()
            return float(p0), float(p1)
        except Exception:
            pass
        return self._direct_phase_p0, self._direct_phase_p1

    def _sync_direct_phase_from_controls(self, refresh_plot=True, sync_sliders=True):
        try:
            p0, p1 = self._get_direct_phase_values()
        except Exception:
            p0, p1 = self._direct_phase_p0, self._direct_phase_p1
        self._set_direct_phase_values(p0, p1, sync_sliders=sync_sliders)
        if refresh_plot:
            try:
                self.draw_figure()
            except Exception:
                pass

    def _ppm_to_point(self, uc, value):
        if value in ('', '*', None):
            return 'default'
        try:
            return int(round(uc.i(float(value), 'ppm')))
        except Exception:
            try:
                return int(round(uc.i(value, 'ppm')))
            except Exception:
                return 'default'

    def _format_phase_value(self, value):
        try:
            return f'{float(value):.2f}'
        except Exception:
            return str(value)

    def _set_direct_phase_values(self, p0=None, p1=None, sync_sliders=False, sync_bound_controls=True):
        previous_sync_flag = getattr(self, '_syncing_direct_phase_controls', False)
        self._syncing_direct_phase_controls = True
        try:
            if p0 is not None:
                self._direct_phase_p0 = round(float(p0), 2)
            if p1 is not None:
                self._direct_phase_p1 = round(float(p1), 2)

            if sync_bound_controls:
                self._set_processing_widget_value('p0', self._format_phase_value(self._direct_phase_p0))
                self._set_processing_widget_value('p1', self._format_phase_value(self._direct_phase_p1))

            if sync_sliders:
                try:
                    self._set_phase_slider_value(self.sld_0, self.sld_0_mode_btn, self._direct_phase_p0)
                except Exception:
                    pass
                try:
                    self._set_phase_slider_value(self.sld_1, self.sld_1_mode_btn, self._direct_phase_p1)
                except Exception:
                    pass
        finally:
            self._syncing_direct_phase_controls = previous_sync_flag

    def _update_phase_readouts(self, p0=None, p1=None, show_fid=False):
        try:
            if show_fid:
                if hasattr(self, 'SetStatusText'):
                    self.SetStatusText('Showing FID')
                return
        except Exception:
            pass

        if p0 is None:
            p0 = getattr(self, '_direct_phase_p0', 0.0)
        if p1 is None:
            p1 = getattr(self, '_direct_phase_p1', 0.0)

        try:
            if hasattr(self, 'phase_annotation') and self.phase_annotation is not None:
                self.phase_annotation.set_text(f'P0 {p0:.2f}\nP1 {p1:.2f}')
        except Exception:
            pass

        try:
            if hasattr(self, 'SetStatusText'):
                self.SetStatusText(f'P0={p0:.2f}  P1={p1:.2f}')
        except Exception:
            pass

    def _render_direct_frequency(self, p0, p1):
        backend = self._direct_phase_backend()
        dic, data = self._build_direct_frequency_data(include_phased=False)

        pass
        pass
        pass

        pass
        try:
            if backend == 'pipe':
                phased_path = self._run_direct_phase_script(p0, p1)
                pass
                dic, data = self._load_direct_spectrum_file(phased_path, label='phased spectrum')
                self.direct_phase_path = phased_path
                self.direct_ft_path = phased_path
            else:
                pipe_proc = self._pipe_proc_module()
                if pipe_proc is None:
                    raise RuntimeError('nmrglue.pipe_proc is not available')
                dic, data = pipe_proc.ps(dic, data, p0=p0, p1=p1)
                pass
                self.direct_phase_path = self.direct_ft_path
            self.direct_phase_dic = dic
            self.direct_phase_data = numpy.array(data, copy=True)
            self.direct_phase_uc = ng.pipe.make_uc(dic, data)
            pass
            pass
            pass
        except Exception as exc:
            pass
            raise

        try:
            dic, data, uc = self._extract_direct_region(dic, data)
            pass
        except Exception as exc:
            pass
            raise
        x = numpy.asarray(uc.ppm_scale())
        y = -numpy.real(numpy.asarray(data))
        pass
        pass
        pass
        pass

        self.direct_phase_dic = dic
        self.direct_phase_data = numpy.array(data, copy=True)
        self.direct_phase_uc = uc
        self.direct_phase_axis = x
        self.direct_phase_real = y
        self.index = x
        self.phase_mask = numpy.ones_like(x, dtype=bool)
        return x, y

    def Read1Dslice(self):

        show_fid = False
        try:
            show_fid = bool(self.cb_show_fid.IsChecked())
        except Exception:
            show_fid = False
        if self._load_direct_fid(show_fid=show_fid):
            return
        self.arrdic = {}
        self.arr = numpy.zeros(256, dtype=complex)
        self.index = numpy.arange(self.arr.shape[0])
        self.READ1D = 0

    def redraw_view(self):
        # This callback belongs to the custom toolbar Redraw button.  It is the
        # one GUI action which deliberately recalculates the vertical limits.
        self.draw_figure(reset_y=True)

    def _phase_bbox_signature(self):
        """Pixel geometry used by the phase-preview blit cache.

        Cache the full figure rather than the axes rectangle.  WXAgg can map
        an axes-sized blit rectangle incorrectly after layout/DPI changes,
        leaving a stale strip at an edge of the displayed spectrum.
        """
        fig = getattr(self, 'fig', None)
        if fig is None:
            return None
        try:
            return tuple(round(float(v), 3) for v in fig.bbox.bounds)
        except Exception:
            return None

    def _on_plot_resize(self, event=None):
        """Never reuse a blit background captured at an old canvas size."""
        self._phase_blit_background = None
        self._phase_blit_bbox = None
        try:
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_plot_draw(self, event=None):
        """Cache a clean, current-size figure background for phase blitting."""
        if not getattr(self, '_phase_blit_enabled', False):
            return
        axes = getattr(self, 'axes', None)
        line = getattr(self, 'phasing', None)
        if axes is None or line is None:
            return
        try:
            self._phase_blit_background = self.canvas.copy_from_bbox(self.fig.bbox)
            self._phase_blit_bbox = self._phase_bbox_signature()
            axes.draw_artist(line)
            annotation = getattr(self, 'phase_annotation', None)
            if annotation is not None:
                axes.draw_artist(annotation)
            self.canvas.blit(self.fig.bbox)
        except Exception:
            self._phase_blit_background = None
            self._phase_blit_bbox = None

    def _blit_phase_preview(self):
        """Fast path for slider-driven phase updates, with safe fallback."""
        axes = getattr(self, 'axes', None)
        line = getattr(self, 'phasing', None)
        background = getattr(self, '_phase_blit_background', None)
        bbox_matches = (getattr(self, '_phase_blit_bbox', None) == self._phase_bbox_signature())
        if getattr(self, '_phase_blit_enabled', False) and axes is not None and line is not None and background is not None and bbox_matches:
            try:
                self.canvas.restore_region(background)
                axes.draw_artist(line)
                annotation = getattr(self, 'phase_annotation', None)
                if annotation is not None:
                    axes.draw_artist(annotation)
                self.canvas.blit(self.fig.bbox)
                return
            except Exception:
                self._phase_blit_background = None
                self._phase_blit_bbox = None
        try:
            self.canvas.draw_idle()
        except Exception:
            pass

    def draw_figure(self, reset_y=False):
        if not hasattr(self, 'fig'):
            return

        # GUI-triggered redraws retain the current vertical view.  Only the
        # toolbar Redraw callback passes reset_y=True.
        old_ylim = None
        if not reset_y:
            try:
                if hasattr(self, 'axes') and self.axes is not None:
                    old_ylim = self.axes.get_ylim()
            except Exception:
                old_ylim = None

        self._phase_blit_background = None
        self._phase_blit_bbox = None
        self.fig.clear()
        try:
            self.Read1Dslice()
        except Exception:
            pass

        self.axes = self.fig.add_subplot(111)
        self.phase_annotation = None
        try:
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
        except Exception:
            pass

        show_fid = False
        try:
            show_fid = bool(self.cb_show_fid.IsChecked())
        except Exception:
            show_fid = False

        if show_fid:
            x = numpy.arange(len(self.arr))
            y_real = numpy.real(numpy.asarray(self.arr))
            y_imag = numpy.imag(numpy.asarray(self.arr))
            self.index = x
            self.phase_mask = numpy.ones_like(x, dtype=bool)
            self.axes.plot(x, y_real, lw=1, label='Real')
            self.phasing, = self.axes.plot(x, y_imag, lw=1, label='Imag')
            self._debug_plot_state('after raw fid plot')
            try:
                self.axes.set_xlabel('Time-domain points')
                self.axes.set_ylabel('FID amplitude')
                self.axes.legend()
            except Exception:
                pass
            y = numpy.concatenate([numpy.asarray(y_real).ravel(), numpy.asarray(y_imag).ravel()])
            if old_ylim is not None:
                try:
                    self.axes.set_ylim(old_ylim)
                except Exception:
                    pass
            elif y.size:
                try:
                    ymin = float(numpy.nanmin(y))
                    ymax = float(numpy.nanmax(y))
                    if ymin == ymax:
                        ymin -= 1.0
                        ymax += 1.0
                    self.axes.set_ylim(ymin, ymax)
                except Exception:
                    pass

        else:
            try:
                p0 = self._phase_slider_value(self.sld_0, self.sld_0_mode_btn)
                p1 = self._phase_slider_value(self.sld_1, self.sld_1_mode_btn)
                pass
                x, y = self._render_direct_frequency(p0, p1)
                pass
                pass
                pass
            except Exception as exc:
                pass
                pass
                x = numpy.arange(len(self.arr))
                y = numpy.real(numpy.asarray(self.arr))
                self.index = x
                self.phase_mask = numpy.ones_like(x, dtype=bool)

            self.phasing, = self.axes.plot(x, y, lw=1)
            try:
                self.phase_annotation = self.axes.text(0.98, 0.98, f'P0 {p0:.2f}\nP1 {p1:.2f}',
                               transform=self.axes.transAxes, ha='right', va='top',
                               fontsize=9, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75, edgecolor='none'))
            except Exception:
                self.phase_annotation = None
                pass
            self._update_phase_readouts(p0, p1, show_fid=False)
            self._debug_plot_state('after spectrum plot')
            try:
                if len(x) > 1:
                    self.axes.set_xlim(x[0], x[-1])
            except Exception:
                pass
            if old_ylim is not None:
                try:
                    self.axes.set_ylim(old_ylim)
                except Exception:
                    pass
            elif y.size:
                try:
                    ymin = float(numpy.nanmin(y))
                    ymax = float(numpy.nanmax(y))
                    if ymin == ymax:
                        ymin -= 1.0
                        ymax += 1.0
                    self.axes.set_ylim(ymin, ymax)
                except Exception:
                    pass

        # Animated artists are excluded from the normal canvas draw so the
        # draw_event callback can cache a clean background for slider blitting.
        try:
            self.phasing.set_animated(bool(self._phase_blit_enabled))
            if getattr(self, 'phase_annotation', None) is not None:
                self.phase_annotation.set_animated(bool(self._phase_blit_enabled))
        except Exception:
            pass

        try:
            self.fig.tight_layout()
        except Exception:
            pass
        try:
            self.canvas.draw()
        except Exception:
            self.canvas.draw_idle()
        try:
            if hasattr(self, 'SetStatusText'):
                if show_fid:
                    status_text = 'Showing FID'
                else:
                    status_text = f'P0={p0:.2f}  P1={p1:.2f}' if 'p0' in locals() and 'p1' in locals() else 'Spectrum'
                self.SetStatusText(status_text)
                self._hover_default_status = status_text
        except Exception:
            pass

    def create_main_panel(self):
        self._build_status_panel()
        self.button_box()
        self.direct_dimension_box()
        self.set_default_values()
        self._install_default_hover_help()

        self.panelSizer = wx.BoxSizer(wx.VERTICAL)
        self.panelSizer.AddSpacer(20)
        if self.statusSizer is not None:
            self.panelSizer.Add(self.statusSizer, flag=wx.GROW)
        self.panelSizer.AddSpacer(10)
        self.panelSizer.Add(self.controls_sizer, flag=wx.GROW)
        self.panelSizer.AddSpacer(10)
        self.panelSizer.AddSpacer(10)
        self.panelSizer.Add(self.dataSizer, flag=wx.GROW)
        self.panelSizer.AddSpacer(20)

        self.splitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.splitSizer.AddSpacer(20)
        self.splitSizer.Add(self.panelSizer, 0)
        self.splitSizer.AddSpacer(20)

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self._phase_blit_enabled = bool(getattr(self.canvas, 'supports_blit', False))
        self._phase_blit_background = None
        self._phase_blit_bbox = None
        self.canvas.mpl_connect('draw_event', self._on_plot_draw)
        self.canvas.mpl_connect('resize_event', self._on_plot_resize)
        self.canvas.mpl_connect('motion_notify_event', self._on_canvas_motion)
        self.canvas.mpl_connect('figure_leave_event', self._on_canvas_leave)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_hover_leave_frame)
        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view)
        self.vbox.Add(self.canvas, 1, wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)

        self.splitSizer.Add(self.vbox, 10, flag=wx.GROW)
        self.SetSizerAndFit(self.splitSizer)
        self.draw_figure()
        self.UpdateLampLights()

    def _normalize_folder_text(self, path: str) -> str:
        path = str(path or '').strip()
        if not path:
            return ''
        try:
            rel = os.path.relpath(path, os.getcwd())
            if not rel.startswith('..'):
                return './' + rel.lstrip('./')
        except Exception:
            pass
        return path

    def OnButtonSetSpecPath(self, event):
        default_path = self._spec_output_dir() or './spec'
        dlg = wx.DirDialog(
            self,
            message='Choose a folder for spectrum output',
            defaultPath=os.path.abspath(default_path),
            style=wx.DD_DEFAULT_STYLE,
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.state.spec_path = self._normalize_folder_text(path)
            if hasattr(self.parent, 'specPathBox'): self.parent.specPathBox.SetValue(self.state.spec_path)
        dlg.Destroy()

    def set_default_values(self):
        """Load only parameters represented by ProcessFrame widgets."""
        p0 = self._parse_param('p0', default='0')
        p1 = self._parse_param('p1', default='0')

        # Direct phasing and FID selection are Process-window state.  Do not
        # load p0_1/p1_1 or any other processing-window settings here.
        self._set_direct_phase_values(p0, p1, sync_sliders=True, sync_bound_controls=False)
        self._set_fid_selection(self._parse_param('FIDsel', default='1'))
        phase_slice_mode = str(self._parse_param('phaseSliceMode', default='First')).strip().lower()
        phase_slice_mode = 'Summed' if phase_slice_mode == 'summed' else 'First'
        if hasattr(self, 'phaseSliceMode'):
            self.phaseSliceMode.SetStringSelection(phase_slice_mode)
        self.phaseSliceModeValue = phase_slice_mode
        if getattr(self, 'state', None) is not None:
            self.state.metadata['phaseSliceMode'] = phase_slice_mode
        projection_type = str(self._parse_param('projectionType', default='skyline')).strip().lower()
        if projection_type not in ('sum', 'skyline'):
            projection_type = 'skyline'
        self.projections.SetValue(projection_type)

        self.xmaxBox.SetValue(str(self._parse_param('xmax', default='')))
        self.xminBox.SetValue(str(self._parse_param('xmin', default='')))

        defaultDir = str(self._parse_param('fiddir', default='')).strip()
        if not defaultDir:
            defaultDir = self.FindData() or './raw'
        if getattr(self, 'state', None) is not None:
            self.state.raw_path = defaultDir
            self.state.spec_path = str(self._parse_param('specPath', default='./spec'))
            # Keep the NMR-tab controls as the visible editors of project paths.
            if hasattr(self.parent, 'outPathBox'):
                self.parent.outPathBox.SetValue(self.state.raw_path)
            if hasattr(self.parent, 'specPathBox'):
                self.parent.specPathBox.SetValue(self.state.spec_path)
        self.UpdateLampLights()

    def OnButtonProcessing(self,event):
        self.sync_current_gui_state()
        frame = getattr(self, 'processing_frame', None)
        if frame is not None:
            try:
                if not frame.IsShown():
                    frame.Show(True)
                try:
                    frame.Layout()
                    frame.Refresh()
                except Exception:
                    pass
                frame.Raise()
                frame.SetFocus()
                return
            except Exception:
                self.processing_frame = None
        self.conv_frame = None
        from spinDecon.gui.dialogs.processing.settings import ProcessingFrame
        self.processing_frame = ProcessingFrame(self)
        self.processing_frame.Show(True)

    def OnButtonProjections(self, event):
        self.sync_current_gui_state()
        if getattr(self, 'projections_frame', None) is not None:
            try:
                self.projections_frame.Raise()
                self.projections_frame.SetFocus()
                return
            except Exception:
                self.projections_frame = None
        from spinDecon.gui.dialogs.processing.projections import ProjectionsFrame
        self.projections_frame = ProjectionsFrame(self)
        self.projections_frame.Show(True)

    def _normalize_nus_schedule(self, value, base_dir=''):
        value = str(value or '').strip()
        if not value:
            return ''
        if not os.path.isabs(value):
            normalized = os.path.normpath(value)
            return '' if normalized in ('.', '') else normalized
        if not base_dir:
            return value
        try:
            return os.path.relpath(os.path.abspath(value), os.path.abspath(base_dir))
        except Exception:
            return value

    def _current_nus_schedule(self):
        base_dir = self._raw_output_dir()
        conv = getattr(self, 'conv_frame', None)
        if conv is not None and hasattr(conv, 'nusFil'):
            try:
                value = self._normalize_nus_schedule(conv.nusFil.GetValue().strip(), base_dir)
                if value:
                    return value
            except Exception:
                pass
        try:
            value = str(self.parent.Parse(self.parent.deconParFile, 'nusFil', default='')).strip()
        except Exception:
            value = ''
        value = self._normalize_nus_schedule(value, base_dir)
        if value:
            return value

        # The conversion stage saves the NUS sampling schedule alongside the
        # converted data under the canonical filename `schedule`.  Processing
        # may be opened without the conversion window, so discover that file
        # directly instead of requiring nusFil to have been persisted first.
        if base_dir:
            for filename in ('schedule', 'nuslist'):
                candidate = os.path.join(base_dir, filename)
                if os.path.isfile(candidate):
                    return self._normalize_nus_schedule(candidate, base_dir)
        return ''

    def _parameter_file_path(self):
        """Resolve the actual parameter file used by ProcessFrame.

        ``deconParFile`` may be a basename or an already-resolved absolute
        path.  Always prefer an existing concrete file so load/save operations
        inspect the current project file rather than a default/template file.
        """
        raw = str(getattr(self.parent, 'deconParFile', '') or '').strip()
        if os.path.isabs(raw):
            return raw
        if raw:
            directory = ''
            try:
                directory = str(self.parent.dirBox.GetValue()).strip()
            except Exception:
                pass
            candidate = os.path.join(directory, raw) if directory else raw
            if os.path.exists(candidate) or not os.path.exists(raw):
                return candidate
        return raw

    def _parse_param(self, key, default=''):
        try:
            return self.parent.Parse(self._parameter_file_path(), key, default=default)
        except Exception:
            return default

    def _parse_allstr(self, key, default=''):
        try:
            value = self.parent.ParseAllStr(self._parameter_file_path(), key)
        except Exception:
            return default
        if value in (0, '', None):
            return default
        return value

    def _parse_bool(self, key, default=False):
        raw = str(self._parse_param(key, default='')).strip().lower()
        if raw in ('y', 'yes', 'true', '1', 't'):
            return True
        if raw in ('n', 'no', 'false', '0', 'f'):
            return False
        return default

    def _conversion_frame(self):
        return getattr(self, 'conv_frame', None)

    def _manual_reference_ppm(self):
        conv = self._conversion_frame()
        if conv is not None and hasattr(conv, 'xcenBox'):
            try:
                value = str(conv.xcenBox.GetValue()).strip()
                if value and value not in ('0', 'None'):
                    return value
            except Exception:
                pass
        # Closed dialogs still leave their newest value in shared live state.
        state = getattr(self, 'state', None)
        if state is not None and 'xcen' in state.gui_settings:
            value = str(state.gui_settings['xcen']).strip()
            if value and value not in ('0', 'None'):
                return value
        try:
            value = str(self._parse_param('xcen', default='')).strip()
            if value and value not in ('0', 'None'):
                return value
        except Exception:
            pass
        return ''

    def _conversion_value(self, name, key=None, default=''):
        conv = self._conversion_frame()
        if conv is not None and hasattr(conv, name):
            try:
                value = getattr(conv, name).GetValue()
                if value is not None:
                    return str(value)
            except Exception:
                pass
        key = key or name
        state = getattr(self, 'state', None)
        if state is not None and key in state.gui_settings:
            value = str(state.gui_settings[key]).strip()
            if key.startswith('label'):
                value = value.replace(' ', '')
        elif key.startswith('label'):
            value = str(self._parse_allstr(key, default=default)).replace(' ', '').strip()
        else:
            value = str(self._parse_param(key, default=default)).strip()
        if value in ('0', 'None'):
            return default
        return value

    def _conversion_checked(self, name, key=None, default=False):
        conv = self._conversion_frame()
        if conv is not None and hasattr(conv, name):
            try:
                return bool(getattr(conv, name).IsChecked())
            except Exception:
                pass
        key = key or name
        state = getattr(self, 'state', None)
        if state is not None and key in state.gui_settings:
            return str(state.gui_settings[key]).strip().lower() in ('y', 'yes', 'true', '1', 't')
        return self._parse_bool(key, default=default)

    def _sync_conversion_dialog(self):
        conv = getattr(self, 'conv_frame', None)
        if conv is not None:
            try:
                conv._copy_to_parent()
            except Exception:
                pass

    def _bind_processing_controls(self, frame):
        """Register the Processing window without aliasing its wx controls.

        Historically ProcessFrame copied every ProcessingFrame widget onto
        itself.  That made two frames appear to own the same controls and made
        script generation depend on those aliases.  ProcessingFrame is now the
        widget owner; ProcessFrame consumers resolve through the frame/shared
        state explicitly.
        """
        self.processing_frame = frame
        try:
            self.ncpus = self._current_ncpus_value(default=self.ncpus)
        except Exception:
            pass
        try:
            p0, p1 = self._get_direct_phase_values()
            self._set_direct_phase_values(p0, p1, sync_sliders=True, sync_bound_controls=False)
        except Exception:
            pass

    def _unbind_processing_controls(self):
        """Forget the Processing window; no ProcessFrame widget aliases exist."""
        self.processing_frame = None

    def _processing_widget(self, name):
        frame = getattr(self, 'processing_frame', None)
        return getattr(frame, name, None) if frame is not None else None

    def _processing_state_key(self, name):
        mapping = {
            'maxIterBox': 'maxIterSMILE', 'ncpusBox': 'ncpus',
            'cb_baseLin': 'lin', 'cb_basepol': 'poly', 'cb_baseSol': 'sol',
            'cb_ft0': 'flip0', 'windowBox0': 'window0',
        }
        if name in mapping:
            return mapping[name]
        for idx in range(1, 4):
            per_dim = {
                f'cb_f{idx}180': f'f{idx}180', f'cb_lp{idx}': f'lp{idx}',
                f'cb_basepol{idx}': f'bl{idx}', f'cb_ft{idx}': f'flip{idx}',
                f'windowBox{idx}': f'window{idx}',
            }
            if name in per_dim:
                return per_dim[name]
        return name

    def _processing_live_value(self, name, default=''):
        ctrl = self._processing_widget(name)
        if ctrl is not None:
            try:
                if hasattr(ctrl, 'IsChecked'):
                    return bool(ctrl.IsChecked())
                return ctrl.GetValue()
            except Exception:
                pass
        state = getattr(self, 'state', None)
        key = self._processing_state_key(name)
        if state is not None and key in state.gui_settings:
            return state.gui_settings[key]
        return default

    def _set_processing_widget_value(self, name, value):
        """Update the owning Processing widget and live state, if available."""
        ctrl = self._processing_widget(name)
        if ctrl is not None:
            try:
                if name.startswith('windowBox') and hasattr(ctrl, 'SetSelection') and isinstance(value, int):
                    ctrl.SetSelection(value)
                else:
                    ctrl.SetValue(value)
            except Exception:
                pass
        state = getattr(self, 'state', None)
        if state is not None:
            state.update_gui_settings({self._processing_state_key(name): value})

    def _discover_dimension_labels(self):
        """Discover raw labels from dataset metadata without opening a child GUI."""
        if getattr(self, 'tp', '') == 'bruk':
            return discover_bruker_labels(self._raw_output_dir(), self.spectral_dim_count)
        return []

    def _hydrate_shared_dimension_labels(self):
        """Initialise the authoritative Process-session raw dimension labels.

        Saved project labels win.  Missing labels are filled from vendor
        metadata (currently Bruker NUC1).  Hydration is deliberately non-dirty.
        """
        state = getattr(self, 'state', None)
        if state is None:
            return {}
        discovered = self._discover_dimension_labels()
        updates = {}
        for idx in range(1, min(int(self.spectral_dim_count), 4) + 1):
            value = clean_dimension_label(self._parse_allstr(f'label{idx}', default=''))
            if not value and idx <= len(discovered):
                value = clean_dimension_label(discovered[idx - 1])
            if value:
                updates[f'label{idx}'] = value
        if updates:
            hydrate = getattr(state, 'hydrate_gui_settings', None)
            if callable(hydrate):
                hydrate(updates, overwrite=True)
            else:
                state.gui_settings.update(updates)
        return updates

    def get_dimension_labels(self):
        """Return authoritative raw spectral labels for this Process session."""
        state = getattr(self, 'state', None)
        live = getattr(state, 'gui_settings', {}) if state is not None else {}
        labels = []
        for idx in range(1, int(self.spectral_dim_count) + 1):
            value = clean_dimension_label(live.get(f'label{idx}', ''))
            if not value:
                value = clean_dimension_label(self._parse_allstr(f'label{idx}', default=''))
            labels.append(value or 'H1')
        return labels

    def set_dimension_labels(self, labels, refresh=True):
        """Publish raw labels to the shared Process store and refresh dependants."""
        updates = {}
        for idx, raw in enumerate(list(labels)[:min(int(self.spectral_dim_count), 4)], 1):
            value = clean_dimension_label(raw)
            if value:
                updates[f'label{idx}'] = value
        state = getattr(self, 'state', None)
        if state is not None and updates:
            state.update_gui_settings(updates)
        self.GetLabs()
        if refresh:
            self.SetLab(refresh=False)
        return updates

    def get_spectral_labels(self):
        """Return canonical (duplicate-disambiguated) Process-session labels."""
        labels = canonical_spectral_labels(self.get_dimension_labels())
        self.labb = labels
        return list(labels)

    def SetLab(self, refresh=True):
        if not hasattr(self, 'lab0') or not hasattr(self, 'labb'):
            return
        pass
        if(type(self.spectral_dim_count)!=str):
            if(self.spectral_dim_count>=2):
                pass

                self.lab0.SetLabel(self.labb[0])
                pass

                self.lab1.SetLabel(self.labb[1])
            if(self.spectral_dim_count>=3):
                pass
                self.lab2.SetLabel(self.labb[2])
            if(self.spectral_dim_count>=4):
                self.lab3.SetLabel(self.labb[3])
        elif(self.has_pseudo_axis and self.spectral_dim_count == 1):
            self.lab0.SetLabel(self.labb[0])
        elif(self.has_pseudo_axis and self.spectral_dim_count == 2):
            self.lab0.SetLabel(self.labb[0])
            self.lab1.SetLabel(self.labb[1])

        if refresh:
            self.draw_figure()

    def OnButtonClean(self,event):
        base = self._spec_output_dir()
        for rel in ('fids', 'XYZA', 'fid_full', 'pdata', 'test.fid', 'test.fid.gz', 'ser_full', 'ser_full.gz', 'slice.fid', 'slice.ft1', 'slice.phased.ft1', 'projections', 'projections1D', 'fid.test.com', 'fid.test.slice.com', 'nmrproc.1D.com', 'nmrproc.1D.phase.com', 'nmrproc.test.com', 'nmrprocLP.com'):
            path = os.path.join(base, rel)
            if os.path.isdir(path):
                os.system('rm -rf '+path)
            elif os.path.exists(path):
                os.system('rm -rf '+path)
    def OnButtonClose(self,event):
        self.Close()

    def OnClose(self, event):
        changed = self.parameters_changed_since_baseline()
        if event.CanVeto() and changed:
            answer = wx.MessageBox(
                "Parameters have been changed, would you like to save?",
                "Before exiting",
                wx.ICON_QUESTION | wx.YES_NO,
            )
            if answer == wx.YES:
                self.save_current_gui_state(reason='process-close')
            else:
                self.discard_parameter_changes()
        try:
            if getattr(self.parent, 'process_frame', None) is self:
                self.parent.process_frame = None
        except Exception:
            pass
        self.Destroy()

    def IsInParse(self,val):
        inny=open(self.parent.deconParFile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(test[0]==val):
                    return 1
        return 0

    def IntToBool(self,booly):
        if(booly==True):
            return 'y'
        else:
            return 'n'

    def _ctrl_value(self, name, default=''):
        ctrl = getattr(self, name, None)
        if ctrl is not None:
            try:
                return ctrl.GetValue()
            except Exception:
                pass
        return self._processing_live_value(name, default)

    def _ctrl_checked(self, name, default=False):
        ctrl = getattr(self, name, None) or self._processing_widget(name)
        if ctrl is not None:
            try:
                return bool(ctrl.IsChecked())
            except Exception:
                pass
        value = self._processing_live_value(name, default)
        if isinstance(value, str):
            return value.strip().lower() in ('y', 'yes', 'true', '1', 't')
        return bool(value)

    def _ctrl_selection(self, name, default=0):
        ctrl = getattr(self, name, None) or self._processing_widget(name)
        if ctrl is not None:
            try:
                return int(ctrl.GetSelection())
            except Exception:
                pass
        state = getattr(self, 'state', None)
        key = self._processing_state_key(name)
        try:
            return int(state.gui_settings[key]) if state is not None and key in state.gui_settings else int(default)
        except Exception:
            return int(default)

    def _write_clean_parameter_line(self, path, key, value):
        lines = []
        target = str(path)
        existing = []
        if os.path.exists(target):
            with open(target, 'r') as fh:
                existing = fh.readlines()
        replaced = False
        for line in existing:
            stripped = line.strip()
            tokens = stripped.split()
            if len(tokens) >= 2 and tokens[0] == key and tokens[1] == '=':
                if not replaced:
                    lines.append(f'{key} = {value}\n')
                    replaced = True
                continue
            lines.append(line)
        if not replaced:
            lines.append(f'{key} = {value}\n')
        with open(target, 'w') as fh:
            fh.writelines(lines)

    def _current_ncpus_value(self, default=1):
        frame = getattr(self, 'processing_frame', None)
        ctrl = getattr(frame, 'ncpusBox', None) if frame is not None else None
        if ctrl is not None:
            try:
                return max(1, int(float(ctrl.GetValue())))
            except Exception:
                pass
        try:
            return max(1, int(float(getattr(self, 'ncpus', default))))
        except Exception:
            pass
        try:
            core_box = getattr(self.parent, 'coreBox', None)
            if core_box is not None:
                return max(1, int(float(core_box.GetValue())))
        except Exception:
            pass
        return max(1, int(default))

    def _current_nmrpipe_input_file(self):
        """Return the authoritative live NMRPipe spectrum selection."""
        parent = getattr(self, 'parent', None)
        box = getattr(parent, 'infileBox', None) if parent is not None else None
        if box is not None:
            try:
                value = str(box.GetValue() or '').strip()
                if value:
                    return value
            except Exception:
                pass
        state = getattr(self, 'state', None)
        return str(getattr(state, 'input_file', '') or '').strip()

    def collect_updates(self):
        """Collect Process-window values without touching the parameter file."""
        self._sync_conversion_dialog()
        p0_value, p1_value = self._saved_direct_phase_values()
        updates = {
            'FIDsel': str(self._current_fid_selection()),
            'phaseSliceMode': (self.phaseSliceMode.GetStringSelection() if hasattr(self, 'phaseSliceMode') else 'First') or 'First',
            'p0': self._format_phase_value(p0_value),
            'p1': self._format_phase_value(p1_value),
            'xmin': self._ctrl_value('xminBox', ''),
            'xmax': self._ctrl_value('xmaxBox', ''),
            'fiddir': getattr(self.state, 'raw_path', './raw'),
            'specPath': getattr(self.state, 'spec_path', './spec'),
            # The main NMR tab owns the selected NMRPipe spectrum name.  Keep
            # it in the Process-family save transaction so spectra produced by
            # Processing/Auto Processing survive the next project session.
            'infile': self._current_nmrpipe_input_file(),
            'projectionType': self.projections.GetValue() or 'skyline',
        }
        if getattr(self, 'state', None) is not None:
            self.state.metadata['FIDsel'] = self._current_fid_selection()
            self.state.metadata['phaseSliceMode'] = updates['phaseSliceMode']
            self.state.metadata['projectionType'] = updates['projectionType']
            self.state.update_gui_settings(updates)
        return updates

    def _collect_current_parameter_updates(self, *, publish=False):
        """Collect persistable Process-family values, optionally publishing live state.

        The non-publishing form is deliberately side-effect free with respect to
        ProjectState and is used by close-time dirty detection.
        """
        state = getattr(self, 'state', None)
        saved_gui = copy.deepcopy(getattr(state, 'gui_settings', {})) if state is not None else None
        saved_metadata = copy.deepcopy(getattr(state, 'metadata', {})) if state is not None else None
        saved_dirty = getattr(state, 'dirty', False) if state is not None else False
        try:
            updates = dict(self.collect_updates())
            for name in ('processing_frame', 'conv_frame', 'projections_frame'):
                frame = getattr(self, name, None)
                collector = getattr(frame, 'collect_updates', None) if frame is not None else None
                if callable(collector):
                    try:
                        updates.update(collector(update_state=publish))
                    except TypeError:
                        updates.update(collector())
            if publish and state is not None:
                state.update_gui_settings(updates)
            return updates
        finally:
            if not publish and state is not None:
                state.gui_settings.clear()
                state.gui_settings.update(saved_gui or {})
                state.metadata.clear()
                state.metadata.update(saved_metadata or {})
                state.dirty = saved_dirty

    @staticmethod
    def _normalise_parameter_value(value):
        return str(value).strip()

    def _read_parameter_file_snapshot(self):
        values = {}
        path = self._parameter_file_path()
        if os.path.exists(path):
            with open(path, 'r') as fh:
                for line in fh:
                    tokens = line.split()
                    if len(tokens) >= 3 and tokens[1] == '=':
                        values[tokens[0]] = tokens[2].strip()
        return values

    def _capture_parameter_baseline(self, updates=None):
        """Record the last accepted Process-family state without writing disk."""
        baseline = self._read_parameter_file_snapshot()
        state = getattr(self, 'state', None)
        if state is not None:
            for key, value in getattr(state, 'gui_settings', {}).items():
                baseline[key] = self._normalise_parameter_value(value)
        if updates:
            for key, value in updates.items():
                baseline[key] = self._normalise_parameter_value(value)
        self._parameter_baseline = baseline
        self._gui_settings_baseline = copy.deepcopy(getattr(state, 'gui_settings', {})) if state is not None else {}
        return baseline

    def parameters_changed_since_baseline(self):
        state = getattr(self, 'state', None)
        current = dict(getattr(state, 'gui_settings', {})) if state is not None else {}
        current.update(self._collect_current_parameter_updates(publish=False))
        baseline = getattr(self, '_parameter_baseline', {})
        for key, value in current.items():
            if self._normalise_parameter_value(value) != self._normalise_parameter_value(baseline.get(key, '')):
                return True
        return False

    def discard_parameter_changes(self):
        """Discard Process-family live changes without touching the parameter file."""
        state = getattr(self, 'state', None)
        if state is None:
            return
        state.gui_settings.clear()
        state.gui_settings.update(copy.deepcopy(getattr(self, '_gui_settings_baseline', {})))
        state.projection_phase_preview = {}
        state.dirty = False

    def sync_current_gui_state(self):
        """Snapshot all currently open GUI windows into shared live state."""
        return self._collect_current_parameter_updates(publish=True)

    def save_current_gui_state(self, *, reason='user'):
        """Atomically commit the current widgets from all open windows."""
        updates = self.sync_current_gui_state()
        savefile = self._parameter_file_path()
        source_path = savefile if os.path.exists(savefile) else None
        update_parameter_file(savefile, updates, source_path=source_path)
        if getattr(self, 'state', None) is not None:
            self.state.dirty = False
            self.state.metadata['last_gui_save_reason'] = str(reason)
        self._capture_parameter_baseline(updates)
        return updates


    #def GetDataPath(self):
    #def GetDataPath(self):
    #    return os.path.join(self.DataStoreBox.GetValue(),self.FidPathBox.GetValue()   )

    def OnCombineBruker(self, event):
        """Combine numbered Bruker child experiments directly into the fid path."""
        from spinDecon.gui.dialogs.processing.bruker_combiner import CombineBrukerFrame
        raw_dir = str(self._raw_output_dir()).strip()
        if not raw_dir:
            errorMessage('Please set the fid path to the folder containing the numbered Bruker experiments.')
            return
        try:
            os.makedirs(raw_dir, exist_ok=True)
        except Exception as exc:
            errorMessage('Cannot use fid path %s: %s' % (raw_dir, exc))
            return
        dlg = CombineBrukerFrame(self, raw_dir)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_OK:
            # The combiner populates the existing fid path.  Re-run format
            # detection now that acqus/acqu2s exist, then refresh status.
            self.tp = None
            self.GetSpectrometerType()
            self.UpdateLampLights()
            try:
                self.draw_figure()
            except Exception:
                pass

    def GetSpectrometerType(self):
        """Detect the spectrometer type for the current raw-data directory.

        The project state is the source of truth for paths; do not cache a
        duplicate storePath/storeExists state on the ProcessFrame.
        """
        raw_dir = str(self._raw_output_dir()).strip()
        if not os.path.exists(raw_dir):
            pass
            return

        if os.path.exists(os.path.join(raw_dir, 'acqus')) or os.path.exists(os.path.join(raw_dir, 'acqu2s')):
            self.tp = 'bruk'
        elif os.path.exists(os.path.join(raw_dir, 'procpar')):
            self.tp = 'var'
            self.parfile = os.path.join(raw_dir, 'procpar')
        else:
            self.tp = 'omega'

    def DoProjections(self, pipefile, output_frame=None):
        """Generate downstream projections after running a processing script."""
        spec_base = self._spec_output_dir()
        pipe_path = os.path.join(spec_base, pipefile)
        projection_type = 'skyline'
        try:
            projection_type = self.projections.GetValue() or 'skyline'
        except Exception:
            pass
        # Pseudo topology must be dispatched before the ordinary spectral-count
        # cases: 2 spectral + 1 pseudo has spectral_dim_count == 2 but its
        # processed array is physically 3D and must first collapse the pseudo
        # axis with MakeProj3P.
        if self.has_pseudo_axis and self.spectral_dim_count == 2:
            if output_frame is not None:
                wx.CallAfter(output_frame.append_text, 'Creating the 2D spectral projection and two 1D spectral projections from the processed pseudo-3D spectrum...\n')
            progress = (lambda text: wx.CallAfter(output_frame.append_text, text)) if output_frame is not None else None
            MakeProj3P(pipe_path, folder=os.path.join(spec_base, 'projections'), progress=progress)
        elif self.spectral_dim_count == 3:
            if output_frame is not None:
                wx.CallAfter(output_frame.append_text, 'Creating 2D and 1D projections from the processed 3D spectrum...\n')
            progress = (lambda text: wx.CallAfter(output_frame.append_text, text)) if output_frame is not None else None
            MakeProj3D(pipe_path, folder=os.path.join(spec_base, 'projections'), OneD=True, projection_type=projection_type, progress=progress)
        elif self.spectral_dim_count == 2:
            if output_frame is not None:
                wx.CallAfter(output_frame.append_text, 'Creating 1D projections from the processed 2D spectrum...\n')
            progress = (lambda text: wx.CallAfter(output_frame.append_text, text)) if output_frame is not None else None
            MakeProj2D(pipe_path, projection_type=projection_type, progress=progress)
        elif self.spectral_dim_count == 4:
            pass

    def RefreshDirectSlice(self, output_frame=None, on_finish=None):
        """Regenerate the direct-dimension preview, streaming NMRPipe output when requested."""
        from spinDecon.gui.dialogs.shell_output import run_command_with_output

        spec_base = self._spec_output_dir()
        script_path = os.path.join(spec_base, 'nmrproc.1D.com') if spec_base else 'nmrproc.1D.com'
        project_root = os.path.abspath(os.path.join(spec_base, os.pardir)) if spec_base else os.getcwd()
        self.nmrPipe.make_proc_script_1d_slice(self, script_path)

        def done(rc=0):
            if rc in (0, None):
                try:
                    self.draw_figure()
                except Exception:
                    pass
                try:
                    self.UpdateLampLights()
                except Exception:
                    pass
            if on_finish is not None:
                try:
                    on_finish()
                except TypeError:
                    on_finish(rc)

        if output_frame is not None:
            output_frame.append_text('\n=== Refreshing direct-dimension preview ===\n')
            return run_command_with_output(
                ['csh', script_path], parent=self, title='Processing Output',
                cwd=project_root or None, output_frame=output_frame, on_finish=done,
                final=False, label='Process 1D preview slice')

        result = subprocess.run(['csh', script_path], cwd=project_root or None, capture_output=True, text=True, check=True)
        done(result.returncode)
        return script_path

    def _maybe_invert_test_ft_after_processing(self, script_path='', pipefile=''):
        """Invert the standard 1D processing test spectrum in place.

        This is only applied to 1D processing outputs named test.ft. It is not
        used when regenerating the direct phasing slice.
        """
        try:
            dim_count = self._spectral_dimension_count()
        except Exception:
            dim_count = 1
        if dim_count != 1 or pipefile != 'test.ft':
            return False

        spec_base = self._spec_output_dir()
        test_path = os.path.join(spec_base, pipefile)
        if not os.path.exists(test_path):
            return False

        try:
            dic, data = ng.pipe.read(test_path)
            data = data * -1
            ng.pipe.write(test_path, dic, data, overwrite=True)
            pass
            return True
        except Exception:
            pass
            raise

    def ResetReads(self):
        """Clear downstream read/peak state before re-running processing."""
        if hasattr(self.parent, 'READ'):
            self.parent.READ = 0
        if hasattr(self.parent, 'PEAK'):
            self.parent.PEAK = 0
        if hasattr(self.parent, 'parent') and hasattr(self.parent.parent, 'KillPage'):
            for page in ('2Dplanes', '1Ddeconv', '2Dslices'):
                try:
                    self.parent.parent.KillPage(page)
                except Exception:
                    pass

    def _processing_script_state(self):
        """Capture current processing values without exposing wx widgets.

        Open widgets have precedence over shared state.  Missing controls are
        supplied from ProjectState.gui_settings, so script generation has one
        explicit, immutable value snapshot.
        """
        frame = getattr(self, 'processing_frame', None)
        state = getattr(self, 'state', None)
        return ProcessingScriptState.capture_current(frame, state)

    def RenderProcessScript(self, lp='n'):
        return self.nmrPipe.RenderProcessScriptState(self, self._processing_script_state(), lp=lp)

    def WriteProcessScript(self, lp='n', outfile=None):
        return self.nmrPipe.WriteProcessScriptState(self, self._processing_script_state(), lp=lp, outfile=outfile)

    def BuildProcessScript(self, lp='n', outfile=None):
        return self.WriteProcessScript(lp=lp, outfile=outfile)

    def ExecuteProcessScript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        return self.nmrPipe.ExecuteProcessScript(self, script_path=script_path, lp=lp, on_finish=on_finish, title=title)

    def RunProcessScript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        return self.nmrPipe.RunProcessScript(self, script_path=script_path, lp=lp, on_finish=on_finish, title=title)

    def renderprocessscript(self, lp='n'):
        return self.RenderProcessScript(lp=lp)

    def writeprocessscript(self, lp='n', outfile=None):
        return self.WriteProcessScript(lp=lp, outfile=outfile)

    def executeprocessscript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        return self.ExecuteProcessScript(script_path, lp=lp, on_finish=on_finish, title=title)

    def runprocessscript(self, script_path, lp='n', on_finish=None, title='Processing Output'):
        return self.RunProcessScript(script_path, lp=lp, on_finish=on_finish, title=title)

    def GetLabs(self):
        labels = canonical_spectral_labels(self.get_dimension_labels())
        if self.has_pseudo_axis:
            real_name = self._conversion_value('RealName', key='RealName', default='')
            if real_name:
                labels.append(real_name)
        self.labb = labels

    def MakeConvScript(self, on_finish=None):
        # Commit current conversion/process widgets before building the script.
        self.save_current_gui_state(reason='conversion-auto')
        self._sync_conversion_dialog()
        self.READ1D = 0
        outfile = self._spec_output_dir() + '/fid.com'
        self.GetLabs()
        self.SetLab(refresh=False)

        if type(self.spectral_dim_count) != str:
            if len(self.labb) != self.spectral_dim_count:
                pass
                pass
                pass
                pass
                pass
                return -1
        else:
            if str(len(self.labb)) + 'p' != self.spectral_dim_count:
                pass
                pass
                pass
                pass
                pass
                return -1

        # Build the vendor-conversion object through the *same* code path as
        # Conversion Script -> Guess.  Keeping a second reconstruction of
        # labels/rk/topology here caused Auto and Guess to diverge for 3p
        # Bruker data (and, after fixing dim='3p', left mismatched list lengths).
        conv = getattr(self, 'conv_frame', None)
        temporary_conv = False
        if conv is None:
            from spinDecon.gui.dialogs.processing.conversion import ConversionFrame
            conv = ConversionFrame(self)
            temporary_conv = True
            try:
                conv.Hide()
            except Exception:
                pass
        try:
            inst = conv._build_vpar()
            self.vpar = inst
            script_path = inst.BuildConversionScript()
        finally:
            if temporary_conv:
                try:
                    conv.Destroy()
                except Exception:
                    pass
        if script_path == -1 or getattr(inst, 'abort', 0) == 1:
            if on_finish is not None:
                on_finish()
            return -1

        needs_extract = False
        try:
            needs_extract = inst._slice_dim_count() > 1 and not (self.has_pseudo_axis and self.spectral_dim_count == 2 and inst.tp == 'bruk')
        except Exception:
            pass
        steps = ['Convert raw data to NMRPipe']
        if needs_extract:
            steps.append('Prepare preview slice')
        steps.extend(['Process preview spectrum', 'Refresh display'])

        def workflow_done(*_args, **_kwargs):
            try:
                output_frame.start_step('Refresh display')
                self.UpdateLampLights()
                # Conversion has produced a new multidimensional preview stack.
                # Pick the strongest quadrature trace before drawing so both
                # the FID selector and any FID-dependent indirect phase are
                # correct on the very first display.
                self.auto_select_strongest_fid()
                # A newly loaded spectrum must establish its own vertical
                # limits instead of inheriting the empty/previous plot limits.
                self.draw_figure(reset_y=True)
                output_frame.append_text('\nConversion workflow complete.\n')
                output_frame.finish_workflow(True)
                output_frame.set_status('Complete')
            except Exception:
                pass
            if on_finish is not None:
                on_finish()

        def conversion_done(rc=0):
            if rc not in (0, None):
                try:
                    output_frame.finish_workflow(False)
                    output_frame.set_status('Conversion failed')
                except Exception:
                    pass
                if on_finish is not None:
                    on_finish()
                return
            try:
                inst.ProcessSlice(output_frame=output_frame, on_finish=workflow_done)
            except Exception as exc:
                output_frame.append_text('\nCould not start preview processing: %s\n' % exc)
                output_frame.finish_workflow(False)
                output_frame.set_status('Conversion failed')
                if on_finish is not None:
                    on_finish()

        output_frame = run_command_with_output(
            ['csh', script_path], parent=self, title='Conversion Output',
            on_finish=conversion_done, final=False,
            label='Convert raw data to NMRPipe', workflow_steps=steps,
            workflow_step=0)
        return output_frame



        # print(key, dic[key])
    # exit()




