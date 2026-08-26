#!/usr/bin/python

###################################################################
# Deconvolve nmr spectrum
###################################################################
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
import numpy,sys,os,wx,platform
import shutil
import nmrglue as ng
from spinDecon.gui.dialogs.processing.process import path_escape
from spinDecon.processing.nmrpipe_scripts import MakeProj4D, MakeProj3D, MakeProj2D, nmrglue_project2D_1D
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from spinDecon.gui.dialogs.errors import errorMessage
from spinDecon.gui.dialogs.shell_output import run_command_with_output
from spinDecon.domain.peaks import peakEntry, diagEntry
from spinDecon.domain.dimensions.viewer_contract import topology_for
import subprocess
from spinDecon.processing.vpar_decon import GetParBrukFile, vpar
from spinDecon.gui.dialogs.processing.nmrpipe_adapter import nmrPipe
import os
import re
# import matplotlib.pyplot as plt
# plt.rcParams['figure.dpi'] = 300


import copy
#import imp
import importlib

from spinDecon.project.data_store import DataStore
from spinDecon.analysis.restricted_3d import (analyse_overlap, classify_overlap, analyse_restricted_3d, classify_restricted_3d)
from spinDecon.project.parameter_store import (
    parse_all_strings as _parse_all_strings,
    parse_float as _parse_float,
    parse_int as _parse_int,
    parse_value as _parse_value,
    update_parameter_file,
    remove_parameter_keys,
)
from spinDecon.project.defaults import UNIDEC_DEFAULTS, available_cpu_count, is_default_value

def ParseFlt(infile,param,default=0.):
    return _parse_float(infile, param, default=default)

def Parse(infile,param,default=''):
    return _parse_value(infile, param, default=default)

def ParseInt(infile,param,default=0):
    return _parse_int(infile, param, default=default)

def ParseAllStr(infile,param):
    return _parse_all_strings(infile, param)

def findnear_index(test,array):
    #array = numpy.asarray(array)
    idx = (numpy.abs(array - test)).argmin()
    return idx

def findmax(array,col):
    test=float(array[0][col])
    imax=0
    for i in range(len(array)):
        if(float(array[i][col])>test):
            test=float(array[i][col])
            imax=i
    return imax

def switchy(Xmin,Xmax):
    return Xmax,Xmin

def readpeaklist(infile):
    peak=[]
    peakfile=open(infile,'r')
    for line in peakfile.readlines():
        test=line.split()
        if(len(test)>0):
            try:
                a=float(test[1])
                peak.append(peakEntry(test))
            except:
                pass
    peakfile.close()
    return peak



class Restricted3DDiagnosticsFrame(wx.Frame):
    """First-pass exception browser for 2D-restricted 3D diagnostics."""
    FILTERS = ("Needs attention", "All", "Heavy overlap", "Missing trace peaks",
               "Needs new XY slice", "Localisation", "OK")

    def __init__(self, decon_frame):
        super().__init__(decon_frame, title="Restricted 3D diagnostics", size=(1000, 620))
        self.decon_frame = decon_frame
        # ProjectState is authoritative at the housekeeping boundary.  Refresh
        # the legacy path controls before mirroring them into this editor.
        decon_frame._apply_state_to_path_controls()
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)

        params = wx.StaticBoxSizer(wx.StaticBox(panel, label="Analysis parameters - grouped by question"), wx.VERTICAL)
        self.controls = {}

        groups = (
            ("1. Heavy overlap",
             "Is a picked 3D component only a small part of the raw intensity?",
             (("Peak/raw fraction <=", "overlap_fraction", "0.05"),)),
            ("2. Missing peaks in the 1D trace",
             "Is coherent raw signal left unexplained by the deconvolution?",
             (("Residual sigma", "residual_sigma_threshold", "3.0"),
              ("Residual fraction", "residual_fraction_warning", "0.05"),
              ("Neighbour J radius", "residual_xy_radius_j", "1"),
              ("Neighbour K radius", "residual_xy_radius_k", "1"),
              ("Minimum supporting traces", "residual_min_support_traces", "3"))),
            ("3. Does this peak need a new XY slice?",
             "Test whether the present XY bore passes through the peak maximum, and estimate where a better XY source lies.",
             (("Z maximum radius", "z_search_radius", "5"),
              ("Minimum maximum fraction", "maximum_fraction_warning", "0.90"),
              ("XY search J radius", "xy_search_radius_j", "3"),
              ("XY search K radius", "xy_search_radius_k", "3"),
              ("XY displacement warning", "xy_displacement_warning", "1.5"))),
        )

        group_row = wx.BoxSizer(wx.HORIZONTAL)
        for title, description, fields in groups:
            box = wx.StaticBoxSizer(wx.StaticBox(panel, label=title), wx.VERTICAL)
            desc = wx.StaticText(panel, label=description)
            desc.Wrap(280)
            box.Add(desc, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 5)
            grid = wx.FlexGridSizer(len(fields), 2, 3, 6)
            grid.AddGrowableCol(0, 1)
            for label, key, value in fields:
                grid.Add(wx.StaticText(panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
                ctrl = wx.TextCtrl(panel, value=value, size=(65, -1))
                self.controls[key] = ctrl
                grid.Add(ctrl, 0, wx.ALIGN_RIGHT)
            box.Add(grid, 0, wx.EXPAND | wx.ALL, 5)
            group_row.Add(box, 1, wx.EXPAND | wx.RIGHT, 5)
        params.Add(group_row, 0, wx.EXPAND | wx.ALL, 5)

        button_box = wx.BoxSizer(wx.HORIZONTAL)
        self.analyse_btn = wx.Button(panel, label="Analyse")
        self.explore_xy_btn = wx.Button(panel, label="Explore proposed XY sources")
        self.help_btn = wx.Button(panel, label="Help: choosing parameters")
        button_box.AddStretchSpacer(1)
        button_box.Add(self.analyse_btn, 0, wx.RIGHT, 5)
        button_box.Add(self.explore_xy_btn, 0, wx.RIGHT, 5)
        button_box.Add(self.help_btn, 0)
        params.Add(button_box, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        outer.Add(params, 0, wx.EXPAND | wx.ALL, 6)

        bar = wx.BoxSizer(wx.HORIZONTAL)
        bar.Add(wx.StaticText(panel, label="Show:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.filter = wx.Choice(panel, choices=list(self.FILTERS)); self.filter.SetSelection(0)
        bar.Add(self.filter, 0, wx.RIGHT, 12)
        self.summary = wx.StaticText(panel, label="Press Analyse to calculate diagnostics.")
        bar.Add(self.summary, 1, wx.ALIGN_CENTER_VERTICAL)
        outer.Add(bar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        self.table = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        columns = (("Source", 100), ("3D peaks", 70), ("Overlap", 85), ("Worst max", 85),
                   ("XY disp", 75), ("Residual", 85), ("Flags", 360))
        for i, (name, width) in enumerate(columns): self.table.InsertColumn(i, name, width=width)
        outer.Add(self.table, 1, wx.EXPAND | wx.ALL, 6)
        panel.SetSizer(outer)
        self.analyse_btn.Bind(wx.EVT_BUTTON, self.OnAnalyse)
        self.explore_xy_btn.Bind(wx.EVT_BUTTON, self.OnExploreXY)
        self.help_btn.Bind(wx.EVT_BUTTON, self.OnHelp)
        self.filter.Bind(wx.EVT_CHOICE, lambda evt: self.RefreshTable())
        self.table.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.RefreshTable()

    def _params(self):
        f = lambda k: float(self.controls[k].GetValue())
        i = lambda k: int(float(self.controls[k].GetValue()))
        return dict(overlap_fraction=f("overlap_fraction"), z_search_radius=i("z_search_radius"),
                    maximum_fraction_warning=f("maximum_fraction_warning"), xy_search_radius_j=i("xy_search_radius_j"),
                    xy_search_radius_k=i("xy_search_radius_k"), xy_displacement_warning=f("xy_displacement_warning"),
                    residual_sigma_threshold=f("residual_sigma_threshold"), residual_fraction_warning=f("residual_fraction_warning"),
                    residual_xy_radius_j=i("residual_xy_radius_j"), residual_xy_radius_k=i("residual_xy_radius_k"),
                    residual_min_support_traces=i("residual_min_support_traces"))

    def OnAnalyse(self, event):
        try:
            result = self.decon_frame.analyse_restricted_3d(**self._params())
            self.summary.SetLabel(f"{result['reference_sources']} sources; {result['attention_sources']} need attention; "
                                  f"{result['overlapped_sources']} overlap; {result['localisation_sources']} localisation; "
                                  f"{result['residual_sources']} unexplained; {result['unmatched_records']} unmatched 3D peaks")
            self.RefreshTable()
        except Exception as exc:
            wx.MessageBox(str(exc), "Diagnostics", wx.OK | wx.ICON_ERROR, self)


    def OnExploreXY(self, event):
        strong = [p for p in self.decon_frame.get_reference_peaks()
                  if getattr(p, "analysis", {}).get("restricted_3d", {}).get("classification", {}).get("strong_missing_xy", False)]
        if not strong:
            wx.MessageBox("The last diagnostics analysis contains no strong missing-XY sources.",
                          "Explore proposed XY sources", wx.OK | wx.ICON_INFORMATION, self)
            return
        frame = getattr(self.decon_frame, 'missingXYExplorerFrame', None)
        try:
            if frame is not None and frame:
                frame.Raise(); frame.Show(True); return
        except Exception:
            pass
        frame = MissingXYExplorerFrame(self.decon_frame, parent=self)
        self.decon_frame.missingXYExplorerFrame = frame
        frame.Show(True)

    def OnHelp(self, event):
        text = """Restricted 3D diagnostics - what question are you trying to answer?

The controls are grouped around three separate scientific problems. They should not be treated as one global quality score. Start with the question you care about, tune only that group if possible, select the resulting sources in the table, and inspect them in the normal 1D/2D slice viewers.

============================================================
1. HEAVY OVERLAP
Question: Is this picked 3D component only a small fraction of the raw signal at its position?
============================================================

Peak/raw fraction <=
The peak-list intensity is divided by the magnitude of the raw intensity at the restricted XY location. A peak is called heavily overlapped when that fraction is at or below this value.

Default: 0.05 means the picked component accounts for only 5% or less of the raw intensity.

Raise it (for example 0.10 or 0.20) to include progressively less severe overlap. Lower it to restrict the list to only the most extreme cases.

Use the 'Heavy overlap' filter to inspect this question on its own. Heavy overlap is important context: localisation measurements on a severely overlapped peak should be interpreted cautiously.

============================================================
2. MISSING PEAKS IN THE 1D TRACE
Question: Did the deconvolution threshold leave a real weak peak unexplained?
============================================================

This test looks at raw-minus-deconvolved intensity and requires the feature to persist at the same Z position in several neighbouring XY traces. That spatial agreement is intended to reject isolated baseline/noise fluctuations.

Residual sigma
Each neighbouring trace must contain residual intensity above this many measured noise sigma. Raise it to be conservative. Lower it (for example from 3.0 toward 2.5 or 2.0) when deliberately hunting weak missed peaks. When lowering it, consider increasing Minimum supporting traces to retain protection against random fluctuations.

Residual fraction
Controls how much coherent unexplained signal is required before the whole source is flagged. Lower it to find subtler missing signal; raise it if too many weak residual features are being reported.

Neighbour J/K radius
Controls how many adjacent XY traces are examined. Radius 1/1 gives up to a 3 x 3 neighbourhood. Increase for features that are broad in XY. Avoid making the neighbourhood so large that unrelated nearby sources routinely enter it.

Minimum supporting traces
The number of neighbouring traces that must independently exceed Residual sigma at the same Z position. Increase this when baseline fluctuations produce false positives. Decrease it for extremely narrow XY features, or when working close to an edge where fewer neighbouring traces exist.

A useful weak-peak search is: lower Residual sigma gradually while increasing Minimum supporting traces. A real peak should remain coherent across neighbouring traces whereas random excursions should not.

Use the 'Missing trace peaks' filter for this question.

============================================================
3. DOES A PEAK NEED A NEW XY SLICE?
Question: Does the existing 2D source pass through the real 3D peak, or should this peak be owned by a different/new XY source? Where should that source be placed?
============================================================

This test combines two observations. First, a clean peak should normally be near a local maximum along its restricted Z trace. Second, at that Z position the program measures the intensity-weighted centre of the surrounding XY signal. This preserves the direction in which a broad/unresolved peak is leaning even when several slices share the same discrete XY maximum. Consistent lean directions across associated 3D peaks are pooled into candidate XY sources.

Z maximum radius
Number of Z points searched on either side of the picked peak. Increase for broad peaks or when picking can be displaced by several digital points. Keep it small enough that the search does not jump to a different neighbouring Z peak.

Minimum maximum fraction
The raw intensity at the picked Z position is compared with the strongest point in that Z window. 0.90 means the picked position should reach at least 90% of the nearby maximum. Raise it to demand more exact centring; lower it to tolerate broad/flat maxima and digital resolution.

XY search J/K radius
Defines how far around the current 2D source the program is allowed to search for the XY maximum at each associated Z peak. These values directly limit how far the analysis can look for a candidate new slice. Increase them when you suspect hidden sources are farther apart. Very large windows risk finding an unrelated source.

XY displacement warning
Distance in XY data points from the current source to the discovered local maximum before the source is considered displaced. Lower it to detect subtle miscentring; raise it to ignore ordinary peak-list/digital-grid error.

Interpretation matters here:
- One or more non-overlapped peaks consistently pointing toward the same nearby XY maximum suggests the existing 2D source is miscentred.
- Different groups of clean 3D peaks leaning in different directions is stronger evidence that the projected 2D peak represented more than one underlying source, even if there are no resolved XY maxima.
- A candidate that coincides with a different existing 2D source is suppressed: it is evidence for that existing source, not a new one.
- The reported centroid-based location is a candidate to inspect in the 2D slice viewer, not an automatic instruction to modify the peak list.

Use 'Needs new XY slice' to concentrate on the strongest current evidence, and 'Localisation' for the broader set of localisation warnings.

============================================================
PRACTICAL WORKFLOW
============================================================

1. Run the defaults and leave the table on Needs attention.
2. Ask one question at a time using Heavy overlap, Missing trace peaks, or Needs new XY slice.
3. Select a row. The shared peak selection should drive the existing 1D and 2D slice windows so you can inspect the evidence in context.
4. Tune the parameters for that problem and Analyse again. Do not feel obliged to tune unrelated groups.
5. Decide whether to accept the result, change the deconvolution threshold, move/add a 2D source, or leave a genuinely complex overlapped case for manual treatment.
6. Recompute the underlying result after making scientific changes, then rerun diagnostics.

The diagnostics are intended to prioritise inspection and provide evidence. They deliberately do not modify the 2D peak list or deconvolution automatically."""
        dlg = wx.Dialog(self, title="Restricted 3D diagnostics help", size=(760, 650), style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        panel = wx.Panel(dlg)
        box = wx.BoxSizer(wx.VERTICAL)
        help_text = wx.TextCtrl(panel, value=text, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_WORDWRAP)
        box.Add(help_text, 1, wx.EXPAND | wx.ALL, 8)
        close = wx.Button(panel, wx.ID_OK, "Close")
        box.Add(close, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        panel.SetSizer(box)
        dlg.ShowModal()
        dlg.Destroy()

    def _include(self, flags):
        choice = self.filter.GetStringSelection() or "Needs attention"
        return {"All": True, "Needs attention": flags.get("needs_attention", False),
                "Heavy overlap": flags.get("overlap", False),
                "Missing trace peaks": flags.get("unexplained", False),
                "Needs new XY slice": flags.get("possible_missing_xy", False),
                "Localisation": flags.get("localisation", False),
                "OK": not flags.get("needs_attention", False)}[choice]

    @staticmethod
    def _fmt(value, digits=3):
        return "-" if value is None else f"{float(value):.{digits}g}"

    def RefreshTable(self):
        self.table.DeleteAllItems(); self._row_names = []
        for peak in self.decon_frame.get_reference_peaks():
            src = getattr(peak, "analysis", {}).get("restricted_3d", {})
            flags = src.get("classification", {})
            if not src or not self._include(flags): continue
            ov, loc, res = src.get("overlap", {}), src.get("localisation", {}), src.get("residual", {})
            labels = []
            if flags.get("overlap"): labels.append("overlap")
            if flags.get("localisation"): labels.append("localisation")
            if flags.get("unexplained"): labels.append("unexplained")
            if flags.get("strong_missing_xy"): labels.append("missing XY (strong)")
            elif flags.get("possible_missing_xy"): labels.append("missing XY?")
            if not labels: labels.append("OK")
            row = self.table.InsertItem(self.table.GetItemCount(), str(getattr(peak, "name", "")))
            vals = (str(src.get("count", 0)), self._fmt(ov.get("worst_fraction")),
                    self._fmt(loc.get("worst_maximum_fraction")), self._fmt(loc.get("max_xy_displacement")),
                    self._fmt(res.get("max_sigma")), ", ".join(labels))
            for col, value in enumerate(vals, 1): self.table.SetItem(row, col, value)
            self._row_names.append(str(getattr(peak, "name", "")))

    def OnSelect(self, event):
        row = event.GetIndex()
        if 0 <= row < len(self._row_names): self.decon_frame.select_reference_peak(self._row_names[row])

    def OnClose(self, event):
        self.decon_frame.restricted3dDiagnosticsFrame = None
        self.Destroy()


class MissingXYExplorerFrame(wx.Frame):
    """Inspect candidate replacement/additional XY sources from the last diagnostics run."""

    def __init__(self, decon_frame, parent=None):
        super().__init__(parent or decon_frame, title="Explore proposed XY sources", size=(820, 720))
        self.decon_frame = decon_frame
        self._sources = self._strong_sources()

        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        top = wx.BoxSizer(wx.HORIZONTAL)
        top.Add(wx.StaticText(panel, label="Strong missing-XY source:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.source_choice = wx.ComboBox(panel, style=wx.CB_READONLY,
                                         choices=[str(getattr(p, "name", "")) for p in self._sources])
        if self._sources:
            self.source_choice.SetSelection(0)
        top.Add(self.source_choice, 1, wx.RIGHT, 8)
        self.status = wx.StaticText(panel, label="")
        top.Add(self.status, 1, wx.ALIGN_CENTER_VERTICAL)
        outer.Add(top, 0, wx.EXPAND | wx.ALL, 8)

        self.fig = Figure(figsize=(7, 6))
        self.canvas = FigCanvas(panel, -1, self.fig)
        outer.Add(self.canvas, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)

        note = wx.StaticText(panel, label=(
            "Blue crosses are the current Reference 2D sources and are labelled with their existing source names. "
            "Orange open circles are proposed XY sources. The proposal nearest the selected source keeps that "
            "source name; any additional proposals receive new source names."))
        note.Wrap(780)
        outer.Add(note, 0, wx.EXPAND | wx.ALL, 8)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.AddStretchSpacer(1)
        self.accept_btn = wx.Button(panel, label="Accept")
        self.accept_btn.SetToolTip("Reserved for the later peak-list editing workflow; no changes are made yet.")
        self.close_btn = wx.Button(panel, label="Close")
        buttons.Add(self.accept_btn, 0, wx.RIGHT, 6)
        buttons.Add(self.close_btn, 0)
        outer.Add(buttons, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        panel.SetSizer(outer)

        self.source_choice.Bind(wx.EVT_COMBOBOX, self.OnSourceChanged)
        self.accept_btn.Bind(wx.EVT_BUTTON, self.OnAccept)
        self.close_btn.Bind(wx.EVT_BUTTON, lambda evt: self.Close())
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Draw()

    def _strong_sources(self):
        out = []
        for peak in self.decon_frame.get_reference_peaks():
            src = getattr(peak, "analysis", {}).get("restricted_3d", {})
            if src.get("classification", {}).get("strong_missing_xy", False):
                out.append(peak)
        return out

    def _selected_source(self):
        idx = self.source_choice.GetSelection()
        return self._sources[idx] if 0 <= idx < len(self._sources) else None

    def _proposed_names(self, count):
        """Generate ``count`` unused names for genuinely additional sources."""
        if count <= 0:
            return []
        current = [str(getattr(p, "name", "")).strip() for p in self.decon_frame.get_reference_peaks()]
        current = [name for name in current if name]
        used = set(current)
        last = current[-1] if current else "Peak0"
        matches = list(re.finditer(r"\d+", last))
        if matches:
            m = matches[-1]
            prefix, suffix = last[:m.start()], last[m.end():]
            width = len(m.group(0))
            number = int(m.group(0))
            def make_name(n):
                return prefix + str(n).zfill(width) + suffix
        else:
            number = 0
            def make_name(n):
                return last + str(n)
        names = []
        while len(names) < count:
            number += 1
            candidate = make_name(number)
            if candidate not in used:
                names.append(candidate)
                used.add(candidate)
        return names

    def _candidate_locations(self, peak):
        """Return directional XY proposals supported by the associated 3D peaks.

        We deliberately use the continuous intensity centroid stored by the
        diagnostic, rather than the discrete plane maximum.  Unresolved sources
        can share one maximum while their individual Z peaks lean in different
        directions.  Similar displacement directions are pooled into one source
        proposal.  A proposal that lands on another existing 2D source is not a
        new source and is suppressed.
        """
        store = self.decon_frame.store
        payload = store.get_peak_list("full")
        records = list(payload.get("records") or payload.get("peaks") or [])
        src = getattr(peak, "analysis", {}).get("restricted_3d", {})
        jscale = numpy.asarray(getattr(store.uc1, "ppms_scale", store.index1), dtype=float)
        kscale = numpy.asarray(getattr(store.uc2, "ppms_scale", store.index2), dtype=float)
        try:
            pj = float(numpy.argmin(numpy.abs(jscale - float(peak.y))))
            pk = float(numpy.argmin(numpy.abs(kscale - float(peak.x))))
        except Exception:
            return []

        observations = []
        for ri in src.get("full_peak_indices", []):
            if not isinstance(ri, int) or not (0 <= ri < len(records)):
                continue
            rec = records[ri]
            diag = rec.get("analysis", {}).get("restricted_3d", {})
            ov, loc = diag.get("overlap", {}), diag.get("localisation", {})
            if ov.get("is_overlapped") or not loc.get("is_warning"):
                continue
            idx = loc.get("xy_centroid_index") or loc.get("xy_max_index")
            if not idx or len(idx) != 2:
                continue
            cj, ck = float(idx[0]), float(idx[1])
            v = numpy.asarray([cj-pj, ck-pk], dtype=float)
            r = float(numpy.linalg.norm(v))
            if not numpy.isfinite(r) or r <= 0.20:
                continue
            observations.append([cj, ck, v, r, str(rec.get("name", ri))])

        # Cluster by direction, not exact grid coordinate.  This captures the
        # scientifically useful pattern of one set of Z peaks leaning one way
        # and another set leaning another way even when neither has a clean XY
        # maximum.  cos(60 deg)=0.5 is deliberately broad within a direction.
        groups = []
        for obs in sorted(observations, key=lambda o: -o[3]):
            unit = obs[2] / obs[3]
            target = None
            best_cos = 0.5
            for group in groups:
                gu = group['vector'] / max(float(numpy.linalg.norm(group['vector'])), 1e-12)
                cosine = float(numpy.dot(unit, gu))
                if cosine > best_cos:
                    best_cos, target = cosine, group
            if target is None:
                groups.append({'items': [obs], 'vector': obs[2].copy()})
            else:
                target['items'].append(obs)
                target['vector'] += obs[2]

        refs = self.decon_frame.get_reference_peaks()
        proposals = []
        for group in groups:
            items = group['items']
            # More strongly displaced slices carry more information about which
            # hidden source dominates, so weight their centroids accordingly.
            w = numpy.asarray([max(o[3], 0.25) for o in items], dtype=float)
            cj = float(numpy.average([o[0] for o in items], weights=w))
            ck = float(numpy.average([o[1] for o in items], weights=w))
            x = float(numpy.interp(ck, numpy.arange(kscale.size), kscale))
            y = float(numpy.interp(cj, numpy.arange(jscale.size), jscale))

            # Never propose a "new" source on top of a different existing 2D
            # peak.  Compare in data-point units so unequal ppm axis scaling does
            # not distort the test.  The selected source itself is exempt: its
            # nearest proposal is allowed to be its replacement/moved position.
            duplicate = False
            for ref in refs:
                if ref is peak:
                    continue
                rj = float(numpy.argmin(numpy.abs(jscale - float(ref.y))))
                rk = float(numpy.argmin(numpy.abs(kscale - float(ref.x))))
                if float(numpy.hypot(cj-rj, ck-rk)) <= 1.0:
                    duplicate = True
                    break
            if duplicate:
                continue
            proposals.append((x, y, [o[4] for o in items], (cj, ck)))
        return proposals

    def _named_candidates(self, peak, candidates):
        """Name proposals so the replacement nearest the old source keeps its identity.

        The nearest proposal is the continuation/replacement of the source under
        investigation and therefore keeps that source name.  Only additional
        XY sources receive newly generated names.
        """
        if not candidates:
            return []
        px, py = float(peak.x), float(peak.y)
        nearest = min(range(len(candidates)),
                      key=lambda i: (candidates[i][0] - px) ** 2 + (candidates[i][1] - py) ** 2)
        new_names = iter(self._proposed_names(max(0, len(candidates) - 1)))
        original_name = str(getattr(peak, 'name', '')).strip() or 'Source'
        return [(original_name if i == nearest else next(new_names), candidate)
                for i, candidate in enumerate(candidates)]

    def _projection(self):
        labels = list(getattr(self.decon_frame, "labb", []) or [])
        if len(labels) < 3:
            raise RuntimeError("A 3D projection requires three axis labels.")
        view = self.decon_frame.get_projection_view(labels[2], labels[1], decon=False, transpose='n')
        if view is None or view.get("ZZ") is None:
            raise RuntimeError("The raw 2D projection used by View Peaks is not available.")
        return view, labels[2], labels[1]

    def OnSourceChanged(self, event):
        self.Draw()

    def OnAccept(self, event):
        wx.MessageBox("Accept is intentionally not connected yet. No peak-list changes have been made.",
                      "Proposed XY sources", wx.OK | wx.ICON_INFORMATION, self)

    def Draw(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        peak = self._selected_source()
        if peak is None:
            ax.text(0.5, 0.5, "No strong missing-XY sources are present in the last analysis.",
                    ha='center', va='center', transform=ax.transAxes)
            self.status.SetLabel("0 targets")
            self.canvas.draw_idle()
            return
        try:
            view, xlabel, ylabel = self._projection()
            XX, YY, ZZ = view.get("XX"), view.get("YY"), numpy.asarray(view.get("ZZ"), dtype=float)
            finite = numpy.abs(ZZ[numpy.isfinite(ZZ)])
            vmax = float(numpy.max(finite)) if finite.size else 0.0
            if vmax > 0:
                # Similar logarithmic contour presentation to View Peaks, but local and read-only.
                low = max(vmax * 0.01, numpy.finfo(float).eps)
                levels = numpy.geomspace(low, vmax, 18)
                ax.contour(XX, YY, ZZ, levels=levels, linewidths=0.65)

            refs = self.decon_frame.get_reference_peaks()
            if refs:
                current_artist = ax.scatter(
                    [float(p.x) for p in refs], [float(p.y) for p in refs],
                    marker='x', s=35, color='tab:blue', label='Current 2D sources', zorder=4)
                # Label the existing sources at their actual XY positions.  Labels use a small
                # screen-space offset so they remain attached correctly when the axes are reversed.
                for ref in refs:
                    ax.annotate(str(getattr(ref, 'name', '')),
                                (float(ref.x), float(ref.y)), xytext=(5, 5),
                                textcoords='offset points', fontsize=8, zorder=7)

            candidates = self._candidate_locations(peak)
            named_candidates = self._named_candidates(peak, candidates)
            if candidates:
                xs = [c[0] for c in candidates]; ys = [c[1] for c in candidates]
                proposed_artist = ax.scatter(
                    xs, ys, marker='o', s=80, linewidths=1.8,
                    facecolors='none', edgecolors='tab:orange',
                    label='Proposed XY sources', zorder=5)
                for proposed_name, (x, y, supporting_peaks, _idx) in named_candidates:
                    ax.annotate(proposed_name, (x, y), xytext=(6, -11),
                                textcoords='offset points', fontsize=9, fontweight='bold',
                                color='tab:orange', zorder=8)

            # Mark the selected current source more prominently, but do not add another legend entry.
            ax.scatter([float(peak.x)], [float(peak.y)], marker='x', s=100, linewidths=2.0, color='tab:blue', zorder=6)

            # Local view includes the source, all proposals, and a margin based on the analysis search window.
            xvals = [float(peak.x)] + [c[0] for c in candidates]
            yvals = [float(peak.y)] + [c[1] for c in candidates]
            xaxis = numpy.asarray(view.get('x_scale') if view.get('x_scale') is not None else numpy.asarray(XX)[0, :], dtype=float)
            yaxis = numpy.asarray(view.get('y_scale') if view.get('y_scale') is not None else numpy.asarray(YY)[:, 0], dtype=float)
            xstep = float(numpy.nanmedian(numpy.abs(numpy.diff(xaxis)))) if xaxis.size > 1 else 0.01
            ystep = float(numpy.nanmedian(numpy.abs(numpy.diff(yaxis)))) if yaxis.size > 1 else 0.01
            params = getattr(self.decon_frame.store, 'analysis', {}).get('restricted_3d', {}).get('parameters', {})
            mx = max(4.0 * xstep, (float(params.get('xy_search_radius_k', 3)) + 2.0) * xstep)
            my = max(4.0 * ystep, (float(params.get('xy_search_radius_j', 3)) + 2.0) * ystep)
            xmin, xmax = min(xvals)-mx, max(xvals)+mx
            ymin, ymax = min(yvals)-my, max(yvals)+my
            # NMR axes conventionally descend in ppm.
            ax.set_xlim(max(xmin, xmax), min(xmin, xmax))
            ax.set_ylim(max(ymin, ymax), min(ymin, ymax))
            ax.set_xlabel(f"{xlabel} (ppm)")
            ax.set_ylabel(f"{ylabel} (ppm)")
            ax.set_title(f"Source {getattr(peak, 'name', '')}: local 2D projection")
            ax.legend(loc='best', fontsize=8)
            self.status.SetLabel(f"{len(candidates)} proposed location(s)")
            self.fig.tight_layout()
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), ha='center', va='center', transform=ax.transAxes, wrap=True)
            self.status.SetLabel("plot unavailable")
        self.canvas.draw_idle()

    def OnClose(self, event):
        self.decon_frame.missingXYExplorerFrame = None
        self.Destroy()


class ReOrganiseFrame(wx.Frame):
    """Transient command window for spectrum reorganisation operations."""
    def __init__(self, decon_frame):
        wx.Frame.__init__(self, decon_frame, title="ReOrganise", style=wx.DEFAULT_FRAME_STYLE & ~wx.RESIZE_BORDER & ~wx.MAXIMIZE_BOX)
        self.decon_frame = decon_frame
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.buttonXY = wx.Button(panel, label="X <-> Y (TP)", size=(150, 28))
        self.buttonXZ = wx.Button(panel, label="X <-> Z (ZTP)", size=(150, 28))
        self.buttonXA = wx.Button(panel, label="X <-> A (ATP)", size=(150, 28))
        self.buttonCirc = wx.Button(panel, label="circ", size=(150, 28))
        self.buttonProject = wx.Button(panel, label="Project", size=(150, 28))
        self.buttonClose = wx.Button(panel, label="Close", size=(150, 28))

        for button in (self.buttonXY, self.buttonXZ, self.buttonXA, self.buttonCirc, self.buttonProject, self.buttonClose):
            sizer.Add(button, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 8)
        sizer.AddSpacer(8)
        panel.SetSizer(sizer)

        self.buttonXY.Bind(wx.EVT_BUTTON, decon_frame.OnButtonXY)
        self.buttonXZ.Bind(wx.EVT_BUTTON, decon_frame.OnButtonXZ)
        self.buttonXA.Bind(wx.EVT_BUTTON, decon_frame.OnButtonXA)
        self.buttonCirc.Bind(wx.EVT_BUTTON, decon_frame.OnButtonCirc)
        self.buttonProject.Bind(wx.EVT_BUTTON, decon_frame.OnButtonProject)
        self.buttonClose.Bind(wx.EVT_BUTTON, self.OnClose)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.update_enabled_state()
        self.Fit()
        self.CentreOnParent()

    def update_enabled_state(self):
        has_pipe = bool(self.decon_frame.has_pipe)
        dim = int(self.decon_frame.dim)
        self.buttonXY.Enable(has_pipe)
        self.buttonXZ.Enable(has_pipe and dim >= 3)
        self.buttonXA.Enable(has_pipe and dim >= 4)
        self.buttonCirc.Enable(has_pipe and dim >= 4)
        self.buttonProject.Enable(has_pipe)

    def OnClose(self, event):
        self.decon_frame.reorganiseFrame = None
        self.Destroy()

class HousekeepingFrame(wx.Frame):
    """Temporary editor for the NMR tab's canonical project path controls.

    The NMR tab remains the owner of all values.  This frame only mirrors them
    while it is open; Save applies the edits to the NMR controls/state and then
    uses the normal project save routine (the single decon parameter file).
    """
    def __init__(self, decon_frame):
        super().__init__(decon_frame, title='Housekeeping', style=wx.DEFAULT_FRAME_STYLE & ~wx.RESIZE_BORDER)
        self.decon_frame = decon_frame
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.FlexGridSizer(rows=4, cols=3, vgap=6, hgap=6)
        grid.AddGrowableCol(1, 1)

        self.dirBox = wx.TextCtrl(panel, value=decon_frame.dirBox.GetValue(), size=(330, -1))
        self.outPathBox = wx.TextCtrl(panel, value=decon_frame.outPathBox.GetValue(), size=(330, -1))
        self.specPathBox = wx.TextCtrl(panel, value=decon_frame.specPathBox.GetValue(), size=(330, -1))
        self.infileBox = wx.TextCtrl(panel, value=decon_frame.infileBox.GetValue(), size=(330, -1))

        rows = (
            ('Working dir:', self.dirBox, self._choose_working_dir),
            ('OutPath:', self.outPathBox, lambda evt: self._set_directory(self.outPathBox)),
            ('SpecPath:', self.specPathBox, lambda evt: self._set_directory(self.specPathBox)),
            ('nmrPipe file:', self.infileBox, self._choose_nmrpipe_file),
        )
        for label, box, handler in rows:
            grid.Add(wx.StaticText(panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(box, 1, wx.EXPAND)
            button = wx.Button(panel, label='...', size=(40, 24))
            button.Bind(wx.EVT_BUTTON, handler)
            grid.Add(button, 0, wx.ALIGN_CENTER_VERTICAL)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.AddStretchSpacer(1)
        save_btn = wx.Button(panel, label='Save')
        close_btn = wx.Button(panel, label='Close')
        save_btn.Bind(wx.EVT_BUTTON, self.OnSave)
        close_btn.Bind(wx.EVT_BUTTON, self.OnClose)
        buttons.Add(save_btn, 0, wx.RIGHT, 6)
        buttons.Add(close_btn, 0)

        outer.Add(grid, 1, wx.EXPAND | wx.ALL, 12)
        outer.Add(buttons, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        panel.SetSizer(outer)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Fit()
        self.CentreOnParent()

    def _choose_working_dir(self, event):
        dlg = wx.DirDialog(self, message='Choose working directory', defaultPath=self.dirBox.GetValue() or os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            self.dirBox.SetValue(dlg.GetPath())
        dlg.Destroy()

    def _choose_directory(self, event):
        # OutPath and SpecPath are project paths rooted at WorkingDir.  Store a
        # relative path when the selected directory is beneath WorkingDir.
        box = self.outPathBox if event.GetEventObject().GetId() == getattr(self, '_out_button_id', -1) else None
        # Determine target from button row by screen position is brittle; bind
        # wrappers below instead.  This method remains for compatibility.

    def _set_directory(self, box):
        working = os.path.abspath(self.dirBox.GetValue() or os.getcwd())
        current = box.GetValue().strip()
        default = current if os.path.isabs(current) else os.path.join(working, current or '.')
        dlg = wx.DirDialog(self, message='Choose folder', defaultPath=os.path.abspath(default))
        if dlg.ShowModal() == wx.ID_OK:
            chosen = os.path.abspath(dlg.GetPath())
            try:
                common = os.path.commonpath([chosen, working])
            except ValueError:
                common = ''
            value = os.path.relpath(chosen, working) if common == working else chosen
            if value == '.': value = './'
            elif not os.path.isabs(value): value = './' + value.replace(os.sep, '/')
            box.SetValue(value)
        dlg.Destroy()

    def _choose_nmrpipe_file(self, event):
        working = os.path.abspath(self.dirBox.GetValue() or os.getcwd())
        spec = self.specPathBox.GetValue().strip() or './spec'
        spec_dir = os.path.abspath(spec if os.path.isabs(spec) else os.path.join(working, spec))
        dlg = wx.FileDialog(self, message='Choose nmrPipe file in SpecPath', defaultDir=spec_dir,
                            wildcard='All files (*.*)|*.*', style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = os.path.abspath(dlg.GetPath())
            try:
                common = os.path.commonpath([path, spec_dir])
            except ValueError:
                common = ''
            if common != spec_dir:
                errorMessage('nmrPipe files must be located inside SpecPath: %s' % spec_dir)
            else:
                self.infileBox.SetValue(os.path.relpath(path, spec_dir).replace(os.sep, '/'))
        dlg.Destroy()

    def OnSave(self, event=None):
        owner = self.decon_frame
        # Apply directory roots first so the SpecPath-relative nmrPipe value is
        # canonicalised against the newly selected roots.
        owner.dirBox.SetValue(self.dirBox.GetValue().strip())
        # Do not mutate process cwd here.  All project resources are resolved
        # explicitly from WorkingDir + RawPath/SpecPath; changing cwd made
        # pseudo3D/Fitting behaviour depend on whether UniDecNMR was launched
        # directly or from SpinHubMain.
        owner.outPathBox.SetValue(self.outPathBox.GetValue().strip())
        owner.specPathBox.SetValue(self.specPathBox.GetValue().strip())
        owner._sync_directory_state_only()
        owner.infileBox.SetValue(owner.state._spec_relative(self.infileBox.GetValue().strip()))
        owner._sync_path_state()
        # There is deliberately one persistence path only: the existing NMR
        # project save routine writes the canonical decon parameter file.
        owner.OnButtonSave(True)
        owner.update_project_lamps()
        self.Destroy()
        owner.housekeepingFrame = None

    def OnClose(self, event=None):
        self.decon_frame.housekeepingFrame = None
        self.Destroy()

class NMRWorkspace(wx.Panel):

    def data_box(self):
        self.uSTA = False
        self.dataLbl = wx.StaticBox(self, -1, 'Data:')

        # File/path controls live together in the Data box.  Keep the
        # historical attribute names so save/load and processing code continue
        # to use exactly the same widgets.
        self.dirLab = wx.StaticText(self.dataLbl, label="Working dir:", size=(83, -1))
        self.dirBox = wx.TextCtrl(self.dataLbl, size=(200, 22))
        self.openDirFileBtn = wx.Button(self.dataLbl, label="...", size=(40, 22))
        self.openDirFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.dirBox))

        self.outPathLab = wx.StaticText(self.dataLbl, label="OutPath:", size=(83, -1))
        self.outPathBox = wx.TextCtrl(self.dataLbl, size=(200, 22))
        self.openOutPathBtn = wx.Button(self.dataLbl, label="...", size=(40, 22))
        self.openOutPathBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.outPathBox, change_cwd=False))

        self.specPathLab = wx.StaticText(self.dataLbl, label="SpecPath:", size=(83, -1))
        self.specPathBox = wx.TextCtrl(self.dataLbl, size=(200, 22))
        self.openSpecPathBtn = wx.Button(self.dataLbl, label="...", size=(40, 22))
        self.openSpecPathBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.specPathBox, change_cwd=False))

        self.infileLab = wx.StaticText(self.dataLbl, label="nmrPipe file:", size=(83, -1))
        self.infileBox = wx.TextCtrl(self.dataLbl, size=(200, 22))
        self.openNMRFileBtn = wx.Button(self.dataLbl, label="...", size=(40, 22))
        self.openNMRFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetSpecFile(evt, self.infileBox))

        self.sizer = wx.GridBagSizer(4, 4)
        self.sizer.Add((10, 0), (0, 0))
        self.sizer.Add((0, 10), (0, 1))
        self.sizer.Add(self.dirLab, (1, 1), border=1,
                       flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP)
        self.sizer.Add(self.dirBox, (1, 2), border=1,
                       flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP | wx.EXPAND)
        self.sizer.Add(self.openDirFileBtn, (1, 3), border=1,
                       flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP)
        for row, lab, box, btn in ((2, self.outPathLab, self.outPathBox, self.openOutPathBtn),
                                   (3, self.specPathLab, self.specPathBox, self.openSpecPathBtn),
                                   (4, self.infileLab, self.infileBox, self.openNMRFileBtn)):
            self.sizer.Add(lab, (row, 1), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT)
            self.sizer.Add(box, (row, 2), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.EXPAND)
            self.sizer.Add(btn, (row, 3), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT)
        self.sizer.AddGrowableCol(2, 1)

        self.dataSizer = wx.StaticBoxSizer(self.dataLbl, wx.VERTICAL)
        self.border = wx.BoxSizer(wx.VERTICAL)
        self.border.Add(self.sizer, 1, wx.EXPAND)
        self.dataSizer.Add(self.border, 1, wx.EXPAND)
        self.dataSizer.AddSpacer(5)

    def project_box(self):
        """Build the NMR project readiness panel without changing button handlers."""
        self.projectLbl = wx.StaticBox(self, -1, 'Project:')
        self.projectSizer = wx.StaticBoxSizer(self.projectLbl, wx.VERTICAL)

        # These are the original live buttons.  Detach them from their old
        # visual sizers and reparent them; their event bindings and attributes
        # are deliberately left unchanged.
        for sizer, button in (
            (self.spectrumButtonRow, self.buttonRead),
            (self.buttonSizer, self.buttonAnalyse),
            (self.sizerPk, self.buttonReadPeak),
            (self.sizerPk, self.buttonReadFullPeak),
        ):
            try:
                sizer.Detach(button)
            except Exception:
                pass
            button.Reparent(self.projectLbl)

        # Peak-list viewer buttons now live beside their corresponding Load
        # buttons in the Project box.  They were deliberately left out of the
        # Peak Lists sizer, so only reparenting is required here.
        self.buttonReferencePeakList.Reparent(self.projectLbl)
        self.buttonFullPeakList.Reparent(self.projectLbl)

        self.projectLamps = {}
        self.projectRowLabels = {}
        self.projectRows = {}
        rows = (
            ('spectrum', 'Spectrum file?', self.buttonRead),
            ('reference', 'Reference peak list?', self.buttonReadPeak),
            ('full', 'Full peak list?', self.buttonReadFullPeak),
            ('decon', 'Deconvolution file?', self.buttonAnalyse),
        )
        for key, label, button in rows:
            row = wx.BoxSizer(wx.HORIZONTAL)
            self.projectRows[key] = row
            lamp = wx.Panel(self.projectLbl, size=(14, 14), style=wx.BORDER_SIMPLE)
            lamp.SetMinSize((14, 14))
            self.projectLamps[key] = lamp
            # Status lamp comes first so readiness can be scanned down the left edge.
            row.Add(lamp, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
            row_label = wx.StaticText(self.projectLbl, label=label)
            self.projectRowLabels[key] = row_label
            row.Add(row_label, 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            if key == 'reference':
                row.Add(self.buttonReferencePeakList, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
            elif key == 'full':
                row.Add(self.buttonFullPeakList, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
            row.Add(button, 0, wx.ALIGN_CENTER_VERTICAL)
            self.projectSizer.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 8)

        projectButtonRow = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonLoadProject = wx.Button(self.projectLbl, label='Load Project')
        self.buttonLoadProject.Bind(wx.EVT_BUTTON, self.OnButtonLoadProject)
        projectButtonRow.Add(self.buttonLoadProject, 1, wx.EXPAND | wx.RIGHT, 4)

        self.buttonSummariseProject = wx.Button(self.projectLbl, label='Summarise Project')
        self.buttonSummariseProject.Bind(wx.EVT_BUTTON, self.OnButtonSummariseProject)
        projectButtonRow.Add(self.buttonSummariseProject, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)

        self.buttonHousekeeping = wx.Button(self.projectLbl, label='Housekeeping')
        self.buttonHousekeeping.Bind(wx.EVT_BUTTON, self.OnButtonHousekeeping)
        projectButtonRow.Add(self.buttonHousekeeping, 1, wx.EXPAND | wx.LEFT, 4)
        self.projectSizer.Add(projectButtonRow, 0, wx.EXPAND | wx.ALL, 8)

        # Refresh immediately whenever any path field is edited or populated.
        for ctrl in (self.infileBox, self.referencePeakBox, self.fullPeakBox, self.dirBox, self.outPathBox, self.specPathBox):
            ctrl.Bind(wx.EVT_TEXT, self._on_project_path_changed)
        self.update_project_lamps()

    def _on_project_path_changed(self, event):
        self.update_project_lamps()
        event.Skip()

    def OnButtonLoadProject(self, event=None):
        """Load the available project components in their dependency order.

        The project-path panel can be populated before the normal system-file
        Load action has synchronised the GUI dimension controls into
        ``ProjectState``.  Spectrum loading consults ``state.topology()``, so
        establish that small piece of canonical state here at the load
        boundary rather than relying on a previous button press.
        """
        pass

        if self.state is not None:
            # ProjectState already owns the declared scientific topology.  In
            # particular, legacy project loads can leave the old selector displaying a
            # physical dimension until makeinp() has inspected the NMR data and
            # canonicalised it.  Do not copy that transitional GUI value back
            # into state here: doing so turns pseudo2D (1 spectral + pseudo)
            # into pseudo3D immediately before the authoritative load boundary.
            self.state.sync_from_values(
                working_dir=self.dirBox.GetValue(),
                raw_path=self.outPathBox.GetValue(),
                spec_path=self.specPathBox.GetValue(),
                input_file=self.infileBox.GetValue(),
                reference_peak_file=self.referencePeakBox.GetValue(),
                full_peak_file=self.fullPeakBox.GetValue(),
                sym_mode=self.cb_grid.IsChecked(),
                decon_bore=self.cb_decon3d.IsChecked(),
            )

        spectrum = self._project_path(self.infileBox.GetValue())
        reference = self._project_path(self.referencePeakBox.GetValue())
        full = self._project_path(self.fullPeakBox.GetValue())
        # For pseudo2D the deconvolution product belongs to the generated 1D
        # projection.  Its exact path can only be resolved after the spectrum
        # has established the canonical topology/axis labels, so defer that
        # lookup until after OnButtonRead below.
        decon = spectrum + '.decon' if spectrum else ''

        # Reading the spectrum establishes the data shape and enables the
        # dependent controls, so it must always be attempted first.
        if spectrum and os.path.isfile(spectrum):
            pass
            self.OnButtonRead(True)
            pass

        if getattr(self, 'READ', 0) != 0:
            # The reference peak list is project state and should be loaded
            # whenever it is available.  cb_decon3d controls how the list is
            # used during deconvolution; it must not gate project loading.
            topology = self._active_topology()
            distinct_reference = not (topology.spectral_dim_count == 2 and not topology.has_pseudo_axis)
            if distinct_reference and reference and os.path.isfile(reference):
                self.OnButtonReadPeak(True)

            if full and os.path.isfile(full):
                self.OnButtonReadFullPeak(True)

            decon = self._active_deconvolution_path(spectrum)
            if decon and os.path.isfile(decon):
                self.OnButtonAnalyse(True)

        self.update_project_lamps()

    def OnButtonHousekeeping(self, event=None):
        """Open the temporary project-path editor."""
        frame = getattr(self, 'housekeepingFrame', None)
        if frame is not None:
            try:
                frame.Raise()
                return
            except Exception:
                self.housekeepingFrame = None
        self.housekeepingFrame = HousekeepingFrame(self)
        self.housekeepingFrame.Show()

    def OnButtonSummariseProject(self, event=None, show_viewer=True):
        """Generate a PDF report with the standard NMR magnet progress strip."""
        progress_frame = None
        try:
            from spinDecon.gui.reporting.project_summary import generate_project_summary, project_summary_stages
            from spinDecon.gui.dialogs.shell_output import ShellOutputFrame
            stages = project_summary_stages(self)
            progress_frame = ShellOutputFrame(self, title='Project Summary Progress')
            progress_frame.set_workflow(stages)
            progress_frame.progressHint.SetLabel('Progress follows the report stages required for this dataset.')
            progress_frame.Show()

            def report_progress(label):
                progress_frame.start_step(label)
                # Report export is intentionally synchronous because it reads
                # live wx/matplotlib windows.  Yield only at coarse stage
                # boundaries so the standard progress control can repaint.
                try:
                    wx.YieldIfNeeded()
                except Exception:
                    pass

            pdf_path, warnings = generate_project_summary(self, progress_callback=report_progress)
            progress_frame.finish_workflow(True)
            try:
                wx.YieldIfNeeded()
            except Exception:
                pass
        except Exception as exc:
            if progress_frame is not None:
                try:
                    progress_frame.finish_workflow(False)
                except Exception:
                    pass
            if event is None and not show_viewer:
                if progress_frame is not None:
                    progress_frame.Destroy()
                raise
            wx.MessageBox(str(exc), 'Project summary', wx.OK | wx.ICON_ERROR, parent=self)
            return None

        if show_viewer:
            try:
                from spinDecon.gui.dialogs.pdf_viewer import PDFViewer
                viewer = PDFViewer(self, size=(900, 700), title='Project Summary')
                viewer.LoadFile(pdf_path)
                viewer.Show()
                self._summary_pdf_viewer = viewer
            except Exception as exc:
                wx.MessageBox('Summary created at:\n%s\n\nViewer error: %s' % (pdf_path, exc),
                              'Project summary', wx.OK | wx.ICON_WARNING, parent=self)
        if progress_frame is not None:
            try:
                wx.CallLater(900, progress_frame.Destroy)
            except Exception:
                progress_frame.Destroy()
        return pdf_path

    def _project_path(self, value):
        """Resolve spectrum-associated files strictly below SpecPath."""
        return self._resolve_spec_file(value)

    def _apply_state_to_path_controls(self):
        """Mirror canonical ProjectState values into legacy UniDecNMR controls.

        These controls remain available to older callbacks, but they no longer
        invent the initial project paths.  ProjectState is the opening-time
        authority and housekeeping is a UI mirror of that state.
        """
        if getattr(self, 'state', None) is None:
            return
        self.dirBox.SetValue(self.state.working_dir or os.getcwd())
        self.outPathBox.SetValue(self.state.raw_path or './raw')
        self.specPathBox.SetValue(self.state.spec_path or './spec')
        self.infileBox.SetValue(self.state.input_file or '')

    def _sync_path_state(self):
        if getattr(self, 'state', None) is None:
            return
        self.state.sync_from_values(
            working_dir=self.dirBox.GetValue(), raw_path=self.outPathBox.GetValue(),
            spec_path=self.specPathBox.GetValue(), input_file=self.infileBox.GetValue(),
            reference_peak_file=self.referencePeakBox.GetValue(),
            full_peak_file=self.fullPeakBox.GetValue())

    def _resolve_spec_file(self, value):
        self._sync_directory_state_only()
        return self.state.resolve_spec_file(value) if getattr(self, 'state', None) else str(value or '')

    def _sync_directory_state_only(self):
        if getattr(self, 'state', None) is not None:
            self.state.working_dir = str(self.dirBox.GetValue() or '').strip()
            self.state.raw_path = str(self.outPathBox.GetValue() or './raw').strip() or './raw'
            self.state.spec_path = str(self.specPathBox.GetValue() or './spec').strip() or './spec'


    def _active_deconvolution_path(self, spectrum=None):
        """Return the deconvolution product appropriate to the active topology.

        Normal spectra own ``<spectrum>.decon``.  Pseudo2D is deliberately
        different: SpinUniDec analyses the materialised 1D spectral projection,
        so the authoritative calculated spectrum is ``<projection>.decon``.
        Keep this distinction at the GUI boundary rather than copying/renaming
        the calculated projection as though it had the shape of the pseudo2D
        source data.
        """
        topology = self._active_topology()
        if topology.spectral_dim_count == 1 and topology.has_pseudo_axis:
            try:
                projection = self.get_pseudo2d_projection_data(ensure_file=False)
            except Exception:
                projection = None
            if projection and projection.get('path'):
                return str(projection['path']) + '.decon'
        spectrum = spectrum or self._project_path(self.infileBox.GetValue())
        return spectrum + '.decon' if spectrum else ''

    def update_project_lamps(self):
        """Refresh file-existence lamps; safe to call from any workflow step."""
        if not hasattr(self, 'projectLamps'):
            return
        spectrum = self._project_path(self.infileBox.GetValue())
        # Lamps are status decoration and must remain safe while a project is
        # incomplete (for example during setup or recovery of a legacy file).
        try:
            decon_path = self._active_deconvolution_path(spectrum)
        except ValueError:
            decon_path = ''
        paths = {
            'spectrum': spectrum,
            'decon': decon_path,
            'reference': self._project_path(self.referencePeakBox.GetValue()),
            'full': self._project_path(self.fullPeakBox.GetValue()),
        }
        for key, path in paths.items():
            exists = bool(path and os.path.isfile(path))
            lamp = self.projectLamps.get(key)
            if lamp is not None:
                lamp.SetBackgroundColour(wx.Colour(46, 160, 67) if exists else wx.Colour(210, 55, 55))
                lamp.SetToolTip(path if path else 'No file selected')
                lamp.Refresh()

        # Viewer buttons are actions on in-memory peak lists, not merely on a
        # selected filename.  Keep them grey until the corresponding list has
        # actually been loaded into the shared data store.
        store = getattr(self, 'store', None)
        peak_lists = getattr(store, 'peak_lists', {}) if store is not None else {}
        if hasattr(self, 'buttonReferencePeakList'):
            self.buttonReferencePeakList.Enable(bool(peak_lists.get('reference')))
        if hasattr(self, 'buttonFullPeakList'):
            self.buttonFullPeakList.Enable(bool(peak_lists.get('full')))

    def spectrum_box(self):
        self.spectrumLbl = wx.StaticBox(self, -1, 'Spectrum:')

        self.buttonProcess = wx.Button(self.spectrumLbl, label="Process", size=(90, 22))
        self.buttonRead = wx.Button(self.spectrumLbl, label="Load")
        self.buttonReOrganise = wx.Button(self.spectrumLbl, label="ReOrganise", size=(90, 22))
        self.reorganiseFrame = None

        # The NMR tab contains actions only; dataset topology is configured in
        # Workflow and stored canonically in ProjectState.
        self.spectrumButtonRow = wx.BoxSizer(wx.HORIZONTAL)
        self.spectrumButtonRow.Add(self.buttonProcess, 1, wx.RIGHT, 4)
        self.spectrumButtonRow.Add(self.buttonRead, 1, wx.RIGHT, 4)
        self.spectrumButtonRow.Add(self.buttonReOrganise, 1)

        self.spectrumSizer = wx.StaticBoxSizer(self.spectrumLbl, wx.VERTICAL)
        self.spectrumSizer.AddSpacer(6)
        self.spectrumSizer.Add(self.spectrumButtonRow, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)
        self.spectrumSizer.AddSpacer(8)

    def pre_read_disabling(self):

        ## 2d Peaklist disable
        self.buttonReadPeak.Enable(False)
        self.openPeakFileBtn.Enable(False)
        self.buttonPeaky.Enable(False)
        self.peakBox.Enable(False)
        self.peakLbl.Enable(False)
        self.peakLab.Enable(False)
        self.fullPeakLab.Enable(False)
        self.fullPeakBox.Enable(False)
        self.openFullPeakFileBtn.Enable(False)
        self.buttonReadFullPeak.Enable(False)
        self.buttonReferencePeakList.Enable(False)
        self.buttonFullPeakList.Enable(False)

        ## Speed buttons
        self.quickButton.Disable()
        self.accurateButton.Disable()
        self.mediumButton.Disable()
        self.speedList.Disable()

        ## Deconvolution Disable
        self.deconLbl.Enable(False)

        self.buttonDecon.Enable(False)
        self.buttonPeakFit.Enable(False)
        self.buttonAnalyse.Enable(False)

        self.coreBox.Enable(False)
        self.threshBox.Enable(False)
        self.facBox.Enable(False)

        self.convBox.Enable(False)
        self.maxiterBox.Enable(False)
        self.initlistBox.Enable(False)
        self.fitRadBox.Enable(False)
        self.fitF1Box.Enable(False)
        self.fitF2Box.Enable(False)

        self.openInitListBtn.Enable(False)
        self.coreLab.Enable(False)
        self.threshLab.Enable(False)
        self.facLab.Enable(False)
        self.convlab.Enable(False)
        self.maxiterLab.Enable(False)
        self.initlistLab.Enable(False)
        self.fitRadLab.Enable(False)
        self.fitF1Lab.Enable(False)
        self.fitF2Lab.Enable(False)

        self.cb_grid_label.Enable(False)
        self.cb_decon3d_label.Enable(False)
        self.cb_decback_label.Enable(False)
        self.cb_fitphases_label.Enable(False)
        self.cb_enhance_label.Enable(False)
        self.cb_initlist_label.Enable(False)
        self.cb_grid.Enable(False)
        self.cb_decon3d.Enable(False)
        self.cb_decback.Enable(False)
        self.cb_fitphases.Enable(False)
        self.cb_enhance.Enable(False)
        self.cb_initlist.Enable(False)

        self.nLab.Enable(False)
        self.sigText.Enable(False)
        self.voigtText.Enable(False)
        self.lorText.Enable(False)

        if self.dim >=1:
            self.sig1Box.Enable(False)
            self.voigt1Box.Enable(False)
            self.lorentz1Box.Enable(False)
            self.n1Lab.Enable(False)
        if self.dim >=2:
            self.sig2Box.Enable(False)
            self.voigt2Box.Enable(False)
            self.lorentz2Box.Enable(False)
            self.n2Lab.Enable(False)
        if self.dim >=3:
            self.sig3Box.Enable(False)
            self.voigt3Box.Enable(False)
            self.lorentz3Box.Enable(False)
            self.n3Lab.Enable(False)
        if self.dim >=4:
            self.sig4Box.Enable(False)
            self.voigt4Box.Enable(False)
            self.lorentz4Box.Enable(False)
            self.n4Lab.Enable(False)

    def pre_read_enabling(self):

        ## 2d Peaklist disable
        self.buttonReadPeak.Enable(True)
        self.openPeakFileBtn.Enable(True)
        self.buttonPeaky.Enable(True)
        self.peakBox.Enable(True)
        self.peakLbl.Enable(True)
        self.peakLab.Enable(True)
        self.fullPeakLab.Enable(True)
        self.fullPeakBox.Enable(True)
        self.openFullPeakFileBtn.Enable(True)
        self.buttonReadFullPeak.Enable(True)
        self.buttonReferencePeakList.Enable(False)
        self.buttonFullPeakList.Enable(False)
        self._update_full_peak_controls()

        ## Speed buttons
        self.quickButton.Enable()
        self.accurateButton.Enable()
        self.mediumButton.Enable()
        self.speedList.Enable()

        ## Deconvolution Disable
        self.deconLbl.Enable(True)

        self.buttonDecon.Enable(True)
        self.buttonPeakFit.Enable(True)
        self.buttonAnalyse.Enable(True)

        self.coreBox.Enable(True)
        self.threshBox.Enable(True)
        self.facBox.Enable(True)

        self.convBox.Enable(True)
        self.maxiterBox.Enable(True)
        self.initlistBox.Enable(True)
        self.fitRadBox.Enable(True)
        self.fitF1Box.Enable(True)
        self.fitF2Box.Enable(True)

        self.openInitListBtn.Enable(True)
        self.coreLab.Enable(True)
        self.threshLab.Enable(True)
        self.facLab.Enable(True)
        self.convlab.Enable(True)
        self.maxiterLab.Enable(True)
        self.initlistLab.Enable(True)
        self.fitRadLab.Enable(True)
        self.fitF1Lab.Enable(True)
        self.fitF2Lab.Enable(True)

        self.cb_grid_label.Enable(True)
        self.cb_decon3d_label.Enable(True)
        self.cb_decback_label.Enable(True)
        self.cb_fitphases_label.Enable(True)
        self.cb_enhance_label.Enable(True)
        self.cb_initlist_label.Enable(True)
        self.cb_grid.Enable(True)
        self.cb_decon3d.Enable(True)
        self.cb_decback.Enable(True)
        self.cb_fitphases.Enable(True)
        self.cb_enhance.Enable(True)
        self.cb_initlist.Enable(True)

        self.nLab.Enable(True)
        self.sigText.Enable(True)
        self.voigtText.Enable(True)
        self.lorText.Enable(True)

        if self.dim >=1:
            self.sig1Box.Enable(True)
            self.voigt1Box.Enable(True)
            self.lorentz1Box.Enable(True)
            self.n1Lab.Enable(True)
        if self.dim >=2:
            self.sig2Box.Enable(True)
            self.voigt2Box.Enable(True)
            self.lorentz2Box.Enable(True)
            self.n2Lab.Enable(True)
        if self.dim >=3:
            self.sig3Box.Enable(True)
            self.voigt3Box.Enable(True)
            self.lorentz3Box.Enable(True)
            self.n3Lab.Enable(True)
        if self.dim >=4:
            self.sig4Box.Enable(True)
            self.voigt4Box.Enable(True)
            self.lorentz4Box.Enable(True)
            self.n4Lab.Enable(True)

    def peaklist_box(self):
        """Build independent controls for the 2D reference and full nD lists."""
        self.peakLbl = wx.StaticBox(self, -1, 'Peak Lists:')

        # Reference list: always 2D and used for projection/slice navigation and
        # for bore-mode deconvolution when "Use 2D peaklist" is selected.
        self.referencePeakLab = wx.StaticText(self.peakLbl, label="Reference 2D:", size=(83,-1))
        self.referencePeakBox = wx.TextCtrl(self.peakLbl, size=(200, 22))
        self.openPeakFileBtn = wx.Button(self.peakLbl, label="...", size=(40,22))
        self.openPeakFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetSpecFile(evt, self.referencePeakBox))
        self.buttonReadPeak = wx.Button(self.peakLbl, label="Load")
        self.buttonReferencePeakList = wx.Button(self.peakLbl, label="Show", size=(-1,22))
        self.buttonPeaky = wx.Button(self.peakLbl, label="Show", size=(-1,22))

        # Full list: dimensionality follows the main spectrum (2D/3D/4D...).
        self.fullPeakLab = wx.StaticText(self.peakLbl, label="Full nD:", size=(83,-1))
        self.fullPeakBox = wx.TextCtrl(self.peakLbl, size=(200, 22))
        self.openFullPeakFileBtn = wx.Button(self.peakLbl, label="...", size=(40,22))
        self.openFullPeakFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetSpecFile(evt, self.fullPeakBox))
        self.buttonReadFullPeak = wx.Button(self.peakLbl, label="Load")
        self.buttonFullPeakList = wx.Button(self.peakLbl, label="Show", size=(-1,22))
        self.buttonRestricted3DDiagnostics = wx.Button(self.peakLbl, label="Diagnostics", size=(-1,22))

        # Compatibility aliases during the migration.  Existing parameter-file
        # and bore-mode code that refers to peakBox continues to mean the
        # reference 2D list, never the full list.
        self.peakLab = self.referencePeakLab
        self.peakBox = self.referencePeakBox

        cnt=0
        self.sizerPk = wx.GridBagSizer(4, 4)
        self.sizerPk.Add(10,0, (cnt,0))
        self.sizerPk.Add(0,10, (cnt,1));cnt+=1
        self.sizerPk.Add(self.referencePeakLab, (cnt, 1), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.referencePeakBox, (cnt, 2), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.openPeakFileBtn, (cnt, 3), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.buttonPeaky, (cnt, 4), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP);cnt+=1
        self.sizerPk.Add(self.buttonReadPeak,(cnt,1));cnt+=1

        self.sizerPk.Add(0,5, (cnt,0));cnt+=1
        self.sizerPk.Add(self.fullPeakLab, (cnt, 1), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.fullPeakBox, (cnt, 2), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.openFullPeakFileBtn, (cnt, 3), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT| wx.TOP)
        self.sizerPk.Add(self.buttonRestricted3DDiagnostics, (cnt, 4), border=1, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP);cnt+=1
        self.sizerPk.Add(self.buttonReadFullPeak,(cnt,1));cnt+=1
        self.sizerPk.Add(0,5, (cnt,0));cnt+=1

        self.peakSizer = wx.StaticBoxSizer(self.peakLbl, wx.VERTICAL)
        self.borderPk = wx.BoxSizer()
        self.borderPk.Add(self.sizerPk, 1, wx.EXPAND)
        self.peakSizer.Add(self.borderPk, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)
        self.peakSizer.AddSpacer(5)

    def decon_box(self, gui_parent=None):
        """Build the UniDec controls as four responsive quadrants.

        The page is deliberately arranged as a 2x2 grid: Convergence and Peak
        shape on the top row, Options and Controls on the bottom row.  The
        Peak shape quadrant owns the dynamic ``sizer3`` grid rebuilt by
        ``DoDim`` so changing dimensionality can expand/contract naturally.
        Historical widget attributes and event bindings are preserved.
        """
        if gui_parent is None:
            gui_parent = self

        self.deconLbl = wx.StaticBox(gui_parent, -1, '')
        self.deconSizer = wx.StaticBoxSizer(self.deconLbl, wx.VERTICAL)

        # wxPython 4.2+ requires every window managed by a wxStaticBoxSizer
        # to be a direct child of that sizer's wxStaticBox.
        self.convergenceBox = wx.StaticBox(self.deconLbl, -1, 'Convergence')
        self.peakShapeBox = wx.StaticBox(self.deconLbl, -1, 'Peak shape')
        self.optionsBox = wx.StaticBox(self.deconLbl, -1, 'Options:')
        self.controlsBox = wx.StaticBox(self.deconLbl, -1, 'Controls:')

        self.buttonDecon = wx.Button(self.controlsBox, label="Decon", size=(-1, 22))
        self.buttonRecon = wx.Button(self.controlsBox, label="Recon", size=(-1, 22))
        self.buttonPeakFit = wx.Button(self.controlsBox, label="Fit Peaks", size=(-1, 22))
        self.buttonAnalyse = wx.Button(self.controlsBox, label="Load")

        self.coreBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.facBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.convBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.maxiterBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.initlistBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.fitRadBox = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        # Optional fixed 2D fitting radii.  These share the persistent keys
        # used by the pseudo-3D Fitting window.  Blank means use automatic
        # FitRad-derived radii.
        self.fitF1Box = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.fitF2Box = wx.TextCtrl(self.convergenceBox, size=(50, 22))
        self.openInitListBtn = wx.Button(self.convergenceBox, label="...", size=(40, 22))

        self.coreLab = wx.StaticText(self.convergenceBox, label="CPUs:")
        self.facLab = wx.StaticText(self.convergenceBox, label="Factor:")
        self.convlab = wx.StaticText(self.convergenceBox, label="Conv:")
        self.maxiterLab = wx.StaticText(self.convergenceBox, label="MaxIter:")
        self.initlistLab = wx.StaticText(self.convergenceBox, label="InitList:")
        self.fitRadLab = wx.StaticText(self.convergenceBox, label="FitRad:")
        self.fitF1Lab = wx.StaticText(self.convergenceBox, label="F1 radius:")
        self.fitF2Lab = wx.StaticText(self.convergenceBox, label="F2 radius:")
        self.openInitListBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.initlistBox))

        # Top-left: convergence/computation settings.
        self.convergenceSizer = wx.StaticBoxSizer(self.convergenceBox, wx.VERTICAL)
        self.sizer2 = wx.GridBagSizer(6, 7)
        rows = (
            (self.coreLab, self.coreBox, self.maxiterLab, self.maxiterBox),
            (self.facLab, self.facBox, self.initlistLab, self.initlistBox),
            (self.convlab, self.convBox, self.fitRadLab, self.fitRadBox),
        )
        for row, (lab1, box1, lab2, box2) in enumerate(rows):
            self.sizer2.Add(lab1, (row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(box1, (row, 1), flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
            self.sizer2.Add(lab2, (row, 2), flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=10)
            self.sizer2.Add(box2, (row, 3), flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.Add(self.openInitListBtn, (1, 4), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.AddGrowableCol(1, 1)
        self.sizer2.AddGrowableCol(3, 1)

        self.speedList = wx.StaticText(self.convergenceBox, label="Speed:")
        self.quickButton = wx.Button(self.convergenceBox, label="Quick", size=(-1, 22))
        self.mediumButton = wx.Button(self.convergenceBox, label="Medium", size=(-1, 22))
        self.accurateButton = wx.Button(self.convergenceBox, label="Accurate", size=(-1, 22))
        self.speedButtons = wx.BoxSizer(wx.HORIZONTAL)
        self.speedButtons.Add(self.speedList, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        self.speedButtons.Add(self.quickButton, 1, wx.RIGHT | wx.EXPAND, 4)
        self.speedButtons.Add(self.mediumButton, 1, wx.RIGHT | wx.EXPAND, 4)
        self.speedButtons.Add(self.accurateButton, 1, wx.EXPAND)
        self.quickButton.Bind(wx.EVT_BUTTON, self.on_quick_button)
        self.mediumButton.Bind(wx.EVT_BUTTON, self.on_medium_button)
        self.accurateButton.Bind(wx.EVT_BUTTON, self.on_accurate_button)
        self.convergenceSizer.Add(self.sizer2, 0, wx.EXPAND | wx.ALL, 8)
        self.convergenceSizer.Add(self.speedButtons, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        # Explicit F1/F2 FIT extraction radii for 2D and 2D+pseudo-axis data.
        self.fitFixedSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fitFixedSizer.Add(self.fitF1Lab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.fitFixedSizer.Add(self.fitF1Box, 1, wx.RIGHT, 12)
        self.fitFixedSizer.Add(self.fitF2Lab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.fitFixedSizer.Add(self.fitF2Box, 1)
        self.convergenceSizer.Add(self.fitFixedSizer, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        # Bottom-left: option checkboxes.
        self.cb_grid_label = wx.StaticText(self.optionsBox, label="Impose symmetry? (eg NOE)")
        self.cb_decon3d_label = wx.StaticText(self.optionsBox, label="Use 2D peaklist")
        self.cb_decback_label = wx.StaticText(self.optionsBox, label="Fit?")
        self.cb_fitphases_label = wx.StaticText(self.optionsBox, label="Phases?")
        self.cb_enhance_label = wx.StaticText(self.optionsBox, label="Enhance?")
        self.cb_initlist_label = wx.StaticText(self.optionsBox, label="Initialise from list?")
        self.cb_grid = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.cb_decon3d = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.cb_decback = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.cb_fitphases = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.cb_fitphases.SetToolTip("Pseudo2D only: fit the shared absorptive + dispersive distortion model during restrained reconstruction.")
        self.cb_enhance = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.cb_enhance.SetToolTip("Write the unclustered single-pass UniDec source distribution (1D-3D only).")
        self.cb_initlist = wx.CheckBox(self.optionsBox, -1, "", style=wx.ALIGN_RIGHT)
        self.sizer4 = wx.GridBagSizer(6, 6)
        options = (
            (self.cb_decon3d, self.cb_decon3d_label),
            (self.cb_grid, self.cb_grid_label),
            (self.cb_initlist, self.cb_initlist_label),
        )
        for row, (checkbox, label) in enumerate(options):
            self.sizer4.Add(checkbox, (row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            self.sizer4.Add(label, (row, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        fit_row = len(options)
        self.sizer4.Add(self.cb_decback, (fit_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer4.Add(self.cb_decback_label, (fit_row, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer4.Add(self.cb_fitphases, (fit_row, 2), flag=wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=12)
        self.sizer4.Add(self.cb_fitphases_label, (fit_row, 3), flag=wx.ALIGN_CENTER_VERTICAL)
        enhance_row = fit_row + 1
        self.sizer4.Add(self.cb_enhance, (enhance_row, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer4.Add(self.cb_enhance_label, (enhance_row, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer4.AddGrowableCol(1, 1)
        self.optionsSizer = wx.StaticBoxSizer(self.optionsBox, wx.VERTICAL)
        self.optionsSizer.Add(self.sizer4, 1, wx.EXPAND | wx.ALL, 8)

        # Top-right: dynamic line-shape controls. DoDim preserves/rebuilds this
        # exact GridBagSizer and its row/column coordinates.
        self.sizer3 = wx.GridBagSizer(6, 8)
        width = 200
        self.nLab = wx.StaticText(self.peakShapeBox, label="Dim:")
        self.nLab.Wrap(width)
        self.sigText = wx.StaticText(self.peakShapeBox, label="Gauss Width(ppm):")
        self.sigText.Wrap(width)
        self.voigtText = wx.StaticText(self.peakShapeBox, label="Voigt (0-1):")
        self.lorText = wx.StaticText(self.peakShapeBox, label="Lorentz Width(ppm):")
        self.sizer3.Add(self.nLab, (2, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer3.Add(self.sigText, (3, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer3.Add(self.voigtText, (4, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer3.Add(self.lorText, (5, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer3.AddGrowableCol(1, 2)
        self.sizer3.SetEmptyCellSize((0, 0))
        self.peakShapeSizer = wx.StaticBoxSizer(self.peakShapeBox, wx.VERTICAL)
        self.peakShapeSizer.Add(self.sizer3, 1, wx.EXPAND | wx.ALL, 8)

        # Compatibility attributes retained for older code paths.
        self.above = wx.BoxSizer(wx.HORIZONTAL)
        self.aboveLine = wx.Panel(self.deconLbl, -1, size=(1, 1)); self.aboveLine.Hide()
        self.below = wx.BoxSizer(wx.HORIZONTAL)
        self.belowLine = wx.Panel(self.deconLbl, -1, size=(1, 1)); self.belowLine.Hide()
        self.below_below = wx.BoxSizer(wx.HORIZONTAL)
        self.below_belowLine = wx.Panel(self.deconLbl, -1, size=(1, 1)); self.below_belowLine.Hide()

        # Bottom-right: action controls.
        self.controlsSizer = wx.StaticBoxSizer(self.controlsBox, wx.VERTICAL)
        self.buttonSizer = wx.GridBagSizer(6, 6)
        self.buttonSizer.Add(self.buttonDecon, (0, 0), flag=wx.EXPAND)
        self.buttonSizer.Add(self.buttonPeakFit, (0, 1), flag=wx.EXPAND)
        self.buttonSizer.Add(self.buttonAnalyse, (0, 2), flag=wx.EXPAND)
        self.buttonSizer.Add(self.buttonRecon, (1, 0), flag=wx.EXPAND)
        for col in range(3):
            self.buttonSizer.AddGrowableCol(col, 1)
        # Keep the action buttons flush with the top of the Controls box.
        self.controlsSizer.Add(self.buttonSizer, 0, wx.EXPAND | wx.ALL, 8)
        self.controlsSizer.AddStretchSpacer(1)

        # Equal, responsive 2x2 quadrants.  Both rows and columns grow, while
        # each child sizer expands to consume its quadrant.
        self.quadrantSizer = wx.GridSizer(rows=2, cols=2, vgap=8, hgap=8)
        self.quadrantSizer.Add(self.convergenceSizer, 1, wx.EXPAND)
        self.quadrantSizer.Add(self.peakShapeSizer, 1, wx.EXPAND)
        self.quadrantSizer.Add(self.optionsSizer, 1, wx.EXPAND)
        self.quadrantSizer.Add(self.controlsSizer, 1, wx.EXPAND)

        self.border2 = wx.BoxSizer(wx.VERTICAL)
        self.border2.Add(self.quadrantSizer, 1, wx.EXPAND | wx.ALL, 10)
        self.deconSizer.Add(self.border2, 1, wx.EXPAND)

    def status_box(self):
        self.statusLbl = wx.StaticBox(self, -1, 'Report:')
        self.statusSizer = wx.StaticBoxSizer(self.statusLbl, wx.VERTICAL)
        self.sizerStat = wx.GridBagSizer(12, 2)
        self.sizerStat.Add((12,0),(0,0))
        self.updateList=[]
        for i in range((15)):
            self.updateList.append(wx.StaticText(self.statusLbl, label="           ",size=(-1,-1)))
            self.sizerStat.Add(self.updateList[-1],(i+1,1))

        self.border4 = wx.BoxSizer(wx.VERTICAL)
        self.border4.Add(self.sizerStat,100, flag=wx.GROW)
        self.statusSizer.Add(self.border4,100, flag=wx.GROW)
        
        
        self.border4.AddSpacer(20)

    def projection_box(self):
        self.projectionLbl = wx.StaticBox(self, -1, 'Projection:')
        self.projectionSizer = wx.StaticBoxSizer(self.projectionLbl, wx.VERTICAL)
        self.fig = Figure()
        self.canvas = FigCanvas(self.projectionLbl, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(10,67))
        self.axes = self.fig.add_subplot(111)
        self.axes.get_yaxis().set_visible(False)
        self.axes.get_xaxis().set_visible(False)
        self.axes.set_frame_on(False)
        # self.axis_title = wx.StaticText(self, label="Projection:",size=(-1,-1))

        self.plotter = wx.BoxSizer(wx.HORIZONTAL)
        # self.plotter.AddSpacer(20)
        # self.plotter.Add(self.axis_title, 0, border=10, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.BOTTOM)
        self.plotter.Add(self.canvas, 100, flag=wx.GROW)
        # self.plotter.AddSpacer(20)
        self.projectionSizer.Add(self.plotter, 1, flag=wx.GROW)

    def noise_box(self):
        # Reserved display area for noise analysis.  Keep this independent of
        # the Projection figure so behaviour can be added without disturbing
        # the existing projection canvas/axes API.
        self.noiseLbl = wx.StaticBox(self, -1, 'Noise:')
        self.noiseSizer = wx.StaticBoxSizer(self.noiseLbl, wx.VERTICAL)
        self.noiseFig = Figure()
        self.noiseCanvas = FigCanvas(self.noiseLbl, -1, self.noiseFig)
        self.noiseCanvas.SetMinSize(wx.Size(10, 200))
        self.noiseAxes = self.noiseFig.add_subplot(111)
        self.noiseAxes.get_yaxis().set_visible(False)
        self.noiseAxes.get_xaxis().set_visible(False)
        self.noiseAxes.set_frame_on(False)
        self.noisePlotter = wx.BoxSizer(wx.HORIZONTAL)
        self.noisePlotter.Add(self.noiseCanvas, 100, flag=wx.GROW)
        self.noiseSizer.Add(self.noisePlotter, 1, flag=wx.GROW)

        # Threshold belongs with the noise display because it directly controls
        # the signal cut-off shown on this plot.  Keep the historical attribute
        # names so parameter save/load and the existing Set handler remain intact.
        self.threshLab = wx.StaticText(self.noiseLbl, label='Threshold:')
        self.threshBox = wx.TextCtrl(self.noiseLbl, size=(50, 22), style=wx.TE_PROCESS_ENTER)
        # Commit threshold changes only when ENTER is pressed.  This avoids
        # repeatedly updating the rest of the GUI while a value is still being typed.
        self.threshBox.Bind(wx.EVT_TEXT_ENTER, self._on_threshold_text_changed)

        self.noiseDetailButton = wx.Button(self.noiseLbl, label='Detail')
        self.noiseDetailButton.SetToolTip('Show detailed noise and intensity statistics')
        self.noiseDetailButton.Bind(wx.EVT_BUTTON, self.on_noise_detail)

        self.noiseControlRow = wx.BoxSizer(wx.HORIZONTAL)
        self.noiseControlRow.AddStretchSpacer(1)
        self.noiseControlRow.Add(self.threshLab, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.noiseControlRow.Add(self.threshBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.noiseControlRow.Add(self.noiseDetailButton, 0, wx.ALIGN_CENTER_VERTICAL)
        self.noiseControlRow.AddStretchSpacer(1)
        self.noiseSizer.Add(self.noiseControlRow, 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 3)

        # Artists used by the lightweight threshold redraw.  The background is
        # captured after a complete noise-plot draw; subsequent threshold edits
        # restore it and redraw only the threshold line/label (Matplotlib blit).
        self._noiseThresholdLine = None
        self._noiseThresholdLabel = None
        self._noiseBlitBackground = None
        self._noiseBlitReady = False
        self.noiseCanvas.mpl_connect('resize_event', self._invalidate_noise_blit)


    def _invalidate_noise_blit(self, event=None):
        self._noiseBlitBackground = None
        self._noiseBlitReady = False

    def _threshold_plot_x(self):
        """Return threshold location in the noise plot's |z| coordinates."""
        stats = getattr(self, 'noiseStats', None)
        if not stats:
            return None
        try:
            fraction = float(self.threshBox.GetValue())
            max_signal = float(stats['max_intensity'])
            mu = float(stats['centre'])
            sigma = float(stats['noise_sigma'])
            if not numpy.isfinite(fraction) or not numpy.isfinite(sigma) or sigma <= 0:
                return None
            # The requested threshold is an absolute spectrum intensity.  The
            # histogram x-axis is |intensity-centre|/sigma, so transform it to
            # the displayed coordinate before drawing the vertical marker.
            threshold_intensity = max_signal * fraction
            return max(abs((threshold_intensity - mu) / sigma), numpy.finfo(float).tiny)
        except (TypeError, ValueError, KeyError, AttributeError):
            return None

    def _create_noise_threshold_artists(self):
        x = self._threshold_plot_x()
        if x is None:
            return
        ax = self.noiseAxes
        self._noiseThresholdLine = ax.axvline(x, color='green', linewidth=1.2,
                                               label='threshold', animated=True)
        self._noiseThresholdLabel = ax.annotate(
            'threshold', xy=(x, 0.98), xycoords=('data', 'axes fraction'),
            xytext=(3, 0), textcoords='offset points', rotation=90,
            va='top', ha='left', fontsize=7, color='green', animated=True)

    def _capture_noise_blit_background(self):
        if self._noiseThresholdLine is None or self._noiseThresholdLabel is None:
            return
        try:
            self._noiseBlitBackground = self.noiseCanvas.copy_from_bbox(self.noiseAxes.bbox)
            self._noiseBlitReady = True
            self._blit_noise_threshold()
        except Exception:
            self._noiseBlitBackground = None
            self._noiseBlitReady = False

    def _blit_noise_threshold(self):
        """Move only the threshold artists, leaving the noise plot untouched."""
        x = self._threshold_plot_x()
        if x is None or self._noiseThresholdLine is None or self._noiseThresholdLabel is None:
            return
        self._noiseThresholdLine.set_xdata([x, x])
        self._noiseThresholdLabel.xy = (x, 0.98)
        if not self._noiseBlitReady or self._noiseBlitBackground is None:
            return
        try:
            self.noiseCanvas.restore_region(self._noiseBlitBackground)
            self.noiseAxes.draw_artist(self._noiseThresholdLine)
            self.noiseAxes.draw_artist(self._noiseThresholdLabel)
            self.noiseCanvas.blit(self.noiseAxes.bbox)
        except Exception:
            self._invalidate_noise_blit()

    def _on_threshold_text_changed(self, event=None):
        # Pressing ENTER in Threshold performs the same update that the former Set
        # button performed.  Ordinary text edits deliberately do nothing until
        # the value is committed with ENTER.
        try:
            float(self.threshBox.GetValue())
        except (TypeError, ValueError):
            if event is not None:
                event.Skip()
            return

        self._blit_noise_threshold()
        self.OnButtonNoise(None)
        if event is not None:
            event.Skip()

    def _noise_detail_text(self):
        """Return a readable explanation of the current noise statistics."""
        stats = getattr(self, 'noiseStats', None)
        if not stats:
            return ('No noise statistics are available yet.\n\n'
                    'Read a spectrum first; the noise analysis is performed automatically.')

        tails = stats.get('tail_counts', {})
        def tail_line(level):
            t = tails.get(level, {})
            pos = int(t.get('positive', 0))
            neg = int(t.get('negative', 0))
            return (f'  Beyond {level} sigma:  positive = {pos:,}, '
                    f'negative = {neg:,}, positive excess = {pos-neg:,}')

        sampled = int(stats.get('sampled_points', 0))
        total = int(stats.get('points', 0))
        centre = stats.get('centre', float('nan'))
        sigma = stats.get('noise_sigma', float('nan'))
        mad_sigma = stats.get('noise_mad_sigma', float('nan'))
        max_intensity = stats.get('max_intensity', float('nan'))
        max_snr = stats.get('max_snr', float('nan'))
        core_fraction = 100.0 * stats.get('core_fraction', 0.0)

        return (
            'NMR NOISE / INTENSITY STATISTICS\n'
            '================================\n\n'
            f'Spectral points: {total:,}\n'
            f'Points used for analysis: {sampled:,}\n\n'
            f'Gaussian centre (mu): {centre:.6g}\n'
            '  The fitted centre of the noise-dominated intensity distribution. '
            'For a well baseline-corrected spectrum this should normally be close to zero.\n\n'
            f'Fitted noise sigma: {sigma:.6g}\n'
            '  The estimated standard deviation of the Gaussian noise. The fit is made '
            'iteratively to the central 2.5-sigma region so that sparse NMR resonances '
            'have much less influence on the estimate. This value is also used as the '
            'absolute noise estimate for signal-to-noise calculations.\n\n'
            f'MAD noise sigma: {mad_sigma:.6g}\n'
            '  A robust independent estimate based on 1.4826 x median absolute deviation. '
            'Agreement between MAD sigma and fitted sigma is a useful indication that the '
            'central distribution is behaving approximately like Gaussian noise.\n\n'
            f'Central fit population: {core_fraction:.2f}% of sampled points\n'
            '  Fraction retained in the final clipped Gaussian core. Resonance-rich, broad, '
            'or strongly non-Gaussian spectra may give a smaller value.\n\n'
            f'Maximum sampled intensity: {max_intensity:.6g}\n'
            f'Maximum positive S/N: {max_snr:.2f}\n'
            '  Maximum positive intensity above the fitted centre, expressed in units of '
            'the fitted noise sigma. It is a point-wise maximum, not a peak-integrated S/N.\n\n'
            'TAIL POPULATIONS\n'
            + '\n'.join(tail_line(level) for level in (2, 3, 4, 5)) + '\n\n'
            '  For symmetric Gaussian noise, positive and negative tail counts should be '
            'similar. An excess of positive points is expected in a conventional, correctly '
            'phased absorption spectrum because real resonances contribute predominantly '
            'positive intensity. Large excesses in both directions can instead indicate '
            'non-Gaussian noise, phase/baseline problems, truncation artefacts, or genuine '
            'signals of both signs.\n\n'
            'HOW TO READ THE HISTOGRAM\n'
            '  The x-axis is absolute intensity relative to the fitted noise sigma. Positive '
            'and negative spectral points are counted separately. Both axes are logarithmic, '
            'so the dense Gaussian noise core and rare high-intensity points can be viewed '
            'together. The dashed Gaussian curve is the expected one-sided population for '
            'ideal Gaussian noise. Systematic departure from that curve, particularly a '
            'positive-tail excess, is evidence for signal or other non-noise structure.\n\n'
            '  These statistics describe the distribution of individual spectral points. '
            'They do not use neighbouring-point correlations, peak shapes, multiplet '
            'structure, or chemical shift, so they are intentionally a simple first-order '
            'description of spectrum sparsity, noise and intensity information.'
        )

    def on_noise_detail(self, event=None):
        """Open a scrollable window explaining the current noise statistics."""
        dlg = wx.Dialog(self, title='Noise Analysis Detail', size=(680, 650),
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        outer = wx.BoxSizer(wx.VERTICAL)
        text = wx.TextCtrl(dlg, value=self._noise_detail_text(),
                           style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH2)
        outer.Add(text, 1, wx.EXPAND | wx.ALL, 10)
        close = wx.Button(dlg, wx.ID_OK, 'Close')
        outer.Add(close, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 10)
        dlg.SetSizer(outer)
        dlg.CentreOnParent()
        dlg.ShowModal()
        dlg.Destroy()


    def analyse_noise_spectrum(self, spectrum=None, max_points=2000000):
        """Analyse and plot the intensity statistics of an NMR spectrum.

        ``spectrum`` may be a numpy-like data array or an object exposing a
        ``data`` attribute.  The analysis is deliberately model-light: a
        robust MAD estimate seeds an iterative Gaussian fit to the central
        (noise-dominated) intensity distribution.  Positive and negative
        absolute intensities are then shown separately on log/log axes.

        Returns a dictionary of useful scalar statistics, making this method
        usable independently of the GUI drawing code.
        """
        source = self.data if spectrum is None else getattr(spectrum, 'data', spectrum)
        try:
            values = numpy.asarray(source, dtype=float).ravel()
        except Exception:
            values = numpy.asarray(source[:], dtype=float).ravel()
        values = values[numpy.isfinite(values)]
        n_total = int(values.size)
        if n_total < 20:
            raise ValueError('At least 20 finite spectral points are required for noise analysis.')

        # Very large multidimensional spectra need not be fully materialised in
        # the histogram.  Uniform striding preserves the global distribution
        # without introducing a random/reproducibility dependency.
        if values.size > max_points:
            step = int(numpy.ceil(values.size / float(max_points)))
            sample = values[::step]
        else:
            sample = values

        median = float(numpy.median(sample))
        mad = float(numpy.median(numpy.abs(sample - median)))
        sigma_mad = 1.482602218505602 * mad
        if not numpy.isfinite(sigma_mad) or sigma_mad <= 0:
            sigma_mad = float(numpy.std(sample))
        if not numpy.isfinite(sigma_mad) or sigma_mad <= 0:
            raise ValueError('Spectrum has zero/undefined intensity spread.')

        # Iterative central Gaussian fit.  The 2.5-sigma clipping strongly
        # suppresses sparse resonances while retaining enough of the Gaussian
        # core to estimate its location and absolute noise scale.
        mu = median
        sigma = sigma_mad
        core = sample
        for _ in range(8):
            mask = numpy.abs(sample - mu) <= 2.5 * sigma
            core = sample[mask]
            if core.size < 20:
                break
            new_mu = float(numpy.mean(core))
            # Correct the variance lost by symmetric 2.5-sigma truncation.
            # For a standard normal Var(X | |X|<2.5) ~= 0.911256.
            new_sigma = float(numpy.std(core)) / numpy.sqrt(0.9112563609)
            if not numpy.isfinite(new_sigma) or new_sigma <= 0:
                break
            if abs(new_sigma - sigma) <= 1e-5 * sigma:
                mu, sigma = new_mu, new_sigma
                break
            mu, sigma = new_mu, new_sigma

        z = (sample - mu) / sigma
        max_abs_z = float(numpy.max(numpy.abs(z)))
        positive = z[z > 0]
        negative = -z[z < 0]
        thresholds = (2, 3, 4, 5)
        tails = {}
        for threshold in thresholds:
            pos_n = int(numpy.count_nonzero(z >= threshold))
            neg_n = int(numpy.count_nonzero(z <= -threshold))
            tails[threshold] = {'positive': pos_n, 'negative': neg_n, 'excess': pos_n - neg_n}

        # Log-spaced bins in |z| expose both the Gaussian core and sparse tails.
        nonzero = numpy.abs(z[numpy.nonzero(z)])
        low = max(0.05, float(numpy.percentile(nonzero, 0.1))) if nonzero.size else 0.05
        high = max(6.0, min(max_abs_z * 1.05, 1.0e6))
        if high <= low:
            high = low * 10.0
        edges = numpy.logspace(numpy.log10(low), numpy.log10(high), 90)
        centres = numpy.sqrt(edges[:-1] * edges[1:])
        pos_counts, _ = numpy.histogram(positive, bins=edges)
        neg_counts, _ = numpy.histogram(negative, bins=edges)

        ax = self.noiseAxes
        ax.clear()
        ax.set_frame_on(True)
        ax.get_xaxis().set_visible(True)
        ax.get_yaxis().set_visible(True)
        # Distinct colours are intentional here: sign is a data category.
        ax.step(centres, pos_counts, where='mid', color='tab:red', linewidth=1.2, label='+')
        ax.step(centres, neg_counts, where='mid', color='tab:blue', linewidth=1.2, label='-')

        # Expected one-sided Gaussian bin counts, normalised to the number of
        # sampled points.  erf is supplied by the Python standard library.
        from math import erf, sqrt
        cdf = numpy.array([0.5 * (1.0 + erf(float(e) / sqrt(2.0))) for e in edges])
        expected = sample.size * numpy.diff(cdf)
        ax.plot(centres, expected, color='0.25', linestyle='--', linewidth=1.0, label='Gaussian')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'Absolute intensity / fitted noise $\sigma$', fontsize=7)
        ax.set_ylabel('Spectral points / bin', fontsize=7)
        ax.tick_params(axis='both', which='both', labelsize=6)
        ax.grid(True, which='both', alpha=0.18)
        # Keep the legend with the compact statistics block in the lower-left.
        # Its order follows the plot creation order: +, -, Gaussian.
        ax.legend(loc='lower left', bbox_to_anchor=(0.025, 0.205),
                  borderaxespad=0.0, fontsize=7, frameon=False)

        max_signal = float(numpy.max(sample))
        snr_max = (max_signal - mu) / sigma
        asym3 = tails[3]['excess']
        stats = {
            'points': n_total,
            'sampled_points': int(sample.size),
            'centre': mu,
            'noise_sigma': sigma,
            'noise_mad_sigma': sigma_mad,
            'max_intensity': max_signal,
            'max_snr': float(snr_max),
            'core_fraction': float(core.size) / float(sample.size),
            'tail_counts': tails,
            'positive_excess_3sigma': asym3,
        }
        self.noiseStats = stats
        self.noiseVal = sigma

        # Compact at-a-glance summary.  The complement of the clipped
        # Gaussian-core population is a deliberately rough measure of spectral
        # occupancy: it includes genuine signal as well as any non-Gaussian
        # artefacts, so the Detail window gives the fuller interpretation.
        signal_fraction = 100.0 * (1.0 - stats['core_fraction'])
        summary = (f"Noise: {sigma:.2e}\n"
                   f"Max S/N: {snr_max:.2f}\n"
                   f"Spectral points: {signal_fraction:.2f}%")
        ax.text(0.025, 0.035, summary, transform=ax.transAxes,
                ha='left', va='bottom', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor='0.55', alpha=0.88))

        # Keep very sparse tail bins visible on a common logarithmic scale.
        ax.set_ylim(bottom=1.0e-1)

        self.noiseFig.tight_layout(pad=1.0)
        self._invalidate_noise_blit()
        self._noiseThresholdLine = None
        self._noiseThresholdLabel = None
        self._create_noise_threshold_artists()
        self.noiseCanvas.draw()
        self._capture_noise_blit_background()
        return stats

    # Public alias with a GUI-oriented name for callers that receive a
    # spectrum object after import/read.
    def plot_noise_histogram(self, spectrum=None):
        return self.analyse_noise_spectrum(spectrum)

    @property
    def deconParFile(self):
        return self.state.deconParFile if getattr(self, "state", None) is not None else getattr(self, "_deconParFile", "")

    @deconParFile.setter
    def deconParFile(self, path):
        self._deconParFile = str(path)
        if getattr(self, "state", None) is not None:
            self.state.deconParFile = path

    ####################################
    #is nmrPipe installed?
    def check_nmrPipe(self):
        try:
            result = subprocess.Popen(["csh -c 'which nmrPipe'"], stdout=subprocess.PIPE, shell=True).communicate()[0]
            
        except Exception:
            
            return False 
        if 'Command not found' in str(result):
            print(result)
            return False

        return True

    _DATASTORE_FIELDS = {
        "dic", "data", "uc0", "uc1", "uc2", "uc3", "index0", "index1", "index2", "index3",
        "XX", "YY", "ZZ", "XX2", "YY2",
        "uc0min", "uc0max", "uc1min", "uc1max", "uc2min", "uc2max", "uc3min", "uc3max",
        "dmax", "noiseVal", "projectedData", "datadec", "dicdec", "peak",
        "spectrumfile", "STD_raw_path", "STD_std_path", "STD", "mixingTimes", "labb", "dim",
        "pkIdx", "pkSlice1D", "pkSlice1Ddec", "Grps", "noeTags", "READ", "DECON", "PEAK",
        "pseudo", "uSTA"
    }

    def __setattr__(self, name, value):
        if name in self._DATASTORE_FIELDS and "store" in self.__dict__:
            setattr(self.store, name, value)
            return
        super().__setattr__(name, value)

    def __getattr__(self, name):
        if name in self._DATASTORE_FIELDS and "store" in self.__dict__:
            return getattr(self.store, name)
        raise AttributeError(f"{type(self).__name__!s} object has no attribute {name!r}")

    def _reset_store(self):
        if getattr(self, "store", None) is not None:
            pass
            self.store.reset()
            pass

    def _spectral_physical_axes(self):
        """Return ``(physical_index, label)`` for spectral axes only.

        ``self.dim`` is the canonical spectral count.  Peak coordinates must
        never include the sampled real/pseudo axis, even though ``labb`` and
        the NumPy array describe physical axes.
        """
        spectral_count = int(getattr(self, 'dim', 0) or 0)
        labels = list(getattr(self, 'labb', []) or [])
        pseudo = bool(self.state.pseudo_axis)
        if not pseudo:
            return [(i, str(labels[i]) if i < len(labels) else 'f%d' % (i + 1))
                    for i in range(spectral_count)]

        # For pseudo-dimensional data, physical-axis identity has already been
        # canonicalised at load time and is stored in ProjectState metadata.
        # Use that topology directly rather than trying to rediscover the pseudo
        # axis from a hard-coded label list.  Real-axis labels are user/data
        # dependent (for example ``usta``), so label heuristics can reject a
        # perfectly valid canonical pseudo2D dataset.
        topology = self._active_topology()
        axes = []
        for axis in topology.spectral_axes:
            physical_index = int(axis.physical_index)
            label = (str(labels[physical_index])
                     if physical_index < len(labels)
                     else str(axis.label or ('f%d' % (int(axis.spectral_index) + 1))))
            axes.append((physical_index, label))
        if len(axes) == spectral_count:
            return axes

        raise ValueError('Cannot identify spectral axes for pseudo-dimensional peak data.')

    def _spectral_axis_labels(self):
        return [label for _, label in self._spectral_physical_axes()]

    def _peak_list_suffix_for_dim(self, dim=None):
        dim = int(dim if dim is not None else getattr(self, 'dim', 0) or 0)
        return {1: '.1D.list', 2: '.2D.list', 3: '.3D.list', 4: '.4D.list'}.get(dim, f'.{dim}D.list')

    def _peak_list_path_for_spectrum(self, spectrum_path=None, dim=None):
        spectrum_path = str(spectrum_path or '').strip()
        if not spectrum_path:
            spectrum_path = self._resolve_input_path(self.infileBox.GetValue()) if hasattr(self, 'infileBox') else ''
        if not spectrum_path:
            return ''
        return spectrum_path + self._peak_list_suffix_for_dim(dim)

    """
    def _sync_peakfile_box(self, spectrum_path=None, dim=None):
        if hasattr(self, 'peakBox'):
            peak_path = self._peak_list_path_for_spectrum(spectrum_path, dim)
            if peak_path:
                self.peakBox.SetValue(peak_path)
        return getattr(self, 'peakBox', None) and self.peakBox.GetValue() or ''
    """
    

    def _update_full_peak_controls(self):
        """Keep peak-list controls aligned with the active scientific topology.

        Pseudo2D has only one spectral peak list.  A Reference 2D list has no
        meaning for this topology, so hide that row and expose the normal full
        list explicitly as ``Full 1D``.
        """
        dim = int(getattr(self, 'dim', 0) or self.state.spectral_dimensions)
        try:
            topology = self._active_topology()
            pseudo2d = (topology.spectral_dim_count == 1 and topology.has_pseudo_axis
                        and topology.physical_dim_count == 2)
            physical2d = (topology.spectral_dim_count == 2 and not topology.has_pseudo_axis)
            hide_reference = pseudo2d or physical2d
        except Exception:
            pseudo2d = physical2d = hide_reference = False
        if hasattr(self, 'fullPeakLab'):
            self.fullPeakLab.SetLabel("Full 1D:" if pseudo2d else (f"Full {dim}D:" if dim else "Full nD:"))
        # A true physical 2D dataset has no distinct Reference peak list, but
        # the historical Reference-row ``Show`` button is also the entry point
        # to PeakFrame ("Get 2D peaks").  Keep that one action visible while
        # hiding every control that exposes Reference as a separate list/file.
        # PeakFrame itself is topology-aware and edits the authoritative Full
        # 2D list for this case.
        for name in ('referencePeakLab', 'referencePeakBox', 'openPeakFileBtn',
                     'buttonReadPeak', 'buttonReferencePeakList'):
            control = getattr(self, name, None)
            if control is not None:
                control.Show(not hide_reference)
        if hasattr(self, 'buttonPeaky'):
            self.buttonPeaky.Show(physical2d or not hide_reference)
        project_row = getattr(self, 'projectRows', {}).get('reference')
        if project_row is not None:
            try:
                project_row.ShowItems(not hide_reference)
            except AttributeError:
                for item in project_row.GetChildren():
                    window = item.GetWindow()
                    if window is not None:
                        window.Show(not hide_reference)
        if hasattr(self, 'projectSizer'):
            self.projectSizer.Layout()
        if hasattr(self, 'peakSizer'):
            self.peakSizer.Layout()
        if hasattr(self, 'fullPeakBox') and not self.fullPeakBox.GetValue().strip():
            try:
                path = self._peak_list_path_for_spectrum(self._resolve_input_path(self.infileBox.GetValue()), dim)
            except Exception:
                path = ''
            if path:
                self.fullPeakBox.SetValue(self.state._spec_relative(path) if getattr(self, 'state', None) is not None else path)

    def reference_peak_save_destination(self, value):
        """Return (SpecPath-relative value, absolute destination) for a reference list."""
        self._sync_directory_state_only()
        raw = str(value or '').strip()
        if not raw:
            raise ValueError('Enter a file name before saving.')

        spec_dir = os.path.abspath(self.state.spec_dir())
        if os.path.isabs(raw):
            destination = os.path.abspath(raw)
        else:
            try:
                relative = self.state._spec_relative(raw)
            except ValueError:
                raise ValueError('Reference peak lists must be saved inside SpecPath.')
            destination = os.path.abspath(os.path.join(spec_dir, relative))

        try:
            if os.path.commonpath([spec_dir, destination]) != spec_dir:
                raise ValueError('Reference peak lists must be saved inside SpecPath.')
            relative = os.path.relpath(destination, spec_dir).replace(os.sep, '/')
        except (ValueError, OSError):
            raise ValueError('Reference peak lists must be saved inside SpecPath.')
        if relative == '..' or relative.startswith('../'):
            raise ValueError('Reference peak lists must be saved inside SpecPath.')
        return relative, destination

    def save_reference_peak_list(self, value):
        """Persist the current 2D reference peaks and update canonical GUI/project state."""
        relative, destination = self.reference_peak_save_destination(value)
        parent = os.path.dirname(destination)
        if parent and not os.path.isdir(parent):
            os.makedirs(parent)

        # Keep the legacy 2D peak-list format: name, Y, X.
        with open(destination, 'w') as outy:
            for pk in self.get_reference_peaks():
                outy.write('%s\t%f\t%f\n' % (pk.name, pk.y, pk.x))

        self.referencePeakBox.SetValue(relative)
        if getattr(self, 'state', None) is not None:
            self.state.reference_peak_file = relative
            self.state.dirty = True
        self.set_reference_peaks(self.get_reference_peaks(), source_path=destination)
        self.update_project_lamps()
        return destination

    def full_peak_save_destination(self, value):
        """Return (SpecPath-relative value, absolute destination) for a Full list."""
        self._sync_directory_state_only()
        raw = str(value or '').strip()
        if not raw:
            raise ValueError('Enter a file name before saving.')
        spec_dir = os.path.abspath(self.state.spec_dir())
        if os.path.isabs(raw):
            destination = os.path.abspath(raw)
        else:
            try:
                relative = self.state._spec_relative(raw)
            except ValueError:
                raise ValueError('Full peak lists must be saved inside SpecPath.')
            destination = os.path.abspath(os.path.join(spec_dir, relative))
        try:
            if os.path.commonpath([spec_dir, destination]) != spec_dir:
                raise ValueError('Full peak lists must be saved inside SpecPath.')
            relative = os.path.relpath(destination, spec_dir).replace(os.sep, '/')
        except (ValueError, OSError):
            raise ValueError('Full peak lists must be saved inside SpecPath.')
        if relative == '..' or relative.startswith('../'):
            raise ValueError('Full peak lists must be saved inside SpecPath.')
        return relative, destination

    def save_full_peak_list(self, value):
        """Persist the authoritative in-memory Full nD rows in loadable format."""
        relative, destination = self.full_peak_save_destination(value)
        payload = self.get_full_peak_payload()
        rows = list(payload.get('rows') or [])
        if not rows:
            raise ValueError('There is no Full peak list in memory to save.')
        parent = os.path.dirname(destination)
        if parent and not os.path.isdir(parent):
            os.makedirs(parent)
        with open(destination, 'w') as outy:
            for fields in rows:
                outy.write('\t'.join(str(value) for value in fields) + '\n')
        self.fullPeakBox.SetValue(relative)
        if getattr(self, 'state', None) is not None:
            self.state.full_peak_file = relative
            self.state.dirty = True
        payload['source_path'] = destination
        self.corrFile = destination
        self.update_project_lamps()
        return destination

    def get_reference_peaks(self):
        payload = self.store.peak_lists.get('reference', {}) if getattr(self, 'store', None) is not None else {}
        peaks = payload.get('peaks')
        return self.peak if peaks is None else peaks

    def _projection_peak_specs_3d(self):
        """Return canonical projection/display specifications for 3D data.

        ``lookup_axes`` and ``transpose`` describe how the corresponding
        spectrum projection view is retrieved from the projection cache.
        ``display_axes`` describe the *final plotted X/Y axes*.  Peak markers
        are always generated from ``display_axes`` and never transposed.  This
        keeps array-storage orientation separate from physical peak
        coordinates, which is especially important for the third 3D panel.
        """
        labels = list(getattr(self, 'labb', []) or [])
        if len(labels) < 3:
            return []
        return [
            {
                'lookup_axes': (str(labels[2]), str(labels[1])),
                'display_axes': (str(labels[2]), str(labels[1])),
                'transpose': 'n',
            },
            {
                'lookup_axes': (str(labels[2]), str(labels[0])),
                'display_axes': (str(labels[2]), str(labels[0])),
                'transpose': 'n',
            },
            {
                # The third spectrum projection is stored/read in the
                # (labels[1], labels[0]) order but displayed transposed, so
                # its final X/Y axes are (labels[0], labels[1]).
                'lookup_axes': (str(labels[1]), str(labels[0])),
                'display_axes': (str(labels[0]), str(labels[1])),
                'transpose': 'y',
            },
        ]

    def _projection_peak_spec_3d(self, left, right, transpose='n'):
        """Resolve a projection lookup request to its canonical 3D spec."""
        left, right = str(left), str(right)
        transpose = 'y' if str(transpose).lower() == 'y' else 'n'
        for spec in self._projection_peak_specs_3d():
            if spec['lookup_axes'] == (left, right) and spec['transpose'] == transpose:
                return spec
        return None

    def _full_peak_records(self, rows, dim=None):
        """Convert raw nD peak-list rows to dimension-independent records."""
        dim = int(dim if dim is not None else getattr(self, 'dim', 0) or 0)
        labels = list(getattr(self, 'labb', []) or [])
        records = []
        for row_index, fields in enumerate(rows or []):
            if len(fields) < dim + 1:
                continue
            try:
                coords = tuple(float(fields[i]) for i in range(1, dim + 1))
            except (TypeError, ValueError):
                continue
            # Decon peak-list coordinates are f1..fN (fastest to slowest),
            # while labb is indexed in the main-array dimension order.
            axis_values = {}
            spectral_labels = self._spectral_axis_labels()
            if len(spectral_labels) >= dim:
                for coord_index, value in enumerate(coords):
                    axis_values[str(spectral_labels[dim - 1 - coord_index])] = value
            intensity = None
            if len(fields) > dim + 1:
                try:
                    intensity = float(fields[dim + 1])
                except (TypeError, ValueError):
                    intensity = None
            records.append({
                'name': str(fields[0]),
                'coordinates': coords,
                'axis_values': axis_values,
                'intensity': intensity,
                'row_index': row_index,
                'fields': list(fields),
                'analysis': {},
            })
        return records

    def _reference_axis_values_3d(self, peak):
        labels = list(getattr(self, 'labb', []) or [])
        if len(labels) < 3:
            return {}
        # The reference list is always 2D.  In 3D bore mode its x coordinate
        # maps to dimension 2 and y maps to both dimensions 1 and 0.
        return {
            str(labels[2]): float(peak.x),
            str(labels[1]): float(peak.y),
            str(labels[0]): float(peak.y),
        }

    def _rebuild_projected_peak_lists(self):
        """Prepare GUI-ready peak overlays for all 3D projection panels.

        Projected markers are keyed by their final displayed X/Y axes.  The
        spectrum's transpose flag is deliberately *not* part of marker
        generation: transpose is an array-orientation concern only.
        """
        store = getattr(self, 'store', None)
        if store is None:
            return
        for key in list(store.projected_peak_lists):
            if isinstance(key, tuple) and key and key[0] in ('reference', 'full'):
                del store.projected_peak_lists[key]
        # A physical 2D spectrum has no distinct Reference peak list: Full 2D
        # is the sole peak authority.  The Fitting workspace is a single-plane
        # adapter and its ``Peaks`` toolbar action consumes this projected
        # overlay, so publish Full records under the fitting display axes.
        if int(getattr(self, 'dim', 0) or 0) == 2:
            pseudo_view = self.get_pseudo3d_view('raw')
            if pseudo_view is not None:
                display_x = str(pseudo_view['x_label'])
                display_y = str(pseudo_view['y_label'])
                full = self.get_full_peak_payload() or {}
                records = full.get('records') or full.get('peaks') or []
                points = []
                for index, record in enumerate(records):
                    if not isinstance(record, dict):
                        continue
                    axes = record.get('axis_values') or {}
                    try:
                        x = float(axes[display_x])
                        y = float(axes[display_y])
                    except (KeyError, TypeError, ValueError):
                        # Older/restored Full payloads may pre-date axis_values.
                        # Their canonical 2D coordinates are still sufficient.
                        coords = tuple(record.get('coordinates') or ())
                        if len(coords) < 2:
                            continue
                        try:
                            x, y = float(coords[0]), float(coords[1])
                        except (TypeError, ValueError):
                            continue
                    points.append({
                        'x': x,
                        'y': y,
                        'label': str(record.get('name', '')),
                        'source_index': record.get('row_index', index),
                        'source': 'full',
                        'axis_values': {display_x: x, display_y: y},
                    })
                store.projected_peak_lists[('full', display_x, display_y)] = {
                    'peaks': points,
                    'source': 'full',
                    'display_axes': (display_x, display_y),
                }
            return
        if int(getattr(self, 'dim', 0) or 0) != 3:
            return
        specs = self._projection_peak_specs_3d()
        if not specs:
            return

        ref_peaks = self.get_reference_peaks()
        for spec in specs:
            display_x, display_y = spec['display_axes']
            points = []
            for index, peak in enumerate(ref_peaks):
                axes = self._reference_axis_values_3d(peak)
                if display_x in axes and display_y in axes:
                    points.append({
                        'x': float(axes[display_x]),
                        'y': float(axes[display_y]),
                        'label': getattr(peak, 'name', ''),
                        'source_index': index,
                        'source': 'reference',
                        'axis_values': dict(axes),
                    })
            store.projected_peak_lists[('reference', display_x, display_y)] = {
                'peaks': points,
                'source': 'reference',
                'display_axes': (display_x, display_y),
            }

        # Physical 3p displays always use the 2D reference list.  Publish an
        # overlay keyed by the canonical Pseudo3D X/Y labels so the viewer can
        # consume projected_peak_lists directly without knowing physical axis
        # order or rebuilding a private peak representation.
        pseudo_view = self.get_pseudo3d_view('raw') if self._is_pseudo3d_topology() else None
        if pseudo_view is not None:
            display_x = str(pseudo_view['x_label'])
            display_y = str(pseudo_view['y_label'])
            points = []
            for index, peak in enumerate(ref_peaks):
                points.append({
                    'x': float(peak.x),
                    'y': float(peak.y),
                    'label': getattr(peak, 'name', ''),
                    'source_index': index,
                    'source': 'reference',
                    'axis_values': {
                        display_x: float(peak.x),
                        display_y: float(peak.y),
                    },
                })
            store.projected_peak_lists[('reference', display_x, display_y)] = {
                'peaks': points,
                'source': 'reference',
                'display_axes': (display_x, display_y),
            }

        full = self.get_full_peak_payload()
        records = full.get('records') or []
        for spec in specs:
            display_x, display_y = spec['display_axes']
            points = []
            for record in records:
                axes = record.get('axis_values', {})
                if display_x in axes and display_y in axes:
                    points.append({
                        'x': float(axes[display_x]),
                        'y': float(axes[display_y]),
                        'label': record.get('name', ''),
                        'source_index': record.get('row_index'),
                        'source': 'full',
                        'axis_values': dict(axes),
                    })
            store.projected_peak_lists[('full', display_x, display_y)] = {
                'peaks': points,
                'source': 'full',
                'display_axes': (display_x, display_y),
            }

    def get_projected_peak_overlay(self, left, right, transpose='n', source=None):
        """Return centrally prepared markers for a 3D projection display.

        ``left/right/transpose`` identify the same spectrum view requested by
        the GUI.  They are resolved here to the view's final display X/Y axes;
        peak coordinates themselves are never transposed.
        """
        store = getattr(self, 'store', None)
        if store is None:
            return []
        spec = self._projection_peak_spec_3d(left, right, transpose)
        if spec is None:
            return []
        display_x, display_y = spec['display_axes']
        if not store.projected_peak_lists:
            self._rebuild_projected_peak_lists()

        chosen = source
        if chosen is None:
            specs = self._projection_peak_specs_3d()
            use_reference_first = False
            if specs and spec == specs[0]:
                try:
                    use_reference_first = bool(self.cb_decon3d.IsChecked())
                except Exception:
                    use_reference_first = False
            full_payload = store.projected_peak_lists.get(
                ('full', display_x, display_y), {}
            )
            if use_reference_first:
                chosen = 'reference'
            elif full_payload.get('peaks'):
                chosen = 'full'
            else:
                chosen = 'reference'

        payload = store.projected_peak_lists.get(
            (chosen, display_x, display_y), {}
        )
        peaks = [dict(peak) for peak in (payload.get('peaks') or [])]

        selection = self._get_peak_selection()
        if not selection:
            for peak in peaks:
                peak['color'] = '#000000'
                peak['selected'] = False
            return peaks

        selection_source = str(selection.get('source') or '')
        selection_name = str(selection.get('name') or '')
        reference_name = str(selection.get('reference_name') or '')
        ref_axes = dict(selection.get('axis_values') or {})
        labels3d = list(getattr(self, 'labb', []) or [])

        selected_ref_index = selection.get('reference_index')
        for peak in peaks:
            peak['color'] = '#000000'
            peak['selected'] = False
            axes = dict(peak.get('axis_values') or {})

            if selection_source == 'full' and peak.get('source') == 'full':
                if str(peak.get('label', '')) == selection_name:
                    peak['color'] = '#2ca02c'
                    peak['selected'] = True

            elif selection_source == 'reference':
                if peak.get('source') == 'reference':
                    same_ref = False
                    if selected_ref_index is not None and peak.get('source_index') == selected_ref_index:
                        same_ref = True
                    elif reference_name and str(peak.get('label', '')) == reference_name:
                        same_ref = True
                    if same_ref:
                        peak['color'] = '#2ca02c'
                        peak['selected'] = True
                elif peak.get('source') == 'full' and len(labels3d) >= 3 and ref_axes:
                    label_dim2 = str(labels3d[2])
                    label_dim1 = str(labels3d[1])
                    ref_x = ref_axes.get(label_dim2)
                    ref_y = ref_axes.get(label_dim1)
                    if ref_x is not None and ref_y is not None and label_dim2 in axes and label_dim1 in axes:
                        tol_x = self._axis_step_tolerance(label_dim2)
                        tol_y = self._axis_step_tolerance(label_dim1)
                        on_reference_line = (
                            abs(float(axes[label_dim2]) - float(ref_x)) <= tol_x
                            and abs(float(axes[label_dim1]) - float(ref_y)) <= tol_y
                        )
                        if on_reference_line:
                            peak['color'] = '#2ca02c'
                            peak['selected'] = True

        return peaks

    def _notify_analysis_changed(self):
        """Tell the notebook that observable scientific evidence changed.

        Workflow remains read-only: this only requests a re-evaluation of the
        DataStore/project evidence after an existing scientific operation has
        successfully committed its result.
        """
        notebook = getattr(self, 'parent', None)
        notify = getattr(notebook, 'notify_analysis_changed', None)
        if callable(notify):
            notify()

    def set_reference_peaks(self, peaks, source_path=None):
        peaks = list(peaks or [])
        self.peak = peaks  # legacy compatibility: self.peak always means reference peaks
        self.PEAK = 1 if peaks else 0
        if int(getattr(self, 'dim', 0) or 0) == 3 and getattr(self, 'data', None) is not None:
            self.pkIdx = []
            self.pkSlice1D = []
            for pk in peaks:
                if hasattr(pk, 'indexJ') and hasattr(pk, 'indexK'):
                    self.pkIdx.append((pk.indexJ, pk.indexK))
                    self.pkSlice1D.append(self.data[:, pk.indexJ, pk.indexK])
        payload = {
            'peaks': peaks,
            'dimension': 2,
            'source_path': source_path if source_path is not None else self._resolve_spec_file(self.referencePeakBox.GetValue()),
            'pkIdx': self.pkIdx,
            'pkSlice1D': self.pkSlice1D,
        }
        self.store.save_peak_list('reference', **payload)
        self._rebuild_projected_peak_lists()
        self._notify_analysis_changed()
        return peaks

    def classify_reference_peak(self, peak_index, residue_type):
        """Classify one reference 2D peak selected by its list row.

        Classification preserves the legacy peakFrame behaviour: if the peak
        name begins with a digit the residue identifier is prepended; otherwise
        the first character is replaced.  The canonical reference list and all
        open views are then refreshed from the controller.
        """
        peaks = self.get_reference_peaks()
        try:
            peak_index = int(peak_index)
            peak = peaks[peak_index]
        except (TypeError, ValueError, IndexError):
            raise ValueError('The selected reference peak is no longer available.')

        residue_type = str(residue_type or '').strip()
        if len(residue_type) != 1:
            raise ValueError('Residue ID must contain exactly one character.')

        old_name = str(getattr(peak, 'name', '') or '')
        if not old_name:
            new_name = residue_type
        elif old_name[0].isdigit():
            new_name = residue_type + old_name
        else:
            new_name = residue_type + old_name[1:]

        peak.name = new_name
        self.set_reference_peaks(peaks)
        self.refresh_reference_peak_views(selected_name=new_name)

        frame = getattr(self, 'peak_frame', None)
        if frame is not None:
            try:
                frame.draw_figure()
            except (RuntimeError, wx.PyDeadObjectError):
                self.peak_frame = None
        return old_name, new_name

    # ------------------------------------------------------------------
    # Shared Pseudo3D model/controller API
    # ------------------------------------------------------------------
    def get_pseudo3d_groups(self):
        """Return the authoritative Pseudo3D overlap-group mapping."""
        return self.store.Grps

    def replace_pseudo3d_groups(self, groups):
        """Replace overlap groups in the shared DataStore, preserving ownership."""
        self.store.Grps.clear()
        self.store.Grps.update({str(k): list(v) for k, v in (groups or {}).items()})
        return self.store.Grps

    def add_pseudo3d_group(self, name, peaks=None):
        self.store.Grps[str(name)] = list(peaks or [])
        return self.store.Grps[str(name)]

    def remove_pseudo3d_group(self, name):
        return self.store.Grps.pop(str(name), None)

    def add_peak_to_pseudo3d_group(self, name, peak_name):
        group = self.store.Grps.setdefault(str(name), [])
        if peak_name not in group:
            group.append(peak_name)
        return group

    def remove_peak_from_pseudo3d_group(self, name, peak_name):
        group = self.store.Grps.get(str(name), [])
        if peak_name in group:
            group.remove(peak_name)
        return group

    def get_parameter_value(self, name, default=''):
        """Read a project parameter without exposing parameter-file plumbing to views."""
        path = self.state.parameter_file if getattr(self, 'state', None) is not None else self.deconParFile
        return Parse(path, name, default=default)

    def get_parameter_float(self, name, default=0.0):
        path = self.state.parameter_file if getattr(self, 'state', None) is not None else self.deconParFile
        return ParseFlt(path, name, default=default)

    def update_pseudo3d_parameters(self, values):
        """Persist Pseudo3D parameters through the project/controller boundary."""
        path = self.state.parameter_file if getattr(self, 'state', None) is not None else self.deconParFile
        if not path:
            raise RuntimeError('No project parameter file is configured')
        update_parameter_file(path, values, source_path=path)
        if getattr(self, 'state', None) is not None:
            self.state.dirty = True

    def get_threshold_fraction(self):
        """Return the current contour/noise threshold fraction via the controller."""
        try:
            return float(self.threshBox.GetValue())
        except Exception:
            return float(self.get_parameter_float('thresh', default=0.05))

    def get_project_input_file(self):
        """Return the canonical project-relative spectrum filename."""
        if getattr(self, 'state', None) is not None:
            return self.state.input_file
        return self.infileBox.GetValue() if hasattr(self, 'infileBox') else ''

    def get_noise_sigma(self):
        """Return the authoritative spectrum noise sigma, or ``None`` if unavailable."""
        candidates = []
        stats = getattr(self, 'noiseStats', None)
        if isinstance(stats, dict):
            candidates.append(stats.get('noise_sigma'))
        candidates.append(getattr(self, 'noiseVal', None))
        store = getattr(self, 'store', None)
        if store is not None:
            candidates.append(getattr(store, 'noiseVal', None))
        for value in candidates:
            try:
                sigma = float(value)
            except (TypeError, ValueError):
                continue
            if numpy.isfinite(sigma) and sigma > 0.0:
                return sigma
        return None

    def intensity_to_snr(self, intensity):
        """Convert intensity value(s) to S/N without modifying authoritative data."""
        sigma = self.get_noise_sigma()
        if sigma is None:
            return None
        return numpy.asarray(intensity) / sigma

    def get_reference_1d_view(self, peak_index):
        """Return the centrally-owned 1D trace view for a reference peak.

        The trace is the direct dimension-0 cut through the centre of the
        selected reference 2D peak in a 3D spectrum.  GUI consumers should
        render this payload rather than extracting slices from spectrum data
        themselves.
        """
        try:
            peak_index = int(peak_index)
        except (TypeError, ValueError):
            return None

        peaks = list(self.get_reference_peaks() or [])

        # DataStore values may be Python sequences or NumPy arrays.  Never use
        # ``value or []`` here: NumPy deliberately rejects truth-value testing
        # for multi-element arrays.
        raw_slices_value = getattr(self.store, 'pkSlice1D', None)
        raw_slices = [] if raw_slices_value is None else list(raw_slices_value)
        if peak_index < 0 or peak_index >= len(peaks) or peak_index >= len(raw_slices):
            return None

        x_value = getattr(self.store, 'index0', None)
        x = numpy.asarray([] if x_value is None else x_value)
        raw = numpy.asarray(raw_slices[peak_index])
        if x.size == 0 or raw.size == 0:
            return None

        decon = None
        decon_slices_value = getattr(self.store, 'pkSlice1Ddec', None)
        decon_slices = [] if decon_slices_value is None else list(decon_slices_value)
        if peak_index < len(decon_slices):
            candidate = numpy.asarray(decon_slices[peak_index])
            if candidate.size == raw.size:
                decon = candidate

        try:
            threshold = float(self.store.dmax) * float(self.threshBox.GetValue())
        except Exception:
            try:
                threshold = float(self.store.dmax)
            except Exception:
                threshold = 0.0

        # Peak-location stems in Slice1D are derived from the authoritative
        # Full Peak List.  conn_data used to provide these markers, but it is
        # deliberately quarantined and must not be reintroduced as a second
        # peak authority.  Full-list coordinates are already normalised into
        # axis_values by _full_peak_records(), so the Slice1D viewer only
        # needs the coordinate on the displayed (dimension-0) axis and the
        # fitted/picked peak intensity.
        markers = []
        full_payload = self.get_full_peak_payload() or {}
        full_records = list(full_payload.get('records') or [])
        reference_name = str(getattr(peaks[peak_index], 'name', peak_index))
        x_label = str(self.labb[0]) if getattr(self, 'labb', None) else None
        for record in full_records:
            name = str(record.get('name', ''))
            # In higher-dimensional journeys the Full peaks belonging to a
            # reference peak use the canonical <reference>_<number> identity.
            # Keep the matching rule aligned with Slice2D's 1D lower panel.
            if re.match(r'^%s_(\d+)$' % re.escape(reference_name), name) is None:
                continue
            axis_values = record.get('axis_values') or {}
            if x_label not in axis_values:
                continue
            try:
                xpos = float(axis_values[x_label])
                height = float(record.get('intensity'))
            except (TypeError, ValueError):
                continue
            if not numpy.isfinite(xpos) or not numpy.isfinite(height):
                continue
            suffix = re.search(r'_(\d+)$', name)
            markers.append({
                'x': xpos,
                'height': height,
                'label': suffix.group(1) if suffix else name,
                'name': name,
            })

        return {
            'x': x,
            'raw': raw,
            'decon': decon,
            'threshold': threshold,
            'peak': peaks[peak_index],
            'label': str(getattr(peaks[peak_index], 'name', peak_index)),
            'x_label': (str(self.labb[0]) + ' (ppm)') if getattr(self, 'labb', None) else 'ppm',
            'noise_sigma': self.get_noise_sigma(),
            'markers': markers,
        }

    def get_full3d_view_spec(self, bore_dim):
        """Return the authoritative orientation description for the Full3D viewer.

        Array dimensions are kept separate from display axes so every Full3D
        consumer (plane, slider, cross-sections and peak overlays) interprets
        transpose/orientation identically.
        """
        if int(getattr(self, 'dim', 0) or 0) != 3:
            return None
        labels = list(getattr(self.store, 'labb', None) or getattr(self, 'labb', []) or [])
        scales = []
        for name in ('uc0', 'uc1', 'uc2'):
            uc = getattr(self.store, name, None)
            if uc is None or not hasattr(uc, 'ppms_scale'):
                return None
            scales.append(numpy.asarray(uc.ppms_scale))
        if len(labels) < 3:
            labels = [str(i) for i in range(3)]
        bore_dim = int(bore_dim) % 3
        mapping = {
            2: (2, 0, 1),  # data[:, :, n].T -> x dim0, y dim1
            1: (1, 0, 2),  # data[:, n, :].T -> x dim0, y dim2
            0: (0, 1, 2),  # data[n, :, :].T -> x dim1, y dim2
        }
        slice_dim, x_dim, y_dim = mapping[bore_dim]
        return {
            'bore_dim': bore_dim,
            'slice_dim': slice_dim,
            'x_dim': x_dim,
            'y_dim': y_dim,
            'slice_scale': scales[slice_dim],
            'x_scale': scales[x_dim],
            'y_scale': scales[y_dim],
            'slice_label': str(labels[slice_dim]),
            'x_label': str(labels[x_dim]),
            'y_label': str(labels[y_dim]),
        }

    def get_full3d_slice_view(self, bore_dim, slice_index):
        """Return one plotting-ready plane from the main and decon 3D cubes."""
        spec = self.get_full3d_view_spec(bore_dim)
        data = getattr(self.store, 'data', None)
        if spec is None or data is None or numpy.ndim(data) != 3:
            return None
        size = int(data.shape[spec['slice_dim']])
        if size <= 0:
            return None
        n = max(0, min(int(slice_index), size - 1))
        if spec['slice_dim'] == 2:
            raw = numpy.asarray(data[:, :, n]).T
        elif spec['slice_dim'] == 1:
            raw = numpy.asarray(data[:, n, :]).T
        else:
            raw = numpy.asarray(data[n, :, :]).T

        decon = None
        datadec = getattr(self.store, 'datadec', None)
        if datadec is not None and numpy.ndim(datadec) == 3:
            try:
                if spec['slice_dim'] == 2:
                    decon = numpy.asarray(datadec[:, :, n]).T
                elif spec['slice_dim'] == 1:
                    decon = numpy.asarray(datadec[:, n, :]).T
                else:
                    decon = numpy.asarray(datadec[n, :, :]).T
            except Exception:
                decon = None

        view = dict(spec)
        view.update({
            'slice_index': n,
            'slice_value': float(spec['slice_scale'][n]),
            'raw': raw,
            'decon': decon,
        })
        return view

    def get_full3d_intensity_limits(self):
        """Return fixed raw 3D intensity limits for all Full3D slices.

        The limits are derived once from the complete raw 3D cube and cached
        in the shared DataStore metadata. They therefore do not change when
        the displayed slice or orientation changes.
        """
        data = getattr(self.store, 'data', None)
        if data is None or numpy.ndim(data) != 3:
            return None
        cached_id = self.store.metadata.get('full3d_intensity_data_id')
        cached_limits = self.store.metadata.get('full3d_intensity_limits')
        if cached_id == id(data) and cached_limits is not None:
            return tuple(cached_limits)
        arr = numpy.asarray(data)
        finite = arr[numpy.isfinite(arr)]
        if finite.size == 0:
            return None
        lo = float(numpy.min(finite))
        hi = float(numpy.max(finite))
        if lo == hi:
            pad = abs(lo) * 0.05 or 1.0
            lo -= pad
            hi += pad
        limits = (lo, hi)
        self.store.metadata['full3d_intensity_data_id'] = id(data)
        self.store.metadata['full3d_intensity_limits'] = limits
        return limits

    def get_full3d_cross_sections(self, bore_dim, slice_index, x_ppm, y_ppm):
        """Return horizontal/vertical traces through a Full3D display plane."""
        view = self.get_full3d_slice_view(bore_dim, slice_index)
        data = getattr(self.store, 'data', None)
        if view is None or data is None or x_ppm is None or y_ppm is None:
            return None
        x_scale = numpy.asarray(view['x_scale'])
        y_scale = numpy.asarray(view['y_scale'])
        xi = int(numpy.argmin(numpy.abs(x_scale - float(x_ppm))))
        yi = int(numpy.argmin(numpy.abs(y_scale - float(y_ppm))))
        n = int(view['slice_index'])

        if view['slice_dim'] == 2:
            horizontal = numpy.asarray(data[:, yi, n])
            vertical = numpy.asarray(data[xi, :, n])
        elif view['slice_dim'] == 1:
            horizontal = numpy.asarray(data[:, n, yi])
            vertical = numpy.asarray(data[xi, n, :])
        else:
            horizontal = numpy.asarray(data[n, :, yi])
            vertical = numpy.asarray(data[n, xi, :])

        horizontal_decon = vertical_decon = None
        datadec = getattr(self.store, 'datadec', None)
        if datadec is not None and numpy.ndim(datadec) == 3:
            try:
                if view['slice_dim'] == 2:
                    horizontal_decon = numpy.asarray(datadec[:, yi, n])
                    vertical_decon = numpy.asarray(datadec[xi, :, n])
                elif view['slice_dim'] == 1:
                    horizontal_decon = numpy.asarray(datadec[:, n, yi])
                    vertical_decon = numpy.asarray(datadec[xi, n, :])
                else:
                    horizontal_decon = numpy.asarray(datadec[n, :, yi])
                    vertical_decon = numpy.asarray(datadec[n, xi, :])
            except Exception:
                horizontal_decon = vertical_decon = None

        return {
            'x_index': xi,
            'y_index': yi,
            'x_ppm': float(x_scale[xi]),
            'y_ppm': float(y_scale[yi]),
            'horizontal_axis': x_scale,
            'horizontal': horizontal,
            'horizontal_decon': horizontal_decon,
            'vertical_axis': y_scale,
            'vertical': vertical,
            'vertical_decon': vertical_decon,
        }

    @staticmethod
    def full_peak_slice_color(slice_delta):
        """Return the shared greyscale used for Full peaks by slice distance.

        The current plane is darkest; peaks one and two planes away become
        progressively lighter.  Full3D and Slice2D deliberately share this
        rule so their 3D context is visually consistent.
        """
        return str(min(0.15 + 0.25 * abs(int(slice_delta)), 0.85))

    def get_full3d_peak_overlay(self, bore_dim, slice_index, max_slice_distance=2):
        """Return Full-list peaks projected onto the current Full3D plane.

        Selection/highlight state is kept separate from the scientific peak
        list.  Reference selections can highlight their corresponding Full
        peak and, in bore mode with Use 2D peaklist enabled, other Full peaks
        sharing the selected slice receive a secondary highlight.
        """
        view = self.get_full3d_slice_view(bore_dim, slice_index)
        if view is None:
            return []
        payload = self.get_full_peak_payload()
        records = list(payload.get('records') or [])
        if not records:
            return []
        slice_scale = numpy.asarray(view['slice_scale'])
        current = int(view['slice_index'])
        x_label = view['x_label']
        y_label = view['y_label']
        slice_label = view['slice_label']
        selection = self._get_peak_selection()
        selected_name = str(selection.get('name')) if selection and selection.get('source') == 'full' else None
        selected_reference_name = str(selection.get('reference_name') or selection.get('name')) if selection and selection.get('source') == 'reference' else None
        selected_ref_axes = {}
        if selection and selection.get('source') == 'reference':
            selected_ref_axes = dict(selection.get('axis_values') or {})

        # In bore mode with Use 2D peaklist enabled, all Full peaks on the
        # selected reference 1D slice get a light-green secondary highlight.
        bore_reference_slice = False
        ref_slice_value = None
        try:
            bore_reference_slice = bool(self.cb_decon3d.IsChecked()) and selection and selection.get('source') == 'reference'
        except Exception:
            bore_reference_slice = False
        if bore_reference_slice and slice_label in selected_ref_axes:
            ref_slice_value = float(selected_ref_axes[slice_label])

        markers = []
        for record in records:
            axes = record.get('axis_values', {})
            if x_label not in axes or y_label not in axes or slice_label not in axes:
                continue
            peak_slice_index = int(numpy.argmin(numpy.abs(slice_scale - float(axes[slice_label]))))
            delta = peak_slice_index - current
            if abs(delta) > int(max_slice_distance):
                continue
            color = self.full_peak_slice_color(delta)
            selected = False
            same_ref_slice = False
            on_reference_line = False
            if selected_name is not None and str(record.get('name', '')) == selected_name:
                selected = True
                color = '#2ca02c'
            elif selected_reference_name is not None and selected_ref_axes:
                # A Reference 2D peak defines a 1D line through a 3D cube.
                # The existing 1D-slice logic is data[:, indexJ, indexK]:
                # dimensions 1 and 2 are fixed by the reference y/x values,
                # while dimension 0 is free.  Therefore a Full 3D peak is on
                # the selected reference line iff its dim-2 and dim-1
                # coordinates match the reference peak; the current Full3D
                # slice dimension is allowed to vary.
                labels3d = list(getattr(self, 'labb', []) or [])
                if len(labels3d) >= 3:
                    label_dim2 = str(labels3d[2])
                    label_dim1 = str(labels3d[1])
                    ref_x = selected_ref_axes.get(label_dim2)
                    ref_y = selected_ref_axes.get(label_dim1)
                    if ref_x is not None and ref_y is not None:
                        tol_x = self._axis_step_tolerance(label_dim2)
                        tol_y = self._axis_step_tolerance(label_dim1)
                        on_reference_line = (
                            label_dim2 in axes and label_dim1 in axes
                            and abs(float(axes[label_dim2]) - float(ref_x)) <= tol_x
                            and abs(float(axes[label_dim1]) - float(ref_y)) <= tol_y
                        )
                        if on_reference_line:
                            selected = True
                            color = '#2ca02c'
            if not selected and bore_reference_slice and ref_slice_value is not None:
                tolerance = self._axis_step_tolerance(slice_label)
                same_ref_slice = abs(float(axes[slice_label]) - ref_slice_value) <= tolerance
                # This secondary highlight is only meaningful for the old
                # reference-plane indication.  Do not use it to colour the
                # entire plane when a reference selection defines a line.
                if same_ref_slice and delta == 0 and not on_reference_line:
                    color = '#8fd694'
            markers.append({
                'x': float(axes[x_label]),
                'y': float(axes[y_label]),
                'label': str(record.get('name', '')),
                'slice_delta': int(delta),
                'color': color,
                'source_index': record.get('row_index'),
                'selected': selected,
                'same_reference_slice': same_ref_slice,
                'on_reference_line': on_reference_line,
            })
        return markers

    def _get_peak_selection(self):
        selection = self.store.metadata.get('peak_selection') if getattr(self, 'store', None) else None
        return dict(selection) if isinstance(selection, dict) else None

    def clear_peak_selection(self, redraw_full3d=True):
        if getattr(self, 'store', None) is not None:
            self.store.metadata.pop('peak_selection', None)
        if redraw_full3d:
            viewer = getattr(self.parent, 'tabFive', None) if getattr(self, 'parent', None) is not None else None
            if viewer is not None:
                try:
                    viewer.draw_figure(keepaxes=True)
                except Exception:
                    pass

    def _full3d_viewer(self, ensure=False):
        if int(getattr(self, 'dim', 0) or 0) != 3 or getattr(self, 'parent', None) is None:
            return None
        viewer = getattr(self.parent, 'tabFive', None)
        exists = False
        try:
            exists = self.parent.PageExists('Full 3D')
        except Exception:
            exists = viewer is not None
        if viewer is None and ensure and not exists:
            try:
                self.parent.AddTabFive(True, self)
                viewer = getattr(self.parent, 'tabFive', None)
            except Exception:
                viewer = None
        return viewer

    def _reference_record_for_name(self, name):
        for idx, peak in enumerate(self.get_reference_peaks()):
            if str(getattr(peak, 'name', '')) == str(name):
                return idx, peak
        return None, None

    def _full_record_for_name(self, name):
        for record in list(self.get_full_peak_payload().get('records') or []):
            if str(record.get('name', '')) == str(name):
                return record
        return None

    def _axis_step_tolerance(self, label):
        labels = list(getattr(self, 'labb', []) or [])
        for idx, lab in enumerate(labels[:3]):
            if str(lab) != str(label):
                continue
            uc = getattr(self.store, f'uc{idx}', None)
            scale = getattr(uc, 'ppms_scale', None) if uc is not None else None
            if scale is None:
                return 1e-3
            arr = numpy.asarray(scale)
            if arr.size < 2:
                return 1e-3
            diffs = numpy.abs(numpy.diff(arr))
            diffs = diffs[numpy.isfinite(diffs)]
            if diffs.size == 0:
                return 1e-3
            return max(float(numpy.median(diffs)) * 1.5, 1e-5)
        return 1e-3

    def _reference_for_full_record(self, record):
        name = str(record.get('name', ''))
        exact_idx, exact = self._reference_record_for_name(name)
        if exact is not None:
            return exact_idx, exact
        axes = record.get('axis_values', {})
        if int(getattr(self, 'dim', 0) or 0) != 3 or not axes:
            return None, None
        best = None
        best_score = float('inf')
        labels = list(getattr(self, 'labb', []) or [])
        if len(labels) < 3:
            return None, None
        for idx, peak in enumerate(self.get_reference_peaks()):
            ref_axes = self._reference_axis_values_3d(peak)
            score = 0.0
            valid = True
            # A reference peak defines the restricted XY bore only.  The
            # full peak's third coordinate is free to vary along that bore and
            # must not participate in 2D<->3D source association.
            for label in (str(labels[2]), str(labels[1])):
                if label not in axes or label not in ref_axes:
                    valid = False
                    break
                tol = self._axis_step_tolerance(label)
                score += abs(float(axes[label]) - float(ref_axes[label])) / tol
            if valid and score < best_score:
                best = (idx, peak)
                best_score = score
        if best is not None and best_score <= 6.0:
            return best
        return None, None

    def _selection_axis_values(self, source, name=None, record=None, peak=None):
        if source == 'reference':
            if peak is None:
                _, peak = self._reference_record_for_name(name)
            if peak is None:
                return {}
            return self._reference_axis_values_3d(peak)
        if record is None and name is not None:
            record = self._full_record_for_name(name)
        return dict(record.get('axis_values') or {}) if record is not None else {}

    def _reference_selection_for_full_record(self, record):
        idx, peak = self._reference_for_full_record(record)
        return idx, peak

    def _physical_2d_peak_for_full_record(self, record):
        """Return a peak-like object for an authoritative physical-2D Full record.

        Physical 2D deliberately has no shadow Reference list, so a Full-list
        selection cannot be resolved through ``_reference_for_full_record``.
        Prefer the live Get 2D Peaks object's Peak instance when that editor is
        open; otherwise build the small x/y/name view required by Projection.
        """
        if not isinstance(record, dict):
            return None
        name = str(record.get('name', ''))
        frame = self._live_peak_frame()
        if frame is not None:
            for peak in list(getattr(frame, 'peak', []) or []):
                if str(getattr(peak, 'name', '')) == name:
                    return peak

        axes = record.get('axis_values') or {}
        labels = list(getattr(self, 'labb', []) or [])
        try:
            if len(labels) >= 2 and str(labels[1]) in axes and str(labels[0]) in axes:
                x = float(axes[str(labels[1])])
                y = float(axes[str(labels[0])])
            else:
                # Full-list files are f1,f2 while the displayed 2D plane is
                # X=f2 (direct) and Y=f1 (indirect).
                coords = tuple(record.get('coordinates') or ())
                if len(coords) < 2:
                    return None
                y, x = float(coords[0]), float(coords[1])
        except (KeyError, TypeError, ValueError):
            return None

        from types import SimpleNamespace
        return SimpleNamespace(name=name, x=x, y=y)

    def OnButtonRestricted3DDiagnostics(self, event=None):
        frame = getattr(self, 'restricted3dDiagnosticsFrame', None)
        if frame is None:
            frame = Restricted3DDiagnosticsFrame(self)
            self.restricted3dDiagnosticsFrame = frame
        frame.Show(); frame.Raise()

    def analyse_restricted_3d(self, **parameters):
        """Run the complete first-pass restricted-3D diagnostic suite."""
        return analyse_restricted_3d(self.store, self._reference_selection_for_full_record, **parameters)

    def analyse_restricted_3d_overlap(self, overlap_fraction=0.05):
        """Annotate reference/full peak records with the first diagnostic.

        This is intentionally GUI-independent apart from reusing the existing
        reference/full association rule.  Later diagnostic panels can call it
        after project load or from an explicit Analyse button.
        """
        return analyse_overlap(
            self.store,
            self._reference_selection_for_full_record,
            overlap_fraction=float(overlap_fraction),
        )

    def reclassify_restricted_3d_overlap(self, overlap_fraction=0.05):
        """Change the overlap threshold without rescanning the 3D data."""
        return classify_overlap(self.store, overlap_fraction=float(overlap_fraction))

    def _update_slice_combo(self, combo, name):
        if combo is None or name is None:
            return False
        idx = combo.FindString(str(name))
        if idx == wx.NOT_FOUND:
            return False
        combo.SetSelection(idx)
        return True

    def _update_slice_views_for_selection(self, selection, source_view=None, source_pane=None):
        ref_name = selection.get('reference_name')
        if not ref_name:
            return
        for attr, combos in (('tabThree', ('ComboBox1',)), ('tabFour', ('ComboBox1', 'ComboBox2'))):
            view = getattr(self.parent, attr, None)
            if view is None:
                continue
            # A selection made inside Slice2D must not navigate its independent
            # top and bottom reference selectors together.  In the normal 2D/2D
            # view only the pane that originated the click is eligible for
            # navigation.  Decon/Orth/1D are derived from the top/left slice, so
            # their lower pane intentionally continues to follow ComboBox1.
            view_combos = combos
            if attr == 'tabFour' and source_view is view and source_pane in ('top', 'bottom'):
                dependent = False
                try:
                    dependent = bool(view.checkDecon.GetValue() or
                                     view.checkOrth.GetValue() or
                                     view.check1D.GetValue())
                except Exception:
                    pass
                if dependent:
                    view_combos = ('ComboBox1',)
                elif source_pane == 'top':
                    view_combos = ('ComboBox1',)
                else:
                    view_combos = ('ComboBox2',)
            changed = False
            for combo_name in view_combos:
                combo = getattr(view, combo_name, None)
                changed = self._update_slice_combo(combo, ref_name) or changed
            if changed:
                # SetSelection/SetValue performed programmatically does not emit
                # wx.EVT_COMBOBOX.  The 2D slice viewer uses that event to reset
                # its cached axis limits before drawing a different reference
                # peak.  Route programmatic cross-view navigation through the
                # same handler so the data *and* ppm axes move together.
                peak_change_handler = getattr(view, 'on_peak_combo_changed', None)
                if callable(peak_change_handler):
                    try:
                        peak_change_handler(None)
                        continue
                    except (RuntimeError, wx.PyDeadObjectError):
                        continue
                    except Exception:
                        # Older/alternative slice viewers may expose a handler
                        # with different assumptions; retain the legacy redraw
                        # fallback below.
                        pass
                if hasattr(view, 'draw_figure'):
                    try:
                        view.draw_figure()
                    except TypeError:
                        view.draw_figure(False)

    def _set_full3d_selection(self, selection, ensure=False):
        viewer = self._full3d_viewer(ensure=ensure)
        if viewer is None:
            return
        if selection is None:
            try:
                viewer.draw_figure(keepaxes=True)
            except Exception:
                pass
            return
        axes = dict(selection.get('axis_values') or {})
        if not axes:
            return
        try:
            spec = viewer._spec()
            scale = numpy.asarray(spec['slice_scale']) if spec else None
            label = spec.get('slice_label') if spec else None
            if scale is not None and scale.size and label in axes:
                idx = int(numpy.argmin(numpy.abs(scale - float(axes[label]))))
                viewer.set_slice_index(idx, keepaxes=True)
            else:
                viewer.draw_figure(keepaxes=True)
        except Exception:
            try:
                viewer.draw_figure(keepaxes=True)
            except Exception:
                pass

    def _store_peak_selection(self, source, name, axis_values, reference_name=None, reference_index=None, record=None):
        self.store.metadata['peak_selection'] = {
            'source': str(source),
            'name': str(name),
            'axis_values': dict(axis_values or {}),
            'reference_name': None if reference_name is None else str(reference_name),
            'reference_index': reference_index,
            'full_record_index': None if record is None else record.get('row_index'),
        }
        return self.store.metadata['peak_selection']

    def select_reference_peak(self, peak_name):
        """Synchronise a 2D reference selection across slice views and Full3D."""
        peaks = self.get_reference_peaks()
        index, peak = self._reference_record_for_name(peak_name)
        if peak is None:
            return False
        selection = self._store_peak_selection(
            'reference', peak.name, self._reference_axis_values_3d(peak),
            reference_name=peak.name, reference_index=index,
        )
        self._update_slice_views_for_selection(selection)
        # A physical 2D Projection is another view of this same reference plane.
        # Repaint it from the shared selection so its peak ornament follows the
        # Reference Peak List viewer just like the Get 2D Peak List frame does.
        try:
            topology = topology_for(self)
            projection = getattr(self.parent, 'tabTwo', None)
            if topology.spectral_dim_count == 2 and not topology.has_pseudo_axis and projection is not None:
                focus_projection = getattr(projection, 'focus_2d_peak', None)
                if callable(focus_projection):
                    focus_projection(peak)
        except Exception as exc:
            print('2D reference selection -> Projection failed: %s' % exc)
        viewer = self._full3d_viewer(ensure=False)
        if viewer is not None:
            self._set_full3d_selection(selection, ensure=False)
        frame = getattr(self, 'peak_frame', None)
        if frame is not None:
            try:
                frame.focus_peak(peak)
            except (RuntimeError, wx.PyDeadObjectError):
                self.peak_frame = None
        return True

    def _select_pseudo2d_full_peak_in_projection(self, peak_name, record=None, zoom=True):
        """Route a Full 1D list selection directly to the pseudo2D Projection view.

        The generic Full-nD selection machinery is reference/slice oriented.  In
        pseudo2D the Projection panel *is* the reference view, so list ``Show``
        must address it directly rather than relying on the nD cross-view path.
        """
        try:
            topology = topology_for(self)
            if not (topology.spectral_dim_count == 1 and topology.has_pseudo_axis):
                return False
            notebook = self.parent
            projection = getattr(notebook, 'tabTwo', None)
            if projection is None:
                notebook.AddTabTwo(True, self)
                projection = getattr(notebook, 'tabTwo', None)
            if projection is None or not hasattr(projection, 'select_full_peak_from_list'):
                return False
            row_index = record.get('row_index') if isinstance(record, dict) else None
            selected = bool(projection.select_full_peak_from_list(
                peak_name, zoom=zoom, row_index=row_index))
            if not selected:
                return False
            # AddPage/_replace_tab does not necessarily select the new/existing
            # page.  Make Show visibly navigate to the Projection reference view.
            try:
                for idx in range(notebook.GetPageCount()):
                    if notebook.GetPage(idx) is projection:
                        notebook.SetSelection(idx)
                        break
            except Exception:
                pass
            try:
                projection.SetFocus()
            except Exception:
                pass
            return True
        except Exception as exc:
            # Do not silently swallow this route: a console diagnostic makes a
            # future platform-specific wx failure actionable without disturbing
            # normal GUI use.
            print('Pseudo2D Full Peak Show -> Projection failed: %s' % exc)
            return False

    def select_full_peak(self, peak_name, source_view=None, source_pane=None):
        """Synchronise a Full nD peak selection without touching peakFrame.

        ``source_view``/``source_pane`` preserve independent Slice2D top/bottom
        navigation when the selection originated in that viewer.
        """
        record = self._full_record_for_name(peak_name)
        if record is None:
            return False

        # Pseudo2D Full peaks are one-dimensional and the Projection panel is
        # their reference frame.  Route there first and return: running the nD
        # reference/slice synchronisation before this was the source of the
        # unreliable Full-list Show behaviour.
        try:
            topology = topology_for(self)
            if topology.spectral_dim_count == 1 and topology.has_pseudo_axis:
                return self._select_pseudo2d_full_peak_in_projection(
                    record.get('name', peak_name), record=record, zoom=True)
        except Exception as exc:
            print('Pseudo2D Full Peak Show topology/route failed: %s' % exc)
        # Physical 2D Full-list Show synchronises both open views of the same plane.
        try:
            topology = topology_for(self)
            if topology.spectral_dim_count == 2 and not topology.has_pseudo_axis:
                # Full 2D is the sole peak authority: there is intentionally no
                # Reference list to resolve through here.  Convert the selected
                # Full record directly to the shared displayed peak instead.
                selected_peak = self._physical_2d_peak_for_full_record(record)
                self._store_peak_selection(
                    'full', record.get('name', peak_name), record.get('axis_values') or {},
                    reference_name=getattr(selected_peak, 'name', None), reference_index=None, record=record,
                )
                if selected_peak is not None:
                    projection = getattr(self.parent, 'tabTwo', None)
                    focus_projection = getattr(projection, 'focus_2d_peak', None) if projection is not None else None
                    if callable(focus_projection):
                        focus_projection(selected_peak)
                    frame = self._live_peak_frame()
                    focus_peak = getattr(frame, 'focus_peak', None) if frame is not None else None
                    if callable(focus_peak):
                        focus_peak(selected_peak)
                return selected_peak is not None
        except Exception as exc:
            print('2D Full Peak Show synchronisation failed: %s' % exc)
        ref_index, ref_peak = self._reference_selection_for_full_record(record)
        ref_name = getattr(ref_peak, 'name', None) if ref_peak is not None else None
        selection = self._store_peak_selection(
            'full', record.get('name', peak_name), record.get('axis_values') or {},
            reference_name=ref_name, reference_index=ref_index, record=record,
        )
        self._update_slice_views_for_selection(selection, source_view=source_view, source_pane=source_pane)
        self._set_full3d_selection(selection, ensure=False)
        return True

    def get_full3d_peak_selection_context(self):
        """Return the current cross-view peak selection for Full3D overlays."""
        return self._get_peak_selection()

    def get_full_peak_axis_metadata(self):
        """Return alias controls for canonical Full-list coordinate columns.

        Coordinate columns are f1..fN, while ``labb``/spectrum dimensions run
        in reverse order.  For a 3D list generated from 2D planes (recognised
        from the complete ``nResID_Number`` naming convention), only f3 is
        aliasable; otherwise every coordinate dimension is exposed.
        """
        dim = int(getattr(self, 'dim', 0) or 0)
        spectral_axes = self._spectral_physical_axes()
        labels = [label for _, label in spectral_axes]
        rows = list(self.get_full_peak_payload().get('rows') or [])
        if dim <= 0 or not rows:
            raise ValueError('Full peak-list axis metadata is unavailable.')

        def split_2d_name(name):
            left, sep, right = str(name).rpartition('_')
            return bool(sep and left and re.search(r'\d+', right))

        special_3d = dim == 3 and all(fields and split_2d_name(fields[0]) for fields in rows)
        coord_indices = [2] if special_3d else list(range(dim))
        metadata = []
        for coord_index in coord_indices:
            spectral_dim = dim - 1 - coord_index
            spectrum_dim = spectral_axes[spectral_dim][0]
            vals = getattr(self, 'index%d' % spectrum_dim, None)
            if vals is None or len(vals) < 2:
                raise ValueError('Spectrum axis metadata is unavailable.')
            width = abs(float(max(vals)) - float(min(vals)) + abs(float(vals[1]) - float(vals[0])))
            label = str(labels[spectral_dim]) if 0 <= spectral_dim < len(labels) else 'f%d' % (coord_index + 1)
            metadata.append({
                'axis': coord_index,
                'label': label,
                'width_ppm': width,
                'spectrum_dim': spectrum_dim,
            })
        return metadata

    def alias_full_peak(self, peak_name, coord_index, direction):
        """Alias one authoritative Full-list coordinate and rebuild all views."""
        if direction not in (-1, 1):
            raise ValueError('Alias direction must be -1 or +1.')
        meta = next((m for m in self.get_full_peak_axis_metadata()
                     if int(m['axis']) == int(coord_index)), None)
        if meta is None:
            raise ValueError('This dimension is not available for Full peak aliasing.')
        payload = self.get_full_peak_payload()
        rows = payload.get('rows') or []
        target = None
        for fields in rows:
            if fields and str(fields[0]) == str(peak_name):
                target = fields
                break
        if target is None:
            raise ValueError('The selected Full peak is no longer available.')
        field_index = int(coord_index) + 1
        try:
            before = float(target[field_index])
        except (IndexError, TypeError, ValueError):
            raise ValueError('The selected peak coordinate is not numeric.')
        after = before + int(direction) * abs(float(meta['width_ppm']))
        # Rows are the single authoritative representation.  Derived records
        # and every projection/overlay are rebuilt from these edited values.
        target[field_index] = ('%.10g' % after)
        records = self._full_peak_records(rows, dim=int(getattr(self, 'dim', 0) or 0))
        payload['rows'] = rows
        payload['records'] = records
        payload['peaks'] = records
        self.store.save_peak_list('full', **payload)
        self._rebuild_projected_peak_lists()
        if getattr(self, 'state', None) is not None:
            self.state.dirty = True

        selection = self._get_peak_selection()
        if selection and selection.get('source') == 'full' and str(selection.get('name')) == str(peak_name):
            self.select_full_peak(peak_name)
        # Redraw open consumers. They all obtain peak positions from the Full
        # payload/projected lists rebuilt above.
        viewer = self._full3d_viewer(ensure=False)
        if viewer is not None:
            try:
                viewer.draw_figure(keepaxes=True)
            except Exception:
                pass
        parent = getattr(self, 'parent', None)
        if parent is not None:
            for attr in ('tabThree', 'tabFour'):
                view = getattr(parent, attr, None)
                if view is not None:
                    try:
                        view.draw_figure()
                    except Exception:
                        pass
        return before, after

    def get_reference_peak_axis_metadata(self):
        """Return display-axis metadata used by the shared 2D reference list.

        Each axis reports its label, one spectral-width step in ppm, and the
        underlying main-spectrum dimensions whose cached peak indices should
        be refreshed after a coordinate edit.
        """
        labels = list(getattr(self, 'labb', []) or [])
        dim = int(getattr(self, 'dim', 0) or 0)

        def axis_width(axis_index):
            vals = getattr(self, 'index%d' % axis_index, None)
            if vals is None or len(vals) < 2:
                raise ValueError('Spectrum axis metadata is unavailable.')
            return abs(float(max(vals)) - float(min(vals)) + abs(float(vals[1]) - float(vals[0])))

        if dim == 3 and len(labels) >= 3:
            return [
                {'axis': 'x', 'label': labels[2], 'width_ppm': axis_width(2), 'spectrum_dims': (2,)},
                # Preserve the historical reference-peak bookkeeping: the
                # projected y coordinate is cached against both I and J.
                {'axis': 'y', 'label': labels[1], 'width_ppm': axis_width(1), 'spectrum_dims': (0, 1)},
            ]
        if dim == 4 and len(labels) >= 4:
            return [
                {'axis': 'x', 'label': labels[2], 'width_ppm': axis_width(2), 'spectrum_dims': (2,)},
                {'axis': 'y', 'label': labels[3], 'width_ppm': axis_width(3), 'spectrum_dims': (3,)},
            ]
        if dim == 2 and len(labels) >= 2:
            pseudo = bool(self.state.pseudo_axis)
            if pseudo:
                spectral_axes = self._spectral_physical_axes()
                y_dim, y_label = spectral_axes[0]
                x_dim, x_label = spectral_axes[1]
                return [
                    {'axis': 'x', 'label': x_label, 'width_ppm': axis_width(x_dim), 'spectrum_dims': (x_dim,)},
                    {'axis': 'y', 'label': y_label, 'width_ppm': axis_width(y_dim), 'spectrum_dims': (y_dim,)},
                ]
            return [
                {'axis': 'x', 'label': labels[1], 'width_ppm': axis_width(1), 'spectrum_dims': (1,)},
                {'axis': 'y', 'label': labels[0], 'width_ppm': axis_width(0), 'spectrum_dims': (0,)},
            ]
        raise ValueError('Reference 2D axis metadata is unavailable.')

    def refresh_reference_peak_indices(self, peak):
        """Refresh cached nD indices for an edited reference peak."""
        for meta in self.get_reference_peak_axis_metadata():
            ppm = getattr(peak, meta['axis'])
            for spectrum_dim in meta['spectrum_dims']:
                self.alias(peak, ppm, spectrum_dim)

    def get_reference_peak_headers(self):
        """Column titles for the always-2D reference peak list.

        The coordinate labels follow the same display convention as peakFrame:
        for 3D the reference plane is labb[2] vs labb[1], for 4D it is
        labb[2] vs labb[3], and for native 2D data it is labb[1] vs labb[0].
        """
        labels = list(getattr(self, 'labb', []) or [])
        dim = int(getattr(self, 'dim', 0) or 0)
        if dim == 3 and len(labels) >= 3:
            x_label, y_label = labels[2], labels[1]
        elif dim == 4 and len(labels) >= 4:
            x_label, y_label = labels[2], labels[3]
        elif len(labels) >= 2:
            x_label, y_label = labels[1], labels[0]
        elif len(labels) == 1:
            x_label = y_label = labels[0]
        else:
            x_label, y_label = 'x', 'y'
        return ['ResID', 'Name', f'{x_label} (ppm)', f'{y_label} (ppm)']

    def get_full_peak_headers(self, row_width=None):
        """Column titles for the dimensionality-matched decon peak list.

        Decon list rows are written as ``Name, f1..fN, Intensity``.  The f
        coordinates run from fastest to slowest dimension, so their NMR labels
        are the main-spectrum labels in reverse order.
        """
        dim = int(getattr(self, 'dim', 0) or 0)
        labels = self._spectral_axis_labels()
        coord_labels = []
        for coord_index in range(dim):
            lab_index = dim - 1 - coord_index
            if 0 <= lab_index < len(labels):
                label = labels[lab_index]
            else:
                label = f'f{coord_index + 1}'
            coord_labels.append(f'{label} (ppm)')

        headings = ['Name'] + coord_labels + ['Intensity']
        if row_width is not None:
            row_width = int(row_width)
            if row_width < len(headings):
                headings = headings[:row_width]
            elif row_width > len(headings):
                # Preserve a useful title for unexpected trailing legacy fields
                # without moving the decon intensity column away from dim + 1.
                headings.extend(
                    f'Field {i}' for i in range(len(headings), row_width)
                )
        return headings

    def get_full_peak_payload(self):
        return self.store.peak_lists.get('full', {}) if getattr(self, 'store', None) is not None else {}

    def _read_full_peak_rows(self, infile):
        rows = []
        with open(infile, 'r') as handle:
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                try:
                    float(fields[1])
                except (IndexError, ValueError):
                    continue
                rows.append(fields)
        return rows

    def load_full_peak_list(self, path=None, quiet=False):
        """Load the dimensionality-matched full list without touching reference peaks."""
        path = str(path or self.fullPeakBox.GetValue() or '').strip()
        if not path:
            return False
        if not os.path.isabs(path):
            candidate = self._resolve_spec_file(path)
            if os.path.exists(candidate):
                path = candidate
        if not os.path.exists(path):
            if not quiet:
                errorMessage('Cannot find full peak list.')
            return False
        rows = self._read_full_peak_rows(path)
        spectral_dim_count = self._active_topology().spectral_dim_count
        records = self._full_peak_records(rows, dim=spectral_dim_count)
        self.store.save_peak_list(
            'full', peaks=records, records=records, rows=rows,
            dimension=spectral_dim_count, source_path=path
        )
        self._rebuild_projected_peak_lists()
        rel_path = self.state._spec_relative(path) if getattr(self, 'state', None) is not None else path
        self.fullPeakBox.SetValue(rel_path)
        if getattr(self, 'state', None) is not None:
            self.state.full_peak_file = rel_path
        self.corrFile = path
        self._notify_analysis_changed()
        # Full-list peaks are the only picked peaks in pseudo2D mode, so the
        # Report panel should update immediately when this list is loaded.
        self.Status()
        return True

    def refresh_reference_peak_views(self, selected_name=None):
        """Refresh combo choices in open slice views while preserving selection."""
        names = [pk.name for pk in self.get_reference_peaks()]
        for attr, combos in (('tabThree', ('ComboBox1',)), ('tabFour', ('ComboBox1', 'ComboBox2'))):
            view = getattr(self.parent, attr, None)
            if view is None:
                continue
            for combo_name in combos:
                combo = getattr(view, combo_name, None)
                if combo is None:
                    continue
                current = selected_name or combo.GetValue()
                combo.SetItems(names)
                idx = combo.FindString(current) if current else wx.NOT_FOUND
                if idx == wx.NOT_FOUND and names:
                    idx = 0
                if idx != wx.NOT_FOUND:
                    combo.SetSelection(idx)

    def _projection_view_payloads(self, dic, data, file_labb, source, decon=False):
        """Build immutable plotting payloads for a cached 2D projection.

        Projection files are interpreted once here, in the controller/data
        layer.  The Projection GUI consumes the stored XX/YY/ZZ views and does
        not read NMRPipe files, calculate ppm axes, transpose arrays, or write
        projection cache entries.
        """
        data = numpy.asarray(data)
        if data.ndim != 2:
            raise ValueError("Projection data must be 2D")

        size0, size1 = data.shape
        uc0 = ng.pipe.make_uc(dic, data, dim=0)
        uc1 = ng.pipe.make_uc(dic, data, dim=1)
        index0 = numpy.asarray([uc0.ppm(i) for i in range(size0)])
        index1 = numpy.asarray([uc1.ppm(i) for i in range(size1)])

        xx_n, yy_n = numpy.meshgrid(index1, index0)
        xx_y, yy_y = numpy.meshgrid(index0, index1)
        file_labb = tuple(file_labb)

        common = dict(
            dic=dic,
            data=data,
            index0=index0,
            index1=index1,
            x_axis=index1,
            y_axis=index0,
            uc0=uc0,
            uc1=uc1,
            source=source,
            file_labb=file_labb,
            decon=bool(decon),
        )

        # Preserve the legacy label/transpose contract while keeping the
        # source array canonical.  A transpose creates a display view only;
        # the canonical data object is never replaced or re-cached.
        out = {}
        for a, b in (file_labb, file_labb[::-1]):
            out[(a, b, 'n')] = dict(
                common,
                XX=xx_n,
                YY=yy_n,
                ZZ=data,
                labb=(a, b),
                transpose='n',
            )
            out[(a, b, 'y')] = dict(
                common,
                XX=xx_y,
                YY=yy_y,
                ZZ=data.T,
                labb=(b, a),
                transpose='y',
            )
        return out

    def _spectrum_view_payload(self, dic, data, source='', labb=None, transpose='n'):
        """Build a plotting-ready 2D spectrum view in the controller layer."""
        if not isinstance(data, numpy.ndarray):
            try:
                data = data[:]
            except Exception:
                pass
        data = numpy.asarray(data)
        if data.ndim != 2:
            raise ValueError("Spectrum view data must be 2D")
        uc0 = ng.pipe.make_uc(dic, data, dim=0)
        uc1 = ng.pipe.make_uc(dic, data, dim=1)
        y_axis = numpy.asarray([uc0.ppm(i) for i in range(data.shape[0])])
        x_axis = numpy.asarray([uc1.ppm(i) for i in range(data.shape[1])])
        # Preserve the legacy convenience attribute without making a viewer own it.
        uc0.ppms_scale = y_axis
        uc1.ppms_scale = x_axis
        transpose = 'y' if str(transpose).lower() == 'y' else 'n'
        if transpose == 'y':
            XX, YY = numpy.meshgrid(y_axis, x_axis)
            ZZ = data.T
            display_labb = tuple(reversed(tuple(labb or ()))) if labb else None
        else:
            XX, YY = numpy.meshgrid(x_axis, y_axis)
            ZZ = data
            display_labb = tuple(labb or ()) or None
        return dict(
            dic=dic, data=data, ZZ=ZZ, XX=XX, YY=YY,
            x_axis=x_axis, y_axis=y_axis, uc0=uc0, uc1=uc1,
            source=source, labb=display_labb, transpose=transpose,
        )

    def _cache_spectrum_views(self, key='raw'):
        """Prepare true-2D raw/decon plotting views once in the shared store."""
        if getattr(self, 'store', None) is None:
            return False
        payload = self.store.spectra.get(key) or {}
        data = payload.get('data')
        dic = payload.get('dic')
        if data is None or dic is None or getattr(data, 'ndim', len(getattr(data, 'shape', ()))) != 2:
            return False
        labels = tuple(payload.get('labb') or getattr(self, 'labb', ()) or ())
        # NMR display convention is direct dimension on X, indirect on Y.
        display_labels = (labels[1], labels[0]) if len(labels) >= 2 else labels
        source = payload.get('spectrumfile') or payload.get('source') or ''
        for transpose in ('n', 'y'):
            view = self._spectrum_view_payload(dic, data, source=source, labb=display_labels, transpose=transpose)
            self.store.save_view(('spectrum', key, transpose), **view)
        return True

    def get_pseudo3d_view(self, key='raw'):
        """Return the canonical view for a 2D spectrum plus one real pseudoaxis.

        The stored cube is exposed in logical ``[pseudo, y, x]`` order without
        copying its data.  Physical-axis metadata remains in the payload so
        viewers never need to assume that the pseudoaxis is NMRPipe axis 0.
        """
        store = getattr(self, 'store', None)
        if store is None:
            return None

        # A true 2D spectrum is exposed to the fitting workspace as a logical
        # one-plane pseudo3D cube.  This is deliberately a *view-layer*
        # convention only: the physical spectrum remains 2D and FUDA receives
        # the original .ft2 file with ZCOOR=2D.
        topology = self._active_topology()
        data_ndim = getattr(getattr(self, 'data', None), 'ndim',
                            len(getattr(getattr(self, 'data', None), 'shape', ())))
        is_true_2d = (
            topology.spectral_dim_count == 2 and
            not topology.has_pseudo_axis and
            data_ndim == topology.physical_dim_count == 2
        )
        if not is_true_2d and not self._is_pseudo3d_topology():
            return None

        view_key = ('pseudo3d', key)
        cached = store.get_view(view_key)
        if cached is not None:
            return cached

        spectrum = store.spectra.get(key) or {}
        # Raw data may use the DataStore's authoritative main-spectrum fallback.
        # Decon data are optional and must never fall back to raw: doing that
        # makes the raw 2D plane get drawn again as the green calculated overlay.
        if key == 'raw':
            data = spectrum.get('data', getattr(store, 'data', None))
        else:
            data = spectrum.get('data')
        if data is None:
            return None
        data = numpy.asarray(data)

        labels = tuple(spectrum.get('labb') or getattr(store, 'labb', None)
                       or getattr(self, 'labb', ()) or ())

        if is_true_2d:
            if data.ndim != 2 or len(labels) < 2:
                return None
            # NMRPipe/nmrglue array order is [indirect, direct] == [y, x].
            # expand_dims returns a view where possible, so the fitting panel
            # does not own a duplicate scientific data set.
            y_uc = getattr(store, 'uc0', None) or getattr(self, 'uc0', None)
            x_uc = getattr(store, 'uc1', None) or getattr(self, 'uc1', None)
            if y_uc is None or x_uc is None:
                try:
                    dic = spectrum.get('dic')
                    y_uc = ng.pipe.make_uc(dic, data, dim=0)
                    x_uc = ng.pipe.make_uc(dic, data, dim=1)
                except Exception:
                    return None
            y_axis = numpy.asarray([y_uc.ppm(i) for i in range(data.shape[0])])
            x_axis = numpy.asarray([x_uc.ppm(i) for i in range(data.shape[1])])
            y_uc.ppms_scale = y_axis
            x_uc.ppms_scale = x_axis
            cube = numpy.expand_dims(data, axis=0)
            XX, YY = numpy.meshgrid(x_axis, y_axis)
            payload = dict(
                data=cube, source_data=data,
                pseudo_dim=None, y_dim=0, x_dim=1,
                pseudo_label='2D', y_label=labels[0], x_label=labels[1],
                labb=('2D', labels[0], labels[1]),
                pseudo_uc=None, y_uc=y_uc, x_uc=x_uc,
                pseudo_axis=numpy.asarray([1.0]), y_axis=y_axis, x_axis=x_axis,
                XX=XX, YY=YY,
                source=spectrum.get('spectrumfile') or spectrum.get('source') or '',
                physical_dim=2, is_single_plane=True,
            )
            store.save_view(view_key, **payload)
            return store.get_view(view_key)

        if data.ndim != 3 or len(labels) < 3:
            return None

        # Stage 7: physical axis identity comes from DatasetTopology.
        # Do not independently rediscover the pseudo axis from labels here:
        # doing so can disagree with topology and produce repeated moveaxis
        # source axes.
        try:
            topology = self.state.topology()
            topology.validate_data_ndim(data.ndim)
            pseudo_axis_spec = topology.pseudo_axis
            spectral_axes = self._pseudo3d_spectral_axes()
        except (AttributeError, ValueError):
            return None
        if pseudo_axis_spec is None or len(spectral_axes) != 2:
            return None

        pseudo_dim = pseudo_axis_spec.physical_index
        # Preserve physical spectral-axis order: the first remaining axis is Y
        # and the second is X.  This reproduces the legacy [pseudo, y, x]
        # convention when the pseudoaxis is physical axis 0.
        y_dim, y_label = spectral_axes[0]
        x_dim, x_label = spectral_axes[1]
        if len({pseudo_dim, y_dim, x_dim}) != 3:
            raise ValueError(
                'pseudo3D topology contains repeated physical axis indices: '
                f'pseudo={pseudo_dim}, y={y_dim}, x={x_dim}'
            )

        pseudo_uc = getattr(store, 'uc%d' % pseudo_dim, None)
        y_uc = getattr(store, 'uc%d' % y_dim, None)
        x_uc = getattr(store, 'uc%d' % x_dim, None)
        if pseudo_uc is None or y_uc is None or x_uc is None:
            return None

        def axis_values(axis_dim, uc):
            scale = getattr(uc, 'ppms_scale', None)
            if scale is None:
                scale = getattr(store, 'index%d' % axis_dim, None)
            if scale is None:
                try:
                    scale = [uc.ppm(i) for i in range(data.shape[axis_dim])]
                except Exception:
                    return None
            values = numpy.asarray(scale)
            if values.size != data.shape[axis_dim]:
                return None
            return values

        pseudo_axis = axis_values(pseudo_dim, pseudo_uc)
        y_axis = axis_values(y_dim, y_uc)
        x_axis = axis_values(x_dim, x_uc)
        if pseudo_axis is None or y_axis is None or x_axis is None:
            return None

        # moveaxis normally returns a view: no second scientific cube is owned
        # by Pseudo3D or by the cache.
        cube = numpy.moveaxis(data, (pseudo_dim, y_dim, x_dim), (0, 1, 2))
        XX, YY = numpy.meshgrid(x_axis, y_axis)
        payload = dict(
            data=cube, source_data=data,
            pseudo_dim=pseudo_dim, y_dim=y_dim, x_dim=x_dim,
            pseudo_label=labels[pseudo_dim], y_label=y_label, x_label=x_label,
            labb=(labels[pseudo_dim], y_label, x_label),
            pseudo_uc=pseudo_uc, y_uc=y_uc, x_uc=x_uc,
            pseudo_axis=pseudo_axis, y_axis=y_axis, x_axis=x_axis,
            XX=XX, YY=YY,
            source=spectrum.get('spectrumfile') or spectrum.get('source') or '',
            physical_dim=3, is_single_plane=False,
        )
        store.save_view(view_key, **payload)
        return store.get_view(view_key)

    def get_spectrum_view(self, decon=False, transpose='n'):
        """Return a shared plotting-ready true-2D spectrum view."""
        if getattr(self, 'store', None) is None:
            return None
        key = 'decon' if decon else 'raw'
        transpose = 'y' if str(transpose).lower() == 'y' else 'n'
        view_key = ('spectrum', key, transpose)
        view = self.store.get_view(view_key)
        if view is None:
            self._cache_spectrum_views(key)
            view = self.store.get_view(view_key)
        return view

    def cache_external_2d_view(self, path, namespace='overlay'):
        """Read an explicitly selected 2D spectrum once and store its views centrally."""
        if not path or not os.path.exists(path) or getattr(self, 'store', None) is None:
            return None
        dic, data = ng.pipe.read(path)
        data = numpy.asarray(data)
        if data.ndim != 2:
            raise ValueError('Selected overlay must be a 2D spectrum')
        for transpose in ('n', 'y'):
            view = self._spectrum_view_payload(dic, data, source=path, transpose=transpose)
            self.store.save_view((namespace, path, transpose), **view)
        return self.store.get_view((namespace, path, 'n'))

    def _cache_projection_array(self, dic, data, file_labb, source, decon=False):
        """Cache canonical projection data and all display views centrally."""
        store = getattr(self, "store", None)
        if store is None:
            return False
        views = self._projection_view_payloads(dic, data, file_labb, source, decon=decon)
        prefix = 'decon_projection' if decon else None
        for (a, b, transpose), payload in views.items():
            if prefix is None:
                # Keep legacy raw aliases during the wider GUI migration.
                store.save_projection((a, b, transpose), **payload)
                store.save_projection(('projections', a, b, transpose), **payload)
            else:
                store.save_projection((prefix, a, b, transpose), **payload)
        return True

    def _cache_projection_file(self, infile, file_labb=None, decon=False):
        """Read one projection file and cache its canonical data/views."""
        if not infile or not os.path.exists(infile):
            return False
        try:
            dic, data = ng.pipe.read(infile)
        except Exception as exc:
            print("Failed to cache projection", infile, exc)
            return False
        if file_labb is None:
            filename = os.path.basename(infile)
            stem = filename[:-4] if filename.endswith('.dat') else filename
            if stem.endswith('.decon'):
                stem = stem[:-6]
            parts = stem.split('.')
            if len(parts) < 2:
                return False
            file_labb = (parts[0], parts[1])
        try:
            return self._cache_projection_array(dic, data, file_labb, infile, decon=decon)
        except Exception as exc:
            print("Failed to prepare projection views", infile, exc)
            return False

    def _cache_projection_folder(self, folder, key_prefix=None):
        """Read raw/decon projections and prepare plotting views in DataStore."""
        store = getattr(self, "store", None)
        if store is None or not os.path.isdir(folder):
            return False
        want_decon = key_prefix == "decon_projection"
        cached_any = False
        for filename in sorted(os.listdir(folder)):
            if not filename.endswith(".dat"):
                continue
            is_decon = filename.endswith(".decon.dat")
            if want_decon != is_decon:
                continue
            infile = os.path.join(folder, filename)
            cached_any = self._cache_projection_file(infile, decon=want_decon) or cached_any
        return cached_any

    def get_projection_view(self, left, right, decon=False, transpose='n'):
        """Return a centrally prepared projection plotting payload.

        The GUI never reads files or derives axes.  If an existing projection
        has not yet been cached (for example when opening an older project),
        the controller resolves and caches it here before returning the view.
        """
        left = str(left or '').strip()
        right = str(right or '').strip()
        transpose = 'y' if str(transpose).lower() == 'y' else 'n'
        if not left or not right or getattr(self, 'store', None) is None:
            return None
        key = (('decon_projection', left, right, transpose)
               if decon else (left, right, transpose))
        payload = self.store.projections.get(key)
        if payload is not None and all(name in payload for name in ('XX', 'YY', 'ZZ', 'labb')):
            return payload

        infile = self._resolve_projection_input_path(left, right, decon=decon, transpose=transpose)
        if infile:
            self._cache_projection_file(infile, decon=decon)
            payload = self.store.projections.get(key)
        return payload

    def _cache_decon_projection_folders(self, base_dir):
        """Cache decon projection folders from the current dataset root."""
        candidates = [
            os.path.join(base_dir, "raw", "projection_decon"),
            os.path.join(base_dir, "projection_decon"),
            os.path.join(base_dir, "raw", "projections_decon"),
            os.path.join(base_dir, "out", "projection_decon"),
        ]
        cached_any = False
        for folder in candidates:
            cached_any = self._cache_projection_folder(folder, key_prefix="decon_projection") or cached_any
        return cached_any

    def _cache_3d_decon_projection_files(self, base_dir):
        """Cache 3D decon projections from out/*.decon into the shared store.

        The 3D projection window expects the decon overlays for the three views
        that correspond to the raw projections:
        - yz -> dims (2, 1)
        - xz -> dims (2, 0)
        - xy -> dims (1, 0)
        """
        store = getattr(self, "store", None)
        if store is None:
            return False

        if not hasattr(self, "labb") or len(getattr(self, "labb", [])) < 3:
            return False

        mappings = [
            ("yz", (self.labb[2], self.labb[1], "n")),
            ("xz", (self.labb[2], self.labb[0], "n")),
            ("xy", (self.labb[1], self.labb[0], "y")),
        ]
        cached_any = False
        for plane, (a, b, transpose) in mappings:
            infile = os.path.join(base_dir, "out", f"{plane}.decon")
            if not os.path.exists(infile):
                continue
            try:
                dic, data = ng.pipe.read(infile)
            except Exception as exc:
                print("Failed to cache 3D decon projection", infile, exc)
                continue
            self._cache_projection_array(dic, data, (a, b), infile, decon=True)
            # Also store the canonical plane name for callers that use the
            # generic plane cache directly.
            store.save_projection(("decon_plane", plane), dic=dic, data=data, source=infile, decon=True)
            cached_any = True
        return cached_any

    def _cache_decon_plane_file(self, plane, infile):
        store = getattr(self, "store", None)
        if store is None or not os.path.exists(infile):
            return False
        Xs, Ys, Zs = self.GetData(infile)
        store.save_projection(("decon_plane", plane), Xs=Xs, Ys=Ys, Zs=Zs, source=infile, decon=True)
        return True

    def _decon_shape_matches_main_spectrum(self, decon_data, decon_path=''):
        """Return True only when a deconvolved array matches the loaded spectrum shape."""
        if not hasattr(self, 'data'):
            print('Cannot verify deconvolved spectrum shape: main spectrum is not loaded.')
            return False
        main_shape = tuple(numpy.shape(self.data))
        decon_shape = tuple(numpy.shape(decon_data))
        if decon_shape != main_shape:
            print('Deconvolved spectrum is a different shape from the main spectrum.')
            if decon_path:
                print('Deconvolution file:', decon_path)
            print('Deconvolved shape:', decon_shape)
            print('Main spectrum shape:', main_shape)
            print('Recalculate the deconvolution for the currently loaded spectrum.')
            return False
        return True

    def _load_decon_outputs(self, infile, load_peaks=True):
        decon_path = infile + '.decon'
        if not os.path.exists(decon_path):
            return False

        try:
            self.dicdec, self.datadec = ng.pipe.read(decon_path)
        except Exception as exc:
            print("Failed to read decon output:", decon_path, exc)
            return False

        # Validate every dimensionality before marking the deconvolution as
        # loaded or publishing it to the shared store.  This prevents a stale
        # .decon file from a different spectrum being analysed accidentally.
        if not self._decon_shape_matches_main_spectrum(self.datadec, decon_path):
            self.DECON = 0
            return False

        self.DECON = 1
        topology = self._active_topology()
        spectral_dim_count = topology.spectral_dim_count
        if getattr(self, "store", None) is not None:
            self.store.save_spectrum("decon", dic=self.dicdec, data=self.datadec, spectrumfile=decon_path, dim=spectral_dim_count, labb=getattr(self, "labb", None))
            self._cache_spectrum_views("decon")

        if spectral_dim_count == 1:
            return True

        if spectral_dim_count == 3 and not topology.has_pseudo_axis:
            base_dir = os.path.dirname(os.path.abspath(infile))
            # Current format: raw and calculated projections share the same
            # projections directory but use distinct filenames.  Keep the old
            # readers as fallbacks for historical datasets.
            self._cache_projection_folder(os.path.join(base_dir, "projections"), key_prefix="decon_projection")
            self._cache_decon_projection_folders(base_dir)
            self._cache_3d_decon_projection_files(base_dir)
            self.pkSlice1Ddec = []
            for pkl in range(len(self.peak)):
                ptC = self.pkIdx[pkl][0]
                ptH = self.pkIdx[pkl][1]
                self.pkSlice1Ddec.append(self.datadec[:, ptC, ptH])
            if load_peaks and getattr(self, "store", None) is not None:
                self.store.save_peak_list("decon", peak=self.peak, pkIdx=self.pkIdx, pkSlice1Ddec=self.pkSlice1Ddec)
            return True

        if spectral_dim_count == 4:
            if getattr(self, "store", None) is not None:
                self.store.save_peak_list("decon", peak=self.peak, pkIdx=self.pkIdx, pkSlice1Ddec=self.pkSlice1Ddec)
            for plane in ("za", "xy", "yz", "xz"):
                self._cache_decon_plane_file(plane, os.path.join('out', f'{plane}.decon'))
            return True

        if spectral_dim_count == 2:
            if load_peaks and getattr(self, "store", None) is not None:
                self.store.save_peak_list("decon", peak=self.peak, pkIdx=self.pkIdx, pkSlice1Ddec=self.pkSlice1Ddec)
            return True

        return True

    def __init__(self,parent,deconParFile,state=None,store=None,decon_parent=None):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.build=platform.uname()[0]
        self.deconBin='decon_'+self.build+'_'+os.popen('uname -p').read().rstrip()
        self.paraDeconBin='decon_parallel_'+self.build+'_'+os.popen('uname -p').read().rstrip()

        self.canvas = None
        self.Parse=Parse
        self.ParseFlt=ParseFlt
        self.ParseInt=ParseInt
        self.ParseAllStr=ParseAllStr
        self.pseudo_spectrum = False
        self.has_pipe = self.check_nmrPipe()
        

        self.state = state if state is not None else ProjectState(session_file=deconParFile, parameter_file=deconParFile)
        # Workflow/ProjectState now owns dataset topology.  Preserve the former
        # NMR radio-box default (1 spectral dimension) for new projects.
        if self.state.spectral_dimensions < 1:
            self.state.spectral_dimensions = 1
        self.store = store if store is not None else DataStore()
        self.deconParFile = deconParFile

        # Shared conversion helper instance available throughout the project.
        self.vpar = vpar()
        self.nmrPipe = nmrPipe(self)

        self.dim=0
        self.READ=0  #zero if data not read in, 1 if not
        self.DECON=0 #deconvolution not been performed
        self.PEAK=0
        self.pseudo = False


        self.peak=[]
        self.Grps={}
        self.pkIdx=[] #index of peak 
        self.pkSlice1D=[] #1D slices
        self.pkSlice1Ddec=[] #1D slices

        self.corrFile=''
        self.peak_frame = None

        self.parent=parent

        # Set sizer for the frame, so we can change frame size to match widgets
        self.windowSizer = wx.BoxSizer()
        self.windowSizer.Add(self, 1, wx.ALL | wx.EXPAND)

        self.data_box()
        # Keep the historical controls alive as the canonical NMR-owned values,
        # but do not display the old Data box on the NMR page.
        self.dataLbl.Hide()
        # ProjectState owns project paths; the NMR tab is their sole editor.
        self._apply_state_to_path_controls()
        self.spectrum_box()
        self.peaklist_box()
        self.decon_box(gui_parent=decon_parent)
        self.project_box()
        self.status_box()
        self.projection_box()
        self.noise_box()
        self._install_status_help()








        ################################################
        #Assemble panel
        # self.splitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.panel1=wx.BoxSizer(wx.VERTICAL)
        self.panel1.Add(self.projectSizer, 0, wx.EXPAND)
        self.panel1.AddSpacer(10)
        self.panel1.Add(self.spectrumSizer, 0, wx.EXPAND)
        self.panel1.AddSpacer(10)

        self.panel1.Add(self.peakSizer, 0, wx.EXPAND)

        # Projection and Noise remain the middle column.  Use a compact
        # horizontal box layout rather than a GridBagSizer with growable spacer
        # columns: the latter distributed a surprising amount of width into
        # empty cells and made the three visible columns look misaligned.
        self.panel2 = wx.BoxSizer(wx.VERTICAL)
        self.panel2.Add(self.projectionSizer, 1, flag=wx.EXPAND)
        self.panel2.AddSpacer(8)
        self.panel2.Add(self.noiseSizer, 1, flag=wx.EXPAND)

        self.splitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.splitSizer.Add(self.panel1, 0, wx.EXPAND)
        self.splitSizer.AddSpacer(8)
        # Give plots the spare width: this is the part of the NMR page that
        # benefits most from resizing, while the control/report columns only
        # need enough width to keep their contents readable.
        self.splitSizer.Add(self.panel2, 1, wx.EXPAND)

        # Deconvolution controls may be hosted by the separate UniDec notebook
        # page. Keep the controls logically owned by this frame, but preserve
        # their historical position between plots and Report when shown here.
        if decon_parent is None:
            self.splitSizer.AddSpacer(8)
            self.splitSizer.Add(self.deconSizer, 0, wx.EXPAND)

        self.splitSizer.AddSpacer(8)
        self.splitSizer.Add(self.statusSizer, 0, wx.EXPAND)

        self.fullSizer = wx.BoxSizer(wx.VERTICAL)
        self.fullSizer.Add(self.splitSizer, 1, wx.EXPAND | wx.ALL, 8)


        self.SetSizerAndFit(self.fullSizer)

        #FGA changed- make decon, sparky and nmr sizers the same width
        dataWidth = 0  # Data controls are retained as hidden canonical owners.
        spectrumWidth = self.spectrumSizer.GetSize()[0]
        peakWidth = self.peakSizer.GetSize()[0]
        projectWidth = self.projectSizer.GetSize()[0]
        deconWidth = self.deconSizer.GetSize()[0]
        widthToSetDecn = max(dataWidth, spectrumWidth, peakWidth, projectWidth)

        panelHeight = self.panel1.GetSize()[1]
        deconHeight = self.deconSizer.GetSize()[1]
        statusHeight = self.statusSizer.GetSize()[1]
        heightToSetDecn = max(panelHeight,deconHeight,statusHeight)

        self.projectSizer.SetMinSize((widthToSetDecn, 0))
        self.spectrumSizer.SetMinSize((widthToSetDecn, 0))
        self.peakSizer.SetMinSize((widthToSetDecn, 0))

        self.border2.SetMinSize((widthToSetDecn, 0))
        # Report does not need to inherit the much wider file-control column.
        # A modest minimum keeps every status line visible without reserving
        # unused horizontal space.
        self.statusSizer.SetMinSize((230, heightToSetDecn))

        self.panel1.SetMinSize((widthToSetDecn, heightToSetDecn))
        self.panel2.SetMinSize((195, heightToSetDecn))
        self.deconSizer.SetMinSize((0,heightToSetDecn))

        #self.sparkSizer.SetMinSize((widthToSetDecn, 0))
        #self.nmrSizer.SetMinSize((widthToSetDecn, 0))
        #self.SetSizerAndFit(self.splitSizer)


        #Set the sizers
        self.SetSizerAndFit(self.fullSizer)


        if self.has_pipe == False:
            self.buttonProcess.Disable()
            # ReOrganise remains available so the command window can explain
            # availability through its disabled dimension-dependent buttons.

        # Set event handlers
        self.buttonDecon.Bind(wx.EVT_BUTTON, self.OnButtonDecon)
        self.buttonRecon.Bind(wx.EVT_BUTTON, self.OnButtonRecon)
        self.buttonPeakFit.Bind(wx.EVT_BUTTON, self.OnButtonPeakFit)
        self.buttonReOrganise.Bind(wx.EVT_BUTTON, self.OnButtonReOrganise)
        self.buttonRead.Bind(wx.EVT_BUTTON, self.OnButtonRead)
        self.buttonProcess.Bind(wx.EVT_BUTTON, self.OnButtonProcess)
        self.buttonPeaky.Bind(wx.EVT_BUTTON, self.OnButtonPeaky)
        self.buttonReadPeak.Bind(wx.EVT_BUTTON, self.OnButtonReadPeak)
        self.buttonReadFullPeak.Bind(wx.EVT_BUTTON, self.OnButtonReadFullPeak)
        self.buttonReferencePeakList.Bind(wx.EVT_BUTTON, self.OnButtonReferencePeakList)
        self.buttonFullPeakList.Bind(wx.EVT_BUTTON, self.OnButtonFullPeakList)
        self.buttonRestricted3DDiagnostics.Bind(wx.EVT_BUTTON, self.OnButtonRestricted3DDiagnostics)

        self.buttonAnalyse.Bind(wx.EVT_BUTTON, self.OnButtonAnalyse)

        ##################

        indir=Parse(self.deconParFile,'indir',default='./')

        #if str(indir) != '':
        #    os.chdir(str(indir))
        #    #indir=path_escape(indir)
        #print(indir)
        thresh=ParseFlt(self.deconParFile,'thresh',default=UNIDEC_DEFAULTS['thresh'])
        ncpus=int(ParseInt(self.deconParFile,'ncpus',default=available_cpu_count()))
        fac=ParseFlt(self.deconParFile,'fac',default=UNIDEC_DEFAULTS['fac'])
        conv=ParseFlt(self.deconParFile,'conv',default=UNIDEC_DEFAULTS['conv'])
        maxiter=int(ParseFlt(self.deconParFile,'maxiter',default=UNIDEC_DEFAULTS['maxiter']))
        fitrad=Parse(self.deconParFile,'FitRad',default='')
        fitf1=Parse(self.deconParFile,'3p_radF1',default='')
        fitf2=Parse(self.deconParFile,'3p_radF2',default='')
        # Processed spectra and their reference peak lists live below SpecPath,
        # never below the vendor/acquisition (raw) directory.  The canonical
        # processed filename is determined by the number of *spectral*
        # dimensions; a 2D+pseudo dataset therefore uses test.ft2.
        try:
            _default_spectral_dim = int(Parse(self.deconParFile, 'dim', default=2))
        except (TypeError, ValueError):
            _default_spectral_dim = 2
        _default_spectral_dim = max(1, min(4, _default_spectral_dim))
        _default_spectrum_name = ('test.ft' if _default_spectral_dim == 1
                                  else 'test.ft%d' % _default_spectral_dim)
        infile=Parse(self.deconParFile,'infile',default='./spec/' + _default_spectrum_name)
        pseudo=Parse(self.deconParFile,'pseudo')

        peakfile=Parse(self.deconParFile,'peakfile',default='./spec/' + _default_spectrum_name + '.list')
        voigt1=Parse(self.deconParFile,'voigt1')
        voigt2=Parse(self.deconParFile,'voigt2')
        voigt3=Parse(self.deconParFile,'voigt3')


        self.dirBox.SetValue(str(indir))
        # if(str(indir)!='0'):
        #     os.chdir(str(indir_old))
        self.coreBox.SetValue(str(ncpus))
        self.threshBox.SetValue(str(thresh))
        self.facBox.SetValue(str(fac))
        self.convBox.SetValue(str(conv))
        self.maxiterBox.SetValue(str(maxiter))
        self.fitRadBox.SetValue('' if str(fitrad) == '0' else str(fitrad))
        self.fitF1Box.SetValue('' if str(fitf1) in ('0', 'None') else str(fitf1))
        self.fitF2Box.SetValue('' if str(fitf2) in ('0', 'None') else str(fitf2))
        self.infileBox.SetValue(str(infile))
        self.referencePeakBox.SetValue(str(peakfile))

        """
        try:
            self._sync_peakfile_box(self._resolve_input_path(self.infileBox.GetValue()), int(self.state.spectral_dimensions))
        except Exception:
            pass
        """
        
        
        try:
            loaded_pseudo = bool(int(Parse(self.deconParFile, 'pseudo')))
        except Exception:
            loaded_pseudo = bool(Parse(self.deconParFile, 'pseudo'))
        self.state.pseudo_axis = loaded_pseudo
        self.pseudo = loaded_pseudo
            
        if(Parse(self.deconParFile,'symmode')=='y'):
            self.cb_grid.SetValue(1)
        else:
            self.cb_grid.SetValue(0)

        self.cb_enhance.SetValue(Parse(self.deconParFile, 'enhance', default='n') == 'y')
        self.cb_fitphases.SetValue(Parse(self.deconParFile, 'fitPhases', default='n') == 'y')

        self.Status()
        self._restore_workflow_flags_from_parameter_file()

    def _restore_workflow_flags_from_parameter_file(self):
        """Restore persisted guided-workflow evidence from the system file.

        Older system files simply omit these keys and therefore remain fully
        compatible.  Scientific values remain in their historical fields; the
        flags only record whether those values were explicitly accepted and
        whether dependent pseudo3D fits must be regenerated.
        """
        fitted = str(Parse(self.deconParFile, 'peakShapeFitted', default='0') or '0').strip().lower()
        stale = str(Parse(self.deconParFile, 'pseudoIntensitiesStale', default='0') or '0').strip().lower()
        downstream = str(Parse(self.deconParFile, 'downstreamAnalysis', default='') or '').strip()
        inspected_raw = Parse(self.deconParFile, 'pseudoSeriesInspected', default='0')
        inspected = str(inspected_raw or '0').strip().lower()
        picked_reviewed_raw = Parse(self.deconParFile, 'pickedPeaksChecked', default='0')
        picked_reviewed = str(picked_reviewed_raw or '0').strip().lower()
        fitting_reviewed_raw = Parse(self.deconParFile, 'fittingResultsInspected', default='0')
        fitting_reviewed = str(fitting_reviewed_raw or '0').strip().lower()
        peak_pick_stale_raw = Parse(self.deconParFile, 'peakPickStale', default='0')
        peak_pick_stale = str(peak_pick_stale_raw or '0').strip().lower()
        pass
        fitted = fitted in ('1', 'y', 'yes', 'true')
        stale = stale in ('1', 'y', 'yes', 'true')
        inspected = inspected in ('1', 'y', 'yes', 'true')
        picked_reviewed = picked_reviewed in ('1', 'y', 'yes', 'true')
        fitting_reviewed = fitting_reviewed in ('1', 'y', 'yes', 'true')
        peak_pick_stale = peak_pick_stale in ('1', 'y', 'yes', 'true')
        self.peak_shape_fitted = fitted
        self.pseudo_intensities_stale = stale
        self.peak_pick_stale = peak_pick_stale
        self.downstream_analysis = downstream
        store = getattr(self, 'store', None)
        if store is not None:
            if fitted:
                store.mark_peak_shape_determined(source='system_file')
            if inspected:
                store.mark_pseudo_series_reviewed(source='system_file')
            else:
                store.invalidate_pseudo_series_review()
            if picked_reviewed:
                store.mark_picked_peaks_reviewed(source='system_file')
            else:
                store.invalidate_picked_peaks_review()
            if fitting_reviewed:
                store.analysis['fitting_results_ready'] = True
                store.mark_fitting_results_reviewed(source='system_file')
            else:
                store.invalidate_fitting_results_review()
            if peak_pick_stale:
                store.analysis['peak_pick_stale'] = True
            else:
                store.analysis.pop('peak_pick_stale', None)
            if downstream:
                store.analysis['downstream_analysis'] = downstream
            else:
                store.analysis.pop('downstream_analysis', None)
            if stale:
                store.analysis.pop('pseudo_intensities_ready', None)
                store.analysis.pop('pseudo_series_reviewed', None)
                store.analysis.pop('pseudo_analysis_complete', None)
            # Do not inspect fit/ synchronously here.  This routine runs while
            # the system file is still being restored; at that point the main
            # spectrum/reference-list controls can be populated while READ and
            # the shared peak store are not yet ready.  A synchronous check can
            # therefore report an empty reference list even though a complete
            # set of Protocol3P files exists.  Reconcile disk evidence after
            # project restoration has returned to the wx event loop.
            if fitted and not stale and bool(self.state.pseudo_axis):
                if not getattr(self, '_pseudo_fit_reconcile_pending', False):
                    self._pseudo_fit_reconcile_pending = True
                    wx.CallAfter(self._reconcile_existing_pseudo3d_fit_evidence)
        pass
        self._workflow_debug('restored workflow flags: peakShapeFitted=%s pseudoIntensitiesStale=%s downstreamAnalysis=%r' % (fitted, stale, downstream))

    def _reconcile_existing_pseudo3d_fit_evidence(self):
        """Publish complete on-disk Protocol3P results after a cold load.

        The reference peak list is authoritative: reuse fit/ only when every
        named reference peak has both ``<name>.dat`` and ``<name>.out``.  This
        is deliberately deferred until project/system-file restoration has
        completed so the configured spectrum and reference-list paths are
        stable.  It does not run reconstruction and never overrides the stale
        flag set by a revised peak shape.
        """
        self._pseudo_fit_reconcile_pending = False
        if not bool(getattr(self, 'peak_shape_fitted', False)):
            return False
        if bool(getattr(self, 'pseudo_intensities_stale', False)):
            return False
        if not bool(self.state.pseudo_axis):
            return False
        store = getattr(self, 'store', None)
        if store is None:
            return False
        try:
            # A COMPLETE Establish reference peaks stage means more than the
            # list file existing on disk: the configured list must have been
            # loaded into the shared DataStore before Workflow is allowed to
            # reason about extraction.  On a cold start OnButtonReadPeak()
            # expects the main spectrum to have been initialised first, so do
            # that in the same stage order used interactively:
            # spectrum -> peak shape -> reference list -> extraction evidence.
            if not self.ensure_workflow_reference_stage_loaded():
                self._workflow_debug('cold-start fit evidence: reference stage could not be materialised yet')
                return False
            refs = list(self.get_reference_peaks() or [])
            missing = self.missing_pseudo3d_fit_peaks()
            self._workflow_debug('cold-start fit evidence: reference_count=%d missing=%r fit_dir=%r' %
                                 (len(refs), missing, self.get_fuda_dir()))
            if not refs or missing:
                return False
            store.mark_pseudo_intensities_ready(
                source='existing_protocol3p_fit_cold_start',
                invalidate_review=False,
                fit_directory=self.get_fuda_dir(),
                reference_peak_count=len(refs))
            self._notify_analysis_changed()
            return True
        except Exception as exc:
            self._workflow_debug('cold-start fit evidence check skipped: %r' % (exc,))
            return False

    def peak_shape_saved(self, was_already_fitted=False):
        """Persist peak-shape completion and invalidate dependent pseudo3D fits.

        Re-saving an already fitted shape is treated as a revision: reference
        peaks are reloaded and Protocol3P is forced to run even if fit/ already
        contains a complete set of files.
        """
        self.peak_shape_fitted = True
        if was_already_fitted:
            self.peak_pick_stale = True
            store = getattr(self, 'store', None)
            if store is not None:
                store.analysis['peak_pick_stale'] = True
                store.invalidate_picked_peaks_review()
        if was_already_fitted and self._is_pseudo3d_topology():
            self.pseudo_intensities_stale = True
            store = getattr(self, 'store', None)
            if store is not None:
                store.analysis.pop('pseudo_intensities_ready', None)
                store.analysis.pop('pseudo_intensities', None)
                store.analysis.pop('pseudo_series_reviewed', None)
                store.analysis.pop('pseudo_series_review', None)
                store.analysis.pop('pseudo_analysis_complete', None)
                store.analysis.pop('pseudo_analysis', None)
        self.OnButtonSave(True)
        self._notify_analysis_changed()
        if was_already_fitted and self._is_pseudo3d_topology():
            # A peak-shape/threshold edit invalidates downstream pseudo-axis
            # intensities, but it must not *run* the workflow.  In particular,
            # reference-peak selection is itself an interactive setup step and
            # users must remain free to revise threshold/shape before deciding
            # which peak list is authoritative.  The Workflow page will expose
            # Extract intensities as READY/BLOCKED and enforce its prerequisites
            # only when that action is explicitly requested.
            self._workflow_debug('peak shape revised: pseudo3D intensities marked stale; extraction deferred to workflow')

    def _mark_pseudo3d_recompute_complete(self):
        """Clear the persisted invalidation flag after successful Protocol3P."""
        if not getattr(self, 'pseudo_intensities_stale', False):
            return
        self.pseudo_intensities_stale = False
        # Persist only the workflow flags plus the current normal settings.
        self.OnButtonSave(True)
        self._workflow_debug('pseudo3D recomputation complete; stale flag cleared')

    _REAL_AXIS_LABELS = frozenset((
        'time_T2', 'ID', 'ncyc', 'ncyc_cp', 'gzlvl5', 'gzlvl1'
    ))

    def _is_pseudo3d_topology(self):
        """True exactly for two spectral axes plus one real pseudo axis.

        Stage 7: this is a topology query, not a legacy ``dim`` detector.
        Physical dimensionality and pseudo-axis position come from the
        canonical ProjectState/DatasetTopology contract.
        """
        try:
            topology = self.state.topology()
        except Exception:
            from spinDecon.domain.topology import DatasetTopology
            topology = DatasetTopology.from_counts(
                max(1, int(getattr(self, 'dim', 1))),
                bool(self.state.pseudo_axis),
            )
        return (topology.spectral_dim_count == 2 and
                topology.has_pseudo_axis and
                topology.physical_dim_count == 3)

    def _pseudo3d_spectral_axes(self):
        """Return physical (axis_index, label) pairs for pseudo-3D spectral axes."""
        if not self._is_pseudo3d_topology():
            return []
        try:
            topology = self.state.topology()
            result = []
            for axis in topology.spectral_axes:
                label = axis.label
                if not label and axis.physical_index < len(getattr(self, 'labb', ())):
                    label = self.labb[axis.physical_index]
                result.append((axis.physical_index, label))
            return result
        except Exception:
            # Compatibility only for an incompletely initialised frame.  The
            # scientific dimensionality remains spectral; labels are not used
            # to reinterpret the dimension count.
            return [(i, label) for i, label in enumerate(self.labb[:3])
                    if str(label) not in self._REAL_AXIS_LABELS]

    def _draw_main_1d_projection(self, axis_index):
        """Draw the loaded 1D projection against the requested physical spectral axis."""
        projected = getattr(self, 'projectedData', None)
        if projected is None or numpy.size(projected) == 0:
            return
        uc = getattr(self, 'uc%d' % axis_index, None)
        if uc is None:
            return
        scale = uc.ppms_scale
        if len(scale) != len(projected):
            print('1D projection length does not match %s axis; not drawing projection.' % self.labb[axis_index])
            return
        self.axes.clear()
        self.axes.get_xaxis().set_visible(True)
        self.axes.tick_params(length=1.0)
        self.axes.plot(scale, projected, color='r', lw=0.5)
        xmin, xmax = scale[0], scale[-1]
        ymin, ymax = self.axes.get_yaxis().get_view_interval()
        self.axes.set_xlim(xmin, xmax)
        self.axes.set_xlabel(str(self.labb[axis_index]) + ' (ppm)',
                             fontsize=self.axes.get_xticklabels()[0].get_fontsize()
                             if self.axes.get_xticklabels() else None)
        self.axes.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
        self.fig.subplots_adjust(left=0.05, right=0.95, bottom=0.2, top=0.95)
        self.canvas.draw()

    def _draw_main_pseudo2d_projection(self):
        """Render the canonical pseudo2D 1D projection in the NMR tab."""
        try:
            payload = self.get_pseudo2d_projection_data(ensure_file=False)
        except Exception as exc:
            print('Unable to prepare pseudo2D main projection:', exc)
            return
        if not payload:
            return
        scale = numpy.asarray(payload.get('index', []), dtype=float).squeeze()
        projected = numpy.asarray(payload.get('data', []), dtype=float).squeeze()
        if scale.ndim != 1 or projected.ndim != 1 or scale.size != projected.size or not scale.size:
            return
        self.axes.clear()
        self.axes.get_yaxis().set_visible(False)
        self.axes.get_xaxis().set_visible(True)
        self.axes.set_frame_on(False)
        self.axes.tick_params(length=1.0)
        self.axes.plot(scale, projected, color='r', lw=0.5)
        # NMR convention is high ppm on the left, low ppm on the right.
        self.axes.set_xlim(float(numpy.nanmax(scale)), float(numpy.nanmin(scale)))
        self.axes.set_xlabel(str(payload.get('label') or 'Direct') + ' (ppm)',
                             fontsize=self.axes.get_xticklabels()[0].get_fontsize()
                             if self.axes.get_xticklabels() else None)
        self.fig.subplots_adjust(left=0.05, right=0.95, bottom=0.2, top=0.95)
        self.canvas.draw()

    ###################################
    #Write status menu
    def Status(self):
        #for i in range(10):
        #    self.updateList[i].SetLabel("shit")
        cnt=0
        if(self.READ==0):
            self.updateList[cnt].SetLabel("Folder name: "+os.getcwd().split('/')[-1]);cnt+=1
            self.updateList[cnt].SetLabel("No spectrum in memory")
        else:
            #self.updateList[cnt].SetLabel("Read in %s" % (self.spectrumFile))
            self.updateList[cnt].SetLabel("Folder name: "+os.getcwd().split('/')[-1]);cnt+=1
            self.updateList[cnt].SetLabel("Spectrum Dimensions: %s" % self.dim);cnt+=1

            # Canonical 3p files are physical 3D arrays containing two
            # frequency axes and one real/pseudo axis.  Report each physical
            # axis from its loaded header metadata rather than applying the
            # older assumption that pseudo mode always adds a hidden dimension.
            _real_axis_labels = self._REAL_AXIS_LABELS
            _physical_3p = self._is_pseudo3d_topology()
            if _physical_3p:
                for _axis in range(3):
                    _label = self.labb[_axis]
                    _size = self.specsize[_axis]
                    if str(_label) in _real_axis_labels:
                        _line = "dim%i: %s real (%i pts)" % (_axis + 1, _label, _size)
                    else:
                        _umin = getattr(self, 'uc%dmin' % _axis)
                        _umax = getattr(self, 'uc%dmax' % _axis)
                        _line = "dim%i: %s %.2f to %.2f ppm (%i pts)" % (_axis + 1, _label, _umin, _umax, _size)
                    self.updateList[cnt].SetLabel(_line);cnt+=1
                self.updateList[cnt].SetLabel("will project down %s" % next((x for x in self.labb[:3] if str(x) in _real_axis_labels), self.labb[0]));cnt+=1
                _spectral_axes = self._pseudo3d_spectral_axes()
                if len(_spectral_axes) >= 2:
                    self.updateList[cnt].SetLabel("peak list must be %s:%s" %
                                                  (_spectral_axes[-1][1], _spectral_axes[-2][1]));cnt+=1
                    # The stored 1D projection is the direct (last spectral) axis.
                    self._draw_main_1d_projection(_spectral_axes[-1][0])
                return

            # Defensive guard: a failed/incomplete read should not crash the
            # whole project while Status() is updating labels.
            if getattr(self, 'labb', None) is None or getattr(self, 'specsize', None) is None:
                self.updateList[cnt].SetLabel("Spectrum metadata unavailable")
                return

            if bool(self.state.pseudo_axis):
                liney="%s %s real (%i pts)" % ("dim1:",self.labb[0],self.specsize[0])
            else:
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim1:",self.labb[0],self.uc0min,self.uc0max,self.specsize[0])
            self.updateList[cnt].SetLabel(liney);cnt+=1



            if(self.dim>=2): #if 3D or greater, at least one more.
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim2:",self.labb[1],self.uc1min,self.uc1max,self.specsize[1])
                self.updateList[cnt].SetLabel(liney);cnt+=1


            if(self.dim==2):
                self.axes.clear()
                self.axes.get_xaxis().set_visible(True)

                self.axes.tick_params(length=1.0)
                # For 2D data the projection shown here is the direct
                # dimension, which is the X axis of the 2D display.
                self.axes.plot(self.uc1.ppms_scale, self.projectedData, color='r', lw=0.5)
                xmin, xmax = self.uc1.ppms_scale[0], self.uc1.ppms_scale[-1]
                ymin, ymax = self.axes.get_yaxis().get_view_interval()
                self.axes.set_xlim(xmin, xmax)
                self.axes.set_xlabel(self.labb[1]+" (ppm)", fontsize=self.axes.get_xticklabels()[0].get_fontsize() if self.axes.get_xticklabels() else None)
                self.axes.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
                self.fig.subplots_adjust(left=0.05,right=0.95,bottom=0.2,top=0.95)
                self.canvas.draw()

            if bool(self.state.pseudo_axis):
                if(self.dim>=2):  #if 3D, more more spectral axis.
                    liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim3:",self.labb[2],self.uc2min,self.uc2max,self.specsize[2])
                elif(self.dim==1): #if psuedo2D, write out the spectral axis.
                    liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim2:",self.labb[1],self.uc1min,self.uc1max,self.specsize[1])
                self.updateList[cnt].SetLabel(liney);cnt+=1

                
            if(self.dim>=3):
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim3:",self.labb[2],self.uc2min,self.uc2max,self.specsize[2])
                self.updateList[cnt].SetLabel(liney);cnt+=1
                if(self.dim==3):
                    self.updateList[cnt].SetLabel("will project down %s" % self.labb[0]);cnt+=1
                    self.updateList[cnt].SetLabel("peak list must be %s:%s" % (self.labb[2],self.labb[1]));cnt+=1
                # try:
                #     self.canvas
                # except NameError:


                    self.axes.clear()
                    self.axes.get_xaxis().set_visible(True)

                    self.axes.tick_params(length=1.0)
                    self.axes.plot(self.uc2.ppms_scale, self.projectedData, color='r', lw=0.5)
                    xmin, xmax = self.uc2.ppms_scale[0], self.uc2.ppms_scale[-1]
                    ymin, ymax = self.axes.get_yaxis().get_view_interval()
                    self.axes.set_xlim(xmin, xmax)
                    self.axes.set_xlabel(self.labb[2]+" (ppm)", fontsize=self.axes.get_xticklabels()[0].get_fontsize() if self.axes.get_xticklabels() else None)
                    self.axes.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
                    self.fig.subplots_adjust(left=0.05,right=0.95,bottom=0.2,top=0.95)
                    self.canvas.draw()
                    

            if(self.dim>=4):
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim4:",self.labb[3],self.uc3min,self.uc3max,self.specsize[3])
                self.updateList[cnt].SetLabel(liney);cnt+=1
                if(self.dim==4):
                    self.updateList[cnt].SetLabel("will project down %s:%s" % (self.labb[0],self.labb[1]));cnt+=1
                    self.updateList[cnt].SetLabel("peak list must be %s:%s" % (self.labb[2],self.labb[3]));cnt+=1

        # In pseudo2D the Full 1D list is the authoritative (and only) picked
        # peak collection.  Do not report the legacy reference-peak count.
        _spectral_dimensions = int(getattr(self.state, 'spectral_dimensions', 0) or 0)
        _topology = self._active_topology() if _spectral_dimensions > 0 else None
        if _topology is not None and _topology.spectral_dim_count == 1 and _topology.has_pseudo_axis:
            _full = self.get_full_peak_payload() or {}
            _picked = _full.get('records') or _full.get('peaks') or _full.get('rows') or []
            if _picked:
                self.updateList[cnt].SetLabel("ProjectedPeaks: %d" % len(_picked))
            else:
                self.updateList[cnt].SetLabel("No peaks in projection")
        elif len(self.peak) == 0:
            self.updateList[cnt].SetLabel("No peaks in projection")
        else:
            self.updateList[cnt].SetLabel("ProjectedPeaks: "+str(len(self.peak)))
        cnt+=1


        if(self.DECON==0):
            self.updateList[cnt].SetLabel("No deconvolved spectrum in memory");cnt+=1
        else:
            self.updateList[cnt].SetLabel("Deconvolved spectrum loaded");cnt+=1

        _decon_peaks = self.store.peak_lists.get('decon', {}).get('records', []) if getattr(self, 'store', None) is not None else []
        if not _decon_peaks:
            self.updateList[cnt].SetLabel("No deconvolved peaks");cnt+=1
        else:
            self.updateList[cnt].SetLabel("Deconvolved peaks: "+str(len(_decon_peaks)));cnt+=1

        for i in range(len(self.updateList)-cnt):
            self.updateList[i+cnt].SetLabel("")


        # Status updates change labels only.  This panel is a notebook page, so
        # fitting it here can resize the page to its minimum size while sibling
        # tabs are being added/replaced (notably on the first 3D Read Ref).
        # Re-layout the existing hierarchy instead and leave page sizing to the
        # notebook.
        self.sizerStat.Layout()
        self.statusSizer.Layout()
        self.fullSizer.Layout()
        self.Layout()

    def CleanDim(self):
        #try:
        #    self.splitSizer.Remove(2)
        #    pass
        #except:
        #    pass

        try:
            a=len(self.cln)
        except:
            self.cln=[]
        #print 'elements:',len(self.cln)
        clean=0
        for i in range(len(self.cln)):
            try:

                self.cln[i].Destroy()
                clean+=1
            except:
                pass

        if(clean!=len(self.cln)):
            print('WARNING! error in cleaning')
        self.cln=[]


    ######################################
        
    # ------------------------------------------------------------------
    # Status-bar hover help
    # ------------------------------------------------------------------
    def _status_help_texts(self):
        """Return concise hover help for the controls on the NMR/UniDec pages."""
        return {
            # Data controls (currently hidden on the NMR page, retained for compatibility).
            'dirLab': 'Working directory used as the base location for this project.',
            'dirBox': 'Working directory used as the base location for this project.',
            'openDirFileBtn': 'Choose the project working directory.',
            'outPathLab': 'Directory containing or receiving raw spectrum data.',
            'outPathBox': 'Directory containing or receiving raw spectrum data.',
            'openOutPathBtn': 'Choose the raw-data output directory.',
            'specPathLab': 'Directory containing processed spectrum files.',
            'specPathBox': 'Directory containing processed spectrum files.',
            'openSpecPathBtn': 'Choose the processed-spectrum directory.',
            'infileLab': 'nmrPipe spectrum file used as the main input dataset.',
            'infileBox': 'nmrPipe spectrum file used as the main input dataset.',
            'openNMRFileBtn': 'Choose the nmrPipe spectrum file to use.',

            # Project.
            'buttonLoadProject': 'Load all available spectrum, peak-list and deconvolution data for the project.',
            'buttonSummariseProject': 'Open a summary of the current project files, settings and readiness.',
            'buttonHousekeeping': 'Review and clean project files using the housekeeping tools.',

            # Spectrum.
            'dimLab': 'Select the dimensionality of the NMR spectrum.',
            'buttonProcess': 'Open the spectrum-processing workflow for the selected dataset.',
            'buttonRead': 'Load the selected spectrum and its metadata into the project.',
            'buttonReOrganise': 'Reorder or transpose spectrum axes into the required project orientation.',

            # Peak lists.
            'referencePeakLab': 'Path to the 2D reference peak list used for navigation and optional deconvolution guidance.',
            'referencePeakBox': 'Path to the 2D reference peak list used for navigation and optional deconvolution guidance.',
            'openPeakFileBtn': 'Choose the 2D reference peak-list file.',
            'buttonReadPeak': 'Load the 2D reference peak list into the project.',
            'buttonReferencePeakList': 'Open the loaded 2D reference peak list for inspection or editing.',
            'buttonPeaky': 'Open the peak viewer for the loaded reference peaks.',
            'fullPeakLab': 'Path to the full-dimensional peak list matching the spectrum.',
            'fullPeakBox': 'Path to the full-dimensional peak list matching the spectrum.',
            'openFullPeakFileBtn': 'Choose the full-dimensional peak-list file.',
            'buttonReadFullPeak': 'Load the full-dimensional peak list into the project.',
            'buttonFullPeakList': 'Open the loaded full-dimensional peak list for inspection or editing.',
            'buttonRestricted3DDiagnostics': 'Open diagnostics for checking the full-dimensional peak assignments and restrictions.',

            # UniDec/deconvolution.
            'buttonDecon': 'Run the deconvolution using the current UniDec settings.',
            'buttonPeakFit': 'Fit peaks to the current deconvolution result.',
            'buttonAnalyse': 'Load the deconvolution output back into the project for analysis.',
            'coreLab': 'Number of CPU cores used by the deconvolution calculation.',
            'coreBox': 'Number of CPU cores used by the deconvolution calculation.',
            'facLab': 'Scaling factor used by the deconvolution calculation.',
            'facBox': 'Scaling factor used by the deconvolution calculation.',
            'convlab': 'Convergence criterion used to decide when iterative deconvolution has converged.',
            'convBox': 'Convergence criterion used to decide when iterative deconvolution has converged.',
            'maxiterLab': 'Maximum number of deconvolution iterations allowed.',
            'maxiterBox': 'Maximum number of deconvolution iterations allowed.',
            'initlistLab': 'Peak-list file used to initialise the deconvolution when list initialisation is enabled.',
            'initlistBox': 'Peak-list file used to initialise the deconvolution when list initialisation is enabled.',
            'openInitListBtn': 'Choose the peak-list file used to initialise deconvolution.',
            'fitRadLab': 'Fitting radius supplied to the deconvolution peak-fitting stage when Fit is enabled.',
            'fitRadBox': 'Automatic fitting radius supplied to the deconvolution peak-fitting stage when Fit is enabled.',
            'fitF1Lab': 'F1 extraction radius in ppm for 2D and 2D+pseudo-axis fitting.',
            'fitF1Box': 'F1 extraction radius in ppm for 2D and 2D+pseudo-axis fitting.',
            'fitF2Lab': 'F2 extraction radius in ppm for 2D and 2D+pseudo-axis fitting.',
            'fitF2Box': 'F2 extraction radius in ppm for 2D and 2D+pseudo-axis fitting.',
            'speedList': 'Choose a preset balance between deconvolution speed and accuracy.',
            'quickButton': 'Apply the quick deconvolution preset for faster processing.',
            'mediumButton': 'Apply the medium deconvolution preset for balanced speed and accuracy.',
            'accurateButton': 'Apply the accurate deconvolution preset for more thorough processing.',
            'cb_grid_label': 'Impose symmetry constraints during deconvolution, for example for symmetric NOE data.',
            'cb_grid': 'Impose symmetry constraints during deconvolution, for example for symmetric NOE data.',
            'cb_decon3d_label': 'Use the loaded 2D reference peak list to guide higher-dimensional deconvolution.',
            'cb_decon3d': 'Use the loaded 2D reference peak list to guide higher-dimensional deconvolution.',
            'cb_decback_label': 'Enable fitting of the deconvolution result after deconvolution.',
            'cb_decback': 'Enable fitting of the deconvolution result after deconvolution.',
            'cb_initlist_label': 'Initialise deconvolution from the peak list specified by InitList.',
            'cb_initlist': 'Initialise deconvolution from the peak list specified by InitList.',
            'nLab': 'Dimension number for the line-shape parameters shown in each column.',
            'sigText': 'Gaussian line width in ppm for each spectral dimension.',
            'voigtText': 'Voigt mixing fraction from 0 to 1 for each spectral dimension.',
            'lorText': 'Lorentzian line width in ppm for each spectral dimension.',

            # Noise/report/projection areas.
            'threshLab': 'Noise threshold used for threshold display and downstream peak assessment.',
            'threshBox': 'Set the noise threshold; press Enter to apply the new value.',
            'noiseDetailButton': 'Show detailed noise and intensity statistics.',
            'projectionLbl': 'Displays the current spectrum projection for visual inspection.',
            'noiseLbl': 'Displays the noise distribution and the currently selected threshold.',
            'statusLbl': 'Reports the availability and current state of project inputs and results.',
        }

    def _set_hover_status(self, text):
        frame = wx.GetTopLevelParent(self)
        statusbar = getattr(frame, 'statusbar', None)
        if statusbar is None:
            return
        # Mark the text as hover-owned so leaving one widget cannot erase an
        # operational status message written by another part of the program.
        frame._hover_status_previous = statusbar.GetStatusText()
        frame._hover_status_text = str(text)
        statusbar.SetStatusText(str(text))

    def _clear_hover_status(self, text):
        frame = wx.GetTopLevelParent(self)
        statusbar = getattr(frame, 'statusbar', None)
        if statusbar is None:
            return
        owned = getattr(frame, '_hover_status_text', None)
        if owned == str(text) and statusbar.GetStatusText() == str(text):
            statusbar.SetStatusText(getattr(frame, '_hover_status_previous', '') or '')
            frame._hover_status_text = None
            frame._hover_status_previous = ''

    def _register_status_help(self, widget, text):
        """Attach status-bar help to a leaf/control widget once.

        Container windows deliberately do not own hover help.  In wxPython a
        StaticBox (and similar parents) can receive enter events while the
        pointer moves between its children, which lets generic section help
        overwrite the more useful help belonging to the control under the
        mouse.
        """
        if widget is None or not isinstance(widget, wx.Window) or not text:
            return
        if isinstance(widget, (wx.StaticBox, wx.ScrolledWindow)):
            return
        # Plain panels are layout/event containers rather than actionable
        # controls.  FigureCanvasWxAgg and other specialised wx.Window
        # subclasses are not wx.Panel, so interactive canvases remain eligible.
        if isinstance(widget, wx.Panel):
            return
        if getattr(widget, '_decon_status_help_bound', False):
            widget._decon_status_help_text = str(text)
            return
        widget._decon_status_help_bound = True
        widget._decon_status_help_text = str(text)
        widget.Bind(wx.EVT_ENTER_WINDOW,
                    lambda evt, w=widget: (self._set_hover_status(w._decon_status_help_text), evt.Skip()))
        widget.Bind(wx.EVT_LEAVE_WINDOW,
                    lambda evt, w=widget: (self._clear_hover_status(w._decon_status_help_text), evt.Skip()))

    def _install_status_help(self):
        """Register explicit semantic help for all persistent NMR/UniDec controls."""
        for name, text in self._status_help_texts().items():
            self._register_status_help(getattr(self, name, None), text)

        # Matplotlib canvases are also interactive GUI widgets on the NMR page.
        self._register_status_help(getattr(self, 'canvas', None),
                                   'Spectrum projection display showing the current loaded data projection.')
        self._register_status_help(getattr(self, 'noiseCanvas', None),
                                   'Noise distribution display; the threshold marker follows the Threshold field.')
        for line in getattr(self, 'updateList', []):
            self._register_status_help(line, 'Project report line showing the current state of an input or result.')
        project_help = {
            'spectrum': 'Shows whether the selected spectrum file is available; use Load to read it.',
            'reference': 'Shows whether the selected 2D reference peak list is available; use Load to read it.',
            'full': 'Shows whether the selected full-dimensional peak list is available; use Load to read it.',
            'decon': 'Shows whether a deconvolution result is available for the selected spectrum; use Load to read it.',
        }
        for key, text in project_help.items():
            self._register_status_help(getattr(self, 'projectLamps', {}).get(key), text)
            self._register_status_help(getattr(self, 'projectRowLabels', {}).get(key), text)
        # Do not register section/container StaticBoxes.  Their large hit areas
        # compete with child controls for EVT_ENTER_WINDOW and can overwrite
        # specific help while the pointer is still over a child widget.

    def _register_dimension_status_help(self, dimension, nlab, sigbox, voigtbox, lorbox):
        self._register_status_help(nlab, 'Line-shape parameters for spectral dimension %d.' % dimension)
        self._register_status_help(sigbox, 'Gaussian line width in ppm for spectral dimension %d.' % dimension)
        self._register_status_help(voigtbox, 'Voigt mixing fraction from 0 to 1 for spectral dimension %d.' % dimension)
        self._register_status_help(lorbox, 'Lorentzian line width in ppm for spectral dimension %d.' % dimension)

    def _selected_topology(self):
        """Return the canonical topology selected in Workflow/ProjectState."""
        return self.state.topology()

    def apply_dataset_type(self, spectral_dimensions, pseudo_axis, *, show_error=True):
        """Apply dataset topology without depending on NMR-tab widgets.

        This is the single compatibility boundary used by Workflow and project
        loading.  It preserves the historical pseudo-4D restriction and all
        dimensional-control refresh behaviour previously owned by wx handlers.
        """
        try:
            dim = int(spectral_dimensions)
        except (TypeError, ValueError):
            return False
        if dim < 1 or dim > 4:
            return False
        pseudo = bool(pseudo_axis)
        if dim == 4 and pseudo:
            pseudo = False
            if show_error:
                errorMessage('pseudo4d not yet supported')
        self.state.sync_from_values(spectral_dimensions=dim, pseudo_axis=pseudo)
        self.dim = dim
        self.pseudo = pseudo
        self.SetDim()
        self._update_full_peak_controls()
        return True

    def _active_topology(self):
        """Return canonical topology for scientific decisions.

        ProjectState owns committed topology.  During early GUI construction,
        before state is available, the current selector values are the only
        meaningful source and are converted through the same DatasetTopology
        boundary.
        """
        state = getattr(self, 'state', None)
        if state is not None:
            return state.topology()
        return self._selected_topology()

    def SetDim(self):
        """Update the UniDec per-dimension controls from the NMR dimension selector."""
        pass
        topology = self.state.topology()
        pass
        self.dim = topology.spectral_dim_count
        if getattr(self, 'state', None) is not None:
            self.state.spectral_dimensions = topology.spectral_dim_count
            self.state.pseudo_axis = topology.has_pseudo_axis

        # Rebuild only the active dimensional columns.  DoDim reloads the values
        # from the normal decon parameter/save file, so dimension changes follow
        # the same persistence path as a freshly loaded project.
        self.DoDim(self.dim)

        # Fixed FIT radii describe the two spectral axes.  They are therefore
        # visible for ordinary 2D data and physical 3D data whose final axis
        # is pseudo/real (Protocol3P).
        show_f1_fit = (topology.spectral_dim_count == 2) or (topology.spectral_dim_count == 1 and topology.has_pseudo_axis)
        show_f2_fit = (topology.spectral_dim_count == 2)
        for control in (self.fitF1Lab, self.fitF1Box):
            control.Show(show_f1_fit)
        for control in (self.fitF2Lab, self.fitF2Box):
            control.Show(show_f2_fit)

        # Phase/distortion fitting is defined only for pseudo2D (one spectral
        # dimension plus one real/pseudo dimension).  Keep the control visible
        # beside Fit? but disable it for every other topology.
        phase_ok = topology.spectral_dim_count == 1 and topology.has_pseudo_axis
        self.cb_fitphases.Enable(phase_ok)
        self.cb_fitphases_label.Enable(phase_ok)

        # This panel is a notebook page: lay out the existing hierarchy rather
        # than SetSizerAndFit(), which can leave dynamically-created controls
        # outside the visible page on some wxPython platforms.
        self.sizer3.Layout()
        self.border2.Layout()
        self.deconSizer.Layout()
        self.fullSizer.Layout()
        self.Layout()
        self._update_full_peak_controls()


    #####################################3

    def DoDim(self, dim):
        """Build the UniDec line-shape controls as a 4-row by *dim* grid.

        Row 0 contains dimension numbers; rows 1-3 contain Gaussian width,
        Voigt fraction and Lorentzian width respectively.  The fixed row labels
        live in column 1 and active dimensions occupy columns 2 onward.
        """
        try:
            old_controls = list(self.cln)
        except Exception:
            old_controls = []

        # Explicitly detach dynamic windows before destroying them.  This avoids
        # stale GridBagSizer items when the dimensionality is reduced.
        for control in old_controls:
            try:
                self.sizer3.Detach(control)
            except Exception:
                pass
        self.CleanDim()

        dim = max(1, min(int(dim), 4))
        self.cln = []

        defaults = {
            'sig': 0.2,
            'voigt': 0.2,
            'lor': 0.2,
        }

        for dimension in range(1, dim + 1):
            col = dimension + 1
            nlab = wx.StaticText(self.peakShapeBox, label=str(dimension))
            sigbox = wx.TextCtrl(self.peakShapeBox, size=(50, 22))
            voigtbox = wx.TextCtrl(self.peakShapeBox, size=(50, 22))
            lorbox = wx.TextCtrl(self.peakShapeBox, size=(50, 22))

            sigbox.SetValue(str(ParseFlt(self.deconParFile,
                                         'sig%d' % dimension,
                                         default=defaults['sig'])))
            voigtbox.SetValue(str(ParseFlt(self.deconParFile,
                                           'voigt%d' % dimension,
                                           default=defaults['voigt'])))
            lorbox.SetValue(str(ParseFlt(self.deconParFile,
                                         'lor%d' % dimension,
                                         default=defaults['lor'])))

            # Preserve the historical public attributes used by save/load and
            # the deconvolution routines (sig1Box, voigt1Box, lorentz1Box, ...).
            setattr(self, 'n%dLab' % dimension, nlab)
            setattr(self, 'sig%dBox' % dimension, sigbox)
            setattr(self, 'voigt%dBox' % dimension, voigtbox)
            setattr(self, 'lorentz%dBox' % dimension, lorbox)

            self.sizer3.Add(nlab, (2, col), flag=wx.ALIGN_CENTER_HORIZONTAL)
            self.sizer3.Add(sigbox, (3, col), flag=wx.EXPAND)
            self.sizer3.Add(voigtbox, (4, col), flag=wx.EXPAND)
            self.sizer3.Add(lorbox, (5, col), flag=wx.EXPAND)

            if self.READ != 1:
                nlab.Disable()
                sigbox.Disable()
                voigtbox.Disable()
                lorbox.Disable()

            self.cln.extend([sigbox, voigtbox, lorbox, nlab])
            self._register_dimension_status_help(dimension, nlab, sigbox, voigtbox, lorbox)

        self._update_reorganise_controls()

        # Reflow the complete UniDec quadrant hierarchy whenever dimensionality
        # changes.  This makes the Peak shape box expand/contract immediately
        # and lets the other three boxes share the available space cleanly.
        self.sizer3.Layout()
        for sizer in (getattr(self, 'peakShapeSizer', None),
                      getattr(self, 'quadrantSizer', None),
                      getattr(self, 'border2', None),
                      getattr(self, 'deconSizer', None)):
            if sizer is not None:
                sizer.Layout()
        self.Layout()
        try:
            self.GetParent().Layout()
            self.GetParent().SendSizeEvent()
        except Exception:
            pass


    #################################################
    
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
                print(splitPath)
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

    def onGetSpecFile(self, e, textBox):
        """Choose a spectrum/peak file from SpecPath and store only its name."""
        self._sync_directory_state_only()
        default_dir = self.state.spec_dir() if getattr(self, 'state', None) else os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file in SpecPath", defaultDir=default_dir,
                            defaultFile="", wildcard="All files (*.*)|*.*", style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = os.path.abspath(dlg.GetPath())
            spec_dir = os.path.abspath(default_dir)
            try:
                common = os.path.commonpath([path, spec_dir])
            except ValueError:
                common = ''
            if common != spec_dir:
                errorMessage('Spectrum and peak-list files must be located inside SpecPath: %s' % spec_dir)
            else:
                textBox.SetValue(os.path.relpath(path, spec_dir).replace(os.sep, '/'))
        dlg.Destroy()

    def onGetDir(self, e, textBox, change_cwd=True):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.DirDialog(self, message="Choose a folder",         style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            try:
                textBox.SetValue('.' + splitPath[1])
            except:
                textBox.SetValue(path)
            print("You chose the following file(s):")
            print(path)
            if change_cwd:
                os.chdir(path)
                self.dirBox.SetValue(path)
                print("CWD: ",os.getcwd())

        dlg.Destroy()

    ##################################################
        
    def minfunc(self,x,b):

        k1=int(x[0])
        k2=int(x[1])
        slice2d=b[0]

        print('   ',k1,k2)
        return 1/(slice2d[k1,k2]),1/(slice2d[k1,k2])

    def OptSlice(self,ptC,ptC2,ptC_max,ptC2_max,slice2d):
        k1=ptC-ptC_max
        k2=ptC2-ptC2_max

        temp=slice2d[k1,k2]
        go=0
        while(go==0):
            run=0
            #print '    ',k1,k2
            for i in (-1,0,1):
                for j in (-1,0,1):
                    try:
                        if(numpy.fabs(slice2d[k1+i,k2+j])>temp):
                            temp=numpy.fabs(slice2d[k1+i,k2+j])
                            k1=k1+i
                            k2=k2+j
                            run+=1
                    except:
                        pass
            if(run==0):
                go=1
        return k1,k2

    def _load_decon_peak_list(self, infile):
        """Load a calculated nD peak list into the canonical peak-list schema.

        Reference, full, projection, and deconvolved peak lists are separate
        semantic roles, but all use the same dimension-independent record
        representation.  Legacy conn_data connectivity is deliberately not
        reconstructed here.
        """
        if not infile or not os.path.exists(infile):
            return []
        rows = []
        with open(infile) as handle:
            for line in handle:
                fields = line.split()
                if fields and not fields[0].startswith('#'):
                    rows.append(fields)
        records = self._full_peak_records(rows, dim=self._active_topology().spectral_dim_count)
        if getattr(self, 'store', None) is not None:
            self.store.save_peak_list(
                'decon',
                role='deconvolved',
                dimensionality=self._active_topology().spectral_dim_count,
                rows=rows,
                records=records,
                source=infile,
            )
        return records


    #####################################
    #guess noise levels.
    def OnButtonNoise(self,event):
        self._blit_noise_threshold()
        if(self.READ==0):
            print('Peak lists not read in. Cannot adjust noise')
            return

        thresh=str(self.dmax*float(self.threshBox.GetValue()))

        # Threshold propagation is deliberately capability based.  Not every
        # workflow creates every notebook page (in particular tabFour is not
        # present for ordinary 2D data), so an absent page must never make a
        # threshold commit fail.
        pseudo = getattr(self.parent, 'tabPseudo', None)
        if pseudo is not None and callable(getattr(pseudo, 'sync_main_threshold', None)):
            try:
                pseudo.sync_main_threshold(float(thresh), redraw=True)
            except (RuntimeError, AttributeError, ValueError):
                pass

        tab_two = getattr(self.parent, 'tabTwo', None)
        if tab_two is not None:
            try:
                tab_two.thresh = float(thresh)
                if hasattr(tab_two, 'textbox0'):
                    tab_two.textbox0.SetValue(thresh)
                draw = getattr(tab_two, 'draw_figure', None)
                if callable(draw):
                    draw()
            except (RuntimeError, AttributeError, ValueError):
                pass

        tab_four = getattr(self.parent, 'tabFour', None)
        if tab_four is not None:
            try:
                # Full/4D views expose a slightly richer threshold UI.  Only
                # touch controls that actually exist on the live page.
                tab_four.thresh = float(thresh)
                if hasattr(tab_four, 'textbox0'):
                    tab_four.textbox0.SetValue(thresh)
                if self.dim == 3 and self.PEAK != 0 and hasattr(tab_four, 'textbox_minP'):
                    tab_four.textbox_minP.SetValue(thresh)
                draw = getattr(tab_four, 'draw_figure', None)
                if self.dim == 3 and self.PEAK != 0 and callable(draw):
                    draw()
            except (RuntimeError, AttributeError, ValueError):
                pass

        # PeakFrame is modeless and project-global.  A threshold commit updates
        # the one existing editor in place; it must never create a PeakFrame.
        peak_frame = self._live_peak_frame()
        if peak_frame is not None:
            try:
                sync = getattr(peak_frame, 'sync_main_threshold', None)
                if callable(sync):
                    sync(redraw=True)
            except (RuntimeError, AttributeError, ValueError):
                pass

    ######################################

    def prepare_workflow(self, workflow_key):
        """Prepare this page for an externally requested workflow.

        This is the public workflow boundary used by the notebook.  It keeps
        launchers independent of individual wx controls while preserving the
        existing interactive handlers.  Spectrum-dependent workflows load the
        configured main spectrum when it exists; failures remain non-fatal so
        the user can repair the project from the NMR page.
        """
        if workflow_key == 'prepare':
            if hasattr(self, 'buttonProcess') and self.buttonProcess.IsEnabled():
                self.buttonProcess.SetFocus()
            return True

        if workflow_key in ('decon', 'inspect', 'slices', 'special'):
            try:
                infile = self.infileBox.GetValue()
                spectrum = self._resolve_input_path(infile) if infile else None
                if spectrum and os.path.isfile(spectrum):
                    current = getattr(self, 'spectrumfile', None)
                    if current != spectrum or getattr(self, 'data', None) is None:
                        self.OnButtonRead(None)
            except (OSError, RuntimeError, ValueError, AttributeError) as exc:
                print('Could not prepare workflow spectrum:', exc)

        if workflow_key == 'decon' and hasattr(self, 'buttonDecon'):
            if self.buttonDecon.IsEnabled():
                self.buttonDecon.SetFocus()
        return True

    def OnButtonProcess(self,event): #load process frame
        from spinDecon.gui.dialogs.processing import process as processFrame
        processFrame=importlib.reload(processFrame)
        bool=processFrame.ProcMan(self)

    def _live_peak_frame(self):
        """Return the single live PeakFrame owned by this project, if any."""
        frame = getattr(self, 'peak_frame', None)
        if frame is None:
            return None
        try:
            if frame.IsBeingDeleted():
                self.peak_frame = None
                return None
        except Exception:
            self.peak_frame = None
            return None
        return frame

    def OnButtonPeaky(self,event):  #load/focus the single peak frame
        # PeakFrame is the project-level 2D peak editor.  For physical 2D it edits
        # the authoritative Full list; other topologies retain their normal Reference semantics.
        # Opening a second instance gives two editors for the same state and, in
        # Workflow, threshold/status refreshes could accidentally create another
        # window.  Reuse the existing modeless frame instead.
        existing = self._live_peak_frame()
        if existing is not None:
            refresh = getattr(existing, 'sync_main_threshold', None)
            if callable(refresh):
                refresh(redraw=True)
            try:
                existing.Show(True)
                existing.Raise()
                existing.SetFocus()
            except Exception:
                pass
            return existing

        from spinDecon.gui.workspaces import peak_review as peakFrame
        peakFrame=importlib.reload(peakFrame)
        self.peak_frame = peakFrame.peakFrame(self)
        try:
            self.peak_frame.Bind(wx.EVT_CLOSE, self._on_peak_frame_close)
        except Exception:
            pass
        return self.peak_frame

    def _on_peak_frame_close(self, event):
        self.peak_frame = None
        event.Skip()

    def OnButtonReferencePeakList(self, event):
        from spinDecon.gui.workspaces.full_peak_list import PeakListFrame
        return PeakListFrame(self, mode='reference')

    def OnButtonFullPeakList(self, event):
        from spinDecon.gui.workspaces.full_peak_list import PeakListFrame
        return PeakListFrame(self, mode='full')

    def refresh_full_peak_list_viewers(self):
        """Refresh any already-open Full Peak List views without creating one."""
        viewers = list(getattr(self, '_full_peak_list_viewers', []) or [])
        alive = []
        for viewer in viewers:
            try:
                if viewer and not viewer.IsBeingDeleted():
                    viewer.on_refresh(None)
                    alive.append(viewer)
            except Exception:
                pass
        self._full_peak_list_viewers = alive

    def focus_full_peak_list_viewers(self, peak_name):
        """Highlight ``peak_name`` in any open Full 1D Peak List viewer."""
        for viewer in list(getattr(self, '_full_peak_list_viewers', []) or []):
            try:
                if viewer and not viewer.IsBeingDeleted():
                    viewer.focus_peak_name(peak_name)
            except Exception:
                pass


    def OnButtonReadFullPeak(self, event):
        if self.READ == 0:
            self.OnButtonRead(True)
            if self.READ == 0:
                return False
        result = self.load_full_peak_list()
        self.update_project_lamps()
        return result

    def OnButtonPeakFit(self, event=None, showFlg=True): # load fit frame
        """Open the peak-fit view; also reusable by project-summary export."""
        from spinDecon.gui.workspaces import peak_fit as peakFitFrame
        peakFitFrame=importlib.reload(peakFitFrame)
        frame = peakFitFrame.peakFitFrame(self, showFlg=showFlg)
        return frame


    #######################################
    #Read peak list file
    def OnButtonReadPeak(self,event):
        if(self.READ==0):
            print('No data. Trying to read that in first.')
            self.OnButtonRead(True,cnt=self.cnt)
            if(self.READ==0):
                print('Cannot read in data either.')
                return
            else:
                print('Successfully read in data.')
                print('Continuing...')
        self.ReadPeakListFile()
        self.update_project_lamps()
                
        """
        peak1 = 0
        peak2 = 0
        peak3 = 0
        if self.PEAK == 1:
            peak1 = self.parent.tabFour.ComboBox1.GetSelection()
            peak2 = self.parent.tabFour.ComboBox2.GetSelection()
            peak3 = self.parent.tabThree.ComboBox1.GetSelection()

        self.parent.tabFour.ComboBox1.SetSelection(peak1)
        self.parent.tabFour.ComboBox2.SetSelection(peak2)
        self.parent.tabThree.ComboBox1.SetSelection(peak3)
        self.parent.tabFour.draw_figure()
        self.parent.tabThree.draw_figure()
        """


    ################################################
    # nmrPipe functions to handle transpositions
    def _update_reorganise_controls(self):
        if self.reorganiseFrame is not None:
            try:
                self.reorganiseFrame.update_enabled_state()
            except (RuntimeError, wx.PyDeadObjectError):
                self.reorganiseFrame = None

    def OnButtonReOrganise(self, event):
        if self.reorganiseFrame is None:
            self.reorganiseFrame = ReOrganiseFrame(self)
        else:
            self.reorganiseFrame.update_enabled_state()
        self.reorganiseFrame.Show()
        self.reorganiseFrame.Raise()

    def OnButtonXY(self,event):
        infile=self.infileBox.GetValue()
        self.spectrumfile=self._resolve_input_path(infile)
        outfile=self.spectrumfile+'.tmp'
        outy=open('test.sh','w')
        outy.write('#!/bin/csh -f\n')
        outy.write('echo Transposing XYZA-YXZA ......\n')
        outy.write('cat %s                    \\\n' % self.spectrumfile)
        outy.write('| nmrPipe  -fn TP                                \\\n')
        outy.write('> %s\n' % (outfile))
        outy.close()
        os.system('csh test.sh')
        indir, scriptname = os.path.split(self.infileBox.GetValue())
        #if(self.dim==4):
        #    MakeProj4D(indir)
        self.CheckFiles(self.spectrumfile,outfile)

    def OnButtonXZ(self,event):
        infile=self.infileBox.GetValue()
        self.spectrumfile=self._resolve_input_path(infile)
        outfile=self.spectrumfile+'.tmp'
        outy=open('test.sh','w')
        outy.write('#!/bin/csh -f\n')
        outy.write('echo Transposing XYZA-ZYXA ......\n')
        outy.write('cat %s                    \\\n' % self.spectrumfile)
        outy.write('| nmrPipe  -fn ZTP                                \\\n')
        outy.write('> %s\n' % (outfile))
        outy.close()
        os.system('csh test.sh')
        indir, scriptname = os.path.split(self.infileBox.GetValue())
        #if(self.dim==4):
        #    MakeProj4D(indir)
        self.CheckFiles(self.spectrumfile,outfile)

    def OnButtonXA(self,event):
        if self._active_topology().spectral_dim_count != 4:
            print('Cannot do this - need to be 4D')
            return
        infile=self.infileBox.GetValue()
        self.spectrumfile=self._resolve_input_path(infile)
        outfile=self.spectrumfile+'.tmp'
        if(os.path.exists('tmp')==0):
            os.system('mkdir tmp')
        outy=open('test.sh','w')
        outy.write('#!/bin/csh -f\n')
        outy.write('echo Transposing XYZA-AYZX ......\n')
        outy.write('cat %s                    \\\n' % self.spectrumfile)
        outy.write('| nmrPipe  -fn ATP                                \\\n')
        outy.write('> %s\n' % (outfile))

        outy.close()
        os.system('csh test.sh')
        self.CheckFiles(self.spectrumfile,outfile)

    def OnButtonCirc(self,event):
        if self._active_topology().spectral_dim_count != 4:
            print('Cannot do this - need to be 4D')
            return
        infile=self.infileBox.GetValue()
        self.spectrumfile=self._resolve_input_path(infile)
        outfile=self.spectrumfile+'.tmp'
        if(os.path.exists('tmp')==0):
            os.system('mkdir tmp')
        outy=open('test.sh','w')
        outy.write('#!/bin/csh -f\n')
        outy.write('echo Transposing XYZA-YZAX ......\n')
        outy.write('echo Writing out as seperate files...\n')
        outy.write('nmrPipe -in %s -verb  \\\n' % self.spectrumfile)
        outy.write('|pipe2xyz -out tmp/ft/test%02d%03d.ft4 -x\n')
        outy.write('echo Transposing...\n')
        outy.write('xyz2pipe -in tmp/ft/test%02d%03d.ft4 -x -verb \\\n')
        outy.write('|pipe2xyz -out tmp/lp/test%02d%03d.ft4 -a\n')
        outy.write('echo Bringing files back together...\n')
        outy.write('xyz2pipe -in tmp/lp/test%02d%03d.ft4 -x -verb \\\n')
        outy.write('> %s \n' % outfile)
        outy.write('echo Cleaning up...\n')
        outy.write('rm -rf tmp')
        outy.close()
        os.system('csh test.sh')
        self.CheckFiles(self.spectrumfile,outfile)
        indir, scriptname = os.path.split(self.infileBox.GetValue())

    def OnButtonProject(self,event):
        print('Re-making projections...')
        indir, scriptname = os.path.split(self.infileBox.GetValue())
        spectral_dim_count = self._active_topology().spectral_dim_count
        if spectral_dim_count == 4:
            #MakeProj4D(indir)
            MakeProj4D(self._resolve_input_path(self.infileBox.GetValue()))
        if spectral_dim_count == 3:
            #MakeProj3D(indir)
            MakeProj3D(self._resolve_input_path(self.infileBox.GetValue()),OneD=True)
        if spectral_dim_count == 2:
            MakeProj2D(self._resolve_input_path(self.infileBox.GetValue()))


    def CheckFiles(self,infile,outfile):
        a=os.stat(infile).st_size
        b=os.stat(outfile).st_size
        if(a==b):
            print('Sizes are the same! Saving transposition')
            os.system('mv %s %s' % (outfile,infile))
        else:
            print('An error has occcured - keeping both files')
            print('The two are not the same size, indicating')
            print('some kind of error.')
            print('infile :',infile)
            print('outfile:',outfile)
            os.system('mv %s %s' % (outfile,infile))
        self.OnButtonRead(True)


    #######################################3
    #Remove notebook page
    def DeletePage(self,pageTitle):
        #pageTitle='2Dplanes'
        for index in range(self.parent.GetPageCount()):
            if self.parent.GetPageText(index) == pageTitle:
                self.parent.DeletePage(index)
                self.parent.SendSizeEvent()
                break

    ##################################################
    #Functions for handling STD data.
    def read_STD(self, infile, spec):
        print('Reading:',infile,'in mode',spec)
        dat=[]
        inny=open(infile)
        for line in inny.readlines():#read in data file
            test=line.split()
            if(len(test)>0):
                dat.append((float(test[1]),float(test[3]),float(test[4]),float(test[5])))
        dat=numpy.array(dat)

        times=numpy.unique(dat[:,0])
        print(dat[:,0])
        ppms=numpy.unique(dat[:,1])
        dat=dat.reshape(len(times),len(ppms),4)
        #self.times=times  #save total range of times

        times = times
        print('Making projections...')

        xvals=numpy.average(dat[:,:,1],axis=0)

        I=numpy.sum(dat[:,:,2:],axis=0)
        I=I.transpose()

            # print 'DIFFERENCE'
        STD=I[1,:]-I[0,:] #select spectrum 1

        data=I[1,:] #select spectrum 1
        return xvals, data, STD, times

    def read_STD_pipe(self, spectrumfile):
        dic, data = ng.pipe.read(spectrumfile)
        uc = ng.pipe.make_uc(dic, data, dim=0)
        index0 = uc.ppm_scale()
        print(data.shape)
        return index0, data

    # def read_STD_spectra(self, raw, std):
    #     print('Reading raw spectrum in...')
    #     index0, data = self.read_STD_pipe(raw)
    #     index0_STD, STD = self.read_STD_pipe(std)
    #     if (index0 == index0_STD).all():
    #         # self.index0 = index0
    #         return(data, STD, index0)
    #     else:
    #         dlg = wx.MessageDialog(self, message="Indexes are not the same for STD and raw", style = wx.ICON_ERROR| wx.OK)
    #         if dlg.ShowModal() == wx.ID_OK:
    #             return

    def read_STD_spectra(self, infile):
        
        print('Reading raw spectrum in...')
        pars=(GetParBrukFile(os.path.dirname(os.path.abspath(infile))+'/fq2list'))
        pars2=(GetParBrukFile(os.path.dirname(os.path.abspath(infile))+'/vdlist'))
        excite=[] #excitation ppm values
        if pars[0][0] == 'P':
            pars = pars[1:]
            for par in pars:
                excite.append(float(par[0]))
            
        else:
            print('Excitation list is in the wrong format - currently only PPM is supported')

        dic,data = ng.pipe.read(infile) #read fids
        Size=data.shape

        uc0 = ng.pipe.make_uc(dic,data,dim=0)
        uc1 = ng.pipe.make_uc(dic,data,dim=1)
        index=[] #get ppms.
        for i in range((Size[1])):
            index.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))

        data=data.reshape((len(pars2),len(pars),Size[1])) #reshape the data

        frq=numpy.zeros_like(data)
        mix=numpy.zeros_like(data)

        mixTotal=[]
        for i,p in enumerate(pars2): #go through vdlist, get mixing times
            val=float(re.findall(r"\d+\.?\d*",p[0])[0]) #get numbers from file
            if p[0][-1:] == 'm':  #if there's an 'm', divide by 1000
                val = val/1000.
            mix[i,:,:]=val  #setup 3D array with mixing times
            mixTotal.append(val) #append time
        for i,p in enumerate(pars):
            val=excite[i]
            frq[:,i,:]=val  #setup 3D array with frequencies
            print(val)

        mixTotal=numpy.array(mixTotal) #turn to numpy
        #unique times:
        mix=numpy.unique(mix)  #get unique mixing times.
        print(mix)       #unique mixing times
        print(mixTotal)  #tota mixing times
        datNew=numpy.zeros((len(mix),len(pars),Size[1])) #get new array for data
        #print numpy.sum(data[0,0,:])  #####PROBLEM HERE#####
            
        excite=numpy.array(excite) #turn excitation frequencies into numpy array
        argyExMax=numpy.argmax(numpy.fabs(excite))                           #find maximum number in exciation (off resonance)

        print(data.shape)
        print(datNew.shape)
        print(mixTotal, mix)
        for m,mi in enumerate(mix):  #for each unique mixing time
            for e,ex in enumerate(excite):  #for each unique excitation time
                mask=(numpy.absolute(mixTotal-mi)<0.000001)  #get indicies for mixing times that are aligned with unique mixing times
                #print mask
                #print numpy.sum(data[mask,:,:])
                #print 'a',datNew[m,e,:].shape
                #print 'b',data[mask,e,:].shape
                #print 'c',numpy.sum(data[mask,e,:])
                datNew[m,e,:]=numpy.sum(data[mask,e,:],axis=0)  #sum mixing times to get unique mixing time array only.

        return data[:,0,:], data[:,1,:]-data[:,0,:], index, mix


    def _spec_output_dir(self):
        self._sync_directory_state_only()
        return self.state.spec_dir() if getattr(self, 'state', None) else './spec'

    def _raw_output_dir(self):
        self._sync_directory_state_only()
        return self.state.raw_dir() if getattr(self, 'state', None) else './raw'

    def get_fuda_dir(self):
        """Return the canonical FUDA workspace below SpecPath (fit/)."""
        return os.path.join(self._spec_output_dir(), 'fit')

    def get_fuda_peak_file(self):
        return os.path.join(self.get_fuda_dir(), 'peak.fuda')

    def get_fuda_parameter_file(self):
        return os.path.join(self.get_fuda_dir(), 'param.fuda')

    def missing_pseudo3d_fit_peaks(self):
        """Return reference-peak names that do not yet have fit results.

        Protocol3P writes per-peak results into ``spec/fit``.  A peak is considered available only when both its FUDA-style ``.dat``
        and ``.out`` results exist.  The reference list is authoritative.
        """
        refs = list(self.get_reference_peaks() or [])
        fit_dir = self.get_fuda_dir()
        missing = []
        for peak in refs:
            name = str(getattr(peak, 'name', '') or '').strip()
            if not name:
                continue
            if not all(os.path.isfile(os.path.join(fit_dir, name + ext))
                       for ext in ('.dat', '.out')):
                missing.append(name)
        return missing

    def _workflow_debug(self, message):
        """Legacy pseudo3D debug hook (intentionally silent)."""
        return

    def _pseudo2d_debug(self, message):
        return

    def ensure_workflow_reference_stage_loaded(self):
        """Materialise the completed Establish reference peaks stage.

        Workflow may mark the stage complete from the configured list file on
        disk.  Before any later pseudo-axis stage is evaluated, however, that
        same list must be loaded into the shared DataStore.  Cold-start loading
        requires the main spectrum first because the legacy peak-list reader
        derives indices/slices from it.  Keeping that ordering here avoids
        making later stages repair an incompletely initialised earlier stage.
        """
        infile = self.infileBox.GetValue().strip() if hasattr(self, 'infileBox') else ''
        spectrum = self._resolve_input_path(infile) if infile else ''
        if not spectrum or not os.path.isfile(spectrum):
            self._workflow_debug('reference stage load deferred: processed spectrum is absent: %r' % spectrum)
            return False

        current = getattr(self, 'spectrumfile', None)
        if current != spectrum or getattr(self, 'data', None) is None or not getattr(self, 'READ', 0):
            self._workflow_debug('reference stage: loading spectrum before reference list: %r' % spectrum)
            # Cold-start materialisation of the spectrum for the project whose
            # workflow flags were just restored. OnButtonRead normally resets
            # DataStore (correct for an explicit spectrum change), but that
            # reset must not erase persisted review acceptance here.
            self._preserve_persisted_workflow_review_on_next_read = True
            try:
                self.OnButtonRead(None)
            except Exception as exc:
                self._workflow_debug('reference stage spectrum load failed: %r' % (exc,))
                return False
            finally:
                self._preserve_persisted_workflow_review_on_next_read = False
        if not (getattr(self, 'READ', 0) and getattr(self, 'data', None) is not None):
            self._workflow_debug('reference stage load deferred: spectrum did not initialise')
            return False

        return self.ensure_reference_peak_list_loaded()

    def ensure_workflow_spectrum_loaded(self):
        """Materialise the configured spectrum for a workflow action."""
        infile = self.infileBox.GetValue().strip() if hasattr(self, 'infileBox') else ''
        spectrum = self._resolve_input_path(infile) if infile else ''
        if not spectrum or not os.path.isfile(spectrum):
            return False
        current = getattr(self, 'spectrumfile', None)
        if current != spectrum or getattr(self, 'data', None) is None or not getattr(self, 'READ', 0):
            self._preserve_persisted_workflow_review_on_next_read = True
            try:
                self.OnButtonRead(None)
            finally:
                self._preserve_persisted_workflow_review_on_next_read = False
        return bool(getattr(self, 'READ', 0) and getattr(self, 'data', None) is not None)

    def ensure_full_peak_list_loaded(self):
        """Materialise the configured Full nD list without changing reference peaks."""
        value = self.fullPeakBox.GetValue().strip() if hasattr(self, 'fullPeakBox') else ''
        path = self._resolve_spec_file(value) if value else ''
        if not path or not os.path.isfile(path):
            return False
        payload = self.store.peak_lists.get('full', {}) if getattr(self, 'store', None) is not None else {}
        loaded_path = str(payload.get('source_path') or '')
        rows = payload.get('rows') or []
        try:
            same_file = bool(loaded_path) and os.path.abspath(loaded_path) == os.path.abspath(path)
        except Exception:
            same_file = False
        if not rows or not same_file:
            return bool(self.load_full_peak_list(path, quiet=True))
        return True

    def ensure_deconvolution_loaded(self):
        """Materialise the calculated spectrum after its peak-list dependencies.

        The legacy 3D decon loader builds calculated 1D slices using ``peak``
        and ``pkIdx``.  Consequently reference peaks must already be loaded;
        loading the .decon file first can otherwise fail with an index error.
        """
        if not self.ensure_workflow_spectrum_loaded():
            return False
        topology = self._active_topology()
        if topology.spectral_dim_count >= 3 and not topology.has_pseudo_axis:
            if not self.ensure_reference_peak_list_loaded():
                return False
        spectrum = getattr(self, 'spectrumfile', '') or self._resolve_input_path(self.infileBox.GetValue())
        decon_path = self._active_deconvolution_path(spectrum)
        if not decon_path or not os.path.isfile(decon_path):
            return False
        if getattr(self, 'DECON', 0) and getattr(self, 'datadec', None) is not None:
            return True
        if topology.spectral_dim_count == 1 and topology.has_pseudo_axis:
            projection = self.get_pseudo2d_projection_data(ensure_file=False)
            projection_path = projection.get('path') if projection else ''
            return bool(projection_path and self._load_pseudo2d_projection_decon_outputs(projection_path))
        if topology.spectral_dim_count == 2 and topology.has_pseudo_axis:
            # Pseudo3D has no general full-spectrum .decon product.
            return True
        return bool(self._load_decon_outputs(spectrum, load_peaks=False))

    def ensure_workflow_review_inputs_loaded(self):
        """Load Review picked peaks inputs in strict dependency order.

        Review requires the main spectrum, the reference list used to index
        slices, the authoritative Full nD list, and (where the topology has a
        full calculated spectrum) the deconvolution result.
        """
        if not self.ensure_workflow_spectrum_loaded():
            return False, 'The processed spectrum file could not be loaded.'
        topology = self._active_topology()
        if topology.spectral_dim_count >= 3 and not topology.has_pseudo_axis:
            if not self.ensure_reference_peak_list_loaded():
                return False, 'The reference peak list could not be loaded.'
        if not self.ensure_full_peak_list_loaded():
            return False, 'The full peak list could not be loaded.'
        if not (topology.spectral_dim_count == 2 and topology.has_pseudo_axis):
            if not self.ensure_deconvolution_loaded():
                return False, 'The deconvolution result could not be loaded.'
        return True, ''

    def ensure_reference_peak_list_loaded(self):
        """Load the configured Reference 2D list into the shared DataStore."""
        value = self.referencePeakBox.GetValue().strip()
        path = self._resolve_spec_file(value) if value else ''
        self._workflow_debug('ensure_reference_peak_list_loaded: box=%r resolved=%r exists=%s dim=%r pseudo=%r READ=%r' %
                             (value, path, bool(path and os.path.isfile(path)), getattr(self, 'dim', None),
                              self.state.pseudo_axis if getattr(self, 'state', None) is not None else None, getattr(self, 'READ', None)))
        if not path or not os.path.isfile(path):
            self._workflow_debug('ABORT reference load: configured file is absent.')
            return False
        payload = self.store.peak_lists.get('reference', {}) if getattr(self, 'store', None) is not None else {}
        loaded_path = str(payload.get('source_path') or '')
        peaks = payload.get('peaks') or []
        self._workflow_debug('reference store before load: count=%d source_path=%r' % (len(peaks), loaded_path))
        try:
            same_file = loaded_path and os.path.abspath(loaded_path) == os.path.abspath(path)
        except Exception as exc:
            same_file = False
            self._workflow_debug('path comparison failed: %r' % exc)
        if not peaks or not same_file:
            self._workflow_debug('calling OnButtonReadPeak(None) now')
            try:
                self.OnButtonReadPeak(None)
            except Exception as exc:
                import traceback
                self._workflow_debug('EXCEPTION from OnButtonReadPeak: %r' % exc)
                traceback.print_exc()
                return False
        refs = list(self.get_reference_peaks() or [])
        payload = self.store.peak_lists.get('reference', {}) if getattr(self, 'store', None) is not None else {}
        self._workflow_debug('reference store after load: count=%d source_path=%r PEAK=%r' %
                             (len(refs), payload.get('source_path'), getattr(self, 'PEAK', None)))
        return bool(refs)

    def ensure_pseudo3d_fit_results(self, force_recompute=False):
        """Ensure Protocol3P fit files exist for every reference peak.

        Missing results trigger the established restricted reconstruction with
        Fit and Use 2D peak list selected.  Returns True when all files already
        exist, False when a run was launched (or could not be launched).
        """
        self._workflow_debug('ensure_pseudo3d_fit_results ENTER: pseudo3d_topology=%s fit_dir=%r' % (self._is_pseudo3d_topology(), self.get_fuda_dir()))
        if not self._is_pseudo3d_topology():
            self._workflow_debug('ABORT fit ensure: _is_pseudo3d_topology() is False')
            return True
        # Protocol3P is keyed by the Reference 2D list.  Load it here as a
        # defensive invariant as well as from the workflow transition.
        if not self.ensure_reference_peak_list_loaded():
            return False
        refs = list(self.get_reference_peaks() or [])
        if not refs:
            return False
        missing = self.missing_pseudo3d_fit_peaks()
        force_recompute = bool(force_recompute or getattr(self, 'pseudo_intensities_stale', False))
        self._workflow_debug('fit check: reference_count=%d missing=%r force_recompute=%s' % (len(refs), missing, force_recompute))
        if not missing and not force_recompute:
            # Existing files are just as authoritative as files produced in
            # this session. Publish the evidence so Workflow can advance
            # without requiring the Pseudo3D page to have been opened first.
            if getattr(self, 'store', None) is not None:
                self.store.mark_pseudo_intensities_ready(
                    source='existing_protocol3p_fit',
                    invalidate_review=False,
                    fit_directory=self.get_fuda_dir(),
                    reference_peak_count=len(refs))
                self._notify_analysis_changed()
            return True
        # These are the main-tab controls labelled "Fit?" and
        # "Use 2D peaklist". OnButtonRecon then selects Protocol3P for a
        # physical 2D+pseudo-axis dataset.
        self.cb_decback.SetValue(True)
        self.cb_decon3d.SetValue(True)
        self._workflow_debug('Recon controls set: Fit?=%s Use 2D peaklist=%s; invoking OnButtonRecon' % (self.cb_decback.IsChecked(), self.cb_decon3d.IsChecked()))
        if getattr(self, 'calcy', None) is not None:
            self.calcy.append_text(
                '\nPseudo3D fit results are missing for: %s. '
                'Running Recon with Fit and Use 2D peak list.\n' % (', '.join(missing) if missing else 'peak-shape revision requires all peaks to be refitted'))
        result = self.OnButtonRecon(None)
        self._workflow_debug('OnButtonRecon returned %r (async launch may normally return None)' % (result,))
        return False

    def get_spectrum_path(self):
        """Return the resolved spectrum path for external analysis programs."""
        self._sync_directory_state_only()
        if getattr(self, 'state', None) is not None:
            return self.state.spectrum_path()
        return self._resolve_spec_file(self.get_project_input_file())

    def _resolve_input_path(self, infile):
        # nmrPipe files are always resolved relative to SpecPath.
        return self._resolve_spec_file(infile)


    def get_pseudo2d_projection_data(self, ensure_file=False):
        """Return the canonical 1D spectral projection for 1-spectral+pseudo data.

        The pseudo/real axis is reduced by summation.  When ``ensure_file`` is
        true the same trace is materialised as an NMRPipe 1D file below
        ``SpecPath/projections1D`` so external UniDec programs operate on the
        identical data shown in the Projections window.
        """
        topology = self._active_topology()
        if not (topology.spectral_dim_count == 1 and topology.has_pseudo_axis and
                topology.physical_dim_count == 2):
            return None
        axes = self._spectral_physical_axes()
        if len(axes) != 1:
            return None
        spectral_axis, label = axes[0]
        data = numpy.asarray(self.data)
        pseudo_axis = topology.pseudo_axis.physical_index
        projection = numpy.sum(data, axis=pseudo_axis)
        projection = numpy.asarray(projection).squeeze()
        index = numpy.asarray(getattr(self, 'index%d' % spectral_axis))
        if projection.ndim != 1 or projection.size != index.size:
            raise ValueError('Pseudo2D projection does not match the spectral axis')

        path = os.path.join(self._spec_output_dir(), 'projections1D', str(label) + '.dat')
        if ensure_file:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            # Reuse the established NMRPipe projection writer so the generated
            # 1D header retains the correct spectral calibration.
            # Pseudo2D projects deliberately retain the scientific topology as
            # one spectral dimension even though the stored NMRPipe array has two
            # physical axes.  The legacy projection writer validates FDDIMCOUNT
            # as a *storage* dimension count, so give it a private header copy
            # advertising the two physical axes.  Do not mutate self.dic: the
            # project topology remains 1 spectral + 1 pseudo.
            projection_dic = self.dic.copy()
            projection_dic['FDDIMCOUNT'] = 2
            nmrglue_project2D_1D(projection_dic, data, folder=os.path.dirname(path), projection_type='sum')
            if not os.path.exists(path):
                raise FileNotFoundError('Pseudo2D spectral projection was not created: %s' % path)
            dic1d, file_data = ng.pipe.read(path)
            projection = numpy.asarray(file_data).squeeze()
        else:
            dic1d = None
        return {'path': path, 'dic': dic1d, 'data': projection, 'index': index, 'label': str(label)}

    def _load_pseudo2d_projection_peaks(self, peak_path):
        """Read/publish the compact SpinUniDec pseudo2D projection peak list.

        Adapted SpinUniDec writes exactly: Number, PPM, Intensity.  For testing
        with an older binary, the historical five-column (_N, 0, 0, ppm, I)
        form is also accepted and normalised in memory.
        """
        peaks = []
        if not peak_path or not os.path.exists(peak_path):
            return peaks
        with open(peak_path, 'r') as handle:
            for line in handle:
                fields = line.split()
                if not fields or line.lstrip().startswith('#'):
                    continue
                try:
                    if len(fields) >= 5:
                        number = int(str(fields[0]).lstrip('_') or len(peaks) + 1)
                        ppm, intensity = float(fields[3]), float(fields[4])
                    elif len(fields) >= 3:
                        number = int(float(fields[0]))
                        ppm, intensity = float(fields[1]), float(fields[2])
                    else:
                        continue
                except (TypeError, ValueError):
                    continue
                peaks.append({'number': number, 'ppm': ppm, 'intensity': intensity})
        if getattr(self, 'store', None) is not None:
            self.store.save_peak_list('pseudo2d_projection', peaks=peaks, rows=peaks,
                                      dimension=1, source_path=peak_path)
        self.pseudo2d_projection_peaks = peaks
        return peaks

    def _load_pseudo2d_projection_decon_outputs(self, projection_path):
        """Publish 1D decon output produced from a pseudo2D projection."""
        decon_path = str(projection_path) + '.decon'
        if not os.path.exists(decon_path):
            return False
        dic, data = ng.pipe.read(decon_path)
        data = numpy.asarray(data).squeeze()
        raw = self.get_pseudo2d_projection_data(ensure_file=False)
        if raw is None or data.ndim != 1 or data.size != numpy.asarray(raw['data']).size:
            print('Pseudo2D projection deconvolution shape does not match the raw projection:', decon_path)
            return False
        if getattr(self, 'store', None) is not None:
            self.store.save_spectrum('pseudo2d_projection_decon', dic=dic, data=data,
                                     spectrumfile=decon_path, dim=1, labb=[raw['label']])
        self.pseudo2d_projection_decon = {'dic': dic, 'data': data, 'path': decon_path}
        self.DECON = 1
        return True

    def _projection_search_dirs(self):
        """Projection artefacts live exclusively below SpecPath."""
        base_dir = self._spec_output_dir()
        return [
            os.path.join(base_dir, 'projections'),
            os.path.join(base_dir, 'projection_decon'),
            os.path.join(base_dir, 'projections1D'),
        ]

    def _resolve_projection_input_path(self, left, right, decon=True, transpose='n'):
        """Resolve a raw or calculated projection without crossing namespaces."""
        left = str(left or '').strip()
        right = str(right or '').strip()
        if not left or not right:
            return None

        if decon:
            cache_keys = [
                ('decon_projection', left, right, transpose),
                ('decon_projection', right, left, transpose),
            ]
        else:
            cache_keys = [
                (left, right, transpose),
                (right, left, transpose),
                ('projections', left, right, transpose),
                ('projections', right, left, transpose),
            ]

        raw_spectrum = ''
        try:
            raw_spectrum = os.path.abspath(self.spectrumfile) if getattr(self, 'spectrumfile', '') else ''
        except Exception:
            raw_spectrum = ''

        if getattr(self, 'store', None) is not None:
            for key in cache_keys:
                payload = self.store.projections.get(key)
                if not payload or not payload.get('source'):
                    continue
                source = str(payload['source'])
                if not source:
                    continue
                if decon:
                    if source.endswith('.decon.dat') or os.path.exists(source):
                        return source
                    continue
                try:
                    abs_source = os.path.abspath(source)
                except Exception:
                    abs_source = source
                if abs_source != raw_spectrum and source.lower().endswith('.dat') and not source.lower().endswith('.decon.dat'):
                    return source

        wanted = (
            (f'{left}.{right}.decon.dat', f'{right}.{left}.decon.dat')
            if decon else
            (f'{left}.{right}.dat', f'{right}.{left}.dat')
        )
        for proj_dir in self._projection_search_dirs():
            if not os.path.isdir(proj_dir):
                continue
            for name in wanted:
                path = os.path.join(proj_dir, name)
                if os.path.exists(path):
                    return path
            try:
                for filename in sorted(os.listdir(proj_dir)):
                    if not filename.endswith('.dat'):
                        continue
                    is_decon = filename.endswith('.decon.dat')
                    if is_decon != bool(decon):
                        continue
                    stem = filename[:-4]
                    if is_decon:
                        stem = stem[:-6]
                    parts = stem.split('.')
                    if left in parts and right in parts:
                        return os.path.join(proj_dir, filename)
            except Exception:
                continue
        return None

    ############################################
    # MAIN READ FUNCTION
    def OnButtonRead(self,event):
        self.update_project_lamps()
        indir=self.dirBox.GetValue()
        infile=self.infileBox.GetValue()
        resolved = self._resolve_input_path(infile)
        
        # Dataset topology is selected in Workflow and already lives in state.
        self.dim = int(getattr(self.state, 'spectral_dimensions', 0) or 1)
        self.pseudo = bool(self.state.pseudo_axis)
        pass

        self.spectrumfile=resolved

        if(os.path.exists(self.spectrumfile)==0):
            print('Cannot find input file:')
            print(self.spectrumfile)
            return
        preserve_review = bool(getattr(self, '_preserve_persisted_workflow_review_on_next_read', False))
        preserved_review = None
        if preserve_review and getattr(self, 'store', None) is not None:
            analysis = getattr(self.store, 'analysis', {})
            if analysis.get('pseudo_series_reviewed'):
                preserved_review = {
                    'pseudo_series_reviewed': True,
                    'pseudo_series_review': dict(analysis.get('pseudo_series_review') or {}),
                }
        pass
        self._reset_store()
        if preserved_review and getattr(self, 'store', None) is not None:
            self.store.analysis.update(preserved_review)
            pass
        pass
        if "STD" in infile:
            self.uSTA = True
            if '.ft2' in resolved:
                
                self.data, self.STD, self.index0, self.mixingTimes = self.read_STD_spectra(resolved)
                self.labb = ['1H']
                self.dmax = numpy.max(self.data)
                self.parent.AddTabSTD(True, self)
                # self.parent.AdduSTA_sims_Tab(True, self)

                return
                
        pass
        self.makeinp('', resolved)  #read in and store main data array.
        try:
            _dbg_top = self.state.topology()
            _dbg_top_s = 'spectral=%r pseudo=%r physical=%r pseudo_index=%r' % (
                _dbg_top.spectral_dim_count, _dbg_top.has_pseudo_axis, _dbg_top.physical_dim_count, _dbg_top.pseudo_physical_index)
        except Exception as _dbg_exc:
            _dbg_top_s = 'ERROR %r' % (_dbg_exc,)
        pass
        # ``dmax`` is a derived cache, not project identity.  Older project files
        # and pseudo-dimensional load paths can legitimately reach this point
        # without having restored it.  Establish it from the just-loaded spectrum
        # before any dependent viewer (notably Projections) is constructed.
        if getattr(self, 'dmax', None) is None and getattr(self, 'data', None) is not None:
            self.dmax = float(numpy.max(numpy.fabs(numpy.asarray(self.data))))
        # A project may reach this read boundary before the dimension ComboBox
        # has a selection (notably legacy/pseudo-2D project loads).  At this
        # point the NMRPipe array is authoritative: physical ndim is known and
        # a pseudo axis contributes one non-spectral dimension.  Recover the
        # canonical spectral count rather than allowing dimension=0 into
        # ProjectState.topology().
        if int(self.dim or 0) < 1:
            physical_ndim = int(getattr(self.data, 'ndim', 0) or 0)
            has_pseudo = bool(self.state.pseudo_axis)
            inferred = physical_ndim - int(has_pseudo)
            if inferred < 1:
                inferred = physical_ndim
            if inferred >= 1:
                self.dim = inferred
                self.state.spectral_dimensions = inferred
        try:
            self.analyse_noise_spectrum(self.data)
        except Exception as exc:
            # Noise analysis is diagnostic and must never prevent a spectrum
            # from loading successfully.
            print('Noise analysis failed:', exc)
        # Full nD is a SpecPath-relative project file.  Only derive the
        # conventional <spectrum>.<n>D.list name when the project/user has not
        # explicitly selected a Full nD list.  Decon completion deliberately
        # replaces this value with the newly generated list below.
        if not self.fullPeakBox.GetValue().strip():
            full_path = self._peak_list_path_for_spectrum(self.spectrumfile, self.dim)
            self.fullPeakBox.SetValue(self.state._spec_relative(full_path))
            if getattr(self, 'state', None) is not None:
                self.state.full_peak_file = self.fullPeakBox.GetValue()
        self._update_full_peak_controls()
        if self.state is not None:
            try:
                self.state.sync_from_values(
                    working_dir=self.dirBox.GetValue(),
                    input_file=self.infileBox.GetValue(),
                    peak_file=self.peakBox.GetValue(),
                    spec_path=self.specPathBox.GetValue(),
                    dimension=int(self.dim),
                    pseudo_axis=bool(self.state.pseudo_axis),
                    sym_mode=self.cb_grid.IsChecked(),
                    decon_bore=self.cb_decon3d.IsChecked(),
                )
            except Exception:
                pass
        # Projection products are required by genuine >=3D spectral data and
        # by 2-spectral + pseudo data (three physical axes).  The old
        # ``self.dim >= 3`` test used physical-dimensional semantics and became
        # false for pseudo-3D once ``dim`` was standardised as spectral only.
        topology = self.state.topology() if self.state is not None else None
        needs_projection_cache = bool(
            topology is not None and
            (topology.spectral_dim_count >= 3 or
             (topology.spectral_dim_count == 2 and topology.has_pseudo_axis and
              topology.physical_dim_count == 3))
        )
        if needs_projection_cache:
            base_dir = self._spec_output_dir()
            projection_dir = os.path.join(base_dir, "projections")
            self._cache_projection_folder(projection_dir)
            self._cache_projection_folder(projection_dir, key_prefix="decon_projection")
            self._cache_projection_folder(os.path.join(base_dir, "projections1D"))

        print('Project load complete.')
        self.update_project_lamps()

        self.parent.KillPage('1Ddeconv')
        self.parent.KillPage('2Dslices')


        # Set up tabs from the canonical topology established by makeinp().
        # self.dim and the GUI controls are compatibility/display aliases only;
        # they must not decide which scientific viewer is constructed.
        topology = self._active_topology()
        spectral_dim_count = topology.spectral_dim_count
        pass

        if spectral_dim_count == 1 and topology.has_pseudo_axis:
            self._draw_main_pseudo2d_projection()

        if spectral_dim_count == 1:
            if topology.has_pseudo_axis:
                # Pseudo2D: derived 1D spectral projection plus raw pseudo-series inspector.
                self.parent.AddTabTwo(True, self)
                self.parent.AddTabPseudo2D(True, self)
                try:
                    self.parent.select_page('Projections')
                except Exception:
                    pass
            else:
                self.parent.AddTab1D(True, self)

        elif spectral_dim_count == 2:
            # 2 spectral + pseudo has three physical axes; topology, not self.dim, distinguishes it.
            self.parent.AddTabTwo(True, self)
            # Both true 2D and 2D+pseudo use the shared fitting workspace; the
            # latter additionally carries its real pseudo axis in topology.
            self.parent.AddTabPseudo3D(True, self)

        elif spectral_dim_count == 3:
            self.parent.AddTabTwo(True, self)
            self.parent.AddTabFive(True, self)

        elif spectral_dim_count == 4:
            self.parent.AddTabTwo(True, self)

        self.READ=1
        self.pre_read_enabling()
        
        self.Status()


    ##############################################
        
    def OnButtonExtract(self,event):

        os.system('rm -rf out')
        self.OnButtonRead(True)

    def IntToBool(self,booly):
        if(booly==True):
            return 'y'
        else:
            return 'n'

    def OnButtonQuit(self,event):
        print('exiting')
        sys.exit(100)


    #########################################
    # Set defaults from input file
    def OnButtonLoad(self,event):

        loaded_dim = int(str(Parse(self.deconParFile, 'dim')))
        try:
            loaded_pseudo = bool(int(ParseFlt(self.deconParFile, 'pseudo')))
        except Exception:
            loaded_pseudo = bool(Parse(self.deconParFile, 'pseudo'))
        self.state.sync_from_values(dimension=loaded_dim, pseudo_axis=loaded_pseudo)
        self.dim = loaded_dim
        self.pseudo = loaded_pseudo
            
        #self.dirBox.SetValue(str(indir))
        #if(str(indir)!='0'):
        #    os.chdir(str(indir))
        # Missing UniDec settings mean "use the application default".  They
        # are deliberately not written merely by loading a project.
        ncpus_value = ParseInt(self.deconParFile, 'ncpus', default=available_cpu_count())
        self.coreBox.SetValue(str(int(ncpus_value)))
        self.threshBox.SetValue(str(ParseFlt(self.deconParFile, 'thresh', default=UNIDEC_DEFAULTS['thresh'])))
        self.facBox.SetValue(str(ParseFlt(self.deconParFile, 'fac', default=UNIDEC_DEFAULTS['fac'])))
        self.convBox.SetValue(str(ParseFlt(self.deconParFile, 'conv', default=UNIDEC_DEFAULTS['conv'])))
        self.maxiterBox.SetValue(str(int(ParseFlt(self.deconParFile, 'maxiter', default=UNIDEC_DEFAULTS['maxiter']))))
        fitrad_value = str(Parse(self.deconParFile, 'FitRad', default='') or '').strip()
        self.fitRadBox.SetValue('' if fitrad_value == '0' else fitrad_value)
        fitf1_value = str(Parse(self.deconParFile, '3p_radF1', default='') or '').strip()
        fitf2_value = str(Parse(self.deconParFile, '3p_radF2', default='') or '').strip()
        self.fitF1Box.SetValue('' if fitf1_value == '0' else fitf1_value)
        self.fitF2Box.SetValue('' if fitf2_value == '0' else fitf2_value)
        self.outPathBox.SetValue(str(Parse(self.deconParFile, 'fiddir', default='./raw')))
        self.specPathBox.SetValue(str(Parse(self.deconParFile, 'specPath', default='./spec')))
        # Establish the directory roots before canonicalising any spectrum-
        # associated value.  All three file controls are displayed relative
        # to SpecPath; legacy values that already contain './spec/' are
        # normalised here rather than leaking that prefix into the GUI.
        self._sync_directory_state_only()
        infile_value = str(Parse(self.deconParFile, 'infile', default='') or '').strip()
        peak_value = str(Parse(self.deconParFile, 'peakfile', default='') or '').strip()
        full_peak_value = str(Parse(self.deconParFile, 'fullPeakFile', default='') or '').strip()
        if infile_value == '0': infile_value = ''
        if peak_value == '0': peak_value = ''
        if full_peak_value == '0': full_peak_value = ''
        self.infileBox.SetValue(self.state._spec_relative(infile_value))
        self.peakBox.SetValue(self.state._spec_relative(peak_value))
        # Older project files do not contain fullPeakFile; leave it empty so
        # OnButtonRead derives the conventional spectrum-associated name.
        self.fullPeakBox.SetValue(self.state._spec_relative(full_peak_value))
        
        try:
            self.pseudo = bool(self.state.pseudo_axis)
        except:
            pass

        if(Parse(self.deconParFile,'symmode')=='y'):
            self.cb_grid.SetValue(1)
        else:
            self.cb_grid.SetValue(0)

        if(Parse(self.deconParFile,'deconBore')=='y'):
            self.cb_decon3d.SetValue(1)
        else:
            self.cb_decon3d.SetValue(0)

        self.cb_enhance.SetValue(Parse(self.deconParFile, 'enhance', default='n') == 'y')
        self.cb_fitphases.SetValue(Parse(self.deconParFile, 'fitPhases', default='n') == 'y')

        self.SetDim() #bring in label, noise, peakwidth
        self.Status()
        self.state.set_parameter_file(self.deconParFile)
        self.state.sync_from_values(
            working_dir=self.dirBox.GetValue(),
            raw_path=self.outPathBox.GetValue(),
            input_file=self.infileBox.GetValue(),
            peak_file=self.peakBox.GetValue(),
            full_peak_file=self.fullPeakBox.GetValue(),
            spec_path=self.specPathBox.GetValue(),
            dimension=int(str(Parse(self.deconParFile,'dim'))),
            pseudo_axis=bool(self.state.pseudo_axis),
            sym_mode=self.cb_grid.IsChecked(),
            decon_bore=self.cb_decon3d.IsChecked(),
        )
        self.state.loaded = True
        self._restore_workflow_flags_from_parameter_file()
        self._notify_analysis_changed()

    # Save boxes to file.
    def OnButtonSave(self,event):
        write={}
        write['indir']=self.dirBox.GetValue()
        write['infile']=self.infileBox.GetValue()
        write['peakfile']=self.peakBox.GetValue()
        # Store Full nD exactly like infile/reference: relative to SpecPath,
        # preserving any projection/subfolder component.
        write['fullPeakFile']=self.state._spec_relative(self.fullPeakBox.GetValue()) if getattr(self, 'state', None) is not None else self.fullPeakBox.GetValue()
        write['fiddir']=self.outPathBox.GetValue()
        write['specPath']=self.specPathBox.GetValue()
        write['dim'] = str(self.state.spectral_dimensions)
        write['pseudo'] = int(self.state.pseudo_axis)
        write['peakShapeFitted']=1 if bool(getattr(self, 'peak_shape_fitted', False) or (getattr(self, 'store', None) is not None and self.store.metadata.get('peak_shape_determined'))) else 0
        write['peakFitCount']=int(getattr(self, 'peak_fit_count', ParseInt(self.deconParFile, 'peakFitCount', default=5)))
        write['peakFitLinkWidths']=1 if bool(getattr(self, 'peak_fit_link_widths', ParseInt(self.deconParFile, 'peakFitLinkWidths', default=1))) else 0
        write['pseudoIntensitiesStale']=1 if bool(getattr(self, 'pseudo_intensities_stale', False)) else 0
        write['downstreamAnalysis']=str(getattr(self, 'downstream_analysis', '') or '')
        write['pseudoSeriesInspected']=1 if bool(getattr(self, 'store', None) is not None and self.store.analysis.get('pseudo_series_reviewed')) else 0
        write['pickedPeaksChecked']=1 if bool(getattr(self, 'store', None) is not None and self.store.analysis.get('picked_peaks_reviewed')) else 0
        write['fittingResultsInspected']=1 if bool(getattr(self, 'store', None) is not None and self.store.analysis.get('fitting_results_reviewed')) else 0
        write['peakPickStale']=1 if bool(getattr(self, 'peak_pick_stale', False) or (getattr(self, 'store', None) is not None and self.store.analysis.get('peak_pick_stale'))) else 0
        pass
        
        if(self.dim>=1):
            write['sig1']=self.sig1Box.GetValue()
            write['voigt1'] = self.voigt1Box.GetValue()
            write['lor1'] = self.lorentz1Box.GetValue()
        if (self.dim>=2):
            write['sig2']=self.sig2Box.GetValue()
            write['voigt2'] = self.voigt2Box.GetValue()
            write['lor2'] = self.lorentz2Box.GetValue()
        if(self.dim>=3):
            write['sig3']=self.sig3Box.GetValue()
            write['voigt3'] = self.voigt3Box.GetValue()
            write['lor3'] = self.lorentz3Box.GetValue()
        if(self.dim==4):
            write['sig4']=self.sig4Box.GetValue()
            write['voigt4'] = self.voigt4Box.GetValue()
            write['lor4'] = self.lorentz4Box.GetValue()

        # Default UniDec values are implicit.  Persist only overrides so a
        # future change to an application default is not masked by a value the
        # user never changed.  If an override is changed back to its default,
        # remove the old key from the system file below.
        defaultable = {
            'thresh': self.threshBox.GetValue(),
            'ncpus': self.coreBox.GetValue(),
            'fac': self.facBox.GetValue(),
            'conv': self.convBox.GetValue(),
            'maxiter': self.maxiterBox.GetValue(),
        }
        default_keys = []
        for key, value in defaultable.items():
            if is_default_value(key, value):
                default_keys.append(key)
            else:
                write[key] = value
        write['FitRad']=self.fitRadBox.GetValue()
        # Shared with Frames/Pseudo3D.py Fitting window.
        write['3p_radF1']=self.fitF1Box.GetValue()
        write['3p_radF2']=self.fitF2Box.GetValue()
        write['symmode']=self.IntToBool(self.cb_grid.IsChecked())
        write['deconBore']=self.IntToBool(self.cb_decon3d.IsChecked())
        write['enhance']=self.IntToBool(self.cb_enhance.IsChecked())
        write['fitPhases']=self.IntToBool(self.cb_fitphases.IsChecked())

        self.state.set_parameter_file(self.deconParFile)
        self.state.sync_from_values(
            working_dir=self.dirBox.GetValue(),
            raw_path=self.outPathBox.GetValue(),
            input_file=self.infileBox.GetValue(),
            peak_file=self.peakBox.GetValue(),
            full_peak_file=self.fullPeakBox.GetValue(),
            spec_path=self.specPathBox.GetValue(),
            dimension=self.dim,
            pseudo_axis=bool(self.state.pseudo_axis),
            sym_mode=self.cb_grid.IsChecked(),
            decon_bore=self.cb_decon3d.IsChecked(),
        )

        target_path = os.path.join(self.dirBox.GetValue(), self.deconParFile)
        update_parameter_file(target_path, write, source_path=self.deconParFile)
        remove_parameter_keys(target_path, default_keys)
        # From this point onward save/load refer to this same canonical file.
        self.deconParFile = os.path.abspath(target_path)
        self.state.set_parameter_file(self.deconParFile)

    #####################################################
    # auxilary files for dealing with decon.
    def delfile(self,filey):
        if(os.path.exists(filey)):
            os.system('rm -rf '+filey)

    def cleanUp(self):
        # out/correlate.* peak-list outputs are deprecated; current decon
        # outputs are named from the input spectrum (<spectrum>.nD.list).
        self.delfile('out/diag.1')
        self.delfile('out/diag.2')
        self.delfile('out/diag.3')


    def run_decon(self, dec3dSet=False, dimProj=False, threshFac=1.0, caller='main', recon=False, input_override=None, dimension_override=None, peak_list_override=None, projection_labels=None):
        """Public entry point for external callers like the peaks window.

        Keeps the actual launch logic centralized in OnButtonDecon while giving
        other frames a stable API that does not depend on wx event wiring.
        """
        return self.OnButtonDecon(True, dec3dSet=dec3dSet, dimProj=dimProj, threshFac=threshFac, caller=caller, recon=recon, input_override=input_override, dimension_override=dimension_override, peak_list_override=peak_list_override, projection_labels=projection_labels)

    def _promote_pseudo2d_projection_peak_list(self, projection_peak_path):
        """Copy a pseudo2D projection result to the main spectrum's Full 1D list.

        The projection directory remains a complete record of the SpinUniDec
        calculation; the spectrum-associated ``<spectrum>.1D.list`` is the
        authoritative GUI peak list used by Full 1D and Projections.
        """
        topology = self._active_topology()
        if not (topology.spectral_dim_count == 1 and topology.has_pseudo_axis) or not projection_peak_path:
            return projection_peak_path
        source = os.path.abspath(str(projection_peak_path))
        main_spectrum = getattr(self, 'spectrumfile', '') or self._resolve_input_path(self.infileBox.GetValue())
        main_spectrum = os.path.abspath(str(main_spectrum)) if main_spectrum else ''
        if not main_spectrum:
            return source
        destination = self._peak_list_path_for_spectrum(main_spectrum, 1)
        if source != destination:
            if not os.path.exists(source):
                print('Cannot promote pseudo2D projection peak list; file not found:', source)
                return source
            shutil.copy2(source, destination)
            print('Copied pseudo2D projection peak list:', source, '->', destination)
        rel_destination = self.state._spec_relative(destination) if getattr(self, 'state', None) is not None else destination
        self.fullPeakBox.SetValue(rel_destination)
        if getattr(self, 'state', None) is not None:
            self.state.full_peak_file = rel_destination
            self.state.dirty = True
        return destination

    def _promote_pseudo3d_peakframe_reference_list(self, peak_path):
        """Copy a pseudo3D PeakFrame 2D list to the canonical reference name.

        PeakFrame decon/recon works on a transient 2D projection.  For a
        2-spectral + 1-pseudo project, keep that projection-owned result but
        also copy it into SpecPath/test.ft2.2D.list and make that copy the
        project's Reference 2D peak list.
        """
        topology = self._active_topology()
        if not (topology.spectral_dim_count == 2 and topology.has_pseudo_axis) or not peak_path:
            return peak_path

        source = os.path.abspath(str(peak_path))
        if not os.path.isfile(source):
            print('Cannot promote pseudo3D PeakFrame peak list; file not found:', source)
            return peak_path

        spec_dir = os.path.abspath(self.state.spec_dir())
        os.makedirs(spec_dir, exist_ok=True)
        destination = os.path.join(spec_dir, 'test.ft2.2D.list')
        if source != destination:
            shutil.copy2(source, destination)
            print('Copied pseudo3D reference peak list:', source, '->', destination)

        relative = self.state._spec_relative(destination)
        self.referencePeakBox.SetValue(relative)
        self.state.reference_peak_file = relative
        self.state.dirty = True

        peak_frame = getattr(self, 'peak_frame', None)
        if peak_frame is not None:
            try:
                peak_frame.peakfileBox.SetValue(destination)
            except Exception:
                pass
        return destination

    def _promote_projection_reference_peak_list(self, projection_peak_path):
        """Move a 3D PeakFrame decon peak list out of transient projections."""
        if self._active_topology().spectral_dim_count != 3 or not projection_peak_path:
            return projection_peak_path
        source = os.path.abspath(str(projection_peak_path))
        main_spectrum = getattr(self, 'spectrumfile', '') or self._resolve_input_path(self.infileBox.GetValue())
        main_spectrum = os.path.abspath(str(main_spectrum)) if main_spectrum else ''
        if not main_spectrum:
            return source
        destination = main_spectrum + '.2D.list'
        if source != destination:
            if not os.path.exists(source):
                print('Cannot promote projection peak list; file not found:', source)
                return source
            os.replace(source, destination)
            print('Promoted reference peak list:', source, '->', destination)
        rel_destination = self.state._spec_relative(destination) if getattr(self, 'state', None) is not None else destination
        self.referencePeakBox.SetValue(rel_destination)
        if getattr(self, 'state', None) is not None:
            self.state.reference_peak_file = rel_destination
            self.state.dirty = True
        peak_frame = getattr(self, 'peak_frame', None)
        if peak_frame is not None:
            try:
                peak_frame.peakfileBox.SetValue(destination)
            except Exception:
                pass
        return destination

    def _finish_decon_run(self, run_dim, run_dimProj, run_dec3d, run_ncpus, decset, rc=0, caller='main', expected_peak_path=None, projection_labels=None):
        if rc not in (0, None):
            if int(run_ncpus) > 1:
                wx.MessageBox('Error: parallel might not be supported... try lowering the CPU count to 1', 'Error', wx.OK | wx.ICON_ERROR)
            else:
                print('Decon finished with exit code', rc)
            self.update_project_lamps()
            if getattr(self, 'calcy', None) is not None: self.calcy.complete_decon_progress(False)
            return

        topology = self._active_topology()

        # Protocol2PFit is the pseudo2D restrained fitting/reconstruction
        # analysis.  Like Protocol3P, its products live in fit/*.out and
        # fit/*.dat.  It deliberately does NOT produce a projection .decon or
        # a new peak list: the Full 1D list is an input restraint.  Handle the
        # completion here before the ordinary pseudo2D projection-Decon path,
        # otherwise Recon+Fit incorrectly tries to ingest projection outputs.
        if str(decset.get('pseudo2DFit', '0')) == '1':
            fit_dir = self.get_fuda_dir()
            try:
                fit_outputs = [name for name in os.listdir(fit_dir) if name.endswith('.out')]
            except OSError:
                fit_outputs = []

            # Refresh an already-open pseudo2D fitting viewer if it exposes a
            # refresh hook.  Do not load deconvolution or peak-list products.
            try:
                if hasattr(self.parent, 'tabTwo') and self.parent.PageExists('Projections'):
                    viewer = getattr(self.parent.tabTwo, '_pseudo2d_fitting_window', None)
                    if viewer is not None:
                        refresh = getattr(viewer, 'refresh_results', None)
                        if callable(refresh):
                            refresh()
            except Exception as exc:
                print('Failed to refresh pseudo2D fitting results:', exc)

            self.update_project_lamps()
            if getattr(self, 'calcy', None) is not None:
                if fit_outputs:
                    self.calcy.append_text('\nPseudo2D restrained fitting complete. Results written to fit/.\n')
                    self.calcy.complete_decon_progress(True, 'Complete — pseudo2D fit outputs generated')
                else:
                    self.calcy.append_text('\nPseudo2D restrained fitting finished, but no fit/*.out results were found.\n')
                    self.calcy.complete_decon_progress(False, 'Incomplete — missing pseudo2D fit outputs')
            return

        # Protocol3P is a restrained fitting/reconstruction analysis.  Its
        # primary products are FUDA-compatible fit/*.out and fit/*.dat files;
        # it does not create a new decon peak list for the GUI to ingest.
        if str(decset.get('pseudo3D', '0')) == '1':
            missing = self.missing_pseudo3d_fit_peaks()
            if not missing and getattr(self, 'store', None) is not None:
                self.store.mark_pseudo_intensities_ready(
                    source='protocol3p_recon',
                    fit_directory=self.get_fuda_dir(),
                    reference_peak_count=len(self.get_reference_peaks() or []))
                self._mark_pseudo3d_recompute_complete()
                self._notify_analysis_changed()
            if getattr(self, 'calcy', None) is not None:
                if missing:
                    self.calcy.append_text(
                        '\n3P reconstruction/fitting finished, but fit results are still missing for: %s.\n'
                        % ', '.join(missing))
                    self.calcy.complete_decon_progress(False, 'Incomplete — missing 3P fit outputs')
                else:
                    self.calcy.append_text('\n3P reconstruction/fitting complete. Results written to fit/.\n')
                    self.calcy.complete_decon_progress(True, 'Complete — 3P fit outputs generated')
            self.update_project_lamps()
            return

        # Enhance produces only a calculated spectrum.  There is deliberately
        # no nD peak-list output, so bypass all normal peak-list loading and
        # analysis while still publishing the .decon spectrum for Show Calc.
        if str(decset.get('enhance', '0')) == '1':
            infile = str(decset.get('infile', '') or '')
            enhance_output = infile + '.decon'
            if topology.spectral_dim_count == 3 and run_dimProj == False and not topology.has_pseudo_axis:
                projection_dir = os.path.join(os.path.dirname(os.path.abspath(infile)), 'projections')
                try:
                    MakeProj3D(enhance_output, folder=projection_dir, OneD=False, clean=False, suffix='.decon')
                    self._cache_projection_folder(projection_dir, key_prefix='decon_projection')
                except Exception as exc:
                    print('Failed to project enhanced 3D spectrum:', exc)
            if not self._load_decon_outputs(infile, load_peaks=False):
                wx.MessageBox('Enhanced spectrum was not produced or could not be loaded.', 'Enhance', wx.OK | wx.ICON_ERROR)
                self.update_project_lamps()
                return
            self.parent.AddTabTwo(True, self)
            self.update_project_lamps()
            if getattr(self, 'calcy', None) is not None:
                self.calcy.append_text('\nEnhancement complete. Loaded: %s\n' % enhance_output)
                self.calcy.complete_decon_progress(True, 'Complete — enhanced spectrum generated')
            return

        # The output list dimensionality belongs to the decon job, not necessarily
        # to the main spectrum. PeakFrame can launch a 2D job from a 3D project.
        decon_dim = int(decset.get('dim', run_dim))
        peak_path = expected_peak_path or self._peak_list_path_for_spectrum(
            decset.get('infile', self.spectrumfile), decon_dim
        )
        if (caller == 'main' and run_dimProj == False and topology.spectral_dim_count == 1
                and topology.has_pseudo_axis):
            projection_source = str(decset.get('infile', '') or '')
            if not self._load_pseudo2d_projection_decon_outputs(projection_source):
                wx.MessageBox('Pseudo2D projection deconvolution output could not be loaded.',
                              'Deconvolution', wx.OK | wx.ICON_ERROR)
                self.update_project_lamps()
                return
            # Promote the projection job's list into the normal spectrum-owned
            # Full 1D workflow.  This makes one authoritative list drive both
            # the main-tab field and the Projections peak overlay.
            peak_path = self._promote_pseudo2d_projection_peak_list(peak_path)
            self.corrFile = peak_path
            if not self.load_full_peak_list(peak_path, quiet=True):
                print('Warning: pseudo2D Full 1D peak list is empty or could not be read:', peak_path)
            try:
                if hasattr(self.parent, 'tabTwo') and self.parent.PageExists('Projections'):
                    self.parent.tabTwo.draw_figure(keepaxes=True)
            except Exception as exc:
                print('Failed to refresh pseudo2D Projection after decon:', exc)
            self.parent.AddTabTwo(True, self)
            self.update_project_lamps()
            if getattr(self, 'calcy', None) is not None:
                self.calcy.complete_decon_progress(True)
                self.calcy.append_text('\nPseudo2D projection analysis complete.\n')
                self.calcy.set_status('Complete')
            return

        if caller == 'main' and run_dimProj == False and peak_path:
            # Main decon output is the full dimensionality-matched list.  Do
            # not overwrite the Reference 2D field.
            rel_peak = self.state._spec_relative(peak_path) if getattr(self, 'state', None) is not None else peak_path
            self.fullPeakBox.SetValue(rel_peak)
            if getattr(self, 'state', None) is not None:
                self.state.full_peak_file = rel_peak
            self.load_full_peak_list(peak_path, quiet=True)
            # A fresh main peak-picking run supersedes any prior review. A
            # restrained Recon/Fit consumes the already-checked list and must
            # not send Workflow back to Review picked peaks.
            is_recon = str(decset.get('recon', '0')) == '1'
            self.peak_pick_stale = False
            if getattr(self, 'store', None) is not None:
                self.store.analysis.pop('peak_pick_stale', None)
                if not is_recon:
                    self.store.invalidate_picked_peaks_review()
                if (is_recon and topology.spectral_dim_count == 2 and
                        not topology.has_pseudo_axis and str(decset.get('FIT', '0')).lower() in ('1', 'true')):
                    self.store.mark_fitting_results_ready(
                        source='2d_recon_fit', fit_directory=self.get_fuda_dir())
            self.OnButtonSave(True)
            self._notify_analysis_changed()

        if caller == 'main' and run_dimProj == False:  # main-window nD decon
            if topology.spectral_dim_count == 3:  # project the calculated 3D spectrum
                self.calcy.append_text('\nProjecting deconvolution results...\n') if getattr(self, 'calcy', None) is not None else None
                main_spectrum = getattr(self, 'spectrumfile', '') or self._resolve_input_path(decset['infile'])
                main_spectrum = os.path.abspath(main_spectrum or decset['infile'])
                decon_spectrum = main_spectrum + '.decon'
                projection_dir = os.path.join(
                    os.path.dirname(main_spectrum),
                    'projections',
                )
                MakeProj3D(
                    decon_spectrum,
                    folder=projection_dir,
                    OneD=False,
                    clean=False,
                    suffix='.decon',
                )
                # Cache only the newly calculated set; raw A.B.dat entries in
                # the same folder remain untouched and retain separate keys.
                self._cache_projection_folder(projection_dir, key_prefix='decon_projection')

            # Deprecated out/correlate.* peak-list paths are intentionally not
            # used.  The full list is <input spectrum>.nD.list.
            self.corrFile = peak_path
            self.OnButtonAnalyse(True)
        else:
            # PeakFrame launches its own decon job. For a 3D project this is
            # an independent 2D deconvolution of the displayed projection; for
            # a 2D project it is the displayed spectrum itself.
            # PeakFrame launches an independent 2D deconvolution of the
            # currently displayed projection. For 3D data this generated 2D
            # list becomes the shared reference list, so move it out of the
            # transient projections directory before reading it.
            if caller == 'peakframe' and run_dimProj != False and topology.spectral_dim_count == 3:
                peak_path = self._promote_projection_reference_peak_list(peak_path)
            elif caller == 'peakframe' and topology.spectral_dim_count == 2 and topology.has_pseudo_axis:
                peak_path = self._promote_pseudo3d_peakframe_reference_list(peak_path)
            self.corrFile = peak_path
            projection_source = str(decset.get('infile', '') or '').strip()
            projection_decon = projection_source + '.decon' if projection_source else ''
            if projection_decon and os.path.exists(projection_decon) and getattr(self, 'store', None) is not None:
                try:
                    dic_proj_dec, data_proj_dec = ng.pipe.read(projection_decon)
                    labels = tuple(projection_labels) if projection_labels and len(projection_labels) == 2 else (tuple(self.labb[2:0:-1]) if len(getattr(self, 'labb', [])) >= 3 else None)
                    if labels and len(labels) == 2:
                        for transpose in ('n', 'y'):
                            view = self._spectrum_view_payload(
                                dic_proj_dec, data_proj_dec, source=projection_decon,
                                labb=labels, transpose=transpose
                            )
                            self.store.save_view(('peakframe_decon', labels[0], labels[1], transpose), **view)
                    else:
                        print('Cannot determine 3D peak-frame projection labels for decon output:', projection_decon)
                except Exception as exc:
                    print('Failed to load peak-frame projection decon output:', projection_decon, exc)

            # The launch is asynchronous: consume the new peak list only once
            # decon has completed and any 3D reference-list promotion is done.
            peak_frame = getattr(self, 'peak_frame', None)
            if peak_frame is not None:
                try:
                    peak_frame.peakfileBox.SetValue(peak_path)
                    peak_frame.ReadPeakDecon()
                    peak_frame.draw_figure()
                except Exception as exc:
                    print('Failed to refresh PeakFrame decon peak list:', exc)

        if getattr(self, 'calcy', None) is not None:
            self.calcy.complete_decon_progress(True)

        self.parent.AddTabTwo(True, self)
        self.update_project_lamps()
        if getattr(self, 'calcy', None) is not None:
            self.calcy.append_text('\nAnalysis complete.\n')
            self.calcy.set_status('Complete')

    def OnButtonRecon(self, event):
        """Run a restricted reconstruction using the main Full nD peak list."""
        self._workflow_debug('OnButtonRecon ENTER: pseudo3d_topology=%s Fit?=%s Use2D=%s refs=%d' % (self._is_pseudo3d_topology(), self.cb_decback.IsChecked(), self.cb_decon3d.IsChecked(), len(self.get_reference_peaks() or [])))
        # 3P recon is restrained by the shared 2D reference peak list; ordinary
        # nD recon continues to use the Full nD list.
        peak_value = self.referencePeakBox.GetValue() if self._is_pseudo3d_topology() else self.fullPeakBox.GetValue()
        peak_path = self._resolve_spec_file(peak_value)
        if not peak_path or not os.path.isfile(peak_path):
            label = 'Reference 2D' if self._is_pseudo3d_topology() else 'Full nD'
            wx.MessageBox('A valid %s peak list is required for Recon.' % label, 'Recon', wx.OK | wx.ICON_WARNING)
            return
        return self.run_decon(caller='main', recon=True)

    def OnButtonDecon(self,event,dec3dSet=False,dimProj=False,threshFac=1.0,caller='main',recon=False, input_override=None, dimension_override=None, peak_list_override=None, projection_labels=None):
        if(self.READ==0):
            self.OnButtonRead(True) #try to read in...
            if(self.READ==0):
                wx.MessageBox('No spectrum is loaded. Please load a spectrum before running deconvolution.', 'Deconvolution', wx.OK | wx.ICON_WARNING)
                return


        
        topology = self._selected_topology()
        spectral_dim_count = topology.spectral_dim_count
        pseudo2d_fit = bool(recon and caller == 'main' and input_override is None and
                            topology.spectral_dim_count == 1 and topology.has_pseudo_axis and
                            self.cb_decon3d.IsChecked() and self.cb_decback.IsChecked())
        if getattr(self, 'state', None) is not None:
            self.state.spectral_dimensions = spectral_dim_count
            self.state.pseudo_axis = topology.has_pseudo_axis
        self.thresh=float(self.threshBox.GetValue())*threshFac


        decset={}

        symmy=self.cb_grid.IsChecked()*1.    #symmetric mode?

        if recon:
            # Restricted nD reconstruction is independent of the 2D bore-mode
            # checkbox.  In particular, a checked "Use 2D peaklist" must not
            # alter the restrained full-dimensional calculation.
            dec3d=0
        elif (caller == 'peakframe' and topology.spectral_dim_count == 2 and
              not topology.has_pseudo_axis and dimension_override == 2):
            # PeakFrame on a true physical 2D dataset is an ordinary,
            # unrestrained 2D peak-picking calculation.  Do not inherit the
            # main-window "Use 2D peaklist" checkbox: that control belongs to
            # restrained/reconstruction workflows, not PeakFrame Decon.
            dec3d=0
        elif(dec3dSet==False and dimProj==False):
            dec3d=self.cb_decon3d.IsChecked()*1. #bore mode?
            if dec3d == 1:
                self.OnButtonReadPeak(None)
        else:
            dec3d=1


        if dimension_override is not None:
            # Explicit analysis targets (notably PeakFrame) own their job
            # dimensionality.  Parent-dataset topology must not reinterpret a
            # displayed 2D spectrum/projection as a 3D/4D job.
            dimVal = str(int(dimension_override))
        elif(dimProj!=False):
            # Projection jobs are ordinary spectral jobs: their protocol
            # dimensionality is the number of projected spectral axes.
            dimVal=str(len(dimProj))
        elif pseudo2d_fit:
            # SpinUniDec Protocol2PFit receives the original physical pseudo2D
            # NMRPipe file, so its external protocol dimension is two.
            dimVal=str(topology.physical_dim_count)
        elif (topology.spectral_dim_count == 2 and
              topology.has_pseudo_axis and
              topology.physical_dim_count == 3):
            # spinUnidec dispatches Protocol3P through ``dim == 3`` plus the
            # explicit pseudo3D selector.  ``spectral_dim_count`` remains 2
            # everywhere in the GUI/project model; only this external protocol
            # field uses the physical dimension count.
            dimVal=str(topology.physical_dim_count)
        else:
            dimVal=str(spectral_dim_count)
        decset['dim']=dimVal


        self.ncpus=int(self.coreBox.GetValue())
        #self._sync_peakfile_box(self.spectrumfile if getattr(self, 'spectrumfile', '') else self._resolve_input_path(self.infileBox.GetValue()), self.dim)
        projection_peak_path = None
        if input_override:
            projection_input = os.path.abspath(str(input_override))
            projection_peak_path = projection_input + '.2D.list' if int(dimVal) == 2 else None
        elif dimProj != False:
            projection_input = self._resolve_projection_input_path(self.labb[2], self.labb[1],decon=False) if hasattr(self, 'labb') and len(getattr(self, 'labb', [])) >= 3 else None
            if projection_input:
                projection_peak_path = projection_input + '.2D.list'
        if int(self.ncpus) < 2:
            specstr=self.deconBin
        else:
            specstr=self.paraDeconBin

        decset['ncpus']=self.coreBox.GetValue()
        decset['rand']=str(7)
        decset['maxIter']=self.maxiterBox.GetValue()

        if(dimProj!=False):
            decset['maxIter']=10000 #default value for 2D
        else:
            decset['maxIter']=self.maxiterBox.GetValue()
    
        decset['conv']=self.convBox.GetValue()

        if spectral_dim_count == 1:
            decset['uSTA']=str(self.uSTA)
            decset['baseFile']=str('False')
            





        if input_override:
            target_spec = os.path.abspath(str(input_override))
            dic, data = ng.pipe.read(target_spec)
            max_proj = numpy.fabs(data[numpy.unravel_index(numpy.argmax(numpy.fabs(data)), data.shape)] * self.thresh)
            decset['infile'] = target_spec
            decset['dmax'] = str(max_proj)
            # Preserve the established 3D/4D projection line-shape mapping
            # when dimProj describes that view.  A 2D+pseudo PeakFrame target
            # has the two spectral axes directly in sig1/sig2.
            if dimProj != False:
                decset['sig2'] = self.sig3Box.GetValue()
                decset['voigt2'] = self.voigt3Box.GetValue()
                decset['lor2'] = self.lorentz3Box.GetValue()
                decset['sig1'] = self.sig2Box.GetValue()
                decset['voigt1'] = self.voigt2Box.GetValue()
                decset['lor1'] = self.lorentz2Box.GetValue()
            else:
                decset['sig1'] = self.sig1Box.GetValue()
                decset['voigt1'] = self.voigt1Box.GetValue()
                decset['lor1'] = self.lorentz1Box.GetValue()
                if int(dimVal) >= 2:
                    decset['sig2'] = self.sig2Box.GetValue()
                    decset['voigt2'] = self.voigt2Box.GetValue()
                    decset['lor2'] = self.lorentz2Box.GetValue()

        elif(dimProj!=False):

            decset['sig2']=self.sig3Box.GetValue()            
            decset['voigt2']=self.voigt3Box.GetValue()
            decset['lor2']=self.lorentz3Box.GetValue()

            pseudo_spec = self._resolve_projection_input_path(self.labb[2], self.labb[1],decon=False)
            if not pseudo_spec:
                raise FileNotFoundError(f'Cannot find projection file for {self.labb[2]}.{self.labb[1]}')
            dic, data = ng.pipe.read(pseudo_spec)
            max_proj = numpy.fabs(data[numpy.unravel_index(numpy.argmax(numpy.fabs(data)), data.shape)]*self.thresh)
            decset['infile']=pseudo_spec
            decset['dmax']=str(max_proj)
            if(spectral_dim_count >= 2):
                decset['sig1']=self.sig2Box.GetValue()
                decset['voigt1']=self.voigt2Box.GetValue()
                decset['lor1']=self.lorentz2Box.GetValue()

            
        elif (topology.spectral_dim_count == 1 and topology.has_pseudo_axis):
            decset['sig1'] = self.sig1Box.GetValue()
            decset['voigt1'] = self.voigt1Box.GetValue()
            decset['lor1'] = self.lorentz1Box.GetValue()
            if pseudo2d_fit:
                # Restrained pseudo2D fitting must see every pseudo slice, not
                # the summed projection used by ordinary pseudo2D Decon.
                decset['infile'] = self._resolve_input_path(self.infileBox.GetValue())
                decset['dmax'] = str(self.dmax * self.thresh)
                decset['pseudo2DFit'] = '1'
                decset['FIT'] = '1'
                decset['FitPhase'] = '1' if self.cb_fitphases.IsChecked() else '0'
                decset['dec3d'] = '0'
                decset['peakList'] = self._resolve_spec_file(self.fullPeakBox.GetValue())
                fitf1_value = self.fitF1Box.GetValue().strip()
                if fitf1_value:
                    decset['FitF1'] = fitf1_value
            else:
                projection = self.get_pseudo2d_projection_data(ensure_file=True)
                if projection is None:
                    raise FileNotFoundError('Cannot create pseudo2D spectral projection')
                pseudo_spec = projection['path']
                data = numpy.asarray(projection['data'])
                max_proj = numpy.fabs(data[numpy.argmax(numpy.fabs(data))]) * self.thresh
                decset['infile'] = pseudo_spec
                decset['dmax'] = str(max_proj)
                # Ask the adapted SpinUniDec binary for the compact pure-1D list.
                decset['pseudo2DOutput'] = '1'

        else:  #normal mode

            decset['infile']=self._resolve_input_path(self.infileBox.GetValue())
            decset['dmax']=str(self.dmax*self.thresh)

            decset['sig1']=self.sig1Box.GetValue()            
            decset['voigt1']=self.voigt1Box.GetValue()
            decset['lor1']=self.lorentz1Box.GetValue()            
            if(spectral_dim_count >= 2):
                decset['sig2']=self.sig2Box.GetValue()
                decset['voigt2']=self.voigt2Box.GetValue()
                decset['lor2']=self.lorentz2Box.GetValue()
            if(spectral_dim_count >= 3):
                decset['sig3']=self.sig3Box.GetValue()
                decset['voigt3']=self.voigt3Box.GetValue()
                decset['lor3']=self.lorentz3Box.GetValue()
            if(spectral_dim_count == 4):
                decset['sig4']=self.sig4Box.GetValue()
                decset['voigt4']=self.voigt4Box.GetValue()
                decset['lor4']=self.lorentz4Box.GetValue()
                
        decset['fac']=self.facBox.GetValue()
        # Enhance is intentionally limited to 1D-3D.  4D retains the existing
        # protocol unchanged even if a saved project has the checkbox enabled.
        decset['enhance']='1' if (spectral_dim_count <= 3 and self.cb_enhance.IsChecked() and not recon) else '0'

        if recon:
            decset['recon']='1'
            decset['peakList'] = (os.path.abspath(str(peak_list_override)) if peak_list_override
                                  else self._resolve_spec_file(self.referencePeakBox.GetValue() if self._is_pseudo3d_topology() else self.fullPeakBox.GetValue()))
            decset['dec3d']='0'
            # Protocol3P applies only to a main-spectrum job.  PeakFrame recon
            # is always an ordinary 2D reconstruction of its displayed target.
            if self._is_pseudo3d_topology() and caller == 'main' and input_override is None:
                # Silent protocol selector: physical 3D data with exactly two
                # spectral axes and one real axis are handled by Protocol3P.
                decset['pseudo3D']='1'
                # Protocol3P is intrinsically a restrained pseudo3D fitting
                # run and the C++ implementation requires FIT=1.  Do not make
                # this depend on the separate GUI background-fit checkbox.
                decset['FIT']='1'
                fitf1_value = self.fitF1Box.GetValue().strip()
                fitf2_value = self.fitF2Box.GetValue().strip()
                if fitf1_value:
                    decset['FitF1'] = fitf1_value
                if fitf2_value:
                    decset['FitF2'] = fitf2_value
        elif(dec3d and dec3dSet): #legacy PeakFrame restricted 2D path
            decset['peakList']='out/list.tmp'
            decset['dec3d']=1
            decset['dmax']=-1
        elif(dec3d and dec3dSet==False):
            decset['peakList'] = projection_peak_path if dimProj != False and projection_peak_path else self._resolve_spec_file(self.referencePeakBox.GetValue())
            decset['dec3d']=1


        
        if(spectral_dim_count >= 3 and not topology.has_pseudo_axis):
            decset['symmy']=str(symmy)
            decset['dec3d']='0' if recon else str(dec3d)

            
        physical_2d_peakframe_decon = bool(
            caller == 'peakframe' and not recon and
            topology.spectral_dim_count == 2 and not topology.has_pseudo_axis and
            dimension_override == 2)
        if(self.cb_decback.IsChecked() and decset.get('enhance') != '1' and
           not physical_2d_peakframe_decon):
            # The physical-2D PeakFrame Decon button is a peak-picking action,
            # not the main workflow's Recon/Fit action.  Keep it independent of
            # a stale/checked main-window Fit control.
            decset['FIT']=True
            fitrad_value = self.fitRadBox.GetValue().strip()
            if fitrad_value:
                decset['FitRad']=fitrad_value
            if spectral_dim_count == 2:
                fitf1_value = self.fitF1Box.GetValue().strip()
                fitf2_value = self.fitF2Box.GetValue().strip()
                if fitf1_value:
                    decset['FitF1'] = fitf1_value
                if fitf2_value:
                    decset['FitF2'] = fitf2_value
            
        #specstr2+='\tsquash '+str(self.conv)
        #fields=specstr2.split('\t')
        #print(fields)
        service = getattr(self.parent, "decon_service", None)
        if service is None:
            from spinDecon.project.decon_service import DeconService
            service = DeconService()

        service.write_init_dict(decset)

        self.cleanUp()

        specstr += ' decon.init'
        print('Executing:', specstr, flush=True)
        run_dim = int(decset['dim'])
        run_dec3d = dec3d
        run_dimProj = dimProj
        run_ncpus = self.ncpus
        # Preserve the exact expected output path across the asynchronous run.
        # Projection jobs are always 2D even when the owning project is 3D.
        expected_peak_path = projection_peak_path if dimProj != False else self._peak_list_path_for_spectrum(decset.get('infile', self.spectrumfile), int(decset['dim']))

        # The C++ process writes one convergence value per ApplyIter to this
        # sidecar file.  Remove stale data before launch so the live plot can
        # never briefly show a previous run while the new process starts.
        convergence_file = decset['infile'] + '.conv'
        try:
            os.remove(convergence_file)
        except FileNotFoundError:
            pass
        except OSError as exc:
            print('Warning: could not remove old convergence file:', exc, flush=True)

        self.calcy = run_command_with_output(
            specstr,
            parent=self,
            title='decon stdout',
            on_finish=lambda rc=0: self._finish_decon_run(run_dim, run_dimProj, run_dec3d, run_ncpus, decset, rc, caller=caller, expected_peak_path=expected_peak_path, projection_labels=projection_labels),
            final=False,
            convergence_file=convergence_file,
            decon_profile=decset,
        )


    #Analyse results from spinDecon
    def OnButtonAnalyse(self,event):
        # Loading a calculated spectrum is dependency-aware just like the
        # Workflow actions.  In particular, physical 3D/4D deconvolution
        # loading derives calculated slices from reference peak indices, so a
        # cold project must materialise spectrum -> reference before .decon.
        if not self.ensure_workflow_spectrum_loaded():
            print('Cannot load deconvolution before the main spectrum is loaded.')
            return False
        topology = self._active_topology()
        spectral_dim_count = topology.spectral_dim_count
        if spectral_dim_count >= 3 and not topology.has_pseudo_axis:
            if not self.ensure_reference_peak_list_loaded():
                print('Cannot load deconvolution before the reference peak list is loaded.')
                return False
        self.DECON=0         #decon flag
        self.pkSlice1Ddec=[] #1D slices
        infile=self._resolve_input_path(self.infileBox.GetValue())

        try:
            poll = self.calcy.poll()
            if poll == None:
                print('Calculation still running in background')
        except:
            pass

        #bvbbbbbb
        #self.corrFile = self._sync_peakfile_box(self.spectrumfile if getattr(self, 'spectrumfile', '') else self._resolve_input_path(self.infileBox.GetValue()), self.dim) or self.corrFile
        if getattr(self, 'calcy', None) is not None:
            self.calcy.append_text(f'Loading deconvolution results: {self.corrFile}\n')
        if spectral_dim_count == 1:
            if topology.has_pseudo_axis:
                projection = self.get_pseudo2d_projection_data(ensure_file=False)
                projection_path = projection.get('path') if projection else ''
                if not projection_path or not self._load_pseudo2d_projection_decon_outputs(projection_path):
                    print('Pseudo2D projection deconvolution file does not exist:',
                          (str(projection_path) + '.decon') if projection_path else '')
                    return
            elif not self._load_decon_outputs(infile):
                print('Deconvolution file does not exist.')
                return

        if(spectral_dim_count == 2):
            if topology.has_pseudo_axis:
                # Pseudo3D currently analyses a 2D projection in peakFrame.
                # There is no general full-spectrum deconvolution product for
                # 2 spectral + 1 pseudo datasets, so do not attempt to read
                # <spectrum>.decon here.  Peak/connectivity results can still
                # be loaded below.
                print('Pseudo3D analysis: skipping deconvolved spectrum read.')
            elif not self._load_decon_outputs(infile):
                print('Deconvolution file does not exist.')
                return

        if(spectral_dim_count == 3):
            if not topology.has_pseudo_axis:
                if not self._load_decon_outputs(infile):
                    print('Deconvolution file does not exist.')
                    return
            else:
                infile=self.infileBox.GetValue()
                # if(os.path.exists(infile+'.decon')==1):
                #     self.dicdec,self.datadec=ng.pipe.read(infile+'.decon')
                #     self.DECON=1
                #     if(self.datadec.shape!=self.data.shape):
                #         print('deconvolved spectrum is a different shape.')
                #         print(self.datadec.shape, self.data.shape)
                #         print('recalculate the deconvolution')
                #         self.DECON=0
                #         #numpy.delete(self.datadec)
                #         #numpy.delete(self.dicdec)
                #         return
                #     if(self.DECON==1):
                #         print('Deconvolved spectrum in memory:',infile)
                # else:
                #     print('Deconvolution file does not exist.')
                #     return
                # self.pkSlice1Ddec=[] #1D slices
                # for pkl in range(len(self.peak)):
                #     ptC=self.pkIdx[pkl][0]
                #     ptH=self.pkIdx[pkl][1]
                #     self.pkSlice1Ddec.append(self.datadec[:,ptC,ptH])

        else:
            pass

        if spectral_dim_count == 4:
            if not self._load_decon_outputs(infile):
                print('Deconvolution file does not exist.')
                return

        # Calculated peak output is a peak list, not connectivity data.  Store
        # it using the same canonical record schema as every other peak list.
        if spectral_dim_count == 1 and topology.has_pseudo_axis:
            self._load_pseudo2d_projection_peaks(self.corrFile)
        else:
            self._load_decon_peak_list(self.corrFile)

        if getattr(self, "store", None) is not None:
            payload = self.store.get_peak_list("decon")
            payload.update(pkIdx=self.pkIdx, pkSlice1D=self.pkSlice1D,
                           pkSlice1Ddec=self.pkSlice1Ddec, Grps=self.Grps)

        # Refresh any open views so they can enable peak/decon controls now
        # that the shared datastore has been populated.
        try:
            if hasattr(self.parent, "tabFive") and self.parent.PageExists('Full 3D'):
                self.parent.tabFive._refresh_control_availability()
                self.parent.tabFive.draw_figure(keepaxes=False)
        except Exception as exc:
            print('Failed to refresh Full3D after analyse:', exc)
        try:
            if hasattr(self.parent, "tabTwo") and self.parent.PageExists('Projections'):
                self.parent.tabTwo.draw_figure(keepaxes=False)
        except Exception as exc:
            print('Failed to refresh Projection after analyse:', exc)



    #########################################################
    # Handle aliasing in peak lists to current spectrum
    def alias(self,peak,ppm,dim):
        if(dim==0):
            dd=numpy.fabs(self.index0[1]-self.index0[0] )
            dmax=self.uc0max
            dmin=self.uc0min
            vals=self.index0
        elif(dim==1):
            dd=numpy.fabs(self.index1[1]-self.index1[0] )
            dmax=self.uc1max
            dmin=self.uc1min
            vals=self.index1
        elif(dim==2):
            dd=numpy.fabs(self.index2[1]-self.index2[0] )
            dmax=self.uc2max
            dmin=self.uc2min
            vals=self.index2
        elif(dim==3):
            dd=numpy.fabs(self.index3[1]-self.index3[0])
            dmax=self.uc3max
            dmin=self.uc3min
            vals=self.index3
        else:
            return -1
        #print('alias start:',dim, ppm, dmax, dmin, dd)
        while(ppm>=dmax):
            # print(dim, ppm, dmax, dmin, dd)
            ppm-=(dmax-dmin+dd)
        while(ppm<=dmin):
            # print(dim, ppm, dmax, dmin, dd)
            ppm+=(dmax-dmin+dd)

        i=findnear_index(ppm,vals)
        #print('alias finish:',dim, ppm, dmax, dmin, dd)
        if(dim==0):
            peak.indexI=i
            peak.ppmI=vals[i]
        elif(dim==1):
            peak.indexJ=i
            peak.ppmJ=vals[i]
        elif(dim==2):
            peak.indexK=i
            peak.ppmK=vals[i]
        elif(dim==3):
            peak.indexL=i
            peak.ppmL=vals[i]

    #look for projected data file
    def SetProjectedData(self,inv=False):
        #self.projectedData = numpy.sum(numpy.sum(self.data, axis=0), axis=0)
        projpath=self._spec_output_dir()

        if(inv==False):
            # For physical 3p, choose the direct (last spectral) axis rather
            # than blindly using labb[-1], which may itself be the pseudoaxis.
            spectral_axes = self._pseudo3d_spectral_axes()
            proj_label = spectral_axes[-1][1] if spectral_axes else self.labb[-1]
            projname=os.path.join(projpath,'projections1D',str(proj_label)+'.dat')
        else: #we are 2D
            # For a normal 2D spectrum, data axis 0 is the indirect
            # dimension and data axis 1 is the direct dimension.  The
            # projection plot uses the direct dimension as its X axis, so
            # load the direct-dimension projection here.
            projname=os.path.join(projpath,'projections',self.labb[1]+'.dat')
        print('Reading ',projname)
        if(os.path.exists(projname)):
            doc,self.projectedData=ng.pipe.read(projname)
            return
        print('No 1D projection found. Need to create it!')
        self.projectedData=[]

    ####################################
    # Read in main raw data file.
    def makeinp(self,indir,infile):
        if os.path.isabs(infile) or os.path.exists(infile):
            infile = infile
        else:
            infile = os.path.join(indir, infile)

        if(os.path.exists('out')==0):
           os.mkdir('out')
        if(os.path.exists('out/slice2d')==0):
           os.mkdir('out/slice2d')
        if(os.path.exists('out/fit')==0):
           os.mkdir('out/fit')

        #os.system('ucsf '+indir+'/new.ft2')
        print('Reading spectrum:', infile)

        #quick read to find dimensionality
        if '.ft' in infile:
            self.dic,self.data = ng.pipe.read_lowmem(infile)
        elif '.ucsf' in infile:
            self.dic, self.data = ng.sparky.read(infile)
        else:
            try:
                self.dic,self.data = ng.pipe.read_lowmem(infile)
            except:
                print("Failed read in")
                
        # Loading an ndarray must never redefine scientific dimensionality.
        # ProjectState owns the spectral/pseudo contract; the array contributes
        # only its physical dimensionality and (where available) axis labels.
        physical_dim = self.data.ndim
        _original_physical_dim = physical_dim
        _original_axis_uc = {}
        if '.ft' in infile:
            for _axis in range(_original_physical_dim):
                try:
                    _original_axis_uc[_axis] = ng.pipe.make_uc(self.dic, self.data, dim=_axis)
                except Exception:
                    pass
        _surviving_original_axes = list(range(_original_physical_dim))
        pipe_labels = []
        if isinstance(getattr(self, 'dic', None), dict):
            # FDF1LABEL/FDF2LABEL/... are dimension-number labels, not NumPy
            # physical-axis order.  Reorder them through FDDIMORDER before
            # recording physical-axis identity in ProjectState.
            _dim_labels = []
            for _n in range(1, 5):
                _lab = self.dic.get('FDF%dLABEL' % _n, '')
                _dim_labels.append('' if _lab is None else str(_lab))
            _order = self.dic.get('FDDIMORDER', ())
            try:
                _physical_order = [int(_order[physical_dim - 1 - i]) - 1
                                   for i in range(physical_dim)]
                pipe_labels = [_dim_labels[i] for i in _physical_order]
            except (TypeError, ValueError, IndexError):
                pipe_labels = _dim_labels[:physical_dim]

        real_axis_labels = ('time_T2', 'ID', 'ncyc', 'ncyc_cp', 'gzlvl5', 'gzlvl1', 'usta')
        state = getattr(self, 'state', None)

        # NMRPipe can expose a lower-dimensional spectrum through a file/header
        # with one or more singleton storage axes.  Those axes are storage
        # artefacts, not scientific dimensions.  If the selected canonical
        # topology already tells us the expected physical dimensionality, drop
        # only the exact number of surplus singleton axes.  Never infer a new
        # spectral dimension from ndarray.ndim.
        expected_physical = 0
        if state is not None and state.spectral_dimensions >= 1:
            expected_physical = state.spectral_dimensions + int(bool(state.pseudo_axis))
        if expected_physical and physical_dim > expected_physical:
            surplus = physical_dim - expected_physical
            singleton_axes = [i for i, n in enumerate(self.data.shape) if int(n) == 1]
            if len(singleton_axes) >= surplus:
                drop_axes = singleton_axes[:surplus]
                old_shape = tuple(self.data.shape)
                indexer = tuple(0 if i in drop_axes else slice(None) for i in range(physical_dim))
                self.data = self.data[indexer]
                pipe_labels = [label for i, label in enumerate(pipe_labels) if i not in drop_axes]
                _surviving_original_axes = [i for i in _surviving_original_axes if i not in drop_axes]
                physical_dim = self.data.ndim
                pass

        pass
        if state is not None:
            # Parameter/project state is authoritative.  The one permitted
            # legacy repair is canonicalize_loaded_dimensions(), which converts
            # historical physical-count pseudo projects at this load boundary.
            if state.spectral_dimensions < 1 and int(getattr(self, 'dim', 0) or 0) > 0:
                state.spectral_dimensions = int(self.dim)
            state.pseudo_axis = bool(self.state.pseudo_axis)
            pass
            try:
                migrated = state.canonicalize_loaded_dimensions(
                    physical_dim, pipe_labels, real_axis_labels=real_axis_labels
                )
            except ValueError as exc:
                print('Data shapes do not match. Problem!')
                print(str(exc))
                return
            self.dim = state.spectral_dimensions
            pass
            if migrated:
                print('Migrated legacy pseudo dimensionality during load:',
                      'physical shape:', self.data.shape, 'labels:', pipe_labels,
                      '-> spectral dimension', self.dim)
        else:
            # Isolated legacy callers have no canonical state.  Preserve their
            # selected spectral count, but validate it against the physical data.
            selected_spectral = int(getattr(self, 'dim', 0) or 1)
            expected_physical = selected_spectral + int(bool(getattr(self, 'pseudo', False)))
            if physical_dim != expected_physical:
                print('Data shapes do not match. Problem!')
                print('spectral dimension:', selected_spectral, 'physical shape:', self.data.shape,
                      'labels:', pipe_labels)
                return
            self.dim = selected_spectral
        
        # Enforce the canonical contract at the full-spectrum load boundary.
        # Derived arrays are validated in their own viewers instead.
        state = getattr(self, 'state', None)
        if state is not None and int(getattr(state, 'spectral_dimensions', 0) or 0) > 0:
            from spinDecon.domain.dimensions.guard import assert_full_dataset_contract
            try:
                assert_full_dataset_contract(state.topology(), self.data, where='deconFrame.makeinp')
            except ValueError as exc:
                print('Data shapes do not match. Problem!')
                print(str(exc))
                return

        topology = self._active_topology()
        spectral_dim_count = topology.spectral_dim_count

        # Keep the selector aligned with canonical state; dimensional setup
        # below must branch on topology, not the transitional ``self.dim`` alias.
        pass
        if int(getattr(self, 'dim', 0) or 0) != spectral_dim_count:
            self.dim = spectral_dim_count
            self.SetDim()
        if (os.path.exists(self.deconParFile)==0):
            self.SetDim()
            
        if spectral_dim_count == 1:
            if topology.has_pseudo_axis:
                Size = self.data.shape
                self.specsize = Size
                self.dmax = numpy.max(self.data)
                pseudo_index = int(getattr(topology, 'pseudo_physical_index', 0) or 0)
                spectral_indices = [i for i in range(len(Size)) if i != pseudo_index]
                if len(spectral_indices) != 1:
                    raise RuntimeError('Pseudo2D requires exactly one spectral physical axis')
                spectral_index = spectral_indices[0]
                spectral_original_index = _surviving_original_axes[spectral_index]
                spectral_uc = _original_axis_uc.get(spectral_original_index)
                if spectral_uc is None and '.ucsf' in infile:
                    spectral_uc = ng.sparky.make_uc(self.dic, self.data, dim=spectral_index)
                if spectral_uc is None:
                    raise RuntimeError('Pseudo2D could not construct spectral unit converter')
                spectral_scale = numpy.asarray([spectral_uc.ppm(i) for i in range(Size[spectral_index])])
                pseudo_scale = numpy.arange(Size[pseudo_index], dtype=float)
                physical_scales = [None, None]
                physical_scales[pseudo_index] = pseudo_scale
                physical_scales[spectral_index] = spectral_scale
                self.index0 = numpy.asarray(physical_scales[0])
                self.index1 = numpy.asarray(physical_scales[1])
                self.uc0max, self.uc0min = float(self.index0[0]), float(self.index0[-1])
                self.uc1max, self.uc1min = float(self.index1[0]), float(self.index1[-1])
                self.labb = tuple(pipe_labels)
                self.YY, self.XX = numpy.meshgrid(self.index1, self.index0)
                pass
            else:
                self.dic,self.data=ng.pipe.read(infile)

                uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)

                ord=self.dic['FDDIMORDER']
                lab1=self.dic['FDF2LABEL']

                lab=lab1
                self.labb=lab
                self.dmax=numpy.max(self.data)  #get max intensity in spectrum
                Size=self.data.shape
                self.specsize=Size

                self.uc0max=uc0.ppm(0)
                self.uc0min=uc0.ppm(Size[0]-1)


                self.index0=[]#make index of carbon chemical shifts for index 0
                for i in range((Size[0])):
                    self.index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))

                self.index0=numpy.array(self.index0)


                self.XX = self.index0


        if spectral_dim_count == 2:

            if topology.has_pseudo_axis:
                self.dic,self.data=ng.pipe.read(infile)
                # uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)
                # uc1 = ng.pipe.make_uc(self.dic,self.data,dim=1)
                # uc2 = ng.pipe.make_uc(self.dic,self.data,dim=2)

                self.uc0 = ng.pipe.make_uc(self.dic, self.data, dim=0)
                x0,x1=self.uc0.ppm_limits()
                self.uc0.ppms_scale=numpy.linspace(x0, x1, int(self.uc0._size))
                self.uc1 = ng.pipe.make_uc(self.dic, self.data, dim=1)
                x0,x1=self.uc1.ppm_limits()
                self.uc1.ppms_scale=numpy.linspace(x0, x1, int(self.uc1._size))
                self.uc2 = ng.pipe.make_uc(self.dic, self.data, dim=2)
                x0,x1=self.uc2.ppm_limits()
                self.uc2.ppms_scale=numpy.linspace(x0, x1, int(self.uc2._size))
                
                ord=self.dic['FDDIMORDER']
                lab1=self.dic['FDF1LABEL']
                lab2=self.dic['FDF2LABEL']
                lab3=self.dic['FDF3LABEL']
                lab=lab1,lab2,lab3

  
                self.labb=lab[int(ord[2])-1],lab[int(ord[1])-1],lab[int(ord[0])-1]
                self.dmax=numpy.max(numpy.fabs(self.data))  #get max intensity in spectrum
                Size=self.data.shape
                self.specsize=Size

                self.SetProjectedData()

                self.uc0max=self.uc0.ppm(0)
                self.uc0min=self.uc0.ppm(Size[0]-1)
                self.uc1max=self.uc1.ppm(0)
                self.uc1min=self.uc1.ppm(Size[1]-1)
                self.uc2max=self.uc2.ppm(0)
                self.uc2min=self.uc2.ppm(Size[2]-1)

                #print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
                #print "Labels: ",self.labb
                #print "dimension 0 limits (ppm): ", self.uc0min, self.uc0max  #carbon
                #print "dimension 1 limits (ppm): ", self.uc1min, self.uc1max  #direct
                #print "dimension 2 limits (ppm): ", self.uc2min, self.uc2max  #direct
                #print 'Maximum Intensity:',self.dmax
                self.index0=[]#make index of carbon chemical shifts for index 0
                for i in range((Size[0])):
                    self.index0.append((self.uc0.ppm(0)-i*(-self.uc0.ppm(Size[0]-1)+self.uc0.ppm(0))/(Size[0]-1)))
                self.index1=[]#make index of carbon chemical shifts for index 1
                for i in range((Size[1])):
                    self.index1.append((self.uc1.ppm(0)-i*(-self.uc1.ppm(Size[1]-1)+self.uc1.ppm(0))/(Size[1]-1)))
                self.index2=[]#make index of carbon chemical shifts for index 2
                for i in range((Size[2])):
                    self.index2.append((self.uc2.ppm(0)-i*(-self.uc2.ppm(Size[2]-1)+self.uc2.ppm(0))/(Size[2]-1)))
                self.index0=numpy.array(self.index0)
                self.index1=numpy.array(self.index1)
                self.index2=numpy.array(self.index2)

                self.YY,self.XX,self.ZZ=numpy.meshgrid(self.index1,self.index0,self.index2)


            else:
            
                # self.dic,self.data=ng.pipe.read(infile)
                if '.ft' in infile:
                    self.dic,self.data = ng.pipe.read_lowmem(infile)
                    self.uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)
                    self.uc1 = ng.pipe.make_uc(self.dic,self.data,dim=1)
                elif '.ucsf' in infile:
                    self.dic, self.data = ng.sparky.read(infile)
                    self.uc0 = ng.sparky.make_uc(self.dic,self.data,dim=0)
                    self.uc1 = ng.sparky.make_uc(self.dic,self.data,dim=1)


                

                ord=self.dic['FDDIMORDER']
                lab1=self.dic['FDF1LABEL']
                lab2=self.dic['FDF2LABEL']
                lab3=self.dic['FDF3LABEL']
                lab4=self.dic['FDF4LABEL']
                lab=lab1,lab2,lab3, lab4
                # Pseudo-axis identity is canonical topology, not a label guess.
                self.pseudo = topology.has_pseudo_axis
                # exit()
                self.labb=lab[int(ord[1])-1],lab[int(ord[0])-1]
                
                self.SetProjectedData(inv=True)
                self.dmax=numpy.max(self.data)  #get max intensity in spectrum
                Size=self.data.shape
                self.specsize=Size

                self.uc0max=self.uc0.ppm(0)
                self.uc0min=self.uc0.ppm(Size[0]-1)
                self.uc1max=self.uc1.ppm(0)
                self.uc1min=self.uc1.ppm(Size[1]-1)
            
                self.uc0 = ng.pipe.make_uc(self.dic, self.data, dim=0)
                x0,x1=self.uc0.ppm_limits()
                self.uc0.ppms_scale=numpy.linspace(x0, x1, int(self.uc0._size))
                self.uc1 = ng.pipe.make_uc(self.dic, self.data, dim=1)
                x0,x1=self.uc1.ppm_limits()
                self.uc1.ppms_scale=numpy.linspace(x0, x1, int(self.uc1._size))
                


                self.index0=[]#make index of carbon chemical shifts for index 0
                for i in range((Size[0])):
                    self.index0.append((self.uc0.ppm(0)-i*(-self.uc0.ppm(Size[0]-1)+self.uc0.ppm(0))/(Size[0]-1)))
                self.index1=[]#make index of carbon chemical shifts for index 1
                for i in range((Size[1])):
                    self.index1.append((self.uc1.ppm(0)-i*(-self.uc1.ppm(Size[1]-1)+self.uc1.ppm(0))/(Size[1]-1)))

                self.index0=numpy.array(self.index0)
                self.index1=numpy.array(self.index1)

                self.YY,self.XX=numpy.meshgrid(self.index1,self.index0)

        elif spectral_dim_count == 3: #for 3d data

            self.dic,self.data=ng.pipe.read(infile)
            #print(self.data)
            # uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)
            # uc1 = ng.pipe.make_uc(self.dic,self.data,dim=1)
            # uc2 = ng.pipe.make_uc(self.dic,self.data,dim=2)

            self.uc0 = ng.pipe.make_uc(self.dic, self.data, dim=0)
            x0,x1=self.uc0.ppm_limits()
            self.uc0.ppms_scale=numpy.linspace(x0, x1, int(self.uc0._size))
            self.uc1 = ng.pipe.make_uc(self.dic, self.data, dim=1)
            x0,x1=self.uc1.ppm_limits()
            self.uc1.ppms_scale=numpy.linspace(x0, x1, int(self.uc1._size))
            self.uc2 = ng.pipe.make_uc(self.dic, self.data, dim=2)
            x0,x1=self.uc2.ppm_limits()
            self.uc2.ppms_scale=numpy.linspace(x0, x1, int(self.uc2._size))

            ord=self.dic['FDDIMORDER']
            lab1=self.dic['FDF1LABEL']
            lab2=self.dic['FDF2LABEL']
            lab3=self.dic['FDF3LABEL']
            lab=lab1,lab2,lab3
            # Pseudo-axis identity is canonical topology, not a label guess.
            self.pseudo = topology.has_pseudo_axis

  
            self.labb=lab[int(ord[2])-1],lab[int(ord[1])-1],lab[int(ord[0])-1]
            self.dmax=numpy.max(numpy.fabs(self.data))  #get max intensity in spectrum
            Size=self.data.shape
            self.specsize=Size

            #self.projectedData = numpy.sum(numpy.sum(self.data, axis=0), axis=0)
            #read projection
            self.SetProjectedData()
            

            self.uc0max=self.uc0.ppm(0)
            self.uc0min=self.uc0.ppm(Size[0]-1)
            self.uc1max=self.uc1.ppm(0)
            self.uc1min=self.uc1.ppm(Size[1]-1)
            self.uc2max=self.uc2.ppm(0)
            self.uc2min=self.uc2.ppm(Size[2]-1)



            #print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
            #print "Labels: ",self.labb
            #print "dimension 0 limits (ppm): ", self.uc0min, self.uc0max  #carbon
            #print "dimension 1 limits (ppm): ", self.uc1min, self.uc1max  #direct
            #print "dimension 2 limits (ppm): ", self.uc2min, self.uc2max  #direct
            #print 'Maximum Intensity:',self.dmax
            self.index0=[]#make index of carbon chemical shifts for index 0
            for i in range((Size[0])):
                self.index0.append((self.uc0.ppm(0)-i*(-self.uc0.ppm(Size[0]-1)+self.uc0.ppm(0))/(Size[0]-1)))
            self.index1=[]#make index of carbon chemical shifts for index 1
            for i in range((Size[1])):
                self.index1.append((self.uc1.ppm(0)-i*(-self.uc1.ppm(Size[1]-1)+self.uc1.ppm(0))/(Size[1]-1)))
            self.index2=[]#make index of carbon chemical shifts for index 2
            for i in range((Size[2])):
                self.index2.append((self.uc2.ppm(0)-i*(-self.uc2.ppm(Size[2]-1)+self.uc2.ppm(0))/(Size[2]-1)))
            self.index0=numpy.array(self.index0)
            self.index1=numpy.array(self.index1)
            self.index2=numpy.array(self.index2)

            self.YY,self.XX,self.ZZ=numpy.meshgrid(self.index1,self.index0,self.index2)

            #print self.XX.shape
            #print self.YY.shape
            #print self.ZZ.shape
            #print self.data.shape
            """
            for i in range(len(self.index0)):
                for j in range(len(self.index1)):
                    for k in range(len(self.index2)):
                        #print self.index0[i],self.index1[j],self.index2[k],self.XX[i,j,k],self.YY[i,j,k],self.ZZ[i,j,k]
                        if(numpy.fabs(self.XX[i,j,k]-self.index0[i])>0.01):
                            print 'shit'
                        if(numpy.fabs(self.YY[i,j,k]-self.index1[j])>0.01):
                            print 'shit'
                        if(numpy.fabs(self.ZZ[i,j,k]-self.index2[k])>0.01):
                            print 'shit'
            """


        elif spectral_dim_count == 4: #for 4d data
            #read in spectrum
            self.dic,self.data = ng.pipe.read_lowmem(infile)
            uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)
            uc1 = ng.pipe.make_uc(self.dic,self.data,dim=1)
            uc2 = ng.pipe.make_uc(self.dic,self.data,dim=2)
            uc3 = ng.pipe.make_uc(self.dic,self.data,dim=3)
            Size=self.data.shape

            ord=self.dic['FDDIMORDER']
            lab1=self.dic['FDF1LABEL']
            lab2=self.dic['FDF2LABEL']
            lab3=self.dic['FDF3LABEL']
            lab4=self.dic['FDF4LABEL']
            self.SetDmax(infile) #set dmax either from file or from numpy
            #self.dmax=numpy.max(self.data)
            Size=self.data.shape
            self.specsize=Size


            lab=lab1,lab2,lab3,lab4
            self.labb=lab[int(ord[3])-1],lab[int(ord[2])-1],lab[int(ord[1])-1],lab[int(ord[0])-1]


            self.uc0max=uc0.ppm(0)
            self.uc0min=uc0.ppm(Size[0]-1)
            self.uc1max=uc1.ppm(0)
            self.uc1min=uc1.ppm(Size[1]-1)
            self.uc2max=uc2.ppm(0)
            self.uc2min=uc2.ppm(Size[2]-1)
            self.uc3max=uc3.ppm(0)
            self.uc3min=uc3.ppm(Size[3]-1)

            self.index0=[]#make index of carbon chemical shifts for index 0
            for i in range((Size[0])):
                self.index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
            self.index1=[]#make index of carbon chemical shifts for index 1
            for i in range((Size[1])):
                self.index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
            self.index2=[]#make index of carbon chemical shifts for index 2
            for i in range((Size[2])):
                self.index2.append((uc2.ppm(0)-i*(-uc2.ppm(Size[2]-1)+uc2.ppm(0))/(Size[2]-1)))
            self.index3=[]#make index of carbon chemical shifts for index 2
            for i in range((Size[3])):
                self.index3.append((uc3.ppm(0)-i*(-uc3.ppm(Size[3]-1)+uc3.ppm(0))/(Size[3]-1)))
            self.index0=numpy.array(self.index0)
            self.index1=numpy.array(self.index1)
            self.index2=numpy.array(self.index2)
            self.index3=numpy.array(self.index3)

            self.YY,self.XX=numpy.meshgrid(self.index1,self.index0)
            self.YY2,self.XX2=numpy.meshgrid(self.index3,self.index2)

            #print 'getting max...'
            #dmax=self.data.max() #get max data value

        # Register the authoritative main spectrum after all dimensionality-
        # specific setup is complete.  In particular, a true 2D spectrum is
        # the data PeakFrame displays directly; there is no separate 2D
        # projection to cache first.
        self.store.save_spectrum(
            'raw',
            dic=self.dic,
            data=self.data,
            spectrumfile=infile,
            dim=spectral_dim_count,
            labb=getattr(self, 'labb', None),
        )
        self._cache_spectrum_views("raw")
        # Physical 3p datasets also get one centrally owned logical view in
        # [pseudo, y, x] order for Pseudo3D and future ancillary consumers.
        self.get_pseudo3d_view("raw")
        print('Spectrum read complete.')

    def SetDmax(self,infile):
        if(os.path.exists(infile+'.dmax')==1):
            inny=open(infile+'.dmax')
            for line in inny.readlines():
                test=line.split()
                if(len(test)>1):
                    if(test[0]=='dmax:'):
                        self.dmax=float(test[1])
                        return
        self.dmax=numpy.max(numpy.fabs(self.data))  #get max intensity in spectrum
        outy=open(infile+'.dmax','w')
        outy.write('dmax: %e\n' % (self.dmax))
        outy.close()


    ##########################################################
    #Read peak list file.
    def ReadPeakListFile(self):  #master read-in function for reference peaklists.
        topology = self._active_topology()
        spectral_dim_count = topology.spectral_dim_count
        if spectral_dim_count < 2:
            errorMessage('Reference peak lists are always 2D and require at least 2D main data.')
            return
        self.PEAK=0    #set read in mode to off.
        self.peak=[]   #re-initialise peak list.

        indir=self.dirBox.GetValue()     #get working directory
        peakfile=self.peakBox.GetValue() #get peaklist file

        if spectral_dim_count == 2:  #if 2 dimensional.
            peakListLocation=self._resolve_spec_file(peakfile)
            self.peak=readpeaklist(peakListLocation)

            if topology.has_pseudo_axis: #for pseudo3D analysis.......
                for pk in self.peak:
                    self.alias(pk,pk.x,1)
                    self.alias(pk,pk.y,2)
                # Previously this unconditionally called parent.tabPseudo before
                # set_reference_peaks(). In Workflow the Fitting tab often does
                # not exist yet, so that exception aborts the reference load.
                pseudo = getattr(self.parent, 'tabPseudo', None)
                self._workflow_debug('ReadPeakListFile dim=2 pseudo: parsed %d peaks; tabPseudo=%r. Deferring SetPeaksToFit until after DataStore commit.' % (len(self.peak), type(pseudo).__name__ if pseudo is not None else None))

            else:
                for pk in self.peak:
                    self.alias(pk,pk.x,0)
                    self.alias(pk,pk.y,1)
            

            # Reference peaks remain a reference peak list; no connectivity mirror is built.

        if spectral_dim_count == 3:
            peakListLocation=self._resolve_spec_file(peakfile)
            print("Reading peak list:", peakListLocation)
            if(os.path.exists(peakListLocation)==0):
                errorMessage('Cannot find peak list.')
                return

            self.peak=readpeaklist(peakListLocation)
            for pk in self.peak:
                self.alias(pk,pk.y,0)
                self.alias(pk,pk.y,1)
                self.alias(pk,pk.x,2)
            self.PEAK=1

            self.pkIdx=[] #index of peak positions
            self.pkSlice1D=[] #1D slices
            for p in range(len(self.peak)):
                ptC=self.peak[p].indexJ #carbon
                ptH=self.peak[p].indexK #proton
                self.pkIdx.append((ptC,ptH))
                self.pkSlice1D.append(self.data[:,ptC,ptH])
        elif spectral_dim_count == 4:
            peakListLocation=self._resolve_spec_file(peakfile)

            print("Reading peak list:", peakListLocation)
            if(os.path.exists(peakListLocation)==0):
                errorMessage('Cannot find peak list.')

                return

            self.peak=readpeaklist(peakListLocation)
            for p in range(len(self.peak)):
                self.alias(self.peak[p],self.peak[p].x,0) #i  #H
                self.alias(self.peak[p],self.peak[p].y,1) #j  #C
                self.alias(self.peak[p],self.peak[p].x,2) #k  #H
                self.alias(self.peak[p],self.peak[p].y,3) #l  #C

            #for pk in self.peak:
            #    if(pk.name=='A100C-H'):
            #        print pk.name,pk.x,pk.y,pk.ppmI,pk.ppmJ,pk.ppmK,pk.ppmL
            #        print self.index2[0],self.index2[-1]

            self.PEAK=1


            self.EXTRACT=1
            if(self.EXTRACT==1):
                print('slicing 4D')
                self.thresh=float(self.threshBox.GetValue())
                specstr=''
                #specstr+=' '+str(len(self.specsize))
                specstr+=' '+str(self.coreBox.GetValue())
                specstr+=' '+self._resolve_spec_file(self.peakBox.GetValue())
                specstr+=' '+str(spectral_dim_count)
                specstr+=' '+self._resolve_input_path(self.infileBox.GetValue())
                args=self.deconBin+' '+specstr+' '+str(self.dmax*self.thresh)+' '+str(0)
                print(args)
                os.system(args)
            self.SetDmax(self._resolve_input_path(self.infileBox.GetValue())) #read dmax

        # Commit the newly-read reference list before constructing any slice
        # viewers.  SliceFrame/SliceFrame2D resolve their peak choices from the
        # shared DataStore during __init__; creating them first leaves the
        # reference payload temporarily empty and can make wx.ComboBox selection
        # fail natively on macOS.
        self.set_reference_peaks(self.peak, source_path=peakListLocation)
        self._workflow_debug('ReadPeakListFile committed reference list: count=%d source=%r' % (len(self.get_reference_peaks() or []), peakListLocation))
        if spectral_dim_count == 2 and topology.has_pseudo_axis:
            pseudo = getattr(self.parent, 'tabPseudo', None)
            setter = getattr(pseudo, 'SetPeaksToFit', None) if pseudo is not None else None
            if callable(setter):
                try:
                    setter()
                    self._workflow_debug('Refreshed existing pseudo panel after reference DataStore commit.')
                except Exception as exc:
                    self._workflow_debug('Existing pseudo panel refresh failed: %r' % exc)
            else:
                self._workflow_debug('No existing pseudo panel to refresh; this is expected before Inspect fitting results.')

        if spectral_dim_count == 3 and not topology.has_pseudo_axis:
            self.parent.AddTabThree(True,self)
            self.parent.AddTabFour(True,self)

        if spectral_dim_count == 4:
            # The historical 4D slice viewer is retired: it is tightly coupled
            # to MAGMA and legacy conn_data.  Preserve the reference peak list
            # and project state, but do not launch the legacy viewer.
            self._workflow_debug(
                '4D reference list loaded; legacy Slice4D/MAGMA viewer is retired.'
            )

        self.refresh_reference_peak_views()
        self.Status()



    ###############################################
    #Automatic setting of unidec
        
    def on_quick_button(self, event):
        self.convBox.SetValue('1e-5')
        self.maxiterBox.SetValue('25')

    def on_medium_button(self, event):
        self.convBox.SetValue('1e-7')
        self.maxiterBox.SetValue('50')

    def on_accurate_button(self, event):
        self.convBox.SetValue('1e-8')
        self.maxiterBox.SetValue('100')


# Historical public name retained for external launchers/importers.
deconFrame = NMRWorkspace
