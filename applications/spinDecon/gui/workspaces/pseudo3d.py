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
import wx
from spinDecon.domain.peaks import peakEntry
from spinDecon.gui.context import context_for, project_for, data_for
from spinDecon.domain.dimensions.viewer_contract import topology_for
import os
import numpy
import re
import subprocess
#import imp
import importlib
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar

class _ToolbarToggleState:
    """Small wx.CheckBox-compatible state holder for toolbar-owned toggles."""
    def __init__(self, value=False):
        self._value = bool(value)
        self._enabled = True
    def GetValue(self): return self._value
    def IsChecked(self): return self._value
    def SetValue(self, value): self._value = bool(value)
    def Enable(self, enabled=True): self._enabled = bool(enabled)
    def IsEnabled(self): return self._enabled

import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm
from matplotlib.widgets import RectangleSelector
import nmrglue as ng
import scipy.optimize as opt
import threading,copy
from spinDecon.processing.nmrpipe_scripts import MakeProj4D, MakeProj3D
from spinDecon.gui.plotting.array_utils import ensure_xy_points, scatter_xy_points
from spinDecon.project.parameter_store import update_parameter_file
from spinDecon.line_fitting.line_fitting import Unidec_line_fitting

from wx.lib.mixins.listctrl import ColumnSorterMixin

############################################################################
# Pseudo 3D frame including interfacing with line_fitting
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def sine_function(x, amp, J,c):
    return amp*numpy.sin(numpy.pi*J*x)

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()


class FudaProgressFrame(wx.Frame):
    """Live, non-blocking view of FUDA stdout and fitting progress."""
    def __init__(self, parent, units):
        wx.Frame.__init__(self, parent, title='FUDA fitting progress', size=(760, 430),
                          style=wx.DEFAULT_FRAME_STYLE | wx.RESIZE_BORDER)
        self.units = units
        self.completed = set()
        self.current = None
        self.stage = 'Starting FUDA'
        self.iteration = 0
        self._seen_outputs = set()

        panel = wx.Panel(self)
        root = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(panel, label='FUDA fitting')
        font = title.GetFont(); font.SetPointSize(font.GetPointSize()+3); font.SetWeight(wx.FONTWEIGHT_BOLD)
        title.SetFont(font)
        root.Add(title, 0, wx.LEFT|wx.RIGHT|wx.TOP, 12)
        self.summary = wx.StaticText(panel, label=self._summary_text())
        root.Add(self.summary, 0, wx.EXPAND|wx.ALL, 12)
        self.gauge = wx.Gauge(panel, range=max(1, len(units)), style=wx.GA_HORIZONTAL)
        root.Add(self.gauge, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 12)

        grid = wx.FlexGridSizer(2, 4, 6, 12); grid.AddGrowableCol(1, 1); grid.AddGrowableCol(3, 1)
        grid.Add(wx.StaticText(panel,label='Current:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.current_text=wx.StaticText(panel,label='Preparing input')
        grid.Add(self.current_text,1,wx.EXPAND)
        grid.Add(wx.StaticText(panel,label='Stage:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.stage_text=wx.StaticText(panel,label=self.stage)
        grid.Add(self.stage_text,1,wx.EXPAND)
        grid.Add(wx.StaticText(panel,label='Iteration:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.iter_text=wx.StaticText(panel,label='-')
        grid.Add(self.iter_text,1,wx.EXPAND)
        grid.Add(wx.StaticText(panel,label='Fit metric:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.metric_text=wx.StaticText(panel,label='-')
        grid.Add(self.metric_text,1,wx.EXPAND)
        root.Add(grid,0,wx.EXPAND|wx.ALL,12)

        root.Add(wx.StaticText(panel,label='FUDA output'),0,wx.LEFT|wx.RIGHT|wx.TOP,12)
        self.output=wx.TextCtrl(panel, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH2|wx.HSCROLL)
        # Keep the diagnostic shell useful but compact so the progress window
        # does not dominate the Pseudo3D workspace. It still expands when the
        # user manually enlarges the frame.
        self.output.SetMinSize((-1, 115))
        root.Add(self.output,1,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP|wx.BOTTOM,12)
        self.close_btn=wx.Button(panel,label='Close'); self.close_btn.Enable(False)
        self.close_btn.Bind(wx.EVT_BUTTON, lambda evt:self.Close())
        root.Add(self.close_btn,0,wx.ALIGN_RIGHT|wx.LEFT|wx.RIGHT|wx.BOTTOM,12)
        panel.SetSizer(root); self.SetMinSize((620,330))

    def _summary_text(self):
        peaks=sum(len(u['peaks']) for u in self.units)
        groups=sum(1 for u in self.units if u['group'] is not None)
        singles=sum(1 for u in self.units if u['group'] is None)
        return '%d fitting units: %d groups + %d individual peaks (%d peaks total)' % (len(self.units),groups,singles,peaks)

    def append_line(self, line):
        self.output.AppendText(line)
        if not line.endswith('\n'): self.output.AppendText('\n')
        self.output.ShowPosition(self.output.GetLastPosition())
        text=line.strip()
        m=re.search(r'Name of Group-of-peaks:\s*(\S+)', text)
        if m:
            self.current=m.group(1); self.stage='Initial fit'; self.iteration=0
            self.current_text.SetLabel(self._unit_label(self.current)); self.stage_text.SetLabel(self.stage); self.iter_text.SetLabel('-')
        if '*** Step One ***' in text: self.stage='Step 1 - intensities / linewidths'
        elif '*** Step Two ***' in text: self.stage='Step 2 - spectral parameters'
        elif '*** Final Step ***' in text: self.stage='Final optimisation'
        elif 'The fit converged successfully' in text: self.stage='Converged - writing results'
        self.stage_text.SetLabel(self.stage)
        mi=re.search(r'Iter:\s*sd\s*=\s*([^ ]+)\s+enorm\s*=\s*([^ ]+)', text)
        if mi:
            self.iteration += 1; self.iter_text.SetLabel(str(self.iteration)); self.metric_text.SetLabel('sd %s   enorm %s' % mi.groups())
        mo=re.search(r'Gnuplot ScriptFile has been written to the file .*[/\\]([^/\\]+)\.gnu', text)
        if mo: self._mark_peak_output(mo.group(1))
        ma=re.search(r'Peak:\s*(\S+)\s+has already been fitted', text)
        if ma: self._mark_peak_output(ma.group(1))

    def _unit_label(self, peak):
        for u in self.units:
            if peak in u['peaks']:
                return ('Group %s (%d peaks)' % (u['group'],len(u['peaks']))) if u['group'] is not None else peak
        return peak

    def _mark_peak_output(self, peak):
        self._seen_outputs.add(peak)
        for i,u in enumerate(self.units):
            if i not in self.completed and all(x in self._seen_outputs for x in u['peaks']): self.completed.add(i)
        self.gauge.SetValue(min(len(self.completed),len(self.units)))
        self.summary.SetLabel('%s   |   %d / %d complete' % (self._summary_text(),len(self.completed),len(self.units)))

    def finish(self, returncode):
        if returncode == 0:
            self.completed=set(range(len(self.units))); self.gauge.SetValue(len(self.units)); self.stage='Complete'
            self.summary.SetLabel('%s   |   Complete' % self._summary_text())
        else:
            self.stage='Failed (exit code %d)' % returncode
            self.summary.SetLabel('%s   |   Calculation failed' % self._summary_text())
        self.stage_text.SetLabel(self.stage); self.close_btn.Enable(True); self.Raise()

class Pseudo3D(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.topology = topology_for(tabOne)
        self.spectral_dim_count = self.topology.spectral_dim_count
        self.physical_dim_count = self.topology.physical_dim_count
        self.dim = self.spectral_dim_count  # compatibility alias: spectral only
        self.app_context = context_for(tabOne, parent)
        self.nmr_workspace = (self.app_context.nmr_workspace
                              if self.app_context is not None and self.app_context.nmr_workspace is not None
                              else tabOne)
        self.pseudo_service = self.app_context.pseudo if self.app_context is not None else None
        self.peak_service = self.app_context.peaks if self.app_context is not None else None
        if self.pseudo_service is None:
            from spinDecon.analysis.pseudo_service import PseudoAxisService
            self.pseudo_service = PseudoAxisService(tabOne)
        if self.peak_service is None:
            from spinDecon.analysis.peak_service import PeakService
            self.peak_service = PeakService(tabOne)
        self.state = project_for(tabOne, parent)
        self.store = data_for(tabOne, parent)
        self.sum=(0.,2.)
        # Contour minimum starts at the same absolute intensity threshold
        # used by the main NMR page.  Subsequent ENTER commits on the main
        # threshold control are propagated by sync_main_threshold().
        self.thresh = self.peak_service.threshold()
        self.offset=0
        self.peaks_drawn=False
        self.plotting_atom_results = False

        dmin, dmax = self.pseudo_service.unit_conversion_bounds


        
        self.create_main_panel()
        self.draw_figure()
        self.Fit()

    def _pseudo3d_view(self, key='raw'):
        """Return the canonical shared [pseudo, y, x] view for raw/decon data."""
        view = self.store.get_view(('pseudo3d', key)) if self.store is not None else None
        if view is None:
            view = self.pseudo_service.view(key)
        if view is None:
            raise RuntimeError('Pseudo3D requires its canonical view in the shared data store')
        return view

    def _is_single_plane_2d(self):
        """True when this fitting workspace is adapting a physical 2D spectrum."""
        view = self._pseudo3d_view()
        return bool(view.get('is_single_plane', False) and view.get('physical_dim') == 2)

    def _is_physical_2d_adapter(self):
        """Compatibility predicate for the shared Fitting workspace.

        A physical 2D spectrum reuses the Pseudo3D/Fitting UI as a single-plane
        adapter, but its restraint/fit peaks belong to the Full 2D list rather
        than the pseudo3D Reference 2D list.
        """
        return self._is_single_plane_2d()

    def _decon_pseudo3d_view(self):
        """Return the optional decon view without making its absence an error.

        Raw data are mandatory in this workspace; calculated data are not.
        Keeping those contracts separate lets a true 2D spectrum open before
        deconvolution and lets the toolbar become available as soon as a
        matching ``.decon`` spectrum has been published to the shared store.
        """
        view = self.store.get_view(('pseudo3d', 'decon')) if self.store is not None else None
        if view is None:
            view = self.pseudo_service.view('decon')
        return view

    def _reference_peaks(self):
        """Return the fitting/grouping peak objects for the active topology.

        Pseudo3D owns an independent Reference 2D list.  A physical 2D
        spectrum does not: its Full 2D list is the same scientific collection
        and is therefore the sole authority.  Full-list records are adapted to
        the small ``peakEntry`` interface expected by this legacy fitting UI.
        """
        if self.store is None:
            return []
        if not self._is_physical_2d_adapter():
            return self.store.peak_lists.get('reference', {}).get('peaks') or []

        records = self.store.peak_lists.get('full', {}).get('peaks') or []
        peaks = []
        for record in records:
            if not isinstance(record, dict):
                peaks.append(record)
                continue
            coords = tuple(record.get('coordinates') or ())
            if len(coords) < 2:
                fields = list(record.get('fields') or [])
                if len(fields) >= 3:
                    coords = (fields[1], fields[2])
            if len(coords) < 2:
                continue
            peak = peakEntry([record.get('name', ''), coords[0], coords[1]])
            peak.analysis = record.get('analysis', {})
            peaks.append(peak)
        return peaks

    def _reference_peak_overlay(self):
        """Return centrally projected reference peaks for the canonical view."""
        view = self._pseudo3d_view()
        list_key = 'full' if self._is_physical_2d_adapter() else 'reference'
        key = (list_key, str(view['x_label']), str(view['y_label']))
        payload = self.store.projected_peak_lists.get(key, {}) if self.store is not None else {}
        if not payload:
            self.pseudo_service.rebuild_projected_peaks()
            payload = self.store.projected_peak_lists.get(key, {}) if self.store is not None else {}
        return payload.get('peaks') or []

    def _reference_peak_by_name(self, name):
        for peak in self._reference_peaks():
            if peak.name == name:
                return peak
        return None

    def _groups(self):
        """Return authoritative overlap groups from the shared DataStore."""
        return self.store.Grps if self.store is not None else {}

    def _fuda_dir(self):
        """Return the controller-owned FUDA workspace below SpecPath."""
        return self.pseudo_service.fuda_dir()

    def _fuda_peak_file(self):
        return self.pseudo_service.fuda_peak_file()

    def _fuda_parameter_file(self):
        return self.pseudo_service.fuda_parameter_file()

    def _notify_analysis_changed(self):
        if self.pseudo_service is not None:
            self.pseudo_service.notify_changed()
            return
        self.pseudo_service.notify_changed()

    def _mark_pseudo_intensities_ready(self, **details):
        if self.store is None:
            return
        self.store.mark_pseudo_intensities_ready(**details)
        self._notify_analysis_changed()

    def _mark_pseudo_series_reviewed(self, **details):
        if self.store is None:
            return
        self.store.mark_pseudo_series_reviewed(**details)
        self._notify_analysis_changed()

    def _mark_pseudo_analysis_complete(self, **details):
        if self.store is None:
            return
        self.store.mark_pseudo_analysis_complete(**details)
        self._notify_analysis_changed()

    def _replace_groups(self, groups):
        return self.pseudo_service.replace_groups(groups)

    def _add_group(self, name, peaks=None):
        return self.pseudo_service.add_group(name, peaks)

    def _remove_group(self, name):
        return self.pseudo_service.remove_group(name)

    def _add_peak_to_group(self, name, peak_name):
        return self.pseudo_service.add_peak_to_group(name, peak_name)

    def _remove_peak_from_group(self, name, peak_name):
        return self.pseudo_service.remove_peak_from_group(name, peak_name)

    def _ensure_line_fitting(self, reset=False):
        """Create the fitter from canonical shared Pseudo3D/reference-peak data.

        The fitter receives a transient numeric adapter because the legacy fitting
        engine works in point coordinates.  Pseudo3D does not retain a second peak
        model; authoritative peaks remain in DataStore.peak_lists['reference'].
        """
        if not reset and getattr(self, 'line_fitting', None) is not None:
            return self.line_fitting

        view = self._pseudo3d_view()
        rows = []
        names = []
        for point in self._reference_peak_overlay():
            try:
                # Legacy fitter convention: [data-y index, data-x index, y ppm, x ppm].
                y_index = int(round(float(view['y_uc'].f(str(point['y']) + ' ppm'))))
                x_index = int(round(float(view['x_uc'].f(str(point['x']) + ' ppm'))))
            except (TypeError, ValueError):
                continue
            rows.append([y_index, x_index, float(point['y']), float(point['x'])])
            names.append(str(point.get('label', '')))

        peak_points = numpy.asarray(rows, dtype=float)
        if peak_points.size == 0:
            peak_points = numpy.empty((0, 4), dtype=float)

        self.line_fitting = Unidec_line_fitting(
            view['data'], peak_points, names,
            self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y,
            self.nu1, self.nu2, self.pseudo_service.threshold_fraction(),
            view['y_uc'], view['x_uc'])
        return self.line_fitting

    def saveBox(self,event):
        write = {
            '3p_min': self.textbox0.GetValue(),
            '3p_fac': self.textbox1.GetValue(),
            '3p_num': self.textbox2.GetValue(),
            '3p_radF1': self.radF1.GetValue(),
            '3p_radF2': self.radF2.GetValue(),
        }
        if self._groups():
            write['3p_grps'] = self.MakeGrpStr()
        self.pseudo_service.update_parameters(write)

    def set_default_values(self): # unpack group save
        if self.pseudo_service.parameter_value('3p_min', default='') not in ('', '0', 0):
            self.textbox0.SetValue(str(self.pseudo_service.parameter_float('3p_min')))
            self.textbox1.SetValue(str(self.pseudo_service.parameter_float('3p_fac')))
            self.textbox2.SetValue(str(self.pseudo_service.parameter_float('3p_num')))

        if self.pseudo_service.parameter_value('3p_radF1', default='') not in ('', '0', 0):
            self.radF1.SetValue(str(self.pseudo_service.parameter_float('3p_radF1')))
            self.radF2.SetValue(str(self.pseudo_service.parameter_float('3p_radF2')))

        self.RestoreGrpStr(str(self.pseudo_service.parameter_value('3p_grps', default='0')))
        self.SetPeaksToFit()

    def _peak_group_name(self, peak_name):
        """Return the group containing peak_name, or None for an orphan."""
        for grp, members in self._groups().items():
            if peak_name in members:
                return str(grp)
        return None

    def SetPeaksToFit(self):
        """Refresh fitting membership and the combined fitting-results list."""
        self.orph = {peak.name: 1 for peak in self._reference_peaks()}
        for grp, pks in self._groups().items():
            for pk in pks:
                if pk in self.orph:
                    self.orph[pk] = 0
        if hasattr(self, 'fittingList'):
            self.RefreshFittingList()

    def _fitting_columns(self):
        return ['Peak', 'grp', '%err', 'f01(ppm)', 'w1(Hz)', 'g1',
                'f02(ppm)', 'w2(Hz)', 'g2']

    def _read_overlap_group_from_out(self, out_file):
        """Return the Overlap_group integer recorded in a FUDA-style .out file.

        Both FUDA and spinUnidec write a two-line header where the value is on
        the line immediately following the ``Peak Name  Overlap_group`` label.
        ``None`` means that the header was absent or malformed.
        """
        try:
            with open(out_file) as handle:
                lines = iter(handle)
                for line in lines:
                    if 'Peak Name' not in line or 'Overlap_group' not in line:
                        continue
                    try:
                        value_line = next(lines)
                    except StopIteration:
                        return None
                    fields = value_line.split()
                    if fields and fields[0] == '#':
                        fields = fields[1:]
                    if len(fields) < 2:
                        return None
                    try:
                        return int(fields[1])
                    except ValueError:
                        return None
        except OSError:
            return None
        return None

    def _update_groups_from_fit_outputs(self):
        """Import spinUnidec overlap groups from the current fit directory.

        FUDA writes ``Overlap_group = -1`` and therefore must not alter the
        overlap groups restored from the project's SAVE/LOAD state.  A
        non-negative value identifies spinUnidec output.  If at least one such
        value is found, the fit output becomes authoritative: all non-negative
        values are grouped by their recorded number and replace ``store.Grps``.
        Peaks whose output records -1 remain ungrouped.
        """
        fit_dir = self._fuda_dir()
        reference_names = [peak.name for peak in self._reference_peaks()]
        imported = {}
        saw_spinunidec_group = False

        for peak_name in reference_names:
            out_file = os.path.join(fit_dir, peak_name + '.out')
            if not os.path.isfile(out_file):
                continue
            group_no = self._read_overlap_group_from_out(out_file)
            if group_no is None or group_no == -1:
                continue
            if group_no < -1:
                continue
            saw_spinunidec_group = True
            imported.setdefault(str(group_no), []).append(peak_name)

        if not saw_spinunidec_group:
            # FUDA (all -1), no results, or malformed/legacy output: retain the
            # groups loaded from the normal system SAVE/LOAD file unchanged.
            return False

        self._replace_groups(imported)
        return True

    def _parse_fitting_result(self, peak_name):
        """Return FUDA result columns for pseudo3D *and* native 2D output.

        nmrPipeFit prefixes spectral parameter rows with ``#`` for arrayed
        pseudo3D fits, but deliberately writes the same rows without ``#`` in
        its ZCOOR=2D path.  Parse by parameter name rather than by comment
        prefix/order so both formats populate the same GUI columns.
        """
        values = [''] * 7
        out_file = os.path.join(self._fuda_dir(), peak_name + '.out')
        if os.path.isfile(out_file):
            params = {}
            in_results = False
            with open(out_file) as handle:
                for line in handle:
                    test = line.split()
                    if not test:
                        continue
                    if 'Results' in test and 'fit' in test:
                        in_results = True
                        continue
                    if in_results and (test[0].startswith('##') or
                                       ('Z-coordinate' in line)):
                        break
                    if not in_results:
                        continue
                    # Arrayed output: '# f01(ppm) value esd'
                    # Native 2D:      '  f01(ppm) value esd'
                    fields = test[1:] if test[0] == '#' else test
                    if len(fields) < 3 or fields[0] == 'Parameter':
                        continue
                    name = fields[0].split('(')[0]
                    if name in ('f01','w1','g1','f02','w2','g2'):
                        try:
                            params[name] = '%.3f' % float(fields[1])
                        except ValueError:
                            params[name] = fields[1]
            for col, name in enumerate(('f01','w1','g1','f02','w2','g2'), start=1):
                values[col] = params.get(name, '')

        dat_file = os.path.join(self._fuda_dir(), peak_name + '.dat')
        if os.path.isfile(dat_file):
            yc, yd = [], []
            with open(dat_file) as handle:
                for line in handle:
                    test = line.split()
                    if len(test) == 4 and test[0] != '#':
                        try:
                            yc.append(float(test[2])); yd.append(float(test[3]))
                        except ValueError:
                            pass
            if yc:
                yc, yd = numpy.asarray(yc), numpy.asarray(yd)
                denom = numpy.max(numpy.abs(yc))
                if denom:
                    values[0] = '%.3f' % (numpy.sqrt(numpy.average((yc-yd)**2.)) / denom * 100.)
        return values

    def _ordered_fitting_peak_names(self):
        """Return reference peaks with every overlap group kept contiguous.

        Group order follows the system group dictionary and member order follows
        each group's stored member list.  Remaining (ungrouped) reference peaks
        follow in reference-list order.
        """
        reference_names = [peak.name for peak in self._reference_peaks()]
        reference_set = set(reference_names)
        ordered, seen = [], set()
        for grp, members in self._groups().items():
            for name in members:
                if name in reference_set and name not in seen:
                    ordered.append(name)
                    seen.add(name)
        ordered.extend(name for name in reference_names if name not in seen)
        return ordered

    def fitting_summary_rows(self):
        """Return the exact data shown in the Fitting results list.

        This is the shared report/UI provider: report generation must not parse
        the list control, so the GUI and PDF cannot drift apart.
        """
        rows = []
        for name in self._ordered_fitting_peak_names():
            grp = self._peak_group_name(name)
            rows.append([name, grp if grp is not None else '-'] +
                        list(self._parse_fitting_result(name)))
        return self._fitting_columns(), rows

    def fitting_report_units(self):
        """Return overlap groups plus ungrouped peaks as report fitting units."""
        reference = set(self._ordered_fitting_peak_names())
        units, grouped = [], set()
        for grp, members in self._groups().items():
            members = [name for name in members if name in reference]
            if members:
                units.append({'group': str(grp), 'peaks': members})
                grouped.update(members)
        for name in self._ordered_fitting_peak_names():
            if name not in grouped:
                units.append({'group': None, 'peaks': [name]})
        return units

    def export_fitting_report_figures(self, report_dir, units=None):
        """Export the live two-pane Fitting view once per fitting unit.

        The left pane remains the GUI's 2D projection and ``ReadFuda`` renders
        the same rotated 3D data/fit view used by the Fitting window.
        """
        paths = []
        report_dir = os.fspath(report_dir)
        old_peak = getattr(self, 'selected_fitting_peak', None)
        old_slice = self.PeakComboSlice.GetSelection() if hasattr(self, 'PeakComboSlice') else -1
        old_accum = self.cb_accum.GetValue() if hasattr(self, 'cb_accum') else False
        try:
            if hasattr(self, 'cb_accum'):
                self.cb_accum.SetValue(False)
            if hasattr(self, 'PeakComboSlice') and self.PeakComboSlice.GetCount():
                self.PeakComboSlice.SetSelection(0)
            for number, unit in enumerate(units if units is not None else self.fitting_report_units(), 1):
                representative = unit['peaks'][0]
                dat = os.path.join(self._fuda_dir(), representative + '.dat')
                if not os.path.isfile(dat):
                    continue
                self.selected_fitting_peak = representative
                self.SetCurrentPeaks()
                self.ReadFuda()
                # Mirror an actual selection in the Fitting window: ReadFuda
                # redraws the canvas, so the group markers must be applied
                # afterwards before the report snapshot is taken.
                self.show_fitting_peak_markers(self.fudapeaks)
                filename = 'pseudo3d_fit_%03d.pdf' % number
                self.canvas.print_figure(os.path.join(report_dir, filename))

                # The main two-pane figure deliberately shows slice zero for a
                # stable overview.  Also export the 3D data/fit comparison for
                # every pseudo-axis slice so the report exposes the full fit.
                slice_files = []
                elev = getattr(self.axesFuda, 'elev', 30)
                azim = getattr(self.axesFuda, 'azim', -60)
                # Use one common Z scale for every slice in this fitting
                # unit.  This makes changes in amplitude across the pseudo
                # axis visually meaningful instead of autoscaling each panel.
                z_limits = self._fitting_slice_z_limits()
                for slice_index in range(len(getattr(self, 'dat', []))):
                    slice_fig = Figure(figsize=(3.0, 3.0))
                    slice_ax = slice_fig.add_subplot(111, projection='3d')
                    self._draw_fitting_slice(slice_ax, slice_index)
                    if z_limits is not None:
                        slice_ax.set_zlim(*z_limits)
                    slice_ax.view_init(elev=elev, azim=azim)
                    slice_ax.set_title(self._fitting_slice_label(slice_index), fontsize=8)
                    slice_fig.tight_layout(pad=0.5)
                    slice_name = 'pseudo3d_fit_%03d_slice_%03d.pdf' % (number, slice_index + 1)
                    slice_fig.savefig(os.path.join(report_dir, slice_name), bbox_inches='tight')
                    plt.close(slice_fig)
                    slice_files.append(slice_name)
                unit['slice_figures'] = slice_files
                paths.append((filename, unit))
        finally:
            self.selected_fitting_peak = old_peak
            if hasattr(self, 'PeakComboSlice') and old_slice >= 0:
                self.PeakComboSlice.SetSelection(old_slice)
            if hasattr(self, 'cb_accum'):
                self.cb_accum.SetValue(old_accum)
            if old_peak:
                try:
                    self.SetCurrentPeaks(); self.ReadFuda()
                except Exception:
                    pass
        return paths

    def RefreshFittingList(self, preserve_selection=True):
        """Scan SpecPath/fit and rebuild the complete fitting peak list.

        Every reference peak is listed even when no FUDA output exists; in that
        case its result columns remain blank.  Group members are deliberately
        adjacent so the table mirrors the current overlap configuration.
        """
        if not hasattr(self, 'fittingList'):
            return
        selected = getattr(self, 'selected_fitting_peak', None) if preserve_selection else None
        self.fittingList.Freeze()
        try:
            self.fittingList.DeleteAllItems()
            _columns, fitting_rows = self.fitting_summary_rows()
            for values in fitting_rows:
                name = values[0]
                row = self.fittingList.InsertItem(self.fittingList.GetItemCount(), name)
                for col, value in enumerate(values[1:], start=1):
                    self.fittingList.SetItem(row, col, str(value))
                if name == selected:
                    self.fittingList.Select(row)
                    self.fittingList.Focus(row)
            if selected is None and self.fittingList.GetItemCount():
                self.selected_fitting_peak = self.fittingList.GetItem(0, 0).GetText()
        finally:
            self.fittingList.Thaw()

    def show_fitting_window(self, silent=False):
        """Refresh the SpinUniDec fitting palette and optionally show it.

        Report generation uses ``silent=True`` so the exact same fitting-window
        refresh path is executed without flashing a modeless window on screen.
        """
        self._update_groups_from_fit_outputs()
        self.SetPeaksToFit()
        if not silent:
            self._show_tool_window(self.fittingFrame)

    def fitting_window_report_data(self):
        """Return fitting rows and units from the refreshed Fitting window.

        SpinUniDec output is imported by ``show_fitting_window`` before the list
        is rebuilt.  The report deliberately reads the resulting list control,
        rather than using project-memory groups as an independent source.
        """
        self.show_fitting_window(silent=True)
        columns = self._fitting_columns()
        rows = [[self.fittingList.GetItem(row, col).GetText()
                 for col in range(len(columns))]
                for row in range(self.fittingList.GetItemCount())]
        units, by_group, order = [], {}, []
        for row in rows:
            peak, group = row[0], row[1]
            if group and group != '-':
                if group not in by_group:
                    by_group[group] = []
                    order.append(group)
                by_group[group].append(peak)
            else:
                units.append({'group': None, 'peaks': [peak]})
        grouped = [{'group': group, 'peaks': by_group[group]} for group in order]
        # Preserve fitting-list order: grouped blocks first because the list is
        # itself ordered that way, then orphan peaks in their displayed order.
        return columns, rows, grouped + units

    def _on_fitting_list_selected(self, event):
        self.selected_fitting_peak = self.fittingList.GetItem(event.GetIndex(), 0).GetText()
        self.on_fitting_peak_selected(event)

    def invalidate_fits_for_peaks(self, peak_names):
        """Remove FUDA outputs made invalid by an overlap-group edit."""
        fuda_dir = self._fuda_dir()
        for peak_name in set(peak_names):
            for ext in ('.dat', '.out'):
                path = os.path.join(fuda_dir, peak_name + ext)
                try:
                    if os.path.isfile(path):
                        os.remove(path)
                except OSError:
                    pass
        if not self.cb_accum.GetValue():
            self.axesFuda.clear()
            self.canvas.draw_idle()
        self.SetPeaksToFit()

    #the string: groups delimited by pipe |, and grp and entires delimited by :, and entired deliminted by ,
    def MakeGrpStr(self):
        grpStr=''
        for j,(grp,vals) in enumerate(self._groups().items()):
            if(len(vals)==0): #don't save if empty.
                continue
            
            if(j!=0):
                grpStr+='|'
            grpStr+=grp+':'
            
            for i,val in enumerate(vals):
                if(i!=0):
                    grpStr+=','
                grpStr+=val

        return grpStr
    
    def RestoreGrpStr(self,grpStr):
        if(grpStr=='0'):
            return

        grps=grpStr.split('|')

        for grpLine in grps:
            grp=grpLine.split(':')[0]
            vals=grpLine.split(':')[1]

            self._add_group(grp, vals.split(','))
        
        
    def onFocus(self, event):
        event.Skip()

    def _make_modeless_window(self, title):
        """Create a modeless tool window whose controls remain owned by Pseudo3D."""
        frame = wx.Frame(self.GetTopLevelParent(), title=title,
                         style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(frame)
        frame.Bind(wx.EVT_CLOSE, lambda evt, f=frame: (f.Hide(), evt.Veto()))
        return frame, panel

    def _show_tool_window(self, frame):
        if not frame.IsShown():
            frame.Show()
        frame.Raise()

    def _close_tool_window(self, frame):
        frame.Hide()

    def _on_contour_enter(self, event):
        """Persist contour settings and redraw after Enter in any contour field."""
        self.saveBox(None)
        self.draw_figureGO()

    def sync_main_threshold(self, threshold, redraw=True):
        """Apply the main NMR absolute threshold as this view's contour minimum.

        This is intentionally dimension-independent: physical pseudo3D and
        the one-plane 2D adapter should respond identically when Threshold is
        committed on the NMR page.
        """
        try:
            value = float(threshold)
        except (TypeError, ValueError):
            return
        if not numpy.isfinite(value) or value <= 0:
            return
        self.thresh = value
        if hasattr(self, 'textbox0'):
            self.textbox0.SetValue(str(value))
        if redraw and hasattr(self, 'axes'):
            self.draw_figureGO()

    def contour_box(self):
        self.contourFrame, panel = self._make_modeless_window('Contours')
        row = wx.BoxSizer(wx.HORIZONTAL)
        self.text1 = wx.StaticText(panel, -1, 'Min:')
        self.text2 = wx.StaticText(panel, -1, 'Factor:')
        self.text3 = wx.StaticText(panel, -1, 'Number:')
        self.textbox0 = wx.TextCtrl(panel, size=(100,22), style=wx.TE_PROCESS_ENTER)
        self.textbox1 = wx.TextCtrl(panel, size=(50,22), style=wx.TE_PROCESS_ENTER)
        self.textbox2 = wx.TextCtrl(panel, size=(50,22), style=wx.TE_PROCESS_ENTER)
        self.textbox0.SetValue(str(self.thresh))
        self.textbox1.SetValue(str(1.2))
        self.textbox2.SetValue(str(15))
        for ctrl in (self.textbox0, self.textbox1, self.textbox2):
            ctrl.Bind(wx.EVT_TEXT_ENTER, self._on_contour_enter)
        close = wx.Button(panel, -1, 'Close', size=(-1,22))
        close.Bind(wx.EVT_BUTTON, lambda evt: self._close_tool_window(self.contourFrame))
        for widget in (self.text1, self.textbox0, self.text2, self.textbox1,
                       self.text3, self.textbox2, close):
            row.Add(widget, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        panel.SetSizer(row)
        row.Fit(panel)
        self.contourFrame.Fit()

    def drawing_box(self):
        # These controls are created on the panel, then reparented into the
        # native Matplotlib toolbar in create_main_panel().
        self.cb_grid = wx.CheckBox(self, -1, 'Peaks', style=wx.ALIGN_RIGHT)
        self.cb_grid.Hide()  # State-only; visible toggle is the Matplotlib Peaks tool.
        self.cb_calc = _ToolbarToggleState(False)
        self.up_button2d = wx.Button(self, -1, '+', size=(20,22))
        self.down_button2d = wx.Button(self, -1, '-', size=(20,22))
        self.slice_box = wx.TextCtrl(self, size=(35,22), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.Bind(wx.EVT_BUTTON, self.on_up_button, self.up_button2d)
        self.Bind(wx.EVT_BUTTON, self.on_down_button, self.down_button2d)

    def on_up_button(self, event):
        if self.current_slice < self._pseudo3d_view()['data'].shape[0]-1:
            self.current_slice+=1
        self.thresh=float(numpy.max(self._pseudo3d_view()['data'][0])*self.pseudo_service.threshold_fraction())

        self.draw_figureGO()

    def on_down_button(self, event):
        if self.current_slice > 0:
            self.current_slice-=1

        self.thresh=float(numpy.max(self._pseudo3d_view()['data'][0])*self.pseudo_service.threshold_fraction())


        self.draw_figureGO()

    """
def convex_hull_graham(points):
    '''
    Returns points on convex hull in CCW order according to Graham's scan algorithm. 
    By Tom Switzer <thomas.switzer@gmail.com>.
    '''
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    def cmp(a, b):
        return (a > b) - (a < b)

    def turn(p, q, r):
        return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)

    def _keep_left(hull, r):
        while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
            hull.pop()
        if not len(hull) or hull[-1] != r:
            hull.append(r)
        return hull

    points = sorted(points)
    l = reduce(_keep_left, points, [])
    u = reduce(_keep_left, reversed(points), [])
    return l.extend(u[i] for i in range(1, len(u) - 1)) or l
    """


    def anal_box(self):
        self.analysisFrame, panel = self._make_modeless_window('Analysis')
        row = wx.BoxSizer(wx.HORIZONTAL)
        self.AnalCombo = wx.ComboBox(panel, -1, size=(80,22), style=wx.CB_READONLY)
        # These are the downstream analyses currently implemented for a
        # physical 2D+pseudo intensity series.  Keep the choice here (rather
        # than duplicating it in Workflow) so the specialist panel remains the
        # authoritative analysis selector.
        for label in self.available_downstream_analyses():
            self.AnalCombo.Append(label)
        saved = self.selected_downstream_analysis()
        if saved and self.AnalCombo.FindString(saved) != wx.NOT_FOUND:
            self.AnalCombo.SetStringSelection(saved)
        elif self.AnalCombo.GetCount():
            self.AnalCombo.SetSelection(0)
        self.Analysebutton = wx.Button(panel, -1, 'Confirm / open', size=(-1,22))
        self.Analysebutton.Bind(wx.EVT_BUTTON, self.OnAnalyseButton)
        close = wx.Button(panel, -1, 'Close', size=(-1,22))
        close.Bind(wx.EVT_BUTTON, lambda evt: self._close_tool_window(self.analysisFrame))
        for widget in (self.AnalCombo, self.Analysebutton, close):
            row.Add(widget, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        panel.SetSizer(row)
        row.Fit(panel)
        self.analysisFrame.Fit()

    def available_downstream_analyses(self):
        """Return analysis types implemented for this pseudo3D datatype."""
        if self._is_single_plane_2d():
            return []
        return ['CPMG', 'Decay']

    def selected_downstream_analysis(self):
        if self.pseudo_service is not None:
            return self.pseudo_service.downstream_analysis
        return self.pseudo_service.downstream_analysis

    def save_downstream_analysis(self, selection):
        """Persist the confirmed analysis type in the normal system file."""
        selection = str(selection or '').strip()
        if selection not in self.available_downstream_analyses():
            return False
        self.pseudo_service.set_downstream_analysis(selection)
        return True

    def show_analysis_selector(self):
        saved = self.selected_downstream_analysis()
        if saved and self.AnalCombo.FindString(saved) != wx.NOT_FOUND:
            self.AnalCombo.SetStringSelection(saved)
        self._show_tool_window(self.analysisFrame)

    def open_saved_analysis(self):
        saved = self.selected_downstream_analysis()
        if not saved or saved not in self.available_downstream_analyses():
            self.show_analysis_selector()
            return False
        self.AnalCombo.SetStringSelection(saved)
        return self.OnAnalyseButton(None, workflow_entry=True)

    def OnAnalyseButton(self,event, workflow_entry=False):
        sele=self.AnalCombo.GetValue()
        if not self.save_downstream_analysis(sele):
            wx.MessageBox('Choose an analysis type first.', 'Analysis', wx.OK | wx.ICON_WARNING)
            return False

        self.flash_status_message('Opening %s analysis' % sele)

        if(sele=='CPMG'):
            from spinDecon.gui.workspaces import cpmg as CPMGframe
            cpmgResults=importlib.reload(CPMGframe)
            bool=cpmgResults.CPMGMan(self, auto_prepare=workflow_entry)
        elif(sele=='Decay'):
            from spinDecon.gui.workspaces import decay as DecayFrame
            DecayResults=importlib.reload(DecayFrame)
            bool=DecayResults.DecayMan(self, auto_prepare=workflow_entry)
        return True
        
    def fitting_box(self):
        """Build the combined fitting/results window."""
        self.fittingFrame, panel = self._make_modeless_window('Fitting')
        self.fittingStatusBar = self.fittingFrame.CreateStatusBar(1)
        self.fittingStatusBar.SetStatusText('Ready')
        self.selected_fitting_peak = None

        # Left: the former Overview information is now the primary peak selector.
        self.fittingList = wx.ListCtrl(panel, -1, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for col, label in enumerate(self._fitting_columns()):
            self.fittingList.InsertColumn(col, label)
        self.fittingList.SetColumnWidth(0, 110)
        self.fittingList.SetColumnWidth(1, 55)
        for col in range(2, len(self._fitting_columns())):
            self.fittingList.SetColumnWidth(col, 85)
        self.fittingList.Bind(wx.EVT_LIST_ITEM_SELECTED, self._on_fitting_list_selected)

        self.peak_fit_button = wx.Button(panel, -1, 'Peak Fit')
        self.peak_draw_button = wx.Button(panel, -1, 'Draw!')
        self.cb_accum = wx.CheckBox(panel, -1, 'Accumulate')
        self.cb_err = wx.CheckBox(panel, -1, 'Difference')
        self.Upbutton = wx.Button(panel, -1, '+', size=(32,-1))
        self.Downbutton = wx.Button(panel, -1, '-', size=(32,-1))
        self.PeakComboSlice = wx.ComboBox(panel, -1, style=wx.CB_READONLY)
        for slice_number in range(int(self._pseudo3d_view()['data'].shape[0])):
            self.PeakComboSlice.Append(str(slice_number + 1))
        if self.PeakComboSlice.GetCount(): self.PeakComboSlice.SetSelection(0)
        self.radF1 = wx.TextCtrl(panel, value='0.1', size=(70,-1), style=wx.TE_PROCESS_ENTER)
        self.radF2 = wx.TextCtrl(panel, value='0.4', size=(70,-1), style=wx.TE_PROCESS_ENTER)
        self.Groupbutton = wx.Button(panel, -1, 'Groups')
        self.Savebutton = wx.Button(panel, -1, 'Save')
        self.SaveResultsbutton = wx.Button(panel, -1, 'SaveResults')
        self.fitallbutton = wx.Button(panel, -1, 'Fit all')
        self.Cleanbutton = wx.Button(panel, -1, 'Clean fuda')
        close = wx.Button(panel, -1, 'Close')

        self.peak_fit_button.Bind(wx.EVT_BUTTON, self.peak_fit)
        self.peak_draw_button.Bind(wx.EVT_BUTTON, self.peak_draw)
        self.Groupbutton.Bind(wx.EVT_BUTTON, self.groupBox)
        self.Savebutton.Bind(wx.EVT_BUTTON, self.saveBox)
        self.SaveResultsbutton.Bind(wx.EVT_BUTTON, self.on_save_results)
        self.fitallbutton.Bind(wx.EVT_BUTTON, self.on_fitall)
        self.Cleanbutton.Bind(wx.EVT_BUTTON, self.on_clean)
        self.Upbutton.Bind(wx.EVT_BUTTON, self.OnUpbutton)
        self.Downbutton.Bind(wx.EVT_BUTTON, self.OnDownbutton)
        self.PeakComboSlice.Bind(wx.EVT_COMBOBOX, self.peak_draw)
        close.Bind(wx.EVT_BUTTON, lambda evt: self._close_tool_window(self.fittingFrame))

        right = wx.BoxSizer(wx.VERTICAL)
        slice_row = wx.BoxSizer(wx.HORIZONTAL)
        slice_row.Add(wx.StaticText(panel, -1, 'Slice:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        slice_row.Add(self.PeakComboSlice, 1, wx.RIGHT, 5)
        slice_row.Add(self.Upbutton, 0, wx.RIGHT, 3); slice_row.Add(self.Downbutton, 0)
        right.Add(slice_row, 0, wx.EXPAND | wx.BOTTOM, 8)
        grid = wx.FlexGridSizer(2, 2, 6, 6); grid.AddGrowableCol(1, 1)
        grid.Add(wx.StaticText(panel, -1, 'F1 radius:'), 0, wx.ALIGN_CENTER_VERTICAL); grid.Add(self.radF1, 1, wx.EXPAND)
        grid.Add(wx.StaticText(panel, -1, 'F2 radius:'), 0, wx.ALIGN_CENTER_VERTICAL); grid.Add(self.radF2, 1, wx.EXPAND)
        right.Add(grid, 0, wx.EXPAND | wx.BOTTOM, 8)
        right.Add(self.cb_accum, 0, wx.BOTTOM, 4); right.Add(self.cb_err, 0, wx.BOTTOM, 10)
        for widget in (self.peak_fit_button, self.peak_draw_button, self.Groupbutton,
                       self.Savebutton, self.SaveResultsbutton, self.fitallbutton, self.Cleanbutton):
            right.Add(widget, 0, wx.EXPAND | wx.BOTTOM, 5)
        right.AddStretchSpacer(1); right.Add(close, 0, wx.EXPAND)

        main = wx.BoxSizer(wx.HORIZONTAL)
        main.Add(self.fittingList, 3, wx.EXPAND | wx.ALL, 8)
        main.Add(right, 1, wx.EXPAND | wx.TOP | wx.RIGHT | wx.BOTTOM, 8)
        panel.SetSizer(main)
        self.fittingFrame.SetMinSize((820, 360))
        self.fittingFrame.SetSize((1100, 520))
        panel.Layout()
        self.RefreshFittingList(False)

    def _set_fitting_status(self, message):
        """Show concise fitting feedback without writing legacy debug output."""
        bar = getattr(self, 'fittingStatusBar', None)
        if bar is not None:
            bar.SetStatusText(str(message))

    def on_save_results(self, event=None):
        """Save the currently displayed fitting-results table as SpecPath/out.tab."""
        # Refresh first so the exported file is an exact snapshot of the table
        # the user sees, including current overlap groups and calculated %err.
        self.RefreshFittingList()
        output_file = os.path.join(self.pseudo_service.output_dir(), 'out.tab')
        columns = self._fitting_columns()
        try:
            with open(output_file, 'w') as handle:
                # Metadata header.  Keep every header line comment-prefixed so
                # the numeric/results body remains straightforward to parse.
                handle.write('# Fitting results\n')
                handle.write('# F1 radius: %s ppm\n' % self.radF1.GetValue().strip())
                handle.write('# F2 radius: %s ppm\n' % self.radF2.GetValue().strip())
                handle.write('# Overlap groups:\n')
                for grp in sorted(self._groups(), key=lambda value: (
                        0, int(value)) if str(value).lstrip('-').isdigit() else (1, str(value))):
                    members = self._groups().get(grp, [])
                    # Only report genuine overlap groups in the header.
                    # Single-member groups are orphans and are intentionally omitted.
                    if len(members) > 1:
                        handle.write('# grp%s %s\n' % (grp, ', '.join(str(member) for member in members)))
                handle.write('# ' + '\t'.join(columns) + '\n')
                for row in range(self.fittingList.GetItemCount()):
                    values = [self.fittingList.GetItem(row, col).GetText()
                              for col in range(len(columns))]
                    handle.write('\t'.join(values) + '\n')
        except OSError as exc:
            self._show_fitting_message('Could not save fitting results to %s\n\n%s' %
                                       (output_file, exc), 'SaveResults')
            return
        self._set_fitting_status('Saved fitting results: %s' % os.path.basename(output_file))
        self._show_fitting_message('Fitting results saved to:\n%s' % output_file,
                                   'SaveResults')

    def _show_fitting_message(self, message, title='Fitting'):
        wx.MessageBox(message, title, wx.OK | wx.ICON_INFORMATION,
                      parent=self.GetTopLevelParent())

    def _has_fuda_results(self):
        """Return True when the current SpecPath/fit contains fit output."""
        fuda_dir = self._fuda_dir()
        if not os.path.isdir(fuda_dir):
            return False
        try:
            return any(name.endswith(('.out', '.dat'))
                       for name in os.listdir(fuda_dir))
        except OSError:
            return False

    def overviewBox(self, event):
        """Open a fresh overview of the current FUDA calculation."""
        if not self._has_fuda_results():
            self._show_fitting_message(
                'No FUDA fitting results were found in:\n%s\n\nRun a peak fit first, then open Overview.'
                % self._fuda_dir(), 'Overview')
            return
        old = getattr(self, '_overview_summary_frame', None)
        if old is not None:
            try:
                old.Destroy()
            except RuntimeError:
                pass
        try:
            self._overview_summary_frame = overManFrame(
                self.GetTopLevelParent(), 10, 'Overview', self)
            self._overview_summary_frame.Show(True)
            self._overview_summary_frame.Raise()
        except Exception as exc:
            self._show_fitting_message(
                'The FUDA results were found, but the Overview window could not be opened.\n\n%s'
                % exc, 'Overview')

    def groupBox(self, event):
        """Open a fresh group summary/editor using the current store state."""
        if not self._reference_peaks():
            self._show_fitting_message(
                'No 2D peaks are currently available.\n\nLoad or create the Full 2D peak list before opening Groups.',
                'Groups')
            return
        old = getattr(self, '_group_summary_frame', None)
        if old is not None:
            try:
                old.Destroy()
            except RuntimeError:
                pass
        try:
            self._group_summary_frame = groupManFrame(
                self.GetTopLevelParent(), 10, 'Groups', self)
            self._group_summary_frame.Show(True)
            self._group_summary_frame.Raise()
        except Exception as exc:
            self._show_fitting_message(
                'The Groups window could not be opened.\n\n%s' % exc, 'Groups')


    def _selected_fuda_result_exists(self):
        """Return True if the current peak/group has readable FUDA output."""
        self.SetCurrentPeaks()
        if not self.fudapeaks:
            return False
        return os.path.isfile(os.path.join(self._fuda_dir(), self.fudapeaks[0] + '.dat'))

    def on_fitting_peak_selected(self, event):
        """Display the selected fit and mark its complete group on both spectra."""
        # SetCurrentPeaks resolves a table row to all members of its overlap
        # group.  Keep the existing marker on the main Projection spectrum.
        self.SetCurrentPeaks()
        if self._is_single_plane_2d():
            main_view = getattr(self.parent, 'tabTwo', None)
            marker_fn = getattr(main_view, 'show_fitting_peak_markers', None)
            if callable(marker_fn):
                marker_fn(self.fudapeaks)

        # ReadFuda()/the empty-result path can redraw this canvas, so add the
        # left-hand spectrum markers afterwards.  This keeps them visible even
        # when selecting a peak that already has fitted results.
        if self._selected_fuda_result_exists():
            self.ReadFuda()
        else:
            if not self.cb_accum.GetValue():
                self.axesFuda.clear()
                self.canvas.draw()
        self.show_fitting_peak_markers(self.fudapeaks)

    def show_fitting_peak_markers(self, peak_names):
        """Mark the selected peak/overlap group on the left-hand 2D spectrum.

        The fitting table resolves a selected row to ``peak_names`` before this
        method is called, so an overlap group is shown in its entirety.  The
        authoritative reference-peak coordinates are used (X=direct,
        Y=indirect).  A clean axes background is captured and the X artists are
        then blitted where supported; a normal draw is used as a safe fallback.
        """
        if not hasattr(self, 'axes') or not hasattr(self, 'canvas'):
            return

        names = {str(name) for name in (peak_names or [])}
        refs = [pk for pk in self._reference_peaks()
                if str(getattr(pk, 'name', '')) in names]

        # Remove the previous transient selection before capturing the clean
        # contour background.  A full draw here also makes this robust to a
        # contour/slice redraw that occurred since the previous selection.
        for artist in getattr(self, '_fitting_selection_2d_artists', []):
            try:
                artist.remove()
            except (ValueError, RuntimeError):
                pass
        self._fitting_selection_2d_artists = []

        try:
            self.canvas.draw()
            background = self.canvas.copy_from_bbox(self.axes.bbox)
            for pk in refs:
                artist, = self.axes.plot(
                    [float(pk.x)], [float(pk.y)],
                    marker='x', linestyle='None', markersize=9,
                    markeredgecolor='black', color='black',
                    markeredgewidth=1.5, zorder=30, animated=True)
                self._fitting_selection_2d_artists.append(artist)
                self.axes.draw_artist(artist)
            self.canvas.blit(self.axes.bbox)
            self._fitting_selection_2d_background = background
        except (AttributeError, NotImplementedError, RuntimeError):
            # Some wx/Matplotlib backends do not expose usable blitting.
            for artist in self._fitting_selection_2d_artists:
                try:
                    artist.set_animated(False)
                except AttributeError:
                    pass
            self.canvas.draw_idle()


    def on_fitall(self,event):
        self._set_fitting_status('Fitting all peaks...')
        self.SetAllPeaks()
        # Fit All is one FUDA job: peak.fuda contains every reference peak and
        # param.fuda contains every applicable overlap group.
        self.DoFuda(backGrd='y')
        # Completion refresh/drawing is performed by _on_fuda_finished().




    def on_clean(self,event):
        self._set_fitting_status('Cleaning fitting results...')
        import shutil
        shutil.rmtree(self._fuda_dir(), ignore_errors=True)
        self._set_fitting_status('Fitting results cleaned')

        
    def on_norm_button(self, event):
        self.plot_scatters()

    

    def plot_scatters(self, axis='default'):
        if axis == 'default':
            axis = self.axes_proj
        axis_font = {'fontname':'Arial', 'size':'14'}
        for x in self.scatters:
            x.remove()
        axis.cla()
        fitter = self._ensure_line_fitting()
        self.intensities = fitter.intensities
        self.scatters=[]
        self.Js={}
        data = []
        names = []
        number_scatters = 0
        norm_button = getattr(self, 'norm_button', None)
        if norm_button is not None and norm_button.GetValue():
            data = self.scatter_data_norm
        else:
            for key in self.intensities.keys():
                data.append(self.intensities[key][0])
                names.append(key)

        xs= range(self._pseudo3d_view()['data'].shape[0])
        for num, line in enumerate(data):
            if len(self.gzlvl1)==self._pseudo3d_view()['data'].shape[0]:
                xs=numpy.array(self.gzlvl1)
                line = numpy.log(line)
                self.scatters.append(axis.scatter(xs, line,color='C'+str(number_scatters), marker='x'))
                m,b = numpy.polyfit(xs, line, 1)
                axis.plot(xs, m*xs+b, ls='--', color='C'+str(number_scatters))
                axis.text(0.1,0.9-float(number_scatters)*0.05, '$Deff = $%.2e $cm^2 s^-1$' % m, transform=axis.transAxes, color='C'+str(number_scatters),  **axis_font)
            elif len(self.T1s) == self._pseudo3d_view()['data'].shape[0]:
                xs=numpy.array(self.T1s)
                sorted = numpy.argsort(xs)
                xs = xs[sorted]
                line=line[sorted]
                line = numpy.log(line)
                r,i0 = numpy.polyfit(xs, line, 1)
                i0=numpy.exp(i0)
                axis.plot(xs, i0*numpy.exp(r*xs), ls='--', color='C'+str(number_scatters))
                self.scatters.append(axis.scatter(xs, numpy.exp(line),color='C'+str(number_scatters), marker='x'))
                axis.text(0.1,0.9-float(number_scatters)*0.05, r'$R_1 = %.2f s^{-1}$' % (-r), transform=axis.transAxes, color='C'+str(number_scatters),  **axis_font)
            elif len(self.T2s) == self._pseudo3d_view()['data'].shape[0]:
                xs=numpy.array(self.T2s)
                sorted = numpy.argsort(xs)
                xs = xs[sorted]

                line=line[sorted]
                line = numpy.log(line)
                r,i0 = numpy.polyfit(xs, line, 1)
                i0=numpy.exp(i0)
                
                self.selected_array = [self.selected]
               

                if str(names[num]) not in self.selected_array:
                    self.scatters.append(axis.scatter(xs, numpy.exp(line), marker='x', alpha=0.5,color='gray'))
                    axis.plot(xs, i0*numpy.exp(r*xs), ls='--', color='gray', alpha=0.5)
                    
                    if self.plotting_atom_results==True and str(names[num]) != 'prelim':
                        resi=int(re.findall(r'[0-9]+',str(names[num]))[0])

                        self.axes_atom_results.bar(resi, -r, color='gray', edgecolor='k', width=1.0, linewidth=0.5, picker=False)

                    if len(names) ==1:
                      axis.text(0.1,0.9-float(number_scatters)*0.05, r'$R_2 = %.2f s^{-1}$' % (-r), transform=axis.transAxes, color='lightgray', alpha=0.5,  **axis_font)
                else:
                    self.scatters.append(axis.scatter(xs, numpy.exp(line), marker='x', color='r', zorder=10000))
                    axis.plot(xs, i0*numpy.exp(r*xs), ls='--', color='r', zorder=10000)
                    if self.plotting_atom_results==True and str(names[num]) != 'prelim':
                        resi=int(re.findall(r'[0-9]+',str(names[num]))[0])

                        self.axes_atom_results.bar(resi, -r, color='r', edgecolor='k', width=1.0, linewidth=0.5, picker=False)

                    axis.text(0.1,0.9-float(number_scatters)*0.05, r'$R_2 = %.2f s^{-1}$' % (-r), transform=axis.transAxes, color='r',  **axis_font, zorder=10000)

            elif len(self.taus) == self._pseudo3d_view()['data'].shape[0]:
                xs=numpy.array(self.taus)
                sorted = numpy.argsort(xs)
                xs = xs[sorted]
                line=line[sorted]


                ys, J = self.fit_sine_curve(xs, line)
                self.Js[names[num]] = J

                self.selected_array = [self.selected]


                if str(names[num]) not in self.selected_array:
                    self.scatters.append(axis.scatter(xs, line, marker='x', alpha=0.5,color='gray'))
                    axis.plot(xs, ys, ls='--', color='gray',alpha=0.5)
                    if self.plotting_atom_results==True and str(names[num]) != 'prelim':
                        resi=int(re.findall(r'[0-9]+',str(names[num]))[0])

                        self.axes_atom_results.bar(resi, J, color='gray', edgecolor='k', width=1.0, linewidth=0.5, picker=False)

                    if len(names) ==1:
                      axis.text(0.75,0.9-float(number_scatters)*0.05, str(names[num])+r' J = %.2f Hz' % (J), transform=axis.transAxes, color='lightgray', alpha=0.5,  **axis_font)
                else:
                    self.scatters.append(axis.scatter(xs, line, marker='x', color='r', zorder=10000))
                    axis.plot(xs, ys, ls='--', color='r', zorder=10000)
                    if self.plotting_atom_results==True and str(names[num]) != 'prelim':
                        resi=int(re.findall(r'[0-9]+',str(names[num]))[0])

                        self.axes_atom_results.bar(resi, J, color='r', edgecolor='k', width=1.0, linewidth=0.5, picker=False)

                    axis.text(0.75,0.9-float(1)*0.05, str(names[num])+r' J = %.2f Hz' % (J), transform=axis.transAxes, color='r',  **axis_font, zorder=10000)

                


                self.print_Js()
            elif len(self.nu_CPMG) == self._pseudo3d_view()['data'].shape[0]:
                xs=numpy.array(self.nu_CPMG)
                mask0 = numpy.argwhere(xs==0.0)
                intensity_zero = numpy.average(line[mask0])
                mask = numpy.argwhere(xs!=0.0)
                ys = -self.time_T2*(numpy.log(line/intensity_zero))
                self.scatters.append(axis.scatter(xs[mask], ys[mask],color='C'+str(number_scatters), marker='x'))
                

            else:
                self.scatters.append(axis.scatter(xs, line,color='C'+str(number_scatters), marker='x'))
            number_scatters +=1
        self.canvas.draw()

    def print_Js(self):
        filename = 'out/J_values.out'
        outy = open(filename, 'w')
        outy.write("Peak Name\tJ (Hz)\n")

        for key in self.Js.keys():
            outy.write("%s\t%f\n" % (key, self.Js[key]))
            # outy.write("\n")
        outy.close()

    def fit_sine_curve(self, xs, line):
        initial_guess = (numpy.max(numpy.fabs(line))*1.3, 15., 0)
        popt, pcov, infodict, mesg, ier = opt.curve_fit(sine_function, xs, line, p0 = initial_guess, maxfev=100000, ftol=1.49012e-12, xtol=1e-12,gtol=1e-12, full_output=1)
        amp, J, c = popt
        return sine_function(xs, amp, J, c), J

    def on_pick(self,event):
        print(event.mouseevent.inaxes)
        if event.mouseevent.inaxes==self.axes:
            ind = event.ind
            print('picked:', self.peaks_text[ind[0]].get_text())
            self.selected = self.peaks_text[ind[0]].get_text()
            self.fuda_number = 0
            self._ensure_line_fitting()

            for key in self.line_fitting.plotting_resim_data.keys():
                if self.selected == str(key):
                    self.line_fitting.plot_fuda_fit(str(self.selected), 0, self.axes_proj, self.canvas)
            self.plot_scatters(self.axes_scatter)
            
        # if event.mouseevent.inaxes==self.axes_atom_results:
        #     print('picked:', event.artist)
            
    def on_click(self, event):
        if event.inaxes==self.axes_atom_results:
            # print(numpy.round(event.xdata))
            self.selected=self.residues[numpy.round(event.xdata)]
            self.plot_scatters(self.axes_scatter)
            self.fuda_number = 0
            self._ensure_line_fitting()
            self.line_fitting.plot_fuda_fit(str(self.selected), 0, self.axes_proj, self.canvas)



    def on_show_calc(self, event):
        """Blit the matching deconvolved plane without rebuilding contours."""
        visible = bool(self.cb_calc.GetValue())
        artists = getattr(self, 'calc_artists', [])
        if visible and not artists:
            self.cb_calc.SetValue(False)
            if hasattr(self, 'toolbar'):
                self.toolbar.set_decon_active(False)
            return
        for artist in artists:
            artist.set_visible(visible)
        background = getattr(self, '_calc_background', None)
        if background is not None:
            try:
                self.canvas.restore_region(background)
                if visible:
                    for artist in artists:
                        self.axes.draw_artist(artist)
                self.canvas.blit(self.axes.bbox)
                return
            except (AttributeError, NotImplementedError, RuntimeError):
                pass
        self.canvas.draw_idle()

    def on_cb_grid(self, event):

        if self.peaks_drawn == False:
            self.canvas.mpl_connect('pick_event', self.on_pick)
            peak_locs = []
            for point in self._reference_peak_overlay():
                loc1 = float(point['x'])
                loc2 = float(point['y'])
                lab = point.get('label', '')
                self.peaks_text.append(self.axes.text(loc1-0.04, loc2-0.04, lab, fontsize=12))
                peak_locs.append([loc1, loc2])
            peak_locs = ensure_xy_points(peak_locs)
            self.peaks_scatter.append(scatter_xy_points(self.axes, peak_locs, c='k',s=50,zorder=2,marker='x', picker=True, pickradius=5))
            self.peaks_drawn = True
        else:
            self.peaks_drawn=False
            for x in self.peaks_scatter:
                x.set_visible(False)
                x.remove()
            for x in self.peaks_text:
                x.set_visible(False)
                x.remove()

        self.canvas.draw()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure()
        # print(self.fig.dpi)
        # exit()
        self.canvas = FigCanvas(self, -1, self.fig)
        self._key_press_cid = self.canvas.mpl_connect('key_press_event', self.on_key)
        self.axes = self.fig.add_subplot(121)
        self.axesFuda = self.fig.add_subplot(122,projection='3d')
        #self.axesFuda = fig.gca(projection='3d')
        
        # self.axes_proj = self.fig.add_subplot(122)
        self.cursor_shown = False
        self.number_scatters=0
        self.pressed = False
        self.moved = False
        self.rectangles = []
        self.scatters = []
        self.scatter_data = []
        self.scatter_data_norm = []
        self.verticals = []
        self.not_yet_drawn = True
        self.gzlvl1 = []
        self.T1s = []
        self.T2s = []
        self.nu_CPMG = []
        self.taus = []
        self.peak_fitted_input=False
        self.peak_fitting_input=False
        self.fitted_unoverlapped=0
        self.current_slice=0
        self.unoverlapped=[]
        self.Gamma_x = 0.02
        self.Gamma_y = 0.2
        self.sigma_y = 0.2
        self.sigma_x = 0.02
        self.nu1 = 0.2
        self.nu2 = 0.2
        self.peaks_text=[]
        self.peaks_scatter=[]
        self.selected = ''
        self.line_fitting = None
        self._fitting_selection_2d_artists = []
        self._fitting_selection_2d_background = None

        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        # Build the modeless palettes first; their controls remain attributes of
        # this Pseudo3D panel and all handlers remain methods of this panel.
        self.contour_box()
        self.fitting_box()
        self.anal_box()
        self.drawing_box()

        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view, peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon, contour_callback=self._toolbar_contours, coordinates=False)
        self.toolbar.AddSeparator()
        # A single contiguous horizontal row: native Matplotlib tools, drawing
        # controls, separator, then the three modeless palette buttons.
        for widget in (self.up_button2d, self.down_button2d, self.slice_box):
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)
        self.toolbar.AddSeparator()
        self.fittingToolButton = wx.Button(self.toolbar, -1, 'Fitting', size=(-1,22))
        self.analysisToolButton = wx.Button(self.toolbar, -1, 'Analysis', size=(-1,22))
        self.fittingToolButton.Bind(wx.EVT_BUTTON, lambda evt: self.show_fitting_window())
        self.analysisToolButton.Bind(wx.EVT_BUTTON, lambda evt: self.show_analysis_selector())
        for widget in (self.fittingToolButton, self.analysisToolButton):
            self.toolbar.AddControl(widget)
        self.toolbar.bind_control_status_help(self.up_button2d, 'Next pseudo-3D slice')
        self.toolbar.bind_control_status_help(self.down_button2d, 'Previous pseudo-3D slice')
        self.toolbar.bind_control_status_help(self.slice_box, 'Enter pseudo-3D slice number')
        self.toolbar.bind_control_status_help(self.fittingToolButton, 'Open fitting controls')
        self.toolbar.bind_control_status_help(self.analysisToolButton, 'Open analysis controls')
        if self._is_single_plane_2d():
            # There is no experimental pseudo axis in a true 2D spectrum.
            # Keep contouring/peak fitting available, but suppress controls
            # whose meaning requires multiple pseudo planes.
            self.up_button2d.Enable(False)
            self.down_button2d.Enable(False)
            self.slice_box.Enable(False)
            self.analysisToolButton.Enable(False)
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 2)
        self.SetSizer(self.vbox)
        self.Layout()
        # Establish the compact geometry before the first spectral draw.
        # Normal redraws deliberately do not recompute subplot positions.
        self._apply_plot_layout()

        self.set_default_values() #load defaults.

    def _projection_view(self, l1, l2, transpose='n'):
        """Retrieve a plotting-ready projection from the shared controller/store."""
        getter = getattr(self.pseudo_service, 'projection_view', None)
        if getter is None:
            return None
        return getter(l1, l2, decon=False, transpose=transpose)

    def _apply_plot_layout(self, fitting=False):
        """Apply a stable compact plot geometry.

        ``tight_layout`` was previously run during ordinary redraws.  Because
        its result depends on the current artists and the realised canvas
        size, the axes could jump noticeably when Draw, slice navigation, or
        a contour update caused a repaint.  Keep the compact geometry as an
        explicit property of the Pseudo3D view instead.
        """
        if fitting:
            # The fitting view contains the main spectrum plus three result
            # axes, so retain a little more inter-axis spacing while keeping
            # the same compact outer margins.
            self.fig.subplots_adjust(
                left=0.065, right=0.985, bottom=0.085, top=0.975,
                wspace=0.34, hspace=0.28,
            )
        else:
            self.fig.subplots_adjust(
                left=0.065, right=0.985, bottom=0.085, top=0.975,
                wspace=0.16, hspace=0.16,
            )

    def draw_figure(self):
        self.draw_figureGO()

    def setup_fitting_view(self):
        self.axes.remove()
        try:
            self.axes_proj.remove()
            self.axes_atom_results.remove()
            self.axes_scatter.remove()
        except:
            pass
        self.axes = self.fig.add_subplot(121)
        self.draw_figure()
        gs = self.fig.add_gridspec(ncols = 4, nrows = 2)
        self.axes_proj = self.fig.add_subplot(gs[0,2], projection='3d')
        self.axes_atom_results = self.fig.add_subplot(gs[0,3])
        self.axes_atom_results.spines['right'].set_visible(False)
        self.axes_atom_results.spines['top'].set_visible(False)
        self.canvas.mpl_connect('button_release_event', self.on_click)

        

        self.setup_atom_results_view()
        self.axes_scatter = self.fig.add_subplot(224)
        self.axes_scatter.spines['right'].set_visible(False)
        self.axes_scatter.spines['top'].set_visible(False)
        self._apply_plot_layout(fitting=True)


    def setup_atom_results_view(self):
        self.residues = {}
        for point in self._reference_peak_overlay():
            name = point.get('label', '')
            try:
                resi=int(re.findall(r'[0-9]+',name)[0])
                self.residues[resi] = name
            except (IndexError, ValueError):
                print("Cannot extract residue number:", name)
                continue

        if not self.residues:
            self.plotting_atom_results = False
            return
        self.axes_atom_results.set_xlim(min(self.residues.keys()), max(self.residues.keys())+0.5)
        self.axes_atom_results.set_xlabel('Peak Number')
        self.plotting_atom_results = True
        


    def line_select_callback(self, eclick, erelease):
        """
        Callback for line selection.

        *eclick* and *erelease* are the press and release events.
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        x_a = int(self._pseudo3d_view()['x_uc'].f(str(eclick.xdata)+' ppm'))
        x_b = int(self._pseudo3d_view()['x_uc'].f(str(erelease.xdata)+' ppm'))
        y_a = int(self._pseudo3d_view()['y_uc'].f(str(eclick.ydata)+' ppm'))
        y_b = int(self._pseudo3d_view()['y_uc'].f(str(erelease.ydata)+' ppm'))
        y_1 = min(y_a, y_b)
        y_2 = max(y_a, y_b)
        x_1 = min(x_a, x_b)
        x_2 = max(x_a, x_b)

        count = 0

        lower_x = min(x1, x2)
        lower_y = min(y1, y2)
        higher_x = max(x1, x2)
        higher_y = max(y1, y2)

        bottom_right = (lower_x, lower_y)
        top_left = (higher_x, higher_y)
        

        for point in self._reference_peak_overlay():
            x_peak, y_peak = float(point['x']), float(point['y'])
            if x_peak < top_left[0] and x_peak > bottom_right[0] and y_peak < top_left[1] and y_peak > bottom_right[1]:
                if count != 0:
                    print('This box contains more than one peak!')
                    return
                count += 1
                self.selected = point.get('label', '')

        # exit()

        self.number_scatters+=1

        # self.axes_proj.set_axis_off()
        # self.axes_proj.remove()
        self.setup_fitting_view()
        data_coords = (y_1,y_2,x_1,x_2)

        self.fuda_number = 0
        fitter = self._ensure_line_fitting()
        result = fitter.prelim_fuda_thread(data_coords, self.selected)
        self.scatter_data.append(result)
        # self.selected = 'prelim'
        # self.intensities['prelim'] = [result, 0,0]
        fitter.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)
        self.fuda_number += 1
        self.plot_scatters(self.axes_scatter)

        if self.peak_fitting_input == True:
            self.peak_fitted_input = True
            #self.info_text.set_text("Thanks!  How does the fit look?")

    def draw_figureGO(self):
        """ Redraws the figure
        """
        self.axes.clear()
        view = self._pseudo3d_view()
        self.thresh=float(numpy.max(numpy.abs(view['data'][0]))*self.pseudo_service.threshold_fraction())


        #levels = [self.thresh]
        #levels=self.

        #self.text1=wx.StaticText(self, -1, 'Min:')
        #self.text2=wx.StaticText(self, -1, 'Factor:')
        #self.text3=wx.StaticText(self, -1, 'Number:')
        #self.textbox0 = wx.TextCtrl(self, size=(100,22), style=wx.TE_PROCESS_ENTER)
        #self.textbox1 = wx.TextCtrl(self, size=(40,22), style=wx.TE_PROCESS_ENTER)
        #self.textbox2 = wx.TextCtrl(self, size=(40,22), style=wx.TE_PROCESS_ENTER)


        contour_min = float(self.textbox0.GetValue())
        # 2D fitting should initially show the authoritative main spectrum.
        # A project may contain contour settings saved for an older pseudo3D
        # dataset; if that minimum is outside this spectrum, reset it to the
        # normal threshold rather than producing an apparently empty pane.
        slice_max = float(numpy.max(numpy.abs(view['data'][self.current_slice])))
        if self._is_single_plane_2d() and (contour_min <= 0 or contour_min >= slice_max):
            contour_min = max(slice_max * self.pseudo_service.threshold_fraction(), numpy.finfo(float).eps)
            self.textbox0.SetValue(str(contour_min))
        levels = [contour_min]
        for x in range(int(float(self.textbox2.GetValue()))):
            levels.append(levels[-1]*float(self.textbox1.GetValue()))

        levels = numpy.array(levels)




        cmap = plt.get_cmap('Oranges')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        self.axes.set_xlabel(view['x_label'], fontsize=8)
        self.axes.set_ylabel(view['y_label'], fontsize=8)
        XX1, YY1 = view['XX'], view['YY']

        # The shared projection store now returns a canonical plotting view:
        # XX/YY have the same orientation as a physical 3p spectral slice
        # (pseudo, indirect, direct).  The legacy Pseudo3D viewer transposed
        # the slice because the old projection path built the mesh in the
        # opposite orientation.  Applying that transpose to the canonical
        # view changes (480, 539) -> (539, 480) and makes Matplotlib reject
        # the contour.
        zs = numpy.asarray(view['data'][self.current_slice])

        # Keep this failure diagnostic explicit: if a non-canonical 3p file
        # reaches this viewer we want to report the orientation problem rather
        # than silently transpose scientific data.
        if zs.shape != XX1.shape or zs.shape != YY1.shape:
            raise RuntimeError(
                'Pseudo3D spectral slice/projection orientation mismatch: '
                'slice=%s XX=%s YY=%s labels=%s' %
                (zs.shape, XX1.shape, YY1.shape, view['labb'])
            )

        self.axes.contour(XX1,YY1,zs, levels=levels, colors='red', linewidths=0.5)
        self.axes.contour(XX1,YY1,-zs, levels=levels, colors='blue', linewidths=0.5)

        # Build the matching calculated plane once per spectral redraw.  Mark
        # these collections animated so the normal draw produces a clean
        # background that can be reused when the toolbar toggle changes.
        self.calc_artists = []
        decon_view = self._decon_pseudo3d_view()
        if decon_view is not None and self.current_slice < decon_view['data'].shape[0]:
            dz = numpy.asarray(decon_view['data'][self.current_slice])
            if dz.shape == XX1.shape:
                for values in (dz, -dz):
                    cs = self.axes.contour(XX1, YY1, values, levels=levels,
                                           colors='green', linewidths=0.7)
                    # Matplotlib <=3.7 exposed individual contour
                    # collections; newer releases make the QuadContourSet
                    # itself the Artist.  Support both APIs so the toolbar has
                    # a real object to show/hide in either case.
                    collections = getattr(cs, 'collections', None)
                    artists = list(collections) if collections is not None else [cs]
                    for artist in artists:
                        artist.set_animated(True)
                        artist.set_visible(bool(self.cb_calc.GetValue()))
                        self.calc_artists.append(artist)
        if hasattr(self, 'toolbar'):
            # Availability follows the shared decon spectrum, not a stale set
            # of artists from an earlier draw.  This matters for the one-plane
            # 2D adapter, where deconvolution may be run after the window opens.
            decon_available = decon_view is not None
            self.toolbar.enable_decon(decon_available)
            if not decon_available:
                self.cb_calc.SetValue(False)
                self.toolbar.set_decon_active(False)

        """
        self.selector = RectangleSelector(self.axes, self.line_select_callback,
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # disable middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=False)
        """

        
        self.axes.set_ylim(YY1[0][0], YY1[-1][0])
        self.axes.set_xlim(XX1[0][0], XX1[0][-1])
        #self.info_text = self.axes.text(0.1,0.9,'Hello', transform=self.axes.transAxes)
        self.canvas.draw()
        try:
            self._calc_background = self.canvas.copy_from_bbox(self.axes.bbox)
            if self.cb_calc.GetValue():
                for artist in self.calc_artists:
                    if artist.get_visible():
                        self.axes.draw_artist(artist)
                self.canvas.blit(self.axes.bbox)
        except (AttributeError, NotImplementedError, RuntimeError):
            self._calc_background = None


    def DoFuda(self,backGrd='n'):

        fuda_dir = self._fuda_dir()
        os.makedirs(fuda_dir, exist_ok=True)
        fudapeakfile = self._fuda_peak_file()
        fudafile = self._fuda_parameter_file()
        
        outy=open(fudapeakfile,'w')
        for f in self.fudapeaks:
            for peak in self._reference_peaks():
                #print('fff',peak.name,f)
                if(peak.name==f):
                    # FUDA spectral-axis convention is transposed relative to
                    # the system reference peak object: write name, F1(y), F2(x).
                    outy.write('%s %f %f\n' % (peak.name, peak.y, peak.x))
        outy.close()
        
        outy=open(fudafile,'w')
        outy.write('PEAKLIST=%s\n' % fudapeakfile)
        #outy.write('SPECFILE=raw/test.ft2\n') #bbbbbb
        outy.write('SPECFILE=%s\n' % self.pseudo_service.spectrum_path()) 
        outy.write('NOISE=3263.0\n')

        view = self._pseudo3d_view()
        # nmrPipeFit has a native non-arrayed 2D path.  Do not fabricate a
        # pseudo coordinate or temporary 3D file: pass the original .ft2 and
        # select that path explicitly with ZCOOR=2D.
        zcoor = '2D' if view.get('physical_dim') == 2 else view['pseudo_label']
        outy.write('ZCOOR=%s\n' % zcoor)
        #outy.write('DELAYFACTOR=1.000\n')
        #outy.write('BASELINE=N\n')
        outy.write('VERBOSELEVEL=5\n')
        outy.write('PRINTDATA=Y\n')

        self.MaxIter=2000
        self.FudaTol=1E-9
        outy.write('LM=(MAXFEV=%i;TOL=%e)\n' % (self.MaxIter,self.FudaTol))

        #outy.write('#DISCARD_SLICES=(1)\n')
        #outy.write('#BASELINE=Y\n')
        #
        #Specify the default values. All values are in ppm:
        #
        outy.write('DEF_LINEWIDTH_F2=%f\n' % (float(self.radF2.GetValue())/2.))
        outy.write('DEF_LINEWIDTH_F1=%f\n' % (float(self.radF1.GetValue())/2.))
        outy.write('DEF_RADIUS_F2=%f\n' % float(self.radF2.GetValue()))
        outy.write('DEF_RADIUS_F1=%f\n'% float(self.radF1.GetValue()))
        outy.write('SHAPE=GLORE\n')
        outy.write('ISOTOPESHIFT=N\n')
        #
        ##

        GrpsCpy=copy.deepcopy(self._groups())

        GrpInc=[]
        for i,pk in enumerate(self.fudapeaks):
            for grp,vals in GrpsCpy.items():
                if(pk in vals):
                    if(grp not in GrpInc):
                        GrpInc.append(grp)
                    break

        self._set_fitting_status('Preparing fit for %d peak%s...' % (len(self.fudapeaks), '' if len(self.fudapeaks) == 1 else 's'))
        for grp in GrpInc:
            outy.write('OVERLAP_PEAKS=(')
            for i,pk in enumerate(self._groups()[grp]):
                if(i!=0):
                    outy.write(';')
                outy.write(pk)
            outy.write(')\n')

        outy.close()

        # Run FUDA asynchronously and stream its combined stdout/stderr to a
        # progress window. Fit All remains one FUDA process and one parameter file.
        specargs = ['nmrPipeFit.py', fudafile, fuda_dir, 'unidecNMR']
        units = self._fuda_progress_units()
        progress = FudaProgressFrame(self.GetTopLevelParent(), units)
        self._fuda_progress_frame = progress
        progress.Show(True); progress.Raise()
        self._set_fitting_status('FUDA fitting in progress...')

        def worker():
            rc = -1
            try:
                proc = subprocess.Popen(specargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        universal_newlines=True, bufsize=1)
                for line in iter(proc.stdout.readline, ''):
                    if not line and proc.poll() is not None: break
                    wx.CallAfter(progress.append_line, line)
                proc.stdout.close(); rc = proc.wait()
            except Exception as exc:
                wx.CallAfter(progress.append_line, 'FUDA launch error: %s\n' % exc)
            wx.CallAfter(self._on_fuda_finished, rc, progress)
        threading.Thread(target=worker, name='FUDA-worker', daemon=True).start()
        return progress

    def _fuda_progress_units(self):
        """Build the same group/orphan work manifest represented by param.fuda."""
        selected=set(self.fudapeaks); units=[]; grouped=set()
        for grp, vals in self._groups().items():
            members=[p for p in vals if p in selected]
            if members:
                units.append({'group': str(grp), 'peaks': members}); grouped.update(members)
        for pk in self.fudapeaks:
            if pk not in grouped: units.append({'group': None, 'peaks': [pk]})
        return units

    def _on_fuda_finished(self, returncode, progress):
        progress.finish(returncode)
        if returncode != 0:
            self._set_fitting_status('FUDA fitting failed')
            return
        self._set_fitting_status('FUDA fitting complete')
        # Results now exist: refresh the combined table and draw the user's
        # current peak/group without blocking the wx event loop during fitting.
        self.RefreshFittingList()
        self.SetCurrentPeaks()
        result_files = []
        try:
            result_files = [name for name in os.listdir(self._fuda_dir())
                            if name.endswith(('.out', '.dat'))]
        except OSError:
            pass
        if result_files:
            self._mark_pseudo_intensities_ready(
                source='fuda', result_count=len(result_files), fit_directory=self._fuda_dir())
        if self.fudapeaks and os.path.isfile(os.path.join(self._fuda_dir(), self.fudapeaks[0] + '.dat')):
            try: self.ReadFuda()
            except Exception as exc: progress.append_line('Result display warning: %s\n' % exc)

    def _draw_fitting_slice(self, axes, slice_index):
        """Draw one parsed SpinUniDec pseudo slice using the GUI fit styling."""
        if slice_index < 0 or slice_index >= len(getattr(self, 'dat', [])):
            return
        for row in self.dat[slice_index]:
            x = numpy.array(row[:, 0])
            y = numpy.array(row[:, 1])
            zd = numpy.array(row[:, 2])
            zc = numpy.array(row[:, 3])
            axes.plot(x, y, zd, color='red')
            axes.plot(x, y, zc, color='green')
            if hasattr(self, 'cb_err') and self.cb_err.GetValue():
                axes.plot(x, y, zd-zc, color='blue')

    def _fitting_slice_z_limits(self):
        """Return common Z limits spanning every plotted fitting slice.

        Data and fitted intensities are always included.  When the fitting
        window is displaying the residual/error trace it is included too, so
        no report slice is clipped by the shared scale.
        """
        values = []
        show_error = hasattr(self, 'cb_err') and self.cb_err.GetValue()
        for pseudo_slice in getattr(self, 'dat', []):
            for row in pseudo_slice:
                try:
                    zd = numpy.asarray(row[:, 2], dtype=float)
                    zc = numpy.asarray(row[:, 3], dtype=float)
                except (IndexError, TypeError, ValueError):
                    continue
                values.extend((zd, zc))
                if show_error:
                    values.append(zd - zc)
        if not values:
            return None
        finite = numpy.concatenate([v[numpy.isfinite(v)] for v in values
                                    if numpy.asarray(v).size])
        if finite.size == 0:
            return None
        zmin = float(numpy.min(finite))
        zmax = float(numpy.max(finite))
        if zmin == zmax:
            pad = max(abs(zmin) * 0.05, 1.0)
            return zmin - pad, zmax + pad
        # Match Matplotlib's usual breathing room while keeping every slice
        # on exactly the same numerical scale.
        pad = 0.05 * (zmax - zmin)
        return zmin - pad, zmax + pad

    def _fitting_slice_label(self, slice_index):
        """Return the most informative available label for a pseudo-axis slice."""
        try:
            view = self._pseudo3d_view()
            for key in ('pseudo_axis', 'z_axis', 'pseudo_values'):
                values = view.get(key)
                if values is not None and len(values) > slice_index:
                    return 'Slice %d (%s)' % (slice_index + 1, values[slice_index])
        except Exception:
            pass
        return 'Slice %d' % (slice_index + 1)

    def ReadFuda(self):

        if(self.cb_accum.GetValue()==False):
            self.axesFuda.clear()

        for i in range(1):
            #for f in self.fudapeaks:
            f=self.fudapeaks[0]

            dat=[]
            r=[]
            d=[]
            cnt=0

            inny=open(os.path.join(self._fuda_dir(), f + '.dat'))
                
            for line in inny.readlines(): #for each line in file....
                test=line.split()  #split by whitespace....
                #print (test,cnt,len(test))
                if(len(test)>0):   #if there's at least one entry...
                    if(test[0][0]!='#'):   #if we're not a line that starts with a hash...
                        if(len(test)==4):   #if we have exactly 4 entries split by whitespace, this is our guy.
                            
                            if(cnt==1): # one blank line separates rows within a slice
                                if len(d)>0:
                                    r.append(numpy.array(d))
                                    d=[]
                            elif(cnt>1): # two or more blank lines separate Z slices
                                # Do not manufacture an empty row/slice when a writer
                                # emits an extra separator at a slice boundary.
                                if len(d)>0:
                                    r.append(numpy.array(d))
                                    d=[]
                                if len(r)>0:
                                    dat.append(r)
                                    r=[]
                            cnt=0                                

                                
                            x=float(test[0])
                            y=float(test[1])
                            yd=float(test[2])
                            yc=float(test[3])

                            row=numpy.array((x,y,yd,yc))
                            d.append(row)
                        else: #we have a line that is empty.
                            cnt+=1
                else:
                    cnt+=1
            if(len(d)>0): # add any trailing row
                r.append(numpy.array(d))
            if len(r)>0: # add the final slice even without a trailing separator
                dat.append(r)

            # Defensive cleanup for legacy/current writers: only non-empty
            # spectral slices count toward the pseudo/Z stack.
            dat = [plane for plane in dat if any(getattr(row, 'size', 0) > 0 for row in plane)]
            self.dat=dat

            sl=0
            #print('dat:',len(self.dat))
            #print('ffffffff',self.PeakComboSlice.Count)
            
            # SliceCombo represents the loaded pseudo/Z dimension.  FUDA should
            # produce the same number of slices; report a mismatch explicitly
            # rather than silently redefining the GUI from result-file length.
            expected_slices = self.PeakComboSlice.GetCount()
            if len(self.dat) != expected_slices:
                raise RuntimeError(
                    'FUDA result Z-stack length mismatch: result=%d loaded=%d' %
                    (len(self.dat), expected_slices))

            sl = self.PeakComboSlice.GetSelection()
            if sl < 0:
                sl = 0
                self.PeakComboSlice.SetSelection(0)

            
            self._draw_fitting_slice(self.axesFuda, sl)
            self._set_fitting_status('Showing %d fitted slice%s' % (len(self.dat), '' if len(self.dat) == 1 else 's'))
        self.canvas.draw()


    def SetCurrentPeaks(self):
        """Resolve the selected table peak to its complete overlap fit set."""
        currPeak = getattr(self, 'selected_fitting_peak', None)
        if not currPeak:
            self.fudapeaks = []
            return
        grp = self._peak_group_name(currPeak)
        self.fudapeaks = list(self._groups()[grp]) if grp is not None else [currPeak]
        self._set_fitting_status('Selected %d peak%s for fitting' % (len(self.fudapeaks), '' if len(self.fudapeaks) == 1 else 's'))


    def SetAllPeaks(self):
        #get current selection
        #first, clean up.
        self.on_clean(True)
        #now put all peaks into the list.
        self.fudapeaks=[]
        for pk in self._reference_peaks():
            self.fudapeaks.append(pk.name)
        self._set_fitting_status('Selected %d peak%s for fitting' % (len(self.fudapeaks), '' if len(self.fudapeaks) == 1 else 's'))

        
        
    def peak_fit(self, event):
        #print('peak_fitting not implemented yet')

        #######FUDA implementation#########
        #1. need list of peaks.
        #2.  need to be able to group peaks into overlap sets
        #3. for each set, get to set individual settings.
        #4. run fuda for each group, or run over everythin.
        #5. view results, versus Z slice.


        self.SetCurrentPeaks()
        self.DoFuda() # execute FUDA asynchronously; completion callback refreshes/draws
        
        
        #write fuda file.
        
        
        # Canonical fitting converters are view['y_uc'] and view['x_uc']

        # if self.peak_fitted_input == False:
        #     self.info_text.set_text('Please drag around an isolated peak!')
        #     self.canvas.draw()
        #     self.peak_fitting_input = True


        #             self.on_cb_grid(event)
        #             print("No peaks read in!")
        #             return

        #     print('aslkdjgf')
        #     return

        # thread = threading.Thread(target=self.fuda_thread)
        # thread.setDaemon(True)
        # thread.start()

        # self.canvas.draw()

    def peak_draw(self, event):  #redraw fuda slice
        #currPeak=self.PeakCombo.GetValue()
        #self.fudapeaks=[]
        #self.fudapeaks.append(currPeak)
        self.SetCurrentPeaks()        
        self.ReadFuda()
        
    def OnUpbutton(self, event): #increment fuda slice
        selection = self.PeakComboSlice.GetSelection()
        if selection < self.PeakComboSlice.GetCount() - 1:
            self.PeakComboSlice.SetSelection(selection + 1)
        self.peak_draw(event)
        
    def OnDownbutton(self, event): #decrement fuda slice
        if(self.PeakComboSlice.GetSelection()>0):
            self.PeakComboSlice.SetSelection(self.PeakComboSlice.GetSelection()-1)
        self.peak_draw(event)


        
    def fuda_thread(self):
        self._ensure_line_fitting().finding_overlaps()
        for point in self._reference_peak_overlay():
            peak = point.get('label', '')
            self.line_fitting.fit_unoverlapped_peak(peak)
            wx.CallAfter(self.plot_scatters, self.axes_scatter)


        self.overlap_thread()

    def overlap_thread(self):
        self._ensure_line_fitting()
        for x in range(len(self.line_fitting.final_overlaps)):
            self.line_fitting.fit_overlapped_peaks(x)
            wx.CallAfter(self.plot_scatters, self.axes_scatter)

    def on_key(self, event):
        if event.key in ('n', 'r') and self.selected != '':
            self._ensure_line_fitting()
        if event.key=='n':
            if self.selected != '':
                if self.fuda_number < len(self._pseudo3d_view()['data']) and self.fuda_number>-1:
                    try:
                        for key in self.line_fitting.plotting_resim_data.keys():
                            if str(self.selected) == str(int(key)):
                                # print('fuda_fit', self.selected)
                                self.line_fitting.plot_fuda_fit(key, self.fuda_number, self.axes_proj, self.canvas)
                                self.fuda_number+=1

                    except:
                        for key in self.line_fitting.plotting_resim_data.keys():
                            if str(self.selected) == str(key):
                                # print('fuda_fit', self.selected)
                                self.line_fitting.plot_fuda_fit(key, self.fuda_number, self.axes_proj, self.canvas)
                                self.fuda_number+=1
        if event.key=='p':
            self.axes_proj.cla()
            self.axes_proj.set_axis_off()
            self.axes_proj.remove()
            self.axes_proj = self.fig.add_subplot(122, projection=None)
            self.plot_scatters()

        if event.key =='r' and self.selected != '':
           

            overlapped = self.line_fitting.is_peak_overlapped(self.selected)
            if overlapped == None:
                return
            elif overlapped == False:  
                self.line_fitting.fit_unoverlapped_peak(self.selected)
                self.line_fitting.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)

            else:
                self.line_fitting.fit_overlapped_peaks(overlapped)
                self.line_fitting.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)

        if event.key=='enter' and self.peak_fitted_input == True:
            self.peak_fitting_input = False
            #self.info_text.set_text('Thanks: now going to iterate through')
            # self.peak_fit(event)
            # self.canvas.draw()
            # self.peak_fitted_input = False
        elif event.key=='enter':
            self.fit_unoverlapped_peak()

    def fit_unoverlapped_peak(self):
        if self.fitted_unoverlapped < len(self.unoverlapped):
            unover = self.unoverlapped[self.fitted_unoverlapped]

            fuda_data = self._pseudo3d_view()['data'][:,unover[0]:unover[1], unover[2]:unover[3]]
            if fuda_data.shape[1] == 0 or fuda_data.shape[2] == 0:
                self._set_fitting_status('Skipped an empty fitting region')
                self.fitted_unoverlapped += 1
                return
            fitter = self._ensure_line_fitting()
            peak_name = unover[8] if len(unover) > 8 else 'prelim'
            peak_loc = numpy.unravel_index(numpy.fabs(fuda_data[0]).argmax(), fuda_data[0].shape)
            fitter.fuda(fuda_data, peak_loc, name=str(peak_name))
            self.fitted_unoverlapped +=1
        else:
            self._set_fitting_status('All selected peaks fitted')

    def fit_ST_equation(self, event):
        ##############################################################################
        # Constants
        #
        # STEJSKAL & TANNER:
        # I/I0 = exp(-gamma^2 G^2 delta^2 [BigT-delta/3] Diff)
        # with
        # I                  signal intensity with diffusion weighting
        # I0                 signal intensity without diffusion weighting
        # gamma              gyromagnetic ratio of protons (rad s-1 G-1)
        gamma=2.675222E4
        # Gmax		     strength of the gradient pulse (G cm-1)
        Gmax=60.
        # delta              duration of the gradient pulse (s) (read in below)
        # BigT               time between the two gradient pulses (read in below)
        # Diff               diffusion constant (cm^2 s-1) (calculated below)
        #
        ##############################################################################
        # print(os.popen('cd raw && vpar seqfil && cd ..').read().split('\n')[2].replace('\"', ''))
        # print('hmqc_c13_500_methyl_diffusion_lek')
        if(os.popen('cd raw && vpar seqfil && cd ..').read().split('\n')[2].replace('\"', '').replace(' ','')=='hmqc_c13_500_methyl_diffusion_lek'):
            gzread = (os.popen('cd raw && vpar gzlvl5 && cd ..').read())
            delta = float(os.popen('cd raw && vpar gt5 && cd ..').read().split('\n')[2])
            BigT = float(os.popen('cd raw && vpar bigT && cd ..').read().split('\n')[2])
            # gamma = gamma*(5./6.)
            fields = gzread.split('\n')
            gzlvl1 = []
            if len(fields) > 2:
                gzs = fields[2].split(' ')
                for gz in gzs:
                    try:
                        if float(gz)>2:
                            G = Gmax*float(gz)/30000.


                            gzlvl1.append(-delta**2*G**2*gamma**2*(BigT-(1/3)*(delta)))

                    except:
                        pass
        else:
            gzread = (os.popen('cd raw && vpar gzlvl1 && cd ..').read())
            delta = float(os.popen('cd raw && vpar gt1 && cd ..').read().split('\n')[2])
            try:
                BigT = float(os.popen('cd raw && vpar BigT && cd ..').read().split('\n')[2])
            except:
                BigT = float(os.popen('cd raw && vpar bigT && cd ..').read().split('\n')[2])

            fields = gzread.split('\n')
            gzlvl1 = []
            if len(fields) > 2:
                gzs = fields[2].split(' ')
                for gz in gzs:
                    try:
                        if float(gz)>2:
                            G = Gmax*float(gz)/30000.

                            gzlvl1.append(-delta**2*G**2*gamma**2*(BigT-(1/3)*(delta)))
                    except:
                        pass


        self.axes.set_yticks(numpy.arange(len(gzlvl1)))
        self.axes.set_yticklabels(gzlvl1)
        self.gzlvl1 = gzlvl1

        self.plot_scatters()

    def fit_cpmg(self, event):
        ncyc_cp = (os.popen('cd raw && vpar ncyc_cp && cd ..').read())
        time_T2 = float(os.popen('cd raw && vpar time_T2 && cd ..').read().split('\n')[2])
        fields = ncyc_cp.split('\n')
        nu_CPMG = []
        if len(fields) > 2:
            ncycles = fields[2].split(' ')
            for r in ncycles:
                try:
                    # if float(r)>2:
                        ncyc = float(r)

                        # print(delta)
                        nu_CPMG.append(ncyc/time_T2 )
                except:
                    pass
        # self.axes.set_xticks(T1s)
        # self.axes.set_yticklabels(T1s)
        self.time_T2 = time_T2
        self.nu_CPMG = nu_CPMG 

        self.plot_scatters(self.axes_scatter)

    def fit_t1(self, event):

        gzread = (os.popen('cd raw && vpar ncyc && cd ..').read())
        fields = gzread.split('\n')
        T1s = []
        if len(fields) > 2:
            rs = fields[2].split(' ')
            for r in rs:
                try:
                    # if float(r)>2:
                        ncyc = float(r)

                        # print(delta)
                        T1s.append(ncyc*(2.0*12.5e-3))
                except:
                    pass
        # self.axes.set_xticks(T1s)
        # self.axes.set_yticklabels(T1s)
        self.T1s = T1s

        self.plot_scatters()

    def fit_j(self, event):
        msg = "What's the name of the variable inept delay?"
        dlg = wx.TextEntryDialog(None, msg)
        res = dlg.ShowModal()
        if res == wx.ID_CANCEL:
            return False
        j_name = dlg.GetValue()
        # self.first_drop = False
        try:
            gzread = (os.popen('cd raw && vpar '+j_name+' && cd ..').read())
        except:
            print('Variable not found')
            return False
        fields = gzread.split('\n')
        taus = []
        if len(fields) > 2:
            ts = fields[2].split(' ')
            for t in ts:
                try:
                        tau = float(t)
                        taus.append(tau)
                except:
                    pass
        # self.axes.set_xticks(T1s)
        # self.axes.set_yticklabels(T1s)
        self.taus = taus


        self.plot_scatters(self.axes_scatter)

    def fit_t2(self, event):

        gzread = (os.popen('cd raw && vpar time_T2 && cd ..').read())
        pwn = float(os.popen('cd raw && vpar pwn && cd ..').read().split('\n')[2])
        fields = gzread.split('\n')
        T2s = []
        if len(fields) > 2:
            rs = fields[2].split(' ')
            for r in rs:
                try:
                    ncyc = float(r)
                    # T2s.append(ncyc*(32.0*pwn*1E-6 + 32.0*450.0e-6))
                    # T2s.append(ncyc*(32.0*pwn*1E-6 + 32.0*450.0e-6))
                    T2s.append(ncyc)
                except:
                    pass
        self.T2s = T2s

        self.plot_scatters(self.axes_scatter)

    def on_scroll(self, event):
        self.ymin,self.ymax=self.axes_h.get_ylim()
        self.axes_h.set_ylim(self.ymin+(self.ymin*0.05*event.step), self.ymax+(self.ymax*0.05*event.step))
        self.axes_h.draw_artist(self.h_line)

    def on_mouse_move(self,event):
        if self.not_yet_drawn == True:

            self.background = self.canvas.copy_from_bbox(self.axes.bbox)

            self.current_h = 0
            self.v_line = self.axes_h.axvline(self._pseudo3d_view()['x_axis'][0], color = 'r', linewidth=2)
            self.h_line, = self.axes_h.plot(self._pseudo3d_view()['x_axis'], numpy.zeros_like(self._pseudo3d_view()['x_axis']), color='k', linewidth = 0.5)
            self.axes_h.set_ylim(numpy.min(self._pseudo3d_view()['data']), numpy.max(self._pseudo3d_view()['data']))


            self.canvas.draw()
            self.not_yet_drawn = False
        if event.inaxes == None:
            self.h_line.set_visible(False)
            # self.canvas.draw()
            self.axes_h.draw_artist(self.h_line)
        else:
            self.h_line.set_visible(True)
            # self.canvas.draw()
            self.axes_h.draw_artist(self.h_line)
        if self.pressed == True:
            self.moved = True
        if self.axes != event.inaxes:


            inv = self.axes.transData.inverted()
            new_dataPoint = int(inv.transform(numpy.array((event.x, event.y)).reshape(1, 2)).ravel()[1])
            # new_dataPoint = int(event.ydata) #(int(numpy.floor(self.combinedTransform.transform(pt_data2)[1])))
            self.canvas.restore_region(self.background)
            self.v_line.set_xdata(event.xdata)
            self.h_line.set_ydata(self._pseudo3d_view()['data'][self.current_slice, new_dataPoint, :])
            self.axes_h.draw_artist(self.h_line)
            self.axes_h.draw_artist(self.v_line)
            self.current_h = new_dataPoint
            # self.canvas.blit(self.axes.bbox)
            self.canvas.draw()

            # if new_dataPoint != self.current_h:

    def _toolbar_decon(self, active):
        self.cb_calc.SetValue(bool(active))
        # Rebuild the contour artists when necessary.  In particular a true
        # 2D fitting window can pre-date the deconvolution run, so it may have
        # no calculated artists even though the shared decon spectrum is now
        # available.  A redraw also guarantees raw contours remain red/blue
        # and the optional calculated overlay alone is green.
        decon_view = self._decon_pseudo3d_view()
        if active and decon_view is not None and not getattr(self, 'calc_artists', []):
            self.draw_figureGO()
            return
        if active and decon_view is None:
            self.cb_calc.SetValue(False)
            self.toolbar.set_decon_active(False)
            self.toolbar.enable_decon(False)
            return
        self.on_show_calc(None)

    def _toolbar_peaks(self, active):
        # Pseudo3D historically toggled artists rather than reading the checkbox.
        # Keep that fast overlay path while synchronising the retained state widget.
        self.cb_grid.SetValue(bool(active))
        if bool(active) != bool(self.peaks_drawn):
            self.on_cb_grid(None)

    def _toolbar_contours(self):
        self._show_tool_window(self.contourFrame)

    def redraw_view(self):
        self.on_draw_button(None)

    def on_draw_button(self, event):
        self.draw_figure()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()

    def on_N_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.draw_figure()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)

    def on_exit(self, event):
        self.Destroy()

    def on_about(self, event):
        msg = """ A demo using wxPython with matplotlib:

         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw!")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def flash_status_message(self, msg, flash_len_ms=1500):
        statusbar = getattr(self, 'statusbar', None) or getattr(self.nmr_workspace, 'statusbar', None)
        if statusbar is None:
            return
        self.statusbar = statusbar
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER,
            self.on_flash_status_off,
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')



class groupMan(wx.App):
    def __init__(self,inherit):
        self.frame_groupManFrame=groupManFrame(None,10,'Groups',inherit)
        self.frame_groupManFrame.Show(True)
#        return Frame1(parent)


# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class groupManFrame(wx.Frame):
#    title = 'AssBox'
    def __init__(self,parent, id, title,inherit):
        self.parent=inherit
        self._init_ctrls(parent,inherit)



    def _init_ctrls(self,prnt,parent):
        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=prnt,pos=wx.Point(358, 184), size=wx.Size(250, 20),
              style=wx.DEFAULT_FRAME_STYLE, title='Groups ...')
        self.SetClientSize(wx.Size(900, 280))

        panel=wx.Panel(self,-1)

        #self.prnt=prnt
        #print(self.prnt.peak)
        #self.parent=parent
        self.corrDict={}

        self.SelGrp=-1
        

        self.lcTxt=wx.StaticText(panel, -1, 'Orphans:')
        self.lc=SortedListCtrl(panel,self.corrDict)

        cnt=0
        self.lc.InsertColumn(cnt, 'Name');cnt+=1

        #self.lc.InsertColumn(cnt, self.parent.iLab.GetLabel());cnt+=1
        #self.lc.InsertColumn(cnt, self.parent.jLab.GetLabel());cnt+=1
        #self.lc.InsertColumn(cnt, 'ppmI(ppm)');cnt+=1
        #self.lc.InsertColumn(cnt, 'ppmJ(ppm)');cnt+=1


        #self.lc.SetColumnWidth(0, 140)
        #self.lc.SetColumnWidth(1, 153)

        self.lcInTxt=wx.StaticText(panel, -1, 'Group Contents:')
        self.lcIn=SortedListCtrl(panel,self.corrDict)

        cnt=0
        self.lcIn.InsertColumn(cnt, 'Name');cnt+=1
        #self.lc.InsertColumn(cnt, self.parent.iLab.GetLabel());cnt+=1
        #self.lc.InsertColumn(cnt, self.parent.jLab.GetLabel());cnt+=1
        #self.lcIn.InsertColumn(cnt, 'ppmI(ppm)');cnt+=1
        #self.lcIn.InsertColumn(cnt, 'ppmJ(ppm)');cnt+=1


        self.lcGrpTxt=wx.StaticText(panel, -1, 'Groups:')        
        self.lcGrp=SortedListCtrl(panel,self.corrDict)

        cnt=0
        self.lcGrp.InsertColumn(cnt, 'Name');cnt+=1

        #self.lc.InsertColumn(cnt, self.parent.iLab.GetLabel());cnt+=1
        #self.lc.InsertColumn(cnt, self.parent.jLab.GetLabel());cnt+=1
        self.lcGrp.InsertColumn(cnt, 'Members');cnt+=1

        self.NewGrpButton =  wx.Button(panel, -1, 'New Group',(710,10))
        self.RemGrpButton= wx.Button(panel, -1, 'Remove Group',(710,60))

        
        
        self.RefreshButton= wx.Button(panel, -1, 'Refresh',(710,10))
        self.ShowButton =  wx.Button(panel, -1, 'Show',(710,10))
        self.AddToGrpButton= wx.Button(panel, -1, '+',(710,60))
        self.RemFromGrpButton = wx.Button(panel, -1, '-',(710,160))

        self.CloseButton = wx.Button(panel, -1, 'Close',(710,160))

        self.GuessButton = wx.Button(panel, -1, 'Guess',(710,160))


        self.lcGrp.Bind(wx.EVT_LIST_ITEM_SELECTED, self.SelectGrp)
                

        self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnShow,self.lc)


        self.Bind (wx.EVT_BUTTON, self.OnNewGrp, self.NewGrpButton)
        self.Bind (wx.EVT_BUTTON, self.OnRemGrp, self.RemGrpButton)
        
        self.Bind (wx.EVT_BUTTON, self.OnRefresh, self.RefreshButton)
        self.Bind (wx.EVT_BUTTON, self.OnShow, self.ShowButton)
        self.Bind (wx.EVT_BUTTON, self.OnClose, self.CloseButton)
        self.Bind (wx.EVT_BUTTON, self.OnAddToGrpButton, self.AddToGrpButton)
        self.Bind (wx.EVT_BUTTON, self.OnRemFromGrpButton, self.RemFromGrpButton)

        self.Bind (wx.EVT_BUTTON, self.OnGuessButton, self.GuessButton)


        
        #self.vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        vbox1 = wx.BoxSizer(wx.VERTICAL)

        vbox1.Add(self.lcGrpTxt)
        vbox1.Add(self.lcGrp)



        
        hbox.Add(vbox1,1,wx.EXPAND)
        
        #hbox.Add(self.lcGrp, 1, wx.EXPAND)

        hbox.AddSpacer(10)

        vbox2 = wx.BoxSizer(wx.VERTICAL)
        vbox2.Add(self.lcInTxt)
        vbox2.Add(self.lcIn)

        hbox.Add(vbox2)
        
        hbox.AddSpacer(10)

        vbox3 = wx.BoxSizer(wx.VERTICAL)
        vbox3.Add(self.lcTxt)
        vbox3.Add(self.lc)

        hbox.Add(vbox3)
        

        hbox.AddSpacer(10)
        
        vbox=wx.BoxSizer(wx.VERTICAL)

        vbox.Add(self.NewGrpButton)
        vbox.Add(self.RemGrpButton)
        
        vbox.AddSpacer(10)


        vbox.Add(self.RefreshButton,  wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.ShowButton,  wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.CloseButton, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.GuessButton, wx.ALIGN_CENTER| wx.TOP)
        vbox.AddSpacer(10)
        vbox.Add(self.AddToGrpButton, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.RemFromGrpButton, wx.ALIGN_CENTER| wx.TOP)
        hbox.Add(vbox)
        panel.SetSizer(hbox)

        self.OnRefresh(True)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        #hbox  = wx.BoxSizer(wx.HORIZONTAL)
        #hbox.Add(self.Addbutton, 1, wx.EXPAND)
        #hbox.Add(self.Removebutton, 1, wx.EXPAND)
        #hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        #hbox.Add(self.Closebutton, 1, wx.EXPAND)
        #hbox.Add(self.Savebutton, 1, wx.EXPAND)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        self.Centre()

        hbox.Fit(self)
        self.Show(True)


        """
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        self.drawing_box()
        self.contour_box()
        self.fit_box = self.fitting_box()

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vbox2)
        self.hbox.AddSpacer(5)

        self.hbox.Add(self.cntrSizer)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.fit_box)

        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        """
        
        #hbox.Add(vbox2, 1, wx.EXPAND)
        #self.SetSizer(hbox)
        #        self.SetSizer(self.vbox)
        #        self.vbox.Fit(self)

    def GetNameGrp(self):
        go=1
        while(go==1):
            for n in range(1000):
                tig=0
                for k in self.parent._groups().keys():
                    if(str(n+1)==k):
                        tig=1
                        break
                if(tig==0):
                    return str(n+1)
        print('Could not allocate name')
            
                    
    def OnNewGrp(self,event):
        newGrp=self.GetNameGrp()

        #print('ttt',newGrp,self.parent.Grps.keys())
        #while(newGrp not in list(self.parent.Grps.keys())):
        #    n+=1
        #    newGrp=str(len(self.parent.Grps.keys())+n)
        #    print('ttt',newGrp,list(self.parent.Grps.keys()))
        
        print('Adding:',newGrp)
        self.parent._add_group(newGrp)
        self.OnRefresh(event)
    def OnRemGrp(self,event):
        c = self.lcGrp.GetFocusedItem()
        grp=self.lcGrp.GetItem(c, col=0).GetText()

        self.SelGrpInd=-1
        self.SelGrp=-1
        
        #put the contents back in the orphan drawer
        for pk in self.parent._groups()[grp]:
            self.orph[pk]=1
        
        affected = list(self.parent._groups().get(grp, []))
        self.parent._remove_group(grp)
        self.parent.invalidate_fits_for_peaks(affected)
        #del self.parent.Grp(str(c+1))

        #index = self.lc.GetFocusedItem()
        #self.lc.DeleteItem(index)
        #print(index)
        #self.parent.draw_figure()
        self.OnRefresh(True)



        
        pass

        
        
        
    def OnAddToGrpButton(self,event):

        c = self.lc.GetFocusedItem()
        ToAdd=self.lc.GetItem(c, col=0).GetText()


        self.orph[ToAdd]=0 #remove from orphans

        g = self.lcGrp.GetFocusedItem()
        grp =self.lcGrp.GetItem(g, col=0).GetText()


        if ToAdd not in self.parent._groups()[grp]:
            affected = list(self.parent._groups().get(grp, [])) + [ToAdd]
            self.parent._add_peak_to_group(grp, ToAdd)
            self.parent.invalidate_fits_for_peaks(affected)
            print('Adding',ToAdd,'to',grp)
        else:
            print('Already added',ToAdd,'to',grp)

        print (self.parent._groups())
        if(self.SelGrpInd!=-1):  
            self.lcGrp.SetItem(self.SelGrpInd,1,'('+str(len(self.parent._groups()[self.SelGrp]))+')') #add atom                

        self.OnRefreshRest(True)
        #
        pass

    def OnRemFromGrpButton(self,event):
        ToRem = self.lcIn.GetItem(self.lcIn.GetFocusedItem(), col=0).GetText()
        if ToRem in self.parent._groups().get(self.SelGrp, []):
            print('Removing',ToRem,'from grp',self.SelGrp)
            affected = list(self.parent._groups().get(self.SelGrp, []))
            self.parent._remove_peak_from_group(self.SelGrp, ToRem)
            self.parent.invalidate_fits_for_peaks(affected)

        if(self.SelGrpInd!=-1):
            self.lcGrp.SetItem(self.SelGrpInd,1,'('+str(len(self.parent._groups()[self.SelGrp]))+')') #add atom                

        self.OnRefreshRest(True)
                

    def onItemSelected(self, event):
        """"""

        #self.parent.select = []
        currentItem = event.GetIndex()
        print('cyrrIt:',currentItem)
        #car = self.corrDict[currentItem]
        #print(int(car[1]))
        #self.parent.select.append(currentItem)
        #self.parent.draw_figure()
        #self.parent.SELECT=0

        #count = self.lc.GetItemCount()
        #self.sorted_artists = [self.list.GetItem(itemId=row, col=0).GetText() for row in xrange(count)]
        #print self.sorted_artists
        #print self.sorted_artists[currentItem]


    def AtoI(self,val):
        for i, peak in enumerate(self.parent._reference_peaks()):
            if(val==peak.name):
                return i

    """
    def OnAdd(self, event):
        sele=self.lc.GetFirstSelected()
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row, col=0).GetText() for row in xrange(count)][sele]
        col2 = [self.lc.GetItem(row, col=1).GetText() for row in xrange(count)][sele]
        #print col1,col2,self.AtoI(col1)
        self.parent.ComboBox1.SetSelection(self.AtoI(col1))
        self.parent.ComboBox2.SetSelection(self.AtoI(col2))
        #self.parent.NOE=1

        self.parent.on_draw_button(True)
    """

    def OnRemove(self, event):
        print('Removing item')
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)


        print(index)
        #self.parent.draw_figure()
        self.OnRefresh(True)


    def OnClose(self, event):
        self.Close()



    def SelectGrp(self,event):
        print()
        try:
            print ('index of selection:',self.lcGrp.GetFocusedItem())
            self.SelGrpInd=self.lcGrp.GetFocusedItem()
            self.SelGrp=self.lcGrp.GetItem(self.lcGrp.GetFocusedItem(), col=0).GetText()

        except:
            self.SelGrpInd=-1
            self.SelGrp=-1
        print('Selected group:',self.SelGrp)
        self.OnRefreshRest(True)        

    def OnRefreshRest(self,event):
        print()
        print('Refreshing peak lists')
        print('Grps:',self.parent._groups())
        
        self.orph={}
        for peak in self.parent._reference_peaks():
            n=peak.name
            self.orph[n]=1


        self.lcIn.DeleteAllItems()
        num_items = self.lcIn.GetItemCount()
        for grp,pks in self.parent._groups().items():
            for pk in pks:
                self.orph[pk]=0
                if(self.SelGrp!=-1 and grp==self.SelGrp): #insert into current group contents list if group selection is right
                    self.lcIn.InsertItem(num_items,(pk))

        num_items = self.lc.GetItemCount()
        self.lc.DeleteAllItems()
        for key,vals in self.orph.items():
            if(vals==1):
                self.lc.InsertItem(num_items,key)  #add assignment
                #self.lc.SetItem(num_items, 0,n) #add atom

                
        self.UpdateFitList()


    def UpdateFitList(self):

        self.parent.SetPeaksToFit()

        

    
    def OnRefresh(self, event):

        print()
        print('Refreshing groups list')
        print('Grps:',self.parent._groups())
        self.lcGrp.DeleteAllItems()
        num_items = self.lcGrp.GetItemCount()
        cnt=0
        for grp,pks in self.parent._groups().items():
            ind=self.lcGrp.InsertItem(num_items,(grp))  #add group
            #self.lcGrp.SetItem(num_items, 0,grp) #add group
            self.lcGrp.SetItem(ind,1,'('+str(len(self.parent._groups()[grp]))+')') #add atom                
            cnt+=1
        self.OnRefreshRest(True)


    ###############################################################

    def OnGuessButton(self,event):
        radF1=float(self.parent.radF1.GetValue())
        radF2=float(self.parent.radF2.GetValue())
        self.FindOverlapGroups(radF1,radF2)
        print ("overlaps:",self.overlap)

        old_affected = [pk for vals in self.parent._groups().values() for pk in vals]
        groups = {str(i + 1): self.overlap[i] for i in range(len(self.overlap))}
        new_affected = [pk for vals in groups.values() for pk in vals]
        self.parent._replace_groups(groups)
        self.parent.invalidate_fits_for_peaks(old_affected + new_affected)
        self.OnRefresh(True)

        
    def AddIfNew(self,test,array):
        for i in range(len(array)):
            if(array[i]==test):
                return 0
        array.append(test)
        return 1


    def GetPeak(self,pk_name):
        for pk in self.parent._reference_peaks():
            if(pk_name==pk.name):
                return pk
        print("Failed to find a peak.")
        sys.exit(100)
        
    def GetNeighbours(self,pk_name,radF1,radF2):
        # Reference peaks are now dimension-neutral peakEntry objects with
        # canonical spectral coordinates ``y`` (F1) and ``x`` (F2).  The old
        # pseudo3D implementation reached into legacy 3D aliases ppmJ/ppmK,
        # which do not exist for ordinary 2D peak lists.
        pk=self.GetPeak(pk_name)
        f1=float(pk.y); f2=float(pk.x)
        min_f1=f1-radF1; max_f1=f1+radF1
        min_f2=f2-radF2; max_f2=f2+radF2

        for pkn in self.parent._reference_peaks():
            if pkn.name == pk.name:
                continue
            nf1=float(pkn.y); nf2=float(pkn.x)
            if min_f1 < nf1 < max_f1 and min_f2 < nf2 < max_f2:
                self.AddIfNew(pkn.name,self.grpTest)

    def FindOverlapGroups(self,radF1,radF2):
        self.overlap=[]
        for i,pk in enumerate(self.parent._reference_peaks()):
            self.grpTest=[]
            self.grpTest.append(pk.name)
            #grp=self.GetNeighbours(pk.name,radiusC,radiusH,grp)
            go=0
            while(go==0):
                start=len(self.grpTest)
                for p in self.grpTest:
                    self.GetNeighbours(p,radF1,radF2)
                if(len(self.grpTest)==start):#finished
                    go=1
            if(len(self.grpTest)>1):#if the group is complete
                self.AddIfNew(sorted(self.grpTest),self.overlap)


    ##################################################################




    def OnShow(self,event):
        #index = self.lc.GetFocusedItem()
        index=self.lc.GetFirstSelected()

        if index < 0:
            return
        peak_name = self.lc.GetItem(index, 0).GetText()
        point = next((p for p in self.parent._reference_peak_overlay()
                      if str(p.get('label', '')) == peak_name), None)
        if point is None:
            return
        print('Showing ', point.get('label', ''))
        axes = self.parent.axes
        view = self.parent._pseudo3d_view()

        x_axis = numpy.asarray(view['x_axis'])
        y_axis = numpy.asarray(view['y_axis'])
        widX = abs(float(x_axis[-1]) - float(x_axis[0])) / 10.
        widY = abs(float(y_axis[-1]) - float(y_axis[0])) / 10.
        x_peak = float(point['x'])
        y_peak = float(point['y'])
        axes.set_xlim(x_peak-widX, x_peak+widX)
        axes.set_ylim(y_peak-widY, y_peak+widY)
        self.parent.canvas.draw_idle()


        return
    def OnListBox1Listbox(self, event):
        '''
        click list item and display the selected string in frame's title
        '''
#        selName = self.listBox1.GetStringSelection()
#        self.SetTitle(selName)
        return

    def OnButton2Button(self, event):
        '''
        click button to clear the listbox items
        '''
        self.listBox1.Clear()

class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

    def Update(self,dicty):
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty
        #print(dicty[0])

    def CustColumnSorter(self, key1, key2):
        col = self._col
        ascending = self._colSortFlag[col]
        ascending=1
        item1 = self.itemDataMap[key1][col]
        item2 = self.itemDataMap[key2][col]

        self.num_cols=[0,2,3,]
        if col in self.num_cols:
            #just convert them to float, cmp do comparing float well
            item1 = float(item1)
            item2 = float(item2)

        cmpVal = cmp(item1, item2)

        # If the items are equal then pick something else to make the sort value unique
        if cmpVal == 0:
            cmpVal = cmp(*self.GetSecondarySortValues(col, key1, key2))

        if ascending:
            return cmpVal
        else:
            return -cmpVal

    def GetColumnSorter(self):
        return self.CustColumnSorter

    







class overMan(wx.App):
    def __init__(self,inherit):
        self.frame_groupManFrame=overManFrame(None,10,'Overview',inherit)
        self.frame_groupManFrame.Show(True)
#        return Frame1(parent)




class overManFrame(wx.Frame):
#    title = 'AssBox'
    def __init__(self,parent, id, title,inherit):
        self.parent=inherit
        self._init_ctrls(parent,inherit)



    def _init_ctrls(self,prnt,parent):
        # BOA generated methods
        wx.Frame.__init__(self, name='', parent=prnt,pos=wx.Point(358, 184), size=wx.Size(250, 20),
              style=wx.DEFAULT_FRAME_STYLE, title='Overview ...')
        self.SetClientSize(wx.Size(900, 280))

        panel=wx.Panel(self,-1)


        self.corrDict={}

        

        self.lc=SortedListCtrl(panel,self.corrDict)

        self.DolcCols()

        self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.Select)        
        #self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnShow,self.lc)

        self.lc.Bind(wx.EVT_LIST_COL_CLICK, self.OnButtonSort)
        

        self.CloseButton = wx.Button(panel, -1, 'Close',(710,160))
        self.CloseButton.Bind (wx.EVT_BUTTON, self.OnCloseButton)
        self.RefreshButton = wx.Button(panel, -1, 'Refresh',(710,160))
        self.RefreshButton.Bind (wx.EVT_BUTTON, self.OnRefreshButton)

        
        #self.vbox = wx.BoxSizer(wx.VERTICAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.lc)
        vbox.AddSpacer(10)

        hbox=wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.RefreshButton)
        hbox.Add(self.CloseButton)
        vbox.Add(hbox)


        vbox.Fit(self)
        panel.SetSizer(vbox)

        #self.OnRefresh(True)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        #hbox  = wx.BoxSizer(wx.HORIZONTAL)
        #hbox.Add(self.Addbutton, 1, wx.EXPAND)
        #hbox.Add(self.Removebutton, 1, wx.EXPAND)
        #hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        #hbox.Add(self.Closebutton, 1, wx.EXPAND)
        #hbox.Add(self.Savebutton, 1, wx.EXPAND)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)



        self.Centre()

        self.Populate()
        
        self.Show(True)

        """





        """
        #self.SetSizer(self.vbox)
        #self.vbox.Fit(self)

    def DolcCols(self):
        cnt=0
        self.lc.InsertColumn(cnt,'Peak');cnt+=1
        self.lc.InsertColumn(cnt,'%serr' % '%');cnt+=1
        self.lc.InsertColumn(cnt,'f01(ppm)');cnt+=1
        self.lc.InsertColumn(cnt,'w1(Hz)');cnt+=1
        self.lc.InsertColumn(cnt,'g1');cnt+=1
        self.lc.InsertColumn(cnt,'f02(ppm)');cnt+=1
        self.lc.InsertColumn(cnt,'w2(Hz)');cnt+=1
        self.lc.InsertColumn(cnt,'g2');cnt+=1

        print ('total columns:',cnt)


    def GetCol(self,rw,count):

        try:
            return numpy.array([float(self.lc.GetItem(row, rw).GetText()) for row in range(count)])
        except:
            #print(rw,count)
            #print([self.lc.GetItem(row, rw).GetText() for row in range(count)])
            return numpy.array([self.lc.GetItem(row, rw).GetText() for row in range(count)])

        
    def OnButtonSort(self,event):
        col=event.GetColumn()
        print('Getting column',col)
        column_count = self.lc.GetColumnCount()
        row_count = self.lc.GetItemCount()
        print('col count;',column_count)
        cols=[] #get cols
        for i in range(column_count):
            cols.append(self.GetCol(i,row_count))
        self.lc.ClearAll()
        self.DolcCols()
        s=numpy.flip(numpy.argsort(cols[col]))
        for ii,arg in enumerate(s):
            #print (arg,col1[arg],col2[arg],col3[arg],col4[arg])
            self.lc.InsertItem(ii,ii)
            for i in range(column_count):
                self.lc.SetItem(ii,i,str(cols[i][arg]))
                
        
        
    def Populate(self):
        self.orph={}

        self.lc.DeleteAllItems()
        for peak in self.parent._reference_peaks():
            n=peak.name
            self.orph[n]=1

        for grp,pks in self.parent._groups().items():
            for pk in pks:
                self.orph[pk]=0
        for pk,val in self.orph.items():
            num_items = self.lc.GetItemCount()
            ind=self.lc.InsertItem(num_items,pk)  #add assignment                    
            self.ParseFile(ind,pk)

            



    def ParseFile(self,ind,pk):

        f=os.path.join(self.parent._fuda_dir(), pk + '.out')
        if(os.path.exists(f)==False):
            return

        print('Reading ',f)
        inny=open(f)
        cnt=2
        tag=0
        for line in inny.readlines():
            if(len(line)>0 and line[0]=='#'):
                test=line.split()
                if(len(test)==4 and test[1]!='Parameter' and tag==1):
                    #print(cnt, test)
                    val='%.3f' % float(test[2])
                    self.lc.SetItem(ind,cnt,val) #add atom
                    cnt+=1
                if(len(test)==7 and test[2]=='Results'):
                    tag=1
                if(test[0][:2]=='##'):
                    tag=2


        f=os.path.join(self.parent._fuda_dir(), pk + '.dat')                    
        if(self.orph[pk]==0): #in a group
            tick=0
            for grp,vals in self.parent._groups().items():
                if(len(vals)>0):
                    if(pk==vals[0]):
                        tick=1
                        break
            if(tick==0):
                return

        if(os.path.exists(f)==False):
            return
        print('Reading ',f)

        yc=[]
        yd=[]
        
        inny=open(f)
        for line in inny.readlines():
            test=line.split()
            if(len(test)==4 and test[0]!='#'):
                yc.append(float(test[2]))
                yd.append(float(test[3]))

        yc=numpy.array(yc)
        yd=numpy.array(yd)

        chi2=numpy.sqrt(numpy.average((yc-yd)**2.))
        
        val='%.3f' % (chi2/numpy.max(yc)*100)
        self.lc.SetItem(ind,1,val) #add atom


    def Select(self,event):
        print()
        print ('index of selection:',self.lc.GetFocusedItem())
        self.SelGrpInd=self.lc.GetFocusedItem()
        self.SelGrp=self.lc.GetItem(self.lc.GetFocusedItem(), col=0).GetText()

        print('Selected group:',self.SelGrp)


        #is this in a group:
        sval=self.SelGrp
        for grp,vals in self.parent._groups().items():
            if(self.SelGrp in vals):
                sval='grp '+str(grp)

        ind=self.parent.PeakCombo.FindString(sval)
        self.parent.PeakCombo.SetSelection(ind)
        self.parent.peak_draw(True)
        

    def OnRefreshButton(self, event):
        rows= self.lc.GetItemCount()
        for i in range(rows): #redo parsing to set values.
            pk=self.lc.GetItem(i,col=0).GetText()
            self.ParseFile(i,pk)

        
        
    def OnCloseButton(self, event):
        self.Close()


        
