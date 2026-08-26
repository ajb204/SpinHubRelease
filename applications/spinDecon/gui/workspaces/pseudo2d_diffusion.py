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
#from attr import s
import wx
from spinDecon.domain.dimensions.viewer_contract import topology_for
from spinDecon.gui.context import context_for
from spinDecon.analysis.diffusion_service import DiffusionService
import string
import os
import numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm
# from matplotlib.widgets import Cursor
import matplotlib.patches as patches
from spinDecon.project.parameter_store import update_parameter_file
from scipy.optimize import curve_fit
# from .vpar_decon import vpar
############################################################################
# Frame for 1d slices
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()

class DiffusionROIFrame(wx.Frame):
    """Modeless table kept in sync with diffusion ROI selections."""

    COLUMNS = ("Trace", "ROI", "Selection", "Min", "Max", "D (fit)", "Fit error", "D (Gaussian)", "Gaussian error")

    def __init__(self, diffusion):
        super().__init__(diffusion, title="Diffusion ROIs", size=(900, 360))
        self.diffusion = diffusion
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(panel, label="Selected diffusion regions"), 0, wx.ALL, 8)
        self.roi_list = wx.ListCtrl(panel, style=wx.LC_REPORT)
        for i, label in enumerate(self.COLUMNS):
            self.roi_list.InsertColumn(i, label)
        widths = (55, 55, 80, 105, 105, 130, 110, 130, 125)
        for i, width in enumerate(widths):
            self.roi_list.SetColumnWidth(i, width)
        self.trace_images = wx.ImageList(42, 14)
        self.roi_list.AssignImageList(self.trace_images, wx.IMAGE_LIST_SMALL)
        sizer.Add(self.roi_list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        merge = wx.Button(panel, label="Merge selected")
        remove = wx.Button(panel, label="Remove selected")
        close = wx.Button(panel, label="Close")
        buttons.Add(merge, 0, wx.ALL, 8)
        buttons.Add(remove, 0, wx.ALL, 8)
        buttons.AddStretchSpacer()
        buttons.Add(close, 0, wx.ALL, 8)
        sizer.Add(buttons, 0, wx.EXPAND)
        panel.SetSizer(sizer)
        merge.Bind(wx.EVT_BUTTON, self._merge)
        remove.Bind(wx.EVT_BUTTON, self._remove)
        close.Bind(wx.EVT_BUTTON, lambda evt: self.Close())
        self.roi_list.Bind(wx.EVT_LIST_ITEM_SELECTED, self._select)
        self.roi_list.Bind(wx.EVT_LIST_ITEM_DESELECTED, self._select)
        self.roi_list.Bind(wx.EVT_LEFT_DCLICK, self._begin_cell_edit)
        self._cell_editor = None
        self.Bind(wx.EVT_CLOSE, self._close)
        self.refresh()

    def _selected(self):
        result = []
        item = -1
        while True:
            item = self.roi_list.GetNextItem(item, state=wx.LIST_STATE_SELECTED)
            if item == -1:
                return result
            result.append(item)

    @staticmethod
    def _number(value):
        return "—" if value is None or not numpy.isfinite(value) else "%.3e" % value

    def refresh(self):
        selected = set(self._selected())
        self.roi_list.DeleteAllItems()
        self.trace_images.RemoveAll()
        for i, stats in enumerate(self.diffusion.roi_stats):
            is_point = stats.get('kind') == 'point'
            row = [
                "", str(i + 1), "Point" if is_point else "Region",
                "%.5g" % stats['minimum'],
                "%.5g" % stats['maximum'],
                self._number(stats.get('gradient_D')),
                self._number(stats.get('gradient_error')),
                "—" if is_point else self._number(stats.get('gaussian_D')),
                "—" if is_point else self._number(stats.get('gaussian_sigma')),
            ]
            # A small coloured line is easier to associate with the plotted
            # trace than another textual label in the diffusion panel.
            rgba = matplotlib.colors.to_rgba('C%d' % (i % 10))
            colour = wx.Colour(*[int(255 * c) for c in rgba[:3]])
            bitmap = wx.Bitmap(42, 14)
            dc = wx.MemoryDC(bitmap)
            dc.SetBackground(wx.Brush(wx.Colour(255, 255, 255)))
            dc.Clear()
            dc.SetPen(wx.Pen(colour, 3))
            dc.DrawLine(5, 7, 37, 7)
            dc.SelectObject(wx.NullBitmap)
            image_index = self.trace_images.Add(bitmap)
            idx = self.roi_list.InsertItem(self.roi_list.GetItemCount(), row[0], image_index)
            for col, value in enumerate(row[1:], 1):
                self.roi_list.SetItem(idx, col, value)
        self.select_indices(selected)

    def select_indices(self, indices):
        wanted = set(indices)
        for i in range(self.roi_list.GetItemCount()):
            self.roi_list.SetItemState(i, wx.LIST_STATE_SELECTED if i in wanted else 0, wx.LIST_STATE_SELECTED)


    def _begin_cell_edit(self, event):
        """Edit an ROI Min/Max cell in place; Enter commits and reanalyses it."""
        pos = event.GetPosition()
        try:
            item, _flags, col = self.roi_list.HitTestSubItem(pos)
        except (AttributeError, ValueError):
            item, _flags = self.roi_list.HitTest(pos)
            col = -1
            if item != wx.NOT_FOUND:
                x = pos.x
                left = 0
                for candidate in range(self.roi_list.GetColumnCount()):
                    right = left + self.roi_list.GetColumnWidth(candidate)
                    if left <= x < right:
                        col = candidate
                        break
                    left = right
        if item == wx.NOT_FOUND or col not in (3, 4):
            event.Skip()
            return
        self._finish_cell_edit(commit=False)
        rect = self.roi_list.GetSubItemRect(item, col)
        editor = wx.TextCtrl(self.roi_list, value=self.roi_list.GetItemText(item, col),
                             pos=rect.GetPosition(), size=rect.GetSize(),
                             style=wx.TE_PROCESS_ENTER)
        editor._roi_item = item
        editor._roi_col = col
        editor.Bind(wx.EVT_TEXT_ENTER, self._commit_cell_edit)
        editor.Bind(wx.EVT_KILL_FOCUS, self._cancel_cell_edit)
        self._cell_editor = editor
        editor.SetFocus()
        editor.SelectAll()

    def _commit_cell_edit(self, event):
        editor = self._cell_editor
        if editor is None:
            return
        try:
            value = float(editor.GetValue())
        except (TypeError, ValueError):
            wx.Bell()
            editor.SetFocus()
            editor.SelectAll()
            return
        item, col = editor._roi_item, editor._roi_col
        a, b = self.diffusion.roi_ranges[item]
        if col == 3:
            a = value
        else:
            b = value
        self._finish_cell_edit(commit=False)
        self.diffusion.update_roi_range(item, a, b)

    def _cancel_cell_edit(self, event):
        # Losing focus cancels an unfinished edit; Enter is the explicit commit.
        wx.CallAfter(self._finish_cell_edit, False)
        event.Skip()

    def _finish_cell_edit(self, commit=False):
        editor = self._cell_editor
        self._cell_editor = None
        if editor is not None:
            editor.Destroy()

    def _select(self, event):
        wx.CallAfter(self.diffusion.highlight_rois, self._selected())

    def _remove(self, event):
        self.diffusion.remove_rois(self._selected())

    def _merge(self, event):
        self.diffusion.merge_rois(self._selected())

    def _close(self, event):
        self.diffusion.roi_frame = None
        event.Skip()


class Pseudo2DDiffusion(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.app_context = context_for(parent, tabOne)
        self.diffusion_service = (getattr(self.app_context, 'diffusion', None)
                                  if self.app_context is not None else None)
        if self.diffusion_service is None:
            self.diffusion_service = DiffusionService(tabOne)
        self.topology = topology_for(tabOne)
        self.spectral_dim_count = self.topology.spectral_dim_count
        self.physical_dim_count = self.topology.physical_dim_count
        self.dim = self.spectral_dim_count  # compatibility alias: spectral only
        self.state = (getattr(self.app_context, 'project', None)
                      if self.app_context is not None else getattr(tabOne, 'state', getattr(parent, 'state', None)))
        self.sum=(0.,2.)
        self.peak=self.diffusion_service.peaks
        self.thresh=self.diffusion_service.data_maximum()
        self.offset=0


        
        dmin, dmax = self.diffusion_service.spectral_bounds

        self.create_main_panel()

        self.set_default_values() #upack Grp save
        self.draw_figure()
        self._restore_saved_rois()

    def plot_scatters(self):
        """Plot ROI mean I/I0 traces and their within-ROI standard deviations."""
        for artist in self.scatters:
            try:
                artist.remove()
            except (AttributeError, ValueError):
                pass
        self.axes_grad.cla()
        self.scatters = []

        k = (numpy.asarray(self.gzlvl1, dtype=float)
             if len(self.gzlvl1) == self.diffusion_service.data.shape[0]
             else -numpy.arange(self.diffusion_service.data.shape[0], dtype=float))
        # Internally k is negative in the Stejskal-Tanner fit.  Display -k so
        # attenuation progresses from zero towards positive diffusion weighting.
        display_x = -k
        order = numpy.argsort(display_x)

        for ii, line in enumerate(self.scatter_data):
            mean_trace = numpy.asarray(line, dtype=float)
            fit = numpy.asarray(self.scatter_data_norm[ii], dtype=float)
            errors = numpy.asarray(self.scatter_data_err[ii], dtype=float)
            colour = 'C%d' % (ii % 10)

            valid = numpy.isfinite(display_x) & numpy.isfinite(mean_trace)
            if errors.size == mean_trace.size:
                valid_err = valid & numpy.isfinite(errors)
                idx = order[valid_err[order]]
                artist = self.axes_grad.errorbar(display_x[idx], mean_trace[idx],
                                                 yerr=errors[idx], color=colour, fmt='o')
                self.scatters.append(artist)
            else:
                idx = order[valid[order]]
                self.scatters.append(self.axes_grad.scatter(display_x[idx], mean_trace[idx],
                                                             color=colour, marker='o'))

            valid_curve = numpy.isfinite(display_x) & numpy.isfinite(fit)
            idx = order[valid_curve[order]]
            self.axes_grad.plot(display_x[idx], fit[idx], ls='--', color=colour)

        self.axes_grad.set_ylabel(r'$I/I_0$')
        self.axes_grad.set_xlabel('Diffusion weighting')
        self.canvas.draw_idle()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.axes = self.fig.add_subplot(221)
        self.axes_d = self.fig.add_subplot(223)
        self.axes_dh = self.axes_d.twinx()

        self.axes_h = self.axes.twinx()
        self.axes_proj = self.fig.add_subplot(243)

        self.axes_grad = self.fig.add_subplot(244)

        self.axes_sca = self.fig.add_subplot(247)
        
        self.axes_err = self.fig.add_subplot(248)

        # Use more of the frame for data while retaining enough room for labels.
        self.fig.subplots_adjust(left=0.065, right=0.985, bottom=0.075, top=0.965,
                                 wspace=0.30, hspace=0.28)
        
        self.cursor_shown = False
        self.number_scatters=0
        self.pressed = False
        self.moved = False
        self.rectangles = []
        self.scatters = []
        self.scatter_data = []
        self.scatter_data_norm = []
        self.scatter_data_err = []
        self.scatter_data_val=[]
        self.verticals = []
        self.roi_ranges = []
        self.roi_stats = []
        self.roi_frame = None
        self._highlighted_rois = set()
        self.not_yet_drawn = True
        self.gzlvl1 = []
        # self.drawbutton = wx.Button(self, -1, "Draw!")
        # self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        # self.canvas.mpl_connect('scroll_event', self.on_scroll)

        # Use the project toolbar.  Its Home/Subplots tools are intentionally
        # hidden; for diffusion the project Redraw tool is labelled Draw!.
        class DiffusionToolbar(RedrawNavigationToolbar):
            _redraw_tool = ("Draw!", "Recalculate and redraw diffusion analysis", "redraw_pencil", "redraw_view")

        self.toolbar = DiffusionToolbar(self.canvas, self.on_draw_button, coordinates=False)
        self._build_diffusion_toolbar_controls()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(6)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def _build_diffusion_toolbar_controls(self):
        """Append diffusion controls immediately after Matplotlib's Save tool."""
        tb = self.toolbar
        tb.AddSeparator()
        self.savebutton = wx.Button(tb, label="Save", size=(-1, 26))
        self.savebutton.SetToolTip("Save diffusion parameters")
        self.savebutton.Bind(wx.EVT_BUTTON, self.on_save_button)
        tb.AddControl(self.savebutton)

        self.noiseFacLab = wx.StaticText(tb, label=" NoiseFac: ")
        self.noiseFac = wx.TextCtrl(tb, size=(50, -1), style=wx.TE_PROCESS_ENTER)
        self.noiseFac.SetToolTip("Noise threshold multiplier; press Enter to recalculate diffusion analysis")
        self.noiseFac.Bind(wx.EVT_TEXT_ENTER, self.on_noise_fac_enter)
        for control in (self.noiseFacLab, self.noiseFac):
            tb.AddControl(control)

        self.roi_button = wx.Button(tb, label="ROI", size=(-1, 26))
        self.roi_button.SetToolTip("View and manage selected diffusion regions")
        self.roi_button.Bind(wx.EVT_BUTTON, self.on_roi_button)
        tb.AddControl(self.roi_button)
        tb.Realize()

    def on_noise_fac_enter(self, event=None):
        """Recalculate the complete analysis and all existing ROIs.

        NoiseFac changes the acceptance threshold, so every per-ppm fit and
        every ROI statistic/error/histogram must be rebuilt from the newly
        accepted population.  Preserve only the ROI ppm ranges themselves.
        """
        try:
            value = float(self.noiseFac.GetValue())
        except (TypeError, ValueError):
            wx.Bell()
            return
        if not numpy.isfinite(value) or value < 0:
            wx.Bell()
            return

        ranges = list(self.roi_ranges)
        highlighted = set(self._highlighted_rois)

        # Clear all derived plot artists before rerunning AnalDiff.  Do not use
        # on_draw_button here because that deliberately deletes the ROI ranges.
        for axis in (self.axes_grad, self.axes_proj, self.axes_sca, self.axes_err,
                     self.axes_h, self.axes_d, self.axes_dh):
            axis.cla()
        self.verticals = []
        self.rectangles = []
        self.scatters = []
        self.scatter_data = []
        self.scatter_data_norm = []
        self.scatter_data_err = []
        self.scatter_data_val = []
        self.roi_stats = []

        self.draw_figure()

        # Recalculate each ROI from the new accepted I/I0 population.
        self.roi_ranges = []
        for a, b in ranges:
            self._append_roi(a, b)
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = {i for i in highlighted if i < len(self.roi_ranges)}
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()

    def on_roi_button(self, event=None):
        if self.roi_frame is None:
            self.roi_frame = DiffusionROIFrame(self)
        self.roi_frame.refresh()
        self.roi_frame.Show()
        self.roi_frame.Raise()

    def _refresh_roi_manager(self):
        if self.roi_frame is not None:
            self.roi_frame.refresh()

    def _rebuild_roi_overlays(self):
        for artist in self.verticals + self.rectangles:
            try:
                artist.remove()
            except (ValueError, AttributeError):
                pass
        self.verticals = []
        self.rectangles = []
        for i, (a, b) in enumerate(self.roi_ranges):
            colour = 'C' + str(i)
            self.verticals.append(self.axes_h.axvline(a, color=colour, linewidth=2, ls='--'))
            if a != b:
                self.verticals.append(self.axes_h.axvline(b, color=colour, linewidth=2, ls='--'))
                left, right = max(a, b), min(a, b)
                rgba = list(matplotlib.colors.to_rgba(colour)[:3]) + [0.3]
                rect = patches.Rectangle((left, -0.5), right-left, 4, linewidth=0, facecolor=rgba)
                self.rectangles.append(rect)
                self.axes.add_patch(rect)
        self.highlight_rois(self._highlighted_rois, redraw=False)

    def highlight_rois(self, indices, redraw=True):
        self._highlighted_rois = set(indices)
        rect_i = 0
        for i, (a, b) in enumerate(self.roi_ranges):
            if a != b and rect_i < len(self.rectangles):
                selected = i in self._highlighted_rois
                self.rectangles[rect_i].set_alpha(0.65 if selected else 0.3)
                self.rectangles[rect_i].set_linewidth(2 if selected else 0)
                self.rectangles[rect_i].set_edgecolor('black' if selected else 'none')
                rect_i += 1
        if self.roi_frame is not None:
            self.roi_frame.select_indices(self._highlighted_rois)
        if redraw:
            self.canvas.draw_idle()

    def remove_rois(self, indices):
        remove = set(indices)
        if not remove:
            return
        for attr in ('scatter_data', 'scatter_data_norm', 'scatter_data_err', 'scatter_data_val', 'roi_ranges', 'roi_stats'):
            values = getattr(self, attr)
            setattr(self, attr, [value for i, value in enumerate(values) if i not in remove])
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = set()
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()

    def _roi_indices(self, a, b):
        c1 = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - a)))
        c2 = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - b)))
        lo, hi = sorted((c1, c2))
        return lo, hi + 1

    def _gradient_fit(self, ydat):
        """Fit the top-right attenuation plot: log(I/I0) = D*k + intercept.

        Return the diffusion coefficient (the fitted gradient), its 1-sigma
        standard error, and the fitted attenuation curve.  Invalid/non-positive
        intensities are excluded exactly where the logarithm is taken; this
        avoids the NaNs introduced by logging the complete ROI vector.
        """
        ydat = numpy.asarray(ydat, dtype=float)
        k = numpy.asarray(self.gzlvl1, dtype=float)
        if len(k) != len(ydat):
            return None, None, numpy.asarray(ydat)
        mask = numpy.isfinite(ydat) & numpy.isfinite(k) & (ydat > 0)
        # The slope itself only needs two distinct attenuation points.  Do not
        # make the D value depend on numpy.polyfit(cov=True), which requires
        # additional degrees of freedom and can fail even when the line is
        # perfectly fit-able.
        if mask.sum() < 2:
            return None, None, numpy.full_like(ydat, numpy.nan, dtype=float)
        xfit = k[mask]
        logy = numpy.log(ydat[mask])
        if numpy.unique(xfit).size < 2:
            return None, None, numpy.full_like(ydat, numpy.nan, dtype=float)
        try:
            design = numpy.column_stack((xfit, numpy.ones_like(xfit)))
            coeff, _, _, _ = numpy.linalg.lstsq(design, logy, rcond=None)
            dval, intercept = float(coeff[0]), float(coeff[1])
        except (ValueError, numpy.linalg.LinAlgError):
            return None, None, numpy.full_like(ydat, numpy.nan, dtype=float)

        # 1-sigma standard error of the fitted gradient.  With only two points
        # the slope is defined but there are no residual degrees of freedom, so
        # retain the D value and report no error rather than turning D into NaN.
        derr = None
        if len(xfit) > 2:
            residual = logy - (dval * xfit + intercept)
            sxx = numpy.sum((xfit - numpy.mean(xfit)) ** 2)
            if sxx > 0:
                variance = numpy.sum(residual ** 2) / float(len(xfit) - 2)
                candidate = numpy.sqrt(max(0.0, variance / sxx))
                if numpy.isfinite(candidate):
                    derr = float(candidate)

        curve = numpy.exp(intercept + dval * k)
        return dval, derr, curve

    @staticmethod
    def _histogram_for_values(values):
        """Return an adaptive, equal-width histogram suitable for ROI D values."""
        vals = numpy.asarray(values, dtype=float)
        vals = vals[numpy.isfinite(vals)]
        if vals.size < 3:
            return None

        spread = numpy.ptp(vals)
        if not numpy.isfinite(spread) or spread == 0:
            pad = max(abs(float(vals[0])) * 0.05, numpy.finfo(float).eps)
            edges = numpy.linspace(vals[0] - pad, vals[0] + pad, 4)
        else:
            # Freedman-Diaconis adapts to both ROI size and distribution width.
            edges = numpy.histogram_bin_edges(vals, bins='fd')
            nbins = int(numpy.clip(len(edges) - 1, 3, 50))
            edges = numpy.linspace(float(vals.min()), float(vals.max()), nbins + 1)
        hist, edges = numpy.histogram(vals, bins=edges)
        centres = 0.5 * (edges[:-1] + edges[1:])
        return vals, hist.astype(float), edges, centres

    @staticmethod
    def _gaussian(x, amplitude, mean, sigma):
        return amplitude * numpy.exp(-0.5 * ((x - mean) / sigma) ** 2)

    def _gaussian_fit_for_range(self, a, b):
        """Fit a Gaussian to the adaptive histogram of pixel-wise ROI D values."""
        lo, hi = self._roi_indices(a, b)
        histogram = self._histogram_for_values(self.dsv[lo:hi])
        if histogram is None:
            return None
        vals, hist, edges, centres = histogram
        if not numpy.any(hist):
            return None

        sigma0 = float(numpy.std(vals, ddof=1))
        if not numpy.isfinite(sigma0) or sigma0 <= 0:
            sigma0 = max(float(numpy.ptp(vals)) / 6.0, numpy.finfo(float).eps)
        p0 = (float(hist.max()), float(numpy.mean(vals)), sigma0)
        try:
            params, _ = curve_fit(
                self._gaussian, centres, hist, p0=p0,
                bounds=([0.0, float(vals.min()), numpy.finfo(float).eps],
                        [numpy.inf, float(vals.max()), numpy.inf]),
                maxfev=10000,
            )
            amplitude, mean, sigma = map(float, params)
            sigma = abs(sigma)
        except (RuntimeError, ValueError, FloatingPointError):
            # A sparse ROI may not support a stable nonlinear histogram fit.
            # Its sample moments are a useful, finite fallback.
            amplitude = float(hist.max())
            mean = float(numpy.mean(vals))
            sigma = sigma0

        gx = numpy.linspace(float(edges[0]), float(edges[-1]), 300)
        gy = self._gaussian(gx, amplitude, mean, sigma)
        return {
            'mean': mean, 'sigma': sigma, 'hist': hist, 'edges': edges,
            'centres': centres, 'xfit': gx, 'yfit': gy,
        }

    def _build_roi_data(self, a, b):
        """Build an ROI from accepted, individually normalised I/I0 traces.

        The top-right trace is the mean of the same normalised spectral points
        shown in the bottom-middle diagnostic; its error bars are their standard
        deviation at each diffusion weighting.
        """
        lo, hi = self._roi_indices(a, b)
        is_point = lo + 1 == hi
        region = numpy.asarray(self.normalisedDiffFull[:, lo:hi], dtype=float)
        accepted = numpy.isfinite(self.dsv[lo:hi])
        region = region[:, accepted]

        if region.shape[1]:
            ydat = numpy.nanmean(region, axis=1)
            yerr = (numpy.nanstd(region, axis=1) if region.shape[1] > 1
                    else numpy.zeros(region.shape[0], dtype=float))
        else:
            ydat = numpy.full(self.diffusion_service.data.shape[0], numpy.nan, dtype=float)
            yerr = numpy.full(self.diffusion_service.data.shape[0], numpy.nan, dtype=float)

        dgrad, dgrad_err, yfit = self._gradient_fit(ydat)
        gmean = gsigma = None
        if not is_point:
            gaussian = self._gaussian_fit_for_range(a, b)
            if gaussian is not None:
                gmean, gsigma = gaussian['mean'], gaussian['sigma']
        stats = dict(kind='point' if is_point else 'region', minimum=min(a,b), maximum=max(a,b),
                     gradient_D=dgrad, gradient_error=dgrad_err, gaussian_D=gmean, gaussian_sigma=gsigma)
        return numpy.asarray(ydat), numpy.asarray(yfit), numpy.asarray(yerr), stats

    def _append_roi(self, a, b):
        ydat, yfit, yerr, stats = self._build_roi_data(a, b)
        self.scatter_data.append(ydat)
        self.scatter_data_norm.append(yfit)
        self.scatter_data_err.append(yerr)
        self.scatter_data_val.append(stats['gradient_D'] if stats['gradient_D'] is not None else numpy.nan)
        self.roi_ranges.append((a, b))
        self.roi_stats.append(stats)

    def _plot_roi_histograms(self):
        self.axes_err.cla()
        for i, ((a, b), stats) in enumerate(zip(self.roi_ranges, self.roi_stats)):
            if stats.get('kind') != 'region':
                continue
            gaussian = self._gaussian_fit_for_range(a, b)
            if gaussian is None:
                continue
            edges = gaussian['edges']
            self.axes_err.bar(
                gaussian['centres'], gaussian['hist'], width=numpy.diff(edges),
                align='center', alpha=0.35, color='C'+str(i),
                edgecolor='C'+str(i), linewidth=0.8,
            )
            self.axes_err.plot(gaussian['xfit'], gaussian['yfit'],
                               color='C'+str(i), linewidth=1.5)
        self.axes_err.set_title('ROI diffusion histograms')
        self.axes_err.set_xlabel(r'$D$ / cm$^2$ s$^{-1}$')
        self.axes_err.set_ylabel('Count')

    def update_roi_range(self, index, a, b):
        """Replace one ROI range and recompute every result derived from it."""
        if not (0 <= index < len(self.roi_ranges)):
            return
        try:
            a = float(a)
            b = float(b)
        except (TypeError, ValueError):
            return
        if not (numpy.isfinite(a) and numpy.isfinite(b)):
            return

        # Snap edited ppm values to actual spectral coordinates, matching mouse ROIs.
        ia = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - a)))
        ib = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - b)))
        a = float(self.diffusion_service.spectral_axis[ia])
        b = float(self.diffusion_service.spectral_axis[ib])
        ydat, yfit, yerr, stats = self._build_roi_data(a, b)
        self.scatter_data[index] = ydat
        self.scatter_data_norm[index] = yfit
        self.scatter_data_err[index] = yerr
        self.scatter_data_val[index] = (stats['gradient_D']
                                        if stats['gradient_D'] is not None else numpy.nan)
        self.roi_ranges[index] = (a, b)
        self.roi_stats[index] = stats
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = {index}
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()

    def merge_rois(self, indices):
        indices = sorted(set(indices))
        if len(indices) < 2:
            return
        a = min(min(self.roi_ranges[i]) for i in indices)
        b = max(max(self.roi_ranges[i]) for i in indices)
        # Remove source ROIs, then create one ROI from the complete merged span.
        remove = set(indices)
        for attr in ('scatter_data', 'scatter_data_norm', 'scatter_data_err', 'scatter_data_val', 'roi_ranges', 'roi_stats'):
            values = getattr(self, attr)
            setattr(self, attr, [value for i, value in enumerate(values) if i not in remove])
        self._append_roi(a, b)
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = {len(self.roi_ranges)-1}
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()

    def AnalDiff(self,xs,ys,zs):
        print ('Analysing diffusion coefficients')
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
            print(gzs)
            for gz in gzs:
                try:
                    if float(gz)>2:
                        G = Gmax*float(gz)/30000.

                        # print(delta)
                        gzlvl1.append(-delta**2*G**2*gamma**2*(BigT-(1/3)*(delta)))
                except:
                    pass
        gzlvl1=numpy.array(gzlvl1)


        
        #print (xs)
        #print (ys.shape)
        #print (XX)
        #print (YY)
        #print (zs.shape)
        #print (XX.shape)
        #import sys
        #sys.exit(100)

        #self.kdiff=(copy.deepcopy(ys)*0.002*delta)**2 *(Delta-delta/3)

        # Diffusion signal threshold: use the noise estimate already calculated
        # for the main NMR tab, scaled by the diffusion-specific NoiseFac.
        # This keeps the diffusion analysis consistent with the spectrum noise
        # model and removes the need for a separate ppm noise region here.
        try:
            noise = abs(float(self.diffusion_service.noise_value()))
            nfac = float(self.noiseFac.GetValue())
            if not numpy.isfinite(noise) or not numpy.isfinite(nfac) or nfac < 0:
                raise ValueError("invalid diffusion noise threshold")
            thresh = noise * nfac
        except (AttributeError, TypeError, ValueError):
            thresh = 0.0

        self.gzlvl1 = numpy.asarray(gzlvl1, dtype=float)
        self.ycalcDiffFull = numpy.full_like(zs, numpy.nan, dtype=float)
        self.normalisedDiffFull = numpy.full_like(zs, numpy.nan, dtype=float)
        self.normalisedFitFull = numpy.full_like(zs, numpy.nan, dtype=float)
        self.dsv = numpy.full(len(xs), numpy.nan, dtype=float)
        self.dsv_error = numpy.full(len(xs), numpy.nan, dtype=float)
        self.asv = numpy.full(len(xs), numpy.nan, dtype=float)
        self.extrapolated_i0 = numpy.full(len(xs), numpy.nan, dtype=float)

        if len(self.gzlvl1) != zs.shape[0]:
            print('Diffusion gradient axis does not match pseudo2D data length')
            return

        for i in range(len(xs)):
            trace = numpy.asarray(zs[:, i], dtype=float)
            fit_trace = trace.copy()
            fit_trace[~(numpy.isfinite(trace) & (trace > thresh))] = numpy.nan
            dval, derr, curve = self._gradient_fit(fit_trace)
            # With the internal (negative) diffusion-weighting coordinate, a
            # physically valid attenuation must have a strictly positive
            # gradient.  Zero/negative D corresponds to a flat or increasing
            # trace when displayed against -k, so reject that ppm completely.
            # Leaving every derived array as NaN propagates the rejection to
            # the middle diagnostics, ROI populations and ROI histograms.
            if dval is None or not numpy.isfinite(dval) or dval <= 0.0:
                continue
            self.dsv[i] = dval
            if derr is not None and numpy.isfinite(derr):
                self.dsv_error[i] = derr
            self.ycalcDiffFull[:, i] = curve
            valid = numpy.isfinite(trace) & numpy.isfinite(curve) & (curve > 0)
            if numpy.any(valid):
                amplitude = numpy.nanmedian(trace[valid] / numpy.exp(dval * self.gzlvl1[valid]))
                self.asv[i] = amplitude
                # Express the extrapolated zero-gradient intensity relative to
                # the least diffusion-weighted measured point.  Fits requiring
                # an extrapolation above 2x the measured intensity are treated
                # as unreliable and excluded from maps, histograms and ROIs.
                ref = int(numpy.nanargmin(numpy.abs(self.gzlvl1)))
                reference = trace[ref]
                if numpy.isfinite(reference) and reference > 0:
                    self.extrapolated_i0[i] = amplitude / reference
                if (not numpy.isfinite(self.extrapolated_i0[i]) or
                        self.extrapolated_i0[i] > 2.0):
                    self.dsv[i] = numpy.nan
                    self.dsv_error[i] = numpy.nan
                    self.asv[i] = numpy.nan
                    self.ycalcDiffFull[:, i] = numpy.nan
                elif numpy.isfinite(amplitude) and amplitude > 0:
                    # Persist the exact per-point I/I0 population used later
                    # for ROI means and error bars.
                    self.normalisedDiffFull[:, i] = trace / amplitude
                    self.normalisedFitFull[:, i] = curve / amplitude

        print('done')
        
    def draw_figure(self):
        #try:
        self.draw_figureGO()
        #except:
        #    pass

    def draw_figureGO(self):
        """ Redraws the figure
        """
        self.axes.clear()
        self.axes_h.clear()
        print('drawing')
        #sele1=self.ComboBox1.GetSelection()
        self.thresh=float(self.diffusion_service.threshold())


        xs=self.diffusion_service.spectral_axis
        ys=numpy.arange(self.diffusion_service.data.shape[0])
        zs=numpy.array(self.diffusion_service.data)

        self.ds=self.AnalDiff(xs,ys,zs)
        
        y2s=numpy.zeros_like(ys)
        y2s.fill(self.thresh)
        levels = [self.thresh]
        for x in range(12):
            levels.append(levels[-1]*1.4)

        cmap = plt.get_cmap('Oranges')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        self.axes.set_xlabel(self.diffusion_service.spectral_label,fontsize=8)
        self.axes.pcolormesh(xs,ys,zs,label='data', norm=norm, cmap = cmap) #plot contour
        self.canvas.draw()
        self.axes.set_ylim(ys[0], ys[-1])
        self.axes.set_xlim(xs[0], xs[-1])

        self.axes_d.set_yscale("log")         
        self.axes_d.scatter(xs,self.dsv)
        #self.axes_d.set_ylim(ys[0], ys[-1])
        self.axes_d.set_xlim(xs[0], xs[-1])

        self.axes_dh.plot(xs,self.asv)
        #self.axes_d.set_ylim(ys[0], ys[-1])
        self.axes_dh.set_xlim(xs[0], xs[-1])
        
        print (self.ycalcDiffFull.shape)
        
        valid_diffusion = numpy.isfinite(self.dsv)
        for x in range(self.diffusion_service.data.shape[0]):
            measured = numpy.where(valid_diffusion, self.diffusion_service.data[x, :], numpy.nan)
            fitted = numpy.where(valid_diffusion, self.ycalcDiffFull[x, :], numpy.nan)
            self.axes_proj.plot(xs, measured, lw=0.5, label=str(x), color='k')
            self.axes_proj.plot(xs, fitted, lw=0.5, label=str(x), color='r')

            # Restore the normalised diagnostic.  These are exactly the I/I0
            # values subsequently averaged within an ROI; their spread becomes
            # the top-right error bars.
            self.axes_sca.plot(xs, self.normalisedDiffFull[x, :], lw=0.5,
                               label=str(x), color='k')
            self.axes_sca.plot(xs, self.normalisedFitFull[x, :], lw=0.5,
                               label=str(x), color='r')
        self.axes_sca.set_ylabel(r'$I/I_0$')

            

            

            # self.axes_proj.legend()
        self.axes_proj.set_xlim(xs[0], xs[-1])
        self.axes_sca.set_xlim(xs[0], xs[-1])

        # self.combinedTransform = self.axes_h.transData + self.axes.transData.inverted()
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.fig.canvas.mpl_connect('draw_event', self._on_canvas_draw)
        self.fig.canvas.mpl_connect('button_press_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    def on_key(self, event):
        print(event.key)
        if event.key=='h':
            self.h_line.set_visible(True)
        if event.key =='c' or event.key=='ctrl+c':
            for x in self.verticals:
                x.remove()
            for x in self.scatters:
                x.remove()
            for x in self.rectangles:
                x.remove()
            self.verticals=[]
            self.scatters = []
            self.scatter_data = []
            self.scatter_data_norm = []
            self.scatter_data_err = []
            self.scatter_data_val = []
            
            self.rectangles=[]
            self.roi_ranges = []
            self.roi_stats = []
            self._highlighted_rois = set()
            self.number_scatters = 0
            self._refresh_roi_manager()
            #self.axes_proj.cla()
            self.axes_grad.cla()
            self.axes_err.cla()
            
            self.canvas.draw()

        # if event.key=='v':
            # if not self.cursor_shown:
            #     print(self.cursor_shown)
                # self.cursor = Cursor(self.axes, horizOn = False, vertOn=True, color='r', linewidth=2, useblit=True)
            # else:
            #     self.cursor.clear()
            # self.canvas.draw()

    def on_scroll(self, event):
        # print('scrolling')
        self.ymin,self.ymax=self.axes_h.get_ylim()
        print(self.ymin, self.ymax)
        self.axes_h.set_ylim(self.ymin+(self.ymin*0.05*event.step), self.ymax+(self.ymax*0.05*event.step))
        self.axes_h.draw_artist(self.h_line)
        # self.axes_h.draw()

    def _on_canvas_draw(self, event=None):
        """Refresh the blit background after any full canvas redraw.

        The moving 1-D trace belongs to ``axes_h`` (the twinned top-left
        axes).  Its background must therefore be captured and blitted using
        that axes' current bbox.  This is especially important after changing
        subplot margins, resizing the frame, zooming, or pressing Raw.
        """
        if hasattr(self, 'fig'):
            # Keep a full-canvas background for the animated spectrum.  On
            # some GUI backends (notably macOS with HiDPI scaling), blitting a
            # twinned-axis bbox can leave a strip at the edge of the axes
            # unchanged after subplot geometry/margins have changed.  A figure
            # background uses the canvas' native pixel extent and avoids that
            # bbox conversion mismatch.
            self.background = self.canvas.copy_from_bbox(self.fig.bbox)

    def on_mouse_move(self,event):
        if event.inaxes==self.axes_h:
            if self.not_yet_drawn == True:
                self.current_h = 0
                self.v_line = self.axes_h.axvline(self.diffusion_service.spectral_axis[0], color='r',
                                                  linewidth=2, animated=True)
                self.h_line, = self.axes_h.plot(
                    self.diffusion_service.spectral_axis, numpy.zeros_like(self.diffusion_service.spectral_axis),
                    color='k', linewidth=0.5, animated=True)

                self.axes_h.set_ylim(numpy.min(self.diffusion_service.data), numpy.max(self.diffusion_service.data))

                # Animated artists are omitted from the normal draw.  Capture
                # the complete, newly laid-out top-left axes only afterwards,
                # so the moving trace can be restored across its full width.
                self.canvas.draw()
                self.background = self.canvas.copy_from_bbox(self.fig.bbox)
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
                # print(x,y)
            # if event.inaxes == self.axes:
            #     print(event.ydata)

                # new_dataPoint = int(event.ydata) #(int(numpy.floor(self.combinedTransform.transform(pt_data2)[1])))
                self.canvas.restore_region(self.background)
                self.v_line.set_xdata([event.xdata, event.xdata])

                # Always update both coordinates from the complete spectral
                # dimension used by the contour.  This prevents an animated
                # line from retaining a stale/cropped x array after redraws.
                trace_x = numpy.asarray(self.diffusion_service.spectral_axis)
                trace_y = numpy.asarray(self.diffusion_service.data[new_dataPoint, :])
                npts = min(trace_x.size, trace_y.size)
                self.h_line.set_data(trace_x[:npts], trace_y[:npts])

                self.axes_h.draw_artist(self.h_line)
                self.axes_h.draw_artist(self.v_line)
                self.current_h = new_dataPoint

                # Blit the same full-canvas region used for the saved
                # background.  This is intentionally broader than axes_h.bbox:
                # it fixes partial updates seen with twinned axes on macOS.
                self.canvas.blit(self.fig.bbox)

                # if new_dataPoint != self.current_h:


    def on_draw_button(self, event=None):
        """Reset the plot views and recalculate without deleting ROIs.

        The toolbar redraw/raw action used to emulate the ``c`` key before
        redrawing.  That key handler is the explicit *clear ROI* operation, so
        pressing the toolbar button unintentionally discarded every selected
        region.  Preserve the ROI geometry here, clear only derived artists and
        axes (which also resets zoom/pan limits), rerun the diffusion analysis,
        and finally rebuild each ROI from the new per-ppm results.
        """
        ranges = list(self.roi_ranges)
        highlighted = set(self._highlighted_rois)

        # Clearing the axes resets all zoom/pan views to the limits established
        # by draw_figure().  Do not call on_key('c'): that deliberately deletes
        # ROI geometry.
        for axis in (self.axes_grad, self.axes_proj, self.axes_sca, self.axes_err,
                     self.axes_h, self.axes_d, self.axes_dh):
            axis.cla()

        # axes_h.cla() removes the animated mouse-trace artists.  Recreate
        # them on the next mouse movement and capture a fresh full-width
        # blitting background for the reset layout.
        self.not_yet_drawn = True
        self.background = None
        self.h_line = None
        self.v_line = None

        self.verticals = []
        self.rectangles = []
        self.scatters = []
        self.scatter_data = []
        self.scatter_data_norm = []
        self.scatter_data_err = []
        self.scatter_data_val = []
        self.roi_stats = []

        # Recalculate the complete per-ppm diffusion analysis using the current
        # data, NoiseFac and all current frame settings.
        self.draw_figure()

        # Reanalyse, rather than merely redraw, the regions that existed when
        # the button was pressed.
        self.roi_ranges = []
        for a, b in ranges:
            self._append_roi(a, b)
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = {i for i in highlighted
                                  if i < len(self.roi_ranges)}
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()

    def on_pick(self, event):
        if event.inaxes == self.axes_h and getattr(event, 'dblclick', False) and event.xdata is not None:
            # Double-click an existing plotted ROI to select it in the manager
            # without stealing the normal click/drag gesture used to create ROIs.
            hits = []
            span = abs(self.diffusion_service.spectral_axis[-1] - self.diffusion_service.spectral_axis[0]) if len(self.diffusion_service.spectral_axis) > 1 else 1.0
            tolerance = span * 0.005
            for i, (a, b) in enumerate(self.roi_ranges):
                lo, hi = sorted((a, b))
                if (lo - tolerance) <= event.xdata <= (hi + tolerance):
                    hits.append(i)
            if hits:
                self.highlight_rois(hits)
                if self.roi_frame is None:
                    self.on_roi_button()
                return
        if event.inaxes==self.axes_h:
            self.pressed = True
            self.origin = event.xdata
            print(self.diffusion_service.pseudo_axis)

    def on_release(self, event):
        if event.inaxes != self.axes_h or event.xdata is None:
            self.pressed = False
            self.moved = False
            return

        coord = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - event.xdata)))
        end = float(self.diffusion_service.spectral_axis[coord])
        start = end
        if self.moved and self.origin is not None:
            coord2 = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - self.origin)))
            start = float(self.diffusion_service.spectral_axis[coord2])

        self._append_roi(start, end)
        self.number_scatters = len(self.scatter_data)
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()

        self.pressed = False
        self.moved = False
        if hasattr(self, 'h_line'):
            self.h_line.set_visible(False)
        if hasattr(self, 'v_line'):
            self.v_line.set_visible(False)
        self.canvas.draw_idle()

    def on_save_button(self,event):

        savefile, parameter_name = self.diffusion_service.parameter_file()
        print('Saving to:',savefile)
        
        write={}
        write['noiseFacDiff'] = self.noiseFac.GetValue()
        # Store ROI geometry only; all derived diffusion statistics are rebuilt
        # from the current data when the frame is next opened.  A compact
        # semicolon-separated list keeps the legacy parameter file to one token:
        #     diffusionROIs = min,max;min,max;...
        # A merged ROI is simply saved as its final merged min/max span.
        write['diffusionROIs'] = ';'.join(
            '{:.12g},{:.12g}'.format(float(a), float(b))
            for a, b in self.roi_ranges
        ) if self.roi_ranges else 'none'

        update_parameter_file(savefile, write, source_path=parameter_name)


    def set_default_values(self): #upack Grp save

        #print('TESTING:',t,t==0 )
        #sys.exit(100)
        # NoiseFac is the only diffusion-specific noise parameter.  The base
        # Base noise is supplied by the diffusion application service.
        if self.diffusion_service.parameter('noiseFacDiff') != 0:
            self.noiseFac.SetValue(str(self.diffusion_service.parameter('noiseFacDiff', numeric=True)))
        else:
            self.noiseFac.SetValue('1.0')

        # ROI data cannot be rebuilt until AnalDiff has populated the accepted
        # per-ppm diffusion arrays, so parse the geometry now and restore it
        # immediately after the initial draw_figure().
        self._saved_roi_ranges = []
        raw_rois = str(self.diffusion_service.parameter('diffusionROIs', default='') or '').strip()
        if raw_rois and raw_rois.lower() not in ('0', 'none'):
            for item in raw_rois.split(';'):
                try:
                    a_text, b_text = item.split(',', 1)
                    a, b = float(a_text), float(b_text)
                except (TypeError, ValueError):
                    continue
                if numpy.isfinite(a) and numpy.isfinite(b):
                    self._saved_roi_ranges.append((a, b))

    def _restore_saved_rois(self):
        """Rebuild saved ROI regions and every result derived from them."""
        ranges = list(getattr(self, '_saved_roi_ranges', []))
        if not ranges:
            return
        self.roi_ranges = []
        for a, b in ranges:
            # Snap saved ppm limits to the current spectral grid.  This also
            # makes saved regions robust to harmless ppm precision changes.
            ia = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - a)))
            ib = int(numpy.argmin(numpy.abs(self.diffusion_service.spectral_axis - b)))
            self._append_roi(float(self.diffusion_service.spectral_axis[ia]),
                             float(self.diffusion_service.spectral_axis[ib]))
        self.number_scatters = len(self.scatter_data)
        self._highlighted_rois = set()
        self._rebuild_roi_overlays()
        self.plot_scatters()
        self._plot_roi_histograms()
        self._refresh_roi_manager()
        self.canvas.draw_idle()
