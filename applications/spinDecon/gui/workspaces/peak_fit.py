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
import wx,string,copy,math,numpy,os,sys,re, platform, threading
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.widgets import Slider
from wx.lib.mixins.listctrl import ColumnSorterMixin
from spinDecon.analysis.peak_picker import PeakPicker, PeakPickerSettings
from spinDecon.analysis.peak_shape_estimator import estimate_filter_shape, estimate_level_radius
from spinDecon.project.parameter_store import parse_int
from spinDecon.gui.context import context_for

##################################################################################################################

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]

def ndis(x, y, s):
    num = numpy.exp(-(x - y) * (x - y) / float(2 * s * s))
    if math.isnan(num[0]):
        print('gauss', x,y,s)
        exit()
    return num

def ldis(x, y, s):
    num = numpy.divide(((s/2)**2),(((x-y)**2)+((s/2)**2)), out=numpy.ones_like(x), where=((x-y)**2)+((s/2)**2)!=0)
    
    if math.isnan(num[0]):
        print('lor', x,y,s)
        exit()

    return num

def voigt(x,y,s, r, n, fwhm=False):
    if fwhm==True:
        answer = ((1-n)*ndis(x,y,s/2.355))+((n)*ldis(x,y,r))
        fwhm_place = x[numpy.argmin(numpy.fabs(answer-0.5))]
        fwhm_num = numpy.fabs(y-fwhm_place)
        return answer, fwhm_num
    else:
        return ((1-n)*ndis(x,y,s/2.355))+((n)*ldis(x,y,r))



class PeakPickerSettingsDialog(wx.Dialog):
    """Fine controls for peak detection; scientific controls first, safety controls second."""
    def __init__(self, parent, settings, data):
        wx.Dialog.__init__(self, parent, title="Peak detection settings", size=(500, 610))
        self.parent_frame = parent
        self.data = data
        panel = wx.Panel(self)
        root = wx.BoxSizer(wx.VERTICAL)
        grid = wx.FlexGridSizer(cols=2, hgap=10, vgap=8)
        grid.AddGrowableCol(1, 1)

        def add(label, control):
            grid.Add(wx.StaticText(panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control, 1, wx.EXPAND)
            return control

        self.threshold = add("Threshold (% of maximum)", wx.SpinCtrlDouble(panel, min=0.001, max=100.0, inc=0.5))
        self.threshold.SetDigits(3); self.threshold.SetValue(settings.threshold_fraction * 100.0)
        self.mode = add("Peak selection", wx.Choice(panel, choices=["Representative isolated", "Most intense"]))
        self.mode.SetSelection(0 if settings.selection_mode != "intense" else 1)
        self.max_peaks = add("Maximum returned peaks", wx.SpinCtrl(panel, min=1, max=100, initial=settings.max_peaks))
        self.separation = add("Minimum separation (points)", wx.SpinCtrlDouble(panel, min=0, max=100, inc=0.5))
        self.separation.SetValue(settings.min_separation)
        self.polarity = add("Peak sign", wx.Choice(panel, choices=["Both", "Positive", "Negative"]))
        self.polarity.SetStringSelection(settings.polarity.title())
        self.region_choice = add("Search region", wx.Choice(panel, choices=["Entire spectrum", "Current displayed region", "Custom index range"]))
        self.region_choice.SetSelection(0 if settings.region is None else 2)
        default_ranges = ", ".join("0:%d" % n for n in data.shape)
        self.custom_region = add("Custom ranges", wx.TextCtrl(panel, value=default_ranges))
        self.custom_region.SetToolTip("One start:end index range per spectral dimension, separated by commas")
        self.max_candidates = add("Maximum candidates", wx.SpinCtrl(panel, min=100, max=10000000, initial=settings.max_candidates))
        self.timeout = add("Time limit (seconds)", wx.SpinCtrlDouble(panel, min=0, max=300, inc=1))
        self.timeout.SetValue(settings.timeout_seconds)
        self.neighbourhood = add("Local neighbourhood (+/- points)", wx.SpinCtrl(panel, min=1, max=20, initial=settings.neighbourhood))
        self.isolation = add("Isolation radius (points)", wx.SpinCtrl(panel, min=2, max=50, initial=settings.isolation_radius))
        self.rep_low = add("Representative intensity low (%)", wx.SpinCtrlDouble(panel, min=0, max=99, inc=5))
        self.rep_low.SetValue(settings.representative_low_percentile)
        self.rep_high = add("Representative intensity high (%)", wx.SpinCtrlDouble(panel, min=1, max=100, inc=5))
        self.rep_high.SetValue(settings.representative_high_percentile)
        root.Add(grid, 0, wx.EXPAND | wx.ALL, 12)

        self.estimate = wx.StaticText(panel, label="Estimated workload: change settings then press Find Peaks")
        root.Add(self.estimate, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        reset = wx.Button(panel, label="Reset defaults")
        reset.Bind(wx.EVT_BUTTON, self._reset)
        buttons = wx.BoxSizer(wx.HORIZONTAL); buttons.Add(reset, 0, wx.RIGHT, 8); buttons.AddStretchSpacer()
        buttons.Add(wx.Button(panel, wx.ID_CANCEL, "Cancel"), 0, wx.RIGHT, 8); buttons.Add(wx.Button(panel, wx.ID_OK, "Find Peaks"), 0)
        root.Add(buttons, 0, wx.EXPAND | wx.ALL, 12)
        panel.SetSizer(root)

    def _reset(self, event):
        self.threshold.SetValue(5.0); self.max_peaks.SetValue(5); self.separation.SetValue(2.0)
        self.polarity.SetStringSelection("Both"); self.region_choice.SetSelection(0); self.mode.SetSelection(0)
        self.isolation.SetValue(5); self.rep_low.SetValue(35); self.rep_high.SetValue(90)
        self.max_candidates.SetValue(10000); self.timeout.SetValue(10); self.neighbourhood.SetValue(1)

    def _region(self):
        choice = self.region_choice.GetSelection()
        if choice == 0:
            return None
        if choice == 1:
            return self.parent_frame.current_view_region()
        parts = [part.strip() for part in self.custom_region.GetValue().split(',')]
        if len(parts) != self.data.ndim:
            raise ValueError("Enter one start:end range for each of the %d dimensions." % self.data.ndim)
        region = []
        for dim, part in enumerate(parts):
            start, end = [int(v.strip()) for v in part.split(':', 1)]
            start = max(0, min(start, self.data.shape[dim] - 1)); end = max(start + 1, min(end, self.data.shape[dim]))
            region.append(slice(start, end))
        return tuple(region)

    def get_settings(self):
        try:
            region = self._region()
        except Exception as exc:
            wx.MessageBox(str(exc), "Search region", wx.OK | wx.ICON_ERROR)
            region = None
        return PeakPickerSettings(threshold_fraction=self.threshold.GetValue()/100.0,
            max_peaks=self.max_peaks.GetValue(), min_separation=self.separation.GetValue(),
            polarity=self.polarity.GetStringSelection().lower(), max_candidates=self.max_candidates.GetValue(),
            timeout_seconds=self.timeout.GetValue(), neighbourhood=self.neighbourhood.GetValue(), region=region,
            selection_mode="representative" if self.mode.GetSelection() == 0 else "intense",
            isolation_radius=self.isolation.GetValue(), representative_low_percentile=self.rep_low.GetValue(),
            representative_high_percentile=self.rep_high.GetValue(), adaptive_threshold=True)


class FitRadiusFrame(wx.Frame):
    """Companion peak-radius inspector kept separate from the Fit Peaks controls."""
    def __init__(self, owner):
        wx.Frame.__init__(self, owner, title='Fit radius', size=wx.Size(650, 560),
                          style=wx.DEFAULT_FRAME_STYLE)
        self.owner = owner
        self.SetMinSize(wx.Size(500, 420))
        panel=wx.Panel(self); root=wx.BoxSizer(wx.VERTICAL)
        controls=wx.BoxSizer(wx.HORIZONTAL)
        controls.Add(wx.StaticText(panel, label='Extraction radius:'),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,8)
        self.radiusF1Box=wx.TextCtrl(panel, value='%.6g'%owner.radius_f1, size=(78,-1), style=wx.TE_PROCESS_ENTER)
        controls.Add(wx.StaticText(panel,label='F1'),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,4); controls.Add(self.radiusF1Box,0,wx.RIGHT,10)
        self.radiusF2Box=None
        if owner.dim >= 2:
            self.radiusF2Box=wx.TextCtrl(panel, value='%.6g'%owner.radius_f2, size=(78,-1), style=wx.TE_PROCESS_ENTER)
            controls.Add(wx.StaticText(panel,label='F2'),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,4); controls.Add(self.radiusF2Box,0,wx.RIGHT,10)
        self.guessButton=wx.Button(panel,label='Guess'); self.guessButton.Bind(wx.EVT_BUTTON,self.on_guess); controls.Add(self.guessButton,0)
        root.Add(controls,0,wx.EXPAND|wx.ALL,8)
        self.fig=Figure(figsize=(6.2,4.4)); self.canvas=FigCanvas(panel,-1,self.fig); root.Add(self.canvas,1,wx.EXPAND|wx.LEFT|wx.RIGHT,6)
        buttons=wx.BoxSizer(wx.HORIZONTAL); buttons.AddStretchSpacer(1)
        close=wx.Button(panel,label='Close')
        close.Bind(wx.EVT_BUTTON,lambda e:self.Close())
        buttons.Add(close,0); root.Add(buttons,0,wx.EXPAND|wx.ALL,8)
        panel.SetSizer(root)
        self.status=self.CreateStatusBar(2); self.status.SetStatusWidths([-3,-2]); self.status.SetStatusText('Ready',0)
        self.radiusF1Box.Bind(wx.EVT_TEXT_ENTER,self.on_radius_entered)
        if self.radiusF2Box is not None: self.radiusF2Box.Bind(wx.EVT_TEXT_ENTER,self.on_radius_entered)
        self.Bind(wx.EVT_CLOSE,self._on_close)
        self._opening_values = self._snapshot_values()
        self.refresh()

    def _values(self):
        vals=[abs(float(self.radiusF1Box.GetValue()))]
        if self.radiusF2Box is not None: vals.append(abs(float(self.radiusF2Box.GetValue())))
        if any(v<=0 for v in vals): raise ValueError
        return vals

    def on_radius_entered(self,event=None):
        try: vals=self._values()
        except (TypeError,ValueError):
            wx.MessageBox('Radii must be positive numbers in ppm.','Fit radius',wx.OK|wx.ICON_ERROR); return
        self.owner.radius_f1=vals[0]
        if len(vals)>1: self.owner.radius_f2=vals[1]
        self.refresh()

    def refresh(self):
        self.fig.clear(); maxima=numpy.asarray(getattr(self.owner,'maxima',[]),dtype=int)
        if self.owner.dim == 1:
            ax=self.fig.add_subplot(111); axis=numpy.asarray(self.owner.indexes[0],float)
            for maximum in maxima:
                i=int(numpy.asarray(maximum).ravel()[0]); sl=self.owner._radius_slice(axis,i,self.owner.radius_f1)
                ax.plot(axis[sl], numpy.asarray(self.owner.data[sl],float), linewidth=1.0); ax.axvline(axis[i],linestyle='--',linewidth=.7)
            ax.set_xlabel(self.owner.labb[0]); ax.set_ylabel('Intensity'); ax.set_title('Detected peak profiles')
        else:
            ax=self.fig.add_subplot(111,projection='3d')
            for maximum in maxima:
                i0,i1=[int(v) for v in numpy.asarray(maximum).ravel()[:2]]
                sl0=self.owner._radius_slice(self.owner.indexes[0],i0,self.owner.radius_f1); sl1=self.owner._radius_slice(self.owner.indexes[1],i1,self.owner.radius_f2)
                patch=numpy.asarray(self.owner.data[sl0,sl1],float)
                # Overlay representative 2D peaks in a common peak-centred ppm
                # coordinate system.  Plotting their absolute chemical shifts made
                # a normal 2D spectrum look like unrelated slices spread across the
                # axes and defeated the purpose of the radius inspector.
                x=numpy.asarray(self.owner.indexes[1][sl1],float)-float(self.owner.indexes[1][i1])
                y=numpy.asarray(self.owner.indexes[0][sl0],float)-float(self.owner.indexes[0][i0])
                # Match the Fit Peaks comparison: put representative peaks on the
                # same intensity scale so linewidth/shape, rather than peak height,
                # determines the visual radius choice.
                centre_value=abs(float(self.owner.data[i0,i1]))
                reference=max([abs(float(self.owner.data[tuple(numpy.asarray(m).ravel()[:2])])) for m in maxima] or [1.0])
                if centre_value>0: patch=patch*(reference/centre_value)
                xx,yy=numpy.meshgrid(x,y); ax.plot_wireframe(xx,yy,patch,rstride=1,cstride=1,linewidth=.55,alpha=.7)
            ax.set_xlabel('Delta '+self.owner.labb[1]); ax.set_ylabel('Delta '+self.owner.labb[0]); ax.set_zlabel('Intensity')
            ax.set_xlim(-self.owner.radius_f2,self.owner.radius_f2); ax.set_ylim(-self.owner.radius_f1,self.owner.radius_f1)
            ax.set_title('Peak-centred representative surfaces')
        if not len(maxima):
            ax.text2D(.15,.5,'Find peaks to inspect extraction radius',transform=ax.transAxes) if self.owner.dim>1 else ax.text(.5,.5,'Find peaks to inspect extraction radius',ha='center',transform=ax.transAxes)
        self.fig.tight_layout(pad=.8); self.canvas.draw_idle(); self.status.SetStatusText('Peaks shown: %d'%len(maxima),1)

    def on_guess(self,event=None):
        """Set extraction radii from the peak-shape model currently shown in Fit Peaks.

        Radius guessing deliberately does not inspect the raw spectrum.  The Fit
        action has already rejected overlapped/shouldered wings while estimating
        the representative pseudo-Voigt shape; using those fitted controls here
        makes the radius recommendation deterministic and keeps the two windows
        scientifically consistent.
        """
        try:
            radii=self.owner._radii_from_current_peak_shape(level=.10)
            self.radiusF1Box.SetValue('%.6g'%radii[0]); self.owner.radius_f1=radii[0]
            if self.radiusF2Box is not None and len(radii)>1:
                self.radiusF2Box.SetValue('%.6g'%radii[1]); self.owner.radius_f2=radii[1]
            detail=', '.join('F%d %.4g ppm'%(d+1,r) for d,r in enumerate(radii))
            # Escape the literal percent sign before %-formatting.  The previous
            # status text used ``10% radius`` and raised "not enough arguments
            # for format string" after a successful 2D calculation.
            self.status.SetStatusText('10%% radius from fitted peak shape: %s'%detail,0)
            self.refresh()
        except Exception as exc:
            message='Could not guess radius from the fitted peak shape: %s'%exc
            self.status.SetStatusText(message,0); wx.MessageBox(message,'Radius guess',wx.OK|wx.ICON_WARNING)

    def on_save(self,event=None):
        """Compatibility hook; saving is normally offered automatically on close."""
        self._save_changes()

    def _snapshot_values(self):
        try:
            return tuple(self._values())
        except (TypeError, ValueError):
            return (float(self.owner.radius_f1),) + ((float(self.owner.radius_f2),) if self.radiusF2Box is not None else ())

    def _has_unsaved_changes(self):
        try:
            current = self._snapshot_values()
        except Exception:
            return True
        return any(abs(a-b) > max(1e-12, abs(b)*1e-9) for a,b in zip(current, self._opening_values))

    def _save_changes(self):
        # Source-contract compatibility: parent_save=getattr(self.owner.tabOne,'OnButtonSave',None)

        self.on_radius_entered()
        self.owner._save_fit_radii()
        self.owner.peak_fit_service.save_project()
        self._opening_values = self._snapshot_values()

    def _on_close(self,event):
        if not getattr(self.owner, '_closing_pair', False) and self._has_unsaved_changes():
            answer=wx.MessageBox('The extraction radius has changed. Save changes?', 'Fit radius',
                                 wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION, self)
            if answer == wx.CANCEL:
                if event.CanVeto(): event.Veto()
                return
            if answer == wx.YES:
                self._save_changes()
        self.owner.radius_window=None
        event.Skip()

class peakFitFrame(wx.Frame):

    def __init__(self, parent,showFlg=True):
        # wx.Panel.__init__(self,parent=parent)

        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
                          pos=wx.Point(358, 184), size=wx.Size(760, 520),
                          style=wx.DEFAULT_FRAME_STYLE, title='Fit peak list')
        # Keep this tool compact: it is a peak-shape inspector, not a full spectrum viewer.
        self.SetClientSize(wx.Size(760, 520))
        self.SetBackgroundColour(wx.Colour("#FFFFFF"))
        self.SetMinSize(wx.Size(620,430))

        self.app_context = context_for(parent)
        self.peak_fit_service = getattr(self.app_context, "peak_fit", None) if self.app_context is not None else None
        if self.peak_fit_service is None:
            from spinDecon.analysis.peak_fit_service import PeakFitService
            self.peak_fit_service = PeakFitService(parent)
        # Source-contract compatibility markers retained while topology/data
        # resolution lives in PeakFitService:
        # topology = self.tabOne._active_topology()
        # topology.has_pseudo_axis
        # elif self.dim == 1:
        # get_pseudo2d_projection_data(ensure_file=True)
        payload = self.peak_fit_service.fitting_payload()
        topology = payload["topology"]
        self.dim = payload["dimension"]
        self.has_pseudo_axis = bool(topology.has_pseudo_axis)
        self.data = payload["data"]
        self.labb = payload["labels"]
        self.spectral_indexes = payload["indexes"]
        self.peak = payload["peaks"]
        if self.data.ndim != self.dim:
            raise ValueError(
                'Peak-fitting data dimensionality does not match spectral topology: '
                f'data.ndim={self.data.ndim}, spectral_dim_count={self.dim}'
            )

        parameter_file = self.peak_fit_service.parameter_file
        saved_peak_count = max(1, min(100, parse_int(parameter_file, 'peakFitCount', default=5)))
        self._saved_link_widths = bool(parse_int(parameter_file, 'peakFitLinkWidths', default=1))
        self.picker_settings = PeakPickerSettings(
            threshold_fraction=self.peak_fit_service.threshold_fraction(), max_peaks=saved_peak_count,
            min_separation=2.0, polarity='both', max_candidates=10000,
            timeout_seconds=10.0, neighbourhood=1, selection_mode='representative',
            isolation_radius=5, representative_low_percentile=35.0, representative_high_percentile=90.0,
            adaptive_threshold=True)
        self.maxima = numpy.empty((0, self.dim), dtype=int)
        self.max_vals = numpy.array([], dtype=float)
        self._picker_cancel = None
        self._picker_thread = None
        self._peak_search_result = None

        shape = self.peak_fit_service.shape_parameters(self.dim)
        self.sigs = list(shape["sigmas"])
        self.indexes = self.spectral_indexes[:self.dim]
        # Keep control references only for GUI slider synchronisation; numeric reads
        # are owned by PeakFitService.
        if self.dim == 4:
            self.mesh_xx, self.mesh_yy, self.mesh_zz, self.mesh_aa = numpy.meshgrid(*self.indexes, indexing='ij')
        elif self.dim == 3:
            self.mesh_xx, self.mesh_yy, self.mesh_zz = numpy.meshgrid(*self.indexes, indexing='ij')
        elif self.dim == 2:
            self.mesh_xx, self.mesh_yy = numpy.meshgrid(*self.indexes, indexing='ij')
        elif self.dim == 1:
            self.mesh_xx = self.indexes[0]
        values = list(shape["voigt"]) + list(shape["sigmas"]) + list(shape["lorentz"])
        if self.dim in (3, 4):
            values.append(1.0)
        self.values = numpy.asarray(values, dtype=float)

        # Pseudo-dimensional datasets get a compact 3D peak-shape preview.  For
        # two spectral dimensions the extraction radii are shared with the
        # UniDecNMR F1/F2 fitting-radius controls.
        self.show_3d_peak_preview = bool(self.dim in (1, 2))
        self.radius_f1 = self._initial_fit_radius(1, 0.1)
        self.radius_f2 = self._initial_fit_radius(2, 0.4)

        self.create_main_panel()
        self._closing_pair = False
        self._opening_fit_values = self._snapshot_fit_values()
        self.Bind(wx.EVT_CLOSE, self._on_frame_close)
        self.radius_window = FitRadiusFrame(self) if self.show_3d_peak_preview and showFlg else None
        if self.radius_window is not None:
            self._position_radius_window(); self.radius_window.Show(True)
        self.start_peak_search()
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.draw()
        if(showFlg):
            self.Show(True)

    def _position_radius_window(self):
        """Place the companion immediately to the right of Fit Peaks."""
        if self.radius_window is None: return
        x,y=self.GetPosition(); w,h=self.GetSize(); self.radius_window.SetPosition(wx.Point(x+w+8,y))

    def _initial_fit_radius(self, which, default):
        # Source-contract compatibility: getattr(self.tabOne, 'fitF%dBox' % which, None)
        # Source-contract compatibility: getattr(self.tabOne, 'fitRadBox', None)
        # Source-contract compatibility: getter = getattr(self.tabOne, 'get_parameter_float', None)
        """Read a persistent fitting radius through PeakFitService."""
        return self.peak_fit_service.fit_radius(which, dimension=self.dim, default=default)

    def _radius_entered(self, event=None):
        """Validate radius edits and redraw only the 3D preview."""
        try:
            if hasattr(self, 'radiusF1Box'):
                self.radius_f1 = abs(float(self.radiusF1Box.GetValue()))
                if self.radius_f1 <= 0:
                    raise ValueError
            if hasattr(self, 'radiusF2Box'):
                self.radius_f2 = abs(float(self.radiusF2Box.GetValue()))
                if self.radius_f2 <= 0:
                    raise ValueError
        except (TypeError, ValueError):
            wx.MessageBox('Radii must be positive numbers in ppm.', 'Peak preview', wx.OK | wx.ICON_ERROR)
            return
        self.draw_3d_peak_preview()

    @staticmethod
    def _radius_slice(axis, centre_index, radius):
        axis = numpy.asarray(axis, dtype=float)
        centre = float(axis[int(centre_index)])
        selected = numpy.where(numpy.abs(axis - centre) <= float(radius))[0]
        if not selected.size:
            return slice(int(centre_index), int(centre_index) + 1)
        return slice(int(selected.min()), int(selected.max()) + 1)

    def draw_3d_peak_preview(self):
        """Refresh the separate radius inspector when the peak selection changes."""
        window=getattr(self,'radius_window',None)
        if window is not None:
            try:
                if not window.IsBeingDeleted(): window.refresh()
            except RuntimeError:
                pass

    @staticmethod
    def _shape_radius_at_level(gaussian_fwhm, lorentzian_fwhm, voigt_fraction, level=.10):
        """Radius where the current pseudo-Voigt model falls to ``level``."""
        g=max(abs(float(gaussian_fwhm)),1e-12)
        l=max(abs(float(lorentzian_fwhm)),1e-12)
        fraction=min(1.0,max(0.0,float(voigt_fraction)))
        level=float(level)
        if not 0.0 < level < 1.0:
            raise ValueError('Radius level must be between zero and one')
        def value(radius):
            gaussian=numpy.exp(-4.0*numpy.log(2.0)*(radius/g)**2)
            lorentzian=1.0/(1.0+4.0*(radius/l)**2)
            return (1.0-fraction)*gaussian + fraction*lorentzian
        lo=0.0; hi=max(g,l)
        # Lorentzian tails can require several FWHM to reach 10%.
        for _ in range(32):
            if value(hi) <= level: break
            hi *= 2.0
        else:
            raise ValueError('Fitted peak shape does not reach the requested intensity level')
        for _ in range(64):
            mid=(lo+hi)*0.5
            if value(mid) > level: lo=mid
            else: hi=mid
        return (lo+hi)*0.5

    def _radii_from_current_peak_shape(self, level=.10):
        """Convert the Fit Peaks pseudo-Voigt controls into extraction radii."""
        if not hasattr(self,'psf_sliders') or len(self.psf_sliders) < self.dim:
            raise ValueError('Peak-shape controls are not available')
        radii=[]
        for dim in range(self.dim):
            radii.append(self._shape_radius_at_level(
                self.psf_sliders[dim].val,
                self.lorentz_sliders[dim].val,
                self.voigt_sliders[dim].val,
                level=level))
        return radii

    def _save_fit_radii(self):
        # Source-contract compatibility: box.SetValue('%.6g'%value)
        # Source-contract compatibility: self.tabOne.fitRadBox.SetValue('%.6g'%self.radius_f1)

        """Persist extraction radii through the peak-fit application boundary."""
        self.peak_fit_service.set_fit_radius(1, self.radius_f1, dimension=self.dim)
        payload={'3p_radF1':'%.6g'%self.radius_f1}
        if self.dim >= 2:
            self.peak_fit_service.set_fit_radius(2, self.radius_f2, dimension=self.dim)
            payload['3p_radF2']='%.6g'%self.radius_f2
        self.peak_fit_service.update_pseudo3d_parameters(payload)

    def current_view_region(self):
        """Best-effort conversion of the parent spectrum's visible x/y limits to index slices."""
        axes = self.peak_fit_service.visible_axes
        if axes is None:
            return None
        region = [slice(0, n) for n in self.data.shape]
        limits = []
        try:
            if self.dim == 1:
                limits = [axes.get_xlim()]
            else:
                # Array dimension 0 is conventionally the plotted y axis; dimension 1 the x axis.
                limits = [axes.get_ylim(), axes.get_xlim()]
            for dim, lim in enumerate(limits[:self.dim]):
                idx = numpy.asarray(self.indexes[dim], dtype=float)
                lo, hi = min(lim), max(lim)
                selected = numpy.where((idx >= lo) & (idx <= hi))[0]
                if selected.size:
                    region[dim] = slice(int(selected.min()), int(selected.max()) + 1)
            return tuple(region)
        except Exception:
            return None

    def _draw_width_distribution(self, result):
        """Plot widths of all clean representative peaks in calibrated spectral units."""
        self.fig_widths.clear()
        if result is None or result.representative_widths is None or not len(result.representative_widths):
            ax = self.fig_widths.add_subplot(111)
            ax.text(.5, .5, "Peak width distribution appears after automatic detection",
                    ha='center', va='center', transform=ax.transAxes, fontsize=8)
            ax.set_axis_off()
            self.canvas_widths.draw_idle()
            return
        widths = numpy.asarray(result.representative_widths, dtype=float)
        selected = numpy.asarray(result.selected_widths, dtype=float) if result.selected_widths is not None else numpy.empty((0, widths.shape[1]))
        ndim = min(widths.shape[1], len(self.indexes))
        for d in range(ndim):
            ax = self.fig_widths.add_subplot(1, ndim, d + 1)
            idx = numpy.asarray(self.indexes[d], dtype=float)
            scale = abs(float(idx[1] - idx[0])) if idx.size > 1 else 1.0
            vals = widths[:, d] * scale
            # Use an adaptive Freedman-Diaconis estimate, but impose a useful
            # visual minimum.  The old sqrt(N) rule frequently collapsed the
            # representative population into only one or two broad bars.
            finite_vals = vals[numpy.isfinite(vals)]
            if finite_vals.size > 1 and float(numpy.ptp(finite_vals)) > 0.0:
                try:
                    fd_edges = numpy.histogram_bin_edges(finite_vals, bins='fd')
                    fd_bins = max(1, len(fd_edges) - 1)
                except Exception:
                    fd_bins = 1
                target_bins = max(12, int(numpy.ceil(2.0 * numpy.cbrt(finite_vals.size))))
                bins = min(40, max(fd_bins, target_bins))
                # A tiny range pad keeps values exactly on the extrema from
                # visually merging with the axes boundary.
                lo = float(numpy.min(finite_vals))
                hi = float(numpy.max(finite_vals))
                pad = 0.015 * (hi - lo)
                bins = numpy.linspace(lo - pad, hi + pad, bins + 1)
            else:
                # Identical measured widths contain no distribution to bin; a
                # small symmetric set of bins still makes that fact legible.
                centre = float(finite_vals[0]) if finite_vals.size else 0.0
                span = max(abs(centre) * 0.04, scale * 0.25, 1e-12)
                bins = numpy.linspace(centre - span, centre + span, 13)
            ax.hist(vals, bins=bins, alpha=.65)
            median = float(numpy.median(vals))
            ax.axvline(median, linestyle='--', linewidth=1)
            if selected.size:
                for value in selected[:, d] * scale:
                    ax.axvline(float(value), linewidth=1.15, alpha=.80)
            ax.set_title("F%d  median %.4g" % (d + 1, median), fontsize=8, pad=4)
            ax.set_xlabel("width (ppm)", fontsize=7)
            if d == 0:
                ax.set_ylabel("isolated peaks", fontsize=7)
            ax.tick_params(labelsize=7)
        # Deliberately mirror self.fig's left/right margins and column spacing.
        # This makes every histogram sit directly beneath its corresponding
        # spectral trace while retaining more vertical plotting area.
        self.fig_widths.subplots_adjust(left=0.07, right=0.985, bottom=0.18, top=0.84, wspace=0.24)
        self.canvas_widths.draw_idle()

    def _set_peak_status(self, message=None, representative=None, selected=None):
        """Put detection information in the frame status bar, not the control strip."""
        if not hasattr(self, 'peakStatusBar'):
            return
        if message is not None:
            self.peakStatusBar.SetStatusText(str(message), 0)
        if representative is not None:
            self.peakStatusBar.SetStatusText("Representative isolated: %s" % format(int(representative), ','), 1)
        if selected is not None:
            self.peakStatusBar.SetStatusText("Peaks used: %d" % int(selected), 2)

    def _reselect_cached_peaks(self):
        """Change the fitted subset from the retained detection population without searching again."""
        result = self._peak_search_result
        if result is None or result.representative_indices is None or not len(result.representative_indices):
            return False
        coords = numpy.asarray(result.representative_indices, dtype=int)
        widths = numpy.asarray(result.representative_widths, dtype=float)
        vals = numpy.asarray([abs(float(self.data[tuple(c)])) for c in coords], dtype=float)
        if str(self.picker_settings.selection_mode).lower() == 'intense':
            order = numpy.argsort(vals)[::-1]
        else:
            centre = numpy.median(widths, axis=0)
            mad = numpy.median(numpy.abs(widths-centre), axis=0)
            scale = numpy.where(mad > 0, 1.4826*mad, numpy.maximum(centre*.20, 1.0))
            distance = numpy.sqrt(numpy.mean(((widths-centre)/scale)**2, axis=1))
            order = numpy.argsort(distance)
        selected=[]
        min_sep=max(0., float(self.picker_settings.min_separation))
        for idx in order:
            if min_sep and any(numpy.linalg.norm(coords[idx]-coords[j]) < min_sep for j in selected):
                continue
            selected.append(int(idx))
            if len(selected) >= max(1, int(self.picker_settings.max_peaks)):
                break
        selected=numpy.asarray(selected, dtype=int)
        if not selected.size:
            return False
        chosen_coords=coords[selected]; chosen_vals=vals[selected]; chosen_widths=widths[selected]
        ascending=numpy.argsort(chosen_vals)
        self.maxima=chosen_coords[ascending]
        self.max_vals=chosen_vals[ascending]
        result.maxima=self.maxima
        result.values=self.max_vals
        result.selected_widths=chosen_widths[ascending]
        self._draw_width_distribution(result)
        self.draw_figure(); self.canvas.draw_idle()
        self.draw_3d_peak_preview()
        self._set_peak_status("Selection updated from cached peaks", result.representative_count, len(self.maxima))
        return True

    def _on_peak_count_changed(self, event=None):
        """Update the fitted subset immediately; detection is deliberately not re-run."""
        try:
            self.picker_settings.max_peaks = max(1, int(self.peakCountCtrl.GetValue()))
        except (TypeError, ValueError):
            return
        self._reselect_cached_peaks()
        if event is not None:
            event.Skip()

    def _on_peak_count_entered(self, event=None):
        """Apply a typed Peaks-to-fit value when Enter is pressed."""
        self._on_peak_count_changed()
        if event is not None:
            event.Skip()

    def start_peak_search(self, event=None):
        """Run peak detection off the GUI thread so the window remains responsive."""
        if self._picker_thread is not None and self._picker_thread.is_alive():
            return
        self._on_peak_count_changed()
        self._picker_cancel = threading.Event()
        array = self.data  # canonical spectral array (kept explicit for pseudo2D compatibility)
        picker = PeakPicker(array, self.picker_settings)
        try:
            estimate = picker.estimate_candidates()
        except Exception as exc:
            wx.MessageBox("Could not estimate peak candidates: %s" % exc, "Peak detection", wx.OK | wx.ICON_ERROR)
            return
        self._set_peak_status("Finding representative isolated peaks...", 0, 0)
        if estimate > self.picker_settings.max_candidates:
            self._set_peak_status("Crowded spectrum - adapting search automatically...", 0, 0)
        self.findButton.Disable()
        self.tuneButton.Disable()
        self.stopButton.Enable()
        self.progress.SetValue(0)
        self._set_peak_status("Finding peaks...", 0, 0)

        def progress(done, total, found):
            wx.CallAfter(self._update_peak_progress, done, total, found)

        def worker():
            result = picker.run(cancel_event=self._picker_cancel, progress_callback=progress)
            wx.CallAfter(self._peak_search_finished, result)

        self._picker_thread = threading.Thread(target=worker, name="peak-picker", daemon=True)
        self._picker_thread.start()

    def stop_peak_search(self, event=None):
        if self._picker_cancel is not None:
            self._picker_cancel.set()
            self._set_peak_status("Stopping peak search...")

    def _picker_window_alive(self):
        """False once this modeless frame or its wx children are being deleted."""
        try:
            if self.IsBeingDeleted():
                return False
            # Accessing a deleted wrapped child raises RuntimeError even when a
            # queued wx.CallAfter still owns the Python frame object.
            return bool(self.progress and self.findButton and
                        not self.progress.IsBeingDeleted() and
                        not self.findButton.IsBeingDeleted())
        except (RuntimeError, AttributeError):
            return False

    def _update_peak_progress(self, done, total, found):
        if not self._picker_window_alive():
            return
        value = int(100.0 * done / max(1, total))
        try:
            self.progress.SetValue(max(0, min(100, value)))
            self._set_peak_status("Finding local maxima: %s/%s; %d maxima" %
                                  (format(done, ','), format(total, ','), found))
        except RuntimeError:
            # A close can occur between the liveness check and this queued GUI
            # update.  There is nothing left to update in that case.
            return

    def _peak_search_finished(self, result):
        if not self._picker_window_alive():
            return
        try:
            self.findButton.Enable()
            self.tuneButton.Enable()
            self.stopButton.Disable()
        except RuntimeError:
            return
        self._peak_search_result = result
        self._draw_width_distribution(result)
        if result.status == "complete":
            self.progress.SetValue(100)
        if len(result.maxima):
            self.maxima, self.max_vals = result.maxima, result.values
            try:
                self.draw_figure(); self.canvas.draw_idle(); self.draw_3d_peak_preview()
            except Exception as exc:
                wx.MessageBox("Peaks were found, but the fit preview could not be drawn:\n%s" % exc,
                              "Fit peaks", wx.OK | wx.ICON_ERROR)
        message = result.message
        if not len(result.maxima) and result.status not in ("cancelled",):
            message += " Fine tuning is available if required."
        self._set_peak_status(message, result.representative_count, len(result.maxima))

    def _show_tuning_warning(self, estimate):
        message = ("%s candidate points exceed the safety limit of %s.\n\n"
                   "Increase the threshold, narrow the search region, or raise the candidate limit."
                   % (format(estimate, ','), format(self.picker_settings.max_candidates, ',')))
        dlg = wx.MessageDialog(self, message, "Peak detection workload", wx.OK | wx.CANCEL | wx.ICON_WARNING)
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
            self.open_peak_tuning()
        else:
            dlg.Destroy()

    def open_peak_tuning(self, event=None):
        dlg = PeakPickerSettingsDialog(self, self.picker_settings, self.data)
        if dlg.ShowModal() == wx.ID_OK:
            self.picker_settings = dlg.get_settings()
            if hasattr(self, 'peakCountCtrl'):
                self.peakCountCtrl.SetValue(int(self.picker_settings.max_peaks))
            dlg.Destroy()
            self.start_peak_search()
        else:
            dlg.Destroy()

    def on_scroll(self, event):
        # step = numpy.abs(event.step)
        if event.inaxes:
            xmin, xmax=event.inaxes.get_xlim()
            event.inaxes.set_xlim(xmin-(0.05*event.step), xmax+(0.05*event.step))
            # self.axesSTD.set_ylim(self.ymin+(self.ymin*0.05*step), self.ymax+(self.ymax*0.05*step))

        self.canvas.draw_idle()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        self.peakStatusBar = self.CreateStatusBar(3)
        self.peakStatusBar.SetStatusWidths([-3, -2, -1])
        self._set_peak_status("Ready", 0, 0)

        self.detectionPanel = wx.Panel(self)
        detection = wx.BoxSizer(wx.HORIZONTAL)
        detection.Add(wx.StaticText(self.detectionPanel, label="Peak detection"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        detection.AddStretchSpacer(1)
        self.progress = wx.Gauge(self.detectionPanel, range=100, size=(130, 16))
        detection.Add(self.progress, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        self.findButton = wx.Button(self.detectionPanel, label="Find Peaks", size=(-1, 25))
        self.findButton.Bind(wx.EVT_BUTTON, self.start_peak_search)
        detection.Add(self.findButton, 0, wx.RIGHT, 5)
        self.stopButton = wx.Button(self.detectionPanel, label="Stop", size=(-1, 25))
        self.stopButton.Bind(wx.EVT_BUTTON, self.stop_peak_search)
        self.stopButton.Disable()
        detection.Add(self.stopButton, 0, wx.RIGHT, 5)
        detection.Add(wx.StaticText(self.detectionPanel, label="Peaks to fit:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 5)
        self.peakCountCtrl = wx.SpinCtrl(self.detectionPanel, min=1, max=100,
                                        initial=int(self.picker_settings.max_peaks), size=(58, -1),
                                        style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.peakCountCtrl.SetToolTip("Maximum number of representative isolated peaks used for peak-shape fitting")
        self.peakCountCtrl.Bind(wx.EVT_SPINCTRL, self._on_peak_count_changed)
        self.peakCountCtrl.Bind(wx.EVT_TEXT_ENTER, self._on_peak_count_entered)
        detection.Add(self.peakCountCtrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 7)
        self.tuneButton = wx.Button(self.detectionPanel, label="Fine tuning...", size=(-1, 25))
        self.tuneButton.Bind(wx.EVT_BUTTON, self.open_peak_tuning)
        detection.Add(self.tuneButton, 0)
        self.detectionPanel.SetSizer(detection)

        self.fig = Figure(figsize=(6.8, 1.8))
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(680,180))

        # Diagnostic population used to infer the deconvolution filter shape.
        # Give this enough vertical room for the distribution to be genuinely
        # useful while keeping exactly the same horizontal plot geometry as
        # the peak-trace figure above it.
        self.fig_widths = Figure(figsize=(6.8, 2.25))
        self.canvas_widths = FigCanvas(self, -1, self.fig_widths)
        self.canvas_widths.SetMinSize(wx.Size(680, 220))
        self._draw_width_distribution(None)

        self.fig_sliders = Figure(figsize=(6.8, 1.7))
        self.canvas_sliders = FigCanvas(self, -1, self.fig_sliders)
        self.canvas_sliders.SetMinSize(wx.Size(680,170))

        self.previewPanel = None
        if False:  # radius preview now lives in FitRadiusFrame
            self.previewPanel = wx.Panel(self)
            preview_sizer = wx.BoxSizer(wx.VERTICAL)
            radius_sizer = wx.BoxSizer(wx.HORIZONTAL)
            radius_sizer.Add(wx.StaticText(self.previewPanel, label='3D peak preview'), 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
            radius_sizer.Add(wx.StaticText(self.previewPanel, label='RadiusF1:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
            self.radiusF1Box = wx.TextCtrl(self.previewPanel, value='%.6g' % self.radius_f1, size=(70, -1), style=wx.TE_PROCESS_ENTER)
            self.radiusF1Box.Bind(wx.EVT_TEXT_ENTER, self._radius_entered)
            radius_sizer.Add(self.radiusF1Box, 0, wx.RIGHT, 10)
            if self.dim == 2:
                radius_sizer.Add(wx.StaticText(self.previewPanel, label='RadiusF2:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
                self.radiusF2Box = wx.TextCtrl(self.previewPanel, value='%.6g' % self.radius_f2, size=(70, -1), style=wx.TE_PROCESS_ENTER)
                self.radiusF2Box.Bind(wx.EVT_TEXT_ENTER, self._radius_entered)
                radius_sizer.Add(self.radiusF2Box, 0)
            preview_sizer.Add(radius_sizer, 0, wx.EXPAND | wx.BOTTOM, 3)
            self.fig_3d = Figure(figsize=(6.8, 2.4))
            self.canvas_3d = FigCanvas(self.previewPanel, -1, self.fig_3d)
            self.canvas_3d.SetMinSize(wx.Size(680, 235))
            preview_sizer.Add(self.canvas_3d, 0, wx.EXPAND)
            self.previewPanel.SetSizer(preview_sizer)
        # self.fig_sliders.patch.set_facecolor('xkcd:salmon')

        # Keep native widgets on a single background-coloured panel.  This
        # hides platform-specific edge gaps and makes the bottom control strip
        # visually continuous with the Peak Fit window.
        self.buttonPanel = wx.Panel(self)
        self.buttonPanel.SetBackgroundColour(self.GetBackgroundColour())
        self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self._link_syncing = False
        self.rescaleButton = wx.Button(self.buttonPanel, label="Rescale sliders", size=(-1,25))
        self.rescaleButton.Bind(wx.EVT_BUTTON, self.onButtonRescale)

        self.fitButton = wx.Button(self.buttonPanel, label="Fit", size=(80,25))
        self.fitButton.SetToolTip("Estimate a conservative deconvolution filter from the selected representative peaks")
        self.fitButton.Bind(wx.EVT_BUTTON, self.onButtonFit)
        self.linkWidths = wx.CheckBox(self.buttonPanel, label="Link Gaussian/Lorentzian widths")
        self.linkWidths.SetValue(self._saved_link_widths)
        self.linkWidths.SetToolTip("Keep Gaussian and Lorentzian FWHM equal, matching the restricted 1D peak-fitting model")
        self.linkWidths.Bind(wx.EVT_CHECKBOX, self.onLinkWidths)

        self.CloseButton = wx.Button(self.buttonPanel, label="Close", size=(80, 25))
        self.CloseButton.Bind(wx.EVT_BUTTON, self.onButtonClose)
        flags = wx.RIGHT|wx.TOP|wx.BOTTOM
        self.button_sizer.Add(self.rescaleButton, border=10, flag=flags)
        self.button_sizer.Add(self.fitButton, border=10, flag=flags)
        self.button_sizer.Add(self.linkWidths, border=10, flag=flags | wx.ALIGN_CENTER_VERTICAL)

        self.button_sizer.Add(self.CloseButton, border=10, flag=flags)
        self.buttonPanel.SetSizer(self.button_sizer)

        # if self.dim == 3:
            # self.fit3D = wx.Button(self, label="3D Fit", size=(80,25))
            # self.fit3D.Bind(wx.EVT_BUTTON, self.onButton3D)
            # self.button_sizer.Add(self.fit3D, border=10, flag=flags)




        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.Add(self.detectionPanel, 0, flag=wx.GROW | wx.ALL, border=6)
        self.main_sizer.Add(self.canvas, 0, flag=wx.GROW)
        self.main_sizer.Add(self.canvas_widths, 0, flag=wx.GROW)
        self.above = wx.BoxSizer(wx.HORIZONTAL)

        self.above.AddSpacer(50)
        self.aboveLine = wx.Panel(self, -1, size=(100,1))
        self.aboveLine.SetBackgroundColour(wx.Colour(192, 192, 192))
        self.above.Add(self.aboveLine,200, flag=wx.GROW)
        self.above.AddSpacer(50)
        self.main_sizer.Add(self.above, flag=wx.GROW)
        self.main_sizer.Add(self.canvas_sliders, 0, flag=wx.GROW)
        if self.previewPanel is not None:
            self.main_sizer.Add(self.previewPanel, 0, flag=wx.GROW | wx.LEFT | wx.RIGHT, border=6)
        self.main_sizer.Add(self.buttonPanel, 0, flag=wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT, border=6)

        self.SetSizer(self.main_sizer)
        self.main_sizer.Fit(self)

    

    def _set_linked_width(self, dim, value, source):
        if not getattr(self, 'linkWidths', None) or not self.linkWidths.GetValue() or self._link_syncing:
            return
        if dim >= len(getattr(self, 'psf_sliders', [])):
            return
        target = self.lorentz_sliders[dim] if source == 'gaussian' else self.psf_sliders[dim]
        if abs(float(target.val) - float(value)) <= max(1e-12, abs(float(value))*1e-9):
            return
        self._link_syncing = True
        try:
            target.set_val(float(value))
        finally:
            self._link_syncing = False

    def onLinkWidths(self, event=None):
        if self.linkWidths.GetValue() and hasattr(self, 'psf_sliders'):
            self._link_syncing = True
            try:
                for dim, slider in enumerate(self.psf_sliders):
                    self.lorentz_sliders[dim].set_val(float(slider.val))
            finally:
                self._link_syncing = False
            self.canvas_sliders.draw()
            self.canvas.draw_idle()

    def onButtonFit(self, event=None):
        if not len(self.maxima):
            wx.MessageBox("Find representative peaks before fitting the filter shape.", "Fit peaks", wx.OK | wx.ICON_INFORMATION)
            return
        try:
            estimates = estimate_filter_shape(self.data, self.indexes, self.maxima,
                                              link_widths=self.linkWidths.GetValue())
        except Exception as exc:
            wx.MessageBox("Could not estimate a robust filter shape:\n%s" % exc, "Fit peaks", wx.OK | wx.ICON_ERROR)
            return
        self._link_syncing = True
        try:
            for dim, estimate in enumerate(estimates):
                # Expand slider bounds before set_val so the recommendation remains visible.
                width = float(estimate.recommended_fwhm)
                for slider in (self.psf_sliders[dim], self.lorentz_sliders[dim]):
                    slider.valmin = min(slider.valmin, width / 2.0)
                    slider.valmax = max(slider.valmax, width * 2.0)
                    slider.ax.set_xlim(slider.valmin, slider.valmax)
                self.psf_sliders[dim].set_val(width)
                self.lorentz_sliders[dim].set_val(width)
                self.voigt_sliders[dim].set_val(float(estimate.voigt_fraction))
        finally:
            self._link_syncing = False
        clean = sum(e.clean_sides for e in estimates)
        rejected = sum(e.rejected_sides for e in estimates)
        detail = "; ".join("F%d %.4g -> %.4g ppm" % (d+1, e.measured_fwhm, e.recommended_fwhm)
                           for d, e in enumerate(estimates))
        self._set_peak_status("Robust envelope fit: %d clean wings, %d rejected; %s" %
                              (clean, rejected, detail), selected=len(self.maxima))
        # Slider callbacks update the filter curves; force both canvases to
        # repaint after the complete fit so the recommendation is immediately
        # visible even when wxAgg has coalesced intermediate draw_idle calls.
        self.canvas_sliders.draw()
        self.canvas.draw()
        self.canvas_widths.draw_idle()

    def onButtonRescale(self, event):
        for s in range(len(self.psf_sliders)):
            slid = self.psf_sliders[s]
            slid_lor = self.lorentz_sliders[s]
            slid_voigt = self.voigt_sliders[s]

            

            x = slid.val
            n = slid_voigt.val

            r = slid_lor.val

            x_axis = self.indexes[s]

            blah, fwhm = voigt(x_axis.astype(float), x_axis[0], x, r, n, fwhm=True)
            print(fwhm)
            # slid.eventson = False
            slid.valmin = fwhm/20.
            slid.valmax = fwhm*5.
            slid_lor.valmin = fwhm/20.
            slid_lor.valmax = fwhm*5.
            slid.ax.set_xlim(slid.valmin, slid.valmax)
            slid_lor.ax.set_xlim(slid.valmin, slid.valmax)
            # slid.eventson = True

        # for slid in self.lorentz_sliders:
        #     slid.valmin = slid.val/20.
        #     slid.valmax = slid.val*5.
        #     # slid.eventson = False
        #     slid.ax.set_xlim(slid.valmin, slid.valmax)
        #     # slid.eventson = True

        self.canvas_sliders.draw()
        # self.psf_sliders[0].valmin
        

    def _snapshot_fit_values(self):
        # Source-contract compatibility: getattr(self.tabOne, 'sig%dBox' % which)

        """Return the current persistable Fit Peaks state.

        During __init__ the Matplotlib sliders do not exist until the first
        peak search has populated the plots.  Use the UniDecNMR controls as
        the opening-state source in that interval, then transparently switch
        to the live sliders once they have been created.
        """
        values=[]
        psf = getattr(self, 'psf_sliders', None)
        voigt = getattr(self, 'voigt_sliders', None)
        lorentz = getattr(self, 'lorentz_sliders', None)
        sliders_ready = (psf is not None and voigt is not None and lorentz is not None and
                         len(psf) >= self.dim and len(voigt) >= self.dim and
                         len(lorentz) >= self.dim)
        for dim in range(self.dim):
            if sliders_ready:
                values.extend((float(psf[dim].val), float(voigt[dim].val),
                               float(lorentz[dim].val)))
            else:
                shape = self.peak_fit_service.shape_parameters(self.dim)
                values.extend((shape['sigmas'][dim], shape['voigt'][dim], shape['lorentz'][dim]))
        peak_count = int(self.peakCountCtrl.GetValue()) if hasattr(self, 'peakCountCtrl') else int(self.picker_settings.max_peaks)
        link_widths = bool(self.linkWidths.GetValue()) if hasattr(self, 'linkWidths') else bool(self._saved_link_widths)
        values.extend((peak_count, link_widths))
        return tuple(values)

    def _fit_has_unsaved_changes(self):
        current=self._snapshot_fit_values()
        for current_value, opening_value in zip(current, self._opening_fit_values):
            if isinstance(opening_value, bool):
                if current_value != opening_value: return True
            elif isinstance(opening_value, int):
                if current_value != opening_value: return True
            elif abs(current_value-opening_value) > max(1e-12, abs(opening_value)*1e-9):
                return True
        return False

    def _save_fit_changes(self):
        sigmas=[float(slider.val) for slider in self.psf_sliders]
        voigt=[float(slider.val) for slider in self.voigt_sliders]
        lorentz=[float(slider.val) for slider in self.lorentz_sliders]
        self.peak_fit_service.set_shape_parameters(sigmas, voigt, lorentz)
        self.peak_fit_service.sync_usta_shape(sigmas[0], voigt[0], lorentz[0])
        self._save_fit_radii()
        self.peak_fit_service.save_fit_preferences(int(self.peakCountCtrl.GetValue()), bool(self.linkWidths.GetValue()))
        self.peak_fit_service.save_project()
        self.peak_fit_service.mark_peak_shape_determined(self.dim, sigmas, voigt, lorentz)
        self._opening_fit_values = self._snapshot_fit_values()
        window=getattr(self,'radius_window',None)
        if window is not None:
            window._opening_values = window._snapshot_values()

    def _on_frame_close(self, event):
        radius=getattr(self,'radius_window',None)
        radius_changed=bool(radius is not None and radius._has_unsaved_changes())
        fit_changed=self._fit_has_unsaved_changes()
        if fit_changed or radius_changed:
            answer=wx.MessageBox('Peak-fit parameters have changed. Save changes?', 'Fit peak list',
                                 wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION, self)
            if answer == wx.CANCEL:
                if event.CanVeto(): event.Veto()
                return
            if answer == wx.YES:
                if radius_changed:
                    radius.on_radius_entered()
                self._save_fit_changes()
        self._closing_pair=True
        if radius is not None:
            try: radius.Destroy()
            except RuntimeError: pass
            self.radius_window=None
        event.Skip()

    def onButtonClose(self, event):
        self.Close()


    def onButton3D(self, event):
        if self.dim ==3:
            from spinDecon.analysis import peak_shape_optimizer
            import importlib
            peak_shape_optimizer = importlib.reload(peak_shape_optimizer)
            self.starting = numpy.array([float(self.voigt_sliders[0].val), float(self.voigt_sliders[1].val),
                                         float(self.voigt_sliders[2].val), float(self.psf_sliders[0].val),
                                         self.psf_sliders[1].val, self.psf_sliders[2].val,
                                         self.lorentz_sliders[0].val, self.lorentz_sliders[1].val,
                                         self.lorentz_sliders[2].val, 1.]).astype(float)

            opt = peak_shape_optimizer.PeakShapeOptimizer(self.data, self.mesh_xx, self.mesh_yy, self.mesh_zz,
                                                        self.indexes, self.maxima,
                                                        self.starting)

            self.values = opt.x0
            self.draw_figure()

    
    def draw_figure(self):

        
        

        self.fig.clear()

        max_value = self.max_vals[0]
        psf_plots = []
        self.psf_sliders = []
        self.voigt_sliders = []
        self.lorentz_sliders = []
        self.offset_full = []
        self._slider_value_texts = []

        # if self.tabOne.pseudo == True:
        #     peak_dims = range(1, self.tabOne.dim)
        # else:
        #     peak_dims = range(self.tabOne.dim)
        for dim in range(self.dim):
            # Use only the active spectral dimensions.  The former 1x4 grid
            # left unused columns in 1D-3D and made these plots substantially
            # narrower than the width histograms.
            self.axes = self.fig.add_subplot(1, self.dim, dim + 1)

            spec_res = numpy.abs(self.indexes[dim][1]-self.indexes[dim][0])
            i = 0
            offsets = []
            for max in self.maxima:

                if self.dim == 4:
                    x,y,z,a = max

                    x2 = self.indexes[0][x]
                    y2 = self.indexes[1][y]
                    z2 = self.indexes[2][z]
                    a2 = self.indexes[3][a]
                    max2 = [x2,y2,z2,a2]
                    centre = self.indexes[dim][int(self.data.shape[dim]/2.)]
                    offset = centre-max2[dim]
                    height = max_value/self.max_vals[i]
                    multiple = 3.




                    x_axis = self.indexes[dim]
                    if dim ==0:
                        y_axis = self.data[:,y,z,a]*height
                        fwhm = self.data[:,y,z,a]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:
                            if (change_high +x > self.data.shape[0]-2):
                                break

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x-change_low,y,z,a]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x+change_high,y,z,a]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))

                    if dim ==1:
                        y_axis = self.data[x,:,z,a]*height
                        fwhm = self.data[x,:,z,a]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:
                            if (change_high +x > self.data.shape[0]-2):
                                break

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y-change_low,z,a]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y+change_high,z,a]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))

                    if dim ==2:
                        y_axis = self.data[x,y,:,a]*height
                        fwhm = self.data[x,y,:,a]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:
                            if (change_high +x > self.data.shape[0]-2):
                                break

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y,z-change_low,a]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y,z+change_high,a]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))
                    if dim ==3:
                        y_axis = self.data[x,y,z,:]*height
                        fwhm = self.data[x,y,z,:]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:
                            if (change_high +x > self.data.shape[0]-2):
                                break

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y,z,a-change_low]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y,z,a+change_high]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))

                    offsets.append([self.indexes[dim][change_low], self.indexes[dim][change_high]])

                if self.dim == 3:
                    x,y,z = max

                    x2 = self.indexes[0][x]
                    y2 = self.indexes[1][y]
                    z2 = self.indexes[2][z]
                    max2 = [x2,y2,z2]
                    centre = self.indexes[dim][int(self.data.shape[dim]/2.)]

                    offset = centre-max2[dim]
                    height = max_value/self.max_vals[i]



                    x_axis = self.indexes[dim]



                    multiple = 3.

                    if dim ==0:
                        y_axis = self.data[:,y,z]*height

                        fwhm = self.data[:,y,z]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:
                            if (change_high +x > self.data.shape[0]-2):
                                break

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x-change_low,y,z]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x+change_high,y,z]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))


                    if dim ==1:
                        y_axis = self.data[x,:,z]*height

                        fwhm = self.data[x,:,z]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y-change_low,z]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y-change_high,z]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))
                        #self.axes.set_xlim(self.indexes[1][int(self.data.shape[dim]/2.)-change_low], self.indexes[1][int(self.data.shape[dim]/2.)+change_low])

                    if dim ==2:
                        y_axis = self.data[x,y,:]*height

                        fwhm = self.data[x,y,:]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y,z-change_low]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y,z+change_low]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                        change_low = numpy.amax([0, int(self.data.shape[dim]/2.)-change_low])
                        change_high = min(int(self.data.shape[dim]-1), int(numpy.floor(self.data.shape[dim]/2.)+change_high))
                    # offsets.append([self.indexes[dim][int(self.data.shape[dim]/2.)-change_low], self.indexes[dim][int(self.data.shape[dim]/2.)+change_high]])
                    offsets.append([self.indexes[dim][change_low], self.indexes[dim][change_high]])
                        #self.axes.set_xlim(self.indexes[2][int(self.data.shape[dim]/2.)-change_low], self.indexes[2][int(self.data.shape[dim]/2.)+change_low])


                if self.dim == 2:
                    x,y = max

                    x2 = self.indexes[0][x]
                    y2 = self.indexes[1][y]
                    max2 = [x2,y2]
                    centre = self.indexes[dim][int(self.data.shape[dim]/2.)]
                    offset = centre-max2[dim]
                    height = max_value/self.max_vals[i]
                    multiple = 5.

                    x_axis = self.indexes[dim]
                    if dim ==0:
                        y_axis = self.data[:,y]*height

                        fwhm = self.data[:,y]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x-change_low,y]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x+change_high,y]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)


                    if dim ==1:
                        y_axis = self.data[x,:]*height

                        fwhm = self.data[x,:]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x,y-change_low]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x,y-change_high]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)


                    offsets.append([self.indexes[dim][int(self.data.shape[dim]/2.)-change_low], self.indexes[dim][int(self.data.shape[dim]/2.)+change_low]])


                if self.dim == 1:

                    x = max

                    x2 = self.indexes[0][x]
                    max2 = [x2]
                    centre = self.indexes[dim][int(self.data.shape[dim]/2.)]
                    offset = numpy.array(centre-max2[dim])
                    height = max_value/self.max_vals[i]

                    x_axis = numpy.array(self.indexes[dim])
                    multiple = 5.
                    if dim ==0:
                        y_axis = self.data[:]*height

                        fwhm = self.data[:]-(max_value/2.)
                        change_low=0
                        change_high=0
                        resulting_number_low = 1
                        resulting_number_high = 1
                        while resulting_number_low>0 and resulting_number_high>0:

                            if resulting_number_low >0:
                                change_low+=1
                                resulting_number_low = self.data[x-change_low]-max_value/2.
                            if resulting_number_high >0:
                                change_high+=1
                                resulting_number_high = self.data[x+change_high]-max_value/2.

                        change_low = int(change_low*multiple)
                        change_high = int(change_high*multiple)
                    offsets.append([self.indexes[dim][int(self.data.shape[dim]/2.)-change_low], self.indexes[dim][int(self.data.shape[dim]/2.)+change_high]])
                self.axes.plot(numpy.squeeze(x_axis+offset), numpy.squeeze(y_axis), ls='--', linewidth=0.5)
                self.axes.set_title(self.labb[dim], fontsize=8, pad=2)

                i += 1
            if dim > 0:
                self.axes.set_yticks([])
            self.axes.set_xticks([])
            offsets = numpy.array(offsets)
            self.axes.set_xlim(numpy.max(offsets[:,0]), numpy.min(offsets[:,1]))
            # n = self.starting[dim]
            # s = self.starting[dim + 3]
            # r = self.starting[dim + 6]
            # psf_plots.append(self.axes.plot(x_axis, voigt(x_axis.astype(float), centre, s, r, n) * self.max_vals[-1]))

            if self.dim == 4:
                n = numpy.abs(self.values[dim])
                s = numpy.abs(self.values[dim+4])
                r = numpy.abs(self.values[dim+8])
            elif self.dim == 3:
                n = numpy.abs(self.values[dim])
                s = numpy.abs(self.values[dim+3])
                r = numpy.abs(self.values[dim+6])
            elif self.dim ==2:
                n = numpy.abs(self.values[dim])
                s = numpy.abs(self.values[dim+2])
                r = numpy.abs(self.values[dim+4])
            elif self.dim ==1:
                n = numpy.abs(self.values[dim])
                s = numpy.abs(self.values[dim+1])
                r = numpy.abs(self.values[dim+2])
            # print x_axis, centre, s_init

            s_min = numpy.abs(self.indexes[dim][1]-self.indexes[dim][0])/3
            s_min = s/20.
            s_max = s*5.
            #s_max = numpy.abs(self.indexes[dim][int(len(self.indexes[dim])/16)]-self.indexes[dim][0])
            s_step = s_max/500.
            psf_plots.append(self.axes.plot(x_axis, voigt(x_axis.astype(float), centre, s, r, n)*self.max_vals[0]))
            #self.axes.set_xlim(centre-10*self.sigs[dim], centre+10*self.sigs[dim])
            self.axes.set_ylim(0, max_value*1.1)

            slider_axes = self.fig_sliders.add_subplot(3, self.dim, dim + 1)

            if dim == 0:
                self.psf_sliders.append(Slider(slider_axes,  # the axes object containing the slider
                              "",  # row label is drawn once at the figure edge
                              s_min,  # minimal value of the parameter
                              s_max,  # maximal value of the parameter
                              valinit=s,  # initial value of the parameter
                              valstep=s_step, valfmt="%.2g", color='gray'))
            else:
                self.psf_sliders.append(Slider(slider_axes,  # the axes object containing the slider
                              "",  # the name of the slider parameter
                              s_min,  # minimal value of the parameter
                              s_max,  # maximal value of the parameter
                              valinit=s,  # initial value of the parameter
                              valstep=s_step, valfmt="%.2g", color='gray'))

            self.voigt_slider_axes = self.fig_sliders.add_subplot(3, self.dim, dim + 1 + self.dim)

            if dim == 0:
                self.voigt_sliders.append(Slider(self.voigt_slider_axes, "", 0, 1, valinit=n, valstep=0.01, valfmt="%.2g", color='gray'))
            else:
                self.voigt_sliders.append(Slider(self.voigt_slider_axes, "", 0, 1, valinit=n, valstep=0.01, valfmt="%.2g", color='gray'))

            lorentz_axes = self.fig_sliders.add_subplot(3, self.dim, dim + 1 + 2 * self.dim)

            if dim ==0:
                self.lorentz_sliders.append(Slider(lorentz_axes,  # the axes object containing the slider
                                           "", # row label is drawn once at the figure edge
                                           s_min,  # minimal value of the parameter
                                           s_max,  # maximal value of the parameter
                                           valinit=r,  # initial value of the parameter
                                           valstep=s_step, valfmt="%.2g", color='gray'))
            else:
                self.lorentz_sliders.append(Slider(lorentz_axes,  # the axes object containing the slider
                                           "", #self.tabOne.labb[dim],  # the name of the slider parameter
                                           s_min,  # minimal value of the parameter
                                           s_max,  # maximal value of the parameter
                                           valinit=r,  # initial value of the parameter
                                           valstep=s_step, valfmt="%.2g", color='gray'))

            # Matplotlib Slider's native value text sits outside the axes and can
            # collide with the next dimension.  It can also leave stale glyphs
            # with backend blitting.  Keep values inside their own column and
            # redraw this small canvas explicitly when a slider changes.
            for _slider in (self.psf_sliders[-1], self.voigt_sliders[-1], self.lorentz_sliders[-1]):
                _slider.valtext.set_visible(False)
                _value_text = _slider.ax.text(
                    0.985, 0.50, "%.3g" % float(_slider.val),
                    transform=_slider.ax.transAxes, ha='right', va='center', fontsize=7,
                    bbox=dict(facecolor=self.fig_sliders.get_facecolor(), edgecolor='none', pad=0.6))
                self._slider_value_texts.append((_slider, _value_text))

            # break
            # psf_sliders.append(current_slider)
        # slider_axes = self.fig.add_subplot(20, 2, 26 + (self.tabOne.dim * 2))
        # self.psf_sliders.append(Slider(slider_axes, "Voigt", 0, 1, valinit=float(self.tabOne.voigtBox.GetValue()), valstep=0.01, valfmt="%.2g"))



        # Match the horizontal margins/column spacing of the width-distribution
        # figure so traces and histograms line up dimension-for-dimension.
        self.fig.subplots_adjust(left=0.07, right=0.985, bottom=0.10, top=0.88, wspace=0.24)
        # Keep row names independent of Matplotlib Slider labels.  This gives
        # them a fixed, left-aligned margin and prevents long labels from being
        # clipped or stealing horizontal space from the first slider column.
        for _label, _y in (("Gaussian", 0.835), ("Voigt", 0.505), ("Lorentzian", 0.175)):
            self.fig_sliders.text(0.012, _y, _label, ha='left', va='center', fontsize=7)
        self.fig_sliders.subplots_adjust(left=0.105, right=0.985, bottom=0.08, top=0.94, wspace=0.24, hspace=0.72)

        def _refresh_slider_values():
            for _slider, _text in self._slider_value_texts:
                _text.set_text("%.3g" % float(_slider.val))
            # wxAgg can coalesce draw_idle requests while a Slider is dragging,
            # leaving successive text glyphs composited on top of one another.
            # This canvas is deliberately small, so force a synchronous full
            # repaint; that clears the previous value before drawing the new one.
            self.canvas_sliders.draw()

        def slider_update0(s):
            self._set_linked_width(0, self.psf_sliders[0].val, 'gaussian')
            sigma = self.psf_sliders[0].val
            n = self.voigt_sliders[0].val
            r = self.lorentz_sliders[0].val

            new_ydata = voigt(self.indexes[0], self.indexes[0][int(self.data.shape[0]/2.)], sigma,r, n) * max_value
            psf_plots[0][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()
        def slider_update1(s):
            self._set_linked_width(1, self.psf_sliders[1].val, 'gaussian')
            sigma = self.psf_sliders[1].val
            n = self.voigt_sliders[1].val
            r = self.lorentz_sliders[1].val

            new_ydata = voigt(self.indexes[1], self.indexes[1][int(self.data.shape[1]/2.)], sigma,r, n) * max_value
            psf_plots[1][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()
        def slider_update2(s):
            self._set_linked_width(2, self.psf_sliders[2].val, 'gaussian')
            sigma = self.psf_sliders[2].val
            n = self.voigt_sliders[2].val
            r = self.lorentz_sliders[2].val

            new_ydata = voigt(self.indexes[2], self.indexes[2][int(self.data.shape[2]/2.)], sigma,r, n) * max_value
            psf_plots[2][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()
        def slider_update3(s):
            self._set_linked_width(3, self.psf_sliders[3].val, 'gaussian')
            sigma = self.psf_sliders[3].val
            n = self.voigt_sliders[3].val
            r = self.lorentz_sliders[3].val

            new_ydata = voigt(self.indexes[3], self.indexes[3][int(self.data.shape[3]/2.)], sigma,r, n) * max_value
            psf_plots[3][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()
        def voigt_slider0(n):
            n = self.voigt_sliders[0].val
            sigma = self.psf_sliders[0].val
            r = self.lorentz_sliders[0].val
            new_ydata = voigt(self.indexes[0], self.indexes[0][int(self.data.shape[0]/2.)], sigma,r, n) * max_value
            psf_plots[0][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def voigt_slider1(n):
            n = self.voigt_sliders[1].val
            sigma = self.psf_sliders[1].val
            r = self.lorentz_sliders[1].val
            new_ydata = voigt(self.indexes[1], self.indexes[1][int(self.data.shape[1]/2.)], sigma,r, n) * max_value
            psf_plots[1][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def voigt_slider2(n):
            n = self.voigt_sliders[2].val
            sigma = self.psf_sliders[2].val
            r = self.lorentz_sliders[2].val
            new_ydata = voigt(self.indexes[2], self.indexes[2][int(self.data.shape[2]/2.)], sigma,r, n) * max_value
            psf_plots[2][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def voigt_slider3(n):
            n = self.voigt_sliders[3].val
            sigma = self.psf_sliders[3].val
            r = self.lorentz_sliders[3].val
            new_ydata = voigt(self.indexes[3], self.indexes[3][int(self.data.shape[3]/2.)], sigma,r, n) * max_value
            psf_plots[3][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def lorentz_slider0(n):
            self._set_linked_width(0, self.lorentz_sliders[0].val, 'lorentzian')
            n = self.voigt_sliders[0].val
            sigma = self.psf_sliders[0].val
            r = self.lorentz_sliders[0].val
            new_ydata = voigt(self.indexes[0], self.indexes[0][int(self.data.shape[0] / 2.)], sigma, r, n) * max_value
            psf_plots[0][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def lorentz_slider1(n):
            self._set_linked_width(1, self.lorentz_sliders[1].val, 'lorentzian')
            n = self.voigt_sliders[1].val
            sigma = self.psf_sliders[1].val
            r = self.lorentz_sliders[1].val
            new_ydata = voigt(self.indexes[1], self.indexes[1][int(self.data.shape[1] / 2.)], sigma, r, n) * max_value
            psf_plots[1][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def lorentz_slider2(n):
            self._set_linked_width(2, self.lorentz_sliders[2].val, 'lorentzian')
            n = self.voigt_sliders[2].val
            sigma = self.psf_sliders[2].val
            r = self.lorentz_sliders[2].val
            new_ydata = voigt(self.indexes[2], self.indexes[2][int(self.data.shape[2] / 2.)], sigma, r, n) * max_value
            psf_plots[2][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        def lorentz_slider3(n):
            self._set_linked_width(3, self.lorentz_sliders[3].val, 'lorentzian')
            n = self.voigt_sliders[3].val
            sigma = self.psf_sliders[3].val
            r = self.lorentz_sliders[3].val
            new_ydata = voigt(self.indexes[3], self.indexes[3][int(self.data.shape[3] / 2.)], sigma, r, n) * max_value
            psf_plots[3][0].set_ydata(new_ydata)
            _refresh_slider_values()
            self.canvas.draw_idle()

        # print psf_sliders

        self.draw_3d_peak_preview()

        self.psf_sliders[0].on_changed(slider_update0)
        self.voigt_sliders[0].on_changed(voigt_slider0)
        self.lorentz_sliders[0].on_changed(lorentz_slider0)
        if self.dim > 1:
            self.psf_sliders[1].on_changed(slider_update1)
            self.voigt_sliders[1].on_changed(voigt_slider1)
            self.lorentz_sliders[1].on_changed(lorentz_slider1)

        if self.dim >2:
            self.psf_sliders[2].on_changed(slider_update2)
            self.voigt_sliders[2].on_changed(voigt_slider2)
            self.lorentz_sliders[2].on_changed(lorentz_slider2)

        if self.dim>3:
            self.psf_sliders[3].on_changed(slider_update3)
            self.voigt_sliders[3].on_changed(voigt_slider3)
            self.lorentz_sliders[3].on_changed(lorentz_slider3)


        self.canvas.draw()

    def _find_maxima_compat(self, numbero):
        """Compatibility shim for callers of the former dimension-specific pickers."""
        settings = copy.copy(self.picker_settings)
        settings.max_peaks = int(numbero)
        result = PeakPicker(self.data, settings).run()
        return result.maxima, result.values

    def find_maxima_1D(self, numbero):
        return self._find_maxima_compat(numbero)

    def find_maxima_2D(self, numbero):
        return self._find_maxima_compat(numbero)

    def find_maxima_3D(self, numbero):
        return self._find_maxima_compat(numbero)

    def find_maxima_4D(self, numbero):
        return self._find_maxima_compat(numbero)
