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
from spinDecon.gui.context import context_for, project_for, data_for
from spinDecon.analysis.slice_service import SliceService
from spinDecon.domain.dimensions.viewer_contract import topology_for
import string
import os
import numpy
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
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

############################################################################
# Frame for 1d slices
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

_SUPERSCRIPT_TRANS = str.maketrans('-+0123456789', '⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹')

def _scientific_unicode(value):
    """Format an intensity as e.g. -1.23 x 10⁹ without MathText delimiters."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not numpy.isfinite(value):
        return str(value)
    if value == 0.0:
        return '0'
    exponent = int(numpy.floor(numpy.log10(abs(value))))
    mantissa = value / (10.0 ** exponent)
    mantissa_text = ('%.3g' % mantissa).replace('-', '−')
    return '%s × 10%s' % (mantissa_text, str(exponent).translate(_SUPERSCRIPT_TRANS))

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()

class SliceFrame(wx.Panel):
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
        self.slice_service = (getattr(self.app_context, "slices", None) if self.app_context is not None else None) or SliceService(tabOne)
        self.state = project_for(tabOne, parent)
        self.sum=(0.,2.)
        # Peak data remain in decon_tab/DataStore; this viewer resolves them on demand.
        self.thresh = self.slice_service.max_intensity
        self.offset=0

        dmin, dmax = self.slice_service.axis_limits(1)

        self.create_main_panel()
        self.draw_figure()
        #self.Show()
        self.Fit()

    @property
    def peak(self):
        """Current 2D reference peaks from the shared application store.

        SliceFrame is a viewer and must not own a private peak-list copy.
        Keeping this compatibility property lets the legacy drawing code use
        ``self.peak`` while resolving the authoritative reference list on each
        access.
        """
        return self.slice_service.reference_peaks()

    def drawing_box(self):
        # These controls are part of the toolbar row, so use a plain sizer.
        # An empty StaticBox reserves unnecessary native border/best-size space.
        self.vbox2 = wx.BoxSizer(wx.HORIZONTAL)

        listy=[]
        for i in range(len(self.peak)):
            if(self.peak[i].name not in listy):
                listy.append(self.peak[i].name)
        self.ComboBox1=wx.ComboBox(self, -1, size=(80, 22), choices=listy, style=wx.CB_READONLY)
        if listy:
            self.ComboBox1.SetSelection(0)



        # Compact reference navigator: arrows step the same ComboBox selection
        # previously controlled by the separate Previous/Next buttons.
        self.peakSpin = wx.SpinButton(self, size=(22,22), style=wx.SP_VERTICAL)
        self.peakSpin.SetRange(-99999, 99999)
        self.peakSpin.Bind(wx.EVT_SPIN_UP, self.on_peak_spin_up)
        self.peakSpin.Bind(wx.EVT_SPIN_DOWN, self.on_peak_spin_down)

        self.cb_grid = wx.CheckBox(self, -1,"Peaks",style=wx.ALIGN_RIGHT)
        # State-only compatibility control; the visible Peaks toggle is in the Matplotlib toolbar.
        self.cb_grid.Hide()
        self.cb_calc = _ToolbarToggleState(False)

        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)

        self.vbox2.Add(self.ComboBox1, border=3, flag=self.flags)
        self.vbox2.Add(self.peakSpin, border=3, flag=self.flags)

    def control_box(self):
        # Reference navigation now lives beside ComboBox1 as a compact
        # SpinButton; retain this hook because create_main_panel calls it.
        self.controlSizer = wx.BoxSizer(wx.HORIZONTAL)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        ## Initialising matplotlib plot
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('button_release_event', self.on_pick)
        self.axes = self.fig.add_subplot(111)
        # Keep the toolbar readout readable for NMR intensities: wx.StaticText
        # displays Unicode directly, so avoid Matplotlib's ${...}$ MathText form.
        self.axes.format_coord = lambda x, y: 'x=%.4f, y=%s' % (x, _scientific_unicode(y))
        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view, peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon, coordinates=False)

        ## Adding our control boxes
        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        self.drawing_box()
        self.control_box()

        ## Add application controls to Matplotlib's native wx.ToolBar so its
        ## wx.TB_BOTTOM border spans both the MPL tools and our controls.
        # Keep the complete Matplotlib tool group (Home/Back/Forward/Pan/Zoom/Save)
        # together, then append our controls, then the live coordinate readout.
        self.toolbar.AddSeparator()
        for widget in (self.ComboBox1, self.peakSpin):
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)
        self.toolbar.bind_control_status_help(self.ComboBox1, 'Select reference peak')
        self.toolbar.bind_control_status_help(self.peakSpin, 'Previous or next reference peak')
        self.toolbar.AddSeparator()
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        ## Main vertical sizer: the plot consumes all remaining space.
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 2)

        ## Final layout adjustments
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)


    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    """
    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=string.split(line)
            peak.append(numpy.array(linetosave).astype(numpy.float))
        peakfile.close()
        return numpy.array(peak)
    """

    def draw_figure(self):
        #try:
        self.draw_figureGO()
        #except:
        #    pass

    def draw_figureGO(self):
        """ Redraws the figure
        """
        self.axes.clear()
        self.axes.format_coord = lambda x, y: 'x=%.4f, y=%s' % (x, _scientific_unicode(y))
        # Compact plotting area: minimise unused margins and use an open frame.
        self.fig.subplots_adjust(left=0.075, right=0.985, bottom=0.105, top=0.985)
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.tick_params(top=False, right=False)

        sele1=self.ComboBox1.GetSelection()
        self.thresh = self.slice_service.threshold()

        if(self.spectral_dim_count==3):
            view = self.slice_service.reference_1d_view(sele1)
            if view is None:
                return
            xs=view['x']
            noise_sigma = view.get('noise_sigma')
            if noise_sigma is None or not numpy.isfinite(noise_sigma) or noise_sigma <= 0:
                self.axes.text(0.5, 0.5, 'S/N unavailable: no valid noise estimate',
                               transform=self.axes.transAxes, ha='center', va='center')
                return
            ys=numpy.asarray(view['raw']) / noise_sigma
            self.thresh=float(view['threshold']) / noise_sigma
            y2s=numpy.zeros_like(ys)
            y2s.fill(self.thresh)

            self.axes.set_xlabel(self.slice_service.label(0), fontsize=8)
            self.axes.set_ylabel('S/N', fontsize=8)
            self.axes.format_coord = lambda x, y: 'x=%.4f, S/N=%.3f' % (x, y)
            self.axes.plot(xs,ys,'r',label='Data', lw=1.0)
            self.axes.plot(xs,y2s,'g',label='Noise Threshold', ls='--', lw=0.9)
            self.xmin,self.xmax=xs[0],xs[-1]
            self.axes.set_xlim(self.xmin, self.xmax)
            self.ymin,self.ymax=self.axes.get_ylim()
            self.offset=-1*self.ymin/2


            #self.axes.text(float(self.peak[sele1].ppmI),-float(self.offset),self.peak[sele1].name,fontsize=8,rotation=90)

            #for i in range(len(self.peak)): #write in the peak labels
            #    self.axes.text(float(self.peak[i].ppmI),-float(self.offset),self.peak[i].name,fontsize=8,rotation=90)

            

            if(self.sum[1]==2):
                
                #plt.title(self.ComboBox1.GetValue()+' automated assignment')
                #if self.slice_service.decon_enabled:
                #if(os.path.exists('out/slice2d/'+self.ComboBox1.GetValue()+'.proj.decon')):

                if self.slice_service.decon_enabled:
                    
                    xs_f=view['x']

                    ys_f=view.get('decon')
                    if ys_f is None:
                        ys_f = numpy.zeros_like(ys)
                    else:
                        ys_f = numpy.asarray(ys_f) / noise_sigma
                    ys_d=ys_f-ys-2*float(self.offset)  #deifference
                    if self.cb_calc.IsChecked():
                        self.axes.plot(xs_f,ys_f,'b',label='Deconvolved', lw=1.0)
                        self.axes.plot(xs_f,ys_d,color='#ADD8E6',label='Difference', lw=1.0)
                    current = []
                    

                    
                    if self.cb_grid.IsChecked():
                        for marker in view.get('markers', []):
                            xpos = marker['x']
                            x=(xpos,xpos)
                            y=(0,marker['height'] / noise_sigma)
                            self.axes.plot(x,y,'k', label='Peak locations (%s)' % marker['label'], lw=1.0)
                            self.axes.text(xpos,-float(self.offset),marker['label'],fontsize=8,rotation=90)

                    #
                    #     distance = numpy.abs(numpy.linalg.norm([cn.f1-current[0], cn.f2-current[1]]))
                    #     print(distance)
                    #     if (distance < 5 and cn.v1 !=sele1):
                    #         # print sele1, cn.v1, cn.f1, cn.f2, current[0], current[1]
                    #         x = (cn.f3, cn.f3)
                    #         y = (0, cn.s1)
                    #         self.axes.plot(x, y, 'k', alpha=0.3, lw=1.0)
                    #         self.axes.text(cn.f3, -float(self.offset), cn.tag2, fontsize=8, rotation=90)



            self.axes.legend(fontsize=8)

            self.axes.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


        if(self.spectral_dim_count==4):
            #WILL NEED TO RETURN TO THIS

            input=self.readfile('out/slice2D/'+self.ComboBox1.GetValue()+'.dat.blura')
            print(len(input))
            xs=[];y2s=[]
            LAB=[];DB=[];DBA=[];DS=[];DI=[]
            for i in range(len(input)):
                if(len(input[i])>0):
                    if(numpy.abs(float(input[i][3]))>float(self.thresh)):
                        LAB.append(input[i][0]) #peaks
                    else:
                        LAB.append('')
                    DBA.append(float(input[i][3])-float(self.thresh)) #
                    DB.append(float(input[i][4])) #
                    DS.append(float(input[i][5])) #
                    DI.append(float(input[i][6])+float(self.thresh)) #
                    y2s.append(float(self.thresh))
            xs=numpy.arange(len(LAB))

            #self.axes.set_xlabel(self.slice_service.label(0), fontsize=8)
            self.axes.plot(xs,DBA,'r',label='DBA')
            self.axes.plot(xs,DB,'b',label='DB')
            self.axes.plot(xs,DS,'c',label='DS')
            self.axes.plot(xs,DI,'k',label='DI')
            self.axes.plot(xs,y2s,'y',label='threshold')

            self.LAB=LAB
            #self.axes.set_xticks(rotation='vertical')
            #matplotlib.pyplot.xticks(xs, LAB, rotation='vertical')
            import matplotlib.pyplot as plt
            self.axes.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))

            self.axes.legend(fontsize=8)

            for i in range(len(self.LAB)): #write in the peak labels
                if(numpy.fabs(DB[i])>self.thresh):
                    self.axes.text(xs[i],DB[i],LAB[i],fontsize=8,rotation=90)


            """
            labels = [item.get_text() for item in self.axes.get_xticklabels()]
            for i in range(len(labels)):
                labels[i]=LAB[i]

            #labels[1] = 'Testing'

            self.axes.set_xticklabels(labels)


            """
            #from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
            self.axes.yaxis.set_major_formatter(FormatStrFormatter('${%0.0e}$'))




        self.canvas.draw()

    def format_func(self,value, tick_number):
        # find number of multiples of pi/2
        logging.info(value,tick_number,int(value),len(self.LAB))

        if(numpy.fabs(int(value)-value)>0.01):
            return ''

        try:
            return self.LAB[int(value)]
        except:
            return ''

        """
        N = int(np.round(2 * value / np.pi))
        if N == 0:
            return "0"
        elif N == 1:
            return "$\\pi/2$"
        elif N == 2:
            return "$\\pi$"
        elif N % 2 > 0:
            return "${0}\\pi/2$".format(N)
        else:
            return "${0}\\pi$".format(N // 2)
        """


    def on_cb_grid(self, event):
        self.draw_figure()

    def on_slider_width(self, event):
        self.draw_figure()

    def _toolbar_decon(self, active):
        self.cb_calc.SetValue(bool(active))
        self.draw_figure()

    def _toolbar_peaks(self, active):
        self.cb_grid.SetValue(bool(active))
        self.on_cb_grid(None)

    def redraw_view(self):
        self.draw_figure()

    def on_draw_button(self, event):
        self.redraw_view()

    def on_peak_spin_up(self, event):
        self.on_N_button(event)

    def on_peak_spin_down(self, event):
        self.on_P_button(event)

    def on_P_button(self, event):
        if self.ComboBox1.GetSelection()<1:
            return
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()

    def on_N_button(self, event):
        if self.ComboBox1.GetSelection()>self.ComboBox1.GetCount()-2:
            return
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.draw_figure()


    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        #box_points = event.artist.get_bbox().get_points()

        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        dlg = wx.MessageDialog(
            self,
            msg,
            "Click!",
            wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def on_text_enter(self, event):
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
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER,
            self.on_flash_status_off,
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
