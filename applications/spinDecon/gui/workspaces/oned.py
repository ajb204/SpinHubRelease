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
from spinDecon.analysis.oned_service import OneDService
from spinDecon.domain.dimensions.viewer_contract import topology_for
import string
import os
import numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
import matplotlib.patches as patches
import copy
import nmrglue as ng
from scipy.interpolate import interp1d
############################################################################
# Frame for 1d slices
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()


class FileDrop(wx.FileDropTarget):

    def __init__(self, canvas,axis):

        wx.FileDropTarget.__init__(self)
        self.canvas = canvas
        self.axis = axis
        self.ucs= []
        self.data = []
        self.first_drop = True
        self.extra_plots = []
        # self.data_1d = data_1d

    def OnDropFiles(self, x, y, filenames):

        for name in filenames:
            if '.ft' in name:
                # try:
                    dic, data = ng.pipe.read(name)
                    if len(data.shape) == 1:
                        if self.first_drop:
                            msg = "Entering multiple plot mode: Please enter title of the first dataset"
                            dlg = wx.TextEntryDialog(None, msg)
                            res = dlg.ShowModal()
                            if res == wx.ID_CANCEL:
                                return False
                            self.data_1d.set_label(dlg.GetValue())
                            self.first_drop = False


                        uc= ng.pipe.make_uc(dic,data)
                        self.data.append(data)
                        self.ucs.append(uc)
                        msg = "Please enter title of this data!"
                        dlg = wx.TextEntryDialog(None, msg)
                        res = dlg.ShowModal()
                        if res == wx.ID_CANCEL:
                            self.canvas.draw_idle()
                            return False

                        self.extra_plots.append(self.axis.plot(uc.ppm_scale(), data, label = dlg.GetValue(), linewidth=0.5, picker = 5)[0])
                        self.axis.legend()
                        self.canvas.draw_idle()
                    else:
                        msg = "This is not 1D data - currently more dims are not supported..."
                        dlg = wx.MessageDialog(None, msg)
                        dlg.ShowModal()

                        return False

                # except:
                #
                #     msg = "There was an error reading in the nmrFile - is it pipe format?"
                #     dlg = wx.MessageDialog(None, msg)
                #     dlg.ShowModal()
                #
                #     return False
            else:

                msg = "Can only deal with *.ft* files!"
                dlg = wx.MessageDialog(None, msg)
                dlg.ShowModal()

                return False



        return True

class OneDFrame(wx.Panel):
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
        self.state = project_for(tabOne, parent)
        self.store = data_for(tabOne, parent)
        self.sum=(0.,2.)
        self.one_d = (self.app_context.one_d if self.app_context is not None else None) or OneDService(tabOne)
        self.peak = self.one_d.peaks
        self.thresh = self.one_d.threshold()
        self.offset=0
        self.addition=False
        self.rectangles = []

        dmin, dmax = self.one_d.axis_limits

        self.create_main_panel()
        self.draw_figure()
        #self.Show()
        self.Fit()

    def drawing_box(self):
        """Create the borderless 1D controls used in the Matplotlib toolbar."""
        self.drawbutton = wx.Button(self, -1, "Draw!", size=(-1,22))
        self.cb_grid = wx.CheckBox(self, -1, "Peaks", style=wx.ALIGN_RIGHT)
        self.cb_calc = wx.CheckBox(self, -1, "ShowCalc", style=wx.ALIGN_RIGHT)
        self.integrate_button = wx.ToggleButton(self, -1, "Integrate", size=(-1,22))
        self.addition_button = wx.Button(self, -1, "Addition", size=(-1,22))

        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_calc)
        self.integrate_button.Bind(wx.EVT_TOGGLEBUTTON, self.on_integrate)
        self.addition_button.Bind(wx.EVT_BUTTON, self.on_addition_button)

    def on_addition_button(self, event):
        final_data = []
        uc=[]
        for x in self.dt.extra_plots:
            # factor = x.selected_plots.get_label().split('*')[-1]

            if len(uc) > 0:
                # print(uc[0][0]-uc[0][1])
                # upscale = (uc[0][0]-uc[0][1])/(x.get_xdata()[0]-x.get_xdata()[1])
                # print(x.get_xdata()[0]-x.get_xdata()[1])
                # print(upscale)
                measured_time = x.get_xdata()
                y_data = x.get_ydata()
                print((measured_time[0]-measured_time[1])-(uc[0]-uc[1]))
                if numpy.abs((measured_time[0]-measured_time[1])-(uc[0]-uc[1]))>0.0001:
                    cubic_interp = interp1d(measured_time, y_data, kind='cubic')
                    final_data.append(cubic_interp(uc))
                else:
                    final_data.append(y_data)

                # exit
            else:
                final_data.append(x.get_ydata())
                uc = x.get_xdata()

        final_data = numpy.array(final_data)
        final = numpy.sum(final_data, axis=0)
        if self.addition == False:
            self.addition_plot, = self.axes.plot(uc, final, label='addition', linewidth=0.5, picker = 5)
            self.addition=True
        else:
            self.addition_plot.set_ydata(final)
        self.canvas.draw_idle()

    # def calc_chi_square(self):


    def on_integrate(self, event):
        self.integrating = True
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        print('integrating')

    def on_mouse_move(self,event):
        if self.not_yet_drawn == True:

            self.background = self.canvas.copy_from_bbox(self.axes.bbox)

            self.current_h = 0
            self.v_line = self.axes.axvline(self.one_d.index[0], color = 'r', linewidth=2)
            # print(self.one_d.index.shape, self.one_d.data[0,:].shape)
            # self.h_line, = self.axes_h.plot(self.one_d.index, self.one_d.data[0,:], color='k', linewidth = 0.5, zorder=1000)
            # self.h_line, = self.axes.plot(self.one_d.index, numpy.zeros_like(self.one_d.index), color='k', linewidth = 0.5)
            # self.h_line.set_visible(False)
            # self.axes_h.set_ylim(numpy.min(self.one_d.data), numpy.max(self.one_d.data))


            self.canvas.draw()
            self.not_yet_drawn = False
        # if event.inaxes == None:
            # self.h_line.set_visible(False)
            # self.canvas.draw()
            # self.axes.draw_artist(self.h_line)
        # else:
            # self.h_line.set_visible(True)
            # self.canvas.draw()
            # self.axes.draw_artist(self.h_line)
        if self.pressed == True:
            self.moved = True
        if self.axes != event.inaxes:
            inv = self.axes.transData.inverted()
            new_dataPoint = int(inv.transform(numpy.array((event.x, event.y)).reshape(1, 2)).ravel()[1])
            # print(x,y)
        if event.inaxes == self.axes:
            # print(event.ydata)


            new_dataPoint = int(event.ydata) #(int(numpy.floor(self.combinedTransform.transform(pt_data2)[1])))
            self.canvas.restore_region(self.background)
            self.v_line.set_xdata(event.xdata)
            # self.h_line.set_ydata(self.one_d.data[new_dataPoint,:])
            # self.axes.draw_artist(self.h_line)
            self.axes.draw_artist(self.v_line)
            self.current_h = new_dataPoint
            self.canvas.blit(self.axes.bbox)

    def on_click(self, event):
        self.pressed = True
        self.origin = event.xdata
        print(self.one_d.index)

    def on_release(self, event):
        # print(self.origin, event.x, self.moved)
        if self.moved == False:
            # print('point')
            xs= range(self.one_d.data.shape[0])
            coord = numpy.argmin(numpy.abs(self.one_d.index-event.xdata))
            vline = self.one_d.index[coord]
            # self.scatter_data.append(self.one_d.data[:,coord])
            # self.scatter_data_norm.append(self.one_d.data[:,coord]/numpy.max(self.one_d.data[:,coord]))
            # self.scatters.append(self.axes_proj.scatter(xs, self.one_d.data[:,coord],color='C'+str(self.number_scatters), marker='x'))
            # self.plot_scatters()
            self.verticals.append(self.axes.axvline(vline, color='C'+str(self.number_scatters), linewidth=2, ls='--'))

            # self.axes_proj.autoscale(enable=True, axis='y')
            self.canvas.draw()
            self.number_scatters+=1
        elif self.moved==True and numpy.absolute(self.origin-event.x)>0.2:
            print('box!')
            xs= range(self.one_d.data.shape[0])
            coord = numpy.argmin(numpy.abs(self.one_d.index-event.xdata))
            vline = self.one_d.index[coord]
            coord2 = numpy.argmin(numpy.abs(self.one_d.index-self.origin))
            vline2 = self.one_d.index[coord2]
            left_coord = numpy.min((coord, coord2))
            right_coord = numpy.max((coord, coord2))

            # self.scatter_data.append(numpy.sum(self.one_d.data[:,left_coord:right_coord], axis=1))
            # self.scatter_data_norm.append(numpy.sum(self.one_d.data[:,left_coord:right_coord], axis=1)/numpy.max(numpy.sum(self.one_d.data[:,left_coord:right_coord], axis=1)))
            # self.scatters.append(self.axes_proj.scatter(xs, numpy.sum(self.one_d.data[:,left_coord:right_coord], axis=1),color='C'+str(self.number_scatters), marker='x'))
            # self.plot_scatters()
            # print(vline, numpy.abs(self.one_d.index-event.xdata))
            self.verticals.append(self.axes.axvline(vline, color='C'+str(self.number_scatters), linewidth=0.5, ls='--'))
            self.verticals.append(self.axes.axvline(vline2, color='C'+str(self.number_scatters), linewidth=0.5, ls='--'))
            left = numpy.max((vline, vline2))
            right = numpy.min((vline, vline2))
            # print(right, left)
            color = list(matplotlib.colors.to_rgba('C'+str(self.number_scatters))[:3])
            color.append(0.3)
            print(color)
            print(right, left)
            h0,h1 = self.axes.get_ylim()
            height = numpy.fabs(h1-h0)
            print(height)
            self.rectangles.append(patches.Rectangle((left,h0), right-left, height, linewidth=0, facecolor=color))
            # self.axes.add_patch(self.rectangles[-1])
            print((right+left)/2)
            print(left_coord, right_coord, numpy.sum(self.one_d.data[left_coord:right_coord]))
            self.axes.annotate("%.2f" % numpy.sum(self.one_d.data[left_coord:right_coord]), ((right+left)/2., h1-((self.number_scatters+1)*(height/20.))), xytext=((right+left)/2.0, h1-(self.number_scatters*(height/20.))), color='C'+str(self.number_scatters), ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white', ec='C'+str(self.number_scatters)),
            arrowprops=dict(arrowstyle='-[, widthB='+str(numpy.abs(right-left))+', lengthB=1.5', lw=2.0, color='C'+str(self.number_scatters)))
            # self.axes_proj.autoscale(enable=True, axis='y')
            self.canvas.draw()
            self.number_scatters+=1
        self.pressed=False
        self.moved=False
        # self.h_line.set_visible(False)
        self.v_line.set_visible(False)
        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        # self.h_line.set_visible(True)
        self.v_line.set_visible(True)
        # self.axes.draw_artist(self.h_line)
        self.axes.draw_artist(self.v_line)

    def on_cb_grid(self, event):
        # if(self.tabOne.dim==3 and self.tabOne.DECON==0):
        #     print('No deconvolution data available')
        #     self.cb_calc.SetValue(0)
        self.draw_figure()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.axes = self.fig.add_subplot(111)
        self.dt = FileDrop(self.canvas, self.axes)
        self.canvas.SetDropTarget(self.dt)
        self.selected_plots =None
        self.cursor_shown = False

        self.pressed = False
        self.moved = False

        self.verticals = []
        self.not_yet_drawn = True
        self.number_scatters = 0
        # self.drawbutton = wx.Button(self, -1, "Draw!")
        # self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.toolbar = NavigationToolbar(self.canvas, coordinates=False)
        self.toolbar.Realize()

        # Match the 2D/3D Projection UI: custom controls live directly in
        # Matplotlib's toolbar rather than in separate StaticBox panels.
        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        self.drawing_box()
        self.toolbar.AddSeparator()
        for widget in (self.cb_grid, self.cb_calc, self.drawbutton,
                       self.integrate_button, self.addition_button):
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 4)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)


    def draw_figure(self):
        #try:
        self.draw_figureGO()
        #except:
        #    pass

    def draw_figureGO(self):
        """ Redraws the figure
        """
        self.axes.clear()
        # Open 1D spectrum style: no top/right frame and the x-axis baseline
        # crosses the y-axis at the physical zero-intensity level.
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['bottom'].set_position(('data', 0.0))
        self.axes.xaxis.set_ticks_position('bottom')
        self.axes.yaxis.set_ticks_position('left')
        print('drawing')
        #sele1=self.ComboBox1.GetSelection()
        self.thresh=float(self.one_d.threshold())


        xs=self.one_d.index
        ys=self.one_d.data
        y2s=numpy.zeros_like(ys)
        y2s.fill(self.thresh)

        self.axes.set_xlabel(self.one_d.labels[0],fontsize=8)
        self.data_1d, = self.axes.plot(xs,ys,'r',label='data', linewidth=0.5, picker = 5)
        self.dt.data_1d = self.data_1d
        self.axes.plot(xs,y2s,'g',label='threshold', linewidth=0.5)
        if(self.one_d.deconvolution_enabled):
            print(self.one_d.deconvolution_enabled)
            if(self.cb_calc.GetValue()==1):
                self.axes.plot(xs,self.one_d.datadec,'b',label='deconvolved', linewidth=0.5)
            # Legacy conn_data stick overlays were removed in Stage 148D.
            # The Full Peak List is the authoritative peak collection.

        # Keep y=0 visible on the initial/redraw view so the relocated x-axis
        # is always present.  Normal Matplotlib zooming remains available.
        ymin, ymax = self.axes.get_ylim()
        self.axes.set_ylim(min(ymin, 0.0), max(ymax, 0.0))

        self.xmin,self.xmax=self.axes.get_xlim()
        self.axes.set_xlim(self.xmax, self.xmin)
        self.ymin,self.ymax=self.axes.get_ylim()
        self.offset=-1*self.ymin/2
        self.canvas.mpl_connect('pick_event', self.on_pick)

        self.canvas.draw()



    def on_draw_button(self, event):
        self.draw_figure()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()

    def on_N_button(self, event):
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

        if self.dt.first_drop == False and event.mouseevent.dblclick:
            self.selected_plots=event.artist
            self.selected_plots.set_linewidth(2.)
            self.canvas.draw_idle()

    def on_scroll(self, event):
        # print('scrolling')
        if self.selected_plots == None:
            if event.step != 0.0:
                print(event.step)
                # max_step = max(1.0, event.step)*numpy.sign(event.step)
                max_step = event.step
                self.ymin,self.ymax=self.axes.get_ylim()
                new_ymin = self.ymin+(self.ymin*0.05*max_step)
                new_ymax = self.ymax+(self.ymax*0.05*max_step)
                
                if numpy.sign(new_ymin) == numpy.sign(self.ymin):
                    self.axes.set_ylim(new_ymin, new_ymax)
                else:
                    print(event.step, new_ymin, self.ymin, (self.ymin*0.05*max_step))
        else:
            base_scale=1.02
            # if event.button == 'up':
            # # deal with zoom in
            #     scale_factor = 1./base_scale
            # elif event.button == 'down':
            # # deal with zoom out
            #     scale_factor = base_scale
            scale_factor = numpy.power(base_scale, -event.step) #*event.step
            print(scale_factor)
            self.selected_plots.set_ydata(self.selected_plots.get_ydata()*scale_factor)
            splitted = self.selected_plots.get_label().split('*')
            if len(splitted) >1:
                new = float(splitted[-1])*scale_factor
                self.selected_plots.set_label(splitted[0]+'*%.2f'%new)
            else:
                self.selected_plots.set_label(self.selected_plots.get_label()+' *'+str(scale_factor))
        self.axes.legend()
        self.canvas.draw()

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
