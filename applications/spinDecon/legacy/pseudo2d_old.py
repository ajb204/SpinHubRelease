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
import wx,string,os,numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm
# from matplotlib.widgets import Cursor
import matplotlib.patches as patches
import scipy
# from .vpar_decon import vpar
# from spinDecon.INDIANA.cellDiff import FitDiff
############################################################################
# Frame for 1d slices
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()

class Indiana_dialog(wx.Dialog):
    def __init__(self, parent, little_delta_string, gradient_string, big_delta_string):
        wx.Dialog.__init__(self, parent, wx.ID_ANY, "Name Input", size= (350,220))
        self.p = wx.Panel(self,wx.ID_ANY)
        # self.label_sizer 
        self.lblname = wx.StaticText(self.p, label="Little Delta (s)", pos=(20,15))
        self.little_delta_box = wx.TextCtrl(self.p, value=little_delta_string, pos=(110,15), size=(200,-1))
        self.lblsur = wx.StaticText(self.p, label="Gradient", pos=(20,60))
        self.lblsur = wx.StaticText(self.p, label="Strength (G/cm)", pos=(20,75))
        self.gradient_box = wx.TextCtrl(self.p, value=gradient_string, pos=(110,60), size=(200,-1))
        self.lblnick = wx.StaticText(self.p, label="Big Delta (s)", pos=(20,110))
        self.big_delta_box = wx.TextCtrl(self.p, value=big_delta_string, pos=(110,110), size=(200,-1))
        self.saveButton =wx.Button(self.p, wx.ID_OK, label="OK", pos=(110,160))
        self.closeButton =wx.Button(self.p, label="Cancel", pos=(210,160))
        self.saveButton.Bind(wx.EVT_BUTTON, self.SaveConnString)
        self.closeButton.Bind(wx.EVT_BUTTON, self.OnQuit)
        self.Bind(wx.EVT_CLOSE, self.OnQuit)
        self.Show()

    def OnQuit(self, event):
        self.result_name = None
        self.Destroy()

    def SaveConnString(self, event):
        self.indiana_params={}
        self.indiana_params['delta'] = float(self.little_delta_box.GetValue())
        self.indiana_params['gzlvl1'] = numpy.array(self.gradient_box.GetValue().split(', '), dtype=float)
        self.indiana_params['BigT'] = numpy.array(self.big_delta_box.GetValue().split(', '), dtype=float)
        # print(self.little_delta, self.gradients, self.big_deltas)
        self.EndModal(wx.ID_OK) 
        # self.Destroy()

class Pseudo2D(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.dim=parent.tabOne.dim
        self.tabOne=tabOne
        self.sum=(0.,2.)
        self.peak=tabOne.peak
        self.thresh=tabOne.dmax
        self.offset=0
        self.indiana=False

        dmax=tabOne.uc0max
        dmin=tabOne.uc0min

        self.create_main_panel()
        self.draw_figure()
        #self.Show()
        self.Fit()

    def onFocus(self, event):
        print("Pseudo has focus!")

    def drawing_box(self):
        self.vbox2Lbl = wx.StaticBox(self,-1,'Drawing:')
        self.vbox2=wx.StaticBoxSizer(self.vbox2Lbl,wx.HORIZONTAL)

        self.drawbutton = wx.Button(self, -1, "Draw!", size=(-1,22))
        self.cb_grid = wx.CheckBox(self, -1,"Peaks",style=wx.ALIGN_RIGHT)
        self.cb_calc = wx.CheckBox(self, -1,"ShowCalc",style=wx.ALIGN_RIGHT)

        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_calc)

        self.vbox2.Add(self.cb_grid, border=10, flag=self.flags)
        self.vbox2.Add(self.cb_calc, border=10, flag=self.flags)
        self.vbox2.Add(self.drawbutton, border=10, flag=self.flags)
        self.vbox2.AddSpacer(10)

    def fitting_box(self):
        vbox2Lbl = wx.StaticBox(self,-1,'Fitting:')
        vbox2=wx.StaticBoxSizer(vbox2Lbl,wx.HORIZONTAL)

        self.norm_button = wx.ToggleButton(self, -1, "Normalize", size=(-1,22))
        self.ST_button = wx.ToggleButton(self, -1, "Stejskal Tanner", size=(-1,22))
        self.T1_button = wx.ToggleButton(self, -1, "T1", size=(-1,22))
        self.indiana_button = wx.ToggleButton(self, -1, "INDIANA", size=(-1,22))
        cb_grid = wx.CheckBox(self, -1,"Peaks",style=wx.ALIGN_RIGHT)
        cb_calc = wx.CheckBox(self, -1,"ShowCalc",style=wx.ALIGN_RIGHT)

        self.Bind(wx.EVT_TOGGLEBUTTON, self.on_norm_button, self.norm_button)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.fit_ST_equation, self.ST_button)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.fit_T1, self.T1_button)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.fit_indiana_button, self.indiana_button)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, cb_grid)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, cb_calc)

        vbox2.Add(cb_grid, border=10, flag=self.flags)
        vbox2.Add(cb_calc, border=10, flag=self.flags)
        vbox2.Add(self.norm_button, border=10, flag=self.flags)
        vbox2.Add(self.ST_button, border=10, flag=self.flags)
        vbox2.Add(self.T1_button, border=10, flag=self.flags)
        vbox2.Add(self.indiana_button, border=10, flag=self.flags)
        vbox2.AddSpacer(10)
        return vbox2

    def fit_indiana_button(self, event):
        if self.indiana_button.GetValue() == True:

            self.indiana_input_file = 'input.txt'
            big_delta_string = ''
            gradient_string = ''
            little_delta_string = ''
            # dac = ''
            for line in open(self.indiana_input_file):
                if '#' not in line:
                    fields = line.split()
                    try:
                        if fields[0] == 'dac':
                            dac = float(fields[1])
                    except:
                        pass
            for line in open(self.indiana_input_file):
                if '#' not in line:
                    fields = line.split()
                    
                    try:
                        numbers = numpy.array(fields[1:], dtype=float)
                        name = fields[0]
                        if name == 'BigT':
                            for x in numbers:
                                big_delta_string = big_delta_string+"{:.3f}".format(x)+', '
                            big_delta_string = big_delta_string[:-2]
                        if name == 'delta':
                            for x in numbers:
                                little_delta_string = little_delta_string+"{:.6f}".format(x)+', '
                            little_delta_string = little_delta_string[:-2]
                        if name == 'grad':
                            numbers =numbers*dac
                            for x in numbers:
                                gradient_string = gradient_string+"{:.3f}".format(x)+', '
                            gradient_string = gradient_string[:-2]        

                    except:
                        print()
                
                    # print(fields)
            
            
            dlg = Indiana_dialog(self, little_delta_string, gradient_string, big_delta_string)
            if dlg.ShowModal() == wx.ID_OK:
                self.indiana = True
                self.indiana_bvalues = []
                self.indiana_params = dlg.indiana_params
                grad = self.indiana_params['gzlvl1']
                bigT = self.indiana_params['BigT']
                # for i in range(len(bigT)):
                for j in range(len(grad)):
                    for i in range(len(bigT)):

                # print(BigT_index[i] + gzlvl_index[j]*len(bigT))
                        self.indiana_bvalues.append((bigT[i], grad[j], bigT[i]*(self.indiana_params['delta']*grad[j]*scipy.constants.value('proton gyromag. ratio')*1e-4)**2))
                print(self.indiana_bvalues)
                # sortedLine[i,j] = line[BigT_index[i] + gzlvl_index[j]*len(bigT)]
                # self.plot_scatters()

        else:
            self.indiana = False
            self.plot_scatters()
        
                

    def on_norm_button(self, event):
        self.plot_scatters()

    def plot_scatters(self):
        axis_font = {'fontname':'Arial', 'size':'14'}
        for x in self.scatters:
            x.remove()
        self.axes_proj.cla()
        self.scatters=[]
        number_scatters = 0
        if self.norm_button.GetValue()==True:
            data = self.scatter_data_norm
        else:
            data = self.scatter_data
        if len(self.gzlvl1)==self.tabOne.data.shape[0]:
            xs=numpy.array(self.gzlvl1)
        elif len(self.T1s)==self.tabOne.data.shape[0]:
            xs=numpy.array(self.T1s)
        else:
            xs= range(self.tabOne.data.shape[0])
            
        # self.indiana_bvalues = []

        self.info_text = self.axes.text(0.1,0.9, '', transform=self.axes.transAxes)
        
        for line in data:
            print('line: ', line)
            if self.indiana:
                print('Indiana not implemented yet')
                # grad = self.indiana_params['gzlvl1']
                # bigT = self.indiana_params['BigT']
                # sortedLine = numpy.zeros([len(bigT), len(grad)])
                # BigT_index = numpy.argsort(bigT)
                # gzlvl_index = numpy.argsort(grad)
                # for i in range(len(bigT)):
                #     for j in range(len(grad)):
                #         # print(BigT_index[i] + gzlvl_index[j]*len(bigT))
                #         # self.indiana_bvalues.append((bigT[i], grad[j], bigT[i]*(self.indiana_params['gzlvl1']*grad[j]*scipy.constants.value('proton gyromag. ratio')*1e-4))**2)
                #         sortedLine[i,j] = line[BigT_index[i] + gzlvl_index[j]*len(bigT)]

                
                # self.fitting = FitDiff(sortedLine, self.indiana_params, meth='bi')
                # self.fitting.mode='mono'

                # try:
                #     self.fitting.fitFunc_rest(boot_num = 0, para_flg = 'n', ncpus=0)
                
                #     print(self.fitting.result)
                #     colors = ['red', 'orange', 'yellow', 'green', 'lightblue', 'blue']
                #     if len(data)>1:
                #         self.fitting.plot_results(self.axes_proj, self.fitting.result, self.fitting.BigT, self.fitting.dat, color='C'+str(number_scatters))
                #     else:
                #         self.fitting.plot_results(self.axes_proj, self.fitting.result, self.fitting.BigT, self.fitting.dat)
                # except:
                #     pass
            elif len(self.T1s)!=self.tabOne.data.shape[0]:
                self.scatters.append(self.axes_proj.scatter(xs, line,color='C'+str(number_scatters), marker='x'))
            
            elif len(self.gzlvl1)==self.tabOne.data.shape[0]:
                line = numpy.log(line)
                self.scatters.append(self.axes_proj.scatter(xs, line,color='C'+str(number_scatters), marker='x'))
                m,b = numpy.polyfit(xs, line, 1)
                self.axes_proj.plot(xs, m*xs+b, ls='--', color='C'+str(number_scatters))
                self.axes_proj.text(0.1,0.9-float(number_scatters)*0.05, '$Deff = $%.2e $cm^2 s^-1$' % m, transform=self.axes_proj.transAxes, color='C'+str(number_scatters),  **axis_font)
            else:
                line = numpy.log((-line+1.)/2.)
                self.scatters.append(self.axes_proj.scatter(xs, line,color='C'+str(number_scatters), marker='x'))
                m,b = numpy.polyfit(xs, line, 1)
                self.axes_proj.plot(xs, m*xs+b, ls='--', color='C'+str(number_scatters))
                self.axes_proj.text(0.1,0.9-float(number_scatters)*0.05, 'R1 = $%.2e s^-1$' % m, transform=self.axes_proj.transAxes, color='C'+str(number_scatters),  **axis_font)

            number_scatters +=1
        self.canvas.draw()


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
        self.axes = self.fig.add_subplot(121)
        self.axes_h = self.axes.twinx()
        self.axes_proj = self.fig.add_subplot(122)
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
        # self.drawbutton = wx.Button(self, -1, "Draw!")
        # self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        # self.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.toolbar = NavigationToolbar(self.canvas)

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        self.drawing_box()
        self.fit_box = self.fitting_box()

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vbox2)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.fit_box)

        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.AddSpacer(10)
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
        self.axes_h.clear()
        print('drawing')
        #sele1=self.ComboBox1.GetSelection()
        self.thresh=float(self.tabOne.dmax*float(self.tabOne.threshBox.GetValue()))


        xs=self.tabOne.index1
        ys=numpy.arange(self.tabOne.data.shape[0])
        zs=self.tabOne.data
        y2s=numpy.zeros_like(ys)
        y2s.fill(self.thresh)
        levels = [self.thresh]
        for x in range(12):
            levels.append(levels[-1]*1.4)

        cmap = plt.get_cmap('Oranges')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        self.axes.set_xlabel(self.parent.tabOne.labb[1],fontsize=8)
        self.axes.pcolormesh(xs,ys,zs,label='data', norm=norm, cmap = cmap)
        self.canvas.draw()
        self.axes.set_ylim(ys[0], ys[-1])
        self.axes.set_xlim(xs[0], xs[-1])

        for x in range(self.tabOne.data.shape[0]):
            self.axes_proj.plot(xs, self.tabOne.data[x,:], lw=0.5, label=str(x))
        # self.axes_proj.legend()
        self.axes_proj.set_xlim(xs[0], xs[-1])

        # self.combinedTransform = self.axes_h.transData + self.axes.transData.inverted()
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
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

            self.rectangles=[]
            self.number_scatters = 0
            self.axes_proj.cla()
            self.canvas.draw()

        # if event.key=='v':
            # if not self.cursor_shown:
            #     print(self.cursor_shown)
                # self.cursor = Cursor(self.axes, horizOn = False, vertOn=True, color='r', linewidth=2, useblit=True)
            # else:
            #     self.cursor.clear()
            # self.canvas.draw()

    def fit_T1(self, event):

        gzread = (os.popen('cd raw && vpar d3 && cd ..').read())
        fields = gzread.split('\n')
        T1s = []
        if len(fields) > 2:
            rs = fields[2].split(' ')
            for r in rs:
                try:
                    # if float(r)>2:
                        ncyc = float(r)

                        # print(delta)
                        T1s.append(r)
                except:
                    pass
        # self.axes.set_xticks(T1s)
        # self.axes.set_yticklabels(T1s)
        self.T1s = T1s

        self.plot_scatters()

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
        self.axes.set_yticks(numpy.arange(len(gzlvl1)))
        self.axes.set_yticklabels(gzlvl1)
        self.gzlvl1 = gzlvl1

        self.plot_scatters()

    def on_scroll(self, event):
        # print('scrolling')
        self.ymin,self.ymax=self.axes_h.get_ylim()
        print(self.ymin, self.ymax)
        self.axes_h.set_ylim(self.ymin+(self.ymin*0.05*event.step), self.ymax+(self.ymax*0.05*event.step))
        self.axes_h.draw_artist(self.h_line)
        # self.axes_h.draw()

    def on_mouse_move(self,event):
        if event.inaxes==self.axes_h:
            if self.not_yet_drawn == True:

                self.background = self.canvas.copy_from_bbox(self.axes.bbox)

                self.current_h = 0
                self.v_line = self.axes_h.axvline(self.tabOne.index1[0], color = 'r', linewidth=2)
                # print(self.tabOne.index1.shape, self.tabOne.data[0,:].shape)
                # self.h_line, = self.axes_h.plot(self.tabOne.index1, self.tabOne.data[0,:], color='k', linewidth = 0.5, zorder=1000)
                self.h_line, = self.axes_h.plot(self.tabOne.index1, numpy.zeros_like(self.tabOne.index1), color='k', linewidth = 0.5)
                # self.h_line.set_visible(False)
                self.axes_h.set_ylim(numpy.min(self.tabOne.data), numpy.max(self.tabOne.data))


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
                # print(x,y)
            # if event.inaxes == self.axes:
            #     print(event.ydata)


                # new_dataPoint = int(event.ydata) #(int(numpy.floor(self.combinedTransform.transform(pt_data2)[1])))
                self.canvas.restore_region(self.background)
                self.v_line.set_xdata(event.xdata)
                self.h_line.set_ydata(self.tabOne.data[new_dataPoint,:])
                if self.indiana:
                    # self.info_text.set_label(str(self.indiana_bvalues[new_dataPoint][0])+' '+str(self.indiana_bvalues[new_dataPoint][1]))
                    print(str(self.indiana_bvalues[new_dataPoint][0])+' '+str(self.indiana_bvalues[new_dataPoint][1]))
                self.axes_h.draw_artist(self.h_line)
                self.axes_h.draw_artist(self.v_line)
                self.current_h = new_dataPoint
                self.canvas.blit(self.axes.bbox)

                # if new_dataPoint != self.current_h:


    def on_draw_button(self, event):
        self.draw_figure()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()

    def on_N_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.draw_figure()


    def on_pick(self, event):
        if event.inaxes==self.axes_h:
            self.pressed = True
            self.origin = event.xdata
            print(self.tabOne.index0)


    def on_release(self, event):
        # print(self.origin, event.x, self.moved)
        if event.inaxes==self.axes_h:
            if self.moved == False:
                # print('point')
                xs= range(self.tabOne.data.shape[0])
                coord = numpy.argmin(numpy.abs(self.tabOne.index1-event.xdata))
                vline = self.tabOne.index1[coord]
                self.scatter_data.append(self.tabOne.data[:,coord])
                self.scatter_data_norm.append(self.tabOne.data[:,coord]/numpy.max(self.tabOne.data[:,coord]))
                # self.scatters.append(self.axes_proj.scatter(xs, self.tabOne.data[:,coord],color='C'+str(self.number_scatters), marker='x'))
                self.plot_scatters()
                self.verticals.append(self.axes_h.axvline(vline, color='C'+str(self.number_scatters), linewidth=2, ls='--'))

                # self.axes_proj.autoscale(enable=True, axis='y')
                self.canvas.draw()
                self.number_scatters+=1
            elif self.moved==True and numpy.absolute(self.origin-event.x)>2.:
                # print('integrate')
                xs= range(self.tabOne.data.shape[0])
                coord = numpy.argmin(numpy.abs(self.tabOne.index1-event.xdata))
                vline = self.tabOne.index1[coord]
                coord2 = numpy.argmin(numpy.abs(self.tabOne.index1-self.origin))
                vline2 = self.tabOne.index1[coord2]
                left_coord = numpy.min((coord, coord2))
                right_coord = numpy.max((coord, coord2))
                print(self.tabOne.data[:,left_coord:right_coord], numpy.sum(self.tabOne.data[:,left_coord:right_coord], axis=1))
                self.scatter_data.append(numpy.sum(self.tabOne.data[:,left_coord:right_coord], axis=1))
                self.scatter_data_norm.append(numpy.sum(self.tabOne.data[:,left_coord:right_coord], axis=1)/numpy.max(numpy.sum(self.tabOne.data[:,left_coord:right_coord], axis=1)))
                # self.scatters.append(self.axes_proj.scatter(xs, numpy.sum(self.tabOne.data[:,left_coord:right_coord], axis=1),color='C'+str(self.number_scatters), marker='x'))
                self.plot_scatters()
                # print(vline, numpy.abs(self.tabOne.index1-event.xdata))
                self.verticals.append(self.axes_h.axvline(vline, color='C'+str(self.number_scatters), linewidth=2, ls='--'))
                self.verticals.append(self.axes_h.axvline(vline2, color='C'+str(self.number_scatters), linewidth=2, ls='--'))
                left = numpy.max((vline, vline2))
                right = numpy.min((vline, vline2))
                # print(right, left)
                color = list(matplotlib.colors.to_rgba('C'+str(self.number_scatters))[:3])
                color.append(0.3)
                # print(color)
                self.rectangles.append(patches.Rectangle((left,-0.5), right-left, 4, linewidth=0, facecolor=color))
                self.axes.add_patch(self.rectangles[-1])
                # self.axes_proj.autoscale(enable=True, axis='y')
                self.canvas.draw()
                self.number_scatters+=1
            self.pressed=False
            self.moved=False
            self.h_line.set_visible(False)
            self.v_line.set_visible(False)
            self.canvas.draw()
            self.background = self.canvas.copy_from_bbox(self.axes.bbox)
            self.h_line.set_visible(True)
            self.v_line.set_visible(True)
            self.axes_h.draw_artist(self.h_line)
            self.axes_h.draw_artist(self.v_line)


    # def on_scroll(self, event):
    #     print('scrolling')
    #     self.ymin,self.ymax=self.axes.get_ylim()
    #     self.axes.set_ylim(self.ymin+(self.ymin*0.05*event.step), self.ymax+(self.ymax*0.05*event.step))
    #     self.canvas.draw()

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
