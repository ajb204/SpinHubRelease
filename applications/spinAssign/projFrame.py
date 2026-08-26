#!/usr/bin/python
import wx,string,copy,math,numpy,os,platform,sys
import matplotlib            #import matplotlib
#matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure

from wx.lib.mixins.listctrl import ColumnSorterMixin
import re
###############################################3
#Show 2D projections of the 3Ds with peaklists

def ParseFlt(infile,param):
    #print param
    par=0
    if(os.path.exists(infile)==0):
        #print 'Cannot find input file'
        pass
    else:
        inny=open(infile,'r')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>2):
                if(test[0]==param and test[1]=='='):
                    par=float(test[2])
                    break
    #deal with defaults
    # if(param=='sig1' and par==0):
    #     par= 0.2
    # if(param=='sig2'  and par==0):
    #     par= 0.2
    # if(param=='sig3'  and par==0):
    #     par= 0.2
    # if(param=='sig4'  and par==0):
    #     par= 0.2
    # if(param=='thresh' and par==0):
    #     par= 0.01
    # if(param=='fac'  and par==0):
    #     par= 1.6
    if(param=='squash'  and par==0):
        par= 0.65
    # if(param=='maxiter'  and par==0):
    #     par= 100

    #no default
    return par

def Parse(infile,param):
    if(os.path.exists(infile)==0):
        print('Cannot find input file')
        return 0
    inny=open(infile,'r')
    for line in inny.readlines():
        test=line.split()
        if(len(test)>2):
            if(test[0]==param and test[1]=='='):
                return (' '.join(test[2:]))
    return 0


class Click():
    def __init__(self, canvas, func, button=1):
        self.canvas=canvas
        #print('paddyxy')
        self.func=func
        self.button=button
        self.press=False
        self.move = False
        self.c1=self.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.canvas.mpl_connect('motion_notify_event', self.onmove)
    def onclick(self,event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event, self.ax)
    def onpress(self,event):
        print('Pressed')
        self.press=True
    def onmove(self,event):
        if self.press:
            self.move=True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False


class projFrame(wx.Panel):

    def __init__(self,parent):

        #self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Panel.__init__(self, parent=parent)

        self.parent=parent
        self.create_main_panel()

        #self.makeref()
        self.setup_figure_data()
        
        self.draw_figure()

        
        #self.Show(True)
        self.Fit()



    def runDecon(self):

        self.build=platform.uname()[0]
        self.deconBin='decon_'+self.build

        specs=list(self.parent.tabOne.molecule.spec.keys())


        for spec in specs:
            specDeconParFile=spec+'/deconParFile'
            indir=Parse(specDeconParFile, 'indir')
            thresh=ParseFlt(specDeconParFile,'thresh')
            ncpus=int(Parse(specDeconParFile,'ncpus'))
            fac=ParseFlt(specDeconParFile,'fac')
            squash=ParseFlt(specDeconParFile,'squash')
            maxiter=int(ParseFlt(specDeconParFile,'maxiter'))
            infile=Parse(specDeconParFile,'infile')
            peakfile=Parse(specDeconParFile,'peakfile')
            voigt1=ParseFlt(specDeconParFile,'voigt1')
            voigt2=ParseFlt(specDeconParFile,'voigt2')
            voigt3=ParseFlt(specDeconParFile,'voigt3')
            sig1=ParseFlt(specDeconParFile,'sig1')
            sig2=ParseFlt(specDeconParFile,'sig2')
            sig3=ParseFlt(specDeconParFile,'sig3')
            lor1=ParseFlt(specDeconParFile,'lor1')
            lor2=ParseFlt(specDeconParFile,'lor2')
            lor3=ParseFlt(specDeconParFile,'lor3')
            symmode=Parse(specDeconParFile,'symmode')
            dmax=numpy.max(self.parent.tabOne.molecule.spec[spec].data)
            print(indir)
            if str(indir) != '0':
                os.chdir(str(indir))
            else:
                return
            specstr=self.deconBin
            specstr+=' '+str(ncpus)
            specstr+=' '+str(peakfile)
            specstr+=' '+str(3)
            specstr+=' '+str(infile)
            specstr+=' '+str(dmax*thresh)
            specstr+=' '+str(sig1)
            specstr+=' '+str(sig2)
            specstr+=' '+str(sig3)
            specstr+=' '+str(fac)
            specstr+=' '+str(squash)
            specstr+=' '+str(0)
            specstr+=' '+str(1.0)
            specstr+=' '+str(7)
            specstr+=' '+str(voigt1)
            specstr+=' '+str(voigt2)
            specstr+=' '+str(voigt3)
            specstr+=' '+str(lor1)
            specstr+=' '+str(lor2)
            specstr+=' '+str(lor3)
            # print(specstr)
            os.system(specstr)


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.

        self.first_draw = 1
        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        # click = Click(self.canvas, self.on_pick)
        self.canvas.mpl_connect('button_release_event', self.on_right_click)
        self.canvas.mpl_connect('key_press_event', self.keyboard_press)
        self.canvas.mpl_connect('motion_notify_event', self.draw_bores)
        self.cmaps =  'seismic','bwr','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'

        self.toolbar = NavigationToolbar(self.canvas)

        self.newPeakEntries = []


        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(5)


        self.axes = self.fig.add_subplot(111)

        self.fig2 = Figure()
        self.canvas2 = FigCanvas(self, -1, self.fig2)

        self.canvas2.SetMinSize(wx.Size(100,100))

        specs=self.parent.tabOne.molecule.spec.keys()


        self.axes_bores = self.fig2.subplots(1,5)
        self.axes2 = {}
        self.axes2['hnco'] = self.axes_bores[0]
        self.axes2['hnca'] = self.axes_bores[1]
        
        

        self.axes2['hncaco'] = self.axes_bores[0].twinx()
        self.axes2['hncoca'] = self.axes_bores[1].twinx()
        i = 2
        if 'hncacb' in specs:
            self.axes2['hncacb'] = self.axes_bores[i]

            if('cbcaconh' in specs):
                self.axes2['cbcaconh'] = self.axes_bores[i]
            elif('hncocacb in specs'):
                self.axes2['hncocacb'] = self.axes_bores[i].twinx()
            i+=1

        if 'hncanh' in specs:
            self.axes2['hncanh'] = self.axes_bores[i]
            self.axes2['hncocanh'] = self.axes_bores[i].twinx()
            i+=1

        if 'ctocsy' in specs:
            self.axes2['ctocsy'] = self.axes_bores[i]
            i+=1
        if 'hcconh' in specs:
            self.axes2['hcconh'] = self.axes_bores[i]
            i+=1

        #if 'cbcaconh' in specs:
            
        #    i+=1



        for ax_name in specs:

            # self.axes2[ax_name].set_ylim((0,0))
            self.axes2[ax_name].spines['top'].set_visible(False)
            self.axes2[ax_name].spines['right'].set_visible(False)
            self.axes2[ax_name].spines['left'].set_visible(False)
            self.axes2[ax_name].set_yticks([])
        for ax in self.axes_bores.flatten():
          ax.get_yaxis().set_visible(False)
        self.bores = {}
        self.axes2ymin = {}
        self.axes2ymax = {}

        # self.vbox.Add(self.canvas2, 1, wx.LEFT | wx.TOP | wx.GROW)

        self.selected=[]
        self.proj={}
        self.text1={}
        self.text2={}
        self.text3={}
        self.textbox1={}
        self.textbox2={}
        self.textbox3={}
        self.cblab={}
        self.cbpk={}
        self.cblist={}
        self.cntrLbl={}
        self.cntrSizer={}
        self.comboClbox={}
        self.launch={}
        self.refbox={}
        self.add = 0
        self.select_key = 0
        self.contour_data = {}
        self.peak_data = {}
        self.text_objects = {}

        self.maxy = -1
        flags = wx.ALIGN_LEFT | wx.ALL #| wx.ALIGN_CENTER_VERTICAL
        self.thresh=1E5
        self.vboxCntr=wx.BoxSizer(wx.VERTICAL)


        spec_logic = {
            'hnco':1,
            'hncoca':1,
            'hncocacb':1,
            'hncocanh':1,
            'hnca':2,
            'hncaco':2,
            'hncacb':2,
            'hncanh':2,
            'ctocsy':1,
            'hcconh':1,
            'cbcaconh':1,
        }

        for i,spec in enumerate(specs):
            self.text1[spec]=wx.StaticText(self, -1, 'Min:')
            self.text2[spec]=wx.StaticText(self, -1, 'Fac:')
            self.text3[spec]=wx.StaticText(self, -1, 'Num:')
            self.textbox1[spec]=wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
            self.textbox2[spec]=wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)
            self.textbox3[spec]=wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)
            self.cblist[spec]=wx.CheckBox(self, -1,)
            self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cblist[spec])
            self.cblist[spec].SetValue(True)


            self.comboClbox[spec]=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.cmaps, style=wx.CB_READONLY)


            self.comboClbox[spec].SetSelection(i)


            self.cblab[spec]=wx.CheckBox(self, -1,)
            self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cblab[spec])
            self.cblab[spec].SetValue(True)

            self.cbpk[spec]=wx.CheckBox(self, -1,)
            self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cbpk[spec])
            #####
            if(spec==self.parent.tabOne.molecule.refSpec):
                self.cbpk[spec].SetValue(True)
            else:
                self.cbpk[spec].SetValue(False)
            #self.launch[spec]=wx.Button(self, -1, "go",size=(30,-1))
            #self.Bind(wx.EVT_BUTTON, self.on_launch, self.launch[spec])
            self.textbox1[spec].SetValue(str(self.parent.tabOne.molecule.spec[spec].noise))
            self.textbox2[spec].SetValue(str(1.2))
            self.textbox3[spec].SetValue(str(15))
            self.cntrSizer[spec] = wx.BoxSizer(wx.HORIZONTAL)
            self.cntrSizer[spec].Add(wx.StaticText(self,-1,spec,size=(65,-1)))
            self.cntrSizer[spec].Add(self.cblist[spec])
            self.cntrSizer[spec].Add(self.cblab[spec])
            self.cntrSizer[spec].Add(self.cbpk[spec])
            self.cntrSizer[spec].Add(self.text1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.comboClbox[spec], 0, border=3, flag=flags)
            #self.cntrSizer[spec].Add(self.launch[spec], 0, border=3, flag=flags)

            self.vboxCntr.Add(self.cntrSizer[spec], 0)
            if spec_logic[spec] == 1:
                self.bores[spec], = self.axes2[spec].plot(self.parent.tabOne.molecule.spec[spec].index0-self.parent.tabOne.molecule.spec[spec].ref, self.parent.tabOne.molecule.spec[spec].data[:,10,10], lw=0.7, color='r', label='i-1')
            else:    
                self.bores[spec], = self.axes2[spec].plot(self.parent.tabOne.molecule.spec[spec].index0-self.parent.tabOne.molecule.spec[spec].ref, self.parent.tabOne.molecule.spec[spec].data[:,10,10], lw=0.7, color='darkblue', label='i and i-1')
            if spec != 'ctocsy':
                self.axes2[spec].set_title(spec[-2:].upper())
            else:
                self.axes2[spec].set_title('TOCSY')
            self.axes2ymin[spec] = 0
            self.axes2ymax[spec] = 0

        if('hnco' in specs and 'hncaco' in specs):
            self.axes2[self.parent.tabOne.molecule.refSpec].legend([self.bores[self.parent.tabOne.molecule.refSpec], self.bores['hncaco']], ['i-1', 'i and i-1'], frameon=False, framealpha=0.0)

        if 'hncanh' in specs:
            self.axes2['hncanh'].legend([self.bores['hncanh'], self.bores['hncocanh']], ['i-1, i, i+1', 'i-1'], frameon=False, framealpha=0.0)
        # self.axes2['hncaco'].legend(frameon=False, framealpha=0.0)

        self.fig.tight_layout()
        self.fig2.tight_layout()


        self.hbox=wx.BoxSizer(wx.HORIZONTAL)

        self.hbox_buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.referenced_text = wx.StaticText(self, -1, 'Referenced: ')
        self.referenced = wx.CheckBox(self, -1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.referenced)
        self.hbox_buttons.Add(self.referenced_text, 0, wx.LEFT | wx.TOP)
        self.hbox_buttons.Add(self.referenced, 0, wx.TOP | wx.LEFT)
        self.folded_text = wx.StaticText(self, -1, 'Unfold: ')
        self.folded = wx.CheckBox(self, -1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.folded)

        self.hbox_buttons.Add(self.folded_text, 0, wx.LEFT | wx.TOP)

        self.hbox_buttons.Add(self.folded, 0, wx.TOP | wx.LEFT)
        self.peakListButton = wx.Button(self, -1, "Peak List")
        self.Bind(wx.EVT_BUTTON, self.onPeakList, self.peakListButton)
        self.drawButton = wx.Button(self, -1, "Draw")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawButton)
        self.addButton = wx.ToggleButton(self, -1, "Add")
        self.Bind(wx.EVT_TOGGLEBUTTON, self.on_add_button, self.addButton)
        self.saveButton = wx.Button(self, -1, "SaveNDecon")
        self.Bind(wx.EVT_BUTTON, self.on_save_button, self.saveButton)


        self.deleteButton = wx.Button(self, -1, "Delete")
        self.Bind(wx.EVT_BUTTON, self.on_delete_button, self.deleteButton)
        self.selectButton = wx.ToggleButton(self, -1, "Select")
        self.Bind(wx.EVT_TOGGLEBUTTON, self.on_select_button, self.selectButton)

        self.hbox_buttons.Add(self.peakListButton)
        self.hbox_buttons.Add(self.drawButton)
        self.hbox_buttons.Add(self.addButton)
        self.hbox_buttons.Add(self.saveButton)
        self.vboxCntr.Add(self.hbox_buttons, 0, wx.TOP | wx.LEFT, 3)

        self.hbox_buttons2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox_buttons2.Add(self.deleteButton)
        self.hbox_buttons2.Add(self.selectButton)
        self.vboxCntr.Add(self.hbox_buttons2, 0, wx.TOP | wx.LEFT, 3)





        self.peakfileLab=wx.StaticText(self,-1,"Save to:")
        self.peakfileBox = wx.TextCtrl(self,size=(150,22),style=wx.TE_PROCESS_ENTER)
        self.peakfileBox.SetValue('/raw/test.ft3.list')

        self.hbox_buttons3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox_buttons3.Add(self.peakfileLab)
        self.hbox_buttons3.Add(self.peakfileBox)
        self.vboxCntr.Add(self.hbox_buttons3, 0, wx.TOP | wx.LEFT, 6)

        self.vboxCntr.Add(self.canvas2, 10, wx.LEFT | wx.TOP | wx.EXPAND)




        self.hbox.Add(self.vboxCntr, 1, wx.TOP | wx.GROW)
        self.hbox.Add(self.vbox, 1, wx.TOP | wx.GROW)
        self.SetSizerAndFit(self.hbox)


        self.ccyc=[u'b', u'g', u'r', u'c', u'm', u'y', u'k', u'darkorange', u'pink']

        

    def on_cb_grid(self,event):
        self.draw_figure()

    def on_draw_button(self, event):
        self.first_draw = 1
        self.draw_figure()

    def on_save_button(self, event):

        outfile=self.peakfileBox.GetValue()

        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
            print('Saving to '+spec+outfile)
            os.system('cp '+spec+outfile+' '+spec+outfile+'.backup')
            outy = open(spec+outfile, 'w')
            namesForOut = self.peak_data[spec][4][:, None]
            outWrite = numpy.append(namesForOut, self.peak_data[spec][3], axis = 1)
            for row in outWrite:
                #print(row)
                outy.write(("%s\t%s\t%s\n" % (row[0], row[2], row[1])))




    def on_add_button(self, event):
        print(self.add)
        if self.add == 0:
            self.add = 1
            self.control_text.set_text('Click to add peak')
            self.addButton.SetValue(True)
            print(self.add)
        else:
            self.add = 0
            self.control_text.set_text('')
            self.addButton.SetValue(False)
        self.draw_figure()
        self.canvas.draw()


    def on_select_button(self, event):
        if self.select_key == 0:
            self.select_key = 1
            self.control_text.set_text('Click to select')
            self.selectButton.SetValue(True)
        else:
            self.select = 0
            self.control_text.set_text('')
            self.selectButton.SetValue(False)
        self.draw_figure()
        self.canvas.draw()


    def on_delete_button(self, event):
        if len(self.selected)>0:
            for peak in self.selected:
                peak = str(peak)
                index = self.peak_data['hnco'][4].index(peak)
                for i, spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
                    self.peak_data[spec][0] = numpy.delete(self.peak_data[spec][0], index, 0)
                    self.peak_data[spec][1] = numpy.delete(self.peak_data[spec][1], index, 0)
                    self.peak_data[spec][2] = numpy.delete(self.peak_data[spec][2], index, 0)
                    self.peak_data[spec][3] = numpy.delete(self.peak_data[spec][3], index, 0)
                    self.peak_data[spec][4] = numpy.delete(self.peak_data[spec][4], index, 0)
                    self.peak_data[spec][4] = self.peak_data[spec][4].tolist()
        self.draw_figure()
        self.canvas.draw()





    def keyboard_press(self, event):
        if event.key=='c':
            for i, spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
                self.axes2[spec].clear()
                self.bores[spec], = self.axes2[spec].plot(self.parent.tabOne.molecule.spec[spec].index0-self.parent.tabOne.molecule.spec[spec].ref, self.parent.tabOne.molecule.spec[spec].data[:,10,10])
                self.axes2ymin[spec] = 0
                self.axes2ymax[spec] = 0
        if event.key=='a':
            if self.add == 0:
                self.add = 1
                self.control_text.set_text('Click to add peak')
                self.addButton.SetValue(True)
                print(self.add)
            else:
                self.add = 0
                self.control_text.set_text('')
                self.addButton.SetValue(False)
        if event.key=='s':
            if self.select_key == 0:
                self.select_key = 1
                self.control_text.set_text('Selecting')
            else:
                self.select_key = 0
                self.control_text.set_text('')
                self.selected = []
        self.canvas.draw()

    def on_right_click(self, event):
        if event.button == 3:
            x = event.xdata
            y = event.ydata

            #print(self.parent.tabOne.molecule.spec['hnco'].index2)
            if x and y:
                for i, spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
                  #print(spec)
                  #if spec != 'hnco':
                    x1, y1 = numpy.abs(self.parent.tabOne.molecule.spec[spec].index2-x).argmin(), numpy.abs(self.parent.tabOne.molecule.spec[spec].index1-y).argmin()
                    self.axes2[spec].plot(self.parent.tabOne.molecule.spec[spec].index0-self.parent.tabOne.molecule.spec[spec].ref, self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])
                    max, min = numpy.max(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1]), numpy.min(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])
                    if max > self.axes2ymax[spec]:
                        self.axes2ymax[spec] = max
                    if min < self.axes2ymin[spec]:
                        self.axes2ymin[spec] = min





    def on_pick(self, event):

        print('Click at:', event.xdata, event.ydata)
        if self.add == 1:
            self.add_peak(event.xdata, event.ydata)
            # self.runDecon()
            # exit()
            self.draw_figure()
            return
        if self.select_key == 1:
            self.selected = []
            print('woooo')
            print(len(self.selected))
            # self.select.append((event.xdata,event.ydata))

            x_min, x_max = self.axes.get_xlim()
            y_min, y_max = self.axes.get_ylim()
            xdist = x_max - x_min
            ydist = y_max - y_min
            maxy_spec={}

            for spec in self.text_objects.keys():

                if self.maxy > -1:
                    self.text_objects[spec][self.maxy].set_color('k')
                raddy = []
                for p in self.text_objects[spec]:
                    # print(self.peak[p].name)
                    xval,yval = p.get_position()  # proton
                    # yval = p.y  # nitrogen
                    rad2 = ((xval - event.xdata) / xdist) ** 2. + ((yval - event.ydata) / ydist) ** 2.
                    raddy.append(rad2)


                raddy = numpy.array(raddy)

                self.maxy = numpy.argmin(raddy)

                for peak in self.text_objects[spec]:
                    peak.set_alpha(0.55)

                self.text_objects[spec][self.maxy].set_alpha(1.0)
                self.text_objects[spec][self.maxy].set_color('r')


            self.selected.append(self.maxy)
            # print(self.selected)

            print('woooop')
            print(len(self.selected))
            print(self.selected)
            self.draw_bores(None)
            self.canvas.draw()
            self.canvas.flush_events()
            # self.canvas2.flush_events()




    def findnear_index(self,test,array):
        #array = numpy.asarray(array)
        return (numpy.abs(array - test)).argmin()




    def draw_bores(self, event):
        if event != None:
            x = event.xdata
            y = event.ydata
        else:
            x = 0
            y = 0

        #print(self.parent.tabOne.molecule.spec['hnco'].index2)
        for ax_name in self.parent.tabOne.molecule.spec.keys():
            self.axes2[ax_name].set_ylim((0,0))
            # self.axes2[ax_name].spines['top'].set_visible(False)
            # self.axes2[ax_name].spines['right'].set_visible(False)
            # self.axes2[ax_name].spines['left'].set_visible(False)
        # self.axes2['ca'].set_ylim((0,0))
        
        # self.axes2['cb'].set_ylim((0,0))
        if x and y and len(self.selected) == 0:
            for i, spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
              
                x1, y1 = numpy.abs(self.parent.tabOne.molecule.spec[spec].index2-x).argmin(), numpy.abs(self.parent.tabOne.molecule.spec[spec].index1-y).argmin()
               
                self.bores[spec].set_ydata(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])
                self.axes2[spec].set_ylim(min(min(self.axes2[spec].get_ylim()), numpy.min(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])), max(max(self.axes2[spec].get_ylim()), numpy.max(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])))

        else:
            for i, spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
                # print(self.peak_data[spec][0][self.selected[0]])
                # print(self.peak_data[spec][0][self.selected[0]][0])
                
                x, y = self.peak_data[spec][0][self.selected[0]][0], self.peak_data[spec][0][self.selected[0]][1]
                x1, y1 = numpy.abs(self.parent.tabOne.molecule.spec[spec].index2-x).argmin(), numpy.abs(self.parent.tabOne.molecule.spec[spec].index1-y).argmin()
                # print(x1, y1)
                self.bores[spec].set_ydata(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])
                self.axes2[spec].set_ylim(min(min(self.axes2[spec].get_ylim()), numpy.min(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])), max(max(self.axes2[spec].get_ylim()), numpy.max(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])))


                # print()

                # self.bores[spec].set_ydata(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])
                # self.axes2[spec].set_ylim(min(min(self.axes2[spec].get_ylim()), numpy.min(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])), max(max(self.axes2[spec].get_ylim()), numpy.max(self.parent.tabOne.molecule.spec[spec].data[:,y1,x1])))
            
        self.canvas2.draw()

    def makeref(self):
        return
        self.ref={}
        self.ref[self.parent.tabOne.molecule.refSpec]=0,0

        self.move={}
        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
            if(spec!=self.parent.tabOne.molecule.refSpec):

                hdiff=[]
                ndiff=[]
                for pk in self.parent.tabOne.molecule.spec[self.parent.tabOne.molecule.refSpec].peak2D:
                    Nref=pk.y
                    Href=pk.x
                    tig=0
                    for pk2 in self.parent.tabOne.molecule.spec[spec].peak2D:
                        if pk2.name==pk.name:
                            tig=1
                            break
                    if(tig==1):
                        Nval=pk2.y  #nitrogen
                        Hval=pk2.x  #proton
                        hdiff.append(Hval-Href)
                        ndiff.append(Nval-Nref)

                        if(pk.name not in self.move.keys()):
                            self.move[pk.name]={}
                        if(spec not in self.move[pk.name].keys()):
                            self.move[pk.name][spec]={}

                        self.move[pk.name][spec]['N']=Nval-Nref
                        self.move[pk.name][spec]['H']=Hval-Href

                self.ref[spec]=numpy.median(hdiff),numpy.median(ndiff)
                
                for pk3 in self.move.keys():
                    if spec in self.move[pk3].keys():
                        self.move[pk3][spec]['N']-=self.ref[spec][1]
                        self.move[pk3][spec]['H']-=self.ref[spec][0]
                        if(numpy.fabs(self.move[pk3][spec]['N'])>0.4):
                            print('warning:',pk3,spec,'n',self.move[pk3][spec]['N'])
                            for pk in self.parent.tabOne.molecule.spec[self.parent.tabOne.molecule.refSpec].peak2D:
                                if(pk.name==pk3):
                                    print(pk.x,pk.y)
                        if(numpy.fabs(self.move[pk3][spec]['H'])>0.025):
                            print('warning',pk3,spec,'h',self.move[pk3][spec]['H'])
                            for pk in self.parent.tabOne.molecule.spec[self.parent.tabOne.molecule.refSpec].peak2D:
                                if(pk.name==pk3):
                                    print(pk.x,pk.y)

    def setup_figure_data(self):

        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
            refN=0 #self.ref[spec][1]
            #refH=self.ref[spec][0]
            refH=self.parent.tabOne.molecule.spec[spec].Hmed
            Xs=self.parent.tabOne.molecule.spec[spec].XX_proj-refN
            Ys=self.parent.tabOne.molecule.spec[spec].YY_proj-refH
            Zs=self.parent.tabOne.molecule.spec[spec].data_proj
            self.contour_data[spec] = ([Ys,Xs,Zs])

            ##labels:#
            referenced_peaks = []
            folded_peaks = []
            referenced_folded_peaks = []
            neither_peaks = []
            names = []
            for pk in self.parent.tabOne.molecule.spec[spec].peak2D:
                referenced_peaks.append([pk.x - refH, pk.ppmJ - refN])
                folded_peaks.append([pk.x - refH, pk.y])
                referenced_folded_peaks.append([pk.x - refH, pk.y])
                neither_peaks.append([pk.x - refH,pk.ppmJ])
                names.append(pk.name)

            referenced_peaks = numpy.array(referenced_peaks)
            folded_peaks = numpy.array(folded_peaks)
            referenced_folded_peaks = numpy.array(referenced_folded_peaks)
            neither_peaks = numpy.array(neither_peaks)
            self.peak_data[spec] = [referenced_peaks, folded_peaks, referenced_folded_peaks, neither_peaks, names]
            #print(self.peak_data[spec][1])
            #print('whack')

   


    # def update_text(self):
    def AddCircles(self,X,Y,rH,rN):
        for x,y in zip(X,Y): #for each centre...

            theta=numpy.linspace(0,2*numpy.pi,100)
            xVals=rH*numpy.cos(theta)+x
            yVals=rN*numpy.sin(theta)+y

            self.axes.plot(xVals,yVals,color='r')


    def draw_figure(self, bounds = [-1]):
        
        if self.first_draw == 0:
           x_min, x_max = self.axes.get_xlim()
           y_min, y_max = self.axes.get_ylim()
        self.axes.clear()


        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):  #for each spectrun...

            refN=0 #self.ref[spec][1]
            #refH=self.ref[spec][0]
            refH=self.parent.tabOne.molecule.spec[spec].Hmed

            if math.isnan(refN):
                refN = 0
            if math.isnan(refH):
                refH = 0

            if(self.cblist[spec].IsChecked()):
                Xs=self.parent.tabOne.molecule.spec[spec].XX_proj-refN
                Ys=self.parent.tabOne.molecule.spec[spec].YY_proj-refH
                Zs=self.parent.tabOne.molecule.spec[spec].data_proj


                levels=self.GetLevels(float(self.textbox1[spec].GetValue()),float(self.textbox2[spec].GetValue()),int(self.textbox3[spec].GetValue()))
                # print(Zs, numpy.max(Zs))

                #self.axes.contour(Ys-refH,Xs-refN,Zs,levels,cmap=cm.get_cmap(self.cmaps[self.comboClbox[spec].GetSelection()]),norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
                #print(Xs, Ys, refN, refH)
                # print(Ys-refH)
                self.axes.contour(Ys,Xs,Zs,levels,cmap=cm.get_cmap(self.cmaps[self.comboClbox[spec].GetSelection()]),norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network

                # self.canvas.draw()
                #print(spec,Xs.shape,Ys.shape,Zs.shape)




                #for pk in self.parent.tabOne.molecule.spec[spec].peak2D:
                    # print(pk.x, pk.ppmJ)
            if self.referenced.IsChecked():
                if self.folded.IsChecked():
                    #self.axes.scatter(pk.x - refH, pk.y, color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cblab[spec].IsChecked()):
                        self.axes.scatter(self.peak_data[spec][2][:,0],self.peak_data[spec][2][:,1], color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cbpk[spec].IsChecked()):
                        self.text_objects[spec] = []
                        for pk_number, text in enumerate(self.peak_data[spec][4]):
                            self.text_objects[spec].append(self.axes.text(self.peak_data[spec][2][pk_number,0],self.peak_data[spec][2][pk_number,1], text, color='k', fontsize=12))

                else:
                    #self.axes.scatter(pk.x - refH, pk.ppmJ - refN, color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cblab[spec].IsChecked()):
                        self.axes.scatter(self.peak_data[spec][0][:,0],self.peak_data[spec][0][:,1], color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cbpk[spec].IsChecked()):
                        self.text_objects[spec] = []
                        for pk_number, text in enumerate(self.peak_data[spec][4]):
                            self.text_objects[spec].append(self.axes.text(self.peak_data[spec][0][pk_number,0],self.peak_data[spec][0][pk_number,1], text, color='k', fontsize=12))
            else:
                if self.folded.IsChecked():
                    #self.axes.scatter(pk.x, pk.y, color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cblab[spec].IsChecked()):
                        self.axes.scatter(self.peak_data[spec][1][:,0],self.peak_data[spec][1][:,1], color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cbpk[spec].IsChecked()):
                        self.text_objects[spec] = []
                        for pk_number, text in enumerate(self.peak_data[spec][4]):
                            self.text_objects[spec].append(self.axes.text(self.peak_data[spec][1][pk_number,0],self.peak_data[spec][1][pk_number,1], text, color='k', fontsize=12))
                else:
                    #self.axes.scatter(pk.x,pk.ppmJ,color=self.ccyc[i],zorder=20,marker='x')
                    if(self.cblab[spec].IsChecked()):
                        self.axes.scatter(self.peak_data[spec][3][:,0],self.peak_data[spec][3][:,1], color=self.ccyc[i], zorder=20, marker='x')
                    if(self.cbpk[spec].IsChecked()):
                        self.text_objects[spec] = []
                        for pk_number, text in enumerate(self.peak_data[spec][4]):
                            self.text_objects[spec].append(self.axes.text(self.peak_data[spec][3][pk_number,0],self.peak_data[spec][3][pk_number,1], text, color='k', fontsize=12))
                        # self.axes.text(self.peak_data[spec][3][:,0],self.peak_data[spec][3][:,1], self.peak_data[spec][4], color='r', fontsize=8)





        if self.first_draw == 1:
           y_min,y_max=self.axes.get_ylim()
           x_min,x_max=self.axes.get_xlim()
           self.first_draw = 1

        #self.axes.set_xlim(x_min,x_max)
        #self.axes.set_ylim(y_min,y_max)
        self.axes.set_xlim(x_max,x_min)
        self.axes.set_ylim(y_max,y_min)

        self.control_text = self.axes.text(x_min, y_min, '')
        self.canvas.draw()

        # if self.add == True:

        pass


    # add a new peak to all 4 2d peaklists
    def add_peak(self, H, N):
        print('adding peak at', H, 'ppm in H and', N, 'ppm in N')
        # get spectra
        specs=list(self.parent.tabOne.molecule.spec.keys())
        # initiate counter
        counter=[]

        # search through the 2d lists for the highest peak and add one
        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
            pkNames = self.peak_data[spec][4]
            for pk in pkNames:
                count=re.sub(r'\D', '', pk)
                count=int(count)
                counter.append(count)
        newPkNo=max(counter) + 1
        name = str(newPkNo)+'H-N'


        # add the new peak to all the spectra
        for i,spec in enumerate(self.parent.tabOne.molecule.spec.keys()):
            #refN=self.ref[spec][1]
            #refH=self.ref[spec][0]
            refN=0
            refH=self.parent.tabOne.molecule.spec[spec].Hmed
            if self.referenced.IsChecked():
                refPeak = [[H, N]]
                unRef = [[H + refH, N + refN]]
                folded = refPeak
                refFold = refPeak
            else:
                refPeak = [[H - refH, N - refN]]
                unRef = [[H, N]]
                folded = refPeak
                refFold = refPeak
            refPeak = numpy.array(refPeak)
            unRef = numpy.array(unRef)
            folded = numpy.array(folded)
            refFold = numpy.array(refFold)
            nameArray = numpy.empty((0,1), str)
            print(refFold.shape)
            nameArray = numpy.append(nameArray, numpy.array(name))

            self.peak_data[spec][0] = numpy.append(self.peak_data[spec][0], refPeak, axis = 0)
            self.peak_data[spec][1] = numpy.append(self.peak_data[spec][1], folded, axis = 0)
            self.peak_data[spec][2] = numpy.append(self.peak_data[spec][2], refFold, axis = 0)
            self.peak_data[spec][3] = numpy.append(self.peak_data[spec][3], unRef, axis = 0)
            self.peak_data[spec][4] = numpy.append(self.peak_data[spec][4], nameArray, axis = 0)
            self.peak_data[spec][4] = self.peak_data[spec][4].tolist()










    def GetLevels(self,min_level,fac,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*fac)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels

    def onPeakList(self, event):
        bool=peakManAss(self)



class peakManAss(wx.App):
    def __init__(self,inherit):
        self.frame_peakManAssFrame=peakManAssFrame(None,10,'Peaks',inherit)
        self.frame_peakManAssFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class peakManAssFrame(wx.Frame):
    #    title = 'AssBox'
    def __init__(self, parent, id, title, inherit):
        self._init_ctrls(parent, inherit)
        self.parent = inherit

    def _init_ctrls(self, prnt, parent):
        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=prnt,
                          pos=wx.Point(358, 184), size=wx.Size(250, 20),
                          style=wx.DEFAULT_FRAME_STYLE, title=u'Peaks ...')
        self.SetClientSize(wx.Size(900, 280))

        panel = wx.Panel(self, -1)

        self.parent = parent.parent
        self.corrDict = {}

        self.lc = SortedListCtrl(panel, self.corrDict)

        cnt = 0
        self.lc.InsertColumn(cnt, 'Spectrum');
        cnt += 1
        self.lc.InsertColumn(cnt, 'ResID');
        cnt += 1
        self.lc.InsertColumn(cnt, 'Name');
        cnt += 1

        self.lc.InsertColumn(cnt, 'H');
        cnt += 1
        self.lc.InsertColumn(cnt, 'N');
        cnt += 1
        self.lc.InsertColumn(cnt, 'ppmI(ppm)');
        cnt += 1
        self.lc.InsertColumn(cnt, 'ppmJ(ppm)');
        cnt += 1
        self.lc.InsertColumn(cnt, 'ppmK(ppm)');
        cnt += 1


        # self.lc.SetColumnWidth(0, 140)
        # self.lc.SetColumnWidth(1, 153)

        self.Refreshbutton = wx.Button(panel, -1, 'Refresh', (710, 10))
        self.Showbutton = wx.Button(panel, -1, 'Show', (710, 10))
        self.Removebutton = wx.Button(panel, -1, 'Remove', (710, 60))
        self.Closebutton = wx.Button(panel, -1, 'Close', (710, 160))
        self.NextButton = wx.Button(panel, -1, 'Next', (710, 10))
        self.PrevButton = wx.Button(panel, -1, 'Previous', (710, 10))

        choices = []
        for pk in self.parent.tabOne.molecule.spec[self.parent.tabOne.molecule.refSpec].peak:
            choices.append(pk.name)

        self.spec_combo = wx.ComboBox(panel, choices=choices)

        # self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnShow, self.lc)

        self.Bind(wx.EVT_BUTTON, self.OnRefresh, self.Refreshbutton)
        self.Bind(wx.EVT_BUTTON, self.OnShow, self.Showbutton)
        self.Bind(wx.EVT_BUTTON, self.OnClose, self.Closebutton)
        self.Bind(wx.EVT_BUTTON, self.OnRemove, self.Removebutton)
        self.Bind(wx.EVT_BUTTON, self.OnNext, self.NextButton)
        self.Bind(wx.EVT_BUTTON, self.OnPrev, self.PrevButton)
        self.Bind(wx.EVT_COMBOBOX, self.OnRefresh, self.spec_combo)
        # self.vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.lc, 1, wx.EXPAND)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.Refreshbutton, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.Showbutton, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.Closebutton, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.Removebutton, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.spec_combo, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.PrevButton, wx.ALIGN_CENTER | wx.TOP)
        vbox.Add(self.NextButton, wx.ALIGN_CENTER | wx.TOP)
        hbox.Add(vbox)
        panel.SetSizer(hbox)
        self.spec_combo.SetSelection(0)
        evt = wx.PyCommandEvent(wx.EVT_BUTTON.typeId, self.Refreshbutton.GetId())
        wx.PostEvent(self, evt)
        # self.OnRefresh(True)
        # self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        # panel1.SetSizer(self.vbox)

        # hbox  = wx.BoxSizer(wx.HORIZONTAL)
        # hbox.Add(self.Addbutton, 1, wx.EXPAND)
        # hbox.Add(self.Removebutton, 1, wx.EXPAND)
        # hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        # hbox.Add(self.Closebutton, 1, wx.EXPAND)
        # hbox.Add(self.Savebutton, 1, wx.EXPAND)
        # self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        # panel1.SetSizer(self.vbox)

        self.Centre()
        self.Show(True)

        # hbox.Add(vbox2, 1, wx.EXPAND)
        # self.SetSizer(hbox)

    #        self.SetSizer(self.vbox)
    #        self.vbox.Fit(self)

    def onItemSelected(self, event):
        """"""

        currentItem = event.m_itemIndex

        # car = self.corrDict[currentItem]
        # print(car)

        # count = self.lc.GetItemCount()
        # self.sorted_artists = [self.list.GetItem(itemId=row, col=0).GetText() for row in xrange(count)]
        # print(self.sorted_artists)
        # print(self.sorted_artists[currentItem])

    def AtoI(self, val):
        for i in range(len(self.parent.tabOne.peak)):
            if (val == self.parent.tabOne.peak[i].name):
                return i

    """
    def OnAdd(self, event):
        sele=self.lc.GetFirstSelected()
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row, col=0).GetText() for row in xrange(count)][sele]
        col2 = [self.lc.GetItem(row, col=1).GetText() for row in xrange(count)][sele]
        #print(col1,col2,self.AtoI(col1))
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
        self.parent.peak.pop(index)
        self.parent.draw_figure()
        self.OnRefresh(True)

    def OnClose(self, event):
        self.Close()

    def OnRefresh(self, event):
        print('Refreshing list')
        self.lc.DeleteAllItems()
        corr = []
        self.corrDict = {}
        self.peak_name = self.spec_combo.GetValue()
        z = -1
        for i, spec in enumerate(self.parent.parent.tabOne.molecule.spec.keys()):

            # print(self.parent.parent.tabOne.molecule.spec[spec].peak)
            peak_numbers = []
            for j, pk in enumerate(self.parent.parent.tabOne.molecule.spec[spec].peak):
                if re.findall(r'\d+', pk.name)[0] == re.findall(r'\d+', self.peak_name)[0]:

                    peak_numbers.append(j)


            for j, peak_number in enumerate(peak_numbers):
                add = []
                z+=1
                pk = self.parent.parent.tabOne.molecule.spec[spec].peak[peak_number]
                add.append(spec)
                add.append(re.findall(r'\d+', pk.name)[0])
                add.append(pk.name)
                add.append(pk.f1)
                add.append(pk.f2)
                add.append(pk.f1)
                add.append(pk.f2)
                add.append(pk.f3)
                # if (self.parent.tabOne.dim == 4):
                #     add.append(pk.ppmL)

                self.corrDict[z] = add


        for key, data in self.corrDict.items():
            num_items = self.lc.GetItemCount()
            self.lc.InsertItem(num_items, (data[0]))  # add assignment
            self.lc.SetItem(num_items, 0, (data[0]))  # add atom
            self.lc.SetItem(num_items, 1, str(data[1]))  # add atom

            for i in range(len(data) - 2):
                self.lc.SetItem(num_items, i + 2, str(data[i + 2]))  # add atom
            self.lc.SetItemData(num_items, key)
        self.lc.Update(self.corrDict)

    def OnPrev(self, event):
        i = self.spec_combo.GetSelection()
        self.spec_combo.SetSelection(i-1)
        self.OnRefresh(True)

    def OnNext(self, event):
        i = self.spec_combo.GetSelection()
        self.spec_combo.SetSelection(i + 1)
        self.OnRefresh(True)

    def OnShow(self, event):
        # index = self.lc.GetFocusedItem()
        index = self.lc.GetFirstSelected()

        pk = self.parent.peak[index]
        print('Showing ', pk.name)
        print(pk.ppmI, pk.ppmJ)
        axes = self.parent.axes

        y_min = self.parent.YY[0][0]
        y_max = self.parent.YY[(len(self.parent.YY)) - 1][0]
        x_min = self.parent.XX[0][0]
        x_max = self.parent.XX[0][(len(self.parent.XX[0])) - 1]

        widY = (y_max - y_min) / 10.
        widX = (x_max - x_min) / 10.

        axes.set_xlim(pk.ppmI - widX, pk.ppmI + widX)
        axes.set_ylim(pk.ppmJ - widY, pk.ppmJ + widY)
        self.parent.draw_figure()

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
    def __init__(self, parent, dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        # ColumnSorterMixin.__init__(self, len(dicty.keys()))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

    def Update(self, dicty):
        ColumnSorterMixin.__init__(self, len(dicty.keys()))
        self.itemDataMap = dicty
        # print(dicty[0])

    def CustColumnSorter(self, key1, key2):
        col = self._col
        ascending = self._colSortFlag[col]
        ascending = 1
        item1 = self.itemDataMap[key1][col]
        item2 = self.itemDataMap[key2][col]

        self.num_cols = [0, 2, 3, ]
        if col in self.num_cols:
            # just convert them to float, cmp do comparing float well
            item1 = float(item1)
            item2 = float(item2)

        cmpVal = cmp(item1, item2)

        # If the items are equal then pick something else to make the sort value unique
        if cmpVal == 0:
            cmpVal = apply(cmp, self.GetSecondarySortValues(col, key1, key2))

        if ascending:
            return cmpVal
        else:
            return -cmpVal

    def GetColumnSorter(self):
        return self.CustColumnSorter
