#!/usr/bin/python
import wx,string,copy,math,numpy,os,sys
# import matplotlib.pyplot as plt            #import matplotlib
#matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import nmrglue as ng
from matplotlib.figure import Figure
from assign.assign_main import peakEntry

#from wx.lib.mixins.listctrl import ColumnSorterMixin
#from deconFrame import ParseFlt

class tocsyFrame(wx.Panel):

    def __init__(self,parent):

        #self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        """
        #wx.Frame.__init__(self, None, wx.ID_ANY,
                          "strip plots",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          size=(self.monitorWidth*0.95, self.monitorHeight*0.85),
                          #size=(1300,670)
                          )
        """

        wx.Panel.__init__(self, parent=parent)
        #panel = wx.Panel(self)

        #wx.Panel.__init__(self, parent=parent)
        self.parent=parent.tabOne.molecule

        

        #self.SetBackgroundColour('WHITE')
        #self.panel=wx.Panel(self,-1)
        #self.parent=mo
        #self.tabOne=parent.tabOne
        #self.sym=self.tabOne.cb_grid.IsChecked()
        #copy in the previous variables
        self.index_data=self.index(self.parent.spec['hnco'].peak2D) #inherit the data index
        # print(self.index_data)

        #self.thresh=parent.tabOne.dmax*float(parent.tabOne.threshBox.GetValue())             #inherit noise value
        #self.offset=copy.deepcopy(tabOne.offset)
        #self.offset=0.0                                #
        #self.peak2D=parent.tabOne.peak2D                   #inherit peak list
        #self.spectrumfile=parent.tabOne.spectrumfile   #inherit spectrumfile
        #get 2d strips from 3d data
        #self.GetSlice2d(self.spectrumfile)             #slice up the 2D spectrum

        #
        self.ax_reset0=1
        self.ax_reset1=1       #for keeping the zoom
        self.ax_reset2=1
        self.ax_reset3=1
        self.inc=0            #for incrementing the slices
        self.inc2=0

        self.SELECT=0
        self.ADD=0
        self.MOVE=0
        self.toggle2D=False
        self.toggleBMRB = False
        self.residue = 'A'


        #self.data=self.parent.spec['hnco'].data
        self.create_main_panel()
        # print('22')
        self.draw_figure()
        # print('44')
        #self.Show(True)
        self.SetSizerAndFit(self.main_sizer)

   


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.


        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)


        self.fig.clear()
        
        self.axes = self.fig.add_subplot(111)


        self.wipe() #set selection=0 for all axes

        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)



        self.text1= wx.StaticText(self, -1, 'Min:')
        self.text2= wx.StaticText(self, -1, 'Fac:')
        self.text3= wx.StaticText(self, -1, 'Num:')
        self.textbox1=wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox2=wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox3=wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)

        self.textbox1.SetValue(str(self.parent.spec['ctocsy'].noise))
        self.textbox2.SetValue(str(1.2))
        self.textbox3.SetValue(str(15))


        self.cmaps =  'seismic','bwr','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'

        self.textcolour= wx.StaticText(self, -1, 'Colour:')
        self.comboClbox=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.cmaps, style=wx.CB_READONLY)

        self.comboClbox.SetSelection(0)



        #print(self.parent.spec['hnco'].peak2D #MASTER PEAK LIST)
        #NEED TO ADD SOME ERROR CHECKING: LOOK AT ALL THE PEAK LISTS AND MAKE SURE THEY ARE ALL
        #THE SAME
        self.listy=[]
        for i in range(len(self.parent.spec['hnco'].peak2D)):
            self.listy.append(self.parent.spec['hnco'].peak2D[i].name)

        # Create a new list of peak labels, but highlight ones that require attention 
        # with a '*f' or '*b' if missing a forward or backward connection respectively

        self.checkForwardList = self.parent.checkForward
        self.checkBackwardList = self.parent.checkBackward

        self.list2 = []

        for peakLabel in self.listy:
            if(peakLabel in self.checkForwardList):
                peakLabel = peakLabel + '*f'
            if(peakLabel in self.checkBackwardList):
                peakLabel = peakLabel + '*b'
            self.list2.append(peakLabel)
        
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(150, -1), choices=self.list2, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox1)
        self.peak = self.listy[0]
    





        self.Nbutton = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)



        self.drawbutton = wx.Button(self, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.width1Box = wx.TextCtrl(self, size=(50, -1))
        self.width1Box.SetValue('0.1')


        self.SelectButton = wx.Button(self, -1, "Select")
        self.DeselectButton = wx.Button(self, -1, "DeSelect")
        self.MoveButton = wx.Button(self, -1, "Move")
        self.AddButton = wx.Button(self, -1, "Add")
        self.DeleteButton = wx.Button(self, -1, "Delete")


        self.SelectButton.Bind(wx.EVT_BUTTON,self.onSelectButton)
        self.DeselectButton.Bind(wx.EVT_BUTTON,self.onDeselectButton)
        self.MoveButton.Bind(wx.EVT_BUTTON,self.onMoveButton)
        self.AddButton.Bind(wx.EVT_BUTTON,self.onAddButton)
        self.DeleteButton.Bind(wx.EVT_BUTTON,self.onDeleteButton)



        self.toggle1D_or_2D_label= wx.StaticText(self, -1, 'Plot 2D:')
        self.toggle1D_or_2D = wx.CheckBox(self,-1,)
        self.Bind(wx.EVT_CHECKBOX, self.toggle_2D_checkbox, self.toggle1D_or_2D)
        self.toggle1D_or_2D.SetValue(False)


        self.toggle_bmrb_label = wx.StaticText(self, -1, 'BMRB Plot:')
        self.toggle_bmrb = wx.CheckBox(self,-1,)
        self.Bind(wx.EVT_CHECKBOX, self.toggle_bmrb_checkbox, self.toggle_bmrb)
        self.toggle_bmrb.SetValue(False)


        self.list_of_amino_acids = 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
        
        
        
        self.ComboBoxBMRB=wx.ComboBox(self, -1, pos=(620, 180), size=(150, -1), choices=self.list_of_amino_acids, style=wx.CB_READONLY)
        self.ComboBoxBMRB.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.draw_bmrb, self.ComboBoxBMRB)

        
        self.main_sizer=wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.Add(self.canvas, 10, flag=wx.GROW)

        self.horizontal_sizer1=wx.BoxSizer(wx.HORIZONTAL)

        self.CreateControlBox()
        self.horizontal_sizer1.Add(self.ControlBox)

        self.CreateControlBoxWidth()
        self.horizontal_sizer1.Add(self.ControlBoxWidth)

        self.CreateControlBox2D()
        self.horizontal_sizer1.Add(self.ControlBox2D)

        self.CreateBMRB_Box()
        self.horizontal_sizer1.Add(self.BMRBbox)

        self.main_sizer.Add(self.horizontal_sizer1)


    def toggle_2D_checkbox(self, event):
        if(self.toggle1D_or_2D.IsChecked()==True):
            self.toggle2D = True
        else:
            self.toggle2D = False
        self.draw_figure()
    
    def toggle_bmrb_checkbox(self, event):
        if(self.toggle_bmrb.IsChecked()==True):
            self.toggleBMRB = True
        else:
            self.toggleBMRB = False
        self.draw_figure()

    def draw_bmrb(self, event):
        self.residue = self.ComboBoxBMRB.GetValue()
        print('Plotting BMRB predicted shifts for residue %s' % self.residue)
        self.draw_figure()


    def CreateControlBox(self):
        self.controlLabel = wx.StaticBox(self,-1,'Controls:')
        self.ControlBox = wx.StaticBoxSizer(self.controlLabel, wx.HORIZONTAL)
        self.ControlBox.Add(self.ComboBox1)
        self.ControlBox.Add(self.drawbutton)
        self.ControlBox.Add(self.Nbutton)
        self.ControlBox.Add(self.Pbutton)

    def CreateControlBoxWidth(self):
        self.width1Lab = wx.StaticBox(self, label="Linewidth(ppm):")
        self.ControlBoxWidth = wx.StaticBoxSizer(self.width1Lab, wx.HORIZONTAL)
        self.ControlBoxWidth.Add(self.width1Box)

    
    def CreateControlBox2D(self):
        self.TwoDboxlab = wx.StaticBox(self, label='2D Controls:')
        self.ControlBox2D = wx.StaticBoxSizer(self.TwoDboxlab, wx.HORIZONTAL)
        self.ControlBox2D.Add(self.toggle1D_or_2D_label)
        self.ControlBox2D.Add(self.toggle1D_or_2D)
        self.ControlBox2D.Add(self.text1)
        self.ControlBox2D.Add(self.textbox1)
        self.ControlBox2D.Add(self.text2)
        self.ControlBox2D.Add(self.textbox2)
        self.ControlBox2D.Add(self.text3)
        self.ControlBox2D.Add(self.textbox3)
        self.ControlBox2D.Add(self.textcolour)
        self.ControlBox2D.Add(self.comboClbox)
        self.ControlBox2D.Add(self.AddButton)



    def CreateBMRB_Box(self):
        self.BMRBlabel = wx.StaticBox(self, label='BMRB:')
        self.BMRBbox = wx.StaticBoxSizer(self.BMRBlabel, wx.HORIZONTAL)
        self.BMRBbox.Add(self.toggle_bmrb_label)
        self.BMRBbox.Add(self.toggle_bmrb)
        self.BMRBbox.Add(self.ComboBoxBMRB)
        


    

        




    def onSelectButton(self,event):
        print('Click to select')
        self.SELECT=1
        pass

    def wipe(self):
        self.axes.selection=[]

    def onDeselectButton(self,event):
        self.wipe()
        self.draw_figure()
        pass

    def onMoveButton(self,event):
        self.MOVE=1
        pass
    def onAddButton(self,event):
        print('Click to add')
        self.ADD=1
        pass
    def onDeleteButton(self,event):
        print('Removing selection')
        for axes in self.axes:
            if(len(axes.selection)>0):
                selection=sorted(axes.selection,reverse=True)
                spec=axes.spec
                res=axes.res
                for sele in axes.selection:
                    print('removing:',self.parent.peak[res][spec][sele].name)
                    self.parent.peak[res][spec].pop(sele)
                axes.selection=[]
        self.draw_figure()


    #make an index
    def index(self,array):
        index=[]
        for i in range(len(array)):
            index.append((array[i].name))
        return index

    def AtoI(self,val):
        for i,ind in enumerate(self.index_data):
            if(ind==val):
                return i
        return -1

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=string.split(line)
            peak.append(linetosave)
        peakfile.close()
        return peak

    def findnear_index(self,test,array):
        #array = numpy.asarray(array)
        idx = (numpy.abs(array - test)).argmin()
        return idx
        #return array[idx]


    def GetLevels(self,min_level,fac,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*fac)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels


    def PlotBMRB(self, residue):
        #self.parent.bmrb[atom][typ].append((resn,shift,std))

        self.p1to3 = self.parent.p1to3
        residue = self.p1to3[residue]



        import scipy.stats as stats

        for tp,list_of_entries in self.parent.bmrb['C'].items():
            if(tp!='C'):
                for entry in list_of_entries:
                    if(entry[0]==residue and entry[1]<90.0): # don't plot bmrb statistics for atoms that are out of the ctocsy range
                        sigma = numpy.float(entry[2])
                        mu = numpy.float(entry[1])
                        x = numpy.linspace(mu - 3*sigma, mu + 3*sigma, 1000)
                        gaussian = stats.norm.pdf(x, mu, sigma)
                        self.axes.plot(x,numpy.sqrt(2*numpy.pi*sigma)*gaussian*numpy.max(self.Zs)/2.0, color='black')
                        self.axes.text(mu,0,tp)


                

    #get 2d strips from 3d data
    def ReSlice2d(self,arr,inc,pkl,peak,width,orth=0,lab='hnco'):

        if(orth==0):
            #print(out 2D slice for each peak correlation)
            #print("Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width)
            print(pkl, lab, self.parent.spec[lab].pkIdx, pkl)
            ptC=self.parent.spec[lab].pkIdx[pkl][0]
            ptC=ptC+inc
            ptH=self.parent.spec[lab].pkIdx[pkl][1]
            print(ptH, ptC)


            ptH_max=self.findnear_index(float(peak[pkl].x)+float(width)/2,self.parent.spec[lab].index2)#find the nearest point to desired chemical shift in carbon index
            ptH_min=self.findnear_index(float(peak[pkl].x)-float(width)/2,self.parent.spec[lab].index2)#find the nearest point to desired chemical shift in carbon index

            Xs=self.parent.spec[lab].XX[:,ptC,ptH_max:ptH_min].transpose()
            Ys=self.parent.spec[lab].ZZ[:,ptC,ptH_max:ptH_min].transpose()

            Zs=arr[:,ptC,ptH_max:ptH_min].transpose() #extract the relevant 2d slice
        else:
            #print(out 2D slice for each peak correlation)
            #print("Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width)

            #ptC=self.findnear_index(float(peak[pkl][2]),self.index1)#find the nearest point to desired chemical shift in carbon index
            ptC_max=self.findnear_index(float(peak[pkl].y)+float(width)/2.,self.parent.spec[lab].index1)#find the nearest point to desired chemical shift in carbon index
            ptC_min=self.findnear_index(float(peak[pkl].y)-float(width)/2.,self.parent.spec[lab].index1)#find the nearest point to desired chemical shift in carbon index
            ptH=self.parent.spec[lab].pkIdx[pkl][1]
            #ptH=self.findnear_index(float(peak[pkl][2]),self.index2)#find the nearest point to desired chemical shift in carbon index
            ptH=ptH+inc
            Xs=self.parent.spec[lab].XX[:,ptC_max:ptC_min,ptH].transpose()
            Ys=self.parent.spec[lab].YY[:,ptC_max:ptC_min,ptH].transpose()
            Zs=arr[:,ptC_max:ptC_min,ptH].transpose() #extract the relevant 2d slice

        #print('Done!')
        return Xs,Ys,Zs


        #get 1d strips from 3d data
    def ReSlice1d(self,arr,inc,pkl,lab='hnco'):

        
        #print(out 2D slice for each peak correlation)
        #print("Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width)
      
        ptC=self.parent.spec[lab].pkIdx[pkl][0]
        ptC=ptC+inc
        ptH=self.parent.spec[lab].pkIdx[pkl][1]
        


            
        Xs=self.parent.spec[lab].XX[:,ptC,ptH]
        

        Zs=arr[:,ptC,ptH] #extract the relevant 2d slice
 
        self.Zs = Zs

        return Xs,Zs

    
    def Draw1D(self,sele, spec, xsub, ysub):
        if(spec != 'ctocsy'):
            print('cannot find ',spec,'in', self.specs)
            return

        #ADD THE SLICE
        

        axes=self.axes

        axes.pk=sele
        self.sele=sele
        axes.xsub=xsub
        axes.spec=spec
        axes.peaks=[]
        axes.res=self.index_data[sele]
        try:
            Xs,Zs=self.ReSlice1d(self.parent.spec[spec].data,self.inc2,sele,lab=spec)
        except Exception as e:
            print(e)



        self.axes.plot(Xs,Zs, linewidth = 0.5, color='#d62728')
        self.axes.set_xlim([numpy.min(Xs),numpy.max(Xs)])
        self.axes.invert_xaxis()


        try:
            for peak in self.parent.peak[self.sele2.split('*')[0]]['ctocsy']:
                self.axes.axvline(peak.f3, color='black')
        except:
            pass


    def AddSlice(self,sele,spec,xsub,ysub):
        if(spec != 'ctocsy'):
            print('cannot find ',spec,'in', self.specs)
            return

        #ADD THE SLICE
        
        levels=self.GetLevels(float(self.textbox1.GetValue()),float(self.textbox2.GetValue()),int(self.textbox3.GetValue()))

        axes=self.axes

        axes.pk=sele
        axes.xsub=xsub
        axes.spec=spec
        axes.peaks=[]
        axes.res=self.index_data[sele]
        try:
            Xs,Ys,Zs=self.ReSlice2d(self.parent.spec[spec].data,self.inc2,sele,self.parent.spec[spec].peak2D,self.Width1,lab=spec)
        except Exception as e:
            print(e)
        yave=numpy.average(Ys)


        if(ysub!=0):
            ysa=yave
        else:
            ysa=0


        clbox=self.comboClbox.GetSelection()

        cmap=cm.get_cmap(self.cmaps[clbox])

        self.axes.contour(Xs-xsub, Ys-ysa+ysub, Zs,levels,cmap=cmap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
        
        try:
            for peak in self.parent.peak[self.sele2.split('*')[0]]['ctocsy']:
                self.axes.scatter(peak.f3,peak.f1,marker='x',color='black',linewidth=2.0)
        except:
            pass

        y_max2a=Ys[0][0]
        y_min2a=Ys[(len(Ys))-1][0]
        x_max2a=Xs[0][0]
        x_min2a=Xs[0][(len(Xs[0]))-1]

        if(ysub==0):
            yline=yave
        else:
            yline=ysub
        
        self.grid_cb=True

        if(self.grid_cb):#horizontal line
            xl=(x_max2a-xsub,x_min2a-xsub)
            hl=(yline,yline)
            self.axes.plot(xl,hl,'k',zorder=0, lw=0.5)
        
        self.axes.invert_xaxis()



        return yave




    def draw_figure(self):
        """ Redraws the figure
        """
        if(self.ax_reset0==0):
            self.x_min0,self.x_max0=self.axes0.get_xlim()
            self.y_min0,self.y_max0=self.axes0.get_ylim()
        if(self.ax_reset1==0):
            self.x_min1,self.x_max1=self.axes1.get_xlim()
            self.y_min1,self.y_max1=self.axes1.get_ylim()
        if(self.ax_reset2==0):
            self.x_min2,self.x_max2=self.axes2.get_xlim()
            self.y_min2,self.y_max2=self.axes2.get_ylim()



        
        self.axes.clear()


        if(self.toggle2D==False):
            self.axes.set_title('cccTOCSY 1D Bore')
            self.axes.set_ylabel(r'Intensity',fontsize=12)
            self.axes.set_xlabel(r'$\delta_{C}$/ppm', fontsize=12)
        else:
            self.axes.set_title('cccTOCSY 2D Bore')
            self.axes.set_ylabel(r'$\delta_{H}$/ppm',fontsize=12)
            self.axes.set_xlabel(r'$\delta_{C}$/ppm', fontsize=12)




        
    

                
 
            

        # self.Width1=float(self.width1Box.GetValue())
        # if(self.Width1==0):
        #     self.Width1=1
        #     self.width1Box.SetValue(str(1))

        self.sele1=self.ComboBox1.GetSelection()
        self.sele2=self.ComboBox1.GetStringSelection()


        self.Width1=float(self.width1Box.GetValue())
        if(self.Width1==0):
            self.Width1=1
            self.width1Box.SetValue(str(1))


        # self.grid_cb=self.cb_grid_auto.GetValue()



        if(self.toggle2D==False):
            self.Draw1D(self.sele1, 'ctocsy',0,0)
        else:
            self.AddSlice(self.sele1,'ctocsy',0,0)


        if(self.toggle2D==False and self.toggleBMRB==True):
            self.PlotBMRB(self.residue)




        # self.axes0.set_yticklabels([])
        # self.axes1.set_yticklabels([])
        # self.axes3.set_yticklabels([])
        # self.axes4.set_yticklabels([])
        # self.axes6.set_yticklabels([])
        # self.axes7.set_yticklabels([])
        self.fig.tight_layout()
        self.canvas.draw()

        return





    def on_cb_grid(self, event):
        self.draw_figure()

    def on_cb_decon(self, event):
        if(self.parent.spec['hnco'].DECON==0):
            self.parent.spec['hnco'].OnButtonAnalyse(True)
            if(self.parent.spec['hnco'].DECON==0):
                print('No deconvolved calculation')
                self.cb_decon.SetValue(0)
                return
            else:
                self.parent.spec['hnco'].conn_data=self.parent.spec['hnco'].conn_data
                #print(self.parent.spec['hnco'].conn_data)
        self.datadec=self.parent.spec['hnco'].datadec
        self.inc=0
        self.inc2=0
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.draw_figure()

    def on_cb_grid_auto(self, event):
        self.draw_figure()



    def on_draw_button(self, event):
        self.peak = self.listy[self.ComboBox1.GetSelection()]

        self.inc=0
        self.inc2=0
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.ax_resetCC=1
        self.ax_resetHC=1
        self.draw_figure()

    def on_N_button(self, event):
        self.wipe()
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button(self, event):
        self.wipe()
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.selection=[]
        self.draw_figure()


    def on_bigN_button(self, event):
        self.wipe()
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.ComboBox3.SetSelection(self.ComboBox3.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_bigP_button(self, event):
        self.wipe()
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()-1)
        self.ComboBox3.SetSelection(self.ComboBox3.GetSelection()-1)
        self.selection=[]
        self.draw_figure()



    def on_swap_button(self, event):
        tmp=self.ComboBox1.GetSelection()
        self.ComboBox1.SetSelection(self.ComboBox2.GetSelection())
        self.ComboBox2.SetSelection(tmp)
        self.draw_figure()

    def on_Up_button(self, event):
        self.inc=self.inc+1
        self.draw_figure()

    def on_Down_button(self, event):
        self.inc=self.inc-1
        self.draw_figure()

    def on_NOE_button(self, event):
        bool=AssMan(self)


    def on_N_button2(self, event):
        self.wipe()
        self.ax_reset1=1
        self.ax_reset0=1
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button2(self, event):
        self.wipe()
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()-1)
        self.selection=[]
        self.draw_figure()

    def on_N_button3(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox3.SetSelection(self.ComboBox3.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button3(self, event):
        self.wipe()
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox3.SetSelection(self.ComboBox3.GetSelection()-1)
        self.selection=[]
        self.draw_figure()


    #def on_Up_button2(self, event):
    #    self.xmin2,self.xmax2=plt.xlim()
    #    self.ymin2,self.ymax2=plt.ylim()
    #    self.inc2=self.inc2+1
    #    self.selection=[]
    #    self.draw_figure()

    #def on_Down_button2(self, event):
    #    self.xmin2,self.xmax2=plt.xlim()
    #    self.ymin2,self.ymax2=plt.ylim()
    #    self.inc2=self.inc2-1
    #    self.selection=[]
    #    self.draw_figure()




    def on_deselect_button(self, event):
        print('Clearing selection')
        self.selection=[]
        self.draw_figure()

    def on_delete_button(self, event):
        print('Removing selection')
        for axes in self.axes:
            spec=axes.spec
            res=axes.res
            if(len(axes.selection)>0):
                selection=sorted(axes.selection,reverse=True)
                for sele in self.selection:
                    print('removing:',self.parent.peak[res][spec][sele].name)
                    self.parent.peak[res][spec].pop(sele)
                axes.selection=[]
        self.draw_figure()



    def on_select_button(self,event):
        print('Click on panel 1 to select')
        self.selection=[]
        self.SELECT=1

    def on_add_button(self,event):
        print('Click in panel 1 to initiate search')
        self.SEARCH=1


    def on_search_button(self, event):
        self.thresh=float(self.textbox0.GetValue())
        self.pick_fac=float(self.textbox_pickfac.GetValue())

        self.res=self.ComboBox1.GetSelection()
        self.parent.spec['hnco'].conn_data=analslices1d_spec(self.res,self.selection,self.parent.spec['hnco'].peak2D,self.index_data,self.parent.spec['hnco'].conn_data,self.thresh,self.pick_fac)
        self.draw_figure()

    def on_save_button(self, event):
        self.outfile=self.textbox_savelist.GetValue()
        print('Saving ',len(self.parent.spec['hnco'].conn_data),'entries in connectivity table to ',self.outfile)
        outy=open(self.outfile,'w')
        for i in range(len(self.parent.spec['hnco'].conn_data)):
            outy.write('%s\t%s\n' % (self.index_data[self.parent.spec['hnco'].conn_data[i][0]],self.index_data[self.parent.spec['hnco'].conn_data[i][1]]))
        outy.close()


    def on_launch(self,event):
        # print(event)
        from decon.decon_tab import MyApp
        deconParFile='hnco/deconParFile'
        frame = MyApp(deconParFile)

    def on_load_button(self, event):
        self.outfile=self.textbox_savelist.GetValue()
        if(os.path.exists(self.outfile)==0):
            print('No outputfile of this type is avaiable')
            return

        print('Loading connectivity table from ',self.outfile)
        self.parent.spec['hnco'].corrFile=self.outfile
        self.parent.spec['hnco'].OnButtonAnalyse(True)

        # print(self.parent.spec['hnco'].noeTags)
        self.ComboBoxNoe.SetItems(self.parent.spec['hnco'].noeTags)
        if(len(self.parent.spec['hnco'].noeTags)!=0):
            self.ComboBoxNoe.SetSelection(0)

        if(self.parent.spec['hnco'].DECON==1):
            self.DECON=1
            self.datadec=self.parent.spec['hnco'].datadec
            self.parent.spec['hnco'].conn_data=self.parent.spec['hnco'].conn_data
            self.noeTag=self.parent.spec['hnco'].noeTags
        else:
            self.DECON=0


        """
        from matplotlib.backends.backend_pdf import PdfPages
        pdf=PdfPages('multipage_pdf.pdf')
        for i in range(len(self.noeTags)):

            self.ComboBox1.SetSelection(self.parent.spec['hnco'].conn_data[i][0])
            self.ComboBox2.SetSelection(self.parent.spec['hnco'].conn_data[i][1])
            self.inc=0
            self.ax_reset=1
            self.ax_reset2=1
            self.inc2=0
            self.draw_figure()
            pdf.savefig(self.fig)
            if(i==0):
                break
        pdf.close()
        """



    def on_comboBoxNoe_select(self, event):
        val=self.ComboBoxNoe.GetSelection()

        self.ComboBox1.SetSelection(self.parent.spec['hnco'].conn_data[val].v1)
        self.ComboBox2.SetSelection(self.parent.spec['hnco'].conn_data[val].v2)

        self.inc=0
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.inc2=0
        self.draw_figure()




    #when search button is pressed make selection
    def on_pick(self, event):


        axes=event.inaxes



        if(self.SELECT==1):
            x_min,x_max=axes.get_xlim()
            y_min,y_max=axes.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            argy=[]
            cbs1=self.ComboBox1.GetSelection()


            axes.selection=[]
            for i,peak in enumerate(axes.peaks):
                rad2=((event.xdata-peak[0])/xdist)**2.+((event.ydata-peak[1])/ydist)**2.
                raddy.append(rad2)
                argy.append(i)

            raddy=numpy.array(raddy)
            maxy=argy[numpy.argmin(raddy)]
            axes.selection.append(maxy)

            axes.selectspec=axes.spec
            print('selected:',axes.spec,axes.res,axes.pk,self.parent.peak[axes.res][axes.spec][maxy].name)


            self.SELECT=0
            self.draw_figure()




        if(self.ADD==1):
            pk=axes.pk
            spec=axes.spec
            ptN=self.parent.spec[spec].pkIdx[pk][0] #point value for N
            ptH=self.parent.spec[spec].pkIdx[pk][1] #point value for H
            i=self.findnear_index(self.parent.spec[spec].index0-axes.xsub,event.xdata) #find X location
            inty=self.parent.spec[spec].data[i,ptN,ptH]
            cppm=self.parent.spec[spec].index0[i]
            nppm=self.parent.spec[spec].index1[ptN]
            hppm=self.parent.spec[spec].index2[ptH]

            cppm-=axes.xsub

            print('position selected in original spectrum:',cppm,nppm,hppm,inty)


            if(spec not in self.parent.peak[axes.res].keys()):
                self.parent.peak[axes.res][spec]=[]
            #pktest=self.parent.peak[axes.res][spec][0]
            #print(pktest.name)
            #print(pktest.f1)
            #print(pktest.f2)
            #print(pktest.f3)
            cnt=len(self.parent.peak[axes.res][spec])+1
            while(1==1): #test names to make sure we get a new one
                test=axes.res+'_'+str(cnt)
                tig=0
                for pk2 in self.parent.peak[axes.res][spec]:
                    if(pk2.name==test):
                        tig=1
                        break
                if(tig==0):
                    break
                cnt+=1

            stry=('%s\t%f\t%f\t%f\t%e\t') % (test,hppm,nppm,cppm,inty)
            inst=peakEntry(stry.split())
            print('Added:',inst.name,inst.f1,inst.f2,inst.f3,inst.inty,'to',axes.res,axes.spec)
            self.parent.peak[axes.res][spec].append(inst)

            pk3=self.parent.peak[axes.res][spec][-1]
            if(spec in ('hnco','hncaco','hnca','hncoca')):
                if(pk3.inty<0):
                    self.parent.DoAlias(pk3,spec)
            axes.selection=[]
            axes.selection.append(len(self.parent.peak[axes.res][spec])-1)
            axes.selectspec=spec
            self.draw_figure()
            self.ADD=0



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
            self.canvas.print_figure(path, dpi=2000)
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

"""
class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(dicty.keys()))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

    def Update(self,dicty):
        ColumnSorterMixin.__init__(self,len(dicty.keys()))
        self.itemDataMap = dicty
        print(dicty[0])
"""
