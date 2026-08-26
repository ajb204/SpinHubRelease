#!/usr/bin/python
import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
# matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin


##########################################################################
# 2D plotting of NMR slices
#

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

#def RunFrame(molecule):



class AssFrameMan(wx.App):
    def __init__(self,inherit):
        app = wx.App()
        #frame = AssFrame(molecule)
        self.frame_ProcessFrame=AssFrame(inherit)
        #FGA added
        self.frame_ProcessFrame.Centre(direction=wx.BOTH)
        self.frame_ProcessFrame.Show(True)
        app.MainLoop()

#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class AssFrame(wx.Frame):


    def __init__(self,parent):

        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "strip plots",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          #size=(self.monitorWidth*0.95, self.monitorHeight*0.85),
                          size=(870,870)
                          )
        panel = wx.Panel(self)

        #wx.Panel.__init__(self, parent=parent)
        self.parent=parent
        #self.WXV=int(wx.__version__.split('.')[0])
        """
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title=u'Get 2D peak list ...')
        self.SetClientSize(wx.Size(900, 280))
        """

        #FGA changed
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.Point(258, 184), size=wx.Size(800, 480),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'MAGMA results ...')
        #self.SetClientSize(wx.Size(900, 280))
        #monitorWidth, monitorHeight = wx.GetDisplaySize()
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.DefaultPosition, size=(monitorWidth, monitorHeight),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'Split Plots')

        #self.SetBackgroundColour('WHITE')
        #self.panel=wx.Panel(self,-1)
        #self.parent=mo
        #self.tabOne=parent.tabOne
        #self.sym=self.tabOne.cb_grid.IsChecked()
        #copy in the previous variables
        self.index_data=self.index(self.parent.spec['hnco'].peak2D) #inherit the data index
        #self.thresh=parent.tabOne.dmax*float(parent.tabOne.threshBox.GetValue())             #inherit noise value
        #self.offset=copy.deepcopy(tabOne.offset)
        #self.offset=0.0                                #
        #self.peak2D=parent.tabOne.peak2D                   #inherit peak list
        #self.spectrumfile=parent.tabOne.spectrumfile   #inherit spectrumfile
        #get 2d strips from 3d data
        #self.GetSlice2d(self.spectrumfile)             #slice up the 2D spectrum

        self.GetSlice2d()             #slice up the 2D spectrum

        #
        self.ax_reset0=1
        self.ax_reset1=1       #for keeping the zoom
        self.ax_reset2=1
        self.inc=0            #for incrementing the slices
        self.inc2=0
        self.SELECT=0
        self.SEARCH=0
        self.selection=[]

        #self.data=self.parent.spec['hnco'].data

        self.create_main_panel()
        print('22')
        self.draw_figure()
        print('44')
        self.Show(True)
        self.Fit()


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

        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)


        self.cmaps =  'seismic','bwr','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'


        #[('Perceptually Uniform Sequential' [
        #    'viridis', 'plasma', 'inferno', 'magma']),
        # ('Sequential', [
        #    'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        #    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        #    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
        # ('Sequential (2)', [
        #    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        #    'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        #    'hot', 'afmhot', 'gist_heat', 'copper']),
        # ('Diverging', [
        #    'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
        #    'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
        # ('Qualitative', [
        #    'Pastel1', 'Pastel2', 'Paired', 'Accent',
        #    'Dark2', 'Set1', 'Set2', 'Set3',
        #    'tab10', 'tab20', 'tab20b', 'tab20c']),
        # ('Miscellaneous', [
        #    'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        #    'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
        #    'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]


        self.width1Lab = wx.StaticText(self, label="Width(ppm):")
        self.width1Box = wx.TextCtrl(self, size=(50, -1))
        self.width1Box.SetValue('0.1')

        #self.width2Lab = wx.StaticText(self, label="Width(ppm):")

        #self.width2Box = wx.TextCtrl(self, size=(50, -1))
        #self.width1Box.SetValue(str(ParseFlt(self.tabOne.deconParFile,'widthX')))
        #self.width2Box.SetValue(str(ParseFlt(self.tabOne.deconParFile,'widthY')))



        #min max and lvls for slices
        #self.text_slice=wx.StaticText(self, -1, 'Slices:',size=(80,-1))

        self.proj={}
        self.text1={}
        self.text2={}
        self.text3={}
        self.textbox1={}
        self.textbox2={}
        self.textbox3={}
        self.cblab={}
        self.cblist={}
        self.cntrLbl={}
        self.cntrSizer={}
        self.comboClbox={}
        self.launch={}

        flags = wx.ALIGN_LEFT | wx.ALL #| wx.ALIGN_CENTER_VERTICAL
        self.thresh=1E5
        self.vboxCntr=wx.BoxSizer(wx.VERTICAL)
        for i,spec in enumerate(('hnco','hncaco','hnca','hncoca','hncacb','hncocacb')):
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

            self.comboClbox[spec].SetSelection(i%2)


            self.cblab[spec]=wx.CheckBox(self, -1,)
            self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cblab[spec])
            self.cblab[spec].SetValue(True)

            self.launch[spec]=wx.Button(self, -1, "go",size=(30,-1))
            self.Bind(wx.EVT_BUTTON, self.on_launch, self.launch[spec])


            self.textbox1[spec].SetValue(str(self.parent.spec[spec].noise))
            self.textbox2[spec].SetValue(str(1.2))
            self.textbox3[spec].SetValue(str(15))
            self.cntrSizer[spec] = wx.BoxSizer(wx.HORIZONTAL)
            self.cntrSizer[spec].Add(wx.StaticText(self,-1,spec,size=(65,-1)))
            self.cntrSizer[spec].Add(self.cblist[spec])
            self.cntrSizer[spec].Add(self.cblab[spec])
            self.cntrSizer[spec].Add(self.text1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.comboClbox[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.launch[spec], 0, border=3, flag=flags)
            self.vboxCntr.Add(self.cntrSizer[spec])


        #print(self.parent.spec['hnco'].peak2D #MASTER PEAK LIST)
        #NEED TO ADD SOME ERROR CHECKING: LOOK AT ALL THE PEAK LISTS AND MAKE SURE THEY ARE ALL
        #THE SAME
        listy=[]
        for i in range(len(self.parent.spec['hnco'].peak2D)):
            listy.append(self.parent.spec['hnco'].peak2D[i].name)
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox1)

        self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox2.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox2)

        self.ComboBox3=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox3.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox3)

        #noes=[]
        #for i in range(len(self.tabOne.peak2D)):
        #    listy.append(self.tabOne.peak2D[i][0])
        #self.ComboBoxNoe=wx.ComboBox(self, -1, pos=(620, 180), size=(300, -1), choices=noes, style=wx.CB_READONLY)
        #self.ComboBoxNoe.SetSelection(0)
        #self.Bind(wx.EVT_COMBOBOX, self.on_comboBoxNoe_select, self.ComboBoxNoe)
        #if(self.textbox_savelist.GetValue()!=''):
        #    self.on_load_button(True)

        #self.swapbutton = wx.Button(self, -1,"Swap")
        #self.Bind(wx.EVT_BUTTON, self.on_swap_button, self.swapbutton)


        self.Nbutton = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        self.Nbutton2 = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button2, self.Nbutton2)

        self.Pbutton2 = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button2, self.Pbutton2)

        self.Nbutton3 = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button3, self.Nbutton3)

        self.Pbutton3 = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button3, self.Pbutton3)


        #self.Upbutton = wx.Button(self, -1,"+",size=(30,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_Up_button, self.Upbutton)

        #self.Downbutton = wx.Button(self, -1,"-",size=(30,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_Down_button, self.Downbutton)

        #self.NOEbutton = wx.Button(self, -1,"Peaks")
        #self.Bind(wx.EVT_BUTTON, self.on_NOE_button, self.NOEbutton)

        #self.Upbutton2 = wx.Button(self, -1,"Up")
        #self.Bind(wx.EVT_BUTTON, self.on_Up_button2, self.Upbutton2)
        #
        #self.Downbutton2 = wx.Button(self, -1,"Down")
        #self.Bind(wx.EVT_BUTTON, self.on_Down_button2, self.Downbutton2)


        self.drawbutton = wx.Button(self, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.bigNextbutton = wx.Button(self, -1, "NEXT")
        self.Bind(wx.EVT_BUTTON, self.on_bigN_button, self.bigNextbutton)
        self.bigPrevbutton = wx.Button(self, -1, "PREVIOUS")
        self.Bind(wx.EVT_BUTTON, self.on_bigP_button, self.bigPrevbutton)

        #self.selectbutton = wx.Button(self, -1, "Select")
        #self.deselectbutton = wx.Button(self, -1, "Deselect")
        #self.deletebutton = wx.Button(self, -1, "Delete")
        #self.addbutton = wx.Button(self, -1, "Add")

        #self.Bind(wx.EVT_BUTTON, self.on_select_button, self.selectbutton)
        #self.Bind(wx.EVT_BUTTON, self.on_deselect_button, self.deselectbutton)
        #self.Bind(wx.EVT_BUTTON, self.on_delete_button, self.deletebutton)
        #self.Bind(wx.EVT_BUTTON, self.on_add_button, self.addbutton)

        self.cb_flip = wx.CheckBox(self, -1,
            "Orth",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_decon, self.cb_flip)


        self.cb_decon = wx.CheckBox(self, -1,
            "Decon",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_decon, self.cb_decon)

        self.cb_grid_auto = wx.CheckBox(self, -1,
            "Labels",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid_auto, self.cb_grid_auto)
        self.cb_grid_auto.SetValue(1)




        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)

        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)






        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.vboxC=wx.BoxSizer(wx.VERTICAL)
        self.vboxC.Add(self.drawbutton, 0, border=3, flag=flags)
        #self.vboxC.Add(self.NOEbutton, 0, border=3, flag=flags)
        self.vboxC.Add(self.cb_flip)
        self.vboxC.Add(self.cb_decon)
        self.vboxC.Add(self.cb_grid_auto)
        self.vboxC.Add(self.width1Lab,0,border=3,flag=flags)
        self.vboxC.Add(self.width1Box,0,border=3,flag=flags)
        self.vboxC.Add(self.bigPrevbutton, 0, border=3, flag=flags)
        self.vboxC.Add(self.bigNextbutton, 0, border=3, flag=flags)

        #self.vboxC.Add(self.cb_grid_auto)


        self.hbox.Add(self.vboxC)

        self.vboxCombo=wx.BoxSizer(wx.VERTICAL)

        self.topLbl = wx.StaticBox(self,-1,'Top:')
        self.topSizer = wx.StaticBoxSizer(self.topLbl, wx.VERTICAL)



        self.hbox1=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.ComboBox1, 0, border=3, flag=flags)
        self.hbox1.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox1.Add(self.Nbutton, 0, border=3, flag=flags)
        #self.hbox1a=wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox1a.Add(self.Upbutton, 0, border=3, flag=flags)
        #self.hbox1a.Add(self.Downbutton, 0, border=3, flag=flags)

        self.topSizer.Add(self.hbox1)
        #self.leftSizer.Add(self.hbox1a)

        print('1')


        self.middleLbl = wx.StaticBox(self,-1,'middle:')
        self.middleSizer = wx.StaticBoxSizer(self.topLbl, wx.VERTICAL)
        self.hbox2=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.ComboBox2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Pbutton2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Nbutton2, 0, border=3, flag=flags)
        self.middleSizer.Add(self.hbox2)


        self.bottomLbl = wx.StaticBox(self,-1,'bottom:')
        self.bottomSizer = wx.StaticBoxSizer(self.bottomLbl, wx.VERTICAL)
        self.hbox3=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.ComboBox3, 0, border=3, flag=flags)
        self.hbox3.Add(self.Pbutton3, 0, border=3, flag=flags)
        self.hbox3.Add(self.Nbutton3, 0, border=3, flag=flags)
        self.bottomSizer.Add(self.hbox3)


        print('2')

        self.vboxCombo.Add(self.topSizer)
        self.vboxCombo.Add(self.middleSizer)
        self.vboxCombo.Add(self.bottomSizer)

        self.hbox.Add(self.vboxCombo)
        self.hbox.Add(self.vboxCntr)




        #self.vboxS=wx.BoxSizer(wx.VERTICAL)
        #self.vboxS.Add(self.selectbutton, 0, border=3, flag=flags)
        #self.vboxS.Add(self.deselectbutton, 0, border=3, flag=flags)
        #self.vboxS.Add(self.deletebutton, 0, border=3, flag=flags)
        #self.vboxS.Add(self.addbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.vboxS)


        #self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        #self.hbox.Add(self.width2Lab,0,border=3,flag=flags)
        #self.hbox.Add(self.width2Box,0,border=3,flag=flags)


        self.vboxB=wx.BoxSizer(wx.VERTICAL)


        self.hbox.Add(self.vboxB)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.SetSizerAndFit(self.vbox)


        print('Reading projections...')
        """
        self.Xs_xy,self.Ys_xy,self.Zs_xy=self.parent.tabTwo.Get2D(self.tabOne.labb[1],self.tabOne.labb[0],transpose='n')
        self.Xs_xy_xmin=numpy.min(self.Xs_xy)
        self.Xs_xy_xmax=numpy.max(self.Xs_xy)
        self.Xs_xy_ymin=numpy.min(self.Ys_xy)
        self.Xs_xy_ymax=numpy.max(self.Ys_xy)

        self.Xs_yz,self.Ys_yz,self.Zs_yz=self.parent.tabTwo.Get2D(self.tabOne.labb[2],self.tabOne.labb[1],transpose='n')
        self.Xs_yz_xmin=numpy.min(self.Xs_yz)
        self.Xs_yz_xmax=numpy.max(self.Xs_yz)
        self.Xs_yz_ymin=numpy.min(self.Ys_yz)
        self.Xs_yz_ymax=numpy.max(self.Ys_yz)
        """

        print('Done')


    #make an index
    def index(self,array):
        index=[]
        for i in range(len(array)):
            index.append((array[i].name))
        return index

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


    #get 2d strips from 3d data
    def GetSlice2d(self):
        pass
        #dic=self..dic
        #self.data=self.parent.spec['hnco'].data
        #self.index0=self.parent.spec['hnco'].index0
        #self.index1=self.parent.spec['hnco'].index1
        #self.index2=self.parent.spec['hnco'].index2

        """
        self.DECON=0
        if(self.tabOne.DECON==1):
            self.datadec=self.tabOne.datadec
            if(self.datadec.shape==self.data.shape):
                print('Shape of deconvolved, and raw are different')
                self.DECON=0
            else:
                self.DECON=1
        """

    #get 2d strips from 3d data
    def ReSlice2d(self,arr,inc,pkl,peak,width,orth=0,lab='hnco'):

        if(orth==0):
            #print(out 2D slice for each peak correlation)
            #print("Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width)
            ptC=self.parent.spec[lab].pkIdx[pkl][0]
            ptC=ptC+inc
            ptH=self.parent.spec[lab].pkIdx[pkl][1]

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


    def AddSlice(self,sele,ax,spec,xsub,ysub,cmap):
        #ADD THE SLICE
        if(self.cblist[spec].IsChecked()==False):
            return 0

        levels=self.GetLevels(float(self.textbox1[spec].GetValue()),float(self.textbox2[spec].GetValue()),int(self.textbox3[spec].GetValue()))

        if(ax==0):
            axes=self.axes0
        if(ax==1):
            axes=self.axes1
        if(ax==2):
            axes=self.axes2
        if(ax==3):
            axes=self.axes3
        if(ax==4):
            axes=self.axes4
        if(ax==5):
            axes=self.axes5
        if(ax==6):
            axes=self.axes6
        if(ax==7):
            axes=self.axes7
        if(ax==8):
            axes=self.axes8


        if(self.orth_cb and self.dec_cb):
            Xs,Ys,Zs=self.ReSlice2d(self.datadec,self.inc,sele,self.parent.spec[spec].peak2D,self.Width2,orth=self.orth_cb,lab='hnco')
        elif(self.orth_cb): #orthoganol
            Xs,Ys,Zs=self.ReSlice2d(self.parent.spec[spec].data,self.inc,sele,self.parent.spec[spec].peak2D,self.Width2,orth=self.orth_cb,lab=spec)
        elif(self.dec_cb): #decon
            Xs,Ys,Zs=self.ReSlice2d(self.datadec,self.inc,sele,self.parent.spec[spec].peak2D,self.Width1,orth=self.orth_cb,lab='hnco')
        else: #other combo box
            Xs,Ys,Zs=self.ReSlice2d(self.parent.spec[spec].data,self.inc2,sele,self.parent.spec[spec].peak2D,self.Width1,orth=self.orth_cb,lab=spec)

        yave=numpy.average(Ys)
        if(ysub!=0):
            ysa=yave
        else:
            ysa=0


        clbox=self.comboClbox[spec].GetSelection()

        cmap=cm.get_cmap(self.cmaps[clbox])

        axes.contour(Xs-xsub, Ys-ysa+ysub, Zs,levels,cmap=cmap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network


        """
        if(self.ax_reset[spec]==1):
            y_max2=Ys[0][0]
            y_min2=Ys[(len(Ys))-1][0]
            x_max2=Xs[0][0]
            x_min2=Xs[0][(len(Xs[0]))-1]
            axes.set_xlim(x_min2,x_max2)
            axes.set_ylim(y_min2,y_max2)
            self.ax_reset[spec]=0
        else:
            axes.set_xlim(self.x_min[spec],self.x_max[spec])
            axes.set_ylim(self.y_min[spec],self.y_max[spec])
        """

        y_max2a=Ys[0][0]
        y_min2a=Ys[(len(Ys))-1][0]
        x_max2a=Xs[0][0]
        x_min2a=Xs[0][(len(Xs[0]))-1]

        if(ysub==0):
            yline=yave
        else:
            yline=ysub

        if(self.grid_cb):#horizontal line
            xl=(x_max2a,x_min2a)
            if(self.cb_flip.GetValue()): #if we want orthogonal, get value from first tick box
                hl=(self.dimC,self.dimC)
            elif(self.cb_decon.GetValue()): #otherwise, if want to see decon
                hl=(self.dimH,self.dimH)
            else: #else go into the other combo-box.
                hl=(yline,yline)
            axes.plot(xl,hl,'blue',zorder=0)



            """
            #plot NOE line
            if(self.orth_cb): #if orthoganol
                xl=(self.dimC,self.dimC2)
                hl=(self.dimC,self.dimC)
                axes.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,self.dimC)
                xd=(self.dimC,self.dimC)
                axes.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,self.dimC)
                xd=(self.dimC2,self.dimC2)
                axes.plot(xd,yd,'cyan') #vertical 2
            elif(self.dec_cb):
                xl=(self.dimC,self.dimC2)
                hl=(self.dimH,self.dimH)
                axes.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,self.dimH)
                xd=(self.dimC2,self.dimC2)
                axes.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,self.dimH2)
                xd=(self.dimC,self.dimC)
                axes.plot(xd,yd,'cyan') #vertical 2
            else:
                xl=(self.dimC,self.dimC2)
                hl=(self.dimH2,self.dimH2)
                axes.plot(xl,hl,'cyan') #horizontal
                yd=(y_max2a,self.dimH2)
                xd=(self.dimC2,self.dimC2)
                axes.plot(xd,yd,'cyan') #vertical 1
                yd=(y_max2a,self.dimH2)
                xd=(self.dimC,self.dimC)
                axes.plot(xd,yd,'cyan') #vertical 2
            """


        if(self.cblab[spec].IsChecked()==False):
            return yave


        #do the cross peak labels
        if(self.grid_cb):
            #for cn in self.parent.spec[spec].conn_data:
            #print('shit')
            #print(self.parent.peak[self.index_data[sele]][spec])
            #print(self.parent.peak[self.index_data[sele]].keys())

            for cn in self.parent.peak[self.index_data[sele]][spec]:

                    #axes.scatter(cn.f3-xsub,cn.f1-ysa+ysub,c='k',s=50,zorder=2,marker='x')
                    #axes.text(cn.f3-xsub,cn.f1-ysa+ysub,cn.name,rotation=90,fontsize=8)

                    axes.scatter(cn.f3,cn.f1-ysa+ysub,c='k',s=50,zorder=2,marker='x')
                    axes.text(cn.f3,cn.f1-ysa+ysub,cn.name,rotation=90,fontsize=8)

                    #if(spec=='hncaco'):
                    #xl=(cn.f3-xsub,cn.f3-xsub)
                    xl=(cn.f3,cn.f3)

                    #y_max2a=Ys[0][0]
                    #y_min2a=Ys[(len(Ys))-1][0]



                    if(cn.tp=='main' or cn.tp=='main2'):
                        cl='r'
                        cl2='red'
                        di='u'
                    else:
                        cl='b'
                        cl2='blue'
                        di='d'
                    #axes.arrow(xl[0],yline,10.,self.Width1*10., shape='full', lw=0, length_includes_head=True,zorder=2)
                    #if(spec=='hcaco'):
                    hl=(yline-self.Width1/2.,yline+self.Width1/2.)

                    if(di=='u'):
                        #up arrow

                        axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                        axes.arrow(xl[0],yline+self.Width1/2.-self.Width1/9.,0.,self.Width1/10., shape='full', lw=0, length_includes_head=True,zorder=2,color=cl,head_width=0.5,head_length=self.Width1/10.)
                    else:
                        #down arrow
                        axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                        axes.arrow(xl[0],yline-self.Width1/2.+self.Width1/9.,0.,-self.Width1/10., shape='full', lw=0, length_includes_head=True,zorder=2,color=cl,head_width=0.5,head_length=self.Width1/10.)



            """
            for sele in self.selection:
                cn=self.parent.spec['hnco'].conn_data[sele]
                if(orth_cb and dec_cb):
                    if(cn.v1==sele1):
                        axes.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        axes.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8,color='r')
                if(orth_cb): #if orthoganol
                    if(cn.v1==sele1):
                        axes.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        axes.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8,color='r')
                elif(dec_cb):
                    if(cn.v1==sele1):
                        axes.scatter(cn.f3,cn.f1,c='k',s=50,zorder=2,marker='x')
                        axes.text(cn.f3,cn.f1,cn.tag2,rotation=90,fontsize=8,color='r')
                else: #if not orthoganol
                    if(cn.v1==sele2):
                        axes.scatter(cn.f3,cn.f1,c='k',s=50,zorder=2,marker='x')
                        axes.text(cn.f3,cn.f1,cn.tag2,rotation=90,fontsize=8,color='r')
            """

        """
        #do the main peak labels
        if(self.orth_cb and self.dec_cb):
            axes.text(self.dimC,self.dimC,self.index_data[self.sele2],rotation=90,fontsize=8)
            axes.scatter(self.dimC,self.dimC,c='g',s=100)
        if(self.orth_cb):#if orthoganol
            axes.text(self.dimC,self.dimC,self.index_data[self.sele2],rotation=90,fontsize=8)
            axes.scatter(self.dimC,self.dimC,c='g',s=100)
        elif(self.dec_cb):
            axes.text(self.dimH,self.dimH,self.index_data[self.sele2],rotation=90,fontsize=8)
            axes.scatter(self.dimH,self.dimH,c='g',s=100)
        else:
            axes.text(self.dimC2,self.dimH2,self.index_data[self.sele2],rotation=90,fontsize=8)
            axes.scatter(self.dimC2,self.dimH2,c='g',s=100)
        """

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



        #colormap=cm.seismic

        #colormap=cm.RdYlBu
        from matplotlib.gridspec import GridSpec
        gs1=GridSpec(3,9)

        self.fig.clear()

        self.Width1=float(self.width1Box.GetValue())
        #self.Width2=float(self.width2Box.GetValue())
        if(self.Width1==0):
            self.Width1=1
            self.width1Box.SetValue(str(1))
        #if(self.Width2==0):
        #    self.Width2=1
        #    self.width2Box.SetValue(str(1))





        self.sele1=self.ComboBox1.GetSelection()
        self.sele2=self.ComboBox2.GetSelection()
        self.sele3=self.ComboBox3.GetSelection()

        self.dimC=float(self.parent.spec['hnco'].peak2D[self.sele1].ppmJ)
        self.dimH=float(self.parent.spec['hnco'].peak2D[self.sele1].ppmK)

        label=self.parent.spec['hnco'].peak2D[self.ComboBox1.GetSelection()].name

        self.dimC2=float(self.parent.spec['hnco'].peak2D[self.sele2].ppmJ)
        self.dimH2=float(self.parent.spec['hnco'].peak2D[self.sele2].ppmK)

        self.orth_cb=self.cb_flip.GetValue()
        self.dec_cb=self.cb_decon.GetValue()
        self.grid_cb=self.cb_grid_auto.GetValue()

        self.axes2 = self.fig.add_subplot(gs1[1,:3])
        self.axes1 = self.fig.add_subplot(gs1[1,3:6])
        self.axes0 = self.fig.add_subplot(gs1[1,6:9])

        self.axes5 = self.fig.add_subplot(gs1[0,:3])
        self.axes4 = self.fig.add_subplot(gs1[0,3:6])
        self.axes3 = self.fig.add_subplot(gs1[0,6:9])

        self.axes8 = self.fig.add_subplot(gs1[2,:3])
        self.axes7 = self.fig.add_subplot(gs1[2,3:6])
        self.axes6 = self.fig.add_subplot(gs1[2,6:9])



        if(self.orth_cb):
            self.axes2.set_xlabel(self.parent.spec['hnco'].labb[0]+'(ppm)',fontsize=8)
            self.axes2.set_ylabel(self.parent.spec['hnco'].labb[1]+'(ppm)',fontsize=8)
        else:
            self.axes2.set_xlabel(self.parent.spec['hnco'].labb[0]+'(ppm)',fontsize=8)
            self.axes2.set_ylabel(self.parent.spec['hnco'].labb[2]+'(ppm)',fontsize=8)

        if(self.orth_cb and self.dec_cb):
            self.axes1.set_xlabel(self.parent.spec['hnca'].labb[0]+'(ppm)',fontsize=8)
            #self.axes1.set_ylabel(self.parent.spec['hnca'].labb[1]+'(ppm)',fontsize=8)
        else:
            self.axes1.set_xlabel(self.parent.spec['hnca'].labb[0]+'(ppm)',fontsize=8)
            #self.axes1.set_ylabel(self.parent.spec['hnca'].labb[2]+'(ppm)',fontsize=8)

        if(self.orth_cb and self.dec_cb):
            self.axes0.set_xlabel(self.parent.spec['hncacb'].labb[0]+'(ppm)',fontsize=8)
            #self.axes0.set_ylabel(self.parent.spec['hncacb'].labb[1]+'(ppm)',fontsize=8)
        else:
            self.axes0.set_xlabel(self.parent.spec['hncacb'].labb[0]+'(ppm)',fontsize=8)
            #self.axes0.set_ylabel(self.parent.spec['hncacb'].labb[2]+'(ppm)',fontsize=8)



        cmap1=cm.bwr
        cmap2=cm.RdYlBu



        ysa=self.AddSlice(self.sele2,2,'hnco',0,0,cmap1)
        ysa2=self.AddSlice(self.sele2,2,'hncaco',self.parent.HNCOmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele2,1,'hnca',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele2,1,'hncoca',self.parent.HNCAmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele2,0,'hncacb',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele2,0,'hncocacb',0,ysa,cmap2)


        ysa=self.AddSlice(self.sele1,5,'hnco',0,0,cmap1)
        ysa2=self.AddSlice(self.sele1,5,'hncaco',self.parent.HNCOmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele1,4,'hnca',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele1,4,'hncoca',self.parent.HNCAmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele1,3,'hncacb',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele1,3,'hncocacb',0,ysa,cmap2)

        ysa=self.AddSlice(self.sele3,8,'hnco',0,0,cmap1)
        ysa2=self.AddSlice(self.sele3,8,'hncaco',self.parent.HNCOmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele3,7,'hnca',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele3,7,'hncoca',self.parent.HNCAmed,ysa,cmap2)

        ysa2=self.AddSlice(self.sele3,6,'hncacb',0,ysa,cmap1)
        ysa2=self.AddSlice(self.sele3,6,'hncocacb',0,ysa,cmap2)


        self.canvas.draw()
        return
        ##############################################################333
        #Subplot 3 - the slice



        if(orth_cb and dec_cb):
            self.axes1.set_xlabel(self.parent.spec['hnca'].labb[0]+'(ppm)',fontsize=8)
            self.axes1.set_ylabel(self.parent.spec['hnca'].labb[1]+'(ppm)',fontsize=8)
            Xs,Ys,Zs=self.ReSlice2d(self.data,self.inc,self.sele2,self.parent.spec['hnca'].peak2D,Width2,orth=orth_cb,lab='hnca')
        else:
            self.axes1.set_xlabel(self.parent.spec['hnca'].labb[0]+'(ppm)',fontsize=8)
            self.axes1.set_ylabel(self.parent.spec['hnca'].labb[2]+'(ppm)',fontsize=8)

            Xs,Ys,Zs=self.ReSlice2d(self.parent.spec['hnca'].data,self.inc,self.sele2,self.parent.spec['hnca'].peak2D,Width1,lab='hnca')


        ysa2=numpy.average(Ys)
        print(ysa2)
        self.axes1.contour(Xs, Ys-ysa2+ysa, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network


        Xt,Yt,Zt=self.ReSlice2d(self.parent.spec['hncoca'].data,self.inc,sele2,self.parent.spec['hncoca'].peak2D,Width1,lab='hncoca')


        ytb2=numpy.average(Yt)

        self.axes1.contour(Xt-self.parent.HNCAmed, Yt-ytb2+ysa, Zt,levels,cmap=cm.RdYlBu,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network



        #do the main peak labels
        if(orth_cb and dec_cb):
            #print(dimC,dimC2,dimH,dimH2)
            self.axes1.text(dimC,dimC2,self.index_data[sele2],rotation=90,fontsize=8)
            self.axes1.scatter(dimC,dimC2,c='g',s=100)
        else:
            self.axes1.text(dimC2,dimH2,self.index_data[sele2],rotation=90,fontsize=8)
            self.axes1.scatter(dimC2,dimH2,c='g',s=100)


        y_max2a=Ys[0][0]
        y_min2a=Ys[(len(Ys))-1][0]
        x_max2a=Xs[0][0]
        x_min2a=Xs[0][(len(Xs[0]))-1]

        #print(x_max2a,x_min2a,y_max2a,y_min2a)
        if(grid_cb):#if NOESY list is present
            xl=(x_min2a,x_max2a)
            if(orth_cb and dec_cb):
                hl=(dimC2,dimC2)
            else:
                hl=(dimH2,dimH2)
            self.axes1.plot(xl,hl,'green') #horizontal


            if(orth_cb and dec_cb):
                xl=(dimC2,dimC)
                hl=(dimC2,dimC2)
                self.axes1.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimC)
                xd=(dimC,dimC)
                self.axes1.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimC)
                xd=(dimC2,dimC2)
                self.axes1.plot(xd,yd,'cyan') #vertical 2
            else:
                xl=(dimC2,dimC2)
                hl=(dimH2,dimH2)
                self.axes1.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimH)
                xd=(dimC2,dimC2)
                self.axes1.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimH)
                xd=(dimC2,dimC2)
                self.axes1.plot(xd,yd,'cyan') #vertical 2



        #do the cross peak labels
        if(grid_cb):
            if(len(self.selection)==0):
                sel=True
            else:
                sel=False
            for cn in self.parent.spec['hnca'].conn_data:
                if(orth_cb and dec_cb): #if orthoganol
                    if(cn.v1==sele2):
                        self.axes1.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        if(sel):
                            self.axes1.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8)
                else: #if not orthoganol
                    if(cn.v1==sele2):
                        self.axes1.scatter(cn.f3,cn.f1-ysa2+ysa,c='k',s=50,zorder=2,marker='x')
                        if(sel):
                            self.axes1.text(cn.f3,cn.f1-ysa2+ysa,cn.tag2,rotation=90,fontsize=8)
            for sele in self.selection:
                cn=self.parent.spec['hnca'].conn_data[sele2]
                if(orth_cb and dec_cb): #if orthoganol
                    if(cn.v1==sele2):
                        self.axes1.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        self.axes1.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8,color='r')
                else: #if not orthoganol
                    if(cn.v1==sele2):
                        self.axes1.scatter(cn.f3,cn.f1-ysa2+ysa,c='k',s=50,zorder=2,marker='x')
                        self.axes1.text(cn.f3,cn.f1-ysa2+ysa,cn.tag2,rotation=90,fontsize=8,color='r')



        if(self.ax_reset1==1):
            y_max1=Ys[0][0]-ysa2+ysa
            y_min1=Ys[(len(Ys))-1][0]-ysa2+ysa
            x_max1=Xs[0][0]
            x_min1=Xs[0][(len(Xs[0]))-1]
            self.axes1.set_xlim(x_min1,x_max1)
            self.axes1.set_ylim(y_min1,y_max1)
            self.ax_reset1=0
        else:
            self.axes1.set_xlim(x_min1,x_max1)
            self.axes1.set_ylim(y_min1,y_max1)



        ##############################################################333
        #Subplot 3 - the slice 333



        ysa2=numpy.average(Ys)
        print(ysa2)
        self.axes0.contour(Xs, Ys-ysa2+ysa, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network



        Xt,Yt,Zt=self.ReSlice2d(self.parent.spec['hncocacb'].data,self.inc,sele2,self.parent.spec['hncocacb'].peak2D,Width1,lab='hncocacb')
        ytb2=numpy.average(Yt)
        self.axes0.contour(Xt-self.parent.HNCAmed, Yt-ytb2+ysa, Zt,levels,cmap=cm.RdYlBu,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network



        #do the main peak labels
        if(orth_cb and dec_cb):
            #print(dimC,dimC2,dimH,dimH2)
            self.axes0.text(dimC,dimC2,self.index_data[sele2],rotation=90,fontsize=8)
            self.axes0.scatter(dimC,dimC2,c='g',s=100)
        else:
            self.axes0.text(dimC2,dimH2,self.index_data[sele2],rotation=90,fontsize=8)
            self.axes0.scatter(dimC2,dimH2,c='g',s=100)


        y_max2a=Ys[0][0]
        y_min2a=Ys[(len(Ys))-1][0]
        x_max2a=Xs[0][0]
        x_min2a=Xs[0][(len(Xs[0]))-1]

        #print(x_max2a,x_min2a,y_max2a,y_min2a)
        if(grid_cb):#if NOESY list is present
            xl=(x_min2a,x_max2a)
            if(orth_cb and dec_cb):
                hl=(dimC2,dimC2)
            else:
                hl=(dimH2,dimH2)
            self.axes0.plot(xl,hl,'green') #horizontal


            if(orth_cb and dec_cb):
                xl=(dimC2,dimC)
                hl=(dimC2,dimC2)
                self.axes0.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimC)
                xd=(dimC,dimC)
                self.axes0.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimC)
                xd=(dimC2,dimC2)
                self.axes0.plot(xd,yd,'cyan') #vertical 2
            else:
                xl=(dimC2,dimC2)
                hl=(dimH2,dimH2)
                self.axes0.plot(xl,hl,'cyan') #horizontal
                yd=(y_min2a,dimH)
                xd=(dimC2,dimC2)
                self.axes0.plot(xd,yd,'cyan') #vertical 1
                yd=(y_min2a,dimH)
                xd=(dimC2,dimC2)
                self.axes0.plot(xd,yd,'cyan') #vertical 2



        #do the cross peak labels
        if(grid_cb):
            if(len(self.selection)==0):
                sel=True
            else:
                sel=False
            for cn in self.parent.spec['hncacb'].conn_data:
                if(orth_cb and dec_cb): #if orthoganol
                    if(cn.v1==sele2):
                        self.axes0.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        if(sel):
                            self.axes0.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8)
                else: #if not orthoganol
                    if(cn.v1==sele2):
                        self.axes0.scatter(cn.f3,cn.f1-ysa2+ysa,c='k',s=50,zorder=2,marker='x')
                        if(sel):
                            self.axes0.text(cn.f3,cn.f1-ysa2+ysa,cn.tag2,rotation=90,fontsize=8)
            for sele in self.selection:
                cn=self.parent.spec['hncacb'].conn_data[sele2]
                if(orth_cb and dec_cb): #if orthoganol
                    if(cn.v1==sele2):
                        self.axes0.scatter(cn.f3,cn.f2,c='k',s=50,zorder=2,marker='x')
                        self.axes0.text(cn.f3,cn.f2,cn.tag2,rotation=90,fontsize=8,color='r')
                else: #if not orthoganol
                    if(cn.v1==sele2):
                        self.axes0.scatter(cn.f3,cn.f1-ysa2+ysa,c='k',s=50,zorder=2,marker='x')
                        self.axes0.text(cn.f3,cn.f1-ysa2+ysa,cn.tag2,rotation=90,fontsize=8,color='r')



        if(self.ax_reset0==1):
            y_max1=Ys[0][0]-ysa2+ysa
            y_min1=Ys[(len(Ys))-1][0]-ysa2+ysa
            x_max1=Xs[0][0]
            x_min1=Xs[0][(len(Xs[0]))-1]

            print('hncacb: ',numpy.max(Xs),numpy.min(Xs),numpy.max(Ys),numpy.min(Ys))
            self.axes0.set_xlim(x_min1,x_max1)
            self.axes0.set_ylim(y_min1,y_max1)
            self.ax_reset0=0
        else:
            self.axes0.set_xlim(x_min1,x_max1)
            self.axes0.set_ylim(y_min1,y_max1)





        self.canvas.draw()

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
        self.inc=0
        self.inc2=0
        self.ax_reset0=1
        self.ax_reset1=1
        self.ax_reset2=1
        self.ax_resetCC=1
        self.ax_resetHC=1
        self.draw_figure()

    def on_N_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        if(self.cb_flip.GetValue()):
            self.ax_reset2=1
        if(self.cb_decon.GetValue()):
            self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        if(self.cb_flip.GetValue()):
            self.ax_reset2=1
        if(self.cb_decon.GetValue()):
            self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.selection=[]
        self.draw_figure()


    def on_bigN_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        if(self.cb_flip.GetValue()):
            self.ax_reset2=1
        if(self.cb_decon.GetValue()):
            self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.ComboBox3.SetSelection(self.ComboBox3.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_bigP_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        if(self.cb_flip.GetValue()):
            self.ax_reset2=1
        if(self.cb_decon.GetValue()):
            self.ax_reset2=1

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
        self.ax_reset1=1
        self.ax_reset0=1
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button2(self, event):
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
        print(self.selection)
        self.selection=sorted(self.selection,reverse=True)
        print(self.selection)

        for sele in self.selection:
            print('removing:',self.parent.spec['hnco'].conn_data[sele].tag)
            self.parent.spec['hnco'].conn_data.pop(sele)

        self.selection=[]
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
        print
        print('Saving ',len(self.parent.spec['hnco'].conn_data),'entries in connectivity table to ',self.outfile)
        outy=open(self.outfile,'w')
        for i in range(len(self.parent.spec['hnco'].conn_data)):
            outy.write('%s\t%s\n' % (self.index_data[self.parent.spec['hnco'].conn_data[i][0]],self.index_data[self.parent.spec['hnco'].conn_data[i][1]]))
        outy.close()


    def on_launch(self,event):
        print(event)
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

        print(self.parent.spec['hnco'].noeTags)
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

        print('faaaa')
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
        #print(event.xdata,event.ydata)
        if(self.SELECT==1):
            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            argy=[]
            cbs1=self.ComboBox1.GetSelection()
            print(cbs1)
            for i,cn in enumerate(self.parent.spec['hnco'].conn_data):
                #print(self.parent.spec['hnco'].peak2D[p].name)
                if(cn.v1==cbs1):
                    yval=cn.f3  #carbon
                    #print(cn.v1,cn.v2,event.xdata,event.ydata,cn.f1,cn.f2,cn.f3,cn.f4)
                    rad2=((yval-event.xdata)/ydist)**2.
                    raddy.append(rad2)
                    argy.append(i)
            raddy=numpy.array(raddy)
            maxy=argy[numpy.argmin(raddy)]
            self.selection.append(maxy)

            cn=self.parent.spec['hnco'].conn_data[maxy]
            print('selected:',cn.tag)

            self.sym=self.parent.spec['hnco'].cb_grid.IsChecked()
            if(self.sym==True): #do second frame if needed.
                for i,pk in enumerate(self.parent.spec['hnco'].peak2D):
                    if(cn.p2==pk.name):
                        self.ComboBox2.SetSelection(i)
                for i,cn2 in enumerate(self.parent.spec['hnco'].conn_data):
                    if(cn2.p1==cn.p2):
                        if(cn2.p2==cn.p1):
                            print('Reciprocated:',cn2.tag)
                            self.selection.append(i)

            self.SELECT=0
            print(self.selection,self.SELECT)
            self.draw_figure()
        if(self.SEARCH==1):


            self.selection=[]

            self.sym=self.parent.spec['hnco'].cb_grid.IsChecked()
            if(self.sym==False):
                print('No protocol for adding peaks yet')
                pk=self.parent.spec['hnco'].peak2D[self.ComboBox1.GetSelection()]
                f1=pk.ppmK
                f2=pk.ppmJ
                f3=event.xdata

                from deconFrame import findnear_index
                i=findnear_index(self.parent.spec['hnco'].index0,event.xdata)

                f3=self.parent.spec['hnco'].index0[i]
                inty=self.parent.spec['hnco'].data[pk.indexI,pk.indexJ,i]
                print(f1,f2,f3,inty)
                print('Adding NOE:')
                from deconFrame import connEntry

                cnt=0
                for cn in self.parent.spec['hnco'].conn_data:
                    ref=cn.p1.split('_')[0]
                    if(ref==pk.name):
                        cnt+=1

                stry=('%s_%i\t%f\t%f\t%f\t%e\t') % (pk.name,cnt+1,f1,f2,f3,inty)
                print(stry)
                test=stry.split()
                cnNew=connEntry(test,sym='n',peak=self.parent.spec['hnco'].peak2D,dim=3)
                print('added')
                self.parent.spec['hnco'].conn_data.append(cnNew)
                self.selection.append(len(self.parent.spec['hnco'].conn_data)-1)
                self.draw_figure()
                self.SEARCH=0
                return

            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            argy=[]
            cbs1=self.ComboBox1.GetSelection()
            print(cbs1)
            for i,pk in enumerate(self.parent.spec['hnco'].peak2D):
                xval=pk.ppmI  #proton
                #yval=pk.ppmJ  #carbon
                #print(cn.v1,cn.v2,event.xdata,event.ydata,cn.f1,cn.f2,cn.f3,cn.f4)
                rad2=((xval-event.xdata)/xdist)**2
                raddy.append(rad2)
                argy.append(i)
            raddy=numpy.array(raddy)
            maxy=argy[numpy.argmin(raddy)]



            print('Testing for NOE:')
            add=0

            pk1=self.parent.spec['hnco'].peak2D[self.ComboBox1.GetSelection()]
            pk2=self.parent.spec['hnco'].peak2D[maxy]
            if(pk1.name==pk2.name):
                print('this is a diagonal peak. Not adding')
                add=1

            print(pk1.name,pk2.name)

            #set second combobox to identity of cross peak
            for i,pk in enumerate(self.parent.spec['hnco'].peak2D):
                if(pk2.name==pk.name):
                    self.ComboBox2.SetSelection(i)

            #is this already here?
            for i,cn in enumerate(self.parent.spec['hnco'].conn_data):
                if(cn.p1==pk1.name and cn.p2==pk2.name):
                    print('NOE already present:',cn.tag)
                    add=1
                if(cn.p2==pk1.name and cn.p1==pk2.name):
                    print('NOE (reciprocated) already present:',cn.tag)
                    add=1
            i1=self.parent.spec['hnco'].data[pk1.indexI,pk1.indexJ,pk2.indexK]
            i2=self.parent.spec['hnco'].data[pk2.indexI,pk2.indexJ,pk1.indexK]

            noise=self.parent.spec['hnco'].dmax*float(self.parent.spec['hnco'].threshBox.GetValue())
            print('noise:',noise)
            print('Intensity of selected cross peak:')
            print(i1,i1/noise)
            print('Intensity of reciprocated cross peak:')
            print(i2,i2/noise)
            if(numpy.fabs(i1)<noise):
                print('I1 less than noise. Probably not a cross peak.')
                add=1
            if(numpy.fabs(i2)<noise):
                print('I2 less than noise. Probably not a cross peak.')
                add=1
            print('Ratio (expected to be >50%):')
            print(min(i1,i2)/max(i1,i2)*100.)


            if(add==0):
                print('Adding NOE:')
                from deconFrame import connEntry


                stry=('%s\t%s\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f') % (pk1.name,pk2.name,pk1.ppmK,pk1.ppmJ,pk2.ppmI,i1,i2,i1/noise,i2/noise,0,0,0,0)
                test=stry.split()
                print(stry)
                cnNew=connEntry(test,sym='y',peak=self.parent.spec['hnco'].peak2D,dim=3)
                print('added')
                self.parent.spec['hnco'].conn_data.append(cnNew)
                self.selection.append(len(self.parent.spec['hnco'].conn_data)-1)

                stry=('%s\t%s\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f') % (pk2.name,pk1.name,pk2.ppmK,pk2.ppmJ,pk1.ppmI,i2,i1,i2/noise,i1/noise,0,0,0,0)
                print(stry)
                test=stry.split()
                cnNew=connEntry(test,sym='y',peak=self.parent.spec['hnco'].peak2D,dim=3)
                print('added')
                self.parent.spec['hnco'].conn_data.append(cnNew)
                self.selection.append(len(self.parent.spec['hnco'].conn_data)-1)
                print(self.selection)
                print(len(self.parent.spec['hnco'].conn_data))
            else:
                print('Not adding.')

            #self.selection.append()
            self.draw_figure()
            self.SEARCH=0



    """
    #when search button is pressed make selection
    def on_pick(self, event):
        if(self.cb_grid_select.GetValue()==1):
            if(self.pick_cnt==0):
                self.selection=[]
                self.select=[]
                self.select.append((event.xdata,event.ydata))
                self.pick_cnt=1
                print('Click 1: ',event.xdata,event.ydata)

            else:
                self.selection=[]
                self.select.append((event.xdata,event.ydata))
                self.pick_cnt=0
                print('Click 2: ',event.xdata,event.ydata)
                for i in range(len(self.parent.spec['hnco'].peak2D)):
                    if(float(self.parent.spec['hnco'].peak2D[i][1])>float(self.select[0][0]) and float(self.parent.spec['hnco'].peak2D[i][1])<float(self.select[1][0])):
                        self.selection.append(i)
                print('Peaks in this selection range:',ItoN1(self.selection,self.index_data))


                self.xmin,self.xmax=self.axes.get_xlim()
                self.ymin,self.ymax=self.axes.get_ylim()

                self.draw_figure()
                self.cb_grid_select.SetValue(0)

    """


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
