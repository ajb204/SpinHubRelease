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

class singleFrame(wx.Panel):

    def __init__(self,parent):
        #this is the walk frame

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
        #self.sym=self.tabOne.cb_grid.IsChecked(
        #copy in the previous variables
        self.index_dataRef=self.index(self.parent.spec[self.parent.refSpec].peak2D) #inherit the data index
        # print(self.index_data)f

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


        #self.data=self.parent.spec['hnco'].data
        self.create_main_panel()
        # print('22')
        self.draw_figure()
        # print('44')
        #self.Show(True)
        self.Fit()


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.

        
        self.fig=plt.Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        
        #colormap=cm.RdYlBu
        #self.fig.clear()
        self.SetSlices() #setup indices for the frames
        self.wipe() #set selection=0 for all axes
        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        self.cmaps =  'seismic','bwr','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'
        print('HHHHHHHH')
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
        self.refbox={}
        self.ref_text={}

        flags = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.thresh=1E5
        self.vboxCntr_label = wx.StaticBox(self,-1,'Spectra:')
        self.vboxCntr = wx.StaticBoxSizer(self.vboxCntr_label, wx.VERTICAL)

        #self.parent.spec=list(self.parent.spec.keys())
        for i,spec in enumerate(self.parent.spec.keys()):
            print (i,spec)
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
            self.cntrSizer[spec].Add(wx.StaticText(self,-1,spec.upper(),size=(77,-1)), 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.cblist[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.cblab[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox1[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox2[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.text3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.textbox3[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.comboClbox[spec], 0, border=3, flag=flags)
            self.cntrSizer[spec].Add(self.launch[spec], 0, border=3, flag=flags)


            #if(spec=='hnco' or spec=='hnca' or spec=='hncacb'):
                # self.ref_text[spec]=wx.StaticText(self,-1, '')
                #if(spec=='hnco'):
                #    self.ref_text[spec]=wx.StaticText(self,-1, 'CO Offset:')
                #if(spec=='hnca'):
                #    self.ref_text[spec]=wx.StaticText(self,-1, 'CA Offset:')
                #if(spec=='hncacb'):
                #    self.ref_text[spec]=wx.StaticText(self,-1, 'CACB Offset:')
                #self.cntrSizer[spec].Add(self.ref_text[spec], 0, border=3, flag=flags)
            #    pass
            specList='hncaco','hncoca','hncocacb','hncacb','ctocsy','hcconh'
            if(spec in specList):
                self.refbox[spec]=wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)
                if(spec=='hncaco'):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'CACO Offset:')
                    self.refbox[spec].SetValue(str(self.parent.HNCOmed))
                elif(spec=='hncoca'):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'COCA Offset:')
                    self.refbox[spec].SetValue(str(self.parent.HNCAmed))
                elif(spec=='hncocacb'):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'COCACB Offset:')
                    self.refbox[spec].SetValue(str(self.parent.HNCACBmed))
                elif(spec=='hncacb' ):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'CACB Offset:')
                    self.refbox[spec].SetValue(str(self.parent.HNCAHNCACBmed))
                elif(spec=='ctocsy' ):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'CTOCSY Offset:')
                    self.refbox[spec].SetValue(str(self.parent.CTOCSYmed))
                elif(spec=='hcconh' ):
                    self.ref_text[spec]=wx.StaticText(self,-1, 'HTOCSY Offset:')
                    self.refbox[spec].SetValue(str(self.parent.HTOCSYmed))

                self.cntrSizer[spec].Add(self.ref_text[spec], 0, border=3, flag=flags)
                self.cntrSizer[spec].Add(self.refbox[spec], 0, border=3, flag=flags)
                
            self.vboxCntr.Add(self.cntrSizer[spec])

        print('aya') 

        #print(self.parent.spec['hnco'].peak2D #MASTER PEAK LIST)
        #NEED TO ADD SOME ERROR CHECKING: LOOK AT ALL THE PEAK LISTS AND MAKE SURE THEY ARE ALL
        #THE SAME
        #self.listy=[]
        #for i in range(len(self.parent.spec['hnco'].peak2D)):
        #    self.listy.append(self.parent.spec['hnco'].peak2D[i].name)


       
        self.listy,skip=self.parent.OrderPeaks(self.parent.G1edges,skip=False)
       
        self.index_data=copy.deepcopy(self.listy)
        #self.index_data=self.index(self.parent.spec['hnco'].peak2D) #inherit the data index

        #print(self.listy)
        #sys.exit(100)
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
        
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.list2, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox1)

        self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.list2, style=wx.CB_READONLY)
        self.ComboBox2.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox2)

        self.ComboBox3=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.list2, style=wx.CB_READONLY)
        self.ComboBox3.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox3)





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

        self.TOCSY=False




        if('ctocsy' in self.parent.spec.keys() or 'hcconh' in self.parent.spec.keys()): 
            self.TOCSY=True

        listy1=[]
        for key,vals in self.parent.shiftx2.items():
            listy1.append(str(key))

        listy2=list(self.parent.p1to3.keys())


            
        self.ComboBoxI2=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy1, style=wx.CB_READONLY)
        self.ComboBoxI2.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxI2)

        self.ComboBoxI1=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy1, style=wx.CB_READONLY)
        self.ComboBoxI1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxI1)
            
        self.ComboBoxI3=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy1, style=wx.CB_READONLY)
        self.ComboBoxI3.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxI3)


        self.ComboBoxR2=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy2, style=wx.CB_READONLY)
        self.ComboBoxR2.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxR2)

        self.ComboBoxR1=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy2, style=wx.CB_READONLY)
        self.ComboBoxR1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxR1)
            
        self.ComboBoxR3=wx.ComboBox(self, -1, pos=(620, 180), size=(60, -1), choices=listy2, style=wx.CB_READONLY)
        self.ComboBoxR3.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBoxR3)


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

        #self.cb_flip = wx.CheckBox(self, -1,
        #    "Orth",
        #    style=wx.ALIGN_RIGHT)
        #self.Bind(wx.EVT_CHECKBOX, self.on_cb_decon, self.cb_flip)


        #self.cb_decon = wx.CheckBox(self, -1,
        #    "Decon",
        #    style=wx.ALIGN_RIGHT)
        #self.Bind(wx.EVT_CHECKBOX, self.on_cb_decon, self.cb_decon)

        self.cb_grid_auto = wx.CheckBox(self, -1,
            "Labels",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid_auto, self.cb_grid_auto)
        self.cb_grid_auto.SetValue(1)

        


        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)

        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND |  wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)




        print('JJJJJ')

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.vboxC_label = wx.StaticBox(self,-1,'Controls:')
        self.vboxC = wx.StaticBoxSizer(self.vboxC_label, wx.VERTICAL)
        self.vboxC.Add(self.drawbutton, 0, border=3, flag=flags)
        #self.vboxC.Add(self.contourButton, 0, border=3, flag=flags)
        #self.vboxC.Add(self.NOEbutton, 0, border=3, flag=flags)
        #self.vboxC.Add(self.cb_flip)
        #self.vboxC.Add(self.cb_decon)
        self.vboxC.Add(self.cb_grid_auto,0,border=3,flag=flags)
        self.vboxC.Add(self.width1Lab,0,border=3,flag=flags)
        self.vboxC.Add(self.width1Box,0,border=3,flag=flags)
        self.vboxC.Add(self.bigPrevbutton, 0, border=3, flag=flags)
        self.vboxC.Add(self.bigNextbutton, 0, border=3, flag=flags)

        #self.vboxC.Add(self.cb_grid_auto)


        self.hbox.Add(self.vboxC, 0, border=3, flag=flags)

        self.vboxCombo=wx.BoxSizer(wx.VERTICAL)

        self.topLbl = wx.StaticBox(self,-1,'i+1:')
        self.topSizer = wx.StaticBoxSizer(self.topLbl, wx.VERTICAL)


        self.hbox1=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.ComboBox1, 0, border=3, flag=flags)
        self.hbox1.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox1.Add(self.Nbutton, 0, border=3, flag=flags)
       
        self.hbox1.Add(self.ComboBoxR1, 0, border=3, flag=flags)
        self.hbox1.Add(self.ComboBoxI1, 0, border=3, flag=flags)
        #self.hbox1a=wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox1a.Add(self.Upbutton, 0, border=3, flag=flags)
        #self.hbox1a.Add(self.Downbutton, 0, border=3, flag=flags)
        self.topSizer.Add(self.hbox1)




        self.middleLbl = wx.StaticBox(self,-1,'i:')
        self.middleSizer = wx.StaticBoxSizer(self.middleLbl, wx.VERTICAL)
        self.hbox2=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.ComboBox2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Pbutton2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Nbutton2, 0, border=3, flag=flags)
        
        self.hbox2.Add(self.ComboBoxR2, 0, border=3, flag=flags)
        self.hbox2.Add(self.ComboBoxI2, 0, border=3, flag=flags)


        
        self.middleSizer.Add(self.hbox2)


        self.bottomLbl = wx.StaticBox(self,-1,'i-1:')
        self.bottomSizer = wx.StaticBoxSizer(self.bottomLbl, wx.VERTICAL)
        self.hbox3=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.ComboBox3, 0, border=3, flag=flags)
        self.hbox3.Add(self.Pbutton3, 0, border=3, flag=flags)
        self.hbox3.Add(self.Nbutton3, 0, border=3, flag=flags)
        
        self.hbox3.Add(self.ComboBoxR3, 0, border=3, flag=flags)
        self.hbox3.Add(self.ComboBoxI3, 0, border=3, flag=flags)
        self.bottomSizer.Add(self.hbox3)



        self.vboxCombo.Add(self.topSizer)
        self.vboxCombo.Add(self.middleSizer)
        self.vboxCombo.Add(self.bottomSizer)

        self.hbox.Add(self.vboxCombo, 0, border=3, flag=flags)
        self.hbox.Add(self.vboxCntr, 0, border=3, flag=flags)


        # self.hboxAdj=wx.BoxSizer(wx.VERTICAL)
        self.hboxAdj_label = wx.StaticBox(self,-1,'Peaks:')
        self.hboxAdj = wx.StaticBoxSizer(self.hboxAdj_label, wx.VERTICAL)
        self.hboxAdj.Add(self.SelectButton, 0, border=5, flag=flags)
        self.hboxAdj.Add(self.DeselectButton, 0, border=5, flag=flags)
        self.hboxAdj.Add(self.MoveButton, 0, border=5, flag=flags)
        self.hboxAdj.Add(self.AddButton, 0, border=5, flag=flags)
        self.hboxAdj.Add(self.DeleteButton, 0, border=5, flag=flags)

        self.hbox.Add(self.hboxAdj, 0, border=3, flag=flags)

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

        # for ax in self.axes[3:]:
        # self.axes5.set_yticklabels([])  
        
        self.SetSizerAndFit(self.vbox)

        
        #print('Reading projections...')
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


    def onSelectButton(self,event):
        print('Click to select')
        self.SELECT=1
        pass

    def wipe(self):
        for axes in self.axes:
            axes.selection=[]

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
        for i,ind in enumerate(self.index_dataRef):
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




    def GetLevels(self,min_level,fac,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*fac)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels


    def findnear_index(self,test,array):
        #array = numpy.asarray(array)
        idx = (numpy.abs(array - test)).argmin()
        return idx
        #return array[idx]

    #get 2d strips from 3d data
    def ReSlice2d(self,arr,inc,pkl,peak,width,orth=0,lab='hnco'):

        if(orth==0):
            #print(out 2D slice for each peak correlation)
            #print("Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width)
            #print(pkl, lab, self.parent.spec[lab].pkIdx, pkl)
            # exit()
            ptC=self.parent.spec[lab].pkIdx[pkl][0]
            ptC=ptC+inc
            ptH=self.parent.spec[lab].pkIdx[pkl][1]
            #print(ptH, ptC)


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


    def AddSlice(self,sele,ax,spec,xsub,ysub):
        if(spec not in self.parent.spec.keys()):
            #print('cannot find ',spec,'in', self.parent.spec)
            return
        #print('Adding slice %s' % spec)
        #ADD THE SLICE
        if(self.cblist[spec].IsChecked()==False):
            return 0
        levels=self.GetLevels(float(self.textbox1[spec].GetValue()),float(self.textbox2[spec].GetValue()),int(self.textbox3[spec].GetValue()))

        axes=self.axes[ax] #select axes
        #print("SELE:",sele)
        axes.pk=sele
        axes.xsub=xsub
        axes.spec=spec
        axes.peaks=[]
        axes.res=self.index_data[sele]
        
        #turn sele back into the hnco index.
        try:
            

            Xs,Ys,Zs=self.ReSlice2d(self.parent.spec[spec].data,self.inc2,self.AtoI(self.index_data[sele]),self.parent.spec[spec].peak2D,self.Width1,lab=spec)
        except Exception as e:
            print('Exception:',e)
        yave=numpy.average(Ys)


        if(ysub!=0):
            ysa=yave
        else:
            ysa=0


        clbox=self.comboClbox[spec].GetSelection()

        cmap=cm.get_cmap(self.cmaps[clbox])

        axes.contour(Xs-xsub, Ys-ysa+ysub, Zs,levels,cmap=cmap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
        
        if(spec=='hncocanh' or spec=='hncanh'):
            try:
                axes.axvline(x=self.parent.peak[self.index_data[sele]][spec][0].f2)
            except:
                pass

        
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
            xl=(x_max2a-xsub,x_min2a-xsub)
            hl=(yline,yline)
            axes.plot(xl,hl,'k',zorder=0, lw=0.5)





        if(self.cblab[spec].IsChecked()==False):
            return yave

        #print("got to here")
        #do the cross peak labels
        if(self.grid_cb):
            #print('a')
            #print (axes.res)
            #print (self.parent.peak.keys())
            if(spec in self.parent.peak[axes.res].keys()):
                for i,cn in enumerate(self.parent.peak[axes.res][spec]):

                    #print(cn)
                    color='k'
                    fs=8
                    if(i in axes.selection):
                        if(axes.selectspec==spec):
                            color='r'
                            fs=10
                            #axes.scatter(cn.f3,cn.f1-ysa+ysub,c='r',s=50,zorder=2,marker='x')
                            #axes.text(cn.f3,cn.f1-ysa+ysub,cn.name,rotation=90,fontsize=10,color='r')

                    
                    axes.scatter(cn.f3,yline,c=color,s=50,zorder=2,marker='x')
                    # axes.text(cn.f3,cn.f1-ysa+ysub,cn.name,rotation=90,fontsize=fs,color=color)
                    # axes.scatter(cn.f3,cn.f1,c=color,s=50,zorder=2,marker='x')
                    # axes.text(cn.f3,cn.f1,cn.name,rotation=90,fontsize=fs,color=color)

                    axes.peaks.append((cn.f3,cn.f1-ysa+ysub))


                    xl=(cn.f3,cn.f3) #draw arrows
                    if(cn.tp=='main' or cn.tp=='main2' or cn.tp=='PosMax' or cn.tp=='NegMax' or cn.tp=='plus'):
                        cl='r'
                        cl2='red'
                        di='u'
                    elif(cn.tp=='diag'):
                        di = 'x'
                        cl='r'
                        cl2='red'
                    else:
                        cl='b'
                        cl2='blue'
                        di='d'
                    #axes.arrow(xl[0],yline,10.,self.Width1*10., shape='full', lw=0, length_includes_head=True,zorder=2)
                    #if(spec=='hcaco'):
                    

                    if(di=='u'):
                        #up arrow
                        hl=(yline,yline+self.Width1/2.)

                        axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                        axes.arrow(xl[0],yline+self.Width1/2.-self.Width1/9.,0.,self.Width1/10., shape='full', lw=0, length_includes_head=True,zorder=2,color=cl,head_width=0.03*numpy.fabs(Xs[0][0]-Xs[0][-1]),head_length=self.Width1/10.)
                    elif di=='d':
                        hl=(yline-self.Width1/2.,yline)
                        #down arrow
                        axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                        axes.arrow(xl[0],yline-self.Width1/2.+self.Width1/9.,0.,-self.Width1/10., shape='full', lw=0, length_includes_head=True,zorder=2,color=cl,head_width=0.03*numpy.fabs(Xs[0][0]-Xs[0][-1]),head_length=self.Width1/10.)

                    # if     
                    #     xy1 = (pk3.f3,max(0,pk3.inty*factor))

                    #     xym1 = (pk3.f3, axis_above.get_ylim()[0])
                    #     con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
                    #         axesA=axis_above, axesB=axis, color=color, ls=(0,(5,5)), facecolor=None, zorder = 100000)
                    #     # axis_above.add_artist(con)
                    #     axis.add_artist(con)
                        
                    
            

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
        if(spec=='ctocsy' or spec=='hcconh'):
            
            #get row number
            row=int(ax)//int(self.slicesCNT) #get row nynber
            #print('ax',ax)
            #print('spec',spec)
            #print('row:',row)
            #print('cols',self.slicesCNT)

            ##### FIRST TO BMRB ######

            cl2='green'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxR3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxR1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxR2.GetValue()
            

            #print('selected residue:',selRes)
            if(spec=='ctocsy'):
                dicty=self.parent.bmrbC[self.parent.p1to3[selRes]]
                med=self.parent.CTOCSYmed
            elif(spec=='hcconh'):
                dicty=self.parent.bmrbH[self.parent.p1to3[selRes]]
                med=self.parent.HTOCSYmed

            #print('peaks:',dicty)
             
            for key,val in dicty.items():
                #print(key,val)
                #lab=key       #typ
                
                shift=val[0]  #shift
                err=val[1]    #err
                lab= selRes+':'+key

                if(shift>105):
                    continue
                #print(key,val,shift,lab)
                #re-reference all BMRB values to within tocsy range (try aliasing)
                pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                shiftI=pok.ppmI-med   #and save.
                #print('vols:',vol)

                xl=shiftI,shiftI
                axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                axes.text(shift,hl[1],lab,fontsize=8,color='g')
            

            ##### NOW DO TEMPLATE #####
            cl2='orange'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxI3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxI1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxI2.GetValue()
            
            

            #print('selected residue:',selRes)
            if(spec=='ctocsy'):
                try:
                    dicty=self.parent.shiftx2[int(selRes)]['C'] 
                except:
                    dicty={}
                med=self.parent.CTOCSYmed
            elif(spec=='hcconh'):
                
                try:
                    dicty=self.parent.shiftx2[int(selRes)]['H'] 
                except:
                    dicty={}
                med=self.parent.HTOCSYmed

            #print('peaks:',dicty)
             

            for key,val in dicty.items():
                #print(key,val)
                #lab=key       #typ
                
                shift=val
                lab= selRes+self.parent.shiftx2[int(selRes)]['resn']+':'+key
                
                if(shift>105):
                    continue
                #print(key,val,shift,lab)
                #re-reference all BMRB values to within tocsy range (try aliasing)
                pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                shiftI=pok.ppmI-med   #and save.
                #print('vols:',vol)

                xl=shiftI,shiftI
                axes.plot(xl,hl,cl2,lw=0.5,zorder=2)

                hl=(yline-self.Width1/2.,yline+self.Width1/2.*0.8)

                axes.text(shiftI,hl[1],lab,fontsize=8,color='y')
            
 
        elif(spec=='hnca' or spec=='hncoca'):
            
            #get row number
            row=int(ax)//int(self.slicesCNT) #get row nynber
            #print('ax',ax)
            #print('spec',spec)
            #print('row:',row)
            #print('cols',self.slicesCNT)

            ##### FIRST TO BMRB ######

            cl2='green'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxR3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxR1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxR2.GetValue()
            

            #print('selected residue:',selRes)
            
            dicty=self.parent.bmrbC[self.parent.p1to3[selRes]]
            med=self.slicesAdj[spec]

            #print('peaks:',dicty)
             
            for key,val in dicty.items():
                #print(key,val)
                #lab=key       #typ
                
                shift=val[0]  #shift
                err=val[1]    #err
                lab= selRes+':'+key

                if(key!='CA'):
                    continue
                #print(key,val,shift,lab)
                #re-reference all BMRB values to within tocsy range (try aliasing)
                pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                shiftI=pok.ppmI-med   #and save.
                #print('vols:',vol)

                xl=shiftI,shiftI
                axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                axes.text(shift,hl[1],lab,fontsize=8,color='g')
            

            ##### NOW DO TEMPLATE #####
            cl2='orange'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxI3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxI1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxI2.GetValue()
            
            

            

            med=self.slicesAdj[spec]

            #print('peaks:',dicty)
             
            for kk in (0,1):
                #print('selected residue:',selRes)
                try:
                    dicty=self.parent.shiftx2[int(selRes)+kk]['C'] 
                except:
                    dicty={}


                for key,val in dicty.items():
                    #print(key,val)
                    #lab=key       #typ
                    
                    shift=val
                    lab= str(int(selRes)+kk)+self.parent.shiftx2[int(selRes)+kk]['resn']+':'+key
                    
                    if(key!='CA'):
                        continue
                    #print(key,val,shift,lab)
                    #re-reference all BMRB values to within tocsy range (try aliasing)
                    pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                    self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                    shiftI=pok.ppmI-med   #and save.
                    #print('vols:',vol)

                    xl=shiftI,shiftI
                    axes.plot(xl,hl,cl2,lw=0.5,zorder=2)

                    hl=(yline-self.Width1/2.,yline+self.Width1/2.*0.8-self.Width1/2.*kk*0.2)

                    axes.text(shiftI,hl[1],lab,fontsize=8,color='y')
                


                

        elif(spec=='hnco' or spec=='hncaco'):
            
            #get row number
            row=int(ax)//int(self.slicesCNT) #get row nynber
            #print('ax',ax)
            #print('spec',spec)
            #print('row:',row)
            #print('cols',self.slicesCNT)

            ##### FIRST TO BMRB ######

            cl2='green'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxR3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxR1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxR2.GetValue()
            

            #print('selected residue:',selRes)
            
            dicty=self.parent.bmrbC[self.parent.p1to3[selRes]]
            med=self.slicesAdj[spec]

            #print('peaks:',dicty)
             
            for key,val in dicty.items():
                #print(key,val)
                #lab=key       #typ
                
                shift=val[0]  #shift
                err=val[1]    #err
                lab= selRes+':'+key

                if(key!='C'):
                    continue
                #print(key,val,shift,lab)
                #re-reference all BMRB values to within tocsy range (try aliasing)
                pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                shiftI=pok.ppmI-med   #and save.
                #print('vols:',vol)

                xl=shiftI,shiftI
                axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                axes.text(shift,hl[1],lab,fontsize=8,color='g')
            

            ##### NOW DO TEMPLATE #####
            cl2='orange'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxI3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxI1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxI2.GetValue()
            
            

            
            med=self.slicesAdj[spec]

            #print('peaks:',dicty)
             
            for kk in (0,1):
                #print('selected residue:',selRes)
                try:
                    dicty=self.parent.shiftx2[int(selRes)+kk]['C'] 
                except:
                    dicty={}



                for key,val in dicty.items():
                    #print(key,val)
                    #lab=key       #typ
                    
                    shift=val
                    lab= str(int(selRes)+kk)+self.parent.shiftx2[int(selRes)+kk]['resn']+':'+key
                    
                    if(key!='C'):
                        continue
                    #print(key,val,shift,lab)
                    #re-reference all BMRB values to within tocsy range (try aliasing)
                    pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                    self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                    shiftI=pok.ppmI-med   #and save.
                    #print('vols:',vol)

                    xl=shiftI,shiftI
                    axes.plot(xl,hl,cl2,lw=0.5,zorder=2)

                    hl=(yline-self.Width1/2.,yline+self.Width1/2.*0.8-self.Width1/2.*kk*0.2)

                    axes.text(shiftI,hl[1],lab,fontsize=8,color='y')
                


                

        elif(spec=='hncacb' or spec=='hncocacb'):
            
            #get row number
            row=int(ax)//int(self.slicesCNT) #get row nynber
            #print('ax',ax)
            #print('spec',spec)
            #print('row:',row)
            #print('cols',self.slicesCNT)

            ##### FIRST TO BMRB ######

            cl2='green'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxR3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxR1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxR2.GetValue()
            

            #print('selected residue:',selRes)
            
            dicty=self.parent.bmrbC[self.parent.p1to3[selRes]]
            med=self.slicesAdj[spec]

            #print('peaks:',dicty)
             
            for key,val in dicty.items():
                #print(key,val)
                #lab=key       #typ
                
                shift=val[0]  #shift
                err=val[1]    #err
                lab= selRes+':'+key

                if(key!='CB' and key!='CA'):
                    continue
                #print(key,val,shift,lab)
                #re-reference all BMRB values to within tocsy range (try aliasing)
                pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                shiftI=pok.ppmI-med   #and save.
                #print('vols:',vol)

                xl=shiftI,shiftI
                axes.plot(xl,hl,cl2,lw=0.5,zorder=2)
                axes.text(shift,hl[1],lab,fontsize=8,color='g')
            

            ##### NOW DO TEMPLATE #####
            cl2='orange'

            hl=(yline-self.Width1/2.,yline+self.Width1/2.)
            if(row==2):
                selRes=self.ComboBoxI3.GetValue()
            elif(row==0):
                selRes=self.ComboBoxI1.GetValue()
            elif(row==1):
                selRes=self.ComboBoxI2.GetValue()
            
            

            #print('selected residue:',selRes)
            for kk in (0,1):
                try:
                    dicty=self.parent.shiftx2[int(selRes)+kk]['C'] 
                except:
                    dicty={}


                med=self.slicesAdj[spec]

                #print('peaks:',dicty)
                

                for key,val in dicty.items():
                    #print(key,val)
                    #lab=key       #typ
                    
                    shift=val
                    lab= str(int(selRes)+kk)+self.parent.shiftx2[int(selRes)+kk]['resn']+':'+key
                    
                    
                    if(key!='CB' and key!='CA'):
                        continue
                    #print(key,val,shift,lab)
                    #re-reference all BMRB values to within tocsy range (try aliasing)
                    pok=peakEntry(('test','0','0',shift+med,'1')) #create a fake peak entry...
                    self.parent.spec[spec].alias(pok,shift+med,0)  #and alias it to within the range of the spectrum...
                    shiftI=pok.ppmI-med   #and save.
                    #print('vols:',vol)

                    xl=shiftI,shiftI
                    axes.plot(xl,hl,cl2,lw=0.5,zorder=2)

                    hl=(yline-self.Width1/2.,yline+self.Width1/2.*0.8-self.Width1/2.*kk*0.2)

                    axes.text(shiftI,hl[1],lab,fontsize=8,color='y')
                


            pass




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
        #print('done here')
        return yave

    def AddEntry(self,spec,cnt):
        if(spec not in self.parent.spec.keys()):
            return False
        self.slicesInd[spec]=cnt
        try:
            self.slicesTit[cnt]+=' '+spec
        except:
            self.slicesTit[cnt]=spec

        return True

    def SetSlices(self):
        print('Setting up slices...')
        self.slicesInd={}
        self.slicesTit={}
        

        specOrder=[]
        specOrder.append(('hnco','hncaco'))
        specOrder.append(('hnca','hncoca',))
        specOrder.append(('hncacb','hncocacb','cbcaconh'))
        specOrder.append(('hncanh','hncocanh'))
        specOrder.append(('ctocsy',))
        specOrder.append(('hcconh',))
        cnt=0
        for col in specOrder: #go through each cluster
            go=0 #assume we don't have any entries
            for spec in col:  #go through the entires..
                if(self.AddEntry(spec,cnt)): #try to add an entry...
                    go=1   #note if successful
            if(go==1):  #if successful, increment col
                cnt+=1

        self.slicesCNT=cnt

        print('Columns:',cnt)
        print(self.slicesInd)
        print(self.slicesTit)

        from matplotlib.gridspec import GridSpec
        gs1=GridSpec(3,self.slicesCNT)

        #we need 'cnt' columns
        #and three rows.
        self.axes=[] #will store all axes
        AxX={}
        for j in range(3): #for 3 rows
            for i in range(cnt): #for each column
                if(j==0 and i==0): #bottom left is its own person
                    self.axes.append(self.fig.add_subplot(gs1[j,i:(i+1)]))
                    AxY=self.axes[-1]
                    AxX[i]=self.axes[-1]
                elif(i==0):
                    self.axes.append(self.fig.add_subplot(gs1[j,i:(i+1)],sharex=AxX[i]))
                    AxY=self.axes[-1]
                elif(j==0):
                    self.axes.append(self.fig.add_subplot(gs1[j,i:(i+1)],sharey=AxY))
                    AxX[i]=self.axes[-1]
                else:
                    self.axes.append(self.fig.add_subplot(gs1[j,i:(i+1)],sharey=AxY,sharex=AxX[i]))


        
        """
        self.axes2 = self.fig.add_subplot(gs1[1,:3])
        self.axes1 = self.fig.add_subplot(gs1[1,3:6],sharey=self.axes2)
        

        self.axes0 = self.fig.add_subplot(gs1[1,6:9],sharey=self.axes2)


        self.axes5 = self.fig.add_subplot(gs1[0,:3], sharex=self.axes2) 
        
        self.axes4 = self.fig.add_subplot(gs1[0,3:6],sharey=self.axes5, sharex=self.axes1)
        self.axes3 = self.fig.add_subplot(gs1[0,6:9], sharey=self.axes5, sharex=self.axes0)

        self.axes8 = self.fig.add_subplot(gs1[2,:3], sharex=self.axes2)
        self.axes7 = self.fig.add_subplot(gs1[2,3:6], sharey=self.axes8, sharex=self.axes1)
        self.axes6 = self.fig.add_subplot(gs1[2,6:9], sharey=self.axes8, sharex=self.axes0)

        self.axes11 = self.fig.add_subplot(gs1[2,9:12])
        self.axes10 = self.fig.add_subplot(gs1[1,9:12])
        self.axes9 = self.fig.add_subplot(gs1[0,9:12])
        self.axes=self.axes0,self.axes1,self.axes2,self.axes3,self.axes4,self.axes5,self.axes6,self.axes7,self.axes8,self.axes9,self.axes10,self.axes11

        """
        print('Done')
        

        

    def SetAdjust(self):
        adjust={}
        for spec in self.parent.spec.keys():
            if('hncaco'==spec):
                self.parent.HNCOmed=float(self.refbox['hncaco'].GetValue())
                adjust['hncaco']=self.parent.HNCOmed
            elif('hncoca'==spec):
                self.parent.HNCAmed=float(self.refbox['hncoca'].GetValue())
                adjust['hncoca']=self.parent.HNCAmed
            elif('hncocacb'==spec):
                self.parent.HNCACBmed=float(self.refbox['hncocacb'].GetValue())
                adjust['hncocacb']=self.parent.HNCACBmed
            elif('hncacb' ==spec):
                self.parent.HNCAHNCACBmed=float(self.refbox['hncacb'].GetValue())
                adjust['hncacb']=self.parent.HNCAHNCACBmed
            elif('ctocsy'==spec):
                self.parent.CTOCSYmed=float(self.refbox['ctocsy'].GetValue())
                adjust['ctocsy']=self.parent.CTOCSYmed
            elif('hcconh' ==spec):
                self.parent.HTOCSYmed=float(self.refbox['hcconh'].GetValue())
                adjust['hcconh']=self.parent.HTOCSYmed
            else:
                adjust[spec]=0
        
        self.slicesAdj=adjust
        

    def draw_figure(self):
        """ Redraws the figure
        """
        print('drawing walk figure')
       
        #if(self.ax_reset0==0): #possibly wonky
        #    self.x_min0,self.x_max0=self.axes[self.slicesCNT*0].get_xlim()
        #    self.y_min0,self.y_max0=self.axes[self.slicesCNT*0].get_ylim()
        ##if(self.ax_reset1==0):
        #    self.x_min1,self.x_max1=self.axes[self.slicesCNT*1].get_xlim()
        #    self.y_min1,self.y_max1=self.axes[self.slicesCNT*1].get_ylim()
        #if(self.ax_reset2==0):
        #    self.x_min2,self.x_max2=self.axes[self.slicesCNT*2].get_xlim()
        #    self.y_min2,self.y_max2=self.axes[self.slicesCNT*2].get_ylim()

        for ax in self.axes:
            ax.clear()
        
        #write Y labels for first column
        self.sele1=self.ComboBox1.GetSelection()
        self.sele2=self.ComboBox2.GetSelection()
        self.sele3=self.ComboBox3.GetSelection()

        
        self.axes[0*self.slicesCNT].set_ylabel(self.ComboBox1.GetValue()+"(i+1)",fontsize=10)
        self.axes[1*self.slicesCNT].set_ylabel(self.ComboBox2.GetValue()+"(i)",fontsize=10)
        self.axes[2*self.slicesCNT].set_ylabel(self.ComboBox3.GetValue()+"(i-1)",fontsize=10)        
        
        for i in range(self.slicesCNT):
            self.axes[0*self.slicesCNT+i].set_title(self.slicesTit[i],fontsize=10)

        
        for i,ax in enumerate(self.axes):
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if(i%self.slicesCNT==0): #set all Y tics other than far left to false
                pass
            else:
                plt.setp(ax.get_yticklabels(), visible=False)

            if(i<2*self.slicesCNT): #set x tics other than bottom to false
                plt.setp(ax.get_xticklabels(), visible=False)

 
        self.Width1=float(self.width1Box.GetValue())
        if(self.Width1==0):
            self.Width1=1
            self.width1Box.SetValue(str(1))

        self.grid_cb=self.cb_grid_auto.GetValue()
        
        self.SetAdjust() #setup adjustment values for offset
        #print('adjust slices:',self.slicesAdj)
        #middle row, i
        Tots=0,self.slicesCNT,2*self.slicesCNT
        sele=self.sele1,self.sele2,self.sele3

        
        for tot,sel in zip(Tots,sele):
            ysa=self.AddSlice(sel,self.slicesInd[self.parent.refSpec]+tot,self.parent.refSpec,0,0)

            
            #print(tot,sel,ysa)
            for spec in self.parent.spec.keys():
                print(spec)
                if(spec==self.parent.refSpec):
                    continue
                ysa2=self.AddSlice(sel,self.slicesInd[spec]+tot,spec,self.slicesAdj[spec],ysa)
    
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
