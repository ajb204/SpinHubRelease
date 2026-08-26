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
import wx,string,copy,math,numpy,os,re
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import nmrglue as ng
import matplotlib.colors as colors
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor, Cursor

from spinDecon.gui.workspaces.slice2d import AssMan
from spinDecon.gui.context import context_for, project_for

#########################################################
# Plot planes from 4D deconvolution

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8



class SliceFrame4D(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent,tabOne):
    #def __init__(self,uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile):
        wx.Panel.__init__(self, parent=parent)

        self.parent=parent
        self.tabOne=parent.tabOne
        self.app_context = context_for(self.tabOne, parent)
        self.slice_service = self.app_context.slices if self.app_context is not None else None
        if self.slice_service is None:
            from spinDecon.analysis.slice_service import SliceService
            self.slice_service = SliceService(self.tabOne)
        self.state = project_for(self.tabOne, parent)

        #self.peak=self.slice_service.peaks #inherit the data index
        #wx.Frame.__init__(self, None, -1, self.title)
        #copy in the previous variables

        #self.thresh=parent.tabOne.noiseVal
        self.thresh=0.5
        #self.offset=copy.deepcopy(tabOne.offset)
        self.offset=0.0


        #self.spectrumfile=parent.tabOne.spectrumfile


        self.pick_cnt=0
        self.selection=[]
        self.ax_reset=1
        self.ax_reset2=1
        self.ax_resetH=1

        self.inc1=0
        self.inc2=0
        self.SELECT=0
        self.SEARCH=1
        self.OVER=0
        self.MAGMA=0

        #self.create_status_bar()
        self.create_main_panel()

        # Validate the centrally owned plotting-ready plane.
        if self._projection_view_yz() is None:
            labels = self.slice_service.labels
            raise RuntimeError(
                '4D Slices requires the %s.%s raw projection in the shared data store'
                % (labels[2], labels[1])
            )

        self.draw_figure()
        self.canvas.draw()
        self.Show(True)
        self.Fit()


    def _projection_view_yz(self):
        labels = self.slice_service.labels
        return self.slice_service.projection_view(
            labels[2], labels[1], decon=False, transpose='y'
        )

    @property
    def Xs_yz(self): return self._projection_view_yz()['XX']
    @property
    def Ys_yz(self): return self._projection_view_yz()['YY']
    @property
    def Zs_yz(self): return self._projection_view_yz()['ZZ']
    @property
    def Xs_yz_xmin(self): return numpy.min(self.Xs_yz)
    @property
    def Xs_yz_xmax(self): return numpy.max(self.Xs_yz)
    @property
    def Xs_yz_ymin(self): return numpy.min(self.Ys_yz)
    @property
    def Xs_yz_ymax(self): return numpy.max(self.Ys_yz)

    #make an index
    def index(self,array):
        index=[]
        for i in range(len(array)):
            index.append((array[i][0]))
        return index

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        self.fig=Figure()

        self.fig.clear()
        from matplotlib.gridspec import GridSpec
        gs1=GridSpec(2,8)
        self.axes1 = self.fig.add_subplot(gs1[:,0:3])
        self.axes2 = self.fig.add_subplot(gs1[:,3:-2], sharey=self.axes1)
        self.axesH = self.fig.add_subplot(gs1[0,-1])

        self.canvas = FigCanvas(self, -1, self.fig)
        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)




        listy=self.slice_service.peak_names() #setup first combobox with peak list
        listy2=list(listy) #setup second combobox with peak list

        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        if self.ComboBox1.GetCount():
            self.ComboBox1.SetSelection(0)

        self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy2, style=wx.CB_READONLY)
        if self.ComboBox2.GetCount():
            self.ComboBox2.SetSelection(0)

        self.drawbutton = wx.Button(self, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)


        self.I1pbutton = wx.Button(self, -1, "y+",size=(30,-1))
        self.I1mbutton = wx.Button(self, -1, "y-",size=(30,-1))
        self.I2pbutton = wx.Button(self, -1, "x+",size=(30,-1))
        self.I2mbutton = wx.Button(self, -1, "x-",size=(30,-1))
        self.Bind(wx.EVT_BUTTON, self.I1p, self.I1pbutton)
        self.Bind(wx.EVT_BUTTON, self.I1m, self.I1mbutton)
        self.Bind(wx.EVT_BUTTON, self.I2p, self.I2pbutton)
        self.Bind(wx.EVT_BUTTON, self.I2m, self.I2mbutton)


        self.swapbutton = wx.Button(self, -1, "Swap")
        self.Bind(wx.EVT_BUTTON, self.on_swap_button, self.swapbutton)

        self.setoverbutton = wx.Button(self, -1, "SetOverlay")
        self.Bind(wx.EVT_BUTTON, self.on_setoverlay_button, self.setoverbutton)

        self.fishbutton = wx.Button(self, -1, "Fish")
        self.Bind(wx.EVT_BUTTON, self.FishNOE, self.fishbutton)

        self.Pbutton = wx.Button(self, -1, "Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        self.Nbutton = wx.Button(self, -1, "Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton2 = wx.Button(self, -1, "Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button2, self.Pbutton2)

        self.Nbutton2 = wx.Button(self, -1, "Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button2, self.Nbutton2)

        self.NOEbutton = wx.Button(self, -1,"PeakBox")
        self.Bind(wx.EVT_BUTTON, self.on_NOE_button, self.NOEbutton)

        lblList = ['All','NOEs','None',]
        self.rbox = wx.RadioBox(self,label = 'Labels:',choices = lblList ,majorDimension = 0, style = wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX,self.on_cb_grid,self.rbox)
        self.rbox.SetSelection(2)


        self.cb_peak = wx.CheckBox(self, -1,"ProjPeaks",style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_peak)
        self.cb_peak.SetValue(0)



        self.ComboMagma=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=[], style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_magma, self.ComboMagma)
        self.cb_grid_magma = wx.CheckBox(self, -1,
            "Magma",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid_magma, self.cb_grid_magma)
        self.cb_grid_magma.SetValue(0)



        lblList2 = ['Raw','Decon',]
        self.rbox2 = wx.RadioBox(self,label = 'View:',choices = lblList2 ,majorDimension = 1, style = wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX,self.on_cb_grid,self.rbox2)
        self.rbox2.SetSelection(0)

        self.cb_orth = wx.CheckBox(self, -1,"Orth",style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_orth)
        self.cb_orth.SetValue(0)


        #contour text boxes
        self.text1=wx.StaticText(self, -1, 'Min:')
        self.text2=wx.StaticText(self, -1, 'Fac:')
        self.text3=wx.StaticText(self, -1, 'No:')
        self.text4=wx.StaticText(self, -1, 'ProjectMin:')

        #       self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.textbox0 = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER) #min
        self.textbox1 = wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER) #max
        self.textbox2 = wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER) #number
        self.textbox3 = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER) #number

        self.textbox0.SetValue(str(self.slice_service.threshold()))
        self.textbox1.SetValue(str(1.2))
        self.textbox2.SetValue(str(15))
        self.textbox3.SetValue(str(10. * self.slice_service.threshold()))

        #self.textbox0.Bind(wx.EVT_TEXT_ENTER,self.draw_figure())
        #self.textbox1.Bind(wx.EVT_TEXT_ENTER,self.draw_figure())
        #self.textbox2.Bind(wx.EVT_TEXT_ENTER,self.draw_figure())
        #self.textbox3.Bind(wx.EVT_TEXT_ENTER,self.draw_figure())




        self.selectbutton = wx.Button(self, -1, "Select")
        self.Bind(wx.EVT_BUTTON, self.on_select_button, self.selectbutton)
        self.deletebutton = wx.Button(self, -1, "Delete")
        self.Bind(wx.EVT_BUTTON, self.on_delete_button, self.deletebutton)
        self.deselectbutton = wx.Button(self, -1, "Deselect")
        self.Bind(wx.EVT_BUTTON, self.on_deselect_button, self.deselectbutton)
        self.addbutton = wx.Button(self, -1, "Add")
        self.Bind(wx.EVT_BUTTON, self.on_search_button, self.addbutton)



        #self.savebutton = wx.Button(self, -1, "Save List")
        #self.Bind(wx.EVT_BUTTON, self.on_save_button, self.savebutton)



        #self.textbox_savelist = wx.TextCtrl(self,size=(200,-1),style=wx.TE_PROCESS_ENTER)


        #self.textbox_savelist.SetValue(str('out/cross_man_save.out'))



        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)
        #
        # Layout with box sizers
        #
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)


        self.leftLbl = wx.StaticBox(self,-1,'Left:')
        self.leftSizer = wx.StaticBoxSizer(self.leftLbl, wx.VERTICAL)

        lblList3 = ['Raw','Overlay',]
        self.rbox3 = wx.RadioBox(self,label = 'View:',choices = lblList3 ,majorDimension = 1, style = wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX,self.on_cb_grid,self.rbox3)

        self.hbox1=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.ComboBox1, 0, border=3, flag=flags)
        self.hbox1.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox1.Add(self.Nbutton, 0, border=3, flag=flags)

        self.hbox1a=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1a.Add(self.I1pbutton, 0, border=3, flag=flags)
        self.hbox1a.Add(self.I1mbutton, 0, border=3, flag=flags)
        self.hbox1a.Add(self.I2pbutton, 0, border=3, flag=flags)
        self.hbox1a.Add(self.I2mbutton, 0, border=3, flag=flags)
        self.hbox1a.Add(self.swapbutton, 0, border=3, flag=flags)

        self.hbox1b=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1b.Add(self.rbox3, 0, border=3, flag=flags)
        self.hbox1b.Add(self.setoverbutton, 0, border=3, flag=flags)

        self.rightLbl = wx.StaticBox(self,-1,'Right:')
        self.rightSizer = wx.StaticBoxSizer(self.rightLbl, wx.VERTICAL)
        self.hbox2=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.ComboBox2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Pbutton2, 0, border=3, flag=flags)
        self.hbox2.Add(self.Nbutton2, 0, border=3, flag=flags)



        #self.border4 = wx.BoxSizer()
        #self.border4.Add(self.sizerStat)
        self.leftSizer.Add(self.hbox1)
        self.leftSizer.Add(self.hbox1a)
        self.leftSizer.Add(self.hbox1b)
        self.rightSizer.Add(self.hbox2)
        self.rightSizer.Add(self.rbox2)




        self.vboxC = wx.BoxSizer(wx.VERTICAL)
        flags = wx.ALIGN_LEFT | wx.ALL

        self.vboxC.Add(self.drawbutton, 0, border=3, flag=wx.ALIGN_LEFT)
        self.vboxC.Add(self.NOEbutton, 0, border=3, flag=flags)
        self.vboxC.Add(self.cb_orth, 0, border=3, flag=flags)
        self.vboxC.Add(self.cb_peak, 0, border=3, flag=flags)
        self.vboxC.Add(self.cb_grid_magma, 0, border=3, flag=flags)
        self.vboxC.Add(self.fishbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.vboxC)

        self.vboxD=wx.BoxSizer(wx.VERTICAL)
        self.vboxD.Add(self.rbox)
        self.vboxD.Add(self.ComboMagma)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox.Add(self.vboxD, 0, border=3, flag=flags)

        self.hbox.Add(self.leftSizer)
        self.hbox.Add(self.rightSizer)







        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.cntrLbl = wx.StaticBox(self,-1,'Contours:')
        self.cntrSizer = wx.StaticBoxSizer(self.cntrLbl, wx.VERTICAL)

        self.vboxT=wx.BoxSizer(wx.VERTICAL)



        self.hboxA = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.text0, 0, border=3, flag=flags)
        self.hboxA.Add(self.text1, 0, border=3, flag=flags)
        self.hboxA.Add(self.textbox0, 0, border=3, flag=flags)

        self.hboxB = wx.BoxSizer(wx.HORIZONTAL)
        self.hboxB.Add(self.text2, 0, border=3, flag=flags)
        self.hboxB.Add(self.textbox1, 0, border=3, flag=flags)

        self.hboxB.Add(self.text3, 0, border=3, flag=flags)
        self.hboxB.Add(self.textbox2, 0, border=3, flag=flags)

        self.hboxC = wx.BoxSizer(wx.HORIZONTAL)
        self.hboxC.Add(self.text4, 0, border=3, flag=flags)
        self.hboxC.Add(self.textbox3, 0, border=3, flag=flags)

        self.vboxT.Add(self.hboxA, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vboxT.Add(self.hboxB, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vboxT.Add(self.hboxC, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.cntrSizer.Add(self.hboxA)
        self.cntrSizer.Add(self.hboxB)
        self.cntrSizer.Add(self.hboxC)

        self.hbox.Add(self.cntrSizer)
        flags = wx.ALIGN_LEFT | wx.ALL

        self.vboxB=wx.BoxSizer(wx.VERTICAL)
        self.vboxB.Add(self.selectbutton, 0, border=3, flag=flags)
        self.vboxB.Add(self.deletebutton, 0, border=3, flag=flags)
        self.vboxB.Add(self.deselectbutton, 0, border=3, flag=flags)
        self.vboxB.Add(self.addbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.vboxB)

        self.cursor = MultiCursor(self.canvas, (self.axes1, self.axes2),  color = 'k', lw=0.5, useblit=True, horizOn=True)
        # plt.setp(self.axes2.get_yticklabels(), visible=False)

        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.SetSizerAndFit(self.vbox)

        leftHeight = self.leftSizer.GetSize()[1]
        rightHeight = self.rightSizer.GetSize()[1]
        cntrHeight = self.cntrSizer.GetSize()[1]
        heightToSet = max(leftHeight,rightHeight,cntrHeight)
        self.leftSizer.SetMinSize((0,heightToSet))
        self.rightSizer.SetMinSize((0,heightToSet))
        self.cntrSizer.SetMinSize((0,heightToSet))

        self.SetSizerAndFit(self.vbox)

        #self.vbox.Fit(self)


    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    # def readfile(self,infile):
    #     dic, data = ng.pipe.read(self.infile)
    #     return data

    def readfile(self,infile):
        inny=open(infile,'r')
        input=[]
        for line in inny.readlines():
            input.append(line.split())
        inny.close()
        return input



    def GetSpec2D(self,b,decon=False,orth=False,inc1=0,inc2=0,over=False):

        if(over):
            infile=self.overpath+'/out/slice2D/'+self.ComboBox1.GetValue()+'.dat.out'
            print('reading ',infile)
            if(os.path.exists(infile)==0):
                print('Cannot find file:',infile)
                return 0

            input=self.readfile(infile)
            xs=[]
            ys=[]
            zs=[]
            Xs=[]
            Ys=[]
            Zs=[]
            #if(decon=='y'): #choose which column.
            #    col=3
            #else:
            #    col=2
            col=2

            for i in range(len(input)):
                if(len(input[i])!=0):
                    xs.append(float(input[i][0]))
                    ys.append(float(input[i][1]))
                    zs.append(float(input[i][col]))
                else:
                    Xs.append(xs)
                    Ys.append(ys)
                    Zs.append(zs)
                    zs=[]
                    ys=[]
                    xs=[]
            return numpy.array(Xs),numpy.array(Ys),numpy.array(Zs)


        if(decon):
            if(b==1):
                infile='out/slice2D/'+self.ComboBox1.GetValue()+'.dat.decon'
            else:
                infile='out/slice2D/'+self.ComboBox2.GetValue()+'.dat.decon'
            if(os.path.exists(infile)==0):
                print('Cannot find file:',infile)
            else:
                input=self.readfile(infile)
                xs=[]
                ys=[]
                zs=[]
                Xs=[]
                Ys=[]
                Zs=[]
                # if(decon=='y'): #choose which column.
                #     col=3
                # else:
                #     col=2
                col=3

                for i in range(len(input)):
                    if(len(input[i])!=0):
                        xs.append(float(input[i][0]))
                        ys.append(float(input[i][1]))
                        zs.append(float(input[i][col]))
                    else:
                        Xs.append(xs)
                        Ys.append(ys)
                        Zs.append(zs)
                        zs=[]
                        ys=[]
                        xs=[]
                return numpy.array(Xs),numpy.array(Ys),numpy.array(Zs)
            # return self.parent.tabOne.XX,self.parent.tabOne.YY,self.parent.tabOne.decon_data[:,:,ptC+inc2,ptH+inc1]


        if(b==1 or orth==True):
            p=self.ComboBox1.GetSelection()
        else:
            p=self.ComboBox2.GetSelection()

        if(orth==False or b==1):
            ptC=self.slice_service.peak(p).indexK #H
            ptH=self.slice_service.peak(p).indexL #C
            #print inc1,inc2
            #print self.slice_service.axis(2)[ptC]
            #print self.slice_service.axis(3)[ptH]
            #print self.slice_service.peaks[p].x  #proton
            #print self.slice_service.peaks[p].y  #carbon
            #print self.slice_service.peaks[p].ppmI  #proton
            #print self.slice_service.peaks[p].ppmJ  #carbon
            #print self.slice_service.peaks[p].ppmK  #proton
            #print self.slice_service.peaks[p].ppmL  #carbon
            #print p,ptC,ptH
            #print self.slice_service.data[:,:,ptC+self.inc2,ptH+self.inc1]
            return self.slice_service.mesh('XX'), self.slice_service.mesh('YY'), self.slice_service.data[:,:,ptC+inc2,ptH+inc1]
        else:
            ptC=self.slice_service.peak(p).indexI #H
            ptH=self.slice_service.peak(p).indexJ #C
            return self.slice_service.mesh('XX2'), self.slice_service.mesh('YY2'), self.slice_service.data[ptC+inc2,ptH+inc1,:,:]



    def AddLabel(self,axes,cbs1,cbv1,cbs2,orth=False):
        if(len(self.selection)!=0):

            for sele in self.selection:
                cn=self.slice_service.connections[sele]

                pkindx={}
                for i,pk in enumerate(self.slice_service.peaks):
                    pkindx[pk.name]=i
                if(orth):
                    x=self.slice_service.peaks[pkindx[cn.p2]].ppmK #alias
                    y=self.slice_service.peaks[pkindx[cn.p2]].ppmL #alias
                else:
                    x=self.slice_service.peaks[pkindx[cn.p2]].ppmI #alias
                    y=self.slice_service.peaks[pkindx[cn.p2]].ppmJ #alias

                axes.scatter(x,y,c='k',s=50,marker='x',zorder=2)
                axes.text(x,y,cn.p2,color='r')


        if(orth):
            pkindx={}
            for i,pk in enumerate(self.slice_service.peaks):
                pkindx[pk.name]=i

        if(self.rbox.GetSelection()==1):#if NOESY list is present
            for cn in self.slice_service.connections:
                if(cn.v1==cbs1):
                    if(orth):
                        x=self.slice_service.peaks[pkindx[cn.p2]].ppmK #alias
                        y=self.slice_service.peaks[pkindx[cn.p2]].ppmL #alias
                    else:
                        x=cn.f4
                        y=cn.f3

                    axes.scatter(x,y,c='k',s=50,marker='x',zorder=2)

                    if(len(self.selection)==0):
                        if(cn.v2==cbs2):
                            axes.text(x,y,cn.p2,color='r')
                        else:
                            axes.text(x,y,cn.p2,fontsize=8)
        else:
            if(orth):
                pk=self.slice_service.peaks[cbs1]
                axes.text(pk.ppmK,pk.ppmL,pk.name,color='b')
                axes.scatter(pk.ppmK,pk.ppmL,s=10,c='b',marker='x',zorder=2)
                #print self.cb_noes.IsChecked()
                if(self.rbox.GetSelection()==0):
                    for i,pk in enumerate(self.slice_service.peaks):
                        if(i!=cbs1):
                            axes.scatter(pk.ppmK,pk.ppmL,s=10,c='k',marker='x',zorder=2)
                            self.line=pk.name
                            axes.text(pk.ppmK,pk.ppmL,self.line,rotation=0,fontsize=7)
            else:
                pk=self.slice_service.peaks[cbs1]
                axes.text(pk.ppmI,pk.ppmJ,pk.name,color='b')
                axes.scatter(pk.ppmI,pk.ppmJ,s=10,c='b',marker='x',zorder=2)
                #print self.cb_noes.IsChecked()
                if(self.rbox.GetSelection()==0):
                    for i,pk in enumerate(self.slice_service.peaks):
                        if(i!=cbs1):
                            axes.scatter(pk.ppmI,pk.ppmJ,s=10,c='k',marker='x',zorder=2)
                            self.line=pk.name
                            axes.text(pk.ppmI,pk.ppmJ,self.line,rotation=0,fontsize=7)


        if(self.MAGMA):
            label=self.slice_service.peaks[cbs1].name
            if(label in list(self.peak_dict.keys())): #if this residue can be translated into magmaPDB language...
                print('current peakID:',label)   #peak label
                magID=self.peak_dict[label]     #magma nmr id
                if(magID in list(self.result_dict.keys())): #for each assignment
                    print('magmaIDs:',magID)         #magma nmr id
                    print('nmr:',magID,' pdb_assignments: ',self.result_dict[magID])
                    print(self.ComboMagma.GetSelection())
                    currass=self.result_dict[magID][self.ComboMagma.GetSelection()] #get proposed assignment from combobox
                else:
                    print('Not assigned')
                    currass=self.ComboMagma.GetValue()

                print('curass',currass)
                adjacency=self.pdb[currass] #get PDB assignment adjacency

                labs={}
                for adj in adjacency:       #for each PDB adjacency...
                    if(adj in list(self.result_dictInv.keys())): #if the adjacency PDB has been assigned...
                        #print 'ass:',adj,' assigned to ',self.result_dictInv[adj]
                        newlab=''
                        for resy in self.result_dictInv[adj]:
                            for pk in self.peak_dictInv[resy]:
                                ref=self.indy[pk]
                                if(ref not in list(labs.keys())):
                                    labs[ref]=adj
                                else:
                                    labs[ref]+='/'+adj
                for ref in list(labs.keys()):
                    #ref=self.indy[pk]
                    #xval=self.slice_service.axis(0)[self.slice_service.peaks[ref].indexI]
                    #yval=self.slice_service.axis(1)[self.slice_service.peaks[ref].indexJ]
                    #zval=self.slice_service.axis(2)[self.slice_service.peaks[ref].indexK]
                    #print 'plotting:',adj,ref,xval,yval,zval

                    axes.scatter(self.slice_service.peaks[ref].ppmK,self.slice_service.peaks[ref].ppmL,s=10,c='g',marker='x',zorder=2)
                    axes.text(self.slice_service.peaks[ref].ppmK,self.slice_service.peaks[ref].ppmL,labs[ref],rotation=0,fontsize=8,color='g',va='bottom',ha='right')
                    #self.axes1.scatter(yval,dimH,c='g',s=50,zorder=2,marker='x')
                    #axes.text(yval,dimH+Width1*0.1,labs[ref],rotation=90,fontsize=8,color='g',va='bottom')




    def FishNOE(self,event):
        #for each residue
        if(self.MAGMA==0):
            print('click magma button to use this mode')
            return

        print('Fishing for new peaks')
        addNew=0
        for pk in self.slice_service.peaks:
            label=pk.name
            if(label in list(self.peak_dict.keys())): #if this residue can be translated into magmaPDB language...
                print('current peakID:',label)   #peak label
                magID=self.peak_dict[label]     #magma nmr id
                if(magID in list(self.result_dict.keys())): #for each assignment
                    print('magmaIDs:',magID)         #magma nmr id
                    print('nmr:',magID,' pdb_assignments: ',self.result_dict[magID])
                    if(len(self.result_dict[magID])==1): #is this peak unqiue assigned? confident assignment
                        currass=self.result_dict[magID][0] #get proposed assignment from combobox

                        adjacency=self.pdb[currass] #get PDB assignment adjacency

                        labs={}
                        for adj in adjacency:       #for each PDB adjacency...
                            if(adj in list(self.result_dictInv.keys())): #if the adjacency PDB has been assigned...
                                if(len(self.result_dictInv[adj])==1):
                                    #print 'ass:',adj,' assigned to ',self.result_dictInv[adj]
                                    newlab=''
                                    resy=self.result_dictInv[adj][0]
                                    for pknew in self.peak_dictInv[resy]:


                                        print('testing:',label,pknew)
                                        add,i1,i2=self.NOETest(self.indy[label],self.indy[pknew]) #test possible NOE assignment
                                        if(add):
                                            addNew+=1
                                            self.AddNOE(self.indy[label],self.indy[pknew],i1,i2)
                                        #ref=self.indy[pk]
                                        #if(ref not in labs.keys()):
                                        #    labs[ref]=adj
                                        #else:
                                        #    labs[ref]+='/'+adj
        print()
        print('Fished out:',addNew,'pairs of peaks')
        print()
        self.selection=[]
        self.draw_figure()




    #test whether given selection has reasonable intensity
    #in the reciprocated plane
    def NOETest(self,refy,maxy):
        print('Testing for NOE:')
        add=0

        #pk1=self.slice_service.peaks[self.ComboBox1.GetSelection()]

        pk1=self.slice_service.peaks[refy]
        pk2=self.slice_service.peaks[maxy]
        if(pk1.name==pk2.name):
            print('this is a diagonal peak. Not adding')
            add=1

        print(pk1.name,pk2.name)


        #is this already here?

        nummy=0
        for i,cn in enumerate(self.slice_service.connections):
            if(cn.p1==pk1.name and cn.p2==pk2.name):
                print('NOE already present:',cn.tag)
                nummy+=1
            if(cn.p2==pk1.name and cn.p1==pk2.name):
                print('NOE (reciprocated) already present:',cn.tag)
                nummy+=1
        if(nummy==2):
            print('both peaks are here: not adding.')
            add=1
        #i1=self.slice_service.data[pk1.indexI,pk2.indexJ,pk2.indexK]
        #i2=self.slice_service.data[pk2.indexI,pk1.indexJ,pk1.indexK]

        i1=self.slice_service.data[pk1.indexI,pk1.indexJ,pk2.indexK,pk2.indexL]
        i2=self.slice_service.data[pk2.indexI,pk2.indexJ,pk1.indexK,pk1.indexL]

        print('i1,i2: ',i1,i2)

        #print self.parent.tabOne.XX[pk1.indexI,pk2.indexJ,pk2.indexK]
        #print self.parent.tabOne.YY[pk1.indexI,pk2.indexJ,pk2.indexK]
        #print self.slice_service.data[pk1.indexI,pk2.indexJ,pk2.indexK]
        #print self.slice_service.data[pk2.indexI,pk2.indexJ,pk1.indexK]
        #print self.parent.tabOne.ZZ[pk1.indexI,pk2.indexJ,pk2.indexK]
        #print self.parent.tabOne.XX[pk2.indexI,pk1.indexJ,pk1.indexK]
        #print self.parent.tabOne.YY[pk2.indexI,pk1.indexJ,pk1.indexK]
        #print self.parent.tabOne.ZZ[pk2.indexI,pk1.indexJ,pk1.indexK]
        #print self.slice_service.datadec[pk2.indexI,pk1.indexJ,pk1.indexK]

        #noise=self.parent.tabOne.dmax*float(self.parent.tabOne.threshBox.GetValue())
        self.noise=float(self.textbox0.GetValue())
        print('noise:',self.noise)
        print('Intensity of selected cross peak:')
        print(i1,i1/self.noise)
        print('Intensity of reciprocated cross peak:')
        print(i2,i2/self.noise)
        if(numpy.fabs(i1)<self.noise):
            print('I1 less than noise. Probably not a cross peak.')
            add=1
        if(numpy.fabs(i2)<self.noise):
            print('I2 less than noise. Probably not a cross peak.')
            add=1
        if(add==0): #if the signs are opposite
            if(i1>0 and i2<0):
                add=1
            elif(i1<0 and i2>0):
                add=1

        print('Ratio (expected to be >50%):')
        print(min(i1,i2)/max(i1,i2)*100.)

        if(add):
            return 0,0,0

        #try: #adjust intensity for deconvolved spectrum
        #    i3=self.slice_service.datadec[pk1.indexI,pk2.indexJ,pk2.indexK]
        #    i4=self.slice_service.datadec[pk2.indexI,pk1.indexJ,pk1.indexK]
        #    i1=i1-i3
        #    i2=i2-i4
        #except:
        #    pass
        print('intensity1:',i1)
        print('intensity2:',i2)

        return 1,i1,i2




    #add NOE entry with given intensities
    #check to make sure main and reciprocal are not already there.
    def AddNOE(self,refy,maxy,i1,i2):
        print('Adding NOE:')

        #pk1=self.slice_service.peaks[self.ComboBox1.GetSelection()]
        pk1=self.slice_service.peaks[refy]
        pk2=self.slice_service.peaks[maxy]



        self.selection=[]
        from spinDecon.domain.peaks import connEntry
        skip=0
        for i,cn in enumerate(self.slice_service.connections):
            if(cn.p1==pk1.name and cn.p2==pk2.name):
                print('NOE already present:',cn.tag)
                skip=1
        if(skip==0):
            print('added main')
            #stry=('%s\t%s\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f') % (pk1.name,pk2.name,pk1.ppmK,pk1.ppmJ,pk2.ppmI,i1,i2,i1/self.noise,i2/self.noise,0,0,0,0)
            stry=('%s\t%s\t%f\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f\t%f\t%f') % (pk1.name,pk2.name,pk1.ppmJ,pk1.ppmI,pk2.ppmL,pk2.ppmK,i1,i2,i1/self.noise,i2/self.noise,0,0,0,0,0,0)
            test=stry.split()
            print(stry)
            cnNew=connEntry(test,sym='y',peak=self.slice_service.peaks,dim=4)
            self.slice_service.connections.append(cnNew)
            self.selection.append(len(self.slice_service.connections)-1)


        skip=0
        for i,cn in enumerate(self.slice_service.connections):
            if(cn.p2==pk1.name and cn.p1==pk2.name):
                print('NOE (reciprocated) already present:',cn.tag)
                skip=1
        if(skip==0):
            print('added reciprocal')
            #stry=('%s\t%s\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f') % (pk2.name,pk1.name,pk2.ppmK,pk2.ppmJ,pk1.ppmI,i2,i1,i2/self.noise,i1/self.noise,0,0,0,0)
            stry=('%s\t%s\t%f\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f\t%f\t%f') % (pk2.name,pk1.name,pk2.ppmJ,pk2.ppmI,pk1.ppmL,pk1.ppmK,i2,i1,i2/self.noise,i1/self.noise,0,0,0,0,0,0)
            test=stry.split()
            print(stry)
            cnNew=connEntry(test,sym='y',peak=self.slice_service.peaks,dim=4)
            self.slice_service.connections.append(cnNew)
            self.selection.append(len(self.slice_service.connections)-1)

        print(self.selection)
        print(len(self.slice_service.connections))


    


    def GetLevels(self,max_level,min_level,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels

    def draw_figure(self):
        """ Redraws the figure
        """



        if(self.ax_resetH==0):
            x_minH,x_maxH=self.axesH.get_xlim()
            y_minH,y_maxH=self.axesH.get_ylim()
        if(self.ax_reset==0):
            x_min1,x_max1=self.axes1.get_xlim()
            y_min1,y_max1=self.axes1.get_ylim()
        if(self.ax_reset2==0):
            x_min2,x_max2=self.axes2.get_xlim()
            y_min2,y_max2=self.axes2.get_ylim()

        orth=self.cb_orth.IsChecked()




        max_level=float(self.textbox1.GetValue())
        min_level=float(self.textbox0.GetValue())
        ctr_level=int(self.textbox2.GetValue())
        levels=self.GetLevels(max_level,min_level,ctr_level)

        max_levelP=float(self.textbox1.GetValue())
        min_levelP=float(self.textbox3.GetValue())
        ctr_levelP=int(self.textbox2.GetValue())
        levelsP=self.GetLevels(max_levelP,min_levelP,ctr_levelP)

        self.axes1.clear()
        Xs,Ys,Zs=self.GetSpec2D(1,inc1=self.inc1,inc2=self.inc2)


        colormap=cm.seismic
        #self.axes.title('2D plot for '+self.ComboBox1.GetValue()+' at '+str(self.slice_service.peaks[self.ComboBox1.GetSelection()][2])+' 13C ppm and '+str(self.slice_service.peaks[self.ComboBox1.GetSelection()][1])+' 1H ppm')
        self.axes1.contour(Xs, Ys, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
#        self.fig.set_xlabel('Omega_C (ppm)')
#$        p.set_ylabel('Omega_H (ppm)')
#        self.ax.set_zlabel('Z')

        if(self.rbox3.GetSelection()==1 and self.OVER==1):
            levelsOver=self.GetLevels(max_level,self.minLevelOver,ctr_level)
            print(levelsOver)
            Xo,Yo,Zo=self.GetSpec2D(1,inc1=0,inc2=0,over=True)
            try:
                self.axes1.contour(Xo+self.overX, Yo+self.overY, Zo,levelsOver,cmap=cm.Blues,norm=colors.Normalize(vmin=-numpy.max(levelsOver),vmax=numpy.max(levelsOver))) #plot pdb network
            except:
                pass

        self.axes1.set_xlabel(self.slice_service.labels[0] + r' $\omega$ (ppm)')
        self.axes1.set_ylabel(self.slice_service.labels[1] + r' $\omega$ (ppm)')
        self.AddLabel(self.axes1,self.ComboBox1.GetSelection(),self.ComboBox1.GetValue(),self.ComboBox2.GetSelection())

        if(self.ax_reset==1):
            y_min1=Ys[0][0]
            y_max1=Ys[0][(len(Ys[0]))-1]
            x_min1=Xs[0][0]
            x_max1=Xs[(len(Xs))-1][0]
            self.axes1.set_xlim(x_min1,x_max1)
            self.axes1.set_ylim(y_min1,y_max1)
            self.ax_reset=0
        else:
            self.axes1.set_xlim(x_min1,x_max1)
            self.axes1.set_ylim(y_min1,y_max1)




        self.axes2.clear()
        # self.axes2.set_yticklabels([])

        

        if(orth==True):
            self.axes2.set_xlabel(self.slice_service.labels[2] + r' $\omega$ (ppm)')
            # self.axes2.set_ylabel(self.slice_service.labels[3]+' $\omega$ (ppm)')
        else:
            self.axes2.set_xlabel(self.slice_service.labels[0] + r' $\omega$ (ppm)')
            self.axes2.set_ylabel(self.slice_service.labels[1] + r' $\omega$ (ppm)')

        if(self.rbox2.GetSelection()!=1): #if not dec
            Xt,Yt,Zt=self.GetSpec2D(2,orth=orth, inc1=self.inc1, inc2=self.inc2)
            self.axes2.contour(Xt, Yt, Zt,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
            if(orth):
                self.AddLabel(self.axes2,self.ComboBox1.GetSelection(),self.ComboBox1.GetValue(),self.ComboBox2.GetSelection(),orth=orth)
            else:
                self.AddLabel(self.axes2,self.ComboBox2.GetSelection(),self.ComboBox2.GetValue(),self.ComboBox1.GetSelection())


        else:
            Xt,Yt,Zt=self.GetSpec2D(1,decon=True,orth=orth, inc1=self.inc1,inc2=self.inc2)
            self.axes2.contour(Xt, Yt, Zt,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network
            self.AddLabel(self.axes2,self.ComboBox1.GetSelection(),self.ComboBox1.GetValue(),self.ComboBox2.GetSelection())

        if(self.ax_reset2==1):
            y_min2=Yt[0][0]
            y_max2=Yt[0][(len(Yt[0]))-1]
            x_min2=Xt[0][0]
            x_max2=Xt[(len(Xt))-1][0]
            self.axes2.set_xlim(x_min2,x_max2)
            self.axes2.set_ylim(y_min2,y_max2)

            self.ax_reset2=0
        else:
            self.axes2.set_xlim(x_min2,x_max2)
            self.axes2.set_ylim(y_min2,y_max2)




        ##########################################################################

        #DO HC PROJECTION
        p=self.ComboBox1.GetSelection()
        dimH=self.slice_service.peaks[p].ppmK
        dimC=self.slice_service.peaks[p].ppmL

        dimHs=self.slice_service.axis(2)[self.slice_service.peaks[p].indexK+self.inc2] #H
        dimCs=self.slice_service.axis(3)[self.slice_service.peaks[p].indexL+self.inc1] #C

        label=self.slice_service.peaks[p].name
        #print dimH,dimC
        print(dimHs,dimCs)

        # Subplot 2 - the HC projection
        #self.axes = self.fig.add_subplot(322)


        self.axesH.clear()
        self.axesH.contour( self.Ys_yz, self.Xs_yz, self.Zs_yz,levelsP,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levelsP),vmax=numpy.max(levelsP))) #plot HC plane


        ex=(self.Xs_yz_ymin,self.Xs_yz_ymax) #x horiz
        ey=(dimCs,dimCs) #line across X
        self.axesH.plot( ex, ey, c='g') #plot pdb network
        ey=(dimHs,dimHs) #line along Y
        ex=(self.Xs_yz_xmin,self.Xs_yz_xmax) #x horiz
        self.axesH.plot( ey, ex, c='g') #plot pdb network
        self.axesH.text(dimH,dimC,label,fontsize=8)
        self.axesH.scatter(dimH,dimC,c='k',s=50,marker='x',zorder=2)

        if(self.rbox2.GetSelection()==1): #draw second line
            p2=self.ComboBox2.GetSelection()
            dimH2=self.slice_service.peaks[p2].ppmK
            dimC2=self.slice_service.peaks[p2].ppmL
            label2=self.slice_service.peaks[p2].name
            dimHs2=self.slice_service.axis(2)[self.slice_service.peaks[p2].indexK] #H
            dimCs2=self.slice_service.axis(3)[self.slice_service.peaks[p2].indexL] #C
            ex=(self.Xs_yz_ymin,self.Xs_yz_ymax) #x horiz
            ey=(dimCs2,dimCs2) #line across X
            self.axesH.plot( ex, ey, c='b') #plot pdb network
            ey=(dimHs2,dimHs2) #line along Y
            ex=(self.Xs_yz_xmin,self.Xs_yz_xmax) #x horiz
            self.axesH.plot( ey, ex, c='b') #plot pdb network
            self.axesH.text(dimH2,dimC2,label2,fontsize=8)
            self.axesH.scatter(dimH2,dimC2,c='k',s=50,marker='x',zorder=2)



        for peak in self.slice_service.peaks: #draw peak locations
            self.axesH.scatter(peak.ppmK,peak.ppmL,c='k',s=10,marker='x',zorder=2)
        if(self.cb_peak.IsChecked()): #write on text labels if required
            for peak in self.slice_service.peaks:
                self.axesH.text(peak.ppmK,peak.ppmL,peak.name,fontsize=6)

        self.axesH.set_xlabel(self.slice_service.labels[2]+'(ppm)',fontsize=8)
        self.axesH.set_ylabel(self.slice_service.labels[3]+'(ppm)',fontsize=8)

        if(self.ax_resetH==1):
            #flip axes (NMR)
            #y_minH=float(self.Ys_yz[0][0])
            #y_maxH=float(self.Ys_yz[(len(self.Ys_yz))-1][0])
            #x_minH=float(self.Xs_yz[0][0])
            #x_maxH=float(self.Xs_yz[0][(len(self.Xs_yz[0]))-1])
            x_minH=float(self.Ys_yz[0][0])
            x_maxH=float(self.Ys_yz[(len(self.Ys_yz))-1][0])
            y_minH=float(self.Xs_yz[0][0])
            y_maxH=float(self.Xs_yz[0][(len(self.Xs_yz[0]))-1])
            self.axesH.set_xlim(x_minH,x_maxH)
            self.axesH.set_ylim(y_minH,y_maxH)
            self.ax_resetH=0
        else:
            self.axesH.set_xlim(x_minH,x_maxH)
            self.axesH.set_ylim(y_minH,y_maxH)


#        plt.scatter(float(self.slice_service.peaks[k][2]),float(self.slice_service.peaks[self.ComboBox1.GetSelection()][1]),c='k',s=50)
        """
        for k in range(len(self.slice_service.peaks)):#look through the peak list...
            self.axes.scatter(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),c='k',s=10)
            if(self.cb_grid.IsChecked() or self.cb_grid_auto.IsChecked()):
                tick=0
                for j in range(len(self.slice_service.connections)):#search through the cross peaks
                    if(k==self.slice_service.connections[j][1] and self.ComboBox1.GetSelection()==self.slice_service.connections[j][0]):#to find the reference for the reciprocated cross peak...
                        tick=1
                    if(k==self.slice_service.connections[j][0] and self.ComboBox1.GetSelection()==self.slice_service.connections[j][1]):#to find the reference for the reciprocated cross peak...
                        tick=1
                if(tick==1 and self.cb_grid_auto.IsChecked()):
                    if(len(self.selection)==0):#if there is no selection, print all peak labels
                        self.axes.text(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),self.index_data[k],fontsize=8)
                    else:
                        for l in range(len(self.selection)):#if there is a selection, only show labels for these
                            if(self.selection[l]==k):
                                self.axes.text(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),self.index_data[k],color='r',fontsize=8)
                elif(tick==0 and self.cb_grid.IsChecked() and k!=self.ComboBox1.GetSelection()):
                    self.axes.text(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),self.index_data[k],color='k',fontsize=7)

        self.axes.text(float(self.slice_service.peaks[self.ComboBox1.GetSelection()][1]),float(self.slice_service.peaks[self.ComboBox1.GetSelection()][2]),self.ComboBox1.GetValue())
        self.axes.scatter(float(self.slice_service.peaks[self.ComboBox1.GetSelection()][1]),float(self.slice_service.peaks[self.ComboBox1.GetSelection()][2]),c='r',s=50)
        """


#        if(self.cb_grid_auto.GetValue()==1):#if auto-detect has been run...
#            for j in range(len(self.slice_service.connections)):#search through the cross peaks
#                if(self.index_data[self.slice_service.connections[j][0]]==self.ComboBox1.GetValue()):#find references for the current plane...
#                    for k in range(len(self.slice_service.peaks)):#look through the peak list...
#                        if(self.slice_service.peaks[k][0]==self.index_data[self.slice_service.connections[j][1]]):#to find the reference for the reciprocated cross peak...
#                            #plt.scatter(float(self.slice_service.peaks[k][2]),float(self.slice_service.peaks[self.ComboBox1.GetSelection()][1]),c='k',s=50)
#                            plt.scatter(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),c='k',s=50)
#
#                            if(len(self.selection)==0):#if there is no selection, print all peak labels
#                                plt.text(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),self.index_data[self.slice_service.connections[j][1]])
#                            else:
#                                for l in range(len(self.selection)):#if there is a selection, only show labels for these
#                                    if(self.selection[l]==self.slice_service.connections[j][1]):
#                                        plt.text(float(self.slice_service.peaks[k][1]),float(self.slice_service.peaks[k][2]),self.index_data[self.slice_service.connections[j][1]],color='r')




#str(self.slice_service.peaks[self.sum[0]][0])+'.dat.out')




        # clear the axes and redraw the plot anew
        #
#        self.axes.clear()
#        self.axes.grid(self.cb_grid.IsChecked())
#        self.axes.bar(
#            left=x,
#            height=self.data,
#            width=self.slider_width.GetValue() / 100.0,
#            align='center',
#            alpha=0.44,
#            picker=5)
        # multi = MultiCursor(self.canvas, (self.axes1, self.axes2), color='k', lw=0.5)
        # cursor = Cursor(self.axes1, color = 'k', lw=0.5)
        # plt.setp(self.axes2.get_yticklabels(), visible=False)
        self.canvas.draw()

    def on_setoverlay_button(self,event):
        self.minLevelOver=float(self.textbox0.GetValue())

        self.overpath=self.onGetDir(True)
        #'/Users/futileenterprises/nmrNOESY/Flemming/dataForAndy/new/set_2/LV_double_methyl/150ms_HCCH_NOESY_short'
        print('path to overlay folder:',self.overpath)

        if(self.overpath!=-1):
            self.OVER=1
        else:
            print('Cancelling overlay')
            self.OVER=0
            self.draw_figure()
            return
        dlg=wx.TextEntryDialog(self,"Enter Min Contour level:","Set MinContour")
        dlg.CentreOnParent()
        if(dlg.ShowModal() == wx.ID_OK):
            self.minLevelOver=float(dlg.GetValue())
        else:
            print('Using current level')

        dlg=wx.TextEntryDialog(self,"Enter ppm offset in X:","Set X offset")
        dlg.CentreOnParent()
        if(dlg.ShowModal() == wx.ID_OK):
            self.overX=float(dlg.GetValue())
        else:
            print('Setting to zero')
            self.overX=0

        dlg=wx.TextEntryDialog(self,"Enter ppm offset in Y:","Set Y offset")
        dlg.CentreOnParent()
        if(dlg.ShowModal() == wx.ID_OK):
            self.overY=float(dlg.GetValue())
        else:
            print('Setting to zero')
            self.overY=0




        print('overPath: ',self.overpath)
        print('overMin:  ',self.minLevelOver)
        print('overX:    ',self.overX)
        print('overY:    ',self.overY)
        self.draw_figure()

    #FGA added
    def onGetDir(self, e):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.DirDialog(self, message="Choose a folder",style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print("You chose the following file(s):")
            print(path)
        else:
            print('no file selected')
            path =-1
        dlg.Destroy()
        return path

    def on_cb_grid(self, event):
        self.draw_figure()


    def on_cb_grid_magma(self, event):
        if(self.MAGMA==1):
            self.MAGMA=0
            self.draw_figure()
            return

        from magma.magma import Magma
        from magma.analysis import AnalAll
        self.inst=Magma(self.slice_service.decon_parameter_file,run='n') #get instance of magma
        self.result_dict,nsoln,edges,edgestot,edgesG2=AnalAll(self.inst.P.outdir)
        #result_dict is nmr->pdb
        self.result_dictInv={} #create the inverse dictionary: pdb->nmr
        for key,vals in list(self.result_dict.items()):
            for val in vals:
                if(val not in list(self.result_dictInv.keys())):
                    self.result_dictInv[val]=[]
                if(key not in self.result_dictInv[val]):
                    self.result_dictInv[val].append(key)


        self.peak_dict={}  #self.peak->nmr
        for i in range(len(self.slice_service.peaks)):
            resi=re.findall(r'[0-9]+',self.slice_service.peaks[i].name)[0]
            labadj=resi+self.slice_service.peaks[i].name[0]
            self.peak_dict[self.slice_service.peaks[i].name]=labadj
        self.peak_dictInv={} #nmr -> self.peak
        for key,val in list(self.peak_dict.items()):
            if(val not in list(self.peak_dictInv.keys())):
                self.peak_dictInv[val]=[]
            if(key not in self.peak_dictInv[val]):
                self.peak_dictInv[val].append(key)
        #print self.peak_dictInv
        self.pdb={}
        for i,node in enumerate(self.inst.short_node_list):
            if node not in list(self.pdb.keys()):
                self.pdb[node]=[]
            for adj in self.inst.short_adjacency[i]:
                if(adj not in self.pdb[node]):
                    self.pdb[node].append(adj)

        self.indy={}
        for i,pk in enumerate(self.slice_service.peaks):
            self.indy[pk.name]=i

        self.UpdateMagmaBox()
        self.MAGMA=1
        self.draw_figure()



    def on_cb_decon(self, event):
        self.draw_figure()

    def on_cb_grid_auto(self, event):
        self.draw_figure()

    def I1p(self,event):
        self.inc1+=1
        self.draw_figure()
    def I1m(self,event):
        self.inc1-=1
        self.draw_figure()
    def I2p(self,event):
        self.inc2+=1
        self.draw_figure()
    def I2m(self,event):
        self.inc2-=1
        self.draw_figure()

    def on_swap_button(self,event):
        tmp=self.ComboBox1.GetSelection()
        self.ComboBox1.SetSelection(self.ComboBox2.GetSelection())
        self.ComboBox2.SetSelection(tmp)
        self.draw_figure()

    def on_draw_magma(self, event):
        self.draw_figure()

    def on_draw_button(self, event):
        self.ax_reset=1
        self.ax_reset2=1
        self.ax_resetH=1
        self.inc1=0
        self.inc2=0
        if(self.MAGMA):
            self.UpdateMagmaBox()
        self.draw_figure()

    def on_NOE_button(self, event):
        bool=AssMan(self)

    def getnoise_spec(self,X,Y,Z):
        #self.nxmin=float(self.textbox_xmin.GetValue())
        #self.nxmax=float(self.textbox_xmax.GetValue())
        #self.nymin=float(self.textbox_ymin.GetValue())
        #self.nymax=float(self.textbox_ymax.GetValue())
        #mask1=(X>self.nxmin)*(X<self.nxmax)
        #mask2=(Y>self.nymin)*(Y<self.nymax)
        #maxnoise=numpy.max(Z[(mask1*mask2)])
        #print maxnoise
        #self.textbox0.SetValue(str(maxnoise*float(self.textbox_fac.GetValue())))
        pass

    def on_P_button(self, event):
        self.ax_reset=1
        #self.ax_resetH=1
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        if(self.MAGMA):
            self.UpdateMagmaBox()
        self.selection=[]
        self.draw_figure()


    def UpdateMagmaBox(self):
        label=self.slice_service.peaks[self.ComboBox1.GetSelection()].name
        #print 'current peakID:',label
        magID=self.peak_dict[label]
        #print 'magmaIDs:',magID
        #print 'nmr:',magID,' pdb_assignments: ',self.result_dict[magID]
        self.ComboMagma.Clear()
        if(magID in self.result_dict):
            self.ComboMagma.SetValue('')
            for res in self.result_dict[magID]:
                self.ComboMagma.Append(res)
            self.ComboMagma.SetSelection(0)
        else:
            self.ComboMagma.SetValue('')
            for node in list(self.pdb.keys()):
                self.ComboMagma.Append(node)
            self.ComboMagma.SetSelection(0)


    def on_N_button(self, event):
        self.ax_reset=1
        #self.ax_resetH=1
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)

        if(self.MAGMA):
            self.UpdateMagmaBox()
        self.selection=[]
        self.draw_figure()


    def on_P_button2(self, event):
        self.ax_reset2=1
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()-1)
        self.selection=[]
        self.draw_figure()

    def on_N_button2(self, event):
        self.ax_reset2=1
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.draw_figure()


    def on_select_button(self, event):
        print('Click on panel 1 to select')
        self.selection=[]
        self.SELECT=1

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
            print('removing:',self.slice_service.connections[sele].tag)
            self.slice_service.connections.pop(sele)

        self.selection=[]
        self.draw_figure()




    def on_search_button(self, event):
        print('Click in panel 1 to initiate search')
        self.SEARCH=1



    def on_save_button(self, event):
        self.outfile=self.textbox_savelist.GetValue()
        print()
        print('Saving ',len(self.slice_service.connections),'entries in connectivity table to ',self.outfile)
        outy=open(self.outfile,'w')
        for i in range(len(self.slice_service.connections)):
            outy.write('%s\t%s\n' % (self.index_data[self.slice_service.connections[i][0]],self.index_data[self.slice_service.connections[i][1]]))
        outy.close()


    #when search button is pressed make selection
    def on_pick(self, event):
        #print event.xdata,event.ydata
        if(self.SELECT==1):
            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            argy=[]
            cbs1=self.ComboBox1.GetSelection()
            print(cbs1)
            for i,cn in enumerate(self.slice_service.connections):
                #print self.slice_service.peaks[p].name
                if(cn.v1==cbs1):
                    xval=cn.f4  #proton
                    yval=cn.f3  #carbon
                    #print cn.v1,cn.v2,event.xdata,event.ydata,cn.f1,cn.f2,cn.f3,cn.f4
                    rad2=((xval-event.xdata)/xdist)**2.+((yval-event.ydata)/ydist)**2.
                    raddy.append(rad2)
                    argy.append(i)
            raddy=numpy.array(raddy)
            maxy=argy[numpy.argmin(raddy)]
            self.selection.append(maxy)

            cn=self.slice_service.connections[maxy]
            print('selected:',cn.tag)

            for i,pk in enumerate(self.slice_service.peaks):
                if(cn.p2==pk.name):
                    self.ComboBox2.SetSelection(i)

            for i,cn2 in enumerate(self.slice_service.connections):
                if(cn2.p1==cn.p2):
                    if(cn2.p2==cn.p1):
                        print('Reciprocated:',cn2.tag)

                        self.selection.append(i)

            self.SELECT=0
            print(self.selection,self.SELECT)
            self.draw_figure()
        if(self.SEARCH==1):
            self.selection=[]
            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            raddy=[]
            argy=[]
            cbs1=self.ComboBox1.GetSelection()
            print(cbs1)
            poss=[]

            self.noise=float(self.textbox0.GetValue())

            for i,pk in enumerate(self.slice_service.peaks):
                    xval=pk.ppmI  #proton
                    yval=pk.ppmJ  #carbon
                    #print cn.v1,cn.v2,event.xdata,event.ydata,cn.f1,cn.f2,cn.f3,cn.f4
                    rad2=((xval-event.xdata)/xdist)**2.+((yval-event.ydata)/ydist)**2.
                    raddy.append(rad2)
                    argy.append(i)
                    if(rad2<5E-6): #if selection is less than threshold...
                        print('possible:',pk.name,rad2)
                        poss.append(i)
            if(len(poss)==0):
                raddy=numpy.array(raddy)
                maxy=argy[numpy.argmin(raddy)]
                poss.append(maxy)


           #of the options, test if its reciprocated and return intensities
            scr=[]
            loc=[]
            for i,po in enumerate(poss):
                add,i1,i2=self.NOETest(cbs1,po) #test possible NOE assignment
                if(add==1):
                    scr.append((i1+i2))
                    loc.append((po,i1,i2))

            scrmax=numpy.argmax(scr) #get the one that has the biggest score
            self.AddNOE(cbs1,loc[scrmax][0],loc[scrmax][1],loc[scrmax][2]) #add the one with the biggest combined intensities

            pk2=self.slice_service.peaks[ loc[scrmax][0] ]
            #set second combobox to identity of cross peak
            for i,pk in enumerate(self.slice_service.peaks):
                if(pk2.name==pk.name):
                    self.ComboBox2.SetSelection(i)

            self.draw_figure()
            self.SEARCH=0




            """
            print 'Testing for NOE:'
            add=0



            pk1=self.slice_service.peaks[self.ComboBox1.GetSelection()]
            pk2=self.slice_service.peaks[maxy]
            if(pk1.name==pk2.name):
                print 'this is a diagonal peak. Not adding'
                add=1

            print pk1.name,pk2.name

            #set second combobox to identity of cross peak
            for i,pk in enumerate(self.slice_service.peaks):
                if(pk2.name==pk.name):
                    self.ComboBox2.SetSelection(i)

            #is this already here?
            for i,cn in enumerate(self.slice_service.connections):
                if(cn.p1==pk1.name and cn.p2==pk2.name):
                    print 'NOE already present:',cn.tag
                    add=1
                if(cn.p2==pk1.name and cn.p1==pk2.name):
                    print 'NOE (reciprocated) already present:',cn.tag
                    add=1
            i1=self.slice_service.data[pk1.indexI,pk1.indexJ,pk2.indexK,pk2.indexL]
            i2=self.slice_service.data[pk2.indexI,pk2.indexJ,pk1.indexK,pk1.indexL]

            noise=self.parent.tabOne.dmax*float(self.parent.tabOne.threshBox.GetValue())
            print 'noise:',noise
            print 'Intensity of selected cross peak:'
            print i1,i1/noise
            print 'Intensity of reciprocated cross peak:'
            print i2,i2/noise
            if(i1<noise):
                print 'I1 less than noise. Probably not a cross peak.'
                add=1
            if(i2<noise):
                print 'I2 less than noise. Probably not a cross peak.'
                add=1
            print 'Ratio (expected to be >50%):'
            print min(i1,i2)/max(i1,i2)*100.


            if(add==0):
                print 'Adding NOE:'
                #from deconFrame import connEntry


                stry=('%s\t%s\t%f\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f\t%f\t%f') % (pk1.name,pk2.name,pk1.ppmI,pk1.ppmJ,pk2.ppmK,pk2.ppmL,i1,i2,i1/noise,i2/noise,0,0,0,0,0,0)
                test=stry.split()
                cnNew=self.parent.tabOne.connEntry(test,sym='y',peak=self.slice_service.peaks,dim=4)
                self.slice_service.connections.append(cnNew)
                self.selection.append(len(self.slice_service.connections)-1)

                stry=('%s\t%s\t%f\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f\t%f\t%f') % (pk2.name,pk1.name,pk2.ppmI,pk2.ppmJ,pk1.ppmK,pk1.ppmL,i2,i1,i2/noise,i1/noise,0,0,0,0,0,0)
                test=stry.split()
                cnNew=self.parent.tabOne.connEntry(test,sym='y',peak=self.slice_service.peaks,dim=4)
                self.slice_service.connections.append(cnNew)
                self.selection.append(len(self.slice_service.connections)-1)
            else:
                print 'Not adding.'

            #self.selection.append()
            self.draw_figure()
            self.SEARCH=0
            """



    def on_noise_button(self, event):
        self.spec=self.readfile('out/slice2D/'+self.ComboBox1.GetValue()+'.dat.out')
        #self.nxmin=float(self.textbox_xmin.GetValue())
        #self.nxmax=float(self.textbox_xmax.GetValue())
        #self.nymin=float(self.textbox_ymin.GetValue())
        #self.nymax=float(self.textbox_ymax.GetValue())
        #self.noise_max=float(getnoise_spec(self.spec,self.nxmin,self.nxmax,self.nymin,self.nymax))
        #self.textbox0.SetValue(str(self.noise_max*float(self.textbox_fac.GetValue())))
        self.textbox0.SetValue(str(self.slice_service.max_intensity*float(self.textbox_fac.GetValue())))
        self.draw_figure()


    def on_AutoFit_button(self, event):
        #RUN THE AUTOFIT ALGORITHM
        #self.nxmin=float(self.textbox_xmin.GetValue())
        #self.nxmax=float(self.textbox_xmax.GetValue())
        #self.nymin=float(self.textbox_ymin.GetValue())
        #self.nymax=float(self.textbox_ymax.GetValue())
        #self.noisefac=float(self.textbox_fac.GetValue())
        #self.slice_service.connections=analslices3d(self.slice_service.peaks,self.noisefac,self.nxmin,self.nxmax,self.nymin,self.nymax)
        #print self.plotty
        self.cb_grid_auto.SetValue(1)
        self.draw_figure()




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
