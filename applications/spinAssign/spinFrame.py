#!/usr/bin/python

import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin
import matplotlib.pyplot as plt
from . import textEdit

import assign.reviewFrame


#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

class spinFrame(wx.Frame):


    def __init__(self,parent):

        #print(parent.molecule.G1edges)
        #.molecule.G1edges[pk]
        #sys.exit(10)
        print('xzy')
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "strip plots",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          #size=(self.monitorWidth*0.5, self.monitorHeight*0.5),
                          size=(1300,670)
                          )
        panel = wx.Panel(self)
        self.SetBackgroundColour('WHITE')

        self.parent=parent

        

        self.listy=parent.parent.tabTwo.listy



        self.index={}
        for i,li in enumerate(self.listy):
            self.index[li]=i

        self.create_main_panel()
        self.draw_figure()
        self.Show(True)
        #self.Fit()


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.

        self.checkForwardList = self.parent.molecule.checkForward
        self.checkBackwardList = self.parent.molecule.checkBackward


        # Check for errors for each peak so that can then flag peaks that have errors in the combo box 
        self.parent.molecule.CollateErrorsForComboBox()
        peak_errors = self.parent.molecule.errors_for_combobox
        potential_im1_proline = self.parent.molecule.potential_proline_m1
        potential_glycine = self.parent.molecule.potential_glycines



        # Create a list of 0/1 for each peak depending on if no error/errors found
        peak_errors_list = []
        for i, pk in enumerate(self.listy):
            try:
               peak_errors_list.append(peak_errors[pk])
            except:
               peak_errors_list.append(0)


        # Create a lift of 0/1 for each peak depending on if the i-1 peak is potentially a proline
        potential_im1_proline_list = []
        for i, pk in enumerate(self.listy):
            try:
                potential_im1_proline_list.append(potential_im1_proline[pk])
            except:
                potential_im1_proline_list.append(0)


        potential_glycine_list = []
        for i, pk in enumerate(self.listy):
            try:
                potential_glycine_list.append(potential_glycine[pk])
            except:
                potential_glycine_list.append(0)



        # Create a new list of peak labels, but highlight ones that require attention 
        # with a '*f' or '*b' if missing a forward or backward connection respectively,
        # or with '(error)' if there are any errors
        self.list2 = []

        

        for i, peakLabel in enumerate(self.listy):
            peakLabel_new = peakLabel
            if(peakLabel in self.checkForwardList):
                peakLabel_new = peakLabel_new + '*f'
            if(peakLabel in self.checkBackwardList):
                peakLabel_new = peakLabel_new + '*b'
            if(peak_errors_list[i]==1):
                peakLabel_new = peakLabel_new + '(error)'
            if(potential_im1_proline_list[i]==1):
                peakLabel_new = peakLabel_new + '(i-1=P)'


            self.list2.append(peakLabel_new)



        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(800,250))
        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)

        self.toolbar = NavigationToolbar(self.canvas)

        self.plotLbl = wx.StaticBox(self, -1, 'LocalSpinSystem:')
        self.plotSizer = wx.StaticBoxSizer(self.plotLbl, wx.VERTICAL)
        self.plotSizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND| wx.GROW | wx.ALL )
        self.plotSizer.Add(self.toolbar, 0, wx.EXPAND)


        self.ComboBoxLab = wx.StaticText(self, label="PeakID:")
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(150, -1), choices=self.list2, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox1)

        self.residuesLab=wx.StaticText(self,label="Residue (i): ")
        self.shiftLab=wx.StaticText(self,label="")
        self.ErrorLab=[]
        for i in range(3):
            self.ErrorLab.append(wx.StaticText(self,label=""))

        self.residues2Lab=wx.StaticText(self,label="Residue (i-1): ")
        self.shift2Lab=wx.StaticText(self,label="")



        lblList = ['Network','Shifts','Shiftx2']
        self.rbox = wx.RadioBox(self,label = '',choices = lblList ,majorDimension = 0, style = wx.RA_SPECIFY_COLS)
        self.rbox.Bind(wx.EVT_RADIOBUTTON,self.redraw)

        self.Nbutton = wx.Button(self, -1,"Next",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton = wx.Button(self, -1,"Prev",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        self.loadbutton = wx.Button(self, -1,"Load",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.onGetFile, self.loadbutton)


        self.savebutton = wx.Button(self, -1,"Save",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.OnSaveResults, self.savebutton)


        self.closebutton = wx.Button(self, -1,"Close",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.on_close_button, self.closebutton)

        self.graphbutton = wx.Button(self, -1,"Graph",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.on_graph_button, self.graphbutton)

        self.errorbutton = wx.Button(self, -1,"Error",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.onButtonError, self.errorbutton)

        self.reviewbutton = wx.Button(self, -1,"Review",size=(50,-1))
        self.reviewbutton.bind(self.OnButtonReview)



        
        self.source = AutoWidthListCtrl(self)
        #self.source = wx.ListCtrl(self,style=wx.LC_REPORT,size=(650,300))
        self.source.SetMinSize((800,200))
        self.source.InsertColumn(0, 'Spectrum', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(1, 'Name', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(2, 'F1', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(3, 'F2', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(4, 'F3', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(5, 'Atom', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(6, 'Label', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(7, 'S/N', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(8, 'Assignment', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(9, 'PossibleCorrelations', width = 1200)

        """
        self.source.SetMinSize((600,200))

        self.source.InsertColumn(0, 'Type', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(1, 'Path', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(2, 'Inputfile', width = 200,format=wx.LIST_FORMAT_CENTRE)

        self.source.InsertColumn(3, 'Peaklist',width=200,format=wx.LIST_FORMAT_CENTRE)
        """

        self.PopulateList()




        self.vboxFull=wx.BoxSizer(wx.VERTICAL)
        # Layout with box sizers


        #self.hboxMain=wx.BoxSizer(wx.VERTICAL)
        #self.hboxMain.Add(self.source,wx.ALL | wx.EXPAND)

        self.hbox=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.ComboBoxLab)
        self.hbox.Add(self.ComboBox1)
        self.hbox.Add(self.Nbutton)
        self.hbox.Add(self.Pbutton)
        self.hbox.Add(self.closebutton)
        self.hbox.Add(self.loadbutton)
        self.hbox.Add(self.savebutton)

        self.hbox.Add(self.graphbutton)
        self.hbox.Add(self.errorbutton)
        self.hbox.Add(self.reviewbutton)
        self.hbox.Add(self.rbox)

        self.hbox3=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.residuesLab)
        self.hbox3.Add(self.shiftLab)

        self.hbox4=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox4.Add(self.residues2Lab)
        self.hbox4.Add(self.shift2Lab)




        self.vboxFull.Add(self.hbox)
        self.vboxFull.Add(self.hbox3)
        self.vboxFull.Add(self.hbox4)
        self.vboxFull.Add(self.source,wx.ALL)
        for e in self.ErrorLab:
            self.vboxFull.Add(e)
        self.vboxFull.Add(self.plotSizer)


        #self.draw_figure()
        #self.canvas.draw()
        self.SetSizerAndFit(self.vboxFull)

    def redraw(self,event):
        self.draw_figure()

    def on_pick(self,event):
        if(self.rbox.GetSelection==1):
            return
        #print(event.xdata,event.ydata)
        if(event.xdata>0):
            test=[]
            loc=[]
            for val in self.places['i+1']:
                test.append(numpy.fabs(val[1][1]-event.ydata))
                loc.append(val[0])
            argy=numpy.argmin(test)

            self.parent.parent.tabTwo.ComboBox1.SetSelection(self.index[loc[argy]])
            self.parent.parent.tabTwo.draw_figure()

        else:
            test=[]
            loc=[]
            for val in self.places['i-1']:
                test.append(numpy.fabs(val[1][1]-event.ydata))
                loc.append(val[0])
            #print(loc)
            #print(test)
            argy=numpy.argmin(test)
            #print(loc[argy],test[argy])
            self.parent.parent.tabTwo.ComboBox3.SetSelection(self.index[loc[argy]])
            self.parent.parent.tabTwo.draw_figure()


    def OnSaveResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            defaultFile='peaks',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            outfile=path+'.peak'
            self.parent.molecule.save_peaks(outfile)



            #self.canvas.print_figure(path, dpi=self.dpi)
            #self.parent.parent.flash_status_message("Saved %s" % path)


    def onGetFile(self, e):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
                            wildcard="Peak file (*.peak)|*.peak|", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            infile = dlg.GetPath()
            #print('Reading', infile)
            self.parent.molecule.load_peaks(infile)
        dlg.Destroy()
        self.PopulateList()
        self.draw_figure()


    def PopulateList(self):
        sele=self.ComboBox1.GetSelection()
        pk=self.listy[sele]
        cnt=1
        self.source
        stry=''


        for res in self.parent.molecule.shift[pk]:
            stry+=res[0]+' '
        self.shiftLab.SetLabel(stry)

        stry=''
        for res in self.parent.molecule.shift2[pk]:
            stry+=res[0]+' '
        self.shift2Lab.SetLabel(stry)



        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        #self.parent.molecule.listy()

        print('Keys for %s: %s' %(pk,self.parent.molecule.peak[pk].keys()))

        specord='hnco','hncaco','hnca','hncoca', 'hncacb', 'hncocacb', 'hncanh', 'hncocanh','ctocsy','hcconh'

        
        for spec in specord:
            if(spec in self.parent.molecule.peak[pk].keys()):
                for i,pk3 in enumerate(self.parent.molecule.peak[pk][spec]):
                    #print(pk,spec,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.tp,pk3.inty)
                    #print('asf')
                    num_items = self.source.GetItemCount()
                    self.source.InsertItem(num_items,str(cnt))
                    self.source.SetItem(num_items,0,str(spec))
                    self.source.SetItem(num_items,1,str(pk3.name))
                    self.source.SetItem(num_items,2,'%.2f' % pk3.f1)
                    self.source.SetItem(num_items,3,'%.2f' % pk3.f2)

                    self.source.SetItem(num_items,4,'%.2f' % pk3.f3p)

                    lab=self.parent.molecule.GetLab(spec,pk3.tp)

                    self.source.SetItem(num_items,5,lab)
                    self.source.SetItem(num_items,6,pk3.tp)
                    self.source.SetItem(num_items,7,'%.2f' % (pk3.inty/self.parent.molecule.spec[spec].noise))

                    #write in possible assignment options
                    if(spec=='hncaco'):
                        if(pk3.tp=="main"):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'COf','f'))

                    if(spec=='hnco'):
                        self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                        self.source.SetItem(num_items,9,self.GetAssOptions(pk,'COb','b'))

                    if(spec=='hnca'):
                        if(pk3.tp!="main"):
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CAb','b'))
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                        else:
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CAf','f'))
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))

                    if(spec=='hncacb'):
                        #if(pk3.tp=='PosMin'):
                        #    self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                        #    self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbA','b'))
                        if(pk3.tp=='NegMin'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CBb','b'))
                        #elif(pk3.tp=='PosMax'):
                        #    self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                        #    self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbC','f'))
                        elif(pk3.tp=='NegMax'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CBf','f'))


                    if(spec=='hncocacb'):
                        if(pk3.tp=='Neg'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CBf','f'))
                        elif(pk3.tp=='Pos'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'CBb','b'))
                    
                    if(spec=='hncanh'):
                        if(pk3.tp=='diag'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'NHb1','b'))
                            # self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            # self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncanhD','b'))

                        if(pk3.tp=='minus'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'NHb2','b'))
                            # self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            # self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncanhC','b'))
                    
                        


                    
                    if(spec=='hncocanh'):
                        if(pk3.tp=='plus'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'NHf1','f'))
                        if(pk3.tp=='diag'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'NHf2','f'))




                    # if(spec=='ctocsy'):
                    #     print(pk3)
                    #     print('wjc')
                    #     self.source.SetItem(num_items,10,self.GetAss(pk,'b'))
                    #     self.source.SetItem(num_items,11,self.GetAssOptions(pk,'hncacbA','b'))

                    cnt+=1
        self.errors_for_ComboBox = {}
        ERRORS=self.parent.molecule.GetErrors(pk)
        
        for i,lab in enumerate(self.ErrorLab):
            lab.SetLabel("")
        if(len(ERRORS)==0):
            self.ErrorLab[0].SetLabel("No obvious errors in peak assignment")
        else:
            for i,ERROR in enumerate(ERRORS):
                try:
                    self.ErrorLab[i].SetLabel(ERROR)
                except:
                    pass


        return

        #OLD SOURCE
        cnt=1
        for i,pk in enumerate(self.listy):
            for val in self.parent.molecule.G1edges[pk]: #ed
                resi=self.parent.molecule.resi[i]
                resn=self.parent.molecule.seq[resi]
                num_items = self.source.GetItemCount()
                self.source.InsertItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,str(pk))
                self.source.SetItem(num_items,1,str(val[0]))
                self.source.SetItem(num_items,2,str(val[2]))
                self.source.SetItem(num_items,3,'%.2f' % (val[1]))



            cnt+=1

    def GetAss(self,pk,tp):
        fstr=''
        #print(self.parent.molecule.G1edges)
        try:
            for edge in self.parent.molecule.G1edges[pk]:
                if(edge[2]==tp):
                    val='%.2f%s' % (edge[1],edge[2])
                    fstr+=edge[0]+'('+val+') '
        except:
            pass

        return fstr


    def GetAssOptions(self,pk,spec,tp):
        fstr=''
        if(pk in self.parent.molecule.Optedges.keys()):
            if(spec in self.parent.molecule.Optedges[pk].keys()):
                for edge in self.parent.molecule.Optedges[pk][spec]:
                    #if(edge not in self.parent.molecule.G1edges[pk] and edge[2]==tp):
                    #if(edge[2]==tp):
                        val='%.2f%s' % (edge[1],edge[2])
                        fstr+=edge[0]+'('+val+') '
                        # if(pk == "3H-N"):
                        #     if(spec=='hncocanh'):
                        #         print(fstr)
                        #         exit(100)
        return fstr

    def draw_figure(self):
        self.fig.clear()
        self.sele1=self.ComboBox1.GetSelection()
        self.places={}



        self.axes5 = self.fig.add_subplot(111)
        if(self.rbox.GetSelection()==0):
            self.plotLbl.SetLabel('LocalSpinSystem:')
            self.axes5.set_axis_off()

            cenMain=numpy.array((0,0))
            radMain=1.
            xspace=4

            self.places['i']=[]
            self.places['i'].append((self.listy[self.sele1],cenMain))
            self.AddCircle(self.axes5,cenMain,'r',self.listy[self.sele1],radMain)

            self.axes5.text(cenMain[0],cenMain[1]+radMain+0.5,'i',fontsize=10,color='r',horizontalalignment='center',verticalalignment='center')
            self.axes5.set_xlim(-10,10)
            self.axes5.set_ylim(-6,1)
            #getForwards and score


            cntF=0;cntB=0
            if(self.listy[self.sele1] in self.parent.molecule.G1edges.keys()):
                for val in self.parent.molecule.G1edges[self.listy[self.sele1]]: #edges for selected res
                    if(val[2]=='f'): #for forward NOES..
                        cenMid=cenMain+(xspace,cntF) #step one right
                        if('i+1' not in self.places.keys()):
                            self.places['i+1']=[]
                        self.places['i+1'].append((val[0],cenMid))
                        rad=1.
                        self.AddCircle(self.axes5,cenMid,'b',val[0],rad)
                        if(cntF==0):
                            self.axes5.text(cenMid[0],cenMid[1]+radMain+0.5,'i+1',fontsize=10,color='r',horizontalalignment='center',verticalalignment='center')

                        a=cenMain+(rad,0)    #right of main circle
                        b=cenMid-(rad,0) #left of new circle
                        self.AddArrow(self.axes5,a,b,val[1],'b')



                        cntN=cntF #look at the edges of this guy...
                        for vol in self.parent.molecule.G1edges[val[0]]:
                            if(vol[0]!=self.listy[self.sele1]):

                                cenEnd=cenMid+(xspace,cntN)
                                if(cntN==0):
                                    self.axes5.text(cenEnd[0],cenEnd[1]+radMain+0.5,'i+2',fontsize=10,color='r',horizontalalignment='center',verticalalignment='center')
                                self.AddCircle(self.axes5,cenEnd,'c',vol[0],0.5)
                                a=cenMid+(1,0) #right of middle circle
                                b=cenEnd-(1,0) #left of end circle
                                self.AddArrow(self.axes5,a,b,vol[1],'b')
                                cntN-=1
                            else:
                                if(vol[2]=='b'):
                                    a=cenMid-(1,0.5) #from left of middle circle
                                    b=cenMain+(1,-0.5) #to right of main
                                    self.AddArrow(self.axes5,a,b,vol[1],'g')
                        cntF-=2
                    if(val[2]=='b'):

                        cenMid=cenMain+(-xspace,cntB) #step one right
                        if('i-1' not in self.places.keys()):
                            self.places['i-1']=[]
                        self.places['i-1'].append((val[0],cenMid))
                        rad=1.
                        self.AddCircle(self.axes5,cenMid,'g',val[0],1)
                        if(cntB==0):
                            self.axes5.text(cenMid[0],cenMid[1]+radMain+0.5,'i-1',fontsize=10,color='r',horizontalalignment='center',verticalalignment='center')

                        a=cenMain-(1,0) #left of centre
                        b=cenMid+(1,0) #right of middle

                        self.AddArrow(self.axes5,a,b,val[1],'g')
                        cntN=0 #look at the edges of this guy...
                        for vol in self.parent.molecule.G1edges[val[0]]:
                            if(vol[0]!=self.listy[self.sele1]):
                                cenEnd=cenMid-(xspace,-cntN)
                                self.AddCircle(self.axes5,cenEnd,'c',vol[0],0.5)
                                if(cntN==0):
                                    self.axes5.text(cenEnd[0],cenEnd[1]+radMain+0.5,'i-2',fontsize=10,color='r',horizontalalignment='center',verticalalignment='center')
                                a=cenEnd+(0.5,0)
                                b=cenMid-(1,0)
                                self.AddArrow(self.axes5,a,b,val[1],'g')
                                cntN-=2
                            else:
                                if(vol[2]=='f'):
                                    a=cenMid+(1,-0.5)
                                    b=cenMain+(-1,-0.5)
                                    self.AddArrow(self.axes5,a,b,val[1],'b')

                                    #self.axes5.arrow(-3,-0.5,2,0, color='b',length_includes_head=True,head_length=0.5,head_width=0.5)
                                    #self.axes5.text(-2,-0.5,'%.2f' % vol[1],fontsize=10)

                        cntB-=2

            if(self.listy[self.sele1] in self.parent.molecule.G1edgesFull.keys()):
                for val in self.parent.molecule.G1edgesFull[self.listy[self.sele1]]: #look for rejected edges
                    if(val not in self.parent.molecule.G1edges[self.listy[self.sele1]]):
                        if(val[2]=='f'): #for forward NOES..
                            cenMid=cenMain+(xspace,cntF) #step one right

                            if('i+1' not in self.places.keys()):
                                self.places['i+1']=[]
                            self.places['i+1'].append((val[0],cenMid))

                            rad=0.5
                            self.AddCircle(self.axes5,cenMid,'c',val[0],rad)

                            a=numpy.array((cenMain[0],cenMid[1]))+(1,0)
                            b=cenMid+(-0.5,0)
                            self.AddArrow(self.axes5,a,b,val[1],'c')

                            cntF-=1

                            for vol in self.parent.molecule.G1edgesFull[val[0]]: #look for reciprocity
                            #print(vol[0],self.listy[self.sele1])
                                if(vol[0]==self.listy[self.sele1] and vol[2]=='b'):
                                    b=numpy.array((cenMain[0],cenMid[1]))+(1,0)
                                    a=cenMid+(-0.5,0)
                                    self.AddArrow(self.axes5,a,b,vol[1],'c',down='y')

                        if(val[2]=='b'): #for forward NOES..
                            cenMid=cenMain-(xspace,-cntB) #step one right

                            if('i-1' not in self.places.keys()):
                                self.places['i-1']=[]
                            self.places['i-1'].append((val[0],cenMid))
                            rad=0.5
                            self.AddCircle(self.axes5,cenMid,'c',val[0],rad)

                            a=numpy.array((cenMain[0],cenMid[1]))-(1,0)
                            b=cenMid+(+0.5,0)
                            self.AddArrow(self.axes5,a,b,val[1],'c')

                            for vol in self.parent.molecule.G1edgesFull[val[0]]: #look for reciprocity
                            #print(vol[0],self.listy[self.sele1])
                                if(vol[0]==self.listy[self.sele1] and vol[2]=='f'):
                                    b=numpy.array((cenMain[0],cenMid[1]))-(1,0)
                                    a=cenMid+(+0.5,0)
                                    self.AddArrow(self.axes5,a,b,vol[1],'c',down='y')
                            cntB-=1

        elif(self.rbox.GetSelection()==1): #draw chemical shifts
            self.plotLbl.SetLabel('Residue probabilities:')
            resns,prob,prob2=self.parent.molecule.AnalPeak(self.listy[self.sele1])

            pos=numpy.arange(len(prob))
            self.axes5.bar(pos-0.2,prob,color='orange',width=0.4)
            self.axes5.bar(pos+0.2,prob2,color='cyan',width=0.4)


            pos3=[]
            prob3=[]
            for i in range(len(pos)):
                tig=0
                for ras in self.parent.molecule.shift[self.listy[self.sele1]]:
                    if(resns[i]==ras[0]):
                        tig=1
                        break
                        #if(resns[i] in self.parent.molecule.shift[self.listy[self.sele1]]):
                if(tig==1):
                    pos3.append(i)
                    prob3.append(prob[i])
            pos3=numpy.array(pos3)
            self.axes5.bar(pos3-0.2,prob3,color='r',width=0.4)
            pos4=[]
            prob4=[]
            for i in range(len(pos)):
                tig=0
                for ras in self.parent.molecule.shift2[self.listy[self.sele1]]:
                    if(resns[i]==ras[0]):
                        tig=1
                        break
                if(tig==1):
                    #if(resns[i] in self.parent.molecule.shift2[self.listy[self.sele1]]):
                    pos4.append(i)
                    prob4.append(prob2[i])
            pos4=numpy.array(pos4)
            self.axes5.bar(pos4+0.2,prob4,color='b',width=0.4)



            self.axes5.set_xticks(pos)
            self.axes5.set_xticklabels(resns)
            self.axes5.set_ylabel("Probability",fontsize=8)


            yl=self.parent.molecule.tolMax,self.parent.molecule.tolMax
            xl=0,len(pos)
            self.axes5.plot(xl,yl)
            yl=self.parent.molecule.tolMin,self.parent.molecule.tolMin
            xl=0,len(pos)
            self.axes5.plot(xl,yl)

            Gender=['maxLim','minLim','i(ex)','i-1(ex)','i','i-1']
            self.axes5.legend(Gender,loc=2)

            cen=(0,0)
            self.places['i']=[]
            self.places['i'].append((self.listy[self.sele1],cen))

            for val in self.parent.molecule.G1edges[self.listy[self.sele1]]: #edges for selected res
                if(val[2]=='f'): #for forward NOES..
                    if('i+1' not in self.places.keys()):
                        self.places['i+1']=[]
                    self.places['i+1'].append((val[0],cen))
                if(val[2]=='b'):
                    if('i-1' not in self.places.keys()):
                        self.places['i-1']=[]
                    self.places['i-1'].append((val[0],cen))




        elif(self.rbox.GetSelection()==2): #draw chemical shifts
            resns,probs=self.parent.molecule.CompareShiftx2(self.listy[self.sele1])
            pos=numpy.arange(len(probs))

            self.axes5.bar(pos,probs,color='orange',width=1)

            self.axes5.set_xticks(pos)
            self.axes5.set_xticklabels(resns)
            self.axes5.set_ylabel("Probability",fontsize=8)

        #print('shlad')
        #plt.savefig('TheFig.png')

        self.canvas.draw()

        #update parent
        if('i' in self.places.keys()):
            self.parent.parent.tabTwo.ComboBox2.SetSelection(self.index[self.places['i'][0][0]])
        if('i+1' in self.places.keys()):
            if(len(self.places['i+1'])>0):
                self.parent.parent.tabTwo.ComboBox1.SetSelection(self.index[self.places['i+1'][0][0]])
        if('i-1' in self.places.keys()):
            if(len(self.places['i-1'])>0):
                self.parent.parent.tabTwo.ComboBox3.SetSelection(self.index[self.places['i-1'][0][0]])
        # self.parent.parent.tabTwo.draw_figure()




        return


    def AddArrow(self,axes,a,b,lab,col,down='n'):
        axes.arrow(a[0],a[1],b[0]-a[0],b[1]-a[1], color=col,length_includes_head=True,head_length=0.5,head_width=0.5)
        if(down=='y'):
            axes.text((a[0]+b[0])*0.5,(b[1]+a[1])*0.5-0.2,'%.2f' % lab,fontsize=8,color=col,horizontalalignment='center',verticalalignment='center')
        else:
            axes.text((a[0]+b[0])*0.5,(b[1]+a[1])*0.5+0.2,'%.2f' % lab,fontsize=8,color=col,horizontalalignment='center',verticalalignment='center')

    def AddCircle(self,axes,cen,col,lab,rad):
        circle1 = plt.Circle(cen, rad, color=col)
        axes.add_artist(circle1)
        axes.text(cen[0],cen[1],lab,fontsize=10,horizontalalignment='center',verticalalignment='center')

    def on_draw_button(self, event):
        self.PopulateList()
        self.draw_figure()

    def onButtonError(self,event):
        ERRORS=[]
        for pk in self.parent.molecule.peak.keys():
            entry=self.parent.molecule.GetErrors(pk)
            if(len(entry)!=0):
                ERRORS.append(('',pk,))
                ERRORS.append(entry)
        # print(ERRORS)
        textEdit.MyFrame(ERRORS,stream='y')

    def OnButtonReview(self,event):
        import reviewFrame
        assFrame=reload(reviewFrame)
        bool=reviewFrame.reviewFrameMan(self)

        
    def on_N_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.selection=[]
        self.PopulateList()
        self.draw_figure()

    def on_P_button(self, event):
        self.ax_reset1=1
        self.ax_reset0=1
        self.inc=0
        #if(self.cb_flip.GetValue()):
        #    self.ax_reset2=1
        #if(self.cb_decon.GetValue()):
        #    self.ax_reset2=1

        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.selection=[]
        self.PopulateList()
        self.draw_figure()

    def on_close_button(self, event):
        self.Close()

    def on_graph_button(self, event):
        self.parent.molecule.normSpec()
        #self.parent.molecule.assSpec()
        self.parent.molecule.AssCACB()
        #self.parent.molecule.EdgeScreen()
        #self.parent.molecule.WriteInit()
