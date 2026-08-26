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




#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

class peakAssFrame(wx.Frame):


    def __init__(self,parent):

        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "strip plots",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          size=(self.monitorWidth*0.5, self.monitorHeight*0.5),
                          #size=(1300,670)
                          )
        panel = wx.Panel(self)
        self.SetBackgroundColour('WHITE')

        self.Nbutton = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)
        self.Pbutton = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)
        self.closebutton = wx.Button(self, -1,"Close")
        self.Bind(wx.EVT_BUTTON, self.on_close_button, self.closebutton)



        self.parent=parent
        self.molecule=parent.molecule
        self.peak=parent.molecule.peak

        self.listy=[]
        for i in range(len(self.molecule.spec['hnco'].peak2D)):
            self.listy.append(self.molecule.spec['hnco'].peak2D[i].name)

        for pk in self.listy:
            ERRORS=self.GetErrors(pk)
            if(len(ERRORS)==0):
                print(pk,'Looking good!')
            else:
                print(pk)
                for ERROR in ERRORS:
                    print('   ',ERROR)


        self.create_main_panel()
        print('22')
        #self.draw_figure()
        print('44')
        self.Show(True)
        #self.Fit()


    def create_main_panel(self):
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.listy, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.on_draw_button, self.ComboBox1) 

        self.index={}
        for i,li in enumerate(self.listy):
            self.index[li]=i

        self.shiftLab=wx.StaticText(self,label="")

        self.ErrorLab=[]
        for i in range(4):
            self.ErrorLab.append(wx.StaticText(self,label=""))
        
        self.residuesLab=wx.StaticText(self,label="Residue Options: ")

        self.source = AutoWidthListCtrl(self)
        self.source.SetMinSize((1000,400))
        self.source.InsertColumn(0, 'spectrum', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(1, 'name', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(2, 'f1', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(3, 'f2', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(4, 'f3', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(5, 'Atom', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(6, 'Label', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(7, 'S/N', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(8, 'Assignment', width = 100,format=wx.LIST_FORMAT_CENTRE) 
        self.source.InsertColumn(9, 'AssignmentOptions', width = 100,format=wx.LIST_FORMAT_CENTRE) 

        self.hbox=wx.BoxSizer(wx.VERTICAL)

        self.hbox2=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.ComboBox1)
        self.hbox2.Add(self.Nbutton)
        self.hbox2.Add(self.Pbutton)
        self.hbox2.Add(self.closebutton)

        self.hbox.Add(self.hbox2)
        
        self.hbox3=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.residuesLab)
        self.hbox3.Add(self.shiftLab)
        self.hbox.Add(self.hbox3)

        self.hbox.Add(self.source)
        
        for E in self.ErrorLab:
            self.hbox.Add(E)

        self.SetSizerAndFit(self.hbox)

        self.PopulateList()
    
    def on_draw_button(self,event):
        self.PopulateList()

    def on_N_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.PopulateList()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.PopulateList()

    def on_close_button(self, event):
        self.Close()


    def GetErrors(self,pk):
        ERRORS=[]
        specord='hnco','hncaco','hnca','hncoca','hncacb','hncocacb'
        for spec in specord:
            if(spec in self.peak[pk].keys()):
                if(spec=='hnco'):
                    if(len(self.peak[pk][spec])>1):
                        ERRORS.append("too many peaks in the HNCO (1 expected) %s %s" % pk,spec)
                if(spec=='hncaco'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCACO (2 expected) %s %s" % pk,spec)
                if(spec=='hnca'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCA (max 2 expected) %s %s" % pk,spec)
                if(spec=='hncoca'):
                    if(len(self.peak[pk][spec])>1):
                        ERRORS.append("too many peaks in the HNCOCA (max 1 expected) %s %s" % pk,spec)
                if(spec=='hncacb'):
                    if(len(self.peak[pk][spec])>4):
                        ERRORS.append("too many peaks in the HNCOCA (max 4 expected) %s %s" % pk,spec)
                    nmax=0;nmin=0;pmax=0;pmin=0;
                    for i,pk2 in enumerate(self.peak[pk][spec]):
                        if(pk2.tp=='PosMax'):
                            pmax+=1
                        if(pk2.tp=='PosMin'):
                            pmin+=1
                        if(pk2.tp=='NegMax'):
                            nmax+=1
                        if(pk2.tp=='NegMin'):
                            nmin+=1
                    if(pmin>1 or pmax>1 or nmax>1 or nmin>1):
                        ERRORS.append("MISS-ASSIGNED HNCACB: pmax %i pmin %i nmax %i nmin %i %s %s" %(pmax,pmin,nmax,nmin,pk,spec))

                if(spec=='hncocacb'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCOCA (max 2 expected) %s %s" % pk,spec)
        self.GetShift(pk)
        for key,vals in self.shufty.items():
            if(len(vals)>1):
                vols=[] #extract shifts
                for val in vals:
                    vols.append(val[0])
                if(numpy.std(vols)>1):
                    ERRORS.append('ERROR: inconsistent shifts: %s %.2f %.2f %s' % (key,numpy.average(vols),numpy.std(vols),vols))
                    
                    diffs=[]
                    for i in range(len(vols)):
                        diff=0
                        for j in range(len(vols)):
                            if(i!=j):
                                diff+=vols[i]-vols[j]
                        diffs.append(diff)
                    argy=numpy.argmax(diffs)
                    spoc=vals[argy][1]
                    ii=vals[argy][2]
                    pkErr=self.peak[pk][spoc][ii]
                    ERRORS.append('%s %s %.2f %s %.2f' % (spoc,pkErr.name,pkErr.f3,pkErr.tp,pkErr.inty))

        return ERRORS
                        
    def GetShift(self,pk):
        self.shufty={}
        for spec in self.peak[pk].keys():
            for i,pk3 in enumerate(self.peak[pk][spec]):
                lab=self.GetLab(spec,pk3.tp)
                if(lab not in self.shufty.keys()):
                    self.shufty[lab]=[]
                self.shufty[lab].append((pk3.f3,spec,i))
            

    def GetLab(self,spec,tp):
        lab=''
        if(spec=='hnco'):
            lab='CO(i-1)'
        if(spec=='hncaco'):
            if(tp=="main"):
                lab='CO(i)'                    
            else:
                lab='CO(i-1)'                    
        if(spec=='hnca'):
            if(tp=="main"):
                lab='CA(i)'                    
            else:
                lab='CA(i+1)'   
        if(spec=='hncoca'):
            lab='CA(i+1)'                    
                 
        if(spec=='hncacb'):
            if(tp=="PosMax"):
                lab='CA(i)'                    
            elif(tp=="NegMax"):
                lab='CB(i)'                    
            elif(tp=="PosMin"):
                lab='CA(i+1)'                    
            elif(tp=="NegMin"):
                lab='CB(i+1)'                    
        if(spec=='hncocacb'):
            if(tp=="NegMin"):
                lab='CB(i+1)'                    
            else:
                lab='CA(i+1)'                    
        return lab

    def PopulateList(self):
        sele=self.ComboBox1.GetSelection()
        pk=self.listy[sele]
        cnt=1
        self.source
        stry=''
        for res in self.parent.molecule.shift[pk]:
            stry+=res+' '
        self.shiftLab.SetLabel(stry)
        print('uijk')


        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

                                  

        specord='hnco','hncaco','hnca','hncoca','hncocacb'
        for spec in specord:
            if(spec in self.peak[pk].keys()):
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    #print(pk,spec,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.tp,pk3.inty)
                    num_items = self.source.GetItemCount()
                    self.source.InsertItem(num_items,str(cnt))
                    self.source.SetItem(num_items,0,str(spec))
                    self.source.SetItem(num_items,1,str(pk3.name))
                    self.source.SetItem(num_items,2,'%.2f' % pk3.f1)            
                    self.source.SetItem(num_items,3,'%.2f' % pk3.f2)            
                    self.source.SetItem(num_items,4,'%.2f' % pk3.f3)            

                    lab=self.GetLab(spec,pk3.tp)

                    self.source.SetItem(num_items,5,lab)            
                    self.source.SetItem(num_items,6,pk3.tp)            
                    self.source.SetItem(num_items,7,'%.2f' % (pk3.inty/self.parent.molecule.spec[spec].noise))            
                    if(spec=='hnco'):
                        #if(pk3.tp=="main"):
                            fstr=''
                            for edge in self.parent.molecule.G1edges[pk]:
                                if(edge[2]=='f'):
                                    val='%.2f%s' % (edge[1],edge[2])
                                    fstr+=edge[0]+'('+val+') '
                                    self.source.SetItem(num_items,8,fstr)

                            fstr=''
                            for edge in self.parent.molecule.G1edgesFull[pk]:
                                if(edge not in self.parent.molecule.G1edges[pk]):
                                    if(edge[2]=='f'):
                                        val='%.2f%s' % (edge[1],edge[2])
                                        fstr+=edge[0]+'('+val+') '
                                        self.source.SetItem(num_items,9,fstr)
                    if(spec=='hnca'):
                        if(pk3.tp!="main"):
                            fstr=''
                            for edge in self.parent.molecule.G1edges[pk]:
                                if(edge[2]=='b'):
                                    val='%.2f%s' % (edge[1],edge[2])
                                    fstr+=edge[0]+'('+val+') '
                                    self.source.SetItem(num_items,9,fstr)
                            self.source.SetItem(num_items,8,fstr)

                            fstr=''
                            for edge in self.parent.molecule.G1edgesFull[pk]:
                                if(edge not in self.parent.molecule.G1edges[pk]):
                                    if(edge[2]=='b'):
                                        val='%.2f%s' % (edge[1],edge[2])
                                        fstr+=edge[0]+'('+val+') '
                                        self.source.SetItem(num_items,9,fstr)



                    print()
                    #self.source.SetItem(num_items,3,'%.2f' % (val[1]))    
                    cnt+=1

        ERRORS=self.GetErrors(pk)
        for i,lab in enumerate(self.ErrorLab):
            lab.SetLabel("")
        if(len(ERRORS)==0):
            self.ErrorLab[0].SetLabel("No obvious errors in peak assignment")
        else:
            for i,ERROR in enumerate(ERRORS):
                self.ErrorLab[i].SetLabel(ERROR)


        pass
