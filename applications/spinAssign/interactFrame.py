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
import assign.textEdit

from SettingsUnidec import Parse


from magma.magma import Magma

#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

class interactFrame(wx.Frame):


    def __init__(self,parent):

        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "Interactive Assignment",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          size=(self.monitorWidth*0.5, self.monitorHeight*0.5),
                          #size=(1300,670)
                          )
        panel = wx.Panel(self)
        self.SetBackgroundColour('WHITE')

        self.parent=parent

        self.listy=parent.parent.tabTwo.listy
        self.G1_nodes=parent.parent.tabOne.molecule.G1_nodes

        self.inst=Magma(self.parent.magmaParFile,run='n') #get instance of magma


        self.SELECT=0
        self.G1=[]
        self.G2=[]
        self.results={}
        self.assSeq={}
        self.place1=[]
        self.node1=[]
        self.place2=[]
        self.node2=[]
        self.PICK=1

        self.index={}
        for i,li in enumerate(self.listy):
            self.index[li]=i
        self.indexSeq={}
        for i,li in enumerate(self.G1_nodes):
            self.indexSeq[li]=i

        print(self.indexSeq)
        self.create_main_panel()
        print('1')
        self.draw_figure()
        print('2')
        self.Show(True)
        print('3')

    def onGuessButton(self,event):
        self.Reccomend()

    def Reccomend(self):

        nodeList=[]
        noeAdj=[]
        for i,node in enumerate(self.inst.noe_node_list):
            if node not in list(self.results.keys()):
                nodeList.append(node)
                row=[]
                for adj in self.inst.noe_adjacency[i]:
                    if(adj not in self.results.keys()):
                        row.append(adj)
                noeAdj.append(row)

        connected_subgraphs = self.inst.G.SplitConnectionsOverSubgraphs(nodeList,noeAdj)
        subgraphRef={}
        print('Breaking data graph into disconnected subgraphs')
        print('la')
        for i,subgraph in enumerate(connected_subgraphs):
            try:

                sub_conn_dict = self.inst.G.GetConnectionDict(subgraph)

                noe_node_list,noe_adjacency= self.inst.G.GetNodesAdjacency(sub_conn_dict)
                subgraphRef[i]={}
                subgraphRef[i]['subgraph']=subgraph
                subgraphRef[i]['nodes']=noe_node_list
                subgraphRef[i]['adj']=noe_adjacency
                subgraphRef[i]['noes']=self.inst.CountNOEs(noe_adjacency)
                subgraphRef[i]['noesScr']=self.inst.CountNOEs(noe_adjacency,weight=self.inst.P.weightEdgesG1)
                print('DisconnectedDataGraph:  ',i+1,'Nodes: ',len(subgraphRef[i]['nodes']),' Restraints: ',subgraphRef[i]['noes'],' MaxScr: ',subgraphRef[i]['noesScr'])
            except:
                pass


        iscr=[]
        ikey=[]
       
        for key,vals in subgraphRef.items(): #for each subgraph...
            sub_noe_node_list=vals['nodes']
            ncnt=len(sub_noe_node_list)
            #if(ncnt>nmax):
            #    nmax=ncnt
            #    ikey=key
            iscr.append(ncnt)
            ikey.append(int(key))


        #sub_noe_node_list=subgraphRef[ikey]['nodes']
        #sub_noe_adjacency=subgraphRef[ikey]['adj']

        scrs=[]
        nodes=[]
        seqs=[]
        for key,vals in subgraphRef.items(): #for each subgraph...
            node,seq,scr=self.DoMove(vals['nodes'],vals['adj'])

            #if(len(scr)==1):
            scrs.append(scr)
            nodes.append(node)
            seqs.append(seq)
            print(key,':',len(vals['nodes']),node,seq,scr)
            print
        for i,node in enumerate(nodes):
            print(node,seqs[i],scrs[i])

        print(nodes)
        argy=numpy.argmax(scrs)
        node=nodes[argy]
        seq=seqs[argy]
        scr=scrs[argy]
        print('reccomended:',node,seq,'gives:',scr)

        if(node==-1):
            print('No good options!')
        else:
            self.ComboBox1.SetSelection(self.index[node])
            self.ComboBox2.SetSelection(self.indexSeq[seq])
        self.onCombo(True)

    def DoMove(self,nodeList,noeAdj):

        node,nodeVal=self.WalkBack(nodeList,noeAdj)

        while(1==1):
            seq,scr,probs=self.GetSuggestion(node)
            #if(len(seq)!=1):
            #    break
            print('GGreccomended:',node,seq,'gives:',scr)
            if((scr<3 and len(nodeList)>3) or scr==1):
                print('am in the if statement')
                node,nodeVal=self.StepForward(node,nodeList,noeAdj)
                if(node==-1):
                    print('break 1')
                    break
            else:
                print('break 2')
                break
        print('GQreccomended:',node,seq,'gives:',scr)
        return node,seq,scr


    def GetSuggestion(self,node):
        print('node = %s' % node)
        print('candidates = %s' % self.parent.parent.tabOne.molecule.candidates)
        opts=self.parent.parent.tabOne.molecule.candidates[node]
        self.resultsTmp=copy.deepcopy(self.results)
        self.assSeqTmp=copy.deepcopy(self.assSeq)
        extras=[]
        assigns=[]
        opty=[]
        print('opts = %s' % opts)
        for opt in opts:
            print('opt = %s' % opt)
            print('extras for each opt= %s' % extras)
            if(opt not in self.assSeq.keys()):
                print('extras at start of if statement = %s' % extras)
                self.results=copy.deepcopy(self.resultsTmp)
                self.assSeq=copy.deepcopy(self.assSeqTmp)
                assign=self.DoAssign(node,opt)
                extras.append(len(assign))
                assigns.append(assign)
                opty.append(opt)
                print('extras end of if statement = %s' % extras)
        self.results=copy.deepcopy(self.resultsTmp)
        self.assSeq=copy.deepcopy(self.assSeqTmp)
        #self.assigns=numpy.array(assigns)
        self.assigns=assigns
        print('extras after = %s' % extras)
        opty=numpy.array(opty)
        argy=numpy.argmax(extras)
        extras=numpy.array(extras)
        mask=(extras==extras[argy])
        if(len(extras[mask])==1):
            return opty[mask][0],extras[mask][0],1
        print('found',len(opty[mask]),'solutions of size',extras[argy])
        probs=[]
        for i in range(len(assigns)):
            if(extras[i]==extras[argy]):
                #print('testing',opts[i],extras[i])
                prob=1
                for ass in assigns[i]:
                    pk=ass[0]
                    se=ass[1]
                    #print('examining ',pk,se)
                    #print(self.parent.parent.tabOne.molecule.assRef.keys())
                    for vols in self.parent.parent.tabOne.molecule.assRef[se]:
                        #print(vols)
                        if(vols[0]==pk):
                            prob*=vols[1]
                            break
                probs.append(prob)
                #print(prob)
        argy=numpy.argmax(probs)
        #print('result:',opts[mask][argy],extras[mask][argy],probs[argy])
        return opty[mask][argy],extras[mask][argy],probs[argy]

        """
        if(len(extras[mask])>1):
            print('WARNING: multiple ways of getting this assignment')
            opts=numpy.array(opts)
            print(opts[mask])

        return opts[argy],extras[argy]
        """


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.


        self.fig=Figure()
        self.fig2=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(400,250))
        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        self.toolbar = NavigationToolbar(self.canvas)

        self.canvas2 = FigCanvas(self, -1, self.fig2)
        self.canvas2.SetMinSize(wx.Size(400,250))
        # Bind the 'pick' event for selection
        self.canvas2.mpl_connect('button_press_event', self.on_pick2)
        self.toolbar2 = NavigationToolbar(self.canvas2)

        self.canvasSizer=wx.BoxSizer(wx.HORIZONTAL)

        self.plotLbl = wx.StaticBox(self, -1, 'SpinSystems:')
        self.plotSizer = wx.StaticBoxSizer(self.plotLbl, wx.VERTICAL)
        self.plotSizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND| wx.GROW | wx.ALL )
        self.plotSizer.Add(self.toolbar, 0, wx.EXPAND)

        self.plotLbl2 = wx.StaticBox(self, -1, 'Sequence:')
        self.plotSizer2 = wx.StaticBoxSizer(self.plotLbl2, wx.VERTICAL)
        self.plotSizer2.Add(self.canvas2, 1, wx.LEFT | wx.TOP | wx.EXPAND| wx.GROW | wx.ALL )
        self.plotSizer2.Add(self.toolbar2, 0, wx.EXPAND)


        self.canvasSizer.Add(self.plotSizer,1,wx.GROW)
        self.canvasSizer.Add(self.plotSizer2,1,wx.GROW)

        self.ComboBoxLab = wx.StaticText(self, label="Peak:")
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.listy, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.onCombo, self.ComboBox1)

        self.ComboBoxLab2 = wx.StaticText(self, label="Seq:")
        self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=self.G1_nodes, style=wx.CB_READONLY)
        self.ComboBox2.SetSelection(0)
        self.Bind(wx.EVT_COMBOBOX, self.onCombo, self.ComboBox2)

        self.residuesLab=wx.StaticText(self,label="Residue (i): ")
        self.shiftLab=wx.StaticText(self,label="")
        self.ErrorLab=[]
        for i in range(3):
            self.ErrorLab.append(wx.StaticText(self,label=""))

        self.residues2Lab=wx.StaticText(self,label="Residue (i-1): ")
        self.shift2Lab=wx.StaticText(self,label="")



        lblList = ['Network','Shifts','Shiftx2','Mask']
        self.rbox = wx.RadioBox(self,label = '',choices = lblList ,majorDimension = 0, style = wx.RA_SPECIFY_COLS)
        self.rbox.SetSelection(3)

        self.rbox.Bind(wx.EVT_RADIOBUTTON,self.redraw)

        self.Assbutton = wx.Button(self, -1,"Assign",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.onAssButton, self.Assbutton)

        #self.Nbutton = wx.Button(self, -1,"Next",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        #self.Pbutton = wx.Button(self, -1,"Prev",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        #self.loadbutton = wx.Button(self, -1,"Load",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.onGetFile, self.loadbutton)


        self.savebutton = wx.Button(self, -1,"Save",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.OnSaveResults, self.savebutton)


        self.closebutton = wx.Button(self, -1,"Close",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.on_close_button, self.closebutton)

        self.clearbutton = wx.Button(self, -1,"Clear",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.onClearButton, self.clearbutton)

        self.Guessbutton = wx.Button(self, -1,"Guess",size=(50,-1))
        self.Bind(wx.EVT_BUTTON, self.onGuessButton, self.Guessbutton)



        #self.removebutton = wx.Button(self, -1,"Remove",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_graph_button, self.graphbutton)

        #self.graphbutton = wx.Button(self, -1,"Graph",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.on_graph_button, self.graphbutton)

        #self.errorbutton = wx.Button(self, -1,"Error",size=(50,-1))
        #self.Bind(wx.EVT_BUTTON, self.onButtonError, self.errorbutton)

        self.source = AutoWidthListCtrl(self)
        #self.source = wx.ListCtrl(self,style=wx.LC_REPORT,size=(650,300))
        self.source.SetMinSize((200,200))
        self.source.InsertColumn(0, 'Peak', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(1, 'Seq', width = 150,format=wx.LIST_FORMAT_LEFT)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick, self.source)

        self.source2 = AutoWidthListCtrl(self)
        self.source2.SetMinSize((200,200))
        self.source2.InsertColumn(0, 'Seq', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.source2.InsertColumn(1, 'Peak', width = 150,format=wx.LIST_FORMAT_LEFT)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick2, self.source2)

        self.source3 = AutoWidthListCtrl(self)
        self.source3.SetMinSize((200,200))
        self.source3.InsertColumn(0, 'Peak', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source3.InsertColumn(1, 'Seq', width = 100,format=wx.LIST_FORMAT_LEFT)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick3, self.source3)

        self.source4 = AutoWidthListCtrl(self)
        self.source4.SetMinSize((200,200))
        self.source4.InsertColumn(0, 'Peak', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source4.InsertColumn(1, 'Seq', width = 100,format=wx.LIST_FORMAT_LEFT)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick4, self.source4)

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

        self.hbox.Add(self.ComboBoxLab2)
        self.hbox.Add(self.ComboBox2)
        #self.hbox.Add(self.Nbutton)
        #self.hbox.Add(self.Pbutton)
        self.hbox.Add(self.Assbutton)
        self.hbox.Add(self.closebutton)
        self.hbox.Add(self.clearbutton)
        self.hbox.Add(self.Guessbutton)

        self.hbox.Add(self.savebutton)

        #self.hbox.Add(self.graphbutton)
        #self.hbox.Add(self.errorbutton)
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

        self.sources=wx.BoxSizer(wx.HORIZONTAL)
        self.sources.Add(self.source)
        self.sources.Add(self.source2)
        self.sources.Add(self.source3)
        self.sources.Add(self.source4)

        self.vboxFull.Add(self.sources,wx.ALL)
        for e in self.ErrorLab:
            self.vboxFull.Add(e)
        self.vboxFull.Add(self.canvasSizer)


        #self.draw_figure()
        #self.canvas.draw()
        self.SetSizerAndFit(self.vboxFull)



    def redraw(self,event):
        self.draw_figure()

    def on_pick(self,event):
        if(self.rbox.GetSelection()!=0):
            return
        xx=event.xdata
        yy=event.ydata
        self.place1=numpy.array(self.place1)
        vals=(xx-self.place1[:,0])**2.+(yy-self.place1[:,1])**2.
        argy=numpy.argmin(vals)
        print(argy)
        print(self.node1[argy])
        self.ComboBox1.SetSelection(self.index[self.node1[argy]])

        if(self.node1[argy] not in self.results.keys()):
            self.PICK=1
            count = self.source.GetItemCount()
            vals=[self.source.GetItem(row, 0).GetText() for row in range(count)]
            for i,val in enumerate(vals):
                if(self.source.IsSelected(i)):
                    self.source.Select(i, on=False)
            for i,val in enumerate(vals):
                if(self.node1[argy]==val and self.source.IsSelected(i)==False):
                    self.source.Focus(i)
                    self.source.Select(i,on=True)
        else:
            self.PICK=2
            count = self.source4.GetItemCount()
            vals=[self.source4.GetItem(row, 0).GetText() for row in range(count)]
            for i,val in enumerate(vals):
                if(self.source4.IsSelected(i)):
                    self.source4.Select(i, on=False)
            for i,val in enumerate(vals):
                if(self.node1[argy]==val and self.source4.IsSelected(i)==False):
                    self.source4.Focus(i)
                    self.source4.Select(i,on=True)

        self.OnDoubleClick(True)

    def on_pick2(self,event):
        if(self.rbox.GetSelection()!=0):
            return
        xx=event.xdata
        yy=event.ydata
        self.place2=numpy.array(self.place2)
        vals=(xx-self.place2[:,0])**2.+(yy-self.place2[:,1])**2.
        argy=numpy.argmin(vals)
        self.ComboBox2.SetSelection(self.indexSeq[self.node2[argy]])


        if(self.node2[argy] not in self.assSeq.keys()):
            self.PICK=3
            count = self.source2.GetItemCount()
            vals=[self.source2.GetItem(row, 0).GetText() for row in range(count)]
            for i,val in enumerate(vals):
                if(self.source2.IsSelected(i)):
                    self.source2.Select(i, on=False)
            for i,val in enumerate(vals):
                if(self.node2[argy]==val and self.source2.IsSelected(i)==False):
                    self.source2.Focus(i)
                    self.source2.Select(i,on=True)
        else:
            self.PICK=4
            count = self.source4.GetItemCount()
            vals=[self.source4.GetItem(row, 1).GetText() for row in range(count)]
            for i,val in enumerate(vals):
                if(self.source4.IsSelected(i)):
                    self.source4.Select(i, on=False)
            for i,val in enumerate(vals):
                if(self.node2[argy]==val and self.source4.IsSelected(i)==False):
                    self.source4.Focus(i)
                    self.source4.Select(i,on=True)

        self.OnDoubleClick2(True)


    def OnSaveResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save Results...",
            defaultDir=os.path.join(os.getcwd(),'dat'),
            defaultFile='Assigned_interactive.list',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()


            outy=open(path,'w')
            for key,val in self.results.items():
                outy.write('%s : ' % key)
                outy.write('%s ' % val)
                outy.write('\n')
            for key,vals in self.parent.parent.tabOne.molecule.candidates.items():
                if(key not in self.results.keys()):
                    outy.write('%s : ' % key)
                    for val in vals:
                        outy.write('%s ' % val)
                    outy.write('\n')
            outy.close()

            # for key, val in 


            #self.canvas.print_figure(path, dpi=self.dpi)
            #self.parent.parent.flash_status_message("Saved %s" % path)


    def onGetFile(self, e):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
                            wildcard="Peak file (*.peak)|*.peak|", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            infile = dlg.GetPath()
            print('Reading', infile)
            self.parent.parent.tabOne.molecule.load_peaks(infile)
        dlg.Destroy()
        self.PopulateList()
        self.draw_figure()


    def PopulateList(self):

        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)
        num_items = self.source2.GetItemCount()
        for i in range(num_items):
            self.source2.DeleteItem(0)

        cnt=1
        col1=[]
        col2=[]
        for key in self.parent.parent.tabOne.molecule.G1edges.keys():
            col1.append(key)
            cc=0
            for opt in self.parent.parent.tabOne.molecule.candidates[key]: #residues
                if(opt not in self.assSeq.keys()):
                    cc+=1
            col2.append(cc) #get indices

        args=numpy.argsort(col2)
        cnt=0
        for arg in args:
            key=col1[arg]
            #print(self.parent.parent.tabOne.molecule.candidates[key])

            if(key not in self.results.keys()):
                num_items = self.source.GetItemCount()
                self.source.InsertItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,key)

                stry=''
                for opt in self.parent.parent.tabOne.molecule.candidates[key]:
                    if(opt not in self.assSeq.keys()):
                        stry+=opt+' '
                self.source.SetItem(num_items,1,stry)
                cnt+=1


        cnt=0
        for shiftxresi in self.parent.parent.tabOne.molecule.shiftxresis:
            lab=str(shiftxresi)+self.parent.parent.tabOne.molecule.shiftx2[shiftxresi]['resn']
            if(lab not in self.assSeq.keys() and lab[-1]!='P' and int(lab[:-1])!=self.parent.parent.tabOne.molecule.FirstResidue):
                stry=''
                if(lab in self.parent.parent.tabOne.molecule.assRef.keys()):
                    for lob in self.parent.parent.tabOne.molecule.assRef[lab]:
                        if(lob[0] not in self.results.keys()):
                            stry+=lob[0]+' '

                num_items = self.source2.GetItemCount()
                self.source2.InsertItem(num_items,str(cnt))
                self.source2.SetItem(num_items,0,lab)
                self.source2.SetItem(num_items,1,stry)


        pk=self.listy[self.ComboBox1.GetSelection()]
        stry=''
        
        for res in self.parent.parent.tabOne.molecule.shift[pk]:
            stry+=res[0]+' '
        self.shiftLab.SetLabel(stry)

        stry=''
        for res in self.parent.parent.tabOne.molecule.shift2[pk]:
            stry+=res[0]+' '
        self.shift2Lab.SetLabel(stry)


        """
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        #self.parent.parent.tabOne.molecule.listy()

        specord='hnco','hncaco','hnca','hncoca','hncacb','hncocacb'
        for spec in specord:
            if(spec in self.parent.parent.tabOne.molecule.peak[pk].keys()):
                for i,pk3 in enumerate(self.parent.parent.tabOne.molecule.peak[pk][spec]):
                    #print(pk,spec,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.tp,pk3.inty)
                    num_items = self.source.GetItemCount()
                    self.source.InsertItem(num_items,str(cnt))
                    self.source.SetItem(num_items,0,str(spec))
                    self.source.SetItem(num_items,1,str(pk3.name))
                    self.source.SetItem(num_items,2,'%.2f' % pk3.f1)
                    self.source.SetItem(num_items,3,'%.2f' % pk3.f2)

                    self.source.SetItem(num_items,4,'%.2f' % pk3.f3p)

                    lab=self.parent.parent.tabOne.molecule.GetLab(spec,pk3.tp)

                    self.source.SetItem(num_items,5,lab)
                    self.source.SetItem(num_items,6,pk3.tp)
                    self.source.SetItem(num_items,7,'%.2f' % (pk3.inty/self.parent.parent.tabOne.molecule.spec[spec].noise))

                    #write in possible assignment options
                    if(spec=='hncaco'):
                        if(pk3.tp=="main"):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncaco','f'))

                    if(spec=='hnco'):
                        self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                        self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hnco','b'))

                    if(spec=='hnca'):
                        if(pk3.tp!="main"):
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncoca','b'))
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                        else:
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hnca','f'))
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))

                    if(spec=='hncacb'):
                        if(pk3.tp=='PosMin'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbA','b'))
                        elif(pk3.tp=='NegMin'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'b'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbB','b'))
                        elif(pk3.tp=='PosMax'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbC','f'))
                        elif(pk3.tp=='NegMax'):
                            self.source.SetItem(num_items,8,self.GetAss(pk,'f'))
                            self.source.SetItem(num_items,9,self.GetAssOptions(pk,'hncacbD','f'))

                    cnt+=1

        ERRORS=self.parent.parent.tabOne.molecule.GetErrors(pk)
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

        """
        return

        #OLD SOURCE
        cnt=1
        for i,pk in enumerate(self.listy):
            for val in self.parent.parent.tabOne.molecule.G1edges[pk]: #ed
                resi=self.parent.parent.tabOne.molecule.resi[i]
                resn=self.parent.parent.tabOne.molecule.seq[resi]
                num_items = self.source.GetItemCount()
                self.source.InsertItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,str(pk))
                self.source.SetItem(num_items,1,str(val[0]))
                self.source.SetItem(num_items,2,str(val[2]))
                self.source.SetItem(num_items,3,'%.2f' % (val[1]))



            cnt+=1

    def GetAss(self,pk,tp):
        fstr=''
        try:
            for edge in self.parent.parent.tabOne.molecule.G1edges[pk]:
                if(edge[2]==tp):
                    val='%.2f%s' % (edge[1],edge[2])
                    fstr+=edge[0]+'('+val+') '
        except:
            pass

        return fstr


    def GetAssOptions(self,pk,spec,tp):
        fstr=''
        if(pk in self.parent.parent.tabOne.molecule.Optedges.keys()):
            if(spec in self.parent.parent.tabOne.molecule.Optedges[pk].keys()):
                for edge in self.parent.parent.tabOne.molecule.Optedges[pk][spec]:
                    #if(edge not in self.parent.parent.tabOne.molecule.G1edges[pk] and edge[2]==tp):
                    #if(edge[2]==tp):
                        val='%.2f%s' % (edge[1],edge[2])
                        fstr+=edge[0]+'('+val+') '
        return fstr

    def GetNodeNumber(self,test,nodeList):
        for j,node in enumerate(nodeList):
            if(node==test):
                return j
        print('could not find node')
        return -1


    def WalkBack(self,sub_noe_node_list,sub_noe_adjacency,Ass=False):
        node=sub_noe_node_list[0] #take this as first guess
        nodeVal=0                 #take this as first guess
        cnt=0
        while(1==1):
            cnt+=1
            tick=0
            print(cnt,node,nodeVal)
            if(cnt==50):
                sys.exit(100)
            for adj in sub_noe_adjacency[nodeVal]:
                if(self.inst.NMR.noes[node][adj][2]=='b'):
                    node=adj
                    nodeVal=self.GetNodeNumber(node,sub_noe_node_list)
                    tick=1
                    break
            if(tick==0):
                break
        return node,nodeVal

    def StepForward(self,node,sub_noe_node_list,sub_noe_adjacency):
        nodeVal=self.GetNodeNumber(node,sub_noe_node_list)
        for adj in sub_noe_adjacency[nodeVal]:
            if(self.inst.NMR.noes[node][adj][2]=='f'):
                node=adj
                nodeVal=self.GetNodeNumber(node,sub_noe_node_list)
                return node,nodeVal
        return -1,-1

    def onCombo(self,event):
        self.PopulateResults()
        self.PopulateList()
        self.draw_figure()

    def draw_figure(self):
        self.fig.clear()
        self.fig2.clear()

        self.sele1=self.ComboBox1.GetSelection()
        self.places={}

        self.axes5 = self.fig.add_subplot(111)
        if(self.rbox.GetSelection()==0): # NETWORK PLOT
            #self.plotLbl.SetLabel('LocalSpinSystem:')
            self.axes5.set_axis_off()

            yrun=0
            xMax=0

            self.place1=[]
            self.node1=[]
            self.node2=[]
            self.place2=[]

            cnt=0
            for key,vals in self.inst.subgraphRef.items(): #for each subgraph...
                sub_noe_node_list=vals['nodes']
                sub_noe_adjacency=vals['adj']


                # print('a')
                node,nodeVal=self.WalkBack(sub_noe_node_list,sub_noe_adjacency)

                print('interact start node:',node)
                cenMain=numpy.array((0,cnt))

                #wire up the circles.
                jobs={}
                jobs[0]=[]
                jobs[0].append((node,nodeVal,0))
                ylevel=0
                xlevel={}
                xlevel[0]=0
                radMain=1
                done=0
                bigCNT=0
                while(1==1):
                    bigCNT+=1
                    if(bigCNT==1000):
                        break
                    #print(jobs, ylevel)
                    node,nodeVal,xpos=jobs[ylevel][0] #unpack current place in the queue
                    #print('drawing',node,'at',xpos,ylevel,'(',nodeVal,')')
                    if(xpos>xMax):
                        xMax=xpos
                    jobs[ylevel].pop(0) #remove entry zero from front of queue

                    xlevel[ylevel]=xpos

                    cenMain=(xpos,ylevel+yrun)
                    col='r'
                    lob=''
                    if(self.SELECT):
                        if(node in self.G1):
                            col='g'
                            lob=node

                    if(node in self.results.keys()):
                        col='y'
                    self.AddCircle(self.axes5,cenMain,col,lob,radMain) #draw circle and label
                    self.place1.append(cenMain)
                    self.node1.append(node)


                    cnt=0
                    for adj in sub_noe_adjacency[nodeVal]:
                        #print(node,adj)
                        if(self.inst.NMR.noes[node][adj][2]=='f'):
                            #node=adj
                            for j in range(len(sub_noe_node_list)):
                                if(adj==sub_noe_node_list[j]):
                                    newVal=j
                                    break
                            if(cnt not in jobs.keys()):
                                jobs[cnt]=[]

                            if(cnt==ylevel):
                                jobs[cnt].append((adj,newVal,xlevel[ylevel]+2))
                            else:
                                jobs[cnt].append((adj,newVal,xlevel[ylevel]))
                            cnt+=1

                    xlevel[ylevel]+=2
                    jump=0
                    while(1==1):
                        jump+=1

                        if(len(jobs.keys())==0): #if no jobs, we're done.
                            done=1
                            break

                        if(ylevel in list(jobs.keys())): #is current level in jobs?
                            if(len(jobs[ylevel])==0):
                                del jobs[ylevel]

                        if(len(jobs.keys())==0): #if no jobs, we're done.
                            done=1
                            break
                        else:
                            ylevel=list(numpy.min(jobs.keys()))[0]
                        if(jump==100):
                            break

                    #print(jump,done)
                    if(done==1):
                        break

                yrun+=4

            self.axes5.set_xlim(-2,xMax+2)
            self.axes5.set_ylim(-2,yrun+2)


            self.axes6 = self.fig2.add_subplot(111)
            ylevel=0
            xpos=0

            self.axes6.set_axis_off()

            xMax=0
            pcnt=0
            for resi in self.parent.parent.tabOne.molecule.resi:

                resn=self.parent.parent.tabOne.molecule.seq[resi]
                if(resn=='P'):
                    if(pcnt==0):
                        ylevel+=4
                        xpos=0
                    pcnt=1
                else:
                    pcnt=0
                    lab=str(resi)+resn
                    cenMain=(xpos,ylevel)

                    col='b'
                    lob=''
                    if(self.SELECT):
                        #print('tg',lab,self.G2,lab in self.G2)
                        if(lab in self.G2):
                            col='g'
                            lob=lab

                    #print(col,lob)
                    if(lab in self.assSeq.keys()):
                        col='y'
                    self.AddCircle(self.axes6,cenMain,col,lob,radMain) #draw circle and label
                    self.place2.append(cenMain)
                    self.node2.append(lab)

                    xpos+=2
                    if(xpos>xMax):
                        xMax=xpos

            self.axes6.set_xlim(-2,xMax+2)
            self.axes6.set_ylim(-2,ylevel+2)

            """
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
            if(self.listy[self.sele1] in self.parent.parent.tabOne.molecule.G1edges.keys()):
                for val in self.parent.parent.tabOne.molecule.G1edges[self.listy[self.sele1]]: #edges for selected res
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
                        for vol in self.parent.parent.tabOne.molecule.G1edges[val[0]]:
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
                        for vol in self.parent.parent.tabOne.molecule.G1edges[val[0]]:
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

            if(self.listy[self.sele1] in self.parent.parent.tabOne.molecule.G1edgesFull.keys()):
                for val in self.parent.parent.tabOne.molecule.G1edgesFull[self.listy[self.sele1]]: #look for rejected edges
                    if(val not in self.parent.parent.tabOne.molecule.G1edges[self.listy[self.sele1]]):
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

                            for vol in self.parent.parent.tabOne.molecule.G1edgesFull[val[0]]: #look for reciprocity
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

                            for vol in self.parent.parent.tabOne.molecule.G1edgesFull[val[0]]: #look for reciprocity
                            #print(vol[0],self.listy[self.sele1])
                                if(vol[0]==self.listy[self.sele1] and vol[2]=='f'):
                                    b=numpy.array((cenMain[0],cenMid[1]))-(1,0)
                                    a=cenMid+(+0.5,0)
                                    self.AddArrow(self.axes5,a,b,vol[1],'c',down='y')
                            cntB-=1
            """
        elif(self.rbox.GetSelection()==1): #draw chemical shifts
            self.plotLbl.SetLabel('Residue probabilities:')
            resns,prob,prob2=self.parent.parent.tabOne.molecule.AnalPeak(self.listy[self.sele1])
            pos=numpy.arange(len(prob))

            self.axes5.bar(pos-0.2,prob,color='orange',width=0.4)
            self.axes5.bar(pos+0.2,prob2,color='cyan',width=0.4)


            pos3=[]
            prob3=[]
            for i in range(len(pos)):
                for vol in self.parent.parent.tabOne.molecule.shift[self.listy[self.sele1]]:
                    if(vol[0]==resns[i]):
                        pos3.append(i)
                        prob3.append(prob[i])
            pos3=numpy.array(pos3)
            self.axes5.bar(pos3-0.2,prob3,color='r',width=0.4)
            pos4=[]
            prob4=[]
            for i in range(len(pos)):
                for vol in self.parent.parent.tabOne.molecule.shift2[self.listy[self.sele1]]:
                    if(vol[0]==resns[i]):
                        pos4.append(i)
                        prob4.append(prob2[i])
            pos4=numpy.array(pos4)
            self.axes5.bar(pos4+0.2,prob4,color='b',width=0.4)



            self.axes5.set_xticks(pos)
            self.axes5.set_xticklabels(resns)
            self.axes5.set_ylabel("Probability",fontsize=8)


            yl=self.parent.parent.tabOne.molecule.tolMax,self.parent.parent.tabOne.molecule.tolMax
            xl=0,len(pos)
            self.axes5.plot(xl,yl)
            yl=self.parent.parent.tabOne.molecule.tolMin,self.parent.parent.tabOne.molecule.tolMin
            xl=0,len(pos)
            self.axes5.plot(xl,yl)

            Gender=['maxLim','minLim','i(ex)','i-1(ex)','i','i-1']
            self.axes5.legend(Gender,loc=2)

            cen=(0,0)
            self.places['i']=[]
            self.places['i'].append((self.listy[self.sele1],cen))

            for val in self.parent.parent.tabOne.molecule.G1edges[self.listy[self.sele1]]: #edges for selected res
                if(val[2]=='f'): #for forward NOES..
                    if('i+1' not in self.places.keys()):
                        self.places['i+1']=[]
                    self.places['i+1'].append((val[0],cen))
                if(val[2]=='b'):
                    if('i-1' not in self.places.keys()):
                        self.places['i-1']=[]
                    self.places['i-1'].append((val[0],cen))

        elif(self.rbox.GetSelection()==2): #draw chemical shifts
            resns,probs=self.parent.parent.tabOne.molecule.CompareShiftx2(self.listy[self.sele1])

            newprobs=[]
            newresns=[]
            for i in range(len(resns)):
                if(resns[i] not in self.assSeq.keys()):
                    newresns.append(resns[i])
                    newprobs.append(probs[i])
            probs=newprobs
            resns=newresns
            pos=numpy.arange(len(probs))

            self.axes5.bar(pos,probs,color='orange',width=1)

            self.axes5.set_xticks(pos)
            self.axes5.set_xticklabels(resns)
            self.axes5.set_ylabel("Probability",fontsize=8)

        elif(self.rbox.GetSelection()==3):
            #first, invert the assRef dictionary.
            mol=self.parent.parent.tabOne.molecule
           
            def RestoreBadNOEs(bad,noes):
                for pk,vals in noes.items():
                    if(pk in bad):
                        #print()
                        #print (pk,vals)
                        #print (pk,mol.G1edges[pk])
                        noes[pk]=mol.G1edges[pk]
            def CleanNOEs(assGreat,noes):
                for pk,at in assGreat.items():
                    #print(pk,at)
                    #print(assOpt[at])
                    #print(noes[pk])

                    for pok,noelist in noes.items():
                        noeListNew=[]
                        for noe in noelist:
                            keep=1  #assume noe is fine...
                            if(noe[0]==pk and pok!=pk): #if it links to one of our assignments...
                                tig=0
                                for pk2,val2,d2 in noes[pk]: #look in the reference for current peak...
                                    if(d2=='f'):
                                        if(pok==pk2 and noe[2]=='b'):  #matches one.
                                            #print('   ','good!')
                                            tig=1
                                    if(d2=='b'):
                                        if(pok==pk2 and noe[2]=='f'):  #mathches one.
                                            #print('   ','good!')
                                            tig=1
                                if(tig==0): #if no matches, it's bad. kill it!
                                    print ('   ','bad noe:',pok,noe)
                                    keep=0
                            if(keep==1):
                                noeListNew.append(noe)
                        noes[pok]=noeListNew



            ################################
            ####### Lets go pistachio!######
            print('-----------------------')
            print("Running assignment.....")

            assOpt=copy.deepcopy(mol.assRef)  #current assignment options
            InvAss=self.MakeInvAss(assOpt)    #inverted assignment options
            noes=copy.deepcopy(mol.G1edges)   #all confirmed edges
            #noes=PurgeNOES(0.2,noes)
            self.assStatus={}
            atoms=self.CreateAtomList(orph=True)  #get list of atoms
            atoms=atoms[1:] #cannot assign N terminus?.

            guess=self.ReadGuess()       #read in a list if guesses if available.
            #print(atoms)
            #return
            #prepare a fix list if needed
            #fix=[]
            #fix.append(('10','91T'))
            fix=[]
            #for key,vals in guess.items():
            #    fix.append((key,vals))

            #peaks=mol.OrderPeaks(noes,skip=False)  #order the peaks according to detectable segments
            peaks,skip=mol.OrderPeaks(noes,skip=[])  #remove orphans skip=[], keep orphans skip=False
            
            print('peaks (%i):' % len(peaks))
            print('atoms (%i):' % (len(atoms)))
            self.CountGraphs(noes)
            
            self.CheckOptionsForGuess(guess,InvAss)
            #return
            #ScreenForCTOCSY(InvAss)
            #CheckOptionsForGuess(guess,InvAss)

            from datetime import datetime
            t1=datetime.now()
            
            #print('AAA', assOpt['158'+mol.seq[158]])
            #do the job!
            assNew={}

            peaks,assNew,assGreat,cost,noes,assOpt,assInv,bad=self.RunCycles(peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=fix)
            #print(self.assStatus['112'])
            #return
            

            #look for permutations in the networks.

            #self.DoConnectionCheck(peaks,atoms,noes,assNew,assOpt,InvAss,assGreat)
            
            peaks,noes=self.RefineGraph(peaks,noes,assNew)  #look for duplicate NOEs and purge. Can create new orphans.
            peaks,assNew,assGreat,cost,noes,assOpt,assInv,bad=self.RunCycles(peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=fix)

            
            

            #self.DoConnectionCheck(peaks,atoms,noes,assNew,assOpt,InvAss,assGreat)

            #redefine all assignments that aren't bad as great.
            #assOpt,InvAss,assGreat=self.MakeGreat(peaks,assNew,assOpt,assGreat,bad)
            #restore options for the bad peaks
            #print('bad',bad)

            #restore assignment options and NOEs for the bad ones.
            assOpt=copy.deepcopy(mol.assRef)  #restore assignment options.
            assOpt,InvAss,assGreat=self.MakeGreat(peaks,assNew,assOpt,assGreat,bad)  #adjust options, fixing greats.
            RestoreBadNOEs(bad,noes)  #restore NOEs to non-great assignments.
            CleanNOEs(assGreat,noes)  #remove NOE options from
            

            peaks,assNew,assGreat,cost,noes,assOpt,assInv,bad=self.RunCycles(peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=fix,permissive=False)
            #print('BBB', assOpt['158'+mol.seq[158]],noes['158'],noes['157'])
 
            assOpt,InvAss,assGreat=self.MakeGreat(peaks,assNew,assOpt,assGreat,bad)
            
            
            self.DoConnectionCheck(peaks,atoms,noes,assNew,assOpt,InvAss,assGreat)

            peaks,assNew,assGreat,cost,noes,assOpt,assInv,bad=self.RunCycles(peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=fix,permissive=True)
            assOpt,InvAss,assGreat=self.MakeGreat(peaks,assNew,assOpt,assGreat,bad)

            if('ctocsy' in mol.spec.keys()):
                #sets=mol.spec['ctocsy'].GetSets(peaks,CSPlim=0.1) #go over peaks and find those that are close in HN (within CSP)
                cspH=Parse('assignParFile','cspH')
                if(cspH==False):
                    cspH=0.1
                cspN=Parse('assignParFile','cspN')
                if(cspN==False):
                    cspN=0.5
                sets=mol.spec['ctocsy'].GetSets(peaks,cspH=float(cspH),cspN=float(cspN)) #go over peaks and find those that are close in HN (within CSP)
                
                #try to assign unassigned peaks within 'match' of another peak in the overlap set.
                self.CheckExcludeOverlap(sets,'ctocsy',0.2)   #find peaks that can be killed.
                self.CheckExcludeOverlap(sets,'hcconh',0.02)
            

            

            #redefine all assignments that aren't bad as great.
            #assOpt,InvAss,assGreat=self.MakeGreat(peaks,assNew,assOpt,assGreat,bad)
            #to save
            self.assFin=assNew

            self.assFinInv={}
            for key,val in assNew.items():
                self.assFinInv[val]=key

            self.assConf=assGreat 
            self.assNOES=noes
            self.assPeaks=peaks
            self.assAtoms=atoms
            self.assOpt=assOpt
            self.assCost=cost
 
            stryLine=[]
            
            
            
            print()
            print('-----------------------')
            print("Assignments: (conf=%i)" % (len(assGreat.keys())))
            print('-----------------------')



            hdr='%5s %5s  %15s %7s %8s ' % ('Peak','Atom','Category','SCR','Confident')
            if('ctocsy' in mol.spec.keys()):
                hdr+='%3s %20s %6s %3s ' % ('i-1','cTOCSY','SCR','Ex')
            if('hcconh' in mol.spec.keys()):
                hdr+='%3s %32s %6s %3s ' % ('i-1','hcconh','SCR','Ex')
            print(hdr)
            
            

            cols='cccccc'
            hdr='\\textbf{%5s} & \\textbf{%5s} & & \\textbf{%15s} & \\textbf{%7s} & \\textbf{%8s}' % ('Peak','Atom','Category','SCR','Confident')
            if('ctocsy' in mol.spec.keys()):
                cols+='clcc'
                hdr+='& \\textbf{%3s} & \\textbf{%20s} & \\textbf{%6s} & \\textbf{%3s} ' % ('i-1','cTOCSY','SCR','Ex')
            if('hcconh' in mol.spec.keys()):
                cols+='clcc'
                hdr+='& \\textbf{%3s} & \\textbf{%32s} & \\textbf{%6s} & \\textbf{%3s} ' % ('i-1','hcconh','SCR','Ex')
            hdr+='\\\\\n'

            self.latexTable=[]
            self.latexTable.append('\\begin{tiny}\n')
            self.latexTable.append('\\noindent \\begin{tabular}{%s}' % cols)
            self.latexTable.append(hdr)
            self.latexTable.append('\\hline\\\\\n')


            for i,peak in enumerate(peaks):
                self.DoReport(i,peak,peaks,atoms,assNew,noes,assOpt,assGreat=assGreat,great='report',stryLine=stryLine,latex=True)


            self.latexTable.append('\\end{tabular}\n\n')
            self.latexTable.append('\\end{tiny}\n\n')

            if(len(stryLine)>0):
                for stry in stryLine:
                    print(stry)

            atomRemain=[]
            for atom in atoms:
                if(atom not in assGreat.values()):
                    atomRemain.append(atom)

            peakRemain=[]
            for pk in peaks:
                if(pk not in assNew.keys()):
                    peakRemain.append(pk)
                     

            poks,skop=mol.OrderPeaks(noes,skip=[])  
            unass=[]
            for pk in poks:
                if(pk not in peaks):
                    unass.append(pk)
            #print("Unassigned:")
            #print(peakRemain) 
            print("Unassignable:")
            print(unass)
            print('Final orphans:')
            print(skop)
            self.assOrph=skop
            self.assUnAss=unass
            self.assRemain=atomRemain



            #exclude='15','40','22'
            exclude=[]
            print("Unassigned atoms:")
            print(atomRemain)   
            pkO=[]
            for pk in skop:
                if(pk in exclude):
                    continue
                pkO.append(pk)
            """
            cost=numpy.zeros((len(pkO),len(atomRemain)))
            cost=cost.transpose()
            for j,at in enumerate(atomRemain):
                for i,pk in enumerate(pkO):
                    prob=1
                    if(at in assOpt.keys()):
                        for v in assOpt[at]: 
                                #print (v,pk)
                                if(v[0]==pk):
                                    prob=v[1]
                    #print(prob)
                    probNew,seen,expect,found,rem,c=self.DoTOCSYscore('C',pk,at,'ctocsy')
                    cost[j,i]=probNew*prob
            from scipy.optimize import linear_sum_assignment
            #ass=linear_sum_assignment((1-cost)**0.5)
            #ass=linear_sum_assignment(numpy.fabs(1-cost))
            #print(cost)
            ass=linear_sum_assignment(cost,maximize=True)
            #assNew={}
                
            
            for b,a in zip(ass[0],ass[1]):
                 print(pkO[a],atomRemain[b])
            """
            
            
            #get CSP overlap sets.

            
            #for pk,specs in mol.peak.items():
            #    for spec,pk3s in specs.items():
            #        for pk3 in pk3s:
            #            print(pk,spec,pk3.name,pk3.tp,pk3.f3p)
            #to handle the excess peaks:
            #39, we see two value CGs, but shiftx predicts 1 in range.
            #1. same for L, eg 48
            #2. sort out the Valine/Leucine overlap sets. Ie if we have extra peaks and we are L or V, re-run without merging in C.
            #3. high cost scores for some residue such as 18. why? some folding values need attention.
            #4. The doubly folded T peaks? Are they really threonines not alanines? check.
            #5. Move some functions back into assign_main where appropriate.

            #todo. First extend the Htocsy print list.

            #return
            #print (peaks)
            #print('guess:',len(guess))
            t2=datetime.now()  
            print('Time taken:',t2-t1)
            #print((t2-t1)*graphs)
            
            #re-order peaks according to assignment

            mol.DoG1histograms('FinalG1edges',noes)
            mol.WriteTolHist()
            Xs,Ys,Zs=self.ExpandSquare(cost)

            #print(Xs,Ys,Zs)
            print('done')
            #return

         

            
            #normalise
            #for i,peak in enumerate(atoms):
            #    rowsum=numpy.sum(Zs[2*i,:])
            #    if(rowsum!=0):
            #        Zs[2*i,:]=Zs[2*i,:]/(rowsum/2)
            #    rowsum=numpy.sum(Zs[2*i+1,:])
            #    if(rowsum!=0):
            #        Zs[2*i+1,:]=Zs[2*i+1,:]/(rowsum/2)

            #print (Zs)
            #from matplotlib.ticker import MultipleLocator
            from mpl_toolkits.axes_grid1 import SubplotDivider
            from mpl_toolkits.axes_grid1.mpl_axes import Axes as LocatableAxes

            #contour levels
            max_level=1.0
            min_level=0.0
            ctr_level=100
            levels=[]
            for i in range(ctr_level):
                levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))


            #make colorscheme
            cdict = {'red':   [(0.0,  1.0, 1.0),
                        (0.5,  1.0, 1.0),
                        (1.0,  0.0, 0.0)],

                'green': [(0.0,  1.0, 1.0),
                        (0.25, 1.0, 1.0),
                        (0.75, 0.0, 0.0),
                        (1.0,  0.0, 0.0)],

                'blue':  [(0.0,  1.0, 1.0),
                        (0.5,  0.0, 0.0),
                        (1.0,  0.0, 0.0)]}

            #setup matplotlib
            #pdf=PdfPages(analdir.split('/')[0]+'.pdf')
            #if(os.path.exists('report')==0):
            #    os.system('mkdir report')
            #pdf=PdfPages('report/correct.pdf')
            #fig = py.figure()
            #ax=fig.add_subplot(111)

            my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
            #plt.grid(True)
            self.axes5.set_xlabel('peakID')
            self.axes5.set_ylabel('atomID')
            #self.axes5.set_title(analdir)
            self.axes5.contourf(Ys, Xs, Zs,levels,cmap=my_cmap) #plot
            #for sq in line: #plot squares if atom refernesces are there, but no peak
            #    self.axes5.plot(sq[:,0],sq[:,1],c='k')


            #add confident and correct text labels
            #Xpos=numpy.max(Xs)*0.85
            #Ypos=numpy.max(Ys)
            #self.axes5.text(Xpos,Ypos*0.08,'confident: '+str(conf),fontsize=8,ha='left',va='center')
             #if(corr!='n'):
            #    self.axes5.text(Xpos,Ypos*0.05,'correct: '+str(corr),fontsize=6,ha='left',va='center')


            #sort out the axes
            self.axes5.tick_params(axis='x',which='both',length=2)
            self.axes5.tick_params(axis='y',which='both',length=2)
            xtics=numpy.arange(len(atoms))
            plt.setp(self.axes5.xaxis.get_majorticklabels(),rotation=90,fontsize=6)
            self.axes5.set_xticks(xtics,labels=atoms)
            ytics=numpy.arange(len(peaks))
            plt.setp(self.axes5.yaxis.get_majorticklabels(),fontsize=6)
            self.axes5.set_yticks(ytics,labels=peaks)


         
            #make the colorbar
            a=numpy.linspace(levels[0],100.,256).reshape(1,-1)
            a=numpy.vstack((a,a))
            run=10
            cols=5
            targ=45
            divider=SubplotDivider(self.fig,run,cols,targ)
            ax_cb=LocatableAxes(self.fig,divider.get_position())
            self.fig.add_axes(ax_cb)
            im=ax_cb.imshow(a,cmap=my_cmap,aspect=10)
            cb=plt.colorbar(im,cax=ax_cb)
            plt.axis("off")
            plt.setp(cb.ax.get_yticklabels(),visible=False)
            plt.setp(cb.ax.get_xticklabels(),visible=False)
            #self.Fit()
            #py.show()
            #fig.savefig(pdf,format='pdf')
            #pdf.close()


            """
            self.axes6 = self.fig2.add_subplot(111)

            infile='chsqc/raw/test.ft2'
            self.dic_proj,self.data_proj=ng.pipe.read(infile)

            Size=self.data_proj.shape
            
            uc0 = ng.pipe.make_uc(self.dic_proj,self.data_proj,dim=0)
            uc1 = ng.pipe.make_uc(self.dic_proj,self.data_proj,dim=1)
            ord_proj=self.dic_proj['FDDIMORDER']
            lab1_proj=self.dic_proj['FDF1LABEL']
            lab2_proj=self.dic_proj['FDF2LABEL']
            #lab3_proj=self.dic['FDF3LABEL']
            lab_proj=lab1_proj,lab2_proj, #lab3_proj
            #print(ord_proj)
            #print(lab_proj)
            self.labb_proj=lab_proj[int(ord_proj[1])-1],lab_proj[int(ord_proj[0])-1]

            self.uc0max=uc0.ppm(0)
            self.uc0min=uc0.ppm(Size[0]-1)
            self.uc1max=uc1.ppm(0)
            self.uc1min=uc1.ppm(Size[1]-1)

            

            print('--------------------------------------------------')
            print('Reading:',infile)
            print("Spectrum dimensions (pts): ",Size)   #print(the spectral dimensions)
            print("Labels: ",self.labb_proj)
            print("dimension 0 limits (ppm): ", self.uc0min, self.uc0max)  #carbon)
            print("dimension 1 limits (ppm): ", self.uc1min, self.uc1max)  #direct)
            
            #print('Maximum Intensity:',self.dmax)
            self.index0=[]#make index of carbon chemical shifts for index 0
            for i in range((Size[0])):
                self.index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
            self.index1=[]#make index of carbon chemical shifts for index 1
            for i in range((Size[1])):
                self.index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
            #self.index2=[]#make index of carbon chemical shifts for index 2
            #for i in range((Size[2])):
            #    self.index2.append((uc2.ppm(0)-i*(-uc2.ppm(Size[2]-1)+uc2.ppm(0))/(Size[2]-1)))
            self.index0=numpy.array(self.index0)
            self.index1=numpy.array(self.index1)
            #self.index2=numpy.array(self.index2)
            self.alias0=numpy.max((self.uc0min,self.uc0max))-numpy.min((self.uc0min,self.uc0max))+numpy.fabs(self.index0[0]-self.index0[1])
            self.alias1=numpy.max((self.uc1min,self.uc1max))-numpy.min((self.uc1min,self.uc1max))+numpy.fabs(self.index1[0]-self.index1[1])
            #self.alias2=numpy.max((self.uc2min,self.uc2max))-numpy.min((self.uc2min,self.uc2max))+numpy.fabs(self.index2[0]-self.index2[1])
            self.YY,self.XX=numpy.meshgrid(self.index1,self.index0)

            def GetLevels(min_level,fac,ctr_level):
                    levels=[]
                    levels.append(min_level)
                    for i in range(ctr_level-1):
                        levels.append(levels[i]*fac)
                    levels=numpy.array(levels)
                    levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
                    return levels


            levels=GetLevels(2E6,1.2,10)

            self.XX+=5.537-5.6494
            print(self.XX)
            self.YY+=55.04-54.4

            print(self.YY)
            self.axes6.contour(self.XX,self.YY,self.data_proj,levels,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #plot pdb network

            for pk,specs in mol.peak.items():
                if('ctocsy' not in specs.keys()):
                    continue
                if('hcconh' not in specs.keys()):
                    continue

                for l in 'A',:#'B','G','D','E':
                    y=False
                    x=False
                    for pk3 in mol.peak[pk]['ctocsy']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            y=pk3.f3p
                            break
                    for pk3 in mol.peak[pk]['hcconh']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            x=pk3.f3p
                            break
                    if(x!=False and y!=False):
                        self.axes6.scatter(x,y,c='k',marker='x',zorder=2)

                for l in 'B',:#'B','G','D','E':
                    y=False
                    x=False
                    for pk3 in mol.peak[pk]['ctocsy']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            y=pk3.f3p
                            break
                    for pk3 in mol.peak[pk]['hcconh']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            x=pk3.f3p
                            break
                    if(x!=False and y!=False):
                        self.axes6.scatter(x,y,c='b',marker='x',zorder=2)

                for l in 'G',:#'B','G','D','E':
                    y=False
                    x=False
                    for pk3 in mol.peak[pk]['ctocsy']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            if(y==False):
                                y=[]
                            y.append(pk3.f3p)
                            #break
                    for pk3 in mol.peak[pk]['hcconh']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            if(x==False):
                                x=[]
                            x.append(pk3.f3p)
                            #break
                    if(x!=False and y!=False):
                        for a in x:
                            for b in y:
                               self.axes6.scatter(a,b,c='r',marker='x',zorder=2)

                for l in 'D',:#'B','G','D','E':
                    y=False
                    x=False
                    for pk3 in mol.peak[pk]['ctocsy']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            if(y==False):
                                y=[]
                            y.append(pk3.f3p)
                            break
                    for pk3 in mol.peak[pk]['hcconh']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            if(x==False):
                                x=[]
                            x.append(pk3.f3p)
                            break
                    if(x!=False and y!=False):
                        for a in x:
                            for b in y:
                                self.axes6.scatter(a,b,c='g',marker='x',zorder=2)
                    elif(y!=False):
                        pass
                        #print('cannot match:',l,pk,y,y-(mol.spec['ctocsy'].uc0max-mol.spec['ctocsy'].uc0min))

                for l in 'E',:#'B','G','D','E':
                    y=False
                    x=False
                    for pk3 in mol.peak[pk]['ctocsy']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            y=pk3.f3p
                            #break
                    for pk3 in mol.peak[pk]['hcconh']:
                        if(len(pk3.tp)<2):
                            continue
                        tp=pk3.tp[1]
                        if(tp==l):
                            x=pk3.f3p
                            #break
                    if(x!=False and y!=False):
                        self.axes6.scatter(x,y,c='k',marker='x',zorder=2)
                    elif(y!=False):
                        print('cannot match:',l,pk,y,y-(mol.spec['ctocsy'].uc0max-mol.spec['ctocsy'].uc0min))
            """



            print('hello!')
            #print (self.parent.parent.tabOne.molecule.assRef)

        self.canvas.draw()
        self.canvas2.draw()

        """
        #update parent
        self.parent.parent.tabTwo.ComboBox2.SetSelection(self.index[self.places['i'][0][0]])
        if('i+1' in self.places.keys()):
            if(len(self.places['i+1'])>0):
                self.parent.parent.tabTwo.ComboBox1.SetSelection(self.index[self.places['i+1'][0][0]])
        if('i-1' in self.places.keys()):
            if(len(self.places['i-1'])>0):
                self.parent.parent.tabTwo.ComboBox3.SetSelection(self.index[self.places['i-1'][0][0]])
        self.parent.parent.tabTwo.draw_figure()
        """

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

    def onAssButton(self,event):

        sele1=self.listy[self.ComboBox1.GetSelection()]
        sele2=self.G1_nodes[self.ComboBox2.GetSelection()]


        if(sele1 in self.results.keys()):
            print(sele1,'has already been assigned')
            print(self.results.items())
            #import sys
            #sys.exit(100)
            return
        if(sele2 in self.assSeq.keys()):
            print(sele2,'has already been assigned')
            return

        self.assSeq[sele2]=1

        num_items = self.source3.GetItemCount()
        self.source3.InsertItem(num_items,str(num_items))
        self.source3.SetItem(num_items,0,sele1)
        self.source3.SetItem(num_items,1,sele2)


        assign=self.DoAssign(sele1,sele2)
        self.PopulateResults()
        self.PopulateList()
        self.draw_figure()


    #Find unassigned peaks.
    #If that peak is in an overlap set...
    #look for assigned peaks in the other peaks
    #within MATCH of the set.
    #if so, remove the peak from the first set as this is almost certainly leakage.
    def CheckExcludeOverlap(self,sets,spec,match):
        mol=self.parent.parent.tabOne.molecule
        for set in sets: #for each set...
            print(set)   #print the set.
            for pk in set: #for each peak in the set...
                if(spec in mol.peak[pk].keys()):
                    kill=[]
                    for pk3 in mol.peak[pk][spec]: #for each peak in the set, see if we can find matches.
                        if(len(pk3.tp)==0): #if the entry has no typ label
                            #is this unexplained peak assigned within one of the others?
                            print('unexplained:',spec,pk,pk3.name,pk3.f3p,pk3.tp)
                            f3=pk3.f3p
                            for pok in set: #go through the set again...
                                if(pok==pk):   #skip if we have the self-match
                                    continue
                                if(spec in mol.peak[pok].keys()):  #now go through peaks in the new slice...
                                    vals=[]
                                    for pk4 in mol.peak[pok][spec]:
                                        vals.append(numpy.fabs(pk4.f3p-f3))
                                    argy=numpy.argmin(vals)   #get the peak that is closest in chemical shift to target f3
                                    if(vals[argy]<match and len(mol.peak[pok][spec][argy].tp)>0):  #if the match is within tolerance and labelled...
                                        print('MATCH:',vals[argy],pok,mol.peak[pok][spec][argy].name,mol.peak[pok][spec][argy].f3p,mol.peak[pok][spec][argy].tp)
                                        kill.append(pk3.name) #kill the original unlabelled peak.
                                        break
                    test=[] #now go through the peaks again...
                    for pk3 in mol.peak[pk][spec]:
                        if(pk3.name not in kill):  #do not kill this one.
                            test.append(pk3)
                        else:   #but if name in kill list, get rid.
                            print('killing',spec,'peak:',pk,pk3.name)
                    mol.peak[pk][spec]=test   #set new tuple for the peak.



    #take 2D slice and expand to make a square for pretty plotting.
    def ExpandSquare(self,cost):
                #turn assignment into plottable square
                ii=cost.shape[0]
                jj=cost.shape[1]
                #Zs=numpy.zeros((len(peaks)*2,len(atoms)*2))
                #Xs=numpy.zeros((len(peaks)*2,len(atoms)*2))
                #Ys=numpy.zeros((len(peaks)*2,len(atoms)*2))
                Zs=numpy.zeros((ii*2,jj*2))
                Xs=numpy.zeros((ii*2,jj*2))
                Ys=numpy.zeros((ii*2,jj*2))

                for i in range(ii): #count up to length of atoms
                    for j in range(jj): #count up to lenght of atoms
                        
                        Xs[2*i+0,2*j+0]=i-0.48
                        Xs[2*i+1,2*j+0]=i+0.48
                        Ys[2*i+0,2*j+0]=j-0.48
                        Ys[2*i+1,2*j+0]=j-0.48
                            
                        Xs[2*i+0,2*j+1]=i-0.48
                        Xs[2*i+1,2*j+1]=i+0.48
                        Ys[2*i+0,2*j+1]=j+0.48
                        Ys[2*i+1,2*j+1]=j+0.48

                        prob=cost[i,j]

                        Zs[2*i+0,2*j+0]=prob
                        Zs[2*i+1,2*j+0]=prob
                        Zs[2*i+0,2*j+1]=prob
                        Zs[2*i+1,2*j+1]=prob
                return Xs,Ys,Zs

    #make cost matrix from chemical shift probabilities      
    def MakeCostSquare(self,peaks,atoms,noes,assOpt,scale,depth):
                cost=numpy.zeros((len(peaks),len(atoms)))
                cnty=numpy.zeros((len(peaks),len(atoms)))
                for i,peak in enumerate(peaks): #count up to length of atoms
                    for j,atom in enumerate(atoms): #count up to lenght of atoms
                        if(atom in assOpt.keys()):
                            for v in assOpt[atom]:
                                
                                if(v[0]==peak):
                                    #if(peak=='158' and atom=='158S'):
                                    #    print('  ','unadjusted cost:',peak,atom,v)

                                    prob=v[1]
                                    #print('bbbbbb',peak,atom,prob)
                                    #if(peak=='11'):
                                    #    if(atom=='96R'):
                                    #        print ('AAAAAA',peak,atom,prob)

                                    #if(prob>1E-3):
                                    prob,cnt=self.AdjustProbability(peak,j,prob,noes,atoms,len(atoms),assOpt,scale,depth)

                                    #if(peak=='158' and atom=='158S'):
                                    #    print('  ','adjusted cost:',peak,atom,v,prob)
                                    cost[i,j]=prob
                                    cnty[i,j]=cnt
                                    break
                #for all the probabilities avaiable for a peak.
                #if some assignment possibilities see more links, these are probably more probable.
                for i,peak in enumerate(peaks):
                    maxy=numpy.max(cnty[i,:])
                    mask=cnty[i,:]==maxy
                    cost[i,:]*=mask
                return cost

    def UpdateAss(self,assNew,a,b,s):                    
        
        tick=0
        if(a not in assNew.keys()):
            tick=1
        elif(a not in self.assStatus.keys()):
            tick=1
        elif(assNew[a]!=b):
            tick=1
        if(tick==1):
            self.assStatus[a]=s
        assNew[a]=b

    #take cost square, optimise, re-order peak list/cost matrix.
    def OptimiseCost(self,peaks,atoms,noes,assOpt,InvAss,assNew,guess,scale=1.0,fix=[],depth=1):
                cost=self.MakeCostSquare(peaks,atoms,noes,assOpt,scale,depth)  #make a cost square from assOpt
                #print('COST',cost[numpy.where(peaks=='11')[0][0],numpy.where(atoms=='96R')[0][0]])
                if(len(fix)>0):
                    for f in fix:
                        pkF=f[0]
                        atF=f[1]

                        pkFi=numpy.where(pkF==peaks)[0][0]
                        atFi=numpy.where(atF==atoms)[0][0]
                        cost[pkFi,:]=0.0
                        cost[pkFi,atFi]=1.0

                #peaks=numpy.array(peaks)
                #print(numpy.where(peaks=='158')[0],numpy.where(atoms=='158S')[0])
                #print('COST',cost[numpy.where(peaks=='158')[0][0],numpy.where(atoms=='158S')[0][0]])

                from scipy.optimize import linear_sum_assignment
                #ass=linear_sum_assignment((1-cost)**0.5)
                #ass=linear_sum_assignment(numpy.fabs(1-cost))
                #print(cost)
                #print('COST',cost.shape)
                #peaks=numpy.array(peaks)
                #print('COST',cost[numpy.where(peaks=='86')[0][0],numpy.where(atoms=='151T')[0][0]])
                #print('COST',cost[numpy.where(peaks=='86')[0][0],numpy.where(atoms=='85V')[0][0]])

                cost=cost.transpose()
                ass=linear_sum_assignment(cost,maximize=True)
                #assNew={}
                
                discrep=[]
                newpunt=0
                good=0
                for b,a in zip(ass[0],ass[1]):
                    #print(a,b,len(peaks),len(atoms))
                    
                    try:
                        #print(peaks[a],atoms[b],guess[peaks[a]])
                        if(atoms[b]!=guess[peaks[a]]):
                            discrep.append((peaks[a],atoms[b],guess[peaks[a]]))
                        else:
                            good+=1
                    except:
                        #print(peaks[a],atoms[b])
                        newpunt+=1

                    self.UpdateAss(assNew,peaks[a],atoms[b],'o')
    
                #print('COST',cost[numpy.where(peaks=='11')[0][0],numpy.where(atoms=='96R')[0][0]])

                print('Optimisation run results:')
                print("     Matches:",good)
                print("     New punts:",newpunt)
                print('     Discrepancies:',len(discrep),discrep)
                print("     Total:",len(discrep)+good+newpunt)
                
                SHOWDISC=False
                if(SHOWDISC):
                    for disc in discrep:
                        print()
                        print(disc)
                        #try:
                        #    print(assOpt[disc[1]])
                        #    print(assOpt[disc[2]])
                        #except:
                        #    pass
                        if(disc[0] not in InvAss.keys()):
                            print('wtf?')
                            print(sorted(list(InvAss.keys())))
                            continue
                        if(disc[1] not in InvAss[disc[0]].keys()):
                            print('  ',disc[1],'not an option')
                        else:
                            print('  ',InvAss[disc[0]][disc[1]])
                        if(disc[2] not in InvAss[disc[0]].keys()):
                            print('  ',disc[2],'not an option')
                        else:
                            print('  ',InvAss[disc[0]][disc[2]])

                        ii=numpy.where(peaks==disc[0])[0][0]
                        jj=numpy.where(atoms==disc[1])[0][0]
                        print('placed cost:',cost[ii,jj],cost.shape,peaks[ii],atoms[jj])
                        jj=numpy.where(atoms==disc[2])[0][0]
                        print('our cost:',cost[ii,jj],cost.shape,peaks[ii],atoms[jj])

                        i=numpy.where(disc[0]==peaks)[0][0]
                        self.DoReport(i,disc[0],peaks,atoms,assNew,noes,assOpt)

                    #print(InvAss[disc[0]])
                
                #print('COST',cost[numpy.where(peaks=='11')[0][0],numpy.where(atoms=='96R')[0][0]])

                peaksNew=[]
                for atom in atoms:
                    for b,a in zip(ass[0],ass[1]):
                        if(atoms[b]==atom):
                            peaksNew.append(peaks[a])    
                peaks=numpy.array(peaksNew)
                cost=self.MakeCostSquare(peaks,atoms,noes,assOpt,scale=scale,depth=depth)
                cost=cost.transpose()
                return peaks,assNew,cost

    #create a list of atoms. specify whether to include orphans and prolines
    def CreateAtomList(self,P=False,orph=True):
        atoms=[]
        seq=self.parent.parent.tabOne.molecule.seq
        for key,vals in seq.items():
            if(P==False and vals=='P'):
                continue
            if(orph==False): #if we don't want the orphans...
                tick=0
                #print (type(key))
                #print(key,str(int(key)+1),str(int(key)+1) in seq.keys(),seq[str(int(key)+1)]!='P')
                if((int(key)+1) in seq.keys() and seq[(int(key)+1)]!='P'):
                    tick=1
                if((int(key)-1) in seq.keys() and seq[(int(key)-1)]!='P'):
                    tick=1
                if(tick==0):
                        continue
            atoms.append((str(key)+vals))
        return numpy.array(atoms)

    #find something? not sure!      
    def GetBMRB(self,atom,typ,res1):   
        if(len(res1)==1):
            res=self.parent.parent.tabOne.molecule.p1to3[res1]
        else:
            res=res1
        for r in self.parent.parent.tabOne.molecule.bmrb[atom][typ]:
            if(r[0]==res):
                return r[1:]
        return False  
    """                  
    def MatchTOCSY(self,nuc,resn,peak,shifts):
                if(nuc=='C'):
                    spec='ctocsy'
                    bmrb=self.parent.parent.tabOne.molecule.bmrbC
                    med=self.parent.parent.tabOne.molecule.CTOCSYmed
                elif(nuc=='H'):
                    spec='hcconh'
                    bmrb=self.parent.parent.tabOne.molecule.bmrbH   
                    med=self.parent.parent.tabOne.molecule.HTOCSYmed 
                
                
                if(spec not in shifts.keys()):
                    return {}

                shifts=shifts[spec]

                valR=[]
                vilR=[]
                if(nuc=='C'):
                    for i,pk in enumerate(shifts): #take all but the already classified CA
                        if(pk.tp!='CA(i-1)'):
                            # print(pk.tp)
                            valR.append(pk.f3p)
                            vilR.append(i)
                else:
                    for i,pk in enumerate(shifts): #take all but the already classified CA
                            valR.append(pk.f3p)
                            vilR.append(i)
                #print('shifts',valR)
                typs=bmrb[resn].keys() # these are the Cas, Cbs etc of the given amino acid
                # print(typs)
                vol=[]
                vop=[]
                ti=[]
                if(nuc=='C'):
                    for t in typs:
                        if(t!='CA' and t!='C'): #take all shifts but CA.
                            if(bmrb[resn][t][0]<110):
                                if(resn=='LYS' and t=='CE'): #HNs on side chain not seen
                                    continue
                                vol.append(bmrb[resn][t][0]) #exclude carbonyls from TOCSY
                                vop.append((bmrb[resn][t][1]*1.5)**2.) # fudge factor for tolerance.
                                ti.append(t) #save type.
                else:
                    for t in typs:
                        if(t!='H' ): #take all shifts but CA.
                            if(bmrb[resn][t][0]<120):
                                vol.append(bmrb[resn][t][0]) #exclude carbonyls from TOCSY
                                vop.append((bmrb[resn][t][1]*1.5)**2.) # fudge factor for tolerance.
                                ti.append(t) #save type.
                ti=numpy.array(ti)
                # print(resn)
                # print(valR) # these are where the peaks are in the tocsy
                # print(vol) # these are where we expect to see the peaks in the tocsy for a given residue

                #reference all BMRB values to within tocsy range (try aliasing)
                #print('pre:',vol)
                from .assign_main import peakEntry
                for i,vo in enumerate(vol): #for each list from the BMRB...
                    pok=peakEntry(('test','0','0',vo+med,'1')) #create a fake peak entry...
                    self.parent.parent.tabOne.molecule.spec[spec].alias(pok,vo+med,0)  #and alias it to within the range of the spectrum...
                    vol[i]=pok.ppmI-med   #and save.
                #print('vols:',vol)
                
                res={}

                if(len(valR)>0 and len(vol)>0): #make sure residue has at least enough shifts to explain data
                    XX,YY=numpy.meshgrid(vol,valR) #create grids of constant row = bmrb and col = expt
                    diff=numpy.fabs(XX-YY)**2. #get the absolute value of their differences
                    from scipy.optimize import linear_sum_assignment
                    ass=linear_sum_assignment(diff)
                    #print('vol',vol)
                    #print('val',valR)
                    #if(len(ass)<len(valR)):
                    #    print('Cannot be this')
                    #print(ass)
                    for a,b in zip(ass[0],ass[1]):
                        #print (a,b)
                        prob=numpy.exp(-1.*(diff[a,b])/(2*vop[b])) #save probability
                        #print('dd',diff[a,b]**0.5,vop[b])
                        #print(ti[b],vol[b],shifts[vil[a]].f3p,shifts[vil[a]].tp,prob)
                        
                        res[ti[b]]=valR[a],vol[b],vop[b],prob #save probability
                        self.parent.parent.tabOne.molecule.peak[peak][spec][vilR[a]].tp=ti[b]+'(i-1)' #classify
                        #print(ti[b]+'(i-1)')
                    for ii in range(len(valR)): #force unassigned peaks to zero
                    #    #print('looking at ',i,val[i],ass[0])
                        if(ii not in ass[0]):
                    #        #print('yipyip') #forcing this to zero (crude)
                            res[nuc+'X'+str(ii)]=valR[ii]
                
                return res
    """
    def GetKeyFromVal(self,test,dicty):
        for key,val in dicty.items():
                if(val==test):
                    return key
        return False



            
    def ReadGuess(self,):
        guess={}
        if(os.path.exists('dat/guess.txt')==False):
            return {}
        inny=open('dat/guess.txt')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>=2):
                atom=test[1]
                for key,vals in self.parent.parent.tabOne.molecule.seq.items():
                    if(atom==str(key)):
                        lab=str(key)+vals
            
                guess[test[0]]=lab
        return guess

    def CheckOptionsForGuess(self,guess,InvAss):
                print('Making sure the expected assignment options are possible...')
                bad=0
                for pk,at in guess.items():
                    
                    tig=0
                    if(pk in InvAss.keys()):
                        if(at in InvAss[pk].keys()):
                            #print ('Good!')
                            tig=1
                            
                    if(tig==0):
                        bad+=1
                        print()
                        print(' Guess assignment:',pk,at)
                        print(" not in assOpt:",InvAss[pk])
                        print(' not good :(')
                    #print(InvAss[pk])   
                if(bad==0):
                    print('All expected assignments are possible.')
                else:
                    print('Ah! Not all assignments are possible!')    
                return
    """
    def ScreenForCTOCSY(self,InvAss):
                for peak,val in InvAss.items():
                    resnsA,probsA,probs2A=self.parent.parent.tabOne.molecule.AnalPeak(peak) 
                    resnsA=numpy.array(resnsA)
                    probNew=[]
                    for v in val.keys():
                        print(peak,v)

                        resn_i=v[-1]
                        resi_i=int(v[:-1])
                        
                        arg=numpy.where(resnsA==resn_i)[0][0]
                        print('I  :',resn_i,resi_i,probsA[arg])

                        resi_im1=str(resi_i-1)
                        resn_im1=''
                        for a in atoms:
                            if(a[:-1]==resi_im1):
                                resn_im1=a[-1]
                                break
                        if(resn_im1==''):
                            probNew.append(probsA[arg])
                            continue
                        
                        arg=numpy.where(resnsA==resn_im1)[0][0]
                        print('IM1:',resn_im1,resi_im1,probs2A[arg])
                        probN=probs2A[arg]*probsA[arg]
                        probNew.append(probN)

                    probNew=numpy.array(probNew)
                    if(numpy.sum(probNew)<1E-5):
                        print('all are bad!')
                        continue

                    valNew=[]
                    kill=[]
                    for i,v in enumerate(val.keys()):
                        if(probNew[i]>1E-4):
                            valNew.append((v[0],probNew[i]))
                        else:
                            kill.append((peak,v))
                    for k,v in kill:
                        print('killing:',k,v)
                        del InvAss[k][v]
    """
    #take the NOES. If there are repeated NOEs
    #and if one maps to two assignments that are +/-1, kill the rest.
    #should probably check these are both assGreats.
    def RefineGraph(self,peaks,noes,assNew):
            print('Pruning NOEs....')
            listy={}
            for pk,vals in noes.items():
                 for pk2,val,dir in vals:
                    lab=pk2+dir
                    if(lab not in listy.keys()):
                         listy[lab]=[]
                    listy[lab].append((pk,pk2,val,dir))
            kill=[]
            good=[]
            for pk,vals in listy.items():
                 if(len(vals)>1):
                      #print (pk,vals)
                      test=[]
                      
                      bad=[]
                      for v in vals:
                        pk,pk2,val,dir=v
                        if(pk not in assNew.keys()):
                             continue    
                        if(pk2 not in assNew.keys()):
                             continue
                            
                        atom1i=int(assNew[pk][:-1])
                        atom2i=int(assNew[pk2][:-1])
                        tick=0
                        if dir=='b':
                             if(atom2i+1==atom1i):
                                  tick=1
                        elif(dir=='f'):
                            if(atom2i-1==atom1i):
                                    tick=1
                        if(tick==1):
                             #print('good',atom1i,atom2i)
                             good.append((pk,pk2))
                        else:
                             #print('bad',atom1i,atom2i)
                             bad.append((pk,pk2))
                        test.append(tick)
                      if(numpy.sum(test)==1):
                           #print('one is good. ')
                           #print('killing',bad)
                           for b in bad:
                                kill.append(b)
            #print(len(noes.keys()))ßƒ
            for g in good:
                 print('Good links:',g)
            for pk,pk2 in kill:
                vals=noes[pk]
                valsNew=[]
                for (pk3,val,dir) in vals:
                     if(pk3==pk2): #kill
                          print('Removing redundant link: ',pk,pk2,val)

                     else:
                        valsNew.append((pk3,val,dir))
                noes[pk]=valsNew
                if(len(valsNew)==0):
                     print("Created a new orphan:",pk)


            peakNew=[]
            for peak in peaks:
                if(len(noes[peak])==0):
                    print('removing:',peak,assNew[peak])
                    del assNew[peak]
                    del self.assStatus[peak]
                    pass
                else:
                    peakNew.append(peak)
            peaks=peakNew
            return peaks,noes

    #take NOE list and where we have multiple strong connections between 
    #two nodes, work out total graphs
    def CountGraphs(self,noes,verb=True):
        print("Checking for multiple entries in graphs....")
        multi={}
        graphs=1
        for pk,vals in noes.items():
            cntb=0
            cntf=0
            for (pk2,val,dir) in vals:
                if(dir=='b'):
                    cntb+=1
                if(dir=='f'):
                    cntf+=1

            if(cntb>1):
                graphs*=cntb

                for (pk2,val,dir) in vals:
                    if(dir=='b'):
                        if(pk not in multi.keys()):
                            multi[pk]=[]
                        multi[pk].append((pk2,dir,val))
                        if(verb):
                            print (pk,pk2,dir,val)
                
            if(cntf>1):
                graphs*=cntf
                for (pk2,val,dir) in vals:
                    if(dir=='f'):
                        if(pk not in multi.keys()):
                            multi[pk]=[]
                        multi[pk].append((pk2,dir,val))
                        if(verb):
                            print (pk,pk2,dir,val)
        print('total graphs:',graphs)
        return graphs,multi

    #make a new assignment options dictionary
    #from a InvAss (peak ordered options)
    def MakeAssOptFromInv(self,peaks,atoms,InvAss):
        assOpt={}
        for atom in atoms:
            pklist=[]
            for peak in peaks:
                try:
                    scr=InvAss[peak][atom]
                    pklist.append((peak,scr))
                except:
                    pass
            if(len(pklist)!=0):
                assOpt[atom]=pklist
        return assOpt
            
    def MakeInvAss(self,assOpt):
        #invert assignment options square
        InvAss={}
        for key,vals in assOpt.items():
            for v,scr in vals:
                if(v not in InvAss.keys()):
                    InvAss[v]={}
                InvAss[v][key]=scr
        return InvAss

    #screen all NOEs below specified num value
    def PurgeNOES(self,num,noes):
        noeNew={}
        for key,vals in noes.items():
            v=[]
            for (pk2,val,dir) in vals:
                if(val<num):
                    v.append((pk2,val,dir))
            noeNew[key]=v
        return noeNew


    


    ########################################################
    #take another look at the NOEs and assignment options.
    #if we have an assignment that maps to the residues being
    #one appart, then keep this NOE and reject the rest.
    def CheckGraphs(self,noes,assNew):
        REFINE=False
        graphs=1
        for pk,vals in noes.items():
            cntb=0
            cntf=0
            for (pk2,val,dir) in vals:
                if(dir=='b'):
                    cntb+=1
                if(dir=='f'):
                    cntf+=1

            if(cntb>1):
                graphs*=cntb
                #print()
                valNew=[]
                tick=0
                for (pk2,val,dir) in vals: #unpack the NOES..
                    #print (pk,pk2,val,dir)
                    if(dir=='b'):  #if there are backward guys...
                        #print (pk,pk2,dir,val)
                        resi1=int(assNew[pk][:-1])  #get atom number for i
                        if(pk2 not in assNew.keys()):
                                continue
                        resi2=int(assNew[pk2][:-1]) #get atom number for connection
                        #print(resi1,resi2)
                        if(resi2+1==resi1): #if they are separated appropriately...
                            REFINE=True
                            #print('yay')
                            #this is the connection to keep.
                            print('REFINING GRAPH:',pk,pk2,dir,val)
                            valNew.append((pk2,val,dir))
                            tick=1
                if(tick==1):
                    for (pk2,val,dir) in vals: #keep the 'fs
                        if(dir=='f'):
                            valNew.append((pk2,val,dir))
                    noes[pk]=valNew

            if(cntf>1):
                #print()
                graphs*=cntf
                tick=0
                valNew=[]
                for (pk2,val,dir) in vals:
                    if(dir=='f'):
                        #print (pk,pk2,dir,val)
                        resi1=int(assNew[pk][:-1])
                        if(pk2 not in assNew.keys()):
                            continue
                        resi2=int(assNew[pk2][:-1])
                        if(resi2-1==resi1):
                            REFINE=True
                            #print('yay')
                            #this is the connection to keep.
                            print('REFINING GRAPH:',pk,pk2,dir,val)
                            valNew.append((pk2,val,dir))
                            tick=1
                if(tick==1):
                    for (pk2,val,dir) in vals: #keep the bs
                        if(dir=='b'):
                            valNew.append((pk2,val,dir))
                    noes[pk]=valNew


                
        #if(REFINE):                
        #    print('New total graphs:',graphs)
        return noes

    #fix a specific peak in the assOpt array.
    def PurgeAss(self,assOpt,peak):
        for key,vals in assOpt.items():
            vnew=[]
            for v in vals:
                if(v[0]!=peak):
                    vnew.append(v)
            assOpt[key]=vnew
        return assOpt

    #get the assignments that cleanly link two ends.
    #then try to elongate the two ends. keep going
    #until the assignment list stops elongating.
    def DoGreatClean(self,peaks,atoms,assNew,noes,assOpt):
        assGreat={}
        for i,peak in enumerate(peaks):
                #if(peak=='37'):
                self.DoReport(i,peak,peaks,atoms,assNew,noes,assOpt,assGreat=assGreat,great=True)
        #print(assGreat)
        #print(len(assGreat.keys()))
        RUN=False
        curr=len(assGreat.keys())
        while(1==1):
            #print('curr',curr)
            for i,peak in enumerate(peaks):
                success=self.DoReport(i,peak,peaks,atoms,assNew,noes,assOpt,assGreat=assGreat,great='extend')
                if(success):
                        RUN=True
            if(len(assGreat.keys())==curr):
                break
            curr=len(assGreat.keys())
        #print(assGreat)
        #print(len(assGreat.keys()))
        return assGreat,RUN
    
    #set all non-bad assignments as great
    def MakeGreat(self,peaks,assNew,assOpt,assGreat,bad):    
        for peak in peaks:
            if(peak not in bad):
                assGreat[peak]=assNew[peak]
        assOpt,InvAss=self.PurgeGreat(assOpt,assGreat)    
        return assOpt,InvAss,assGreat

    #Fix all assignments in assGreast within assOpt
    def PurgeGreat(self,assOpt,assGreat):
        #now fix all the assignments in assGreat and redo.
        for peak,atom in assGreat.items():
            #print(peak,atom)
            assOpt=self.PurgeAss(assOpt,peak)
            assOpt[atom]=[]
            assOpt[atom].append((peak,1.)) #make only 1 option for this atom.
        InvAss=self.MakeInvAss(assOpt)
        return assOpt,InvAss

    def DoCommonConnect(self,peaks,atoms,noes,assNew,assOpt,InvAss,assGreat,permissive=False):
        linky={}
        bad=[]
        RUN=False
        for i,peak in enumerate(peaks):
            stry,prob,great=self.DoReport(i,peak,peaks,atoms,assNew,noes,assOpt,assGreat=assGreat,great='report',verb=False)
        
            if(permissive==False and stry!='GREAT!'):  #only great is great. strict on what is great.
                bad.append(peak)
            if(permissive==True and stry=='BAD'):  #forward and backward great are great.
                bad.append(peak)
            
                    
            linky[int(assNew[peak][:-1])]=stry
        mol=self.parent.parent.tabOne.molecule
        for resi,stry in linky.items():
            if(stry=='BACKWARD GREAT!'):
                if(resi+2 in linky.keys()):
                        if(linky[resi+2]=='FORWARD GREAT!'):

                            pkAssF=self.GetKeyFromVal(str(resi)+mol.seq[resi],assNew)
                            pkAssB=self.GetKeyFromVal(str(resi+2)+mol.seq[resi+2],assNew)  
                            segF=self.GetLinks(pkAssF,'f',noes)
                            segB=self.GetLinks(pkAssB,'b',noes)
                            common=numpy.intersect1d(segF,segB)
                            if(len(common)==1 and common[0] in bad):
                                #this should be assigned here.
                                print('ASSIGNING THROUGH COMMON CONNECTION:',common[0],str(resi+1)+mol.seq[resi+1])
                                assGreat[common[0]]=str(resi+1)+mol.seq[resi+1]
                                
                                self.UpdateAss(assNew,common[0],str(resi+1)+mol.seq[resi+1],'c')
                                #assNew[common[0]]=str(resi+1)+mol.seq[resi+1]
                                #self.assStatus[common[0]]='c'
                                RUN=True

                                badNew=[]
                                for b in bad:
                                    if(b!=common[0]):
                                        badNew.append(b)
                                bad=badNew
                                #print("COMMON:",pkAssF,pkAssB,common)

        assOpt,InvAss=self.PurgeGreat(assOpt,assGreat)
        
        return assOpt,InvAss,assGreat,bad,RUN

    def DoConnectionCheck(self,peaks,atoms,noes,assNew,assOpt,InvAss,assGreat):

        mol=self.parent.parent.tabOne.molecule       

        def GetConn(pk1,pk2,dir,specs,tols):
            hits=[]
            for i,spec in enumerate(specs): #for connections in each spectrum...
                if(pk1 in spec.keys()): 
                    for j, edge in enumerate(spec[pk1]):  #for all of the edges 
                        if(numpy.fabs(edge[1])<=tols[i]):   #if the score is within the tolerance...
                            if(edge[2]==dir and edge[0]==pk2):
                                print('Found:',dir,pk1,i,edge)
                                hits.append((i,edge))
            return hits

        def CheckConns(pkAssA,pkAssB):
            hitsF=GetConn(pkAssA,pkAssB,'f',mol.spec_edges_forwards,mol.spec_tols_forwards)
            hitsB=GetConn(pkAssB,pkAssA,'b',mol.spec_edges_backwards,mol.spec_tols_backwards)
            check=[]
            for i,sp in enumerate(mol.spec_edges_ref):
                iscr=0
                for h in hitsF:
                    if h[0]==i:
                        iscr+=1
                for h in hitsB:
                    if h[0]==i:
                        iscr+=1
                if(iscr==2):
                    print('Connection found:',sp)
                else:
                    print('Connection not found:',sp)  
                    check.append(sp)
            return check

        """
        def AtoI(val,index_dataRef):
                        for i,ind in enumerate(index_dataRef):
                            if(ind==val):
                                return i
                        return -1                                               
        def index(array):
                    index=[]
                    for i in range(len(array)):
                        index.append((array[i].name))
                    return index
        def findnear_index(test,array):
                    #array = numpy.asarray(array)
                    idx = (numpy.abs(array - test)).argmin()
                    return idx
        """
 
        def DoCOcheck(pkAssA,pkAssB,stryLines):
                print('   checking CO connections...')
                
                pkAF=mol.GetPeakTp('hncaco',pkAssA,'main') #get the main peak from hncaco A (looking forwards)
                pkBB=mol.GetPeakTp('hnco',pkAssB,'')       #get the main peak from hnco B (looking backwards)
                
                if(pkAF!=False):
                    print('from A, A->B (f)',pkAF.name,pkAF.tp,pkAF.f3p)
                    x,z,zN=mol.spec['hnco'].GetIntensity(pkAssB,pkAF.f3p) #A is fine: look back from HNCO for B for A.
                    if(numpy.fabs(zN)>1.2):
                        print('Keeping found peak:',pkAF.f3p,x,z,zN)
                        stry='add,hnco,%s,%.3f,%.2f # missing connection b from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MINOR' % (pkAssB,pkAF.f3p,z,pkAssB,pkAssA,'hnco',zN)
                        stryLines.append(stry)
                        if(pkBB!=False):
                            stry='remove,hnco,%s,%s # missing connection b from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MINOR' % (pkAssB,pkBB.name,pkAssB,pkAssA,'hnco',zN)
                            stryLines.append(stry)
                        return

                else:
                    print('cannot find A->B from A')
                
                #or could take the HNCA minor from A if there are two peaks there.

                if(pkBB!=False):
                    print('from B, B->A (b)',pkBB.name,pkBB.tp,pkBB.f3p)  
                    x,z,zN=mol.spec['hncaco'].GetIntensity(pkAssA,pkBB.f3p)
                    if(zN>1.2):
                        print('Keeping found peak:',pkBB.f3p,x,z,zN)
                        stry='add,hncaco,%s,%.3f,%.2f # missing connection f from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MAIN' % (pkAssA,pkBB.f3p,z,pkAssA,pkAssB,'hncaco',zN)
                        stryLines.append(stry)  #check intensity is sufficient to make it max! this should be the new main.
                        if(pkAF!=False):
                            stry='remove,hncaco,%s,%s # missing connection f from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MAIN' % (pkAssA,pkAF.name,pkAssA,pkAssB,'hncaco',zN)
                            stryLines.append(stry)  #check intensity is sufficient to make it max! this should be the new main.
                        return

                else:
                    print('cannot find B->A from B')

                """
                #assume either 'main' in HNCACO from A looking forward or 'main' in HNCO in B looking backwards
                if(pkBB==False and pkAF!=False): #if one other is false. #No A in HNCO of B
                    x,z,zN=mol.spec['hnco'].GetIntensity(pkAssB,pkAF.f3p) #A is fine: look back from HNCO for B for A.
                    print('Found peak!',pkAF.f3p,x,z,zN)
                    if(zN>1.2):
                        print('keeping.')
                        stry='add,hnco,%s,%.3f,%.2f # missing connection b from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MINOR' % (pkAssB,pkAF.f3p,z,pkAssB,pkAssA,'hnco',zN)
                        stryLines.append(stry)
                        return
                if(pkAF==False and pkBB!=False):  #B is fine: look in HNCA of A for B
                    x,z,zN=mol.spec['hncaco'].GetIntensity(pkAssA,pkBB.f3p)
                    print('Found peak!',pkBB.f3p,x,z,zN)
                    if(zN>1.2):
                        print('keeping.')
                        stry='add,hncaco,%s,%.3f,%.2f # missing connection f from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MAIN' % (pkAssA,pkBB.f3p,z,pkAssA,pkAssB,'hncaco',zN)
                        stryLines.append(stry)  #check intensity is sufficient to make it max! this should be the new main.
                        return
                """
                print('No connnection found using current classifications')

                print('A')
                mol.ShowPeaks('hnco',pkAssA)
                mol.ShowPeaks('hncaco',pkAssA)
                print('B')
                mol.ShowPeaks('hnco',pkAssB)
                mol.ShowPeaks('hncaco',pkAssB)
                return
                #check the forward view.   
                if(pkBB!=False): # we have backward hncoca main.
                    #loop over HNCA peaks from A, do any match?
                    vals=[]
                    if('hnca' in mol.peak[pkAssA].keys()):
                        for pk3 in mol.peak[pkAssA]['hnca']:
                            vals.append(numpy.fabs(pk3.f3p-pkBB.f3p))
                    argy=numpy.argmin(vals)

                    print('Miss-classification of main in hnca of %s. Adjusting.' % pkAssA)
                    
                    stry='set,hnca,%s,main' % (mol.peak[pkAssA]['hnca'][argy].name)
                    
                    stryLines.append(stry)
                    
                    stry='remove,hnca,%s,%s' % (pkAssA,pkAF.name)
                    stryLines.append(stry)
                    #print(stry)

                    #check intensity 

                    #pkBB2=GetPeak('hnca',pkAssB,'') #get second peak from hnca
                    #if(numpy.fabs(pkBB2.f3p-pkBB.f3p)<1): #if we have a reasonable match between the two CA peaks.

                
                #check the backward view.
                #print (sc)

                #if('hnca' in mol.peak[pk].keys()):
                #    for pk3 in mol.peak[pk][sc]:
                #        print(sc,pk3.name,pk3.tp,pk3.f3p,pk3.inty)


        def DoCAcheck(pkAssA,pkAssB,stryLines):
                #FIRST, assume the main/minor classifications are correct.
                pkAF=mol.GetPeakTp('hnca',pkAssA,'main') #get the main peak from hnca A (looking forwards)
                pkBB=mol.GetPeakTp('hncoca',pkAssB,'') #get the main peak from hncoca main (looking backwards)

                #if we have an HNCA main looking forwards....
                if(pkAF!=False):
                    print('from A, A->B (f)',pkAF.name,pkAF.tp,pkAF.f3p) #look for evidence of point from B->A
                    x,z,zN=mol.spec['hnca'].GetIntensity(pkAssB,pkAF.f3p) #A is fine: look in HNCA for B for A.
                    if(numpy.fabs(zN)>2):
                        print('Keeping found peak:',pkAF.f3p,x,z,zN)
                        stry='add,hncoca,%s,%.3f,%.2f # missing connection b from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MINOR' % (pkAssB,pkAF.f3p,z,pkAssB,pkAssA,'hnca',zN/100)
                        stryLines.append(stry)
                        stry='swap,hnca,%s,%s,%s   # taking most intense as main as matches TOCSY' % (pk,pk3.name,pk4.name)
                        #if there is one, remove the existing backward HNCA peak.
                        if(pkBB!=False):
                            stry='remove,hncoca,%s,%s # existing connection is bad. remove.' % (pkAssB,pkBB.name)
                            stryLines.append(stry)
                            

                        return
                        #print (stry)
                else:
                    print('cannot find A->B from A')
                
                #if we have an HNCOCA looking backwards...
                if(pkBB!=False):
                    print('from B, B->A (b)',pkBB.name,pkBB.tp,pkBB.f3p)  
                    x,z,zN=mol.spec['hnca'].GetIntensity(pkAssA,pkBB.f3p)  #look for evidence of point from A->B
                    
                    if(numpy.fabs(zN)>2):
                        print('Keeping found peak:',pkBB.name,pkBB.tp,pkBB.f3p,x,z,zN)
                        stry='add,hnca,%s,%.3f,%.2f # missing connection f from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MAIN' % (pkAssA,pkBB.f3p,z,pkAssA,pkAssB,'hncoca',zN*100)
                        stryLines.append(stry)  #check intensity is sufficient to make it max! this should be the new main.
                        if(pkAF!=False):
                            stry='remove,hnca,%s,%s # existing connection is bad. remove.' % (pkAssA,pkAF.name)
                            stryLines.append(stry)

                        #print(stry)
                        return
                else:
                    print('cannot find B->A from B')

                """
                #assume either 'main' in HNCA from A looking forward or 'main' in HNCOCA in B looking backwards
                if(pkBB==False and pkAF!=False): #if one other is false.
                    x,z,zN=mol.spec['hnca'].GetIntensity(pkAssB,pkAF.f3p) #A is fine: look in HNCA for B for A.
                    print('Found peak!',pkAF.f3p,x,z,zN)
                    if(zN>2):
                        print('keeping.')
                        stry='add,hnca,%s,%.3f,%.2f # missing connection b from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MINOR' % (pkAssB,pkAF.f3p,z,pkAssB,pkAssA,'hnca',zN)
                        stryLines.append(stry)
                        return
                if(pkAF==False and pkBB!=False):  #B is fine: look in HNCA of A for B
                    x,z,zN=mol.spec['hnca'].GetIntensity(pkAssA,pkBB.f3p)
                    print('Found peak!',pkBB.f3p,x,z,zN)
                    if(zN>2):
                        print('keeping.')
                        stry='add,hncoca,%s,%.3f,%.2f # missing connection f from peak %s to %s in %s. Found peak with s/n %.2f. NEEDS TO BE MAIN' % (pkAssA,pkBB.f3p,z,pkAssA,pkAssB,'hncoca',zN)
                        stryLines.append(stry)  #check intensity is sufficient to make it max! this should be the new main.
                        return
                print('No connnection found using current classifications')
                
                #check the forward view.   
                if(pkBB!=False and pkAF!=False): # we have backward hncoca main.
                    #loop over HNCA peaks from A, do any match?
                    vals=[]
                    if('hnca' in mol.peak[pkAssA].keys()):
                        for pk3 in mol.peak[pkAssA]['hnca']:
                            vals.append(numpy.fabs(pk3.f3p-pkBB.f3p))
                    if(len(vals)>0):
                        argy=numpy.argmin(vals)

                        print('Miss-classification of main in hnca of %s. Adjusting.' % pkAssA)
                        
                        stry='set,hnca,%s,main' % (mol.peak[pkAssA]['hnca'][argy].name)
                        
                        stryLines.append(stry)
                        
                        stry='remove,hnca,%s,%s' % (pkAssA,pkAF.name)
                        stryLines.append(stry)
                        #print(stry)

                        #check intensity 

                        #pkBB2=GetPeak('hnca',pkAssB,'') #get second peak from hnca
                        #if(numpy.fabs(pkBB2.f3p-pkBB.f3p)<1): #if we have a reasonable match between the two CA peaks.
                """
                print("no data for possible connection found. peaks are:")
                print('A')
                mol.ShowPeaks('hnca',pkAssA)
                mol.ShowPeaks('hncoca',pkAssA)
                print('B')
                mol.ShowPeaks('hnca',pkAssB)
                mol.ShowPeaks('hncoca',pkAssB)

                #check the backward view.
                #print (sc)

                #if('hnca' in mol.peak[pk].keys()):
                #    for pk3 in mol.peak[pk][sc]:
                #        print(sc,pk3.name,pk3.tp,pk3.f3p,pk3.inty)


        def ExamineConnection(resi,resiN,noes,stryLines):  
            
            pkAssA=self.GetKeyFromVal(str(resi)+mol.seq[resi],assNew)
            pkAssB=self.GetKeyFromVal(str(resiN)+mol.seq[resiN],assNew)  
            print()
            print('Examining connection between ',pkAssA,str(resi)+mol.seq[resi],'and',pkAssB,str(resiN)+mol.seq[resiN]) 
            segA=self.GetLinks(pkAssA,'f',noes)
            segB=self.GetLinks(pkAssB,'b',noes)
            #print(segA)
            #print(segB)
            if(pkAssB in segA):
                print('Connection seen in NOEs A->B')
            else:
                #print(mol.G1edges[pkAssA])
                print('Connection not seen in NOEs A->B. Searching origianl list:')
                for pk3,val,dir in mol.G1edges[pkAssA]: #first, look in the G1edges. has this been removed by mistake?
                     #print (pk3,val,dir)
                     if(dir=='f'):
                        if(pk3==pkAssB):
                            print('Restoring connection:',pkAssA,pk3,val,dir)
                            #need to restore this impact.
                            #print(noes[pkAssB])
                            noes[pkAssA].append((pk3,val,dir))
                            return
                print('No connection seen in original graph for A->B')


            if(pkAssA in segB):
                print('Connection seen in NOEs B->A')
            else:
                print('Connection NOT seen in NOEs B->A. Searching original list:')
                for pk3,val,dir in mol.G1edges[pkAssB]: #first, look in the G1edges. has this been removed by mistake?
                     #print (pk3,val,dir)
                     if(dir=='b'):
                        if(pk3==pkAssA):
                            print('Restoring connection:',pkAssB,pk3,val,dir)
                            #need to restore this impact.
                            #print(noes[pkAssB])
                            noes[pkAssB].append((pk3,val,dir))
                            return
                print('No connection seen in original graph for B->A')
            

            check=CheckConns(pkAssA,pkAssB) #check the original classifications, return spectra where we lack a connection.
            print ('From original classificaitons, missing connections in:',check)
            if('CA' in check):
                print('Looking for CA connections')
                print('A')
                #ShowPeak('hnca',pkAssA)
                #ShowPeak('hncoca',pkAssA)

                #the i-1 is fine for A
                #the i+1 is fine for B
                #can we see a + A->B from A
                #can we see a - B->A from B
                DoCAcheck(pkAssA,pkAssB,stryLines)
    
                    
                    
                #print('B')
                #ShowPeak('hnca',pkAssB)
                #ShowPeak('hncoca',pkAssB)
                #print(pkAssA,pkAssB)
            
            if('CO' in check):
            
                DoCOcheck(pkAssA,pkAssB,stryLines)
                

            
            #mol.peaks[pkAssF]['hnca']
            #mol.peaks[pkAssB]['hnco']
            #look at HNCA

            #look at HNCO
            #look at HNCACO



        linky={}
        bad=[]
        stryLines=[]
        RUN=False
        for i,peak in enumerate(peaks): #get current state of play.
            stry,prob,great=self.DoReport(i,peak,peaks,atoms,assNew,noes,assOpt,assGreat=assGreat,great='report',verb=False)
            print(peak,assNew[peak],stry)
            if(stry=='BAD'):
                    bad.append(peak)
            linky[int(assNew[peak][:-1])]=stry
        
        for resi,stry in linky.items():
            if(stry=='BACKWARD GREAT!'):
                if(resi+1 in linky.keys()):
                    ExamineConnection(resi,resi+1,noes,stryLines)
            if(stry=='FORWARD GREAT!'):
                if(resi+1 in linky.keys()):
                    if(linky[resi+1]=='FORWARD GREAT!'): 
                        ExamineConnection(resi,resi+1,noes,stryLines)
                    
        #print(noes['58'])                

        print('Manual lines to add:')                   
        for stry in stryLines:
            print(stry)
        #assOpt,InvAss=self.PurgeGreat(assOpt,assGreat)
        
        #return assOpt,InvAss,assGreat,bad,RUN



            #assign, look for confident assignments, elongate the ends
            #clean up the NOE graph.
    def Cycle(self,peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=[],permissive=False):
                peaks,assNew,cost=self.OptimiseCost(peaks,atoms,noes,assOpt,InvAss,assNew,guess,scale=0.1,fix=fix,depth=1)
                #print('bb',assNew['112'],self.assStatus['112'])
                #print('bb',assNew['103'],self.assStatus['103'])

                assGreat,RUN=self.DoGreatClean(peaks,atoms,assNew,noes,assOpt)  #extend if possible
                assOpt,InvAss=self.PurgeGreat(assOpt,assGreat)
                noes=self.CheckGraphs(noes,assNew)
                if(RUN): #repeat if there has  been a change
                    peaks,assNew,cost=self.OptimiseCost(peaks,atoms,noes,assOpt,InvAss,assNew,guess,scale=0.1,fix=fix,depth=1)

                assOpt,InvAss,assGreat,bad,RUN=self.DoCommonConnect(peaks,atoms,noes,assNew,assOpt,InvAss,assGreat,permissive=permissive)
                #print('eb',assNew['112'],self.assStatus['112'])
                #print('eb',assNew['103'],self.assStatus['103'])
                return peaks,assNew,assGreat,cost,noes,assOpt,InvAss,bad
            

            #repeat optimisation/clean cycles until complete
    def RunCycles(self,peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=[],permissive=False):
                curr=0,0
                cnt=0
                while(1==1):
                    cnt+=1
                    print("Cycle starting:",cnt)
                    peaks,assNew,assGreat,cost,noes,assOpt,InvAss,bad=self.Cycle(peaks,atoms,noes,assOpt,InvAss,assNew,guess,fix=fix,permissive=permissive)
                    
                  
                    

                    graphs,multi=self.CountGraphs(noes)
                    print('Cycle complete:',cnt,'Conf:',len(assGreat.keys()),'Graphs:',graphs,'Bad:',len(bad))
                    if(len(assGreat.keys())==curr[0] and graphs==curr[1]):
                        break
                    else:
                        curr=len(assGreat.keys()),graphs
                return peaks,assNew,assGreat,cost,noes,assOpt,InvAss,bad




    def TestSanity(self,peak,atom,dir,i,j,atoms,peaks,noes,assNew,assOpt,verb=True,depth=True):
                    if(verb):
                        print (peak,atom,atoms[j]) #checking assignment.
                    resi=int(atom[:-1])
                    if(dir=='f'): 
                        #this is the C terminus and there is no i+1
                        if(j==len(atoms)-1):
                            #does the assignment have someone strongly connected to it?
                            seg=self.GetLinks(peak,dir,noes)
                            if(len(seg)==0):
                                return False,True
                            else:
                                return False,False

                        nextAtom=atoms[j+1]  #this is the next assignment.
                        resiN=int(nextAtom[:-1])
                        stry='(forwards)'
                    if(dir=='b'):
                        #this is the N terminus, and there's no i-1.
                        if(j==0): 
                            #does the assignment have someone strongly connected backwards?
                            seg=self.GetLinks(peak,dir,noes)
                            if(len(seg)==0):
                                return False,True
                            else:
                                return False,False
                        nextAtom=atoms[j-1]  #this is the next assignment.
                        resiN=int(nextAtom[:-1])
                        stry='(backwards)'

                        mol=self.parent.parent.tabOne.molecule
                        if(mol.seq[resi-1]=='P'): #do a test if P is i-1
                            if(verb):
                                print('Next atom ',stry,' is a proline:' ,nextAtom)

                            #do a TOCSY test if we are going backwards, is i-1 a P?
                            #resnT=mol.shiftx2[shiftxresi-1]['resn']
                            #map TOCSY peaks onto this residue and classify.
                            res=mol.DoScrTOCSY('C',peak,mol.bmrbC['PRO'],cerr=3.,resn='P')
                            prob=1.0
                            cnt=0
                            for tp,vals in res.items():
                                if('X' not in tp):
                                    expt,bmrb,std,probT,aT=vals
                                    prob*=probT
                                    cnt+=1
                            if(cnt>2 and prob>0.5):
                                 return False,True
                            else:
                                 return False,False
                           

                    if(numpy.abs(resiN-resi)!=1):
                            if(verb):
                                print("Jump in atoms: no connection expected.")
                                print(atom,nextAtom)
                            seg=self.GetLinks(peak,dir,noes) #get for connections in current direction...
                            if(len(seg)==0): #no connections in this direction. good match.
                                return False,True
                            else:
                                return False,False
   
                    if(verb):
                        print('next atom ',stry,':',nextAtom)
                    
                    pkAss=self.GetKeyFromVal(nextAtom,assNew)
                    
                    if(pkAss==False):
                        if(verb):
                            print('No assignment placed for next atom ',stry)
                        return False,False
                    if(verb):
                        print("Assigned to:",pkAss)

                    seg=self.GetLinks(peak,dir,noes)
                    
                    t=0
                    for s in seg:
                        if(depth): #go loooking for new assignments...
                            for (pk2,val,d) in noes[peak]: #for all connections to peak...
                                if(d==dir): #in the right direction...
                                    if(verb):                   
                                        print('found ',stry,'link:',s,(pk2,val,d))
                                    #print("  ",pk2,"is assigned to:",assNew[pk2])
                                    for (pk3,vol) in assOpt[nextAtom]: #look for atom in possibilities
                                        if(pk3==pk2): #if the link has a possibility of being placed as the next atom, we're good!
                                            #print("   probability to be placed here:",vol)
                                            #print("   current status of placement:")
                                            #print('bb',assNew.keys())
                                            if(pk2 not in assNew.keys()):
                                                return pkAss,False
                                            jnew=numpy.where(assNew[pk2]==atoms)[0][0]
                                            inew=numpy.where(peaks==pk2)[0][0]
                                            pkAssF,successF=self.TestSanity(pk2,assNew[pk2],'f',inew,jnew,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                                            pkAssB,successB=self.TestSanity(pk2,assNew[pk2],'b',inew,jnew,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                                            if(successF==False and successB==False):
                                                if(verb):
                                                    print ("   NEW ASSIGNMENT:",pk2,nextAtom)
                                                return (pk2,nextAtom,peak,assNew[pk2]),True
                                                #print(assOpt[nextAtom])
                        if(pkAss==s):  #the peak assigned to the next atom has a connection to this one 
                            t=1
                    
                    if(t==1):   #great!
                        return pkAss,True

                    #trawl though options to find nearest connections
                    for key,vals in self.parent.parent.tabOne.molecule.Optedges[peak].items():
                        for v in vals:
                            if(v[0]==pkAss and v[2]==dir):
                                if(verb):
                                    print('Closest connection:',key,v)
                    return pkAss,False
                    
                    #if(t==0 and i<len(peaks)-1):
                    #    print('We should see a f link between',peak,'and',peaks[i+1])
                    #    for key,vals in mol.Optedges[peak].items():
                    #        for v in vals:
                    #            if(v[0]==peaks[i+1] and v[2]=='f'):
                    #                print(key,v)


    #look in NOEs for the peak. 
    #Return matches that hit the specified direction
    def GetLinks(self,peak,d,noes):
                seg=[]
                if(peak in noes.keys()):
                    for (noe,val,dir) in noes[peak]:
                        if(dir==d):
                            seg.append(noe)
                return seg
    def GetMatches(self,jj,seg,atoms,assOpt):
                voptN=[]
                voptS=[]
                for f in seg:
                    if(atoms[jj] in assOpt.keys()): 
                        for v in assOpt[atoms[jj]]:
                            if(v[0]==f):
                                voptN.append(v[0])
                                voptS.append(v[1])
                voptS=numpy.array(voptS)
                return voptN,voptS
            
    #for a given peak, first get new peaks in a given direction (d, forward or backward)
    #then look to see if the correct assignment (jj+1) is there.
    #return the score if so.
    #can also go to depth2.
    def NextProb(self,peak,d,noes,atoms,j,jj,assOpt,scale,depth):
                                
                seg=self.GetLinks(peak,d,noes)
                #2. for each of those, work out the probabilty if the atom is +1
                if(numpy.fabs(int(atoms[j][:-1])-int(atoms[jj][:-1]))!=1.0): #make sure that we are looking for +/-1
                    return scale,False

                voptN,voptS=self.GetMatches(jj,seg,atoms,assOpt)
                if(len(voptS)>0):
                    if(depth==1):
                        arg=numpy.argmax(voptS)
                        return voptS[arg],True
                    """
                    vscores=[]
                    atomID=atoms[jj]
                    if(d=='f'):
                        if(jj>=len(atoms)):
                            arg=numpy.argmax(voptS)
                            return voptS[arg],True

                        atomID_next=atoms[jj+1]

                        if(int(atomID[:-1])+1!=int(atomID_next[:-1]) or atomID_next[-1]=='P'):
                            arg=numpy.argmax(voptS)
                            return voptS[arg],True

                        for i in range(len(voptN)):
                            seg2=self.GetLinks(voptN[i],d,noes)
                            #print(seg,seg2)
                            voptN2,voptS2=self.GetMatches(jj+1,seg2,atoms,assOpt)
                            if(len(voptS2)!=0): 
                                #print('ye')
                                vtmp=voptS2*voptS[i]  
                                vscores=numpy.concatenate((vscores,vtmp))
                            else:
                                vscores=(voptS*scale)
                        vscores=numpy.array(vscores)
                        #print( vscores)
                        return numpy.max(vscores),True
                    
                    elif(d=='b'):
                        if(jj==0):
                            arg=numpy.argmax(voptS)
                            return voptS[arg],True
                        
                        atomID_next=atoms[jj-1]
                        if(int(atomID[:-1])-1!=int(atomID_next[:-1]) ):
                            arg=numpy.argmax(voptS)
                            return voptS[arg],True
                        for i in range(len(voptN)):
                            seg2=self.GetLinks(voptN[i],d,noes)
                            voptN2,voptS2=self.GetMatches(jj+1,seg2,atoms,assOpt)  
                            if(len(voptS2)!=0): 
                                #print('ye')
                                vtmp=voptS2*voptS[i]
                                vscores=numpy.concatenate((vscores,vtmp))
            
                            else:
                                vscores=voptS*scale
                        vscores=numpy.array(vscores)
                        return numpy.max(vscores),True
                    """

                    #if(len(voptS)>1):
                    #    print('best match:',peak,voptS[arg],atoms[jj],voptN[arg],d,voptN)
                return scale,False #return rescaled probability

    def AdjustProbability(self,peak,j,prob,noes,atoms,atomMax,assOpt,scale,depth):
                #1. get the forward NOE links

                
                if(j!=atomMax-1):
                    pf,linkf=self.NextProb(peak,'f',noes,atoms,j,j+1,assOpt,scale,depth) #try to find a good match forwards
                else:
                    linkf=True;pf=1

                #2. get the backward NOE links
                if(j!=0):
                    #if one backwards is a proline, then don't penalise in the cost. 
                    if(self.parent.parent.tabOne.molecule.seq[int(atoms[j][:-1])-1]=='P'):
                        linkb=True;pb=1
                    else:
                        pb,linkb=self.NextProb(peak,'b',noes,atoms,j,j-1,assOpt,scale,depth)  #try to find a good match backwards
                else:
                    linkb=True;pb=1

                #if(peak=='158' and atoms[j]=='158S'):
                #    print('     ',linkb,linkf,prob,pf,pb)
                if(linkb==False and linkf==False):
                    return 0,0
                else:
                    cnt=0
                    if(linkb!=False):
                        cnt+=1
                    if(linkf!=False):
                        cnt+=1
                    #give no penalty if detect zero forward or backward links
                    #if(linkf==False):
                    #    pf=1.0
                    #if(linkb==False):
                    #    pb=1.0
                    return prob*pf*pb,cnt
                #return prob                    

    #get the slice from a spectrum.
    #run the maximising peak picker with current settings (width/StoN)
    #run chemical shift mapping to see if we can map extracted to expected
    #compare to current assignment dict to see if we have anyone new.
    def MatchTOCSYslice(self,peak,spec,resn,nuc,vol,vop,ti,width=1.0,h=1.0,res={},stryLine=[]):
            
                mol=self.parent.parent.tabOne.molecule

                #X,Z=mol.spec[spec].GetSlice(peak)   #extract slice from raw data
                #dx=numpy.fabs(X[1]-X[0])  #work out spacing.
                #dist=width/dx #1ppm divided by dwell space is minimum distance
                #import scipy
                #pks,prop=scipy.signal.find_peaks(Z,distance=dist,width=(None,None),height=mol.spec[spec].noise*h)
                
                X,Z=mol.spec[spec].FindPeakMaxima(peak,width,h)

                
                #self,valR,vol,vop,ti,nuc):
                ras,cost=mol.MatchAtoB(X,vol,vop,ti,nuc)  #map picked peaks onto expected assignments
                if(peak=='28'):
                    print (ti)
                    print(X,Z)
                    print()
                    print(res)
                    print (ras)
                    #analyse the output and return required values
                ros={}
                found=[]
                for ass,vals in ras.items():
                    if(len(vals)==1):  #not assigned.
                         continue
                    if('X' in ass): #not assigned 
                         continue
                    if(ass not in res.keys()):  #assigned and not in reference list.
                        ros[ass]=vals
                        found.append(ass)
                        #print(vals)
                        a=int(vals[4])
                        #valR[a],vol[b],vop[b],prob,a=vals #unpack
                        stry='add,%s,%s,%.2f,%.2f    #  peaks detected in %s! %s %s s/n: %.2f' % (spec,peak,X[a],Z[a],spec,resn,ass,Z[a]/mol.spec[spec].noise)
                        stryLine.append(stry) 

                

                return found,cost,ros


    def CategoriseRes(self,peak,spec,res,ti):
            seen=[]
            excess=[]
            prob=1.0
            cnt=0

            #reset labels.
            if(spec in self.parent.parent.tabOne.molecule.peak[peak].keys()):
                for pk3 in self.parent.parent.tabOne.molecule.peak[peak][spec]:
                    if('CA' not in pk3.tp):
                        pk3.tp=''


            for tp,vals in res.items():
                if('X' in tp): #anybody unassigned? overpick?
                    excess.append(tp)
                    continue
                expt,bmrb,std,probT,aT=vals #unpack.
                
                prob*=probT
                cnt+=1
                seen.append(tp)
                #categorise.
                self.parent.parent.tabOne.molecule.peak[peak][spec][vals[4]].tp=tp+'(i-1)' #adjust label.
      
            rem=[]
            for typ in ti:
                if(typ not in seen): #exclude residues that can't be seen
                    rem.append(typ) 

            return seen,excess,rem,prob



    def DoTOCSYscore(self,nuc,peak,atom,spec,stryLine=[]):
        mol=self.parent.parent.tabOne.molecule
        #spectrum, peak name, specify target ppm and return intensity
        


                
                
        """
        if(ti[b] in res.keys()):
                    continue
            
            found.append(ti[b])
            ros[ti[b]]=X[pks[a]],vol[b],vop[b],prob,Z[pks[a]]
        return found,cost,ros
        """
        
        """
        argy=numpy.flip(numpy.argsort(numpy.fabs(Z[pks]))) #sort intensity big to small.
        vol,vop,ti=mol.ConsolidateShifts(spec,nuc,shifts,cerr)  #get expected intensities from reference.

        XX,YY=numpy.meshgrid(vol,X[pks]) #create grids of constant row = bmrb and col = expt
        diff=numpy.fabs(XX-YY)**2.       #get the absolute value of their differences
        from scipy.optimize import linear_sum_assignment
        ass=linear_sum_assignment(diff)  #get most probable assignment
        found=[]
        cost=0
        for a,b in zip(ass[0],ass[1]):
            cost+=diff[a,b]**0.5
            prob=numpy.exp(-1.*(diff[a,b])/(2*vop[b]))
            if(ti[b] in res.keys()):
                    continue
            stry='add,%s,%s,%.2f,%.2f    #  peaks detected in %s! %s %s s/n: %.2f' % (spec,peak,X[pks[a]],Z[pks[a]],spec,mol.seq[int(atom[:-1])-1],ti[b],Z[pks[a]]/mol.spec[spec].noise)
            stryLine.append(stry) 
            found.append(ti[b])
            ros[ti[b]]=X[pks[a]],vol[b],vop[b],prob,Z[pks[a]]
        return found,cost,ros
        """
        
   
        ######################################################################
        
        if(nuc=='C'):
            cerr=3
            bmrb=mol.bmrbC   
            width=1.0  
            merge=1.5
        elif(nuc=='H'):
            cerr=0.1
            bmrb=mol.bmrbH
            width=0.1
            merge=0.01
            
            
        if( (int(atom[:-1])-1) in mol.shiftx2.keys()):
            shifts=mol.shiftx2[int(atom[:-1])-1][nuc]
        elif( mol.seq[int(atom[:-1])-1]=='P'):
            shifts=bmrb['PRO']
        else:
            return -1,[],[],[],[],0

        resn=mol.seq[int(atom[:-1])-1]  #get residue type of i-1

        
        valR=mol.GetPeaksF3p(peak,spec) #get unidec detected peaks #(not 'A')
        vol,vop,ti=mol.ConsolidateShifts(spec,nuc,shifts,cerr=cerr,merge=merge)  #get peaks from reference
        res,cost=mol.MatchAtoB(valR,vol,vop,ti,nuc)     #map one onto the other
        seen,excess,rem,prob=self.CategoriseRes(peak,spec,res,ti) #categorise results dict (seen, excees and probability)
        

        #if(len(valR)>0):  
        #    res=mol.DoMatchTOCSY(nuc,valR,resn,bmrb[mol.p1to3[resn]],mol.spec[spec].ref,spec,CA=True,SCR=False,cerr=cerr,merge=merge,STRICT=False) #so strict point to point comparison
        #    seen,excess,prob=CategoriseRes(spec,res) #categorise results dict (seen, excees and probability)
        #    if(peak=='28'):
        #            print (excess)
        #            print (res)
        #else:
        #    seen=[];excess=[];prob=1.0;res={}
        #if we have some excess, try going again without the merger.
        #if(len(excess)>1):
        #    merge=False
        #    res=mol.DoScrTOCSY(nuc,peak,shifts,cerr=cerr,resn=resn,merge=merge,IncCA=True) #assign peaks using unidec assignments
        #    seen,excess,prob=CategoriseRes(spec,res) #categorise results dict (seen, excees and probability)
        #    if(peak=='28'):
        #         print (excess)

        #vol,vop,ti=mol.ConsolidateShifts(nuc,shifts,cerr,merge=merge,CA=True) #get shifts, group within threshold
        
        hbottom=0.8  #value to bottom out at
        hinc=0.95    #geometric factor to reduce height by
        h=1.0
        stryTmp=[];ros={};found=[];remNew=[] #initialise
        while(len(rem)>0):
            stryTmp=[]

            found,cost,ros=self.MatchTOCSYslice(peak,spec,resn,nuc,vol,vop,ti,width=width,h=h,res=res,stryLine=stryTmp)
            remNew=[]
            for r in rem: #if we have any remainders...
                if(r not in found):
                    remNew.append(r)

            if(len(remNew)>0): #reduce the s/n threshold
                h*=hinc        #and go again.
            else:
                break
            if(h<hbottom): #if we reach rock bottom, abort.
                break
        
        #if from the final pass we have some peaks to add, add to the list!
        if(len(stryTmp)>0): #add the new lines into the mix.
            for stry in stryTmp:
                stryLine.append(stry)
                #add the peak.
        
        #amend the peak list and set typ
        #Add peak to molecule peak list at f3p, extra letter 'ref' added to the end
        #first unused number taken for indexing.
        #tp set to specified value. x and y taken from peak2D.
        for typ,vals in ros.items():
            if(len(vals)!=5): #classify.
                continue
            mol.AddNewPeak(peak,spec,vals[0],typ+'(i-1)','A')

        #transfer all peaks in seen to found if they have a peak name that ends 'A'
        
        seenNew=[]
        for s in seen:
            tick=0
            for pk3 in mol.peak[peak][spec]:
                if(pk3.tp.split('(i-1)')[0]==s):
                    if(pk3.name[-1]=='A'):
                        found.append(s)
                        tick=1
                        break
            if(tick==0):
                seenNew.append(s)


        return prob,seenNew,excess,found,remNew,cost


    def DoReport(self,i,peak,peaks,atoms,assNew,noes,assOpt,assGreat={},great=True,verb=True,stryLine=[],latex=False):
                #print()
                #if(peak not in assNew.keys()):
                #    print("no entry for:",peak)
                #    return
                mol=self.parent.parent.tabOne.molecule
                def FormatReportStrSide(nuc,atom,seen,found,rem,expect,cost,latex=False):
                        ss=''
                        for s in seen:
                            ss+=' '+s
                        
                        if(len(found)>0):
                            ss+='/' 
                            for f in found:
                                ss+=' '+f
                        
                        if(len(rem)>0):
                            ss+=' ['
                            for r in rem:
                                ss+=' '+r
                        ee=''
                        if(len(expect)>0):
                            ee+=' %i' % len(expect)
                        if(nuc=='C'):
                            stry='%1s %22s %6.2f %3s' % (mol.seq[int(atom[:-1])-1],ss,cost,ee)
                            if(latex!=False):
                                latex+='& %1s & %22s & %6.2f & %3s' % (mol.seq[int(atom[:-1])-1],ss,cost,ee)
                        elif(nuc=='H'):
                            stry='%1s %34s %6.2f %3s' % (mol.seq[int(atom[:-1])-1],ss,cost,ee)
                            if(latex!=False):
                                latex+='& %1s & %34s & %6.2f & %3s' % (mol.seq[int(atom[:-1])-1],ss,cost,ee)
                        return stry,latex
                atom=assNew[peak]  
                j=numpy.where(atoms==atom)[0][0]
                
                if(great==True):
                    pkAssF,successF=self.TestSanity(peak,atom,'f',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    pkAssB,successB=self.TestSanity(peak,atom,'b',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    if(successF and successB):
                        assGreat[peak]=atom
                elif(great=='extend'):
                    #print('extend!')
                    pkAssF,successF=self.TestSanity(peak,atom,'f',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    pkAssB,successB=self.TestSanity(peak,atom,'b',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    if(successF and successB):
                        return False
                    
                    #if one side is great, then add more if there are good options going that way.
                    if(successF ==False):
                        dir='f'
                    elif(successB==False):
                        dir='b'
                    res2,successF2=self.TestSanity(peak,atom,dir,i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=True)

                    if(successF2):
                        a,b,c,d=res2
                        if(a not in assGreat.keys()):
                            print("ADDING THROUGH EXTENSION:",a,b)
                            assGreat[a]=b

                            self.UpdateAss(assNew,a,b,'e')
                            self.UpdateAss(assNew,c,d,'e')
                            #print(self.assStatus['112'])
            
                            #assNew[a]=b
                            #assNew[c]=d
                            #self.assStatus[a]='e'
                            #fself.assStatus[c]='e'
                            #del assNew[c]
                            return True
                    return False

                elif(great=='report'):
                    pkAssF,successF=self.TestSanity(peak,atom,'f',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    pkAssB,successB=self.TestSanity(peak,atom,'b',i,j,atoms,peaks,noes,assNew,assOpt,verb=False,depth=False)
                    stry=''
                    if(successF and successB):
                        stry+='GREAT!'
                    elif(successF):
                        stry+='FORWARD GREAT!'
                    elif(successB):
                        stry+='BACKWARD GREAT!'
                    else: 
                        stry+='BAD'

                    ####GET PROB#####
                    InvAss=self.MakeInvAss(assOpt)
                    #print(peak,InvAss[peak])
                    
                    #print(self.parent.parent.tabOne.molecule.assRef[atom])
                    tig=0
                    for pk,prob in self.parent.parent.tabOne.molecule.assRef[atom]:
                        if(pk==peak):
                            tig=1
                            break
                    if(tig==0):
                        prob=0
                    #try: 
                    #    prob=InvAss[peak][atom]
                    #except:
                    
                    liney='%5s %5s %s %15s %7.3f %8s ' % (peak,atom,self.assStatus[peak],stry,prob,peak in assGreat.keys())
                    
                    lineyTex=False
                    if(latex!=False):
                        lineyTex='%5s & %5s & %s & %15s & %7.3f & %8s ' % (peak,atom,self.assStatus[peak],stry,prob,peak in assGreat.keys())
                        
                    if('ctocsy' in self.parent.parent.tabOne.molecule.spec.keys()):
                        probC,seenC,expectC,foundC,remC,costC=self.DoTOCSYscore('C',peak,atom,'ctocsy',stryLine=stryLine) #get TOCSY probabilty
                        stryC,lineyTex=FormatReportStrSide('C',atom,seenC,foundC,remC,expectC,costC,latex=lineyTex)
                        liney+=' %35s' % (stryC)
                    if('hcconh' in self.parent.parent.tabOne.molecule.spec.keys()):    
                        probH,seenH,expectH,foundH,remH,costH=self.DoTOCSYscore('H',peak,atom,'hcconh',stryLine=stryLine) #get TOCSY probabilty
                        stryH,lineyTex=FormatReportStrSide('H',atom,seenH,foundH,remH,expectH,costH,latex=lineyTex)
                        liney+=' %47s' % (stryH)
                    #to handle the excess peaks:
                    #39, we see two value CGs, but shiftx predicts 1 in range.
                    #same for L, eg 48
                    if(latex!=False):
                        self.latexTable.append(lineyTex+' \\\\\n')
                    if(verb):
                        #%'%5s %5s %15s %7s %8s %19s %11s %5s %3i'


                        #print('%5s %5s %15s %7.3f %8s  %35s %47s' % (peak,atom,stry,probC,peak in assGreat.keys(),stryC,stryH))
                        print(liney)
                        #print(peak,atom,stry,prob,)
                        
                    return stry,prob,peak in assGreat.keys()
                
                #print("GREAT!")
                
                #print(pkAssF,pkAssB)
                #if(pkAssF!=False):
                #    if(pkAssB!=False):
                #        print("LOOK FOR MIDDLE")
                #        print(noes[pkAssB])
                #        print(noes[pkAssF])
                """
                    nextAtom=atoms[j-1]
                    if(nextAtom[-1]=='P'):
                        print('Next atom (backward) is a proline:' ,nextAtom)
                
                        print('next atom (backwards):',nextAtom)
                        seg=GetLinks(peak,'b',noes)
                        t=0
                        for s in seg:
                            if(peaks[i-1]==s):
                                t=1                        
                                print('found backwards link:',s,noes[peak])
                                break
                        if(t==0 and i>0):
                            print('We should see a b link between',peak,'and',peaks[i-1])
                            for key,vals in mol.Optedges[peak].items():
                                for v in vals:
                                    if(v[0]==peaks[i-1] and v[2]=='b'):
                                        print(key,v)
                """
                #print(mol.AnalPeak(peak))

                """
                mol.ShiftScore_i(peak)
                #print(mol.scr)
                print('consistency of shifts with i=',assNew[peak][-1])
                print(mol.scr[assNew[peak][-1]])
                

                #print(mol.scr)
                
                if(j!=0):
                    nextAtom=atoms[j-1]
                    print('consistency of shifts with i-1=',nextAtom[-1])
                    mol.ShiftScore_im1(peak)
                    try:
                        print(mol.scr[nextAtom[-1]])
                    except:
                        pass
                    print('expect 13C:',mol.bmrbC[mol.p1to3[nextAtom[-1]]])
                    print('expect 1H:',mol.bmrbH[mol.p1to3[nextAtom[-1]]])
                """



                

                #specord='hnco','hncaco','hnca','hncoca', 'hncacb', 'hncocacb', 'hncanh', 'hncocanh','ctocsy','hcconh'

                """
                res3=mol.p1to3[assNew[peak][-1]]

                print(' %8s %8s %7s %7s %7s ' % ('Nuc','Name','ppm','ppmBMRB','stdBMRB'))

                def ShowComparison(nuc,name,peaks,spec,bmrb,comp=False):
                    if(spec not in peaks.keys()):
                        return
                    peaks=peaks[spec]
                    for pk3 in peaks:
                        if(comp!=False):
                            if(pk3.tp!=comp):
                                continue
                        if(nuc=='H' and name=='H(i)'):
                            ppm=pk3.f1
                        elif(nuc=='N' and name=='N(i)'):
                            ppm=pk3.f2
                        else:
                            ppm=pk3.f3p
                        print(' %8s %8s %7.3f %7.3f %7.3f ' % (nuc,name,ppm,bmrb[0],bmrb[1]))

                ShowComparison('H','H(i)',mol.peak[peak],'hnco',GetBMRB('H','H',res3))
                ShowComparison('N','N(i)',mol.peak[peak],'hnco',GetBMRB('N','N',res3))
                ShowComparison('C','CO(i)',mol.peak[peak],'hncaco',mol.bmrbC[res3]['C'],comp='main')
                ShowComparison('C','CA(i)',mol.peak[peak],'hnca',mol.bmrbC[res3]['CA'],comp='main')

                if(j!=0):
                    res3=mol.p1to3[atoms[j-1][-1]]
                    ShowComparison('C','CO(i-1)',mol.peak[peak],'hnco',mol.bmrbC[res3]['C'])
                    ShowComparison('C','CA(i-1)',mol.peak[peak],'hnca',mol.bmrbC[res3]['CA'],comp='')

                    res=MatchTOCSY('C',res3,peak,mol.peak[peak])
                    
                    unass=False
                    for key ,vals in res.items():
                        if('X' not in key):
                            ShowComparison('C',key+'(i-1)',mol.peak[peak],'ctocsy',mol.bmrbC[res3][key],comp=key+'(i-1)')
                        else:
                            unass=True
                    if(unass):
                        print('unassigned:')
                        for key ,vals in res.items():
                            if('X' in key):
                                print(key,vals)
                                return
                    missingC=[]
                    for typ in mol.bmrbC[res3].keys():
                        if(typ not in res.keys() and typ!='CA' and typ!='C'):
                            if(mol.bmrbC[res3][typ][0]<110):
                                missingC.append(typ)
                    
                    res=MatchTOCSY('H',res3,peak,mol.peak[peak])
                    #print(res)
                    unass=False
                    for key ,vals in res.items():
                        if('X' not in key):
                            ShowComparison('H',key+'(i-1)',mol.peak[peak],'hcconh',mol.bmrbH[res3][key],comp=key+'(i-1)')
                        else:
                            unass=True
                    if(unass):
                        print('unassigned:')
                        for key ,vals in res.items():
                            if('X' in key):
                                print(key,vals)
                                #return
                    missingH=[]
                    for typ in mol.bmrbH[res3].keys():
                        if(typ not in res.keys() and typ!='H'):
                            missingH.append(typ)

                    #print(missing)
                    if(len(missingC)>0):
                        print('missing C peaks')
                        for miss in missingC:
                            print(miss,mol.bmrbC[res3][miss])

                    if(len(missingH)>0):
                        print('missing H peaks')
                        for miss in missingH:
                            print(miss,mol.bmrbH[res3][miss])

                """
                #for pk3 in mol.peak[peak]['hnco']:
                #    print('%s %s %s' % ()'H',pk3.f1,GetBMRB('H','H',assNew[peak][-1])))
                #for pk3 in mol.peak[peak]['hnco']:
                #    print('N',pk3.f1,GetBMRB('N','N',assNew[peak][-1]))

                #for pk3 in mol.peak[peak]['hncaco']:
                #    if(pk3.tp=='main'):
                #        print('CO(i)',pk3.f3p,mol.bmrbC[res3]['C'])
                #for pk3 in mol.peak[peak]['hnca']:
                #    if(pk3.tp=='main'):
                #        print('CA(i)',pk3.f3p,mol.bmrbC[res3]['CA'])


                #for pk3 in self.parent.molecule.peak[peak]['ctocsy']:
                #    if(pk3.name=='main'):
                #        print('CA(i)',pk3.f3p,self.bmrbC[assNew[peak][-1]]['CA'])
                return
                for spec in specord:
                    if(spec in self.parent.molecule.peak[pk].keys()):
                        for i,pk3 in enumerate(self.parent.molecule.peak[pk][spec]):
                            print(pk,spec,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.tp,pk3.inty)
                            print('asf')
                            num_items = self.source.GetItemCount()
                            self.source.InsertItem(num_items,str(cnt))
                            self.source.SetItem(num_items,0,str(spec))
                            self.source.SetItem(num_items,1,str(pk3.name))
                            self.source.SetItem(num_items,2,'%.2f' % pk3.f1)
                            self.source.SetItem(num_items,3,'%.2f' % pk3.f2)

                            self.source.SetItem(num_items,4,'%.2f' % pk3.f3p)

                            lab=self.parent.molecule.GetLab(spec,pk3.tp)


    def DoAssign(self,pk,seq):
        assign=[]
        for key,vals in self.inst.subgraphRef.items(): #for each subgraph...
            sub_noe_node_list=vals['nodes']
            sub_noe_adjacency=vals['adj']
            if(pk in sub_noe_node_list):
                break

        pkRef=pk
        pkValRef=self.GetNodeNumber(pk,sub_noe_node_list)
        seqRef=seq
        self.results[pk]=seq
        self.assSeq[seq]=1
        assign.append((pk,seq))

        pkVal=pkValRef
        while(1==1):
            #assign current pk/seq
            forward=[]
            for adj in sub_noe_adjacency[pkVal]:
                if(self.inst.NMR.noes[pk][adj][2]=='f'):
                    forward.append(adj)
            bk=0
            if(len(forward)==1 and forward[0] not in self.results.keys()):
                pk=forward[0]
                seqVal=self.GetNodeNumber(seq,self.G1_nodes)

                if(seqVal+1<len(self.G1_nodes)):
                    seq=self.G1_nodes[seqVal+1]
                else:
                    break

                if seq in self.parent.parent.tabOne.molecule.candidates[pk]:
                    if(pk not in self.results.keys() and seq not in self.assSeq.keys()):
                        #print('forward assignment:',pk,seq)
                        self.results[pk]=seq
                        self.assSeq[seq]=1
                        assign.append((pk,seq))
                        pkVal=self.GetNodeNumber(pk,sub_noe_node_list)
                        bk=1
            if(bk==0):
                break

        pkVal=pkValRef
        pk=pkRef
        seq=seqRef

        #print('reset to ',pk,seq)
        while(1==1):
            backward=[]
            for adj in sub_noe_adjacency[pkVal]:
                if(self.inst.NMR.noes[pk][adj][2]=='b'):
                    backward.append(adj)
            #print(backward)
            bk=0
            if(len(backward)==1 and backward[0] not in self.results.keys()):
                pk=backward[0]
                seqVal=self.GetNodeNumber(seq,self.G1_nodes)

                if(seqVal-1>0):
                    seq=self.G1_nodes[seqVal-1]
                else:
                    break

                seq=self.G1_nodes[seqVal-1]
                if seq in self.parent.parent.tabOne.molecule.candidates[pk]:
                    if(pk not in self.results.keys() and seq not in self.assSeq.keys()):
                        #print('backward assignment:',pk,seq)
                        self.results[pk]=seq
                        self.assSeq[seq]=1
                        assign.append((pk,seq))
                        pkVal=self.GetNodeNumber(pk,sub_noe_node_list)
                        bk=1
            if(bk==0):
                break
        return assign
        #move forwards along adjacency

    def PopulateResults(self):
        num_items = self.source4.GetItemCount()
        for i in range(num_items):
            self.source4.DeleteItem(0)


        for key in self.listy:
            if(key in self.results.keys()):
                val=self.results[key]
                num_items = self.source4.GetItemCount()
                self.source4.InsertItem(num_items,str(num_items))
                self.source4.SetItem(num_items,0,key)
                self.source4.SetItem(num_items,1,val)

    def onClearButton(self,event):
        self.results={}
        self.assSeq={}
        num_items = self.source3.GetItemCount()
        for i in range(num_items):
            self.source3.DeleteItem(0)

        num_items = self.source4.GetItemCount()
        for i in range(num_items):
            self.source4.DeleteItem(0)

        self.draw_figure()
        self.PopulateResults()
        self.PopulateList()

    def OnDoubleClick(self,event):
        if(self.PICK==1):
            sele=self.source.GetFirstSelected()
            count = self.source.GetItemCount()
            col1 = [self.source.GetItem(row, 0).GetText() for row in range(count)][sele]
            col2 = [self.source.GetItem(row, 1).GetText() for row in range(count)][sele]

            self.ComboBox1.SetSelection(self.index[col1])
        elif(self.PICK==2):
            sele=self.source4.GetFirstSelected()
            count = self.source4.GetItemCount()
            col1 = [self.source4.GetItem(row, 0).GetText() for row in range(count)][sele]
            col2 = [self.source4.GetItem(row, 1).GetText() for row in range(count)][sele]

            self.ComboBox1.SetSelection(self.index[col1])
        else:
            return

        self.SELECT=1
        self.G1=(col1,)
        self.G2=col2.split()
        print('selected:',self.G1,self.G2)
        self.draw_figure()

    def OnDoubleClick2(self,event):
        if(self.PICK==3):
            sele=self.source2.GetFirstSelected()
            count = self.source2.GetItemCount()
            col1 = [self.source2.GetItem(row, 0).GetText() for row in range(count)][sele]
            col2 = [self.source2.GetItem(row, 1).GetText() for row in range(count)][sele]

            self.ComboBox2.SetSelection(self.indexSeq[col1])

        elif(self.PICK==4):
            sele=self.source4.GetFirstSelected()
            count = self.source4.GetItemCount()
            col1 = [self.source4.GetItem(row, 1).GetText() for row in range(count)][sele]
            col2 = [self.source4.GetItem(row, 0).GetText() for row in range(count)][sele]

            self.ComboBox2.SetSelection(self.indexSeq[col1])

        else:
            sele=self.source2.GetFirstSelected()
            count = self.source2.GetItemCount()
            col1 = [self.source2.GetItem(row, 0).GetText() for row in range(count)][sele]
            col2 = [self.source2.GetItem(row, 1).GetText() for row in range(count)][sele]

            self.ComboBox2.SetSelection(self.indexSeq[col1])


        self.SELECT=1
        self.G2=(col1,)
        self.G1=col2.split()
        print('selected:',self.G1,self.G2)
        self.draw_figure()

    def OnDoubleClick3(self,event):

        sele=self.source2.GetFirstSelected()
        count = self.source2.GetItemCount()
        col1 = [self.source2.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.source2.GetItem(row, 1).GetText() for row in range(count)][sele]

        self.SELECT=1
        self.G2=(col1,)
        self.G1=col2.split()
        print('selected:',self.G1,self.G2)
        self.draw_figure()

    def OnDoubleClick4(self,event):

        sele=self.source2.GetFirstSelected()
        count = self.source2.GetItemCount()
        col1 = [self.source2.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.source2.GetItem(row, 1).GetText() for row in range(count)][sele]

        self.SELECT=1
        self.G2=(col1,)
        self.G1=col2.split()
        print('selected:',self.G1,self.G2)
        self.draw_figure()

    def on_draw_button(self, event):
        self.PopulateList()
        self.draw_figure()

    def onButtonError(self,event):
        ERRORS=[]
        for pk in self.parent.parent.tabOne.molecule.peak.keys():
            entry=self.parent.parent.tabOne.molecule.GetErrors(pk)
            if(len(entry)!=0):
                ERRORS.append(('',pk,))
                ERRORS.append(entry)
        print(ERRORS)
        textEdit.MyFrame(ERRORS,stream='y')


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
        self.parent.parent.tabOne.molecule.normSpec()
        #self.parent.parent.tabOne.molecule.assSpec()
        self.parent.parent.tabOne.molecule.AssCACB()
        #self.parent.parent.tabOne.molecule.EdgeScreen()
        #self.parent.parent.tabOne.molecule.WriteInit()
