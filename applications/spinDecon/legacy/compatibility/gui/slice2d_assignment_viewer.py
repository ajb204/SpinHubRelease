"""Quarantined historical Slice2D assignment/conn_data viewer.

Extracted from gui.workspaces.slice2d during Stage 148B. This source is retained
for possible future recovery and is not part of the active application.
"""

class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

    def Update(self,dicty):
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty
        print(dicty[0])



# 4D NOE COMPATIBILITY BRIDGE
# DEPRECATED: AssMan/AssManFrame are the legacy conn_data viewer.
# conn_data is no longer authoritative; use Frames.peakListFrame.PeakListFrame
# with mode='full'.  Retained temporarily for backwards compatibility only.
class AssMan(wx.App):
    """DEPRECATED launcher for the legacy ``conn_data`` peak-list viewer."""
    def __init__(self,inherit):
        self.frame_AssManFrame=AssManFrame(None,10,'Assignment',inherit)
        self.frame_AssManFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class AssManFrame(wx.Frame):
    """DEPRECATED ``conn_data`` viewer; scheduled for future removal."""
#    title = 'AssBox'
    def __init__(self,parent, id, title,inherit):
        self.parent=inherit
        self._init_ctrls(parent)


    def _init_ctrls(self,prnt):
        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=prnt,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title='ListBox Test ...')
        self.SetClientSize(wx.Size(900, 280))

        panel=wx.Panel(self,-1)
        self.corrDict={}
        self.sym=self.slice_service.connections[0].sym

        self.lc=SortedListCtrl(panel,self.corrDict)

        cnt=0
        self.lc.InsertColumn(cnt, 'Starting resonance');cnt+=1
        if(self.sym=='y'):
            self.lc.InsertColumn(cnt, 'Ending resonance');cnt+=1
        self.lc.InsertColumn(cnt, 'f1(ppm)');cnt+=1
        self.lc.InsertColumn(cnt, 'f2(ppm)');cnt+=1
        self.lc.InsertColumn(cnt, 'f3(ppm)');cnt+=1
        if(self.slice_service.spectral_dimension == 4):
            self.lc.InsertColumn(cnt, 'f4(ppm)');cnt+=1
        self.lc.InsertColumn(cnt, 's/n 1');cnt+=1
        if(self.sym=='y'):
            self.lc.InsertColumn(cnt, 's/n 2');cnt+=1
            self.lc.InsertColumn(cnt, 'diff');cnt+=1
            self.lc.InsertColumn(cnt, 'IntDist');cnt+=1
            self.lc.InsertColumn(cnt, 'Shift');cnt+=1
            self.lc.InsertColumn(cnt, 'Confidence');cnt+=1
            self.lc.InsertColumn(cnt, 'Kilter');cnt+=1

        self.lc.SetColumnWidth(0, 140)
        self.lc.SetColumnWidth(1, 153)


        #outy=open('confNOE.list','w')
        #for key,data in items:
            #outy.write('%s\t' % str(data[0]))
            #outy.write('%s\t' % str(data[1]))
            #for i in range(len(data)-2):
            #    outy.write('%s\t' % str(data[i+2]))  #add atom
            #outy.write('\n')
        #outy.close()

        self.Addbutton =  wx.Button(panel, 10, 'Show',(710,10))
        self.Nextbutton =  wx.Button(panel, -1, 'Next',(710,10))
        self.Previousbutton =  wx.Button(panel, -1, 'Previous',(710,10))

        self.Removebutton= wx.Button(panel, 11, 'Remove',(710,60))
        #self.Clearbutton = wx.Button(panel, 12, 'Clear',(710,110))
        self.Closebutton = wx.Button(panel, 13, 'Close',(710,160))
        self.Savebutton = wx.Button(panel, 14, 'Save',(710,210))
        self.Loadbutton = wx.Button(panel, -1, 'Load',(710,210))

        # self.cb_remAuto = wx.CheckBox(panel, -1,"Auto",style=wx.ALIGN_RIGHT)
        # self.cb_remSing = wx.CheckBox(panel, -1,"Orph",style=wx.ALIGN_RIGHT)
        # self.cb_remReci = wx.CheckBox(panel, -1,"Recip",style=wx.ALIGN_RIGHT)
        # self.Bind(wx.EVT_CHECKBOX, self.OnRefresh, self.cb_remAuto)
        # self.Bind(wx.EVT_CHECKBOX, self.OnRefresh, self.cb_remSing)
        # self.Bind(wx.EVT_CHECKBOX, self.OnRefresh, self.cb_remReci)

        # self.cb_remAuto.SetValue(1)
        # self.cb_remSing.SetValue(1)
        # self.cb_remReci.SetValue(1)

        #wx.StaticText(self, -1, 'Assman', (0,0))

        #self.pdbfile = infile
        self.textbox = wx.TextCtrl(
            panel,
            size=(150,-1),
            style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.textbox.SetValue("out/MyList.out")
        #self.textbox.SetValue(self.pdbfile)




        #self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnAdd,self.lc)

        self.Bind (wx.EVT_BUTTON, self.OnAdd, self.Addbutton)
        self.Bind (wx.EVT_BUTTON, self.OnNext, self.Nextbutton)
        self.Bind (wx.EVT_BUTTON, self.OnPrevious, self.Previousbutton)
        self.Bind (wx.EVT_BUTTON, self.OnRemove, self.Removebutton)
        #self.Bind (wx.EVT_BUTTON, self.OnClear, self.Clearbutton)
        self.Bind (wx.EVT_BUTTON, self.OnClose, self.Closebutton)
        self.Bind (wx.EVT_BUTTON, self.OnSaveResults, self.Savebutton)
        self.Bind (wx.EVT_BUTTON, self.OnLoadResults, self.Loadbutton)

        #self.vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.lc, 1, wx.EXPAND)

        vbox=wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.Addbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Nextbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Previousbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Removebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.Clearbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Closebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Savebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Loadbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.textbox, 0, wx.ALIGN_CENTER| wx.TOP)

        # vbox.Add(self.cb_remAuto, 0, wx.ALIGN_CENTER| wx.TOP)
        # vbox.Add(self.cb_remSing, 0, wx.ALIGN_CENTER| wx.TOP)
        # vbox.Add(self.cb_remReci, 0, wx.ALIGN_CENTER| wx.TOP)


        hbox.Add(vbox)
        panel.SetSizer(hbox)



        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        #hbox  = wx.BoxSizer(wx.HORIZONTAL)
        #hbox.Add(self.Addbutton, 1, wx.EXPAND)
        #hbox.Add(self.Removebutton, 1, wx.EXPAND)
        #hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        #hbox.Add(self.Closebutton, 1, wx.EXPAND)
        #hbox.Add(self.Savebutton, 1, wx.EXPAND)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        self.Centre()
        self.Show(True)
        self.OnRefresh(True)





    def OnRefresh(self, event):
        print('Refreshing list')
        self.lc.DeleteAllItems()
        corr=[]
        self.corrDict={}
        for i,cn in enumerate(self.slice_service.connections):
            add=[]
            add.append(cn.p1)
            if(cn.sym=='y'):
                add.append(cn.p2)
            add.append(cn.f1)
            add.append(cn.f2)
            add.append(cn.f3)
            if(self.slice_service.spectral_dimension == 4):
                add.append(cn.f4)
            add.append(cn.s1)
            if(cn.sym=='y'):
                try:
                    distppm=cn.distppm
                    conf=1
                    if(distppm<0.5):
                        conf+=2
                    if(distppm<0.1):
                        conf+=3
                    add.append(cn.s2)
                    add.append(cn.frac)
                    add.append(cn.distScore)
                    add.append(cn.distppm)
                    add.append(conf)
                    add.append(0.0)
                except:
                    pass
            self.corrDict[i]=add

        # if(self.cb_remAuto.IsChecked()==0):
        if (False):
            for key,vals in list(self.corrDict.items()):
                if(vals[0]==vals[1]):
                    del self.corrDict[key]

        # if(self.cb_remSing.IsChecked()==0):
        if (False):
            rem=[]
            for key,vals in list(self.corrDict.items()):
                tick=0
                for koi,vols in list(self.corrDict.items()):
                    if(vals[0]==vols[1] and vals[1]==vols[0]):
                        tick=1
                        break
                if(tick==0):
                    if(key not in rem):
                        rem.append(key)

            for re in rem:
                del self.corrDict[re]


        # if(self.cb_remReci.IsChecked()==0):
        if (False):
            rem=[]
            for key,vals in list(self.corrDict.items()):
                for koi,vols in list(self.corrDict.items()):
                    if(vals[0]==vols[1] and vals[1]==vols[0]):
                        if(key not in rem):
                            rem.append(key)
                        if(koi not in rem):
                            rem.append(koi)
                        break
            for re in rem:
                del self.corrDict[re]



        for key,data in list(self.corrDict.items()):
            num_items = self.lc.GetItemCount()
            self.lc.InsertItem(num_items,str(data[0]))  #add assignment
            self.lc.SetItem(num_items, 0,str(data[0])) #add atom
            self.lc.SetItem(num_items, 1,str(data[1]))  #add atom

            for i in range(len(data)-2):
                self.lc.SetItem(num_items, i+2,str(data[i+2]))  #add atom
            self.lc.SetItemData(num_items, key)
        self.lc.Update(self.corrDict)



    def OnSaveResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save peaklist...",
            defaultDir=os.path.join(os.getcwd(),'out'),
            #defaultFile=os.path.split(self.deconParFile)[1],
            defaultFile='MyList.out',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.outfile=path
            self.OnSave()
            #self.deconParFile=path
            #if(os.path.exists(self.deconParFile)==0):
            #    outy=open(self.deconParFile,'w');outy.close()
            #self.notebook.tabMagma.deconParFile=path
            #self.notebook.tabMagma.OnButtonSave(True)

            #self.notebook.tabOne.deconParFile=path
            #self.notebook.tabOne.OnButtonSave(True)

            ##self.canvas.print_figure(path, dpi=self.dpi)
            #self.flash_status_message("Saved %s" % path)

    #FGA added
    def OnLoadResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Load session...",
            #defaultDir=os.getcwd(),
            defaultDir=os.path.join(os.getcwd(),'out'),
            defaultFile="",
            wildcard=file_choices,
            style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:

            self.outfile=dlg.GetPath()
            self.OnLoad()
    def onItemSelected(self, event):
        """"""

        #currentItem = event.m_itemIndex
        #car = self.corrDict[currentItem]
        #print car

        #count = self.lc.GetItemCount()
        #self.sorted_artists = [self.list.GetItem(itemId=row, col=0).GetText() for row in xrange(count)]
        #print self.sorted_artists
        #print self.sorted_artists[currentItem]


    def AtoI(self,val):
        for i in range(len(self.slice_service.peaks)):
            if(val==self.slice_service.peaks[i].name):
                return i

    def OnAdd(self, event):
        sele=self.lc.GetFirstSelected()
        count = self.lc.GetItemCount()
        if(self.sym=='y'):
            col1 = [self.lc.GetItem(row, 0).GetText() for row in range(count)][sele]
            col2 = [self.lc.GetItem(row, 1).GetText() for row in range(count)][sele]
            #print col1,col2,self.AtoI(col1)
            self.parent.ComboBox1.SetSelection(self.AtoI(col1))
            self.parent.ComboBox2.SetSelection(self.AtoI(col2))
        else:
            col1 = [self.lc.GetItem(row, 0).GetText() for row in range(count)][sele]
            logging.info(col1)
            resn=col1.split('_')[0]
            logging.info(resn)
            #print col1,col2,self.AtoI(col1)
            self.parent.ComboBox1.SetSelection(self.AtoI(resn))
            #self.parent.ComboBox2.SetSelection(self.AtoI(col2))

        self.parent.on_draw_button(True)


    def OnNext(self, event):
        sele=self.lc.GetFirstSelected()
        self.lc.Select(sele,on=0)
        self.lc.Select(sele+1,on=1)
        self.OnAdd(True)
    def OnPrevious(self, event):
        sele=self.lc.GetFirstSelected()
        self.lc.Select(sele,on=0)
        self.lc.Select(sele-1,on=1)
        self.OnAdd(True)



    def OnRemove(self, event):
        if(self.sym=='y'):
            index = self.lc.GetFocusedItem()
            count=self.lc.GetItemCount()
            item1 = [self.lc.GetItem(row, 0).GetText() for row in range(count)][index]
            item2 = [self.lc.GetItem(row, 1).GetText() for row in range(count)][index]
            logging.info(len(self.slice_service.connections),count)
            logging.info('looking for:',item1,item2)
            for i,cn in enumerate(self.slice_service.connections):
                if(cn.p1==item1 and cn.p2==item2):
                    logging.info('removing:',cn.tag)
                    self.slice_service.connections.pop(i)
                    break
            logging.info('looking for:',item2,item1)
            for i,cn in enumerate(self.slice_service.connections):
                if(cn.p1==item2 and cn.p2==item1):
                    logging.info('removing:',cn.tag)
                    self.slice_service.connections.pop(i)
                    break
        else:
            index = self.lc.GetFocusedItem()
            count=self.lc.GetItemCount()
            item1 = [self.lc.GetItem(row, 0).GetText() for row in range(count)][index]
            logging.info('looking for:',item1)
            for i,cn in enumerate(self.slice_service.connections):
                if(cn.p1==item1):
                    logging.info('removing:',cn.p1)
                    self.slice_service.connections.pop(i)
                    break

        self.OnRefresh(True)
        logging.info(len(self.slice_service.connections))
        self.slice_service.refresh_status()
        self.parent.draw_figure()

        """
        print 'removingL:',item1,item2
        self.lc.DeleteItem(index)
        count=self.lc.GetItemCount()
        for itemIndex in range(count):
            item1f = [self.lc.GetItem(row, 0).GetText() for row in xrange(count)][itemIndex]
            item2f = [self.lc.GetItem(row, 1).GetText() for row in xrange(count)][itemIndex]
            if(item1f==item2 and item2f==item1):
                print 'removingL:',item1f,item2f
                self.lc.DeleteItem(itemIndex)
                break
        count=self.lc.GetItemCount()

        """
        logging.info(len(self.slice_service.connections))
        self.slice_service.refresh_status()
        self.parent.draw_figure()

    def OnClose(self, event):
        self.Close()

    def OnClear(self, event):
        self.lc.DeleteAllItems()

    def OnSave(self):
        logging.info()
        logging.info('Saving list to ',self.outfile)
        outy=open(self.outfile,'w')
        for entry in self.slice_service.connections:
            if(self.sym=='y'):
                if(self.slice_service.spectral_dimension==3):
                    outy.write('%s\t%s\t%f\t%f\t%f\t%e\t%e\t%e\t%e\t%f\t%f\t%f\t%f\n' % (entry.p1,entry.p2,entry.f1,entry.f2,entry.f3,entry.s1,entry.s2,entry.n1,entry.n2,entry.frac,entry.distScore,entry.distppm,entry.Intscore))
                else:
                    outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (entry.p1,entry.p2,entry.f1,entry.f2,entry.f3,entry.f4,entry.s1,entry.s2,entry.n1,entry.n2,entry.frac,entry.distScore,entry.distppm,entry.Intscore))
            else:
                outy.write('%s\t%f\t%f\t%f\t%e\n' % (entry.p1,entry.f1,entry.f2,entry.f3,entry.s1))

        outy.close()
        self.slice_service.refresh_status()


    def OnLoad(self):
        #self.outfile=self.textbox.GetValue()
        logging.info("loading ",self.outfile)
        self.slice_service.load_connections(self.outfile, sym=self.sym)
        self.OnRefresh(True)
        self.parent.draw_figure()


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
