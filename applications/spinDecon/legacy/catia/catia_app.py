#!/usr/bin/python

###################################################################
# Deconvolve nmr spectrum
###################################################################
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

import numpy,sys,os,wx,platform,glob,shutil
import nmrglue as ng
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import subprocess
import re,copy
from spinDecon.parameter_store import read_structured_parameter_file, write_structured_parameter_file




from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin

#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

import wx.lib.mixins.listctrl  as  listmix
class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
    ''' TextEditMixin allows any column to be edited. '''
    
    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """Constructor"""
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)

        

class catiaFrame(wx.Panel):

    def __init__(self,parent,deconParFile):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.state = getattr(parent, "state", None)
        # CATIA workspace follows the main GUI SpecPath rather than a hard-coded
        # project/raw or current-directory location.
        self.raw = self._main_spec_dir()
        self.savefile=os.path.join(self.raw, 'catia.save')
        # Set sizer for the frame, so we can change frame size to match widgets
        self.windowSizer = wx.BoxSizer()
        self.windowSizer.Add(self, 1, wx.ALL | wx.EXPAND)

        self.datasetBox()

    def _main_spec_dir(self):
        """Return the resolved SpecPath directory owned by the main GUI/state."""
        obj = self
        while obj is not None:
            state = getattr(obj, 'state', None)
            if state is not None and hasattr(state, 'spec_dir'):
                try:
                    return os.path.normpath(state.spec_dir())
                except Exception:
                    pass
            obj = getattr(obj, 'parent', None)
        spec_path = self._main_spec_path()
        if os.path.isabs(spec_path):
            return os.path.normpath(spec_path)
        return os.path.normpath(os.path.join(self._main_working_dir(), spec_path))

    def _main_working_dir(self):
        """Return the project working directory used to anchor relative paths."""
        obj = self
        while obj is not None:
            state = getattr(obj, 'state', None)
            if state is not None:
                working_dir = str(getattr(state, 'working_dir', '') or '').strip()
                if working_dir:
                    return os.path.normpath(working_dir)
            obj = getattr(obj, 'parent', None)
        return '.'

    def _main_spec_path(self):
        """Return the SpecPath configured by the main GUI/state."""
        obj = self
        while obj is not None:
            state = getattr(obj, 'state', None)
            if state is not None:
                spec_path = str(getattr(state, 'spec_path', '') or '').strip()
                if spec_path:
                    return spec_path
            ctrl = getattr(obj, 'specPathBox', None)
            if ctrl is not None:
                try:
                    spec_path = str(ctrl.GetValue() or '').strip()
                    if spec_path:
                        return spec_path
                except Exception:
                    pass
            obj = getattr(obj, 'parent', None)
        return './spec'

    def _dataset_spec_dir(self, dataset_root):
        """Resolve a selected dataset's processed-data directory using SpecPath."""
        spec_path = self._main_spec_path()
        if os.path.isabs(spec_path):
            return os.path.normpath(spec_path)
        return os.path.normpath(os.path.join(dataset_root, spec_path))

    def datasetBox(self):


        
        self.sizer = wx.GridBagSizer(1, 6)

        cnt=0
        self.buttonPlus = wx.Button(self, label="+", size = (-1,22))
        self.buttonMinus = wx.Button(self, label="-", size = (-1,22))
        self.buttonSave = wx.Button(self, label="Save", size = (-1,22))
        self.buttonLoad = wx.Button(self, label="Load", size = (-1,22))
        

        self.buttonRefresh = wx.Button(self, label="Refresh", size = (-1,22))
        
        
        self.RunCatia = wx.Button(self, label="Run Catia", size = (-1,22))

        self.buttonPlus.Bind(wx.EVT_BUTTON, self.OnButtonPlus)
        self.buttonMinus.Bind(wx.EVT_BUTTON, self.OnButtonMinus)
        self.buttonSave.Bind(wx.EVT_BUTTON, self.OnButtonSave)
        self.buttonLoad.Bind(wx.EVT_BUTTON, self.OnButtonLoad)
        self.buttonRefresh.Bind(wx.EVT_BUTTON, self.OnButtonRefresh)
        self.RunCatia.Bind(wx.EVT_BUTTON, self.OnRunCatia)
        
        self.sizer.Add(self.buttonPlus, (0, 0), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1
        self.sizer.Add(self.buttonMinus, (0, 1), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1
        self.sizer.Add(self.buttonSave, (0, 2), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1
        self.sizer.Add(self.buttonLoad, (0, 3), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1
        self.sizer.Add(self.buttonRefresh, (0, 4), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1
        self.sizer.Add(self.RunCatia, (0, 5), border=10, flag=wx.ALIGN_CENTER_VERTICAL | wx.RIGHT);cnt+=1


        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.sizer)

        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.dataLbl = wx.StaticBox(self, -1, 'Data:')
        self.dataSizer = wx.StaticBoxSizer(self.dataLbl, wx.VERTICAL)
        #self.dataSizer.Add(self.buttonProcess)
        self.lc=AutoWidthListCtrl(self)
        #self.datasets.SetMinSize((650,300))
        #self.datasets.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)


        self.lc.Bind(wx.EVT_LIST_COL_CLICK, self.OnButtonSort)
        self.lc.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)
        self.lc.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRClick)
        #self.lc.Bind(wx.EVT_LISTBOX_DCLICK, self.OnDClick)
        self.dataSizer.Add(self.lc)

        
        self.peakLbl = wx.StaticBox(self, -1, 'Peaks:')
        self.peakSizer = wx.StaticBoxSizer(self.peakLbl, wx.VERTICAL)
        self.datasets=AutoWidthListCtrl(self)
        self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        self.peakSizer.Add(self.datasets)

        
        hbox1.Add(self.dataSizer)
        hbox1.Add(self.peakSizer)
        vbox.Add(hbox1)
        
        #peaks

        #self.datasets.SetMinSize((650,300))
        #self.datasets.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClickPeak)
        #self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        #self.datasets.Bind(wx.EVT_LIST_COL_CLICK, self.OnButtonSort)

        

        self.parLocal = EditableListCtrl(self,style=wx.LC_REPORT)
        self.parLocal.SetMinSize((300,300))
        locbox = wx.StaticBox(self,-1,'Local Parameters:')
        locboxS=wx.StaticBoxSizer(locbox,wx.VERTICAL)
        locboxS.Add(self.parLocal)


        
        self.fitLocal = AutoWidthListCtrl(self)
        self.fitLocal.SetMinSize((300,300))
        fitbox = wx.StaticBox(self,-1,'Fitted Parameters:')
        fitboxS=wx.StaticBoxSizer(fitbox,wx.VERTICAL)
        fitboxS.Add(self.fitLocal)

        
        
        #self.lc.Bind(wx.EVT_LISTBOX_DCLICK, self.OnAdd)
        
        self.SetupLC()
        self.OnButtonLoad(True)

        
        
        #self.border.Add(self.sizer)
        #self.dataSizer.Add(self.border)
        #self.dataSizer.AddSpacer(5)
        #self.dataSizer.Add(self.lc)
        #self.dataSizer.AddSpacer(5)

        hbox=wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(locboxS)
        hbox.Add(fitboxS)
        
        vbox.Add(hbox)
        

        self.border = wx.BoxSizer()

        self.fullSizer=wx.BoxSizer(wx.HORIZONTAL)
        self.fullSizer.Add(vbox)

        self.create_main_panel()
        self.border.Add(self.fullSizer)
        self.SetSizerAndFit(self.border)

        self.OnButtonLoad(True)
        

    def create_main_panel(self):
        
        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        #self.canvas.mpl_connect('button_press_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)


        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)
        #self.SetSizer(self.vbox)
        #self.vbox.Add(fitbox)
        self.vbox.Fit(self)        
        self.fullSizer.Add(self.vbox)



    def ReadCatia(self,pk):

        self.fitLocal.ClearAll()
        self.fitLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.fitLocal.InsertColumn(1, 'Value', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.fitLocal.InsertColumn(2, 'Status', width = 80,format=wx.LIST_FORMAT_CENTRE)
        num_items = self.fitLocal.GetItemCount()

        infile=self.raw+'/catia/OutPut/'+pk+'.dat'
        if(os.path.exists(infile)==False):
            return False

        self.catiaXr=[]
        self.catiaYr=[]
        self.catiaEr=[]
        self.catiaFr=[]
        self.catiaLr=[]

        self.catiaX=[]
        self.catiaY=[]
        self.catiaE=[]
        self.catiaF=[]
        self.catiaL=''

        
        locFlag=0
        cnt=0
        inny=open(infile)
        for line in inny.readlines():
            if(len(line)>0):
                if(line[0]!='#'):
                    test=line.split()
                    if(len(test)==4):
                        self.catiaX.append(float(test[0]))
                        self.catiaY.append(float(test[1]))
                        self.catiaE.append(float(test[2]))
                        self.catiaF.append(float(test[3]))
                else:
                    test=line.split()
                    if(test[0]=='#DataSet:'):

                        if(len(self.catiaX)>0):
                            self.catiaX=numpy.array(self.catiaX)
                            self.catiaY=numpy.array(self.catiaY)
                            self.catiaE=numpy.array(self.catiaE)
                            self.catiaF=numpy.array(self.catiaF)

                            self.catiaXr.append(self.catiaX)
                            self.catiaYr.append(self.catiaY)
                            self.catiaEr.append(self.catiaE)
                            self.catiaFr.append(self.catiaF)
                            self.catiaLr.append(self.catiaL)

                            self.catiaX=[]
                            self.catiaY=[]
                            self.catiaE=[]
                            self.catiaF=[]

                            
                        self.catiaL=test[1]

                    if(test[0]=='#Atom:'):
                        locFlag=1
                        
                    if(locFlag==0 and len(test)==4):
                        self.fitLocal.InsertItem(num_items,str(cnt))
                        self.fitLocal.SetItem(num_items,0,test[1]) #id

                        try:
                            self.fitLocal.SetItem(num_items,1,'%.3f ' % (float(test[2]))) #A
                        except:
                            self.fitLocal.SetItem(num_items,1,'%s ' % ((test[2]))) #A


                        try:
                            self.fitLocal.SetItem(num_items,2,'%.3f ' % (float(test[3]))) #A
                        except:
                            self.fitLocal.SetItem(num_items,2,'%s ' % ((test[3]))) #A

                            
                        cnt+=1

        if(len(self.catiaX)>0):
            self.catiaX=numpy.array(self.catiaX)
            self.catiaY=numpy.array(self.catiaY)
            self.catiaE=numpy.array(self.catiaE)
            self.catiaF=numpy.array(self.catiaF)

            self.catiaXr.append(self.catiaX)
            self.catiaYr.append(self.catiaY)
            self.catiaEr.append(self.catiaE)
            self.catiaFr.append(self.catiaF)
            self.catiaLr.append(self.catiaL)
                        

        self.catiaXr=numpy.array(self.catiaXr)
        self.catiaYr=numpy.array(self.catiaYr)
        self.catiaEr=numpy.array(self.catiaEr)
        self.catiaFr=numpy.array(self.catiaFr)

        self.chi2Global=numpy.average((self.catiaYr-self.catiaFr)**2.)

        

        self.globout=self.raw+'/catia/OutPut/GlobalParam.fit'
        if(os.path.exists(self.globout)==1):
            inny=open(self.globout)
            for line in inny.readlines():
                #test=line.split('=')
                #if(len(test)==2):
                #        self.fitLocal.InsertStringItem(num_items,str(cnt))
                #        self.fitLocal.SetStringItem(num_items,0,test[1]) #id
                    
                test=line.split()
                if(len(test)==3):
                    self.fitLocal.InsertItem(num_items,str(cnt))
                    self.fitLocal.SetItem(num_items,0,test[0]+' (global)') #id
                    try:
                        self.fitLocal.SetItem(num_items,1,'%.3f ' % (float(test[1]))) #A
                    except:
                        self.fitLocal.SetItem(num_items,1,'%s ' % ((test[1]))) #A

                    try:
                        self.fitLocal.SetItem(num_items,2,'%.3f ' % (float(test[2]))) #A
                    except:
                        self.fitLocal.SetItem(num_items,2,'%s ' % ((test[2]))) #A

                            
                    cnt+=1


        #if(pk not in self.cpmgLocal.keys()):
        #    self.DoFit(pk)
                    
        #for par in 'R0line','pb','kex','R0','dw':
        #    self.fitLocal.InsertStringItem(num_items,str(cnt))
        #    self.fitLocal.SetStringItem(num_items,0,par+' (local)') #id
        #    val=self.cpmgLocal[pk][par]
        #    self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (val)) #A
        #    cnt+=1
            

        #self.fitLocal.InsertStringItem(num_items,str(cnt))
        #self.fitLocal.SetStringItem(num_items,0,'chi2Line') #id
        #self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Line)) #A
        #cnt+=1

        #self.fitLocal.InsertStringItem(num_items,str(cnt))
        #self.fitLocal.SetStringItem(num_items,0,'chi2Local') #id
        #self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Local)) #A
        #cnt+=1

        try:
            self.fitLocal.InsertStringItem(num_items,str(cnt))
            self.fitLocal.SetStringItem(num_items,0,'chi2Global') #id
            self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Global)) #A
            cnt+=1
        except:
            pass

            
            
        return True

        
    def draw_figure(self,event):

        sele=self.datasets.GetFirstSelected()
        #print(sele)
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)][sele]
        print('Selected:',col1,col2)

        #if(self.ReadFuda(col1)==False):
        #    return
        
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.clear()




        #self.DoFit(col1)
        if(self.ReadCatia(col1)==True):

            for i in range(len(self.catiaXr)):
                self.ax.errorbar(self.catiaXr[i],self.catiaYr[i],yerr=self.catiaE,fmt='o')                
                self.ax.plot(self.catiaXr[i],self.catiaFr[i],label=self.catiaLr[i])
                #self.ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.catiaE,fmt='o',label='raw')                
            #self.ax.plot(self.cpmgX,self.cpmgL,label='line')
            #self.ax.plot(self.cpmgX,self.cpmgF,label='local')

        #else:
            #self.ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.cpmgE,fmt='o',label='raw')                
            #self.ax.plot(self.cpmgX,self.cpmgL,label='line')
            #self.ax.plot(self.cpmgX,self.cpmgF,label='local')


        self.ax.legend(loc='upper right')
        self.ax.set_xlabel("nu_cpmg")
        self.ax.set_ylabel("R2_eff")
        self.canvas.draw()
        
        return

        
        
    def UpdatePeaks(self):
        self.datasets.ClearAll()
        self.datasets.InsertColumn(0, 'Peak', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        #self.datasets.InsertColumn(1, 'Sets', width = 80,format=wx.LIST_FORMAT_CENTRE)
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row,0).GetText() for row in range(count)]
        col2 = [self.lc.GetItem(row,1).GetText() for row in range(count)]
        col3 = [self.lc.GetItem(row,2).GetText() for row in range(count)]
        for (c1,c2,c3) in zip(col1,col2,col3):
            self.datasets.InsertColumn(1, c1, width = 80,format=wx.LIST_FORMAT_CENTRE)
            


            

        
    def OnButtonPlus(self, event):
        
        #if self.contentNotSaved:
        #    if wx.MessageBox("Current content has not been saved! Proceed?", "Please confirm",
        #                     wx.ICON_QUESTION | wx.YES_NO, self) == wx.NO:
        #        return

        # otherwise ask the user what new file to open
        with wx.DirDialog(self, "Open dataset folder",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as DirDialog:
            
            if DirDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            # Proceed loading the file chosen by the user
            pathname = DirDialog.GetPath()

            num_items = self.lc.GetItemCount()
            
            self.lc.InsertItem(num_items,num_items)  #add assignment
            
            lab=self.GetLabel(pathname)
            
            self.lc.SetItem(num_items, 0,lab) #add atom        
            self.lc.SetItem(num_items, 1,str(True))  #add atom            
            self.lc.SetItem(num_items, 2,pathname) #add atom        


            #print(pathname)
            #try:
            #    with open(pathname, 'r') as file:
            #        self.doLoadDataOrWhatever(file)
            #except IOError:
            #    wx.LogError("Cannot open file '%s'." % newfile)


    def GetLabel(self,path):
        fle=os.path.join(self._dataset_spec_dir(path), 'catia', 'frame.save')
        
        if(os.path.exists(fle)==False):
            print('cannot find save file:',fle)
            return path.split('/')[-1]


        sfrq=''
        temp=''
        inny=open(fle)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(test[0]=='set' and test[1]=='sfrq'):
                    sfrq=test[2]
                elif(test[0]=='set' and test[1]=='temperature'):
                    temp=test[2]

        if(temp=='' and sfrq==''):
            return path.split('/')[-1]
        else:
            return sfrq+'_'+temp
                    
    def OnButtonMinus(self,event):
        sele=self.lc.GetFirstSelected()
        self.lc.DeleteItem(sele)


    def OnButtonRefresh(self,event):

        #first, loop over the various dataset files, get the paramset files
        #read in the dataset files
        #fix parameters that need to be fixed.
        #crack on!

        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row,0).GetText() for row in range(count)]
        col2 = [self.lc.GetItem(row,1).GetText() for row in range(count)]
        col3 = [self.lc.GetItem(row,2).GetText() for row in range(count)]

        #self.datfile=[]
        #self.globfile=[]
        #self.locfile=[]

        self.sett=[]
        self.local=[]

        self.UpdatePeaks() #clear peak selections
        
        for i,(c1,c2,c3) in enumerate(zip(col1,col2,col3)):
            if(c2=='True'):
                print('Including',c1,c3)

                #1. read in sfrq, temp and hence parset
                self.ReadLocalPars(i,c3)

                #dat=os.path.join(self._dataset_spec_dir(c3), 'catia', 'dataset.inp')
                #if(os.path.exists(dat)):
                #   self.datfile.append(dat)


        self.SetLocal()


        
    def OnRunCatia(self,event):

        self.PathExists((self.raw+'/catia',self.raw+'/catia/dataset',self.raw+'/catia/OutPut'))

        self.datfile=[]
        self.globfile=[]
        self.locfile=[]

        self.KexFit=True
        self.PbFit=True
        self.Conv=1E-3
        self.MaxIter=100
        
        self.WriteCatiaDataset()
        self.WriteCatiaPar()
        self.WriteCatiaFile()

        
        for old_output in glob.glob(os.path.join(self.raw, 'catia', 'OutPut', '*')):
            if os.path.isdir(old_output):
                shutil.rmtree(old_output)
            else:
                os.remove(old_output)
        if(os.uname()[0]=='Darwin'): #run catia
            os.system('catia_Darwin_i386 < '+self.raw+'/catia/Fit.catia')
        else:
            os.system('catia_Linux_x86_64 < '+self.raw+'/catia/Fit.catia')
        
        pass

    def AddIfNew(self,arr,val):
        for a in arr:
            if(val[0]==a[0]):
                return
        arr.append(val)
        
    def ReadLocalPars(self,idd,path):
        infile=os.path.join(self._dataset_spec_dir(path), 'catia', 'frame.save')
        if(os.path.exists(infile)==False):
            print('Cannot find input file.')
            return


        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(test[0]=='set'):
                self.sett.append((test[1],test[2])) #set pars 

            elif(test[0]=='par'):

                if(test[2]=='1'): #val col
                    val=[]
                    val.append(test[1])
                    val.append(test[3])
                elif(test[2]=='2'): #fit col
                    if(test[3]=='fit'):
                        val.append('fit')
                    else:
                        val.append('')
                        
                elif(test[2]=='3'): #file col
                    if(test[3]=='None'):
                        val.append('')
                    else:
                        val.append(test[3])
                        
                    print('adding:',val)
                    
                    self.AddIfNew(self.local,val)
                    
                    #self.local.append(val)   #local pars
                    val=[]
                print (val)
            elif(test[0]=='peak'): #set peaks
                
                if(test[2]=='True'):
                    #if(test[1] not in self.peaks.keys()):
                    #    self.peaks[test[1]]=[]
                    #self.peaks[test[1]].append(idd)

                    count = self.datasets.GetItemCount()
                    col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
                    tick=0
                    for i,c1 in enumerate(col1):
                        if(test[1]==c1):#found.
                            self.datasets.SetItem(i,idd+1,'True')
                            tick=1
                    if(tick==0):
                        num_items = self.datasets.GetItemCount()
                        self.datasets.InsertItem(num_items,num_items)
                        self.datasets.SetItem(num_items,0,test[1])
                        self.datasets.SetItem(num_items,idd+1,'True')

                        
                    #for i,c1 in enumerate(col1):
                #    if(c1==test[1]):
                
            #else:
            #    if(test[0]=='seqfil'):
            #        self.seqfilCombo.SetValue(test[1])
            #    elif(test[0]=='basis'):
            #        self.basisCombo.SetValue(test[1])
            #    elif(test[0]=='RexScreen'):
            #        self.RexScreenBox.SetValue(test[1])



    def SetLocal(self):


        #self.setLocal.ClearAll()
        #self.setLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        #self.setLocal.InsertColumn(1, 'Set', width = 80,format=wx.LIST_FORMAT_CENTRE)

        self.parLocal.ClearAll()        
        self.parLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.parLocal.InsertColumn(1, 'Initial', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.parLocal.InsertColumn(2, 'Fit?', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.parLocal.InsertColumn(3, 'File', width = 80,format=wx.LIST_FORMAT_CENTRE)

        num_items = self.parLocal.GetItemCount()
        cnt=0
        for local in self.local:
            cnt+=1
            self.parLocal.InsertStringItem(num_items,str(cnt))
            self.parLocal.SetStringItem(num_items,0,local[0])
            self.parLocal.SetStringItem(num_items,1,str(local[1]))
            self.parLocal.SetStringItem(num_items,2,str(local[2]))
            if(len(local)==4):
                self.parLocal.SetStringItem(num_items,3,str(local[3]))
            
        #num_items = self.setLocal.GetItemCount()
        #cnt=0
        #for local in self.sett:
        #    cnt+=1
        #    self.setLocal.InsertStringItem(num_items,str(cnt))
        #    self.setLocal.SetStringItem(num_items,0,local[0])
        #    self.setLocal.SetStringItem(num_items,1,str(local[1]))
            
            
        return



    def WriteDatasetHeader(self,outy,path):
        indat=os.path.join(self._dataset_spec_dir(path), 'catia', 'dataset.inp')
        if(os.path.exists(indat)==False):
            print('Cannot find dataset file:',indat)
            return False
        inny=open(indat)
        for line in inny.readlines():
            test=line.split()
            if(test[0]=='DataDirectory'):
                break
            outy.write(line)
        inny.close()
    
    def WriteCatiaDataset(self):

        count = self.lc.GetItemCount()
        col1 = numpy.array([self.lc.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.lc.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.lc.GetItem(row, 2).GetText() for row in range(count)])

        
        for ii,(c1,c2,c3) in enumerate(zip(col1,col2,col3)): #for all datasets in record
            datfil=self.raw+'/catia/dataset.'+str(ii)+'.inp'

            self.datfile.append(datfil)
            outy=open(datfil,'w')

            if(c2!='True'):
                continue
            if(self.WriteDatasetHeader(outy,c3)==False):
                continue

            countP = self.datasets.GetItemCount()
            col  = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(countP)])
            colP = numpy.array([self.datasets.GetItem(row, ii+1).GetText() for row in range(countP)])

            print()
            print('LISTY')
            print (col,colP)
            print()
            outy.write('DataDirectory = %s\n' % (os.path.join(self._dataset_spec_dir(c3), 'fit') + os.sep))
            outy.write('Data = (\n')
            for c,cp in zip(col,colP):
                if(cp=='True'):
                    outy.write('[%s;%s.out.cpmg];\n' % (c,c))
            outy.write(')')
            outy.close()

            
            
    def WriteCatiaPar(self):
        #each dataset requires its own local parameters

        #sele=self.setLocal.GetFirstSelected()
        #print(sele)
        #self.GetActualField() #set dfrq,field and fieldlab
        
        #count = self.setLocal.GetItemCount()
        #col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
        #col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])

        
        #sele=self.parLocal.GetFirstSelected()
        #print(sele)
        #count = self.parLocal.GetItemCount()
        #col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        #col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        #col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])

        parfil=self.raw+'/catia/ParamSet.inp'
        print('Writing parameter file:',parfil)


        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])

        outy=open(self.raw+'/catia/ParamSet.inp','w')
        self.locfile.append(self.raw+'/catia/ParamSet.inp')
        outy.write('format = (')
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % c1)
        outy.write(')\n')
        outy.write('* = (')
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % str(c2))        
        outy.write(')\n')
        outy.close()

        self.fix=[]
        for i,(c1,c3) in enumerate(zip(col1,col3)):
            if(c3=='fit'):
                pass
            else:
                self.fix.append(c1)

        
        """
        self.locfile.append(parfil)
        outy=open(parfil,'w')
        outy.write('format = (')
        for i,(c1,c2,c3,c4) in enumerate(self.local):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % c1)
        outy.write(')\n')
        outy.write('* = (')
        for i,(c1,c2,c3,c4) in enumerate(self.local):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % str(c2))        
        outy.write(')\n')
        outy.close()

        self.fix=[]
        for i,(c1,c2,c3,c4) in enumerate(self.local):
            if(c3=='fit'):
                pass
            else:
                self.fix.append(c1)
        """
        self.globfile.append(self.raw+'/catia/ParamGlobal.inp')
        outy=open(self.raw+'/catia/ParamGlobal.inp','w')
        outy.write('kex=1000.\n')
        outy.write('pb=0.02\n')
        outy.close()

        
    def WriteCatiaFile(self):
        outy=open(self.raw+'/catia/Fit.catia','w')
        for i in range(len(self.datfile)):
            outy.write('ReadDataset(%s)  #Data summary file \n' % self.datfile[i])
        for i in range(len(self.locfile)):
            outy.write('ReadParam_Local(%s)      #Local parameters initial\n' % self.locfile[i])
        for i in range(len(self.globfile)):
            outy.write('ReadParam_Global(%s)   #Global parameters initial\n' % self.globfile[i])


        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])
        col4 = numpy.array([self.parLocal.GetItem(row, 3).GetText() for row in range(count)])
        for i,(c1,c4) in enumerate(zip(col1,col4)):
            if(len(c4)>0):
                fil=c4.split()[0]
                if(os.path.exists(fil)):
                    outy.write('ReadParam(%s;%s;0;1)\n' % (c1,c4))
                else:
                    print('Cannt find fix file:',fil)
            
            
        #outy.write('ReadParam(Omega;<SpecPath>/test.ft2.list;0;1)   #peak list\n')
        #if(os.path.exists(os.path.join(self._dataset_spec_dir('.'), 'catia', 'DeltaOmega.inp'))):
        #    outy.write('ReadParam(DeltaO;<SpecPath>/catia/DeltaOmega.inp;0;1)#delta omegas\n')

        outy.write('# Fix all the static parameters\n')

        
        for i in range(len(self.fix)):
            outy.write('FreeLocalParam(all;%s;false)\n' % self.fix[i])
            

        #    outy.write('SetGlobalParam(kex;250)\n')
        #    outy.write('SetGlobalParam(pb;0.05)\n')

        self.PathExists((self.raw+'/catia/OutPut',))
        outy.write('# Deal with global parameters\n')

        if(self.KexFit):
            outy.write('FreeGlobalParam(kex;true)\n')
        else:
            outy.write('FreeGlobalParam(kex;false)\n')
        if(self.PbFit):
            outy.write('FreeGlobalParam(pb;true)\n')
        else:
            outy.write('FreeGlobalParam(pb;false)\n')
        outy.write('# Minimize\n')
        outy.write('echo(\\n)\n')
        outy.write('#Minimize()\n')
        outy.write('Minimize(print=y;tol=%f;maxiter=%i)\n' % (self.Conv,self.MaxIter))
        outy.write('#\n')
        outy.write('#  // Print some files for plotting\n')
        outy.write('PrintParam('+self.raw+'/catia/OutPut/GlobalParam.fit;global)\n')
        outy.write('PrintParam('+self.raw+'/catia/OutPut/DeltaOmega.fit;DeltaO)\n')
        outy.write('PrintData('+self.raw+'/catia/OutPut/)\n')
        outy.write('echo(\n)\n')
        outy.write('ChiSq(all;all)\n')
        outy.write('exit(0)\n')
        outy.close()


    def PathExists(self,test):
        for t in test:
            if(os.path.exists(t)==False):
                print('Creating:',t)
                os.system('mkdir '+t)


    def OnButtonLoad(self,event):
        self.local = []

        if os.path.exists(self.savefile) == False:
            print('Cannot find savefile:', self.savefile)
            return

        data = read_structured_parameter_file(self.savefile)

        for dataset_id, values in data.get('dataset', {}).items():
            count = self.lc.GetItemCount()
            col1 = numpy.array([self.lc.GetItem(row, 0).GetText() for row in range(count)])

            if dataset_id not in col1:
                num_items = self.lc.GetItemCount()
                self.lc.InsertItem(num_items, num_items)
                self.lc.SetItem(num_items, 0, dataset_id)
                if isinstance(values, (tuple, list)):
                    if len(values) > 0:
                        self.lc.SetItem(num_items, 1, str(values[0]))
                    if len(values) > 1:
                        self.lc.SetItem(num_items, 2, str(values[1]))
                else:
                    self.lc.SetItem(num_items, 1, str(values))

        self.UpdatePeaks()

        for par_name, par_vals in data.get('par', {}).items():
            val = [par_name, par_vals.get('1', '')]
            if par_vals.get('2', '') == 'fit':
                val.append('fit')
            else:
                val.append('')

            if par_vals.get('3', '') == 'None':
                val.append('')
            else:
                val.append(par_vals.get('3', ''))

            print('adding:', val)
            self.local.append(val)
            print(val)

        for peak_id, values in data.get('peak', {}).items():
            count = self.datasets.GetItemCount()
            col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])

            if peak_id not in col1:
                self.datasets.InsertItem(count, count)
                self.datasets.SetItem(count, 0, peak_id)

            count = self.datasets.GetItemCount()
            col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
            vals = list(values) if isinstance(values, (tuple, list)) else [values]
            for i, c1 in enumerate(col1):
                if c1 == peak_id:
                    print(values)
                    countLC = self.lc.GetItemCount()
                    print('cntLC:', countLC)
                    for ii in range(countLC):
                        if len(vals) >= ii + 1:
                            if vals[ii] == 'True':
                                self.datasets.SetItem(i, ii + 1, str(True))
                            elif vals[ii] == 'False':
                                self.datasets.SetItem(i, ii + 1, str(False))

        print('local', self.local)

        self.SetLocal() #transfer pars to arrays

    def OnButtonSave(self,event):
        print('Saving to:',self.savefile)
        write={}
        
        count = self.lc.GetItemCount()
        col1 = numpy.array([self.lc.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.lc.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.lc.GetItem(row, 2).GetText() for row in range(count)])
        write['dataset']={}
        for i,(c1,c2,c3) in enumerate(zip(col1,col2,col3)):
            write['dataset'][c1]=c2,c3


        count = self.datasets.GetItemCount()
        col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.datasets.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.datasets.GetItem(row, 1).GetText() for row in range(count)])
        write['peak']={}
        for i,(c1,c2,c3) in enumerate(zip(col1,col2,col3)):
            write['peak'][c1]=c2,c3

        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])
        col4 = numpy.array([self.parLocal.GetItem(row, 3).GetText() for row in range(count)])



        write['par']={}
        for i,(c1,c2,c3,c4) in enumerate(zip(col1,col2,col3,col4)):
            write['par'][c1]={}
            write['par'][c1][1]=c2

            if(c3==''):
                write['par'][c1][2]='fix'
            else:
                write['par'][c1][2]=c3

            if(c4==''):
                write['par'][c1][3]='None'
            else:
                write['par'][c1][3]=c4
        write_structured_parameter_file(self.savefile, write)


    def SetupLC(self):

        self.lc.InsertColumn(0,'ID')
        self.lc.InsertColumn(1,'Include')
        self.lc.InsertColumn(2,'path')

        
    def OnRClick(self,event):

        sele=self.lc.GetFirstSelected()
        #print(sele)
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row,0).GetText() for row in range(count)][sele]
        col2 = [self.lc.GetItem(row,1).GetText() for row in range(count)][sele]
        col3 = [self.lc.GetItem(row,2).GetText() for row in range(count)][sele]

        print (col1,col2,col3)
        
        #launch window
        from . import CPMGframe
        import importlib
        cpmgResults=importlib.reload(CPMGframe)
        bool=cpmgResults.CPMGMan(self,pth=col3)

        
        
        
    def OnClick(self,event): #when selecting a peak...plot what we have.

        col=event.GetColumn()
        print(event)
        
        sele=self.lc.GetFirstSelected()
        #print(sele)
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.lc.GetItem(row, 1).GetText() for row in range(count)][sele]
        col3 = [self.lc.GetItem(row, 2).GetText() for row in range(count)][sele]
        print('Current:',col1,col2,col3)
        print('Changing inclusion')
        num_items = self.lc.GetItemCount()

        if(col2=='True'):
            self.lc.SetItem(sele,1,str(False))
        else:
            self.lc.SetItem(sele,1,str(True))
            
        pass





    def OnButtonSort(self,event):

        col=event.GetColumn()

        count = self.lc.GetItemCount()
        col1 = numpy.array([self.lc.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.lc.GetItem(row, 1).GetText() for row in range(count)])

        try:
            col3 = numpy.array([float(self.lc.GetItem(row, 2).GetText()) for row in range(count)])
        except:
            col3 = numpy.array([self.lc.GetItem(row, 2).GetText() for row in range(count)])
        """
        try:
            col4 = numpy.array([float(self.lc.GetItem(row, 3).GetText()) for row in range(count)])
        except:
            col4 = numpy.array([self.lc.GetItem(row, 3).GetText() for row in range(count)])
        #print(numpy.argsort(col2))
        """
        self.lc.ClearAll()

        self.SetupLC()
        
        
        if(col==0):
            s=numpy.argsort(col1)
        elif(col==1):
            s=numpy.flip(numpy.argsort(col2))
        elif(col==2):
            s=numpy.flip(numpy.argsort(col3))
            
            
        for i,arg in enumerate(s):
            #print (arg,col1[arg],col2[arg],col3[arg],col4[arg])
            self.lc.InsertItem(i,i)
            self.lc.SetItem(i,0,str(col1[arg]))
            self.lc.SetItem(i,1,str(col2[arg]))
            self.lc.SetItem(i,2,str(col3[arg]))
            #self.datasets.SetStringItem(i,3,str(col4[arg]))
        
