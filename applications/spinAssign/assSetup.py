#!/usr/bin/python

###################################################################
# Deconvolve nmr spectrum
###################################################################

import numpy,sys,os,string,math,wx, platform
from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin
from assign.assign_main import molecule,spectrum
from importlib import reload
import importlib

import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

#import nmrglue as ng
wx.SYS_COLOUR_BTNTEXT=0
#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)


class assSetup(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent,assignParFile):
    #def __init__(self,uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)


        self.molecule=molecule()

        self.build=platform.uname()[0]
        #self.deconBin='decon_'+self.build


        self.assignParFile=assignParFile
        self.parent=parent

        self.READSEQ=0


        #self.result = wx.StaticText(self, label="")
        #self.result.SetForegroundColour(wx.RED)

        self.buttonAdd = wx.Button(self, label="Add")
        self.buttonRemove = wx.Button(self, label="Remove")
        self.buttonAnalyse = wx.Button(self, label="Analyse")
        self.buttonSave = wx.Button(self, label="Save")
        self.buttonLoad = wx.Button(self, label="Load")
        self.buttonQuit = wx.Button(self, label="Quit")
        self.buttonspin = wx.Button(self,label= "Spins")
        self.buttonReadSeq = wx.Button(self, label="ReadSeq")
        self.buttonReadAssignment = wx.Button(self, label="ReadAss")




        self.Bind(wx.EVT_BUTTON, self.on_spin_button, self.buttonspin)




        self.buttonAdd.Bind(wx.EVT_BUTTON, self.OnButtonAdd)
        self.buttonRemove.Bind(wx.EVT_BUTTON, self.OnButtonRemove)
        self.buttonAnalyse.Bind(wx.EVT_BUTTON, self.OnButtonProcess)
        self.buttonSave.Bind(wx.EVT_BUTTON, self.OnButtonSave)
        self.buttonLoad.Bind(wx.EVT_BUTTON, self.OnButtonLoad)
        self.buttonReadSeq.Bind(wx.EVT_BUTTON, self.OnButtonReadSeq)
        self.buttonReadAssignment.Bind(wx.EVT_BUTTON, self.OnButtonReadAss)
        self.buttonQuit.Bind(wx.EVT_BUTTON, self.OnButtonQuit)
        #self.buttonPeak.Bind(wx.EVT_BUTTON, self.OnButtonPeak)

        self.source = AutoWidthListCtrl(self)
        self.source.SetMinSize((600,200))

        self.source.InsertColumn(0, 'Type', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(1, 'Path', width = 100,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(2, 'Inputfile', width = 200,format=wx.LIST_FORMAT_CENTRE)

        self.source.InsertColumn(3, 'Peaklist',width=200,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(4, 'Ref',width=200,format=wx.LIST_FORMAT_CENTRE)
        self.source.InsertColumn(5, 'f1180',width=200,format=wx.LIST_FORMAT_CENTRE)

        self.seq = AutoWidthListCtrl(self)
        self.seq.SetMinSize((120,200))
        self.seq.InsertColumn(0, 'resi', width = 50,format=wx.LIST_FORMAT_CENTRE)
        self.seq.InsertColumn(1, 'resn', width = 50,format=wx.LIST_FORMAT_CENTRE)



        # Set simple sizer for a nice border
        self.dataLbl = wx.StaticBox(self, -1, 'Data sources:')
        self.dataSizer = wx.StaticBoxSizer(self.dataLbl, wx.VERTICAL)
        self.dataSizer.Add(self.source)

        self.statusLbl = wx.StaticBox(self, -1, 'Status:')
        self.statusSizer = wx.StaticBoxSizer(self.statusLbl, wx.VERTICAL)
        self.statusText=[]
        for i in range(10):
            self.statusText.append(wx.StaticText(self,label=""))
            self.statusSizer.Add(self.statusText[-1])


        self.hbox=wx.BoxSizer(wx.HORIZONTAL)


        self.optionLbl = wx.StaticBox(self, -1, 'Options:')
        self.vboxBut = wx.StaticBoxSizer(self.optionLbl, wx.VERTICAL)
        # self.vboxBut=wx.BoxSizer(wx.VERTICAL)
        self.vboxBut.Add(self.buttonAdd, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonRemove, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonAnalyse, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonSave, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonLoad, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonQuit, 1, wx.TOP, border=5)
        self.vboxBut.Add(self.buttonspin, 1, wx.TOP | wx.BOTTOM, border=5)


        self.hbox.Add(self.vboxBut, 0, wx.LEFT | wx.RIGHT, border=5)
        self.hbox.Add(self.dataSizer, 1, wx.LEFT | wx.RIGHT, border=5)
        self.hbox.Add(self.statusSizer, 1, wx.LEFT | wx.RIGHT, border=5)


        self.molLbl = wx.StaticBox(self, -1, 'Molecule:')
        self.molSizer = wx.StaticBoxSizer(self.molLbl, wx.VERTICAL)


        self.molfileLab = wx.StaticText(self, label="SeqFile:")
        self.molfileBox = wx.TextCtrl(self, size=(200, -1))
        self.molfileBut = wx.Button(self, label="...", size=(30,-1))
        self.molfileBut.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.molfileBox))

        self.msizer = wx.GridBagSizer(15, 3)
        cnt=0
        self.msizer.Add(self.molfileLab,(cnt,0))
        self.msizer.Add(self.molfileBox,(cnt,1))
        self.msizer.Add(self.molfileBut,(cnt,2));cnt+=1

        self.molSizer.Add(self.msizer)
        self.molSizer.Add(self.buttonReadSeq)

        self.msizer1 = wx.GridBagSizer(15, 3)
        cnt=0
        self.assLab = wx.StaticText(self, label="SeqFile:")
        self.assBox = wx.TextCtrl(self, size=(200, -1))
        self.assBut = wx.Button(self, label="...", size=(30,-1))
        self.assBut.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.assBox))
        
        self.msizer1.Add(self.assLab,(cnt,0))
        self.msizer1.Add(self.assBox,(cnt,1))
        self.msizer1.Add(self.assBut,(cnt,2));cnt+=1
        
        self.molSizer.Add(self.msizer1)
        self.molSizer.Add(self.buttonReadAssignment)


        self.parLbl = wx.StaticBox(self, -1, 'Parameters:')
        self.parSizer = wx.StaticBoxSizer(self.parLbl, wx.VERTICAL)
        self.parbag = wx.GridBagSizer(15, 2)
        self.parSizer.Add(self.parbag)
    
        self.canvasSizer=wx.BoxSizer(wx.VERTICAL)
        self.refplot = wx.StaticBox(self, -1, 'Referencing:')
        #self.fig = Figure(figsize=(self.max_size/1.3,96./96.), dpi = 96)
        self.figRef = Figure()
        self.canvasRef = FigCanvas(self, -1, self.figRef)
        self.toolbarRef = NavigationToolbar(self.canvasRef)
        
        self.canvasSizer.Add(self.refplot)
        self.canvasSizer.Add(self.canvasRef, 1, wx.LEFT | wx.TOP | wx.EXPAND| wx.GROW | wx.ALL )
        self.canvasSizer.Add(self.toolbarRef, 0, wx.EXPAND)
       

        self.parsLab={}
        self.parsBox={}
        self.parsBut={}
        self.initVal={}
        self.pars=('tolHNCO','tolHNCA','tolHNCACB','tolHNCANH','tolMatch','tolSin','tolMax','tolMin','refSpec','template')
        self.initVal['tolHNCO']=0.02
        self.initVal['tolHNCA']=0.4
        self.initVal['tolHNCACB']=0.3
        self.initVal['tolHNCANH'] = 0.4
        self.initVal['tolMatch']=0.4
        self.initVal['tolSin']=5
        self.initVal['tolMax']=0.6
        self.initVal['tolMin']=0.1
        self.initVal['refSpec']='hnco'
        self.initVal['template']=''
        for i,par in enumerate(self.pars):
            self.parsLab[par]= wx.StaticText(self, label=par+':')
            self.parsBox[par]= wx.TextCtrl(self, size=(200, -1))
            self.parsBox[par].SetValue(str(self.initVal[par]))
            self.parbag.Add(self.parsLab[par],(i,0))
            self.parbag.Add(self.parsBox[par],(i,1))
            if(par=='template'):
                self.parsBut[par]=wx.Button(self, label="...", size=(30,-1))
                self.parsBut[par].Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.parsBox[par]))
                self.parbag.Add(self.parsBut[par],(i,2))
        #self.magmaLbl = wx.StaticBox(self, -1, 'Magma:')
        #self.magmaSizer = wx.StaticBoxSizer(self.magmaLbl, wx.VERTICAL)
        #self.magmaBut = wx.Button(self, label="RunMagma", size=(30,-1))
        #self.magmaBut.Bind(wx.EVT_BUTTON,self.on_run_magma)




        #self.magmaSizer.Add(self.magmaBut)

        self.vboxFull=wx.BoxSizer(wx.VERTICAL)
        self.vboxFull.Add(self.hbox)


        self.vboxM=wx.BoxSizer(wx.VERTICAL)
        self.vboxM.Add(self.molSizer)
        self.vboxM.Add(self.seq)

        self.vboxH=wx.BoxSizer(wx.HORIZONTAL)
        self.vboxH.Add(self.vboxM)
        self.vboxH.Add(self.parSizer)
        self.vboxH.Add(self.canvasSizer)

        self.vboxFull.Add(self.vboxH)

        self.SetSizerAndFit(self.vboxFull)

        self.Status()

    """
    def on_run_magma(self,event):
        #self.molecule.normSpec()
        #self.molecule.assSpec()
        #self.molecule.PeakShifts()
        #self.molecule.EdgeScreen(strict='n')

        #self.molecule.WriteInit('results/1/mces.txt')
        self.molecule.WriteInit()


        os.system('rm results/mces.txt.G')
        self.Nprocessors=1
        args=' results/ full'
        self.mcesPath='mcesCore_Darwin'
        self.mcesPathParallel='mcesCore_parallel_Darwin'
        if(self.Nprocessors>1): #run in parallel
            runline='mpirun -np '+str(self.Nprocessors)+' '+self.mcesPathParallel+' '+args
            #runline='mpirun --mca pml ob1  '+self.mcesPathParallel+' '+args
        else:
            runline=self.mcesPath+' '+args
        print(runline)
        os.system(runline)
    """


    def OnButtonAdd(self,event):
        bool=AddBoxMan(self)

    def OnButtonRemove(self,event):
        sele=self.source.GetFirstSelected()
        print('Removing item',sele)
        self.source.DeleteItem(sele)


        ###############################################


        #self.peakSizer.SetMinSize((widthToSetDecn, 0))
        #self.panel1.SetMinSize((0, heightToSetDecn))
        #self.deconSizer.SetMinSize((0,heightToSetDecn))
        #self.statusSizer.SetMinSize((0,heightToSetDecn))

        #indir=Parse(self.deconParFile,'indir')
        #thresh=ParseFlt(self.deconParFile,'thresh')
        #ncpus=int(Parse(self.deconParFile,'ncpus'))
        #fac=ParseFlt(self.deconParFile,'fac')
        #squash=ParseFlt(self.deconParFile,'squash')
        #maxiter=ParseFlt(self.deconParFile,'maxiter')
        #infile=Parse(self.deconParFile,'infile')
        #peakfile=Parse(self.deconParFile,'peakfile')


        #self.peakBox.SetValue(str(peakfile))



    def Status(self):
        for i in range(10): #reset status
            self.statusText[i].SetLabel("")

        cnt=0
        if(self.READSEQ==0):
            self.statusText[cnt].SetLabel("No sequence loaded")
        else:
            self.statusText[cnt].SetLabel("Residues in sequence: %i" % len(self.molecule.seq.keys()));cnt+=1
            self.statusText[cnt].SetLabel("Expected peaks: %i" % len(self.molecule.seq.keys()))
            ros=0
            for i,resi in enumerate(self.molecule.resi):
                #print(resi,self.molecule.seq[resi])
                if(self.molecule.seq[resi]!='P'):
                    ros+=1
            self.statusText[cnt].SetLabel("Expected peaks (inc N term): %i" % ros)
        cnt+=1


        specload=len(self.molecule.spec.keys())
        if(specload==0):
            self.statusText[cnt].SetLabel("No spectra loaded")
        else:
            self.statusText[cnt].SetLabel("Spectra loaded: %i" % specload)
        cnt+=1

        pcnt=0
        for pk in self.molecule.peak.keys():
            for spec in self.molecule.peak[pk].keys():
                pcnt+=len(self.molecule.peak[pk][spec])
        if(pcnt==0):
            self.statusText[cnt].SetLabel("No 3D peaks")
        else:
            self.statusText[cnt].SetLabel("2D peaks: %i" % len(self.molecule.peak.keys()));cnt+=1
            self.statusText[cnt].SetLabel("3D peaks: %i" % pcnt)
        cnt+=1


        self.SetSizerAndFit(self.vboxFull)

        #for i in range(10):
        #    self.updateList[i].SetLabel("shit")



        """
        cnt=0
        if(self.READ==0):
            self.updateList[cnt].SetLabel("No spectrum in memory")
        else:
            #self.updateList[cnt].SetLabel("Read in %s" % (self.spectrumFile))
            self.updateList[cnt].SetLabel("Spectrum Dimensions: %s" % self.dim);cnt+=1
            print(self.specsize,self.dim)
            liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim1:",self.labb[0],self.uc0min,self.uc0max,self.specsize[0])
            self.updateList[cnt].SetLabel(liney);cnt+=1
            liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim2:",self.labb[1],self.uc1min,self.uc1max,self.specsize[1])
            self.updateList[cnt].SetLabel(liney);cnt+=1

            if(self.dim>=3):
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim3:",self.labb[2],self.uc2min,self.uc2max,self.specsize[2])
                self.updateList[cnt].SetLabel(liney);cnt+=1
                if(self.dim==3):
                    self.updateList[cnt].SetLabel("will project down %s" % self.labb[0]);cnt+=1
                    self.updateList[cnt].SetLabel("peak list must be %s:%s" % (self.labb[1],self.labb[2]));cnt+=1

            if(self.dim>=4):
                liney="%s %s %.2f to %.2f ppm (%i pts)" % ("dim4:",self.labb[3],self.uc3min,self.uc3max,self.specsize[3])
                self.updateList[cnt].SetLabel(liney);cnt+=1
                if(self.dim==4):
                    self.updateList[cnt].SetLabel("will project down %s:%s" % (self.labb[0],self.labb[1]));cnt+=1
                    self.updateList[cnt].SetLabel("peak list must be %s:%s" % (self.labb[2],self.labb[3]));cnt+=1

        #cnt+=1
        if(len(self.peak)==0):
            self.updateList[cnt].SetLabel("No peaks in projection")
        else:
            self.updateList[cnt].SetLabel("ProjectedPeaks: "+str(len(self.peak)))
        cnt+=1


        if(self.DECON==0):
            self.updateList[cnt].SetLabel("No deconvolved spectrum in memory");cnt+=1
        else:
            self.updateList[cnt].SetLabel("Deconvolved spectrum loaded");cnt+=1

        if(len(self.conn_data)==0):
            self.updateList[cnt].SetLabel("No picked peaks");cnt+=1
        else:
            self.updateList[cnt].SetLabel("Total cross peaks: "+str(len(self.conn_data)));cnt+=1
            if(self.cb_grid.IsChecked()==True):
                self.updateList[cnt].SetLabel("Diagonals: "+str(self.diagCnt));cnt+=1
                self.updateList[cnt].SetLabel("Cross:     "+str(self.crossCnt));cnt+=1

        for i in range(len(self.updateList)-cnt):
            self.updateList[i+cnt].SetLabel("")


        self.SetSizerAndFit(self.fullSizer)
        """


    #FGA added
    def onGetFile(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print(path)
            #fu=self.specPathBox.GetValue()
            #print(fu)
            #print(path.split(fu))
            splitPath = path.split(cwd)
            textBox.SetValue('.' + splitPath[1])
            print("You chose the following file(s):")
            print(path)
        dlg.Destroy()

    #FGA added
    def onGetDir(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.DirDialog(self, message="Choose a folder",         style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            try:
                textBox.SetValue('.' + splitPath[1])
            except:
                textBox.SetValue(path)
            print("You chose the following file(s):")
            print(path)
            #os.chdir(path)
            #self.dirBox.SetValue(path)
            #print("CWD: ",os.getcwd())

        dlg.Destroy()




    def OnButtonProcess(self,event):



        for par in self.pars: #set variable name in molecule to value in textbox
            if(par != 'template' and par!='refSpec'):
                exec('self.molecule.'+par+'=float("%s")' % str(self.parsBox[par].GetValue()))
            else:
                exec('self.molecule.'+par+'="%s"' % str(self.parsBox[par].GetValue()))
        print('tolHNCO:',self.molecule.tolHNCO)
        print('template:',self.molecule.template)


        self.OnButtonReadSeq(True)


        num_items = self.source.GetItemCount()
        for cnt in range(num_items):
            a=self.source.GetItem(cnt,0).GetText()
            b=self.source.GetItem(cnt,1).GetText()
            c=self.source.GetItem(cnt,2).GetText()
            d=self.source.GetItem(cnt,3).GetText()
            e=self.source.GetItem(cnt,4).GetText()
            f=self.source.GetItem(cnt,5).GetText()
            #read spectrum
            self.molecule.AddSpec(spectrum(a,b,c,d,e,f))


        # check that all the peaks are in the refSpec!
        specs = self.molecule.spec.keys()  #all the spec names
        self.molecule.reference()  #align spectra according to settings and rules
  
        self.molecule.normSpec()   #classify peaks according to current positions
        self.molecule.assSpec()    #make graph of connections
        self.molecule.EdgeScreen()   #analyse graph for self consistency

        self.molecule.PeakShifts() #analyse chemical shifts to create distribuions
      
   
        #print(self.molecule.G1edges)

        #print(self.molecule.G1edges)
        #sys.exit(100)


        self.molecule.WriteInit()
        self.parent.AddTabTwo()   #draw figures
        self.parent.AddTabThree() #draw figures
        self.parent.AddTabFour()
        ccctocsy_included = False
        inny = open(self.assignParFile, 'r')
        for line in inny.readlines():
            if 'ccctocsy' in line:
                ccctocsy_included=True
                break
        #if(ccctocsy_included==True):
        #    self.parent.AddTabFive()
        self.Status()


    def OnButtonReadSeq(self,event):
        seqfil=self.molfileBox.GetValue()
        if(os.path.exists(seqfil)==0):
            print('Sequence file does not exist:',seqfil)
            return
        self.molecule.GetSeq(seqfil)
        cnt=0
        for i in range(len(self.molecule.resi)):
            cnt+=1
            #print('current num itesm:',num_items)
            resi=self.molecule.resi[i]
            resn=self.molecule.seq[resi]
            num_items = self.seq.GetItemCount()
            self.seq.InsertItem(num_items,str(cnt))
            self.seq.SetItem(num_items,0,str(resi))
            self.seq.SetItem(num_items,1,str(resn))

        self.READSEQ=1
    
    def OnButtonReadAss(self,event):
        assFile=self.assBox.GetValue()
        results={}
        if(os.path.exists(assFile)==0):
            print('Sequence file does not exist:',assFile)
            return
        inny = open(assFile,'r')
        for line in inny.readlines():
            line=line.split('\n')[0]
            if('residue' in line):
                continue
            line=line.split(',')
            if(line[1]==''):
                continue
            else:
                results[line[1]+'H-N']=line[0]+line[2]
        
        self.molecule.assemble_assignment('assignment', results)
        

    def on_spin_button(self, event):
        # import spinFrame
        #from assign.spinFrame import spinFrame
        from . import spinFrame
        spinFrame=importlib.reload(spinFrame)
        # spinFrame=reload(spinFrame)
        bool=spinFrame.spinFrame(self)
        pass


    def OnButtonReadPeak(self,event):
        if(self.READ==0):
            print('No data. Trying to read that in first.')
            self.OnButtonRead(True)
            if(self.READ==0):
                print('Cannot read in data either.')
                return
            else:
                print('Successfully read in data.')
                print('Continuing...')
        self.ReadPeakListFile()



    def DeletePage(self,pageTitle):
        #pageTitle='2Dplanes'
        for index in range(self.parent.GetPageCount()):
            if self.parent.GetPageText(index) == pageTitle:
                self.parent.DeletePage(index)
                self.parent.SendSizeEvent()
                break
        #self.parent.DeletePage('2Dplanes')



    def OnButtonQuit(self,event):
        print('exiting')
        sys.exit(100)

    def OnButtonPeak(self,event):
        import peakAssFrame
        spinFrame=reload(peakAssFrame)
        bool=peakAssFrame.peakAssFrame(self)


    def OnButtonSave(self,event):
        if os.path.exists("assignParFile"):
            dlg = wx.MessageDialog(None, "Do you want to update?",'Updater',wx.YES_NO | wx.ICON_QUESTION)
            result = dlg.ShowModal()

            if result == wx.ID_YES:
                self.performSave()
            else:
                print("No pressed")


    def performSave(self):
        print('Saving to:','assignParFile')
        outy=open('assignParFile','w')
        num_items = self.source.GetItemCount()
        for cnt in range(num_items):
            a=self.source.GetItem(cnt,0).GetText()
            b=self.source.GetItem(cnt,1).GetText()
            c=self.source.GetItem(cnt,2).GetText()
            d=self.source.GetItem(cnt,3).GetText()
            e=self.source.GetItem(cnt,4).GetText()
            f=self.source.GetItem(cnt,5).GetText()
            outy.write('SPECTRUM: %i %s %s %s %s %s %s\n' % (cnt,a,b,c,d,e,f))
        outy.write('MOLFILE: %s\n' % (self.molfileBox.GetValue()))

        for i,par in enumerate(self.pars):
            outy.write('%s %s\n' % (par,self.parsBox[par].GetValue()))
        outy.close()

    def OnButtonLoad(self,event):
        #reset box



        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        if(os.path.exists('assignParfile')==0):
            print('Cannot find ','assignParFile')
            print('Aborting Load')
            return
        cnt=0
        inny=open('assignParFile')
        #self.source.ClearAll()
        for line in inny.readlines():
            test=line.split()
            if(len(test)!=0):
                if(test[0]=='SPECTRUM:' and len(test)==8):
                    cnt+=1
                    num_items = self.source.GetItemCount()
                    print('current num items:',num_items)
                    self.source.InsertItem(num_items,str(cnt))
                    self.source.SetItem(num_items,0,test[2])
                    self.source.SetItem(num_items,1,test[3])
                    self.source.SetItem(num_items,2,test[4])
                    self.source.SetItem(num_items,3,test[5])
                    self.source.SetItem(num_items,4,test[6])
                    self.source.SetItem(num_items,5,test[7])
                if(test[0]=='MOLFILE:' and len(test)==2):
                    self.molfileBox.SetValue(test[1])
                if(test[0] in self.pars):
                    if(len(test)==2):
                        self.parsBox[test[0]].SetValue(test[1])
                    else:
                        self.parsBox[test[0]].SetValue('')
                inny.close()
        print('Loaded:',cnt,'spectra')
        self.OnButtonProcess(True)
        self.DrawFigure()

        """
        write={}
        write['indir']=self.dirBox.GetValue()
        write['infile']=self.infileBox.GetValue()
        write['peakfile']=self.peakBox.GetValue()
        write['dim']=str(self.dimBox.GetSelection()+1)

        write['sig1']=self.sig1Box.GetValue()
        write['sig2']=self.sig2Box.GetValue()
        write['sig3']=self.sig3Box.GetValue()
        if(self.dim==4):
            write['sig4']=self.sig4Box.GetValue()

        write['thresh']=self.threshBox.GetValue()
        write['ncpus']=self.coreBox.GetValue()
        write['fac']=self.facBox.GetValue()
        write['squash']=self.squashBox.GetValue()
        write['maxiter']=self.maxiterBox.GetValue()
        write['symmode']=self.IntToBool(self.cb_grid.IsChecked())

        dec=[]
        print(self.deconParFile)
        if(os.path.exists(self.deconParFile)):
            inny=open(self.deconParFile)
            for line in inny.readlines():
               dec.append(line)
            inny.close()

        outy=open(self.deconParFile,'w')
        for de in dec:
            test=de.split()
            if(len(test)>0):
                if(test[0] in write.keys()):
                    outy.write('%s = %s ' % (test[0],write[test[0]]))
                    for j in range(len(test)-3):
                        outy.write(' %s' % test[j+3])
                    outy.write('\n')
                    del write[test[0]]
                else:
                    outy.write(de)
            else:
                outy.write(de)

        for key,vals in write.items():
            try:
                outy.write('%s = %s\n' % (key,vals[0]))
            except:
                print('problem with ',key)

        outy.close()
        """



    def DrawFigure(self):
        self.molecule.MakeHistogram()


        self.axRef1 = self.figRef.add_subplot(211)
        self.axRef2 = self.figRef.add_subplot(212)
        #self.axRef3 = self.figRef.add_subplot(111)

        maxo=0
        for spec in self.molecule.hist.keys():
            if(spec=='hcconh' or spec=='Hhsqc'):
                continue
            edges,gram=self.molecule.hist[spec]
            #print (spec,edges,gram)
            self.axRef1.plot(edges,gram,label=spec,zorder=2) 
            if(numpy.max(gram)>maxo):
                maxo=numpy.max(gram)

        maxB=0 #get max value in BMRB histograms
        for spec in self.molecule.bmrbHist.keys():
            if(spec[0]!='C'):
                continue
            edges,gram=self.molecule.bmrbHist[spec]
            if(numpy.max(gram)>maxB):
                maxB=numpy.max(gram)
        for spec in self.molecule.bmrbHist.keys():
            if(spec[0]!='C'):
                continue
            edges,gram=self.molecule.bmrbHist[spec]
            #print (spec,edges,gram)
            self.axRef1.plot(edges,gram/maxB*maxo,color='k',zorder=0) 

        #Gender=['maxLim','minLim','i(ex)','i-1(ex)','i','i-1']
        self.axRef1.legend(loc=2)
        self.axRef1.set_xlabel("13C (ppm)",fontsize=10)


        maxo=0
        for spec in self.molecule.hist.keys():
            if(spec!='hcconh' and spec!='Hhsqc'):
                continue
            edges,gram=self.molecule.hist[spec]
            #print (spec,edges,gram)
            self.axRef2.plot(edges,gram,label=spec,zorder=2) 
            if(numpy.max(gram)>maxo):
                maxo=numpy.max(gram)

        maxB=0 #get max value in BMRB histograms
        for spec in self.molecule.bmrbHist.keys():
            if(spec[0]!='H'):
                continue
            edges,gram=self.molecule.bmrbHist[spec]
            if(numpy.max(gram)>maxB):
                maxB=numpy.max(gram)
        for spec in self.molecule.bmrbHist.keys():
            if(spec[0]!='H'):
                continue
            edges,gram=self.molecule.bmrbHist[spec]
            #print (spec,edges,gram)
            self.axRef2.plot(edges,gram/maxB*maxo,color='k',zorder=0) 

        self.axRef2.legend(loc=2)
        self.axRef2.set_xlabel("1H (ppm)",fontsize=10)



        self.figRef.tight_layout()
        self.canvasRef.draw()
        #sys.exit(100)        

        #print(self.molecule.shufty)
        #self.GetShift(pk)
        #sys.exit(100)



    def OnButtonAnalyse(self,event):
        self.dim=(self.dimBox.GetSelection()+1)
        self.DECON=0         #decon flag
        self.pkSlice1Ddec=[] #1D slices
        self.conn_data=[]
        self.noeTags=[]

        try:
            poll = self.calcy.poll()
            if poll == None:
                print('Calculation still running in background')
        except:
            pass

        self.corrFile='out/correlate.3'
        print('Reading outputs from',self.corrFile)

        if(self.dim==3):
            #indir=self.dirBox.GetValue()
            infile=self.infileBox.GetValue()
            if(os.path.exists(infile+'.decon')==1):
                self.dicdec,self.datadec=ng.pipe.read(infile+'.decon')
                self.DECON=1
                if(self.datadec.shape!=self.data.shape):
                    print('deconvolved spectrum is a different shape.')
                    print('recalculate the deconvolution')
                    self.DECON=0
                    #numpy.delete(self.datadec)
                    #numpy.delete(self.dicdec)
                    return
                if(self.DECON==1):
                    print('Deconvolved spectrum in memory:',infile)
            else:
                print('Deconvolution file does not exist.')
                return
            self.pkSlice1Ddec=[] #1D slices
            for pkl in range(len(self.peak)):
                ptC=self.pkIdx[pkl][0]
                ptH=self.pkIdx[pkl][1]
                self.pkSlice1Ddec.append(self.datadec[:,ptC,ptH])

        else:
            pass


        self.conn_data=[]
        if(self.cb_grid.IsChecked()==True):
            sym='y'
        else:
            sym='n'

        print('Reading in ',self.corrFile)

        self.GetConn(self.corrFile,sym=sym)


        print('Done!')






class AddBoxMan(wx.App):
    def __init__(self,inherit):
        self.frame_addFrame=AddBox(inherit)
        self.frame_addFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]

class AddBox(wx.Frame):

    def __init__(self,parent):
        #wx.Panel.__init__(self, parent=parent)
        #self.parent=parent
        #self.tabOne=parent.tabOne
        #self.create_main_panel()
        #self.draw_figure()
        #self.canvas.draw()

        self.parent=parent

        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title=u'Add spectrum...')
        self.SetClientSize(wx.Size(900, 280))


        panel=wx.Panel(self,-1)

        self.typeyLab = wx.StaticText(self, label="SpectrumType:")
        self.listy=[]
        self.listy.append('hnco')
        self.listy.append('hncaco')
        self.listy.append('hnca')
        self.listy.append('hncoca')
        self.listy.append('hncacb')
        self.listy.append('hncocacb')
        self.listy.append('ctocsy')
        self.typey=wx.ComboBox(self, -1,choices=self.listy, style=wx.CB_READONLY)

        self.specPathLab = wx.StaticText(self, label="SpectrumPath:")
        self.specPathBox = wx.TextCtrl(self, size=(200, -1))
        self.specPathBut = wx.Button(self, label="...", size=(30,-1))
        self.specPathBut.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.specPathBox))

        self.deconLab = wx.StaticText(self, label="ParameterFile:")
        self.deconBox = wx.TextCtrl(self, size=(200, -1))
        self.deconBut = wx.Button(self, label="...", size=(30,-1))
        self.deconBut.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.deconBox))

        self.peakLab = wx.StaticText(self, label="PeakList:")
        self.peakBox = wx.TextCtrl(self, size=(200, -1))
        self.peakBut = wx.Button(self, label="...", size=(30,-1))
        self.peakBut.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.peakBox))

        self.addBut = wx.Button(self, label="Add", size=(150,-1))
        self.addBut.Bind(wx.EVT_BUTTON,self.onButtonAdd)

        self.cancelBut = wx.Button(self, label="Cancel", size=(150,-1))
        self.cancelBut.Bind(wx.EVT_BUTTON,self.onButtonCancel)



        self.sizer = wx.GridBagSizer(15, 3)
        cnt=0
        self.sizer.Add(self.typeyLab,(cnt,0))
        self.sizer.Add(self.typey,(cnt,1));cnt+=1

        self.sizer.Add(self.specPathLab,(cnt,0))
        self.sizer.Add(self.specPathBox,(cnt,1))
        self.sizer.Add(self.specPathBut,(cnt,2));cnt+=1

        self.sizer.Add(self.deconLab,(cnt,0))
        self.sizer.Add(self.deconBox,(cnt,1))
        self.sizer.Add(self.deconBut,(cnt,2));cnt+=1

        self.sizer.Add(self.peakLab,(cnt,0))
        self.sizer.Add(self.peakBox,(cnt,1))
        self.sizer.Add(self.peakBut,(cnt,2));cnt+=1

        self.vbox=wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.sizer)
        self.vbox.Add(self.addBut)
        self.vbox.Add(self.cancelBut)
        self.SetSizerAndFit(self.vbox)

        self.Show(True)


    def onButtonAdd(self,event):
        print('Adding....')
        typey=self.listy[self.typey.GetSelection()]
        specpath=self.specPathBox.GetValue()
        parfile=self.deconBox.GetValue()
        peakfile=self.peakBox.GetValue()
        print(specpath,parfile,peakfile)
        for test in (specpath,parfile,peakfile):
            if(os.path.exists(test)==0):
                print('WARNING: cannot find file: ',test)

        num_items = self.parent.source.GetItemCount()
        self.parent.source.InsertItem(num_items,str(num_items))
        self.parent.source.SetItem(num_items,0,typey)
        self.parent.source.SetItem(num_items,1,specpath)
        self.parent.source.SetItem(num_items,2,parfile)
        self.parent.source.SetItem(num_items,3,peakfile)
        self.Close()
        pass

    def onButtonCancel(self,event):
        print('Cancelling')
        self.Close()

    #FGA added
    def onGetFile(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        cwd=self.specPathBox.GetValue()
        #cwd=self.dirBox.GetValue()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=cwd, defaultFile="",
            wildcard="PDB file (decon*)|decon*|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print('p',path)
            fu=self.specPathBox.GetValue().split('.')[1]
            print(fu)
            print(path.split(fu))
            splitPath = path.split(fu)
            print(splitPath)
            textBox.SetValue('.' + splitPath[1])
            print("You chose the following file(s):")
            print(path)

            """
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            textBox.SetValue('.' + splitPath[1])
            print("You chose the following file(s):")
            print(path)
            """
        dlg.Destroy()

    def onGetDir(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.DirDialog(self, message="Choose a folder",         style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            try:
                textBox.SetValue('.' + splitPath[1])
            except:
                textBox.SetValue(path)
            print("You chose the following file(s):")
            print(path)
            #os.chdir(path)
            #self.dirBox.SetValue(path)
            #print("CWD: ",os.getcwd())

        dlg.Destroy()
