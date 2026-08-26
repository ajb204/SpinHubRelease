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

###################################################################
# magmaFrame
###################################################################
import numpy,sys,os,string,math,wx
#from processFrame import ProcMan
import copy,time,re

# from magma.magma import Magma
# from magma.report import Report
# from magma.nofit import Visualise

from subprocess import Popen
import subprocess
import signal
# from spinDecon.integrations.magma import results as magmaResults #bring in the magmaResults class
from spinDecon.gui.dialogs import text_viewer as textEdit
from spinDecon.project.parameter_store import update_parameter_file

#from LoadFilePopUp import LoadFilePopUp

#from threading import Thread
from multiprocessing import Process
#import imp
import importlib


class TestThread(Process):
    def __init__(self,infile):
        """Init Worker Thread Class."""
        Process.__init__(self)
        print(infile)
        self.infile=infile
        self.daemon=True
        self.start()    # start the thread


    #----------------------------------------------------------------------
    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread.
        from magma.magma import Magma
        Magma(self.infile)


        """
        for i in range(6):
            time.sleep(10)
            wx.CallAfter(self.postTime, i)
        time.sleep(1)
        #wx.CallAfter(Publisher.sendMessage, "update", "Thread finished!")
        wx.CallAfter( "update", "Thread finished!")
        """



class magmaFrame(wx.ScrolledWindow):
    """ The main frame of the application
    """
    title = 'MAGMA'

    def OnQuit(self, e):
        self.parent.parent.Close()
    def on_about(self, event):
        msg="fun with methyl assignments"
        dlg = wx.MessageDialog(self, msg, "MAGMA", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()
        m_load = menu_file.Append(-1, "&Open\tCtrl-L", "Open session file")
        #self.Bind(wx.EVT_MENU, self.OnLoadResults, m_load)
        self.parent.parent.Bind(wx.EVT_MENU, self.OnLoadResults, m_load)
        menu_file.AppendSeparator()

        m_save = menu_file.Append(-1, "&Save\tCtrl-S", "Save status")
        self.parent.parent.Bind(wx.EVT_MENU, self.OnSaveResults, m_save)
        #self.Bind(wx.EVT_MENU, self.on_save_plot, m_save)
        menu_file.AppendSeparator()


        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.parent.parent.Bind(wx.EVT_MENU, self.OnQuit, m_exit)
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.parent.parent.Bind(wx.EVT_MENU, self.on_about, m_about)
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.parent.parent.SetMenuBar(self.menubar)



    def __init__(self,parent,deconParFile):
        self.parent=parent
        self.state = getattr(parent, "state", None)
        self.pars={}
        self.WXV=int(wx.__version__.split('.')[0])

        self.DECONFULL=1 #include decon conversion

        """
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        newitem = wx.MenuItem(fileMenu,wx.ID_NEW,'&Quit\tCtrl+Q')
        #qmi.SetBitmap(wx.Bitmap('exit.png'))
        fileMenu.Append(newitem)

        #self.Bind(wx.EVT_MENU, self.OnQuit,1)

        menubar.Append(fileMenu, '&File')
        parent.parent.SetMenuBar(menubar)

        self.SetSize((350, 250))
        self.SetTitle('Icons and shortcuts')
        self.Centre()
        """
        # self.create_menu()



        self.magmaCalc=0
        self.deconParFile=os.path.join(os.getcwd(),deconParFile)
        self.GetPars() #read pars from deconPar

        #FGA changed from wx.Panel
        wx.ScrolledWindow.__init__(self, parent=parent, style=wx.VSCROLL|wx.HSCROLL)
        self.SetScrollRate( 5, 5 )
        #wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.SetBackgroundColour('WHITE')

        #color = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
        #self.SetBackgroundColour(color)




        # Set sizer for the panel content
        self.splitSizer = wx.BoxSizer(wx.HORIZONTAL)


        self.sizerP = wx.GridBagSizer(15, 3)
        self.sizerN = wx.GridBagSizer(15, 3)
        self.sizerM = wx.GridBagSizer(15, 3)

        #self.sizerR = wx.GridBagSizer(15, 3)
        self.sizerR = wx.BoxSizer(wx.VERTICAL)
        ###########################
        #PDB info
        self.dirLbl  = wx.StaticText(self, label="workDir:")
        self.dirVal     = wx.TextCtrl(self,   size=(150, -1))
        self.openDirFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.openDirFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.dirVal))


        self.pdbfileLbl  = wx.StaticText(self, label="PDB file:")
        self.pdbfile     = wx.TextCtrl(self,   size=(150, -1))
        self.pdbfile.SetValue(str(self.pars['pdb_file']))
        self.openPdbFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.openPdbFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.pdbfile))


        self.chainLbl  = wx.StaticText(self, label="Chains:")
        self.chain     = wx.TextCtrl(self,   size=(200, -1))

        #chainTmp=str(self.pars['chains'])
        #print chainTmp
        #chainVals=''
        #for i in range(len(chainTmp)):
        #    chainVals+=chainTmp[i]+' '
        #print chainVals
        #self.chain.SetValue(chainVals)
        self.chain.SetValue(self.pars['chains'])

        self.cb_LV = wx.CheckBox(self, -1,"MergeLVs",style=wx.ALIGN_RIGHT)
        if(self.pars['mergeLV']=="on"):
            self.cb_LV.SetValue(True)


        self.methylLbl  = wx.StaticText(self, label="Methyls:")
        self.methyls=[]
        self.methyls.append('ILE')
        self.methyls.append('LEU')
        self.methyls.append('VAL')
        self.methyls.append('MET')
        self.methyls.append('ALA')
        self.methyls.append('THR')
        self.methylBox = wx.CheckListBox(self,-1, choices=self.methyls,size=(150,127))
        self.SyncBox(self.methylBox,self.methyls,'residues')
        self.methyls=numpy.array(self.methyls)


        #self.LVLbl  = wx.StaticText(self, label="LVenantiomers:")
        self.LVs=[]
        self.LVs.append('off')
        self.LVs.append('R')
        self.LVs.append('S')
        self.LVs.append('merge')

        self.LVBox = wx.RadioBox(self, label = 'LV enantiomers', size=(150,-1),pos = wx.DefaultPosition, choices = self.LVs, majorDimension = 1)
        length = self.LVBox.GetSize()[1]
        # self.LVBox.SetMinSize((200, length))
        self.LVBox.Bind(wx.EVT_RADIOBOX,self.OnRadioBox)
        # self.LVBox.SetSelection(choice)


        #self.LVBox = wx.CheckListBox(self, -1, choices=self.LVs)
        #self.SyncBox(self.LVBox,self.LVs,'proRS')
        self.LVs=numpy.array(self.LVs)

        self.threshLbl  = wx.StaticText(self, label="Distance:")
        self.thresh     = wx.TextCtrl(self,   size=(200, -1))
        self.thresh.SetValue(str(self.pars['short_distance_threshold']))

        self.longthreshLbl  = wx.StaticText(self, label="Long:")
        self.longthresh     = wx.TextCtrl(self,   size=(200, -1))
        self.longthresh.SetValue(str(self.pars['long_distance_threshold']))

        self.shortthreshLbl  = wx.StaticText(self, label="Short:")
        self.shortthresh     = wx.TextCtrl(self,   size=(200, -1))
        self.shortthresh.SetValue(str(self.pars['shortLim']))


        self.buttonShowP = wx.Button(self, label="ShowPDB")
        #self.buttonSave = wx.Button(self, label="Save")
        self.Bind(wx.EVT_BUTTON, self.OnButtonShowP, self.buttonShowP)
        #self.Bind(wx.EVT_BUTTON, self.OnButtonSave, self.buttonSave)

        cnt=0

        self.sizerP.Add(self.dirLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.dirVal   ,(cnt,1),flag=wx.EXPAND);
        self.sizerP.Add(self.openDirFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.pdbfileLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.pdbfile   ,(cnt,1),flag=wx.EXPAND);
        self.sizerP.Add(self.openPdbFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.chainLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.chain   ,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.methylLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.methylBox   ,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.cb_LV   ,(cnt,1),flag=wx.EXPAND);cnt+=1

        #self.sizerP.Add(self.LVLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.LVBox   ,(cnt,1),flag=wx.EXPAND);cnt+=1

        self.sizerP.Add(self.threshLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.thresh      ,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.shortthreshLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.shortthresh      ,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerP.Add(self.longthreshLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerP.Add(self.longthresh      ,(cnt,1),flag=wx.EXPAND);cnt+=1


        self.sizerP.Add(self.buttonShowP   ,(cnt,1),flag=wx.EXPAND);cnt+=1
        #self.sizerP.Add(self.buttonSave   ,(cnt,1),flag=wx.EXPAND);cnt+=1

        #############################


        if(self.DECONFULL==1):
            #Decon info
            self.deconfileLbl  = wx.StaticText(self, label="Decon file:")
            self.deconfile     = wx.TextCtrl(self,   size=(150, -1))
            self.deconfile.SetValue('out/correlate.3')

            #FGA added
            self.selectDeconFileBtn = wx.Button(self, label="...", size=(30,-1))
            self.selectDeconFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.deconfile))


            self.screenLbl  = wx.StaticText(self, label="SN lim:")
            self.screenbox  = wx.TextCtrl(self,   size=(200, -1))
            self.screenLbl2  = wx.StaticText(self, label="Hlim:")
            self.screenbox2  = wx.TextCtrl(self,   size=(200, -1))
            self.screenLbl3  = wx.StaticText(self, label="Clim:")
            self.screenbox3  = wx.TextCtrl(self,   size=(200, -1))
            #self.screenLbl4  = wx.StaticText(self, label="Rat:")
            #self.screenbox4  = wx.TextCtrl(self,   size=(200, -1))
            self.screenbox.SetValue('0.001')
            self.screenbox2.SetValue('0.01')
            self.screenbox3.SetValue('0.1')
            #self.screenbox4.SetValue('0.8')
            self.buttonScreen = wx.Button(self, label="Screen")
            self.cb_recip = wx.CheckBox(self, -1,"Reciprocate",style=wx.ALIGN_RIGHT)


            cnt=0
            self.sizerD = wx.GridBagSizer(15, 3)
            self.sizerD.Add(self.deconfileLbl,(cnt,0),flag=wx.EXPAND)
            self.sizerD.Add(self.deconfile   ,(cnt,1),flag=wx.EXPAND)
            self.sizerD.Add(self.selectDeconFileBtn, (cnt,2), flag=wx.EXPAND);cnt+=1
            self.sizerD.Add(self.screenLbl,(cnt,0),flag=wx.EXPAND)
            self.sizerD.Add(self.screenbox   ,(cnt,1),flag=wx.EXPAND);cnt+=1
            self.sizerD.Add(self.screenLbl2,(cnt,0),flag=wx.EXPAND)
            self.sizerD.Add(self.screenbox2   ,(cnt,1),flag=wx.EXPAND);cnt+=1
            self.sizerD.Add(self.screenLbl3,(cnt,0),flag=wx.EXPAND)
            self.sizerD.Add(self.screenbox3   ,(cnt,1),flag=wx.EXPAND);cnt+=1
            #self.sizerD.Add(self.screenLbl4,(cnt,0),flag=wx.EXPAND)
            #self.sizerD.Add(self.screenbox4   ,(cnt,1),flag=wx.EXPAND);cnt+=1
            self.sizerD.Add(self.cb_recip  ,(cnt,1),flag=wx.EXPAND);cnt+=1

            self.sizerD.Add(self.buttonScreen ,(cnt,1),flag=wx.EXPAND);cnt+=1
            self.Bind(wx.EVT_BUTTON, self.OnButtonScreen, self.buttonScreen)

            #FGA commented out- this was what was making the decon sizer so huge
            #self.decText=[]
            #for i in range(3):
            #    self.decText.append(wx.StaticText(self, label=""))
            #    self.sizerD.Add(self.decText[-1]   ,(cnt,0),flag=wx.EXPAND)
            #    #self.sizerD.Add(self.decText[-1]   ,flag=wx.EXPAND)
            #    cnt+=1


        if(self.DECONFULL==1):
            #############################
            #Sparky info
            self.sparkyLbl  = wx.StaticText(self, label="SparkyList:")
            self.sparkyfile     = wx.TextCtrl(self,   size=(200, -1))
            #FGA added
            self.selectSparkyFileBtn = wx.Button(self, label="...", size=(30,-1))
            self.selectSparkyFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.sparkyfile))
            #self.deconfile.SetValue(str(self.pars['distance_restraints']))
            self.buttonSparky = wx.Button(self, label="SparkyToMagma")

            cnt=0
            self.sizerS = wx.GridBagSizer(15, 3)
            self.sizerS.Add(self.sparkyLbl,(cnt,0),flag=wx.EXPAND)
            self.sizerS.Add(self.sparkyfile   ,(cnt,1),flag=wx.EXPAND)
            self.sizerS.Add(self.selectSparkyFileBtn, (cnt,2), flag=wx.EXPAND);cnt+=1
            self.sizerS.Add(self.buttonSparky  ,(cnt,1),flag=wx.EXPAND);cnt+=1
            self.Bind(wx.EVT_BUTTON, self.OnButtonSparky, self.buttonSparky)

        #self.decText=[]
        #for i in range(3):
        #    self.decText.append(wx.StaticText(self, label=""))
        #    self.sizerD.Add(self.decText[-1]   ,(cnt,0),flag=wx.EXPAND)
        #    #self.sizerD.Add(self.decText[-1]   ,flag=wx.EXPAND)
        #    cnt+=1




        #############################
        #NMR info
        self.noefileLbl  = wx.StaticText(self, label="NOE file:")
        self.noefile     = wx.TextCtrl(self,   size=(150, -1))
        self.noefile.SetValue(str(self.pars['distance_restraints']))
        #FGA changed
        self.chooseNoeFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.chooseNoeFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.noefile))


        self.mergefileLbl  = wx.StaticText(self, label="LVmergeFile:")
        self.mergefile     = wx.TextCtrl(self,   size=(200, -1))
        self.mergefile.SetValue(self.pars['LVmergeFile'])
        self.chooseLVFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.chooseLVFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.mergefile))


        self.ileshiftLbl  = wx.StaticText(self, label="Expt Shifts:")
        self.ileshift     = wx.TextCtrl(self,   size=(150, -1))
        self.selectNmrFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.selectNmrFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.ileshift))

        self.shiftXLbl = wx.StaticText(self, label="Calc Shifts:")
        self.shiftX = wx.TextCtrl(self, size=(150,-1))
        self.selectShiftXFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.selectShiftXFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.shiftX))
        #print self.pars['shiftXFile']
        self.shiftX.SetValue(self.pars['shiftXFile'])



        self.buttonShowN = wx.Button(self, label="ShowNOEs")
        #self.buttonInt = wx.Button(self, label="Intensities")

        self.Bind(wx.EVT_BUTTON, self.OnButtonShowN, self.buttonShowN)

        cnt=0
        self.sizerN.Add(self.noefileLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerN.Add(self.noefile   ,(cnt,1),flag=wx.EXPAND)
        self.sizerN.Add(self.chooseNoeFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerN.Add(self.mergefileLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerN.Add(self.mergefile   ,(cnt,1),flag=wx.EXPAND)
        self.sizerN.Add(self.chooseLVFileBtn ,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerN.Add(self.ileshiftLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerN.Add(self.ileshift,(cnt,1),flag=wx.EXPAND)
        self.sizerN.Add(self.selectNmrFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerN.Add(self.shiftXLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerN.Add(self.shiftX,(cnt,1),flag=wx.EXPAND)
        self.sizerN.Add(self.selectShiftXFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1
        self.sizerN.Add(self.buttonShowN   ,(cnt,1),flag=wx.EXPAND);cnt+=1
        #self.sizerN.Add(self.buttonInt   ,(cnt,1),flag=wx.EXPAND);cnt+=1


        #############################
        #MAGMA info

        self.protoLbl  = wx.StaticText(self, label="Protocol:")
        self.proto=[]
        self.proto.append('subgraphMode')
        self.proto.append('polishMode')
        self.proto.append('nudgeMode')
        self.proto.append('finalMode')
        self.proto.append('longMode')
        #self.proto.append('postDistanceMode')
        self.protoBox = wx.CheckListBox(self, -1, choices=self.proto,size=(150,107))

        check=[]
        for i,pro in enumerate(self.proto):
            if(len(str(self.pars[pro]).split("on"))>1):
                check.append(i)

        if(self.WXV==4):
            self.protoBox.SetCheckedItems(check)
        else:
            self.protoBox.SetChecked(check)

        self.processorsLbl  = wx.StaticText(self, label="Processors:")
        self.processors     = wx.TextCtrl(self,   size=(200, -1))
        self.processors.SetValue(str(self.pars['Nprocessors']))


        self.cb_inty = wx.CheckBox(self, -1,"UseIntensities",style=wx.ALIGN_RIGHT)
        if(self.pars['weightEdgesG1']=="on"):
            self.cb_inty.SetValue(True)

        self.cb_repro = wx.CheckBox(self, -1,"Reprioritise",style=wx.ALIGN_RIGHT)

        #try: #FGA thinkgs this is now fine!
        #    if(self.pars['Reprioritise']=="on"):
        #        self.cb_repro.SetValue(True)
        #except:
        #    pass




        self.buttonGraphs = wx.Button(self, label="GetGraphs")
        self.buttonMAGMA = wx.Button(self, label="RunMAGMA")
        self.buttonResults = wx.Button(self, label="Results")
        self.buttonPDF = wx.Button(self, label="MakePDF")

        self.buttonKill = wx.Button(self, label="KillMAGMA")


        #self.buttonSaveResults = wx.Button(self, label="Save Results")
        #self.buttonLoadResults = wx.Button(self, label="Load Results")
        #self.saveResultsAs = wx.TextCtrl(self, size=(100, -1))


        #self.source = wx.ListCtrl(self, -1, style = wx.LC_REPORT|wx.BORDER_SUNKEN,size=(250,100))
        #self.source.InsertColumn(0, 'protocol', width = 100)
        #self.source.InsertColumn(1, 'status', width = 100)
        #self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick, self.source)


        self.Bind(wx.EVT_BUTTON, self.OnButtonMAGMA, self.buttonMAGMA)
        self.Bind(wx.EVT_BUTTON, self.OnButtonGraphs, self.buttonGraphs)
        self.Bind(wx.EVT_BUTTON, self.OnButtonResults, self.buttonResults)
        self.Bind(wx.EVT_BUTTON, self.OnButtonPDF, self.buttonPDF)
        self.Bind(wx.EVT_BUTTON, self.OnButtonKill, self.buttonKill)
        #self.Bind(wx.EVT_BUTTON, self.OnSaveResults, self.buttonSaveResults)
        #self.Bind(wx.EVT_BUTTON, self.OnLoadResults, self.buttonLoadResults)

        self.approxLbl  = wx.StaticText(self, label="Approx:")
        approx=[]
        approx.append('stripMode')
        approx.append('communityMode')
        approx.append('logicFilter')
        approx.append('ileShiftMode')
        self.approxBox = wx.CheckListBox(self, -1, choices=approx,size=(150,85))
        check=[]
        for i,pro in enumerate(approx):
            if(len(str(self.pars[pro]).split("on"))>1):
                check.append(i)

        if(self.WXV==4):
            self.approxBox.SetCheckedItems(check)
        else:
            self.approxBox.SetChecked(check)



        #self.pars['reprioritise']=Parse(deconParFile,'reprioritise')


        self.preoptLbl  = wx.StaticText(self, label="PreOpt:")
        self.preopt     = wx.TextCtrl(self,   size=(150, -1))
        #print self.pars['priority_iterations'],self.pars['maximum_run_time']
        self.preopt.SetValue(str(self.pars['priority_iterations'])+' '+str(self.pars['maximum_run_time']))

        self.craicLbl  = wx.StaticText(self, label="FixAssigns:")
        self.craic     = wx.TextCtrl(self,   size=(150, -1))
        self.selectCraicFileBtn = wx.Button(self, label="...", size=(30,-1))
        self.selectCraicFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.craic))


        #self.ileshiftLbl  = wx.StaticText(self, label="ileShift:")
        #self.ileshift     = wx.TextCtrl(self,   size=(150, -1))
        #if(self.pars['ileShiftMode']=="on"):
        #    self.ileshift.SetValue(self.pars['ileShiftFile'])


        cnt=0
        self.sizerM.Add(self.protoLbl   ,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.protoBox   ,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerM.Add(self.approxLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.approxBox,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerM.Add(self.cb_inty,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerM.Add(self.cb_repro,(cnt,1),flag=wx.EXPAND);cnt+=1

        self.sizerM.Add(self.preoptLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.preopt,(cnt,1),flag=wx.EXPAND);cnt+=1

        self.sizerM.Add(self.craicLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.craic,(cnt,1),flag=wx.EXPAND);
        self.sizerM.Add(self.selectCraicFileBtn,(cnt,2),flag=wx.EXPAND);cnt+=1

        #self.sizerM.Add(self.ileshiftLbl,(cnt,0),flag=wx.EXPAND)
        #self.sizerM.Add(self.ileshift,(cnt,1),flag=wx.EXPAND);cnt+=1

        self.sizerM.Add(self.processorsLbl,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.processors,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerM.Add(self.buttonGraphs,(cnt,0),flag=wx.EXPAND)
        self.sizerM.Add(self.buttonMAGMA,(cnt,1),flag=wx.EXPAND);cnt+=1
        self.sizerM.Add(self.buttonKill,(cnt,0),flag=wx.EXPAND);

        self.butSize= wx.BoxSizer(wx.HORIZONTAL)
        self.butSize.Add(self.buttonResults)
        self.butSize.Add(self.buttonPDF)
        self.sizerM.Add(self.butSize,(cnt,1),flag=wx.EXPAND);cnt+=1

        #FGA added
        #self.sizerM.Add(self.saveResultsAs,(cnt,0),flag=wx.EXPAND)
        #self.butSize2 = wx.BoxSizer(wx.HORIZONTAL)
        #self.butSize2.Add(self.buttonSaveResults)
        #self.butSize2.Add(self.buttonLoadResults)
        #self.sizerM.Add(self.butSize2,(cnt,1),flag=wx.EXPAND);cnt+=1


        #############################
        #report info

        cnt=0
        self.repText=[]
        for i in range(30):
            self.repText.append(wx.StaticText(self, label=""))
            #self.sizerR.Add(self.repText[-1]   ,(cnt,0),flag=wx.EXPAND)
            self.sizerR.Add(self.repText[-1]   ,flag=wx.EXPAND)
            cnt+=1

        #############################

        self.pdbLbl = wx.StaticBox(self, -1, 'PDB:', size=(250, 140))
        self.pdbSizer = wx.StaticBoxSizer(self.pdbLbl, wx.VERTICAL)
        self.borderP = wx.BoxSizer()
        self.borderP.Add(self.sizerP, 1, wx.ALL | wx.EXPAND, 7)
        self.pdbSizer.Add(self.borderP)

        self.nmrLbl = wx.StaticBox(self, -1, 'NMR:', size=(240, 140))
        self.nmrSizer = wx.StaticBoxSizer(self.nmrLbl, wx.VERTICAL)
        self.borderN = wx.BoxSizer()
        self.borderN.Add(self.sizerN, 1, wx.ALL | wx.EXPAND, 7)
        self.nmrSizer.Add(self.borderN)

        if(self.DECONFULL==1):
            self.decLbl = wx.StaticBox(self, -1, 'Decon:', size=(240, 140))
            self.decSizer = wx.StaticBoxSizer(self.decLbl, wx.VERTICAL)
            self.borderD = wx.BoxSizer()
            self.borderD.Add(self.sizerD, 1, wx.ALL | wx.EXPAND, 7)
            self.decSizer.Add(self.borderD)
            self.sparkLbl = wx.StaticBox(self, -1, 'Sparky:', size=(240, 140))
            self.sparkSizer = wx.StaticBoxSizer(self.sparkLbl, wx.VERTICAL)
            self.borderS = wx.BoxSizer()
            self.borderS.Add(self.sizerS, 1, wx.ALL | wx.EXPAND, 7)
            self.sparkSizer.Add(self.borderS)




        self.magmaLbl = wx.StaticBox(self, -1, 'MAGMA:', size=(240, 140))
        self.magmaSizer = wx.StaticBoxSizer(self.magmaLbl, wx.VERTICAL)
        self.borderM = wx.BoxSizer(wx.VERTICAL)
        self.borderM.Add(self.sizerM, 1, wx.ALL | wx.EXPAND, 7)
        #self.tmpM = wx.BoxSizer(wx.VERTICAL)
        #self.tmpM.Add(self.source,21,wx.ALL|wx.EXPAND,27)
        #self.borderM.Add(self.tmpM)
        #self.SetSizer(self.borderM)


        self.magmaSizer.Add(self.borderM)

        self.repLbl = wx.StaticBox(self, -1, 'Report:', size=(240, 140))
        self.repSizer = wx.StaticBoxSizer(self.repLbl, wx.VERTICAL)
        self.borderR = wx.BoxSizer()
        self.borderR.Add(self.sizerR, 1, wx.ALL | wx.EXPAND, 7)
        self.repSizer.Add(self.borderR)




        self.splitSizer.Add(self.pdbSizer)

        self.decn = wx.BoxSizer(wx.VERTICAL)

        if(self.DECONFULL==1):
            self.decn.Add(self.decSizer)
            self.decn.Add(self.sparkSizer)
        self.decn.Add(self.nmrSizer)

        self.splitSizer.Add(self.decn)
        self.splitSizer.Add(self.magmaSizer)
        self.splitSizer.Add(self.repSizer)

        self.SetSizerAndFit(self.splitSizer)



        #FGA changed- make decon, sparky and nmr sizers the same width
        nmrWidth = self.nmrSizer.GetSize()[0]

        if(self.DECONFULL):
            decWidth = self.decSizer.GetSize()[0]
            sparkWidth = self.sparkSizer.GetSize()[0]
            widthToSetDecn = max(nmrWidth, decWidth, sparkWidth)
        else:
            widthToSetDecn = max(nmrWidth,0)

        if(self.DECONFULL):
            self.decSizer.SetMinSize((widthToSetDecn, 0))
            self.sparkSizer.SetMinSize((widthToSetDecn, 0))
        self.nmrSizer.SetMinSize((widthToSetDecn, 0))
        self.SetSizerAndFit(self.splitSizer)

        #FGA added  CHECLKCHECK
        self.UpdatePars()

    #go through a list, compare to self.par['key']
    def SyncBox(self,box,listy,key):
        test=str(self.pars[key])
        #print test
        check=[]
        for i,li in enumerate(listy):
            if(len(test.split(li))>1):

                #print 'found',li
                check.append(i)
            #else:
            #    box.SetCheckedItems(i)=False

        if(self.WXV==4):
            box.SetCheckedItems(check)
        else:
            box.SetChecked(check)

    def IntToOnOff(self,test):
        if(test):
            return "on"
        else:
            return "off"


    def GetCombo(self,comb,listy):
        mask=numpy.array(comb.GetChecked())
        str=''
        for li in listy[mask]:
            str+=li+' '
        return str

    def WriteMAGMA(self):
        self.SyncPars()

        dec=[]
        if(os.path.exists(self.deconParFile)):
            inny=open(self.deconParFile)
            for line in inny.readlines():
               dec.append(line)
            inny.close()

        #write non-magma pars
        # print(de)
        outy=open(self.deconParFile,'w')
        for de in dec:
            test=de.split()
            if(len(test)>0):
                if(test[0] in list(self.pars.keys())):
                    outy.write('%s %s ' % (test[0],self.pars[test[0]]))
                    #for j in range(len(self.pars[test[0]])-2):
                    #   outy.write('%s' % self.pars[test[0]][j+2])
                    outy.write('\n')
                    del self.pars[test[0]]
                else:
                    outy.write(de)
            else:
                outy.write(de)


        # outy=open(self.deconParFile,'w')
        for key,vals in list(self.pars.items()):
            if(vals!='notSpecified'):
                outy.write('%s %s\n' % (key,vals))
        outy.close()

        self.SyncPars()


    def SyncPars(self):


        self.pars['pdb_file']=self.pdbfile.GetValue()
        self.pars['chains']=self.chain.GetValue()

        try:
            self.pars['residues']=self.GetCombo(self.methylBox,self.methyls)
        except:
            self.pars['residues']=''

        self.pars['mergeLV']=self.IntToOnOff(self.cb_LV.IsChecked())
        self.pars['short_distance_threshold']=self.thresh.GetValue()
        self.pars['long_distance_threshold']=self.longthresh.GetValue()

        self.pars['shortLim']=self.shortthresh.GetValue()

        #self.pars['proRS']=self.GetCombo(self.LVBox,self.LVs)
        self.pars['proRS']=wx.RadioBox.GetStringSelection(self.LVBox)

        self.pars['distance_restraints']=self.noefile.GetValue()

        self.pars['LVmergeFile']=self.mergefile.GetValue()
        #self.mergefile     = wx.TextCtrl(self,   size=(200, -1))

        self.pars['subgraphMode']=self.IntToOnOff(self.protoBox.IsChecked(0))
        self.pars['polishMode']=self.IntToOnOff(self.protoBox.IsChecked(1))
        self.pars['nudgeMode']=self.IntToOnOff(self.protoBox.IsChecked(2))
        self.pars['finalMode']=self.IntToOnOff(self.protoBox.IsChecked(3))
        self.pars['longMode']=self.IntToOnOff(self.protoBox.IsChecked(4))

        self.pars['weightEdgesG1']=self.IntToOnOff(self.cb_inty.IsChecked())
        self.pars['reprioritise']=self.IntToOnOff(self.cb_repro.IsChecked())

        self.pars['Nprocessors']=self.processors.GetValue()
        self.pars['priority_iterations']=self.preopt.GetValue().split()[0]
        self.pars['maximum_run_time']=self.preopt.GetValue().split()[1]
        #self.pars['reprioritise']=self.preopt.GetValue().split()[1]
        self.pars['stripMode']=self.IntToOnOff(self.approxBox.IsChecked(0))
        self.pars['communityMode']=self.IntToOnOff(self.approxBox.IsChecked(1))
        self.pars['logicFilter']=self.IntToOnOff(self.approxBox.IsChecked(2))
        self.pars['ileShiftMode']=self.IntToOnOff(self.approxBox.IsChecked(3))

        if(len(self.ileshift.GetValue())>0 and self.ileshift.GetValue()!='notSpecified'):
            self.pars['ileShiftMode']="on"
            self.pars['ileShiftFile']=self.ileshift.GetValue()
        else:
            self.pars['ileShiftMode']="off"
            try:
                del self.pars['ileShiftFile']
            except:
                pass

        if(len(self.craic.GetValue())>0 and self.craic.GetValue()!='notSpecified'):
            self.pars['craicMode']="on"
            self.pars['craicFile']=self.craic.GetValue()
        else:
            self.pars['craicMode']="off"
            try:
                del self.pars['craicFile']
            except:
                pass

        if(len(self.shiftX.GetValue())>0 and self.shiftX.GetValue()!='notSpecified'):
            self.pars['shiftXFile']=self.shiftX.GetValue()
        else:
            try:
                del self.pars['shiftXFile']
            except:
                pass



    #FGA added
    def UpdatePars(self):
        self.GetPars()

        #PDB
        self.pdbfile.SetValue(self.pars['pdb_file'])
        self.chain.SetValue(self.pars['chains'])
        self.SyncBox(self.methylBox,self.methyls,'residues')

        if(self.pars['mergeLV']=="on"):
            self.cb_LV.SetValue(True)
        else:
            self.cb_LV.SetValue(False)
        self.thresh.SetValue(str(self.pars['short_distance_threshold']))
        self.longthresh.SetValue(str(self.pars['long_distance_threshold']))

        choice=0
        for i in range(len(self.LVs)):
            if(self.LVs[i]==self.pars['proRS']):
                choice=i
        #choice = self.LVs.GetCheckedItems(self.pars['proRS'])
        self.LVBox.SetSelection(choice)

        #Decon
        #self.deconfile.SetValue('out/correlate.2')
        #self.screenbox.SetValue('4')

        #NMR
        self.noefile.SetValue(str(self.pars['distance_restraints']))
        self.mergefile.SetValue(self.pars['LVmergeFile'])


        #MAGMA
        self.proto=[]
        self.proto.append('subgraphMode')
        self.proto.append('polishMode')
        self.proto.append('nudgeMode')
        self.proto.append('finalMode')
        self.proto.append('longMode')
        #self.proto.append('postDistanceMode')
        self.checkBox(self.proto, self.protoBox, 'on')

        self.processors.SetValue(str(self.pars['Nprocessors']))

        if(self.pars['weightEdgesG1']=="on"):
            self.cb_inty.SetValue(True)
        else:
            self.cb_inty.SetValue(False)

        if(self.pars['reprioritise']=="on"):
            self.cb_repro.SetValue(True)
        else:
            self.cb_repro.SetValue(False)
        approx=[]
        approx.append('stripMode')
        approx.append('communityMode')
        approx.append('logicFilter')
        approx.append('ileShiftMode')
        self.checkBox(approx, self.approxBox, 'on')



        if self.pars['priority_iterations'] and self.pars['maximum_run_time']:
            self.preopt.SetValue(str(self.pars['priority_iterations'])+' '+str(self.pars['maximum_run_time']))
        else:
            self.preopt.SetValue('')

        if self.pars['craicMode'] == 'on':
            self.craic.SetValue(self.pars['craicFile'])
        else:
            self.craic.SetValue('')

        if 'ileShiftFile' in self.pars:
            self.ileshift.SetValue(self.pars['ileShiftFile'])
        else:
            self.ileshift.SetValue('')

        if 'shiftXFile' in self.pars:
            self.shiftX.SetValue(self.pars['shiftXFile'])
        else:
            self.shiftX.SetValue('')

    def checkBox(self, list, box, split):
        check=[]
        for i,pro in enumerate(list):
            if(len(str(self.pars[pro]).split(split))>1):
                check.append(i)
        #box.SetCheckedItems(check)
        if(self.WXV==4):
            box.SetCheckedItems(check)
        else:
            box.SetChecked(check)


    #FGA added
    def OnRadioBox(self, e):
        self.pars['proRS'] = self.LVBox.GetStringSelection()

    #FGA added
    def onGetFile(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" \
            "All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
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
            os.chdir(path)
            self.dirVal.SetValue(path)
            print("CWD: ",os.getcwd())

        dlg.Destroy()

    #FGA added
    """
    def OnSaveResults(self, event):
        if os.path.exists('results') == 1:
            if self.saveResultsAs.GetValue():
                newPath = self.saveResultsAs.GetValue()
                os.system('cp -r results/ results'+newPath)
            else:
                msg = "Could not save results as you have not entered a results filename."
                dlg = wx.MessageDialog(self, msg, "Oops", wx.OK)
                dlg.ShowModal()
                dlg.Destroy()
            self.saveResultsAs.SetValue('')
    """

    #FGA added
    def OnLoadResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Load session...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.deconParFile=path
            self.GetPars()
            self.UpdatePars()

            #self.canvas.print_figure(path, dpi=self.dpi)
            self.parent.parent.flash_status_message("Loaded %s" % path)
            os.chdir(os.path.dirname(path))
            print("CWD: ",os.getcwd())
            self.dirVal.SetValue(os.path.dirname(path))


    def OnSaveResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            defaultFile=os.path.split(self.deconParFile)[1],
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            self.deconParFile=path
            if(os.path.exists(self.deconParFile)==0):
                outy=open(self.deconParFile,'w');outy.close()
            self.OnButtonSave(True)



            #self.canvas.print_figure(path, dpi=self.dpi)
            self.parent.parent.flash_status_message("Saved %s" % path)



    def OnButtonKill(self,event):
        print('Tring to Kill Calculation...')


        #print self.magmaCalc.poll()
        if(self.magmaCalc!=0):
            #os.killpg(os.getpgid(self.magmaCalc.pid),signal.SIGTERM)
            self.magmaCalc.terminate()
            #poll = self.magmaCalc.poll()
            #print poll
            #if poll == None:
            #    self.magmaCalc.kill()
            #    os.system('killall -9 mcesCore_parallel_Darwin')
            #    os.system('killall -9 mcesCore_Darwin')
            #    print 'Success.'
            time.sleep(2)
            print('sending killall signal to binaries.')
            os.system('killall -9 mcesCore_parallel_Darwin')
            os.system('killall -9 mcesCore_Darwin')
            os.system('killall -9 mcesCore_parallel_Linux')
            os.system('killall -9 mcesCore_Linux')
            self.magmaCalc=0
        else:
            print('No calculation running.')
            #except:
            #print 'Failed to kill calculation.'
            #pass

        time.sleep(2)
        print('sending killall signal to binaries.')
        os.system('killall -9 mcesCore_parallel_Darwin')
        os.system('killall -9 mcesCore_Darwin')
        os.system('killall -9 mcesCore_parallel_Linux')
        os.system('killall -9 mcesCore_Linux')


    def OnButtonPDF(self,event):


        self.inst=Magma(self.deconParFile,run='n') #get instance of magma
        self.rep=Report(self.inst) #setup report class



        #if(inst.P.subgraphMode=="on" and inst.P.finalMode=="on"):
        #    #rep.result_dict=rep.ReadRes(infile=inst.P.outdir+'_full/combinedResults.res')
        #    rep.result_dict=rep.ReadRes(infile=inst.P.outdir+'/combinedResults.res')
        #else:
        #    rep.ReadRes()

        if(os.path.exists('report/test1.png')==0 or os.path.exists('report/test2.png')==0 or os.path.exists('report/test3.png')==0):
            test=Visualise(self.inst) #make pngs using pymol if not already there

            test.pymol_gen(run='y')
            test.pymol_gen()

            test.monte_min()
            test.pymol_gen_noe(run='y') #make pretty output
            test.pymol_gen_noe() #make pretty output

        if(os.path.exists('report/test1.png')==0 or os.path.exists('report/test2.png')==0 or os.path.exists('report/test3.png')==0):
            print()
            print('Undetermined problem in making the png images using pymol. Trying to carry on.')
            print()

        try:
            self.rep.MakeTexReport()
        except:
            pass

        self.UpdateReport()

        #FGA added
        from spinDecon.gui.dialogs.pdf_viewer import PDFViewer
        pdfV = PDFViewer(None, size=(800, 600))
        pdfV.viewer.UsePrintDirect = False
        if os.path.exists('summary.pdf') == 1:
            pdfV.viewer.LoadFile('summary.pdf')
        pdfV.Show()


        """
        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        #set options in listbox
        try:
            for i in range(len(self.rep.progress)):
                num_items = self.source.GetItemCount()
                self.source.InsertStringItem(num_items,self.rep.progress[i][0])
                self.source.SetStringItem(num_items,0,self.rep.progress[i][0])

                if(os.path.exists(self.rep.progress[i][1]+'/combinedResults.res')): #is this mode done?
                    self.source.SetStringItem(num_items,1,"complete")
                else:
                    self.source.SetStringItem(num_items,1,"incomplete")
        except:
            num_items = self.source.GetItemCount()
            self.source.InsertStringItem(num_items,'subgraphMode')
            self.source.SetStringItem(num_items,0,'subgraphMode')
            self.source.SetStringItem(num_items,1,'incomplete')
        """

    def OnButtonResults(self,event):
        if(os.path.exists(self.deconParFile)==0):
           print('Cannot find magma input file')
           return
        #self.colVal=self.source.GetFirstSelected()
        #count = self.source.GetItemCount()
        #self.col1 = [self.source.GetItem(itemId=row, col=0).GetText() for row in xrange(count)][self.colVal]
        #print self.col1,self.colVal

        self.inst=Magma(self.deconParFile,run='n')
        self.UpdateReport()
        from spinDecon.integrations.magma import results as magmaResults
        magmaResults=imp.reload(magmaResults)
        bool=magmaResults.MagmaResultsMan(self)


    def OnButtonShowN(self,event):
        self.WriteMAGMA()

        self.inst=Magma(self.deconParFile,run='n')
        test=Visualise(self.inst)
        test.monte_min()
        test.pymol_gen_noe()
        self.UpdateReport()

        self.calcy=os.popen('pymol report/pyscript_noe.py')
        # self.PymolExec('pyscript_noe.py')


    def OnButtonShowP(self,event):
        self.WriteMAGMA()
        self.OnButtonSave(True)
        if(os.path.exists(self.pdbfile.GetValue())==0):
            print('cannot find pdbfile')
            return
            #self.pdbfile     = wx.TextCtrl(self,   size=(200, -1))
        self.inst=Magma(self.deconParFile,run='n')
        if(self.inst.success==False):
            print('aborting')
            return
        test=Visualise(self.inst)
        test.pymol_gen()
        self.UpdateReport()

        self.PymolExec('pyscript.py')



    #annoying: need to get X11 permissions to run pymol when in package.
    #doing this via 'open terminal'
    def PymolExec(self,pyscript):
        #print os.getcwd()


        # if(os.uname()[0]=='Linux'):
            self.calcy=os.popen('pymol report/%s' % (pyscript,))
        # else:
        #     files='./report/tmp.sh'
        #     outy=open(files,'w')
        #     outy.write('cd %s\n' % (os.getcwd()))
        #     outy.write('pwd\n')
        #     outy.write('chmod +x %s/report/tmp.sh\n' % os.getcwd())
        #     outy.write('`which pymol` %s/report/%s\n' % (os.getcwd(),pyscript))
        #     outy.close()
        #     self.calcy=os.popen('open -a Terminal %s/report/tmp.sh' % os.getcwd())






    def UpdateReport(self):
        small_graph = self.inst.G.NetworkxGraph(self.inst.noe_node_list,self.inst.noe_adjacency)
        big_graph = self.inst.G.NetworkxGraph(self.inst.short_node_list,self.inst.short_adjacency)
        for i in range(len(self.repText)):
            self.repText[i].SetLabel("")
        cnt=0
        self.repText[cnt].SetLabel("PDB:");cnt+=1
        self.repText[cnt].SetLabel("Nodes: %i Edges: %i" % (len(big_graph.nodes()),len(big_graph.edges())));cnt+=1
        stry=''
        for res in self.inst.P.residues: #print counts of residues from log file
            ref=self.inst.PDB.dicty[res]['ref']
            stry+='%s:%s '% (res,self.inst.PDB.readin[res]/len(self.inst.P.chains))
        self.repText[cnt].SetLabel(stry);cnt+=1
        self.repText[cnt].SetLabel("NMR:");cnt+=1
        self.repText[cnt].SetLabel("Restraints: %s" % (self.inst.CountNOEs(self.inst.noe_adjacency)));cnt+=1
        stry=''
        for res in self.inst.P.residues: #print counts of residues from log file

            try:
                ref=str(len(self.inst.NMR.residues[res[0]]))
                stry+='%s:%s '% (res,ref)
            except:
                pass
        self.repText[cnt].SetLabel(stry);cnt+=1

        self.repText[cnt].SetLabel("Sparsity: %.2f" % (len(small_graph.edges())/float(len(big_graph.edges()))));cnt+=1

        self.repText[cnt].SetLabel('Nodes: %i Edges: %i\n\n' % (len(small_graph.nodes()),len(small_graph.edges())));cnt+=1;
        stry=''
        for res in self.inst.P.residues: #print counts of residues from log file
            ref=self.inst.PDB.dicty[res]['ref']
            stry+='%s:%s '% (res,self.inst.PDB.readin[res]/len(self.inst.P.chains))
        stry=''
        for i in range(len(list(self.inst.subgraphRef.keys()))):
            self.repText[cnt].SetLabel('Subgraph  '+str(i+1)+' Nodes: '+str(len(self.inst.subgraphRef[i]['nodes']))+' Edges: '+str(self.inst.subgraphRef[i]['noes']));cnt+=1

        """
        try:
            #poll = self.magmaCalc.poll()
            #if poll == None:
            #    self.repText[cnt].SetLabel('Calculation running...');cnt+=1
            tack=self.magmaCalc.is_alive()
            if tack == True:
                self.repText[cnt].SetLabel('Calculation running...');cnt+=1

        except:
            #self.repText[cnt].SetLabel('Calculation still running');cnt+=1
            pass
        """

        if(os.path.exists('results')==0):
            os.system('mkdir results')
        if(os.path.exists(self.inst.P.outdir)==0):
            os.system('mkdir '+self.inst.P.outdir)
        if(os.path.exists(self.deconParFile)==0):
            self.WriteMAGMA()
        if(os.path.exists(self.inst.P.outdir+'/input.magma')==0):
            print(self.deconParFile, os.path.normpath(self.deconParFile))
            import shutil
            shutil.copy2(self.deconParFile, self.inst.P.outdir)

        self.GetModes()

        self.repText[cnt].SetLabel('PROGRESS:');cnt+=1
        for i in range(len(self.progress)):
            self.repText[cnt].SetLabel('%s %s' % (self.progress[i][0],self.progress[i][3]));cnt+=1


    def ParseMagma(self):
        inny=open(self.inst.P.outdir+'/input.magma')
        check=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(test[0]=='subgraphMode' and test[1]=="on"):
                    check.append(0)
                if(test[0]=='polishMode' and test[1]=="on"):
                    check.append(1)
                if(test[0]=='nudgeMode' and test[1]=="on"):
                    check.append(2)
                if(test[0]=='finalMode' and test[1]=="on"):
                    check.append(3)
                if(test[0]=='longMode' and test[1]=="on"):
                    check.append(4)
        if(self.WXV==4):
            self.protoBox.SetCheckedItems(check)
        else:
            self.protoBox.SetChecked(check)
    def GetModes(self):
        #self.ParseMagma()
        self.progress=[]
        cnt=0

        if(self.protoBox.IsChecked(0)):
            cnt+=1
            self.progress.append(('subgraphMode',self.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.protoBox.IsChecked(1)):
            cnt+=1
            self.progress.append(('polishMode',self.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.protoBox.IsChecked(2)):
            cnt+=1
            self.progress.append(('nudgeMode',self.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.protoBox.IsChecked(3) or self.protoBox.IsChecked(0)==False):
            cnt+=1
            self.progress.append(('finalMode',self.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.protoBox.IsChecked(4)):
            cnt+=1
            self.progress.append(('longMode',self.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))


        self.progress=numpy.array(self.progress)
        self.progress[-1][1]=self.inst.P.outdir #adjust last folder

        cnt=0
        for i in range(len(self.progress)-1):
            if(os.path.exists(self.progress[i][1])):
                self.progress[i,3]='complete'
                cnt+=1
        if(cnt==len(self.progress)-1):
            if(os.path.exists(self.progress[-1,1]+'/combinedResults.res')):
                self.progress[-1,3]='complete'


    def OnButtonGraphs(self,event):
        self.WriteMAGMA()
        self.inst=Magma(self.deconParFile,run='n')
        self.UpdateReport()

    def OnButtonMAGMA(self,event):

        try:
            tack=self.magmaCalc.is_alive()
            if(tack==True):
                print('Calculation already running.')
                print('Kill it if you are done.')
                return
        except:
            pass

        self.WriteMAGMA()

        self.inst=Magma(self.deconParFile,run='n')
        self.rep=Report(self.inst) #setup report class
        self.UpdateReport()
        #run magam in subprocess, keep track of system pid for easy shutdown
        print('running: magmaRun '+self.deconParFile)
        self.magmaCalc=TestThread(self.deconParFile)
        import multiprocessing
        self.currentProc=multiprocessing.current_process()

        #print self.currentProc.pid
        #print self.currentProc.name
        #self.magmaCalc=subprocess.Popen('magmaRun '+self.deconParFile,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True,preexec_fn=os.setsid)


    def OnButtonSparky(self,event):
        pass

    #take correlate peak name, and make MAGMA peak name
    def FormatPeak(self,raw):
        #raw=raw.split('C-H')[0]
        test1= re.findall(r'[0-9]+',raw)[0]
        return raw[0]+str(test1)



    def MakeMagma(self,peaks,connEntry,exclude,outy):
            s1=connEntry.s1
            s2=connEntry.s2

            sval=s1 #value to save
            #sval=s1/(connEntry.hm+connEntry.cm) #value to save
            #sval=connEntry.Intscore
            #sval=connEntry.distScore
            #sval=connEntry.n2
            #print self.parent.tabOne.dmax*self.sn,s1,s2


            if(numpy.fabs(s1)<self.parent.tabOne.dmax*self.sn):
                print('by',s1)
                return
            #if(numpy.fabs(s2)<self.parent.tabOne.dmax*self.sn):
            #    print 'bye',s2
            #    return


            """
            self.rat=0.5
            norat=0  #if both are intense...
            if(s1>5*self.parent.tabOne.dmax*self.sn):
                if(s2>5*self.parent.tabOne.dmax*self.sn):
                    norat=1 #skip the ratio screen

            if(self.rat!=0 and norat==0):
                rat1=s1/s2
                if(rat1<self.rat):
                    return
                rat2=s2/s1
                if(rat2<self.rat):
                    return
            """

            p1=connEntry.p1
            p2=connEntry.p2
            c1=self.FormatPeak(p1)
            c2=self.FormatPeak(p2)

            #if(numpy.fabs(s1)>snLim and numpy.fabs(s2)>snLim):

            if(c1==c2):
                return 0
            if(p1 in exclude):
                return 0
            if(p2 in exclude):
                return 0

            if(self.cb_recip.IsChecked()):
                recip=0
                for j,connEntry2 in enumerate(self.parent.tabOne.conn_data): #enforce reciprocity
                    d2=connEntry2.p1 #self.FormatPeak(connEntry2.p1)
                    d1=connEntry2.p2 #self.FormatPeak(connEntry2.p2)
                    #if(1==1):
                    if(d1==p1 and d2==p2):
                        recip=1
                        break

                if(recip==0):
                    return 0


            outy.write('%s\t%s\t%e\n' % (c1,c2,numpy.fabs(sval)))
            if(c1 not in peaks):
                peaks.append(c1)
            if(c2 not in peaks):
                peaks.append(c2)
            return 1



    def OnButtonScreen(self,event):

        self.sn=float(self.screenbox.GetValue())



        #read in deconfile
        infile=self.deconfile.GetValue() #get decon output file
        # if(self.parent.tabOne.GetConn(infile)==-1): #read it in
        #     return

        stry="Raw Crosspeaks: %i" % (len(self.parent.tabOne.conn_data))
        #stry+=" Diagonals:  %i" % (len(self.parent.tabOne.diag))

        #snLim=float(self.screenbox.GetValue())

        #cnt=0
        #self.decText[cnt].SetLabel(stry);cnt+=1

        outfile=self.noefile.GetValue()    #write noe file in magma format
        outfile1=self.mergefile.GetValue() #write mergefile
        print('writing:',outfile)

        print(self.screenbox2.GetValue())
        print(self.screenbox3.GetValue())
        try:
            hlim=float(self.screenbox2.GetValue())
            clim=float(self.screenbox3.GetValue())
        except:
            hlim=0
            clim=0

        #try:
        #    self.rat=float(self.screenbox4.GetValue())
        #except:
        #    self.rat=0

        print(hlim,clim)

        exclude=[]

        #self.parent=parent
        dim=self.parent.tabOne.dim
        if(dim==4):
            if(hlim!=0 and clim!=0):
                for pk in self.parent.tabOne.peak:
                    for pk2 in self.parent.tabOne.peak:
                        hdiff=numpy.fabs(pk.ppmK-pk2.ppmK)
                        cdiff=numpy.fabs(pk.ppmL-pk2.ppmL)
                        #print pk.name,pk2.name,hdiff,cdiff
                        if(pk.name!=pk2.name):
                            if(hdiff<hlim and cdiff<clim):
                                addy=pk.name
                                addy2=pk2.name
                                print('Peaks too close - excluding:',addy,addy2)
                                if(addy not in exclude):
                                    exclude.append(addy)
                                if(addy2 not in exclude):
                                    exclude.append(addy2)
        else:
            if(hlim!=0 and clim!=0):
                for pk in self.parent.tabOne.peak:
                    for pk2 in self.parent.tabOne.peak:
                        hdiff=numpy.fabs(pk.ppmK-pk2.ppmK)
                        cdiff=numpy.fabs(pk.ppmJ-pk2.ppmJ)
                        #print pk.name,pk2.name,hdiff,cdiff,pk.ppmI,pk.ppmJ,pk.ppmK
                        if(pk.name!=pk2.name):
                            if(hdiff<hlim and cdiff<clim):
                                addy=pk.name
                                addy2=pk2.name
                                print('Peaks too close - excluding:',addy,addy2)
                                if(addy not in exclude):
                                    exclude.append(addy)
                                if(addy2 not in exclude):
                                    exclude.append(addy2)

        #exclude.append('A114C-H')
        #exclude.append('A91C-H')
        #exclude.append('V107G2C-H')
        #exclude.append('L17D2C-H')
        #exclude.append('V176G2C-H')
        #exclude.append('L123D2C-H')
        print(exclude)


        outy=open(outfile,'w')
        peaks=[]
        for i,connEntry in enumerate(self.parent.tabOne.conn_data):
            self.MakeMagma(peaks,connEntry,exclude,outy)
        outy.close()

        done=[] #write diag
        outy1=open(outfile1,'w')
        for pe in peaks:
            for po in peaks:
                if(pe!=po):
                    test1= re.findall(r'[0-9]+',pe)[0]
                    test2= re.findall(r'[0-9]+',po)[0]
                    if(test1==test2):
                        if(test1 not in done):
                            outy1.write('%s\t%s\n' % (pe,po))
                            done.append(test1)
        outy1.close()





        #screen for whatever parameter

        #write a magma input file
        pass
    def OnButtonSave(self,event):
        print('Saving to ',self.deconParFile)
        self.SyncPars()
        write={}
        for key,vals in list(self.pars.items()):
            write[key]=vals

        update_parameter_file(self.deconParFile, write, source_path=self.deconParFile)

    def Parse(self,infile,key):
        if(os.path.exists(infile)==0):
            return 'notSpecified'
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>1):
                if(test[0]==key):
                    str=''
                    for i in range(len(test)-1):
                        str+=test[i+1]
                    return str
        return 'notSpecified'


    def GetPars(self):
        self.pars={}
        self.pars['pdb_file']=self.Parse(self.deconParFile,'pdb_file').split()[0]
        self.pars['chains']=self.Parse(self.deconParFile,'chains')

        chainTmp=str(self.pars['chains'])

        chainVals=''
        for i in range(len(chainTmp)):
            chainVals+=chainTmp[i]+' '
        self.pars['chains']=chainVals


        self.pars['residues']=self.Parse(self.deconParFile,'residues')
        self.pars['short_distance_threshold']=self.Parse(self.deconParFile,'short_distance_threshold').split()[0]
        self.pars['long_distance_threshold']=self.Parse(self.deconParFile,'long_distance_threshold').split()[0]
        self.pars['shortLim']=self.Parse(self.deconParFile,'shortLim').split()[0]

        self.pars['distance_restraints']=self.Parse(self.deconParFile,'distance_restraints').split()[0]
        self.pars['LVmergeFile']=self.Parse(self.deconParFile,'LVmergeFile').split()[0]
        self.pars['mergeLV']=self.Parse(self.deconParFile,'mergeLV').split()[0]
        self.pars['proRS']=self.Parse(self.deconParFile,'proRS').split()[0]
        self.pars['subgraphMode']=self.Parse(self.deconParFile,'subgraphMode').split()[0]
        self.pars['polishMode']=self.Parse(self.deconParFile,'polishMode').split()[0]
        self.pars['nudgeMode']=self.Parse(self.deconParFile,'nudgeMode').split()[0]
        self.pars['finalMode']=self.Parse(self.deconParFile,'finalMode').split()[0]
        self.pars['longMode']=self.Parse(self.deconParFile,'longMode').split()[0]
        self.pars['Nprocessors']=self.Parse(self.deconParFile,'Nprocessors').split()[0]



        if(self.pars['short_distance_threshold']=='notSpecified'):
            self.pars['short_distance_threshold']=10
        if(self.pars['long_distance_threshold']=='notSpecified'):
            self.pars['long_distance_threshold']=15
        if(self.pars['shortLim']=='notSpecified'):
            self.pars['shortLim']=6
        if(self.pars['Nprocessors']=='notSpecified'):
            self.pars['Nprocessors']=1


        self.pars['weightEdgesG1']=self.Parse(self.deconParFile,'weightEdgesG1').split()[0]
        self.pars['reprioritise']=self.Parse(self.deconParFile,'reprioritise').split()[0]



        self.pars['priority_iterations']=self.Parse(self.deconParFile,'priority_iterations').split()[0]

        self.pars['maximum_run_time']=self.Parse(self.deconParFile,'maximum_run_time').split()[0]
        self.pars['reprioritise']=self.Parse(self.deconParFile,'reprioritise').split()[0]
        self.pars['stripMode']=self.Parse(self.deconParFile,'stripMode').split()[0]
        self.pars['communityMode']=self.Parse(self.deconParFile,'communityMode').split()[0]
        self.pars['logicFilter']=self.Parse(self.deconParFile,'logicFilter').split()[0]
        self.pars['ileShiftMode']=self.Parse(self.deconParFile,'ileShiftMode').split()[0]
        self.pars['ileShiftFile']=self.Parse(self.deconParFile,'ileShiftFile')#.split()[0]

        self.pars['craicMode']=self.Parse(self.deconParFile,'craicMode')#.split()[0]
        self.pars['craicFile']=self.Parse(self.deconParFile,'craicFile')#.split()[0]

        self.pars['shiftXFile']=self.Parse(self.deconParFile,'shiftXFile').split()[0]



        if(len(self.pars['ileShiftFile'])>0):
            self.pars['ileShiftMode']="on"
        else:
            self.pars['ileShiftMode']="off"


        if(len(self.pars['craicFile'])>0):
            self.pars['craicMode']="on"
        else:
            self.pars['craicMode']="off"

        if(self.pars['priority_iterations']=='notSpecified'):
            self.pars['priority_iterations']=1
        if(self.pars['maximum_run_time']=='notSpecified'):
            self.pars['maximum_run_time']=100
