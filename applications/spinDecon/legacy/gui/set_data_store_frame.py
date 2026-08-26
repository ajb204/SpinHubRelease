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

import numpy,sys,os,string,math,wx
import nmrglue as ng
from spinDecon.misc import textEdit
from spinDecon.gui.dialogs.errors import errorMessage

from . import deconFrame as decon
from spinDecon.parameter_store import update_parameter_file

import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import nmrglue as ng
from matplotlib.figure import Figure
import posixpath
from wx.lib.mixins.listctrl import ColumnSorterMixin
from matplotlib.widgets import Cursor
import logging



class SetPathsMan(wx.App):
    def __init__(self,inherit,showFlg=True):
        self.frame_SetPathsFrame=SetPathsFrame(None,40,'Set Datastore',inherit,showFlg=showFlg)
        if(showFlg):
            self.frame_SetPathsFrame.Show(True)
#        return Frame1(parent)

class SetPathsFrame(wx.Frame):
    def __init__(self,parent,id,title,inherit,showFlg=True):
        wx.Frame.__init__(self,  name='', parent=parent,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title='Set Datastore Paths ...')
        self.SetClientSize(wx.Size(900, 280))
        self.SetBackgroundColour(wx.Colour(255, 255, 255))

        self.parent=inherit
        self.state = getattr(inherit, "state", getattr(getattr(inherit, "parent", None), "state", None))
        panel=wx.Panel(self,-1)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.create_main_panel()

        if(showFlg):
            self.Show(True)

        self.Fit()
        self.Show(True)


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.splitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.vbox=wx.BoxSizer(wx.VERTICAL)

        self.DataStoreLab = wx.StaticText(self, label="Data Store:")
        self.DataStoreBox = wx.TextCtrl(self, size=(400, 22))
        self.openDataStoreBtn = wx.Button(self, label="...", size=(40,22))
        self.FidPathLab = wx.StaticText(self, label="Fid Path:")
        self.FidPathBox = wx.TextCtrl(self, size=(400, 22))
        self.openFidPathBtn = wx.Button(self, label="...", size=(40,22))

        self.StoreLab = wx.StaticText(self, label="Data Path:")
        self.StorePath = wx.StaticText(self, label="")
        self.StoreExists = wx.StaticText(self, label="")

        
        self.openDataStoreBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDirL(evt, self.DataStoreBox))
        self.openFidPathBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDirL(evt, self.FidPathBox,default=self.DataStoreBox.GetValue() ))

        self.LoadButton = wx.Button(self, label="Load", size=(100,22))
        self.SaveButton = wx.Button(self, label="Save", size=(100,22))
        self.CloseButton = wx.Button(self, label="Close", size=(100,22))
        self.LocalButton = wx.Button(self, label="Set Local", size=(100,22))

        
        self.LoadButton.Bind(wx.EVT_BUTTON,self.OnLoadButton)
        self.SaveButton.Bind(wx.EVT_BUTTON,self.OnSaveButton)
        self.CloseButton.Bind(wx.EVT_BUTTON,self.OnClose)
        self.LocalButton.Bind(wx.EVT_BUTTON,self.OnLocalButton)
        

        self.sizer = wx.GridBagSizer(2, 3)

        cnt=0
        self.sizer.Add(self.DataStoreLab,(cnt,0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.DataStoreBox,(cnt,1), flag=wx.ALIGN_CENTER_VERTICAL);
        self.sizer.Add(self.openDataStoreBtn,(cnt,2), flag=wx.ALIGN_CENTER_VERTICAL);cnt+=1

        self.sizer.Add(self.FidPathLab,(cnt,0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.FidPathBox,(cnt,1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.openFidPathBtn,(cnt,2), flag=wx.ALIGN_CENTER_VERTICAL);cnt+=1

        self.sizer.Add(self.StoreLab,(cnt,0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.StorePath,(cnt,1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.StoreExists,(cnt,2), flag=wx.ALIGN_CENTER_VERTICAL);cnt+=1
        

        #self.StoreLab = wx.StaticText(self, label="Data Path:")
        #self.StorePath = wx.StaticText(self, label="")
        #self.StoreExists = wx.StaticText(self, label="")

        
        self.buttonSz = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonSz.Add(self.LoadButton)
        self.buttonSz.Add(self.SaveButton)
        self.buttonSz.Add(self.CloseButton)
        self.buttonSz.Add(self.LocalButton)
        

        self.dataLbl = wx.StaticBox(self, -1, 'Set Datastore Paths')
        self.dataSizer = wx.StaticBoxSizer(self.dataLbl, wx.VERTICAL)
        self.dataSizer.Add(self.sizer)

        self.fullSz = wx.BoxSizer(wx.VERTICAL)
        self.fullSz.Add(self.dataSizer)
        self.fullSz.Add(self.buttonSz)
        #self.dataSizer.Add(self.buttonSz)
        
        #self.vbox.Add(self.DataStoreLab)
        #self.splitSizer.Add(self.vbox, 10, flag=wx.GROW)
        self.SetSizerAndFit(self.fullSz)        

        self.OnLoadButton(True)
        
    def OnLoadButton(self,event):
        print("Loading datastore paths...")
        defaultDir=''
        DeconDataStore=str(self.parent.parent.Parse(self.parent.parent.deconParFile,'dataStore',default=defaultDir))
        if(DeconDataStore==''):
           pass
        elif(os.path.exists(DeconDataStore)==True):  #if the path exists, leave it there.
           pass
        elif(os.path.exists(os.path.join(os.path.expanduser("~"),DeconDataStore))==True): #try to  combine path with home
            DeconDataStore=os.path.join(os.path.expanduser("~"),DeconDataStore) #add home to the datapath store
        else:
             print("Datastore not found in deconParFile")
        print('DeconDataStore:',DeconDataStore)
        #print(os.path.expanduser("~"))
        #print(os.path.join(os.path.expanduser("~"),DeconDataStore))
        #print(os.path.exists(os.path.join(os.path.expanduser("~"),DeconDataStore)))
        self.DataStoreBox.SetValue(DeconDataStore)
        self.FidPathBox.SetValue(str(self.parent.parent.Parse(self.parent.parent.deconParFile,'dataLoc',default=defaultDir)))
        print( 'dataLoc', str(self.parent.parent.Parse(self.parent.parent.deconParFile,'dataLoc',default=defaultDir)) )

        self.GetStorePath()

        #self.StorePath.show()
        #self.StorePath.SetName(self.storePath)
        #print(self.StorePath.GetName())

    def GetStorePath(self):
        self.storePath=os.path.join(self.DataStoreBox.GetValue(),self.FidPathBox.GetValue())
        self.storeExist=os.path.exists(self.storePath)
        if(self.storeExist):
            self.parent.GetSpectrometerType()

            if(self.parent.tp=='bruk'):
                sp='Bruker'
            elif(self.parent.tp=='var'):
                sp='Varian'
            else:
                sp='Not sure of spectrometer type'
            
            self.StorePath.SetLabel("Path exists: "+sp)
        else:
            self.StorePath.SetLabel("Path does not exist!")
        print("DataLocation:",self.storePath)
        print("Does the store exist?",self.storeExist)


    def OnLocalButton(self,event):
        dirbox=self.parent.dirBox.GetValue()
        self.DataStoreBox.SetValue('')
        self.FidPathBox.SetValue(dirbox)
        self.GetStorePath()
        pass
        
    def OnSaveButton(self,event):
        savefile=os.path.join(self.parent.parent.dirBox.GetValue(),self.parent.parent.deconParFile)
        print('Saving to:',savefile)
        print("DataStore:",self.DataStoreBox.GetValue())
        print("DataLoc:",self.FidPathBox.GetValue())

        self.GetStorePath()
        if(self.storeExist==False):
            print("File store does not exist! Aborting save!")
            return
        
        write={}
        write['dataStore']=self.DataStoreBox.GetValue()
        write['dataLoc']=self.FidPathBox.GetValue()
        update_parameter_file(savefile, write, source_path=self.parent.parent.deconParFile)

    def OnClose(self, event):
        
        """
        if event.CanVeto() and self.fileNotSaved:

            if wx.MessageBox("Would you like to save?",
                            "Before exiting",
                            wx.ICON_QUESTION | wx.YES_NO) == wx.YES:
                self.OnButtonSave(None)
        """

            
        self.Destroy()  # you may also do:  event.Skip()

