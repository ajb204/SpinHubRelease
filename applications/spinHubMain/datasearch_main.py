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
import os,sys, platform
pathname, scriptname = os.path.split(sys.argv[0])   #get location where this script was executed
if(os.path.exists(os.pathsep+os.path.join(os.getcwd(),'bin') )): #does this path exist?
    os.environ["PATH"]+=os.pathsep+os.path.join(os.getcwd(),'bin')  #if running from an app, this will add bins to the system path
if(len(os.path.dirname(sys.executable).split('deconRun.app'))>1):
    from os.path import expanduser  #if running using the app, change working folder to user directory
    os.chdir(expanduser("~"))

#adding temp location.
#binaries will be copied here, so needs to be in system's path
#only for pyinstaller linux app.
if(platform.uname()[0]=='Linux'):
    try:
        print('MEIPASS:',sys._MEIPASS)
        os.environ["PATH"]+=os.pathsep+sys._MEIPASS
    except:
        pass
    try: #cleanup files in tmp
        files=os.listdir('/tmp')
        for file in files:
            if(len(file.split('MEI'))>1):
                test=os.path.join('/tmp',file)
                if test!=sys._MEIPASS:
                    print('Removing temp file:',test)
                    os.system('rm -rf '+test)
                else:
                    print('this is our guy',test)
    except:
        pass



#cleanup MEIPASS
#removing old temp directories
#import subprocess
#subprocess.call(['ls','-l'])

# Begin importing
import wx
from .Frames.dataSearch import dataSearch



########################################################################
class NotebookDemo(wx.Notebook):
    """
    Notebook class

    """
    def __init__(self, parent,panel):
        wx.Notebook.__init__(self, panel, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )
        self.parent=parent
        self.AddTab_datasearch()


    def AddTab_datasearch(self):
        #try:
        self.KillPage('dataSearch')
        self.tab_datasearch = dataSearch(self)
        #self.tab_datasearch.Bind(wx.EVT_SET_FOCUS, self.tab_datasearch.onFocus)
        self.AddPage(self.tab_datasearch, "dataSearch")

    def OnPageChanged(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        #print 'OnPageChanged,  old:%d, new:%d, sel:%d\n' % (old, new, sel)
        event.Skip()

    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        #print 'OnPageChanging, old:%d, new:%d, sel:%d\n' % (old, new, sel)
        event.Skip()

    def PageExists(self,pageTitle):
        for index in range(self.GetPageCount()):
            if self.GetPageText(index) == pageTitle:
                return 1
        return 0

    def KillPage(self,pageTitle):
        print()
        print('Pages:',self.GetPageCount())
        for index in range(self.GetPageCount()):
            print(self.GetPageText(index))
            if self.GetPageText(index) == pageTitle:
                print('killing page')
                print(self.GetPageCount())
                self.DeletePage(index)
                self.SendSizeEvent()
                print(self.GetPageCount())
                print('done')
                break



########################################################################
class MyApp(wx.Frame):
    """
    Frame that holds all other widgets
    """

    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "dataSearch - "+os.getcwd().split('/')[-1], wx.DefaultPosition,
                        #   size=(self.monitorWidth*0.75, self.monitorHeight*0.85),
                          size=(1370,780)
                          )
        panel = wx.Panel(self)

        self.create_menu()
        self.create_status_bar()

        self.SetBackgroundColour('WHITE')

        self.notebook = NotebookDemo(self,panel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)

        panel.SetSizerAndFit(sizer)
        # self.Bind(wx.EVT_SIZE, self.OnSize)

        self.Layout()

        self.Show()

        #if(os.path.exists(deconParFile)==1):
        #    self.TestPath(deconParFile)
        #    self.DoLoad(deconParFile)
        #    #self.DoLoad(os.path.join(os.getcwd(), deconParFile))

        #self.Maximize(True)

    def OnSize(self, event):
        print('resized!')
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()

        m_new = menu_file.Append(-1, "&New\tCtrl-N", "New session")
        self.Bind(wx.EVT_MENU, self.OnNew, m_new)
        menu_file.AppendSeparator()

        m_load = menu_file.Append(-1, "&Open\tCtrl-L", "Open session file")
        self.Bind(wx.EVT_MENU, self.OnLoadResults, m_load)
        menu_file.AppendSeparator()

        m_save = menu_file.Append(-1, "&Save\tCtrl-S", "Save status")
        self.Bind(wx.EVT_MENU, self.OnSaveResults, m_save)
        menu_file.AppendSeparator()


        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.OnQuit, m_exit)
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def OnQuit(self, e):
        self.Destroy()

    def OnNew(self,event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            #defaultFile=os.path.split(self.deconParFile)[1],
            defaultFile='deconParFile',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.deconParFile=path

            outy=open(self.deconParFile,'w');outy.close()

            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Loaded %s" % path)
            os.chdir(os.path.dirname(path))
            print("CWD: ",os.getcwd())

            self.DoLoad(path)

    def OnSaveResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            #defaultFile=os.path.split(self.deconParFile)[1],
            defaultFile='deconParFile',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            self.deconParFile=path
            if(os.path.exists(self.deconParFile)==0):
                outy=open(self.deconParFile,'w');outy.close()
            # self.notebook.tabMagma.deconParFile=path
            # self.notebook.tabMagma.OnButtonSave(True)

            self.notebook.tabOne.deconParFile=path
            self.notebook.tabOne.OnButtonSave(True)

            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved %s" % path)


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
            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Loaded %s" % path)
            os.chdir(os.path.dirname(path))
            # print("CWD: ",os.getcwd())

            self.DoLoad(path)
    def TestPath(self,deconParFile):
        Parse=self.notebook.tabOne.Parse #messy...

        self.deconParFile=deconParFile
        indir=Parse(self.deconParFile,'indir')
        fiddir=Parse(self.deconParFile,'fiddir')
        # print('start:',indir,fiddir)
        if(indir!=0):
            if(os.path.exists(str(indir))==0):
                indir=self.CheckPath(str(indir))
            if fiddir!=0:
                if(os.path.exists(str(fiddir))==0):
                    fiddir=self.CheckFidPath(str(indir),str(fiddir))

        #self.WriteFID(indir,fiddir)
        # print('finish:',indir,fiddir)

    def CheckPath(self,indir):
        print('cannot find ',indir,'. Trying to update:')
        tast=indir.split("/")
        loop=len(tast)-1
        ref=self.deconParFile

        for i in range(len(tast)): #looping backwards along the files
            ii=loop-i
            test=os.path.join(os.getcwd(),indir.split("/")[ii],ref)
            #print 'testing:',test
            if(os.path.exists(test)==1):
                print('Found new indir:',os.path.join(os.getcwd(),indir.split("/")[ii]))
                #sys.exit(100)
                #os.setcwd(indir)
                # print(os.getcwd())
                return os.path.join(os.getcwd(),indir.split("/")[ii])

        print('Cannot find directory',indir)
        sys.exit(100)

    def CheckFidPath(self,indir,fiddir):
        print('cannot find fiddir: ',fiddir,'. Trying to update:')
        #test=fiddir.split(indir)
        tast=fiddir.split("/")
        print(indir)
        for i in range(len(tast)): #looping backwards along the files
            ii=len(tast)-1-i-1

            splitty=tast[ii] #point to split
            print(splitty)
            click=os.path.join(indir,splitty)
            #try:
            #    click=os.path.join(splitty,fiddir.split(splitty)[-1])
            #except:
            #    click=os.path.join(splitty,fiddir.split(splitty)[-2])
            print(click)
            test=os.path.join(indir,click)
            # print('testing:',test)
            if(os.path.exists(test)==1):
                print('Found new fiddir:',test)
                #sys.exit(100)
                return test
        print('Cannot find directory',indir)
        sys.exit(100)


        print(test)
        print(indir,fiddir)
        fidnew=os.path.join(indir,test[-1])
        if(os.path.exists(fidnew)):
            print('Newfid found:',fidnew)
            return fidnew
        print('Cannot update fidfile.')
        return str(0)


    def WriteFID(self,indir,fiddir):
        try:
            decfile=os.path.join(indir,self.deconParFile)
        except:
            return
        if(os.path.exists(decfile)==0):
            return
        dec=[]
        inny=open(decfile)
        for line in inny.readlines():
            test=line.split()
            tick=0
            if(len(test)>0):
                if(test[0]=='fiddir'):
                    dec.append('fiddir = '+fiddir+'\n')
                    tick=1
                if(test[0]=='indir'):
                    dec.append('indir = '+fiddir+'\n')
                    tick=1
            if(tick==0):
                dec.append(line)
        inny.close()
        outy=open(decfile,'w')
        for de in dec:
            outy.write(de)
        outy.close()



    def DoLoad(self,path):
        #load magmaTab
        self.notebook.deconParFile=path

        #self.notebook.tabMagma.deconParFile=path
        #self.notebook.tabMagma.GetPars()
        #self.notebook.tabMagma.UpdatePars()
        #self.notebook.tabMagma.dirVal.SetValue(os.path.dirname(path))

        #load nmr tab
        #self.notebook.tabOne.deconParFile=path
        #self.notebook.tabOne.dirBox.SetValue(os.path.dirname(path))
        #self.notebook.tabOne.dirBox.SetValue(path)
        self.notebook.tabOne.OnButtonLoad(True)
        self.notebook.tabOne.pre_read_disabling()#load values in input file

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
        msg="UniDecNMR"
        dlg = wx.MessageDialog(self, msg, "UniDecNMR", wx.OK)
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


    def OnClose(self, event):
        self.Close(True)
        sys.exit(100)





#----------------------------------------------------------------------
if __name__ == "__main__":
    #print sys.argv
    if(len(sys.argv)==2):
        deconParFile=sys.argv[1]
    else:
        deconParFile='deconParFile'

    if(os.path.exists(deconParFile)==0):
        outy=open('deconParFile','w');outy.close()

    app = wx.App()
    frame = MyApp(deconParFile)
    app.MainLoop()
