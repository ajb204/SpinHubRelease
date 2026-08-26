#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
##########################################################3
#Front end for assignNMR assignment software
#A.Baldwin 15th April 2019

import os,sys


import wx
from assign.assSetup import assSetup
from assign.singleFrame import singleFrame
from assign.projFrame import projFrame
from assign.magmaAssFrame import magmaFrame
from assign.tocsyFrame import tocsyFrame
from assign.spinSystem import spinSystem

import importlib

#from assignFrame import assFrame


########################################################################
class NotebookAss(wx.Notebook):
    """
    Notebook class
    """
    def __init__(self, parent,panel,assignParFile):
        wx.Notebook.__init__(self, panel, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )

        self.parent=parent
        self.assignParFile=assignParFile

        #try:
        self.tabOne = assSetup(self,self.assignParFile)
        #self.tabOne.SetBackgroundColour(wx.Colour(47,79,79))
        self.tabOne.SetBackgroundColour("WHITE")
        self.AddPage(self.tabOne, "Setup")

        #create dat file if needed
        if(os.path.exists('dat')==False):
            os.system('mkdir dat')

        self.tabMag = magmaFrame(self,'./dat/magmaParFile')
        self.AddPage(self.tabMag, "MAGMA")

        self.test=(1600,1000)
        self.SetSize(self.test)

        #self.tabOne.OnButtonRead(True)
        #self.tabOne.OnButtonCollate(True)

    def AddTabTwo(self):
        try:
            self.KillPage('Walk')
            #print("ba")
            self.tabTwo = singleFrame(self)
            print("ca")
            self.AddPage(self.tabTwo, "Walk")
        except Exception as e:
            print('Failed making tabtwo 2D planes')
            print(e)
            exit()

    def AddTabThree(self):
        # try:
            self.KillPage('Projections')
            self.tabThree = projFrame(self)
            self.AddPage(self.tabThree, "Projections")
        # except Exception as e:
        #     print('Failed making 2D planes')
        #     print(e)
        #     exit()


    def AddTabFour(self):
        return
        self.KillPage('Spin System')

        self.tabFour = spinSystem(self)
        self.AddPage(self.tabFour, "Spin System")
    # def AddTabFive(self):
    #     from shift_assignment import shift_assignment
    #     self.KillPage('Shift Assignment')
    #     self.tabFive = shift_assignment(self)
    #     self.AddPage(self.tabFive, "Shift Assignment")

    def AddTabFive(self):
        self.KillPage('tocsyFrame')
        self.tabFive = tocsyFrame(self)
        self.AddPage(self.tabFive, "tocsyFrame")


    """
    def AddTabThree(self,event,tabOne):
        #try:
        self.KillPage('1Ddeconv')
        self.tabThree = SliceFrame(self,self.tabOne)
        #tabOne.SetBackgroundColour("Gray")
        self.AddPage(self.tabThree, "1Ddeconv")

        #self.SetSize(self.test)
        #self.tabThree.create_main_panel()
        #self.tabThree.draw_figure()
        #except:
        #    print('Failed making projection frame')

    def AddTabFour(self,event,tabOne):
        try:
            #if(tabOne.cb_grid.IsChecked()==True):
            self.KillPage('2Dslices')
            self.tabFour = SliceFrame2D(self,self.tabOne)
            #else:
            #    self.tabFour = SliceFrame2Dnorm(self,self.tabOne)
            self.AddPage(self.tabFour,"2Dslices")

        except:
            print('Failed making 2D slice plots')

    def AddTabFour4D(self,event,tabOne):
        self.tabFour = SliceFrame4D(self,self.tabOne)
        self.AddPage(self.tabFour,"2Dslices")


        #self.tabFive = NOEFrame(self,self.tabOne)
        #self.AddPage(self.tabFive,"NOEstats ")

            #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
            #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)


        #self.tabTwo.create_main_panel()
        #self.tabTwo.draw_figure()
     """

    def OnPageChanged(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        #print('OnPageChanged,  old:%d, new:%d, sel:%d\n' % (old, new, sel))
        event.Skip()

    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        #print('OnPageChanging, old:%d, new:%d, sel:%d\n' % (old, new, sel))
        event.Skip()


    def PageExists(self,pageTitle):
        for index in range(self.GetPageCount()):
            if self.GetPageText(index) == pageTitle:
                return 1
        return 0

    def KillPage(self,pageTitle):
        print
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
class MyAppAss(wx.Frame):
    """
    Frame that holds all other widgets
    """

    #----------------------------------------------------------------------
    def __init__(self,assignParFile,showFlg=True):
        """Constructor"""
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "assignRun",wx.DefaultPosition,
                          #size=(self.monitorWidth*0.95, self.monitorHeight*0.85),
                          size=(1300,670)
                          )
        panel = wx.Panel(self)

        #self.create_menu()
        self.create_status_bar()

        self.SetBackgroundColour('WHITE')

        self.notebook = NotebookAss(self,panel,assignParFile)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)

        panel.SetSizerAndFit(sizer)

        self.Layout()

        if(showFlg):
            self.Show()
        #self.Maximize(True)


    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()
        m_load = menu_file.Append(-1, "&Open\tCtrl-L", "Open session file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_load)
        menu_file.AppendSeparator()

        m_save = menu_file.Append(-1, "&Save\tCtrl-S", "Save status")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_save)
        menu_file.AppendSeparator()


        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)


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
        #msg = """ A demo using wxPython with matplotlib:
        #
        # * Use the matplotlib navigation bar
        # * Add values to the text box and press Enter (or click "Draw!")
        # * Show or hide the grid
        # * Drag the slider to modify the width of the bars
        # * Save the plot to a file using the File menu
        # * Click on a bar to receive an informative message
        #"""
        msg="fun with methyl assignments"
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


    def OnClose(self, event):
        self.Close(True)
        sys.exit(100)





#----------------------------------------------------------------------
if __name__ == "__main__":
    #print(sys.argv)
    if(len(sys.argv)==2):
        assignParFile=sys.argv[1]
    else:
        assignParFile='assignParFile'

    if(os.path.exists(assignParFile)==0):
        outy=open('assignParFile','w');outy.close()

    app = wx.App()
    frame = MyAppAss(assignParFile)


    app.MainLoop()
