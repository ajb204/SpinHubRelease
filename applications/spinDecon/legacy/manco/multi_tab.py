#!/usr/bin/python
####################################################################
# Front end for manco spectral visualisation and fitting software
#
# A.Baldwin 10th Dec 2010
# A.Baldwin 12th June 2015  #sorted out for data

import wx
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting rsoutines from matplotlib     

from .multiPlot2D import SliceFrame2D
from .multiPlot2Dnorm import SliceFrame2Dnorm
from .multiPlot import SliceFrame
from spinDecon.legacy.noe_workspace import NOEFrame
from spinDecon.gui.workspaces.projection import Projection
from spinDecon.gui.workspaces import nmr as deconFrame
from .multi import multiFrame


########################################################################
class NotebookDemo(wx.Notebook):
    """
    Notebook class
    """
    def __init__(self, parent):



        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )


        self.tabOne = multiFrame(self)
        self.AddPage(self.tabOne, "Setup")

        self.test=(710,670)
        self.SetSize(self.test)

        #self.tabOne.OnButtonRead(True)
        #self.tabOne.OnButtonCollate(True)

    def AddTabTwo(self,event,tabOne):
        self.tabTwo = Projection(self,self.tabOne)
        self.AddPage(self.tabTwo, "2Dplane")
        

    def AddTabThree(self,event,tabOne):
        self.tabThree = SliceFrame(self,self.tabOne)
        #        tabOne.SetBackgroundColour("Gray")
        self.AddPage(self.tabThree, "1Ddeconv")

        self.SetSize(self.test)
        self.tabThree.create_main_panel()
        self.tabThree.draw_figure()
            
    def AddTabFour(self,event,tabOne):
        #if(tabOne.cb_grid.IsChecked()==True):
        self.tabFour = SliceFrame2Dnorm(self,self.tabOne)
        #else:
        #    self.tabFour = SliceFrame2Dnorm(self,self.tabOne)
        self.AddPage(self.tabFour,"2Dslices")
            
        #self.tabFive = NOEFrame(self,self.tabOne)
        #self.AddPage(self.tabFive,"NOEstats ")
            
            #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
            #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)
            

        #self.tabTwo.create_main_panel()
        #self.tabTwo.draw_figure()
        
        



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






########################################################################
class MyApp(wx.Frame):
    """
    Frame that holds all other widgets
    """

    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "Manco",wx.DefaultPosition,
                          size=(710,670)
                          )
        panel = wx.Panel(self)

        self.create_menu()
        self.create_status_bar()

        notebook = NotebookDemo(panel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()

        self.Show()



    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
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
        msg = """ A demo using wxPython with matplotlib:

         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw!")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
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



    app = wx.PySimpleApp()
    frame = MyApp()
    app.MainLoop()


