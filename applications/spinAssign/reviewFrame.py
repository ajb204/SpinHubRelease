#!/usr/bin/python
import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
# matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin


##########################################################################
# 2D plotting of NMR slices
#

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

#def RunFrame(molecule):


class reviewFrameMan(wx.App):
    def __init__(self,inherit):
        app = wx.App()
        #frame = AssFrame(molecule)
        self.frame_ProcessFrame=reviewFrame(inherit)
        #FGA added
        self.frame_ProcessFrame.Centre(direction=wx.BOTH)
        self.frame_ProcessFrame.Show(True)
        app.MainLoop()

class reviewFrame(wx.Frame):


    def __init__(self,parent):

        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "review plots",wx.DefaultPosition,style=wx.DEFAULT_FRAME_STYLE,
                          #size=(self.monitorWidth*0.95, self.monitorHeight*0.85),
                          size=(870,870)
                          )
        panel = wx.Panel(self)
        self.parent=parent

        self.create_main_panel()
        self.draw_figure()
        self.Show(True)
        self.Fit()
        


    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """

        #self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.


        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)

        # Bind the 'pick' event for selection
        self.canvas.mpl_connect('button_press_event', self.on_pick)

        #
        #
        # 1. draw a series of 1D spectra
        # 2. have a button to run scipy on all traces
        # 3. have a listcntrl that shows possible peaks.
        # 4. have the ability to select just one and add.
        # 5. on a second listcntrl, we have our current peaks.
        # 6. need to be able to edit them all.

        
    def draw_figure(self):

        #

        pass
