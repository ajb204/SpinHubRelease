import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from magma.magma import Magma
import wx.lib.scrolledpanel as scrolled
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

class shift_assignment(wx.Panel):

    def __init__(self,parent):

        wx.Panel.__init__(self, parent=parent, name="Shift Assignment")

    #     self.parent=parent.tabMag
    #     self.inst=Magma(self.parent.deconParFile,run='n') #get instance of magma
    #     self.parent=parent.tabOne.molecule
    #     self.place1=[]
    #     self.node1 = []
    #     self.subgraphs = self.inst.subgraphRef.items()
    #     self.system_graph = self.subgrapher()
    #     self.create_main_panel()
    #     self.draw_figure()
    #
    # def create_main_panel(self):
    #
