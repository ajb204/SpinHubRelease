from .slicePlot2D import sliceFrame2D
from .slicePlot import SliceFrame



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
        tabOne = Tab_PDB(self)
        self.AddPage(tabOne, "InputInformation")
        tabTwo=SliceFrame(self,tabOne)
        self.AddPage(tabTwo, "Deconvolve")

#        tabOne.SetBackgroundColour("Gray")

        tabThree = SliceFrame2D(self,tabOne,tabTwo)
        self.AddPage(tabTwo, "2D representation")


        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)

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





