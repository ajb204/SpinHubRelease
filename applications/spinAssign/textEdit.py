import os
import wx

class MyPanel(wx.Panel):

    def __init__(self, parent,run='n'):

        self.parent=parent
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour('GREY')
        self.my_text = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        btn = wx.Button(self, label='Open Text File')
        btn.Bind(wx.EVT_BUTTON, self.onOpen)

        btnC = wx.Button(self, label='Close')
        btnC.Bind(wx.EVT_BUTTON, self.onClose)

        if(run=='y'):
            btnR = wx.Button(self, label='Run')
            btnR.Bind(wx.EVT_BUTTON, self.onRun)
            


        sizer = wx.BoxSizer(wx.VERTICAL)


        sizer.Add(self.my_text, 1, wx.ALL|wx.EXPAND)

        btsizer=wx.BoxSizer(wx.HORIZONTAL)
        btsizer.Add(btn, 0, wx.ALL|wx.CENTER, 5)
        btsizer.Add(btnC, 0, wx.ALL|wx.CENTER, 5)
        if(run=='y'):
            btsizer.Add(btnR, 0, wx.ALL|wx.CENTER, 5)
        sizer.Add(btsizer)

        self.SetSizer(sizer)

    def onOpen(self, event):
        wildcard = "*"
        dialog = wx.FileDialog(self, "Open Text Files", wildcard=wildcard, defaultFile="",
                               style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if dialog.ShowModal() == wx.ID_CANCEL:
            return

        path = dialog.GetPath()

        if os.path.exists(path):
            with open(path) as fobj:
                for line in fobj:
                    self.my_text.WriteText(line)

    def GoOpen(self,infile):
        self.infile=infile
        path = infile
        if os.path.exists(path):
            with open(path) as fobj:
                for line in fobj:
                    self.my_text.WriteText(line)

    def GoStream(self,stream):
        print(stream)
        for line in stream:
            for lon in line:
                self.my_text.WriteText(lon+'\n')


    def onClose(self, event):
        self.parent.Close()

    def onRun(self, event):
        print('Overwriting ',self.infile)
        outy=open(self.infile,'w')
        lines=self.my_text.GetNumberOfLines()
        for i in range(lines):
            outy.write(self.my_text.GetLineText(i)+'\n')
        outy.close()
        os.system('csh '+self.infile)
        pass


class MyFrame(wx.Frame):

    def __init__(self,infile,run='n',stream='n'):
        wx.Frame.__init__(self, None, title='Text File Reader')
        panel = MyPanel(self,run=run)

        print('stream:',stream)
        print(infile)
        if(stream=='n'):
            panel.GoOpen(infile)
        else:
            panel.GoStream(infile)
        self.Show()

if __name__ == '__main__':
    app = wx.App(False)
    frame = MyFrame('results/1D09/log')
    app.MainLoop()
