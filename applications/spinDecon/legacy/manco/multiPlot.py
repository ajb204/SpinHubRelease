#!/usr/bin/python
import wx,string,numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib     
from matplotlib.figure import Figure
import nmrglue as ng

############################################################################
# Frame for 1d slices
#


matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()

class SliceFrame(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)   

        self.parent=parent
        self.sum=(0.,2.)
        self.peak=tabOne.peak

        self.create_main_panel()
        self.draw_figure()
        #self.Show()
        #self.Fit()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        #self.axes = self.fig.add_subplot(111)

        listy=[]
        for i in range(len(self.peak)):
            listy.append(self.peak[i][0])
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
#        self.dpi = 100
#        self.fig = plt.figure()
#        self.canvas = FigCanvas(self, -1, self.fig)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
#        self.axes = self.fig.add_subplot(110)
#        fig = p.figure()
#        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
#        self.axes = self.fig.add_subplot(111)
#        self.canvas = FigCanvas(self, -1, self.fig)
#        self.ax = p3.Axes3D(fig)
#        self.ax.contour3D(X,Y,Z)
#        p.show()
        # Bind the 'pick' event for clicking on one of the bars
        #
#        self.canvas.mpl_connect('pick_event', self.on_pick)
#        self.textbox = wx.TextCtrl(
#            self,
#            size=(200,-1),
#            style=wx.TE_PROCESS_ENTER)
 #       self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)

        self.drawbutton = wx.Button(self, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.Pbutton = wx.Button(self, -1, "Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        self.Nbutton = wx.Button(self, -1, "Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        #self.cb_grid = wx.CheckBox(self, -1,
        #    "Show Grid",
        #    style=wx.ALIGN_RIGHT)
        #self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        #self.slider_label = wx.StaticText(self, -1,
        #    "Bar width (%): ")
        #self.slider_width = wx.Slider(self, -1,
        #    value=20,
        #    minValue=1,
        #    maxValue=100,
        #    style=wx.SL_AUTOTICKS | wx.SL_LABELS)
        #self.slider_width.SetTickFreq(10, 1)
        #self.Bind(wx.EVT_COMMAND_SCROLL_THUMBTRACK, self.on_slider_width, self.slider_width)
        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)
        #
        # Layout with box sizers
        #
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox.Add(self.ComboBox1, 0, border=3, flag=flags)
        self.hbox.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.Nbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
#        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        #self.hbox.Add(self.text0, 0, flag=flags)
        #self.hbox.Add(self.textbox0, 0, flag=flags)
        #self.hbox.Add(self.text1, 0, flag=flags)
        #self.hbox.Add(self.textbox1, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        self.a0,self.a1,self.aa=self.GetRange('hnca/testJig.ft3')
        self.b0,self.b1,self.bb=self.GetRange('hncoca/testJig.ft3')



    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=string.split(line)
            peak.append(linetosave)
        peakfile.close()
        return peak

    def GetRange(self,infile):
        dic,data=ng.pipe.read(infile)
        uc0 = ng.pipe.make_uc(dic,data,dim=0)
        uc1 = ng.pipe.make_uc(dic,data,dim=1)
        uc2 = ng.pipe.make_uc(dic,data,dim=2)

        Size=data.shape
        uc0max=uc0.ppm(0)
        uc0min=uc0.ppm(Size[0]-1)
        uc1max=uc1.ppm(0)
        uc1min=uc1.ppm(Size[1]-1)
        uc2max=uc2.ppm(0)
        uc2min=uc2.ppm(Size[2]-1)
        
        print("Spectrum dimensions (pts): ",Size)   #print the spectral dimensions
        #print "Labels: ",self.labb
        print("dimension 0 limits (ppm): ", uc0min, uc0max)  #carbon 
        print("dimension 1 limits (ppm): ", uc1min, uc1max)  #direct 
        print("dimension 2 limits (ppm): ", uc2min, uc2max)  #direct 
        return uc0min,uc0max,numpy.fabs(uc0.ppm(1)-uc0.ppm(0))

    def draw_figure(self):
        """ Redraws the figure
        """
        self.axes = self.fig.add_subplot(211)
        self.axes.clear()

        #self.thresh=float(self.textbox0.GetValue())
        #self.offset=float(self.textbox1.GetValue())
        
        dirs='hnca','hncoca'
        for dir in dirs:

            input=self.readfile(dir+'/out/slice2D/'+self.ComboBox1.GetValue()+'.proj.out')
            xs=[]
            ys=[]
            y2s=[]
            for i in range(len(input)):
                xs.append(float(input[i][0]))
                ys.append(input[i][1])
                #y2s.append(float(self.thresh))

            xs=numpy.array(xs)
            ys=numpy.array(ys)
            if(dir=='hnca'):
                #print len(xs),xs.shape
                #print xs,self.b0
                #print xs<self.b0
                select=xs<self.b0
                print(xs[xs<self.b0],self.b0)
                print(numpy.sum(xs<self.b0),'bollocks')

                #print xs[xs<self.b0]
                #xs[xs<self.b0]=xs[xs<self.b0]+(xs[xs<self.b0]-self.b0)
                xs[xs<self.b0]=(self.b1)-(self.b0-xs[xs<self.b0])
                argy=numpy.argsort(xs)
                xs=xs[argy]
                ys=ys[argy]


#                print a0,a1,b0,b1



            #self.axes.set_xlabel(self.parent.tabOne.labb[0],fontsize=8)
            self.axes.plot(xs,ys,label='data')
            #self.axes.plot(xs,y2s,'g',label='threshold')

            self.xmin,self.xmax=self.axes.get_xlim()
            self.ymin,self.ymax=self.axes.get_ylim()
            self.offset=-1*self.ymin/2

            for i in range(len(self.peak)): #write in the peak labels
                self.axes.text(float(self.peak[i][2]),-float(self.offset),self.peak[i][0],fontsize=9,rotation=90)


        self.axes = self.fig.add_subplot(212)
        self.axes.clear()

        #self.thresh=float(self.textbox0.GetValue())
        #self.offset=float(self.textbox1.GetValue())
        
        dirs='hnco','hncaco'
        for dir in dirs:

            input=self.readfile(dir+'/out/slice2D/'+self.ComboBox1.GetValue()+'.proj.out')
            xs=[]
            ys=[]
            y2s=[]
            for i in range(len(input)):
                xs.append(float(input[i][0]))
                ys.append(input[i][1])
                #y2s.append(float(self.thresh))

            xs=numpy.array(xs)
            ys=numpy.array(ys)
            if(dir=='hncaco'):
                #print len(xs),xs.shape
                #print xs,self.b0
                #print xs<self.b0
                select=xs<self.b0
                print(xs[xs<self.b0],self.b0)
                print(numpy.sum(xs<self.b0),'bollocks')

                #print xs[xs<self.b0]
                #xs[xs<self.b0]=xs[xs<self.b0]+(xs[xs<self.b0]-self.b0)
                xs[xs<self.b0]=(self.b1)-(self.b0-xs[xs<self.b0])
                argy=numpy.argsort(xs)
                xs=xs[argy]
                ys=ys[argy]


#                print a0,a1,b0,b1

            #print xs

            #self.axes.set_xlabel(self.parent.tabOne.labb[0],fontsize=8)
            self.axes.plot(xs,ys,label='data')
            #self.axes.plot(xs,y2s,'g',label='threshold')

            self.xmin,self.xmax=self.axes.get_xlim()
            self.ymin,self.ymax=self.axes.get_ylim()
            self.offset=-1*self.ymin/2

            for i in range(len(self.peak)): #write in the peak labels
                self.axes.text(float(self.peak[i][2]),-float(self.offset),self.peak[i][0],fontsize=9,rotation=90)



        """
        if(self.sum[1]==2):
            #plt.title(self.ComboBox1.GetValue()+' automated assignment')
            input=self.readfile('out/fit/'+self.ComboBox1.GetValue()+'.fitslice')
            xs_f=[]
            ys_f=[]
            ys_d=[]
            for i in range(len(input)):
                xs_f.append(input[i][0])
                ys_f.append(input[i][2])
                ys_d.append(float(input[i][2])-float(ys[i])-2*float(self.offset))
            self.axes.plot(xs_f,ys_f,'b',label='fit')
            self.axes.plot(xs_f,ys_d,color='#ADD8E6',label='difference')
            self.axes.legend(fontsize=8)
            for i in range(len(input[i])-3):
                xs=[]
                ys=[]
                for j in range(len(input)):
                    xs.append(input[j][0])
                    ys.append(float(input[j][i+3])-float(self.offset))
                self.axes.plot(xs,ys,'y',label='peak'+str(i))


            self.axes.set_title(self.ComboBox1.GetValue()+' slice')
            input=self.readfile('out/fit/'+self.ComboBox1.GetValue()+'.blur')
            xs_f=[]
            ys_f=[]
            ys_d=[]
            for i in range(len(input)):
                xs_f.append(input[i][0])
                ys_f.append(input[i][2])
            self.axes.plot(xs_f,ys_f,'b',label='noes')

            plt.legend(fontsize=8)
            for i in range(len(input[i])-3):
                xs=[]
                ys=[]
                for j in range(len(input)):
                    xs.append(input[j][0])
                    ys.append(float(input[j][i+3])-float(self.offset))
                self.axes.plot(xs,ys,'y',label='peak'+str(i))



        # clear the axes and redraw the plot anew
        #
#        self.axes.clear()
#        self.axes.grid(self.cb_grid.IsChecked())
#        self.axes.bar(
#            left=x,
#            height=self.data,
#            width=self.slider_width.GetValue() / 100.0,
#            align='center',
#            alpha=0.44,
#            picker=5)

        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        self.axes.yaxis.set_major_formatter(FormatStrFormatter('${%0.0e}$'))
        """

        self.canvas.draw()

    def on_cb_grid(self, event):
        self.draw_figure()

    def on_slider_width(self, event):
        self.draw_figure()

    def on_draw_button(self, event):
        self.draw_figure()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()

    def on_N_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.draw_figure()


    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        dlg = wx.MessageDialog(
            self,
            msg,
            "Click!",
            wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def on_text_enter(self, event):
        self.draw_figure()

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



