#!/usr/bin/python
import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import nmrglue as ng
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin

##########################################################################
# 2D plotting of NMR slices
#



def RunFrame(uc1min,uc1max,peak,noiseVal,conn_data,spectrumfile):
    app = wx.PySimpleApp()
    offset=0.0
    thresh=noiseVal
    index_data=index(peak)                         #extract peak index from peak list
    frame = SliceFrame2D(uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile)
    app.MainLoop()


class SliceFrame2Dnorm(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent,tabOne):
    #def __init__(self,uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile):
        wx.Panel.__init__(self, parent=parent)

        self.parent=parent
        self.tabOne=parent.tabOne
        
        #wx.Frame.__init__(self, None, -1, self.title)
        #copy in the previous variables
        self.index_data=self.index(parent.tabOne.peak)
        self.thresh=parent.tabOne.noiseVal*float(parent.tabOne.threshBox.GetValue())
        #self.offset=copy.deepcopy(tabOne.offset)
        self.offset=0.0
        self.peak=parent.tabOne.peak
        self.conn_data=parent.tabOne.conn_data
        self.spectrumfile=parent.tabOne.spectrumfile
        #get 2d strips from 3d data
        self.GetSlice2d(self.spectrumfile)


        self.pick_cnt=0
        self.selection=[]
        self.xmin=0
        self.xmax=0
        self.ymin=0
        self.ymax=0
        self.ax_reset=1       #for keeping the zoom
        self.inc=0            #for incrementing the slices

        self.xmin2=0
        self.xmax2=0
        self.ymin2=0
        self.ymax2=0
        self.ax_reset2=1
        self.inc2=0

        self.NOE=0


        dmax=tabOne.uc1max
        dmin=tabOne.uc1min



        self.create_main_panel()
        self.draw_figure()
        self.canvas.draw()
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

        #min max and lvls for slices
        self.text_slice=wx.StaticText(self, -1, 'Slices:',size=(80,-1))
        self.text1=wx.StaticText(self, -1, 'Min:')
        self.text2=wx.StaticText(self, -1, 'Factor:')
        self.text3=wx.StaticText(self, -1, 'Number:')
        self.textbox0 = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox1 = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox2 = wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)

        #min max and lvls for projections
        self.text_proj=wx.StaticText(self, -1,  'Projections:',size=(80,-1))
        self.text1p=wx.StaticText(self, -1, 'Min:')
        self.text2p=wx.StaticText(self, -1, 'Factor:')
        self.text3p=wx.StaticText(self, -1, 'Number:')
        self.textbox_minP = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox_maxP = wx.TextCtrl(self,size=(100,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox_lvlP = wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)

        #set the default values
        self.textbox0.SetValue(str(self.thresh))
        self.textbox1.SetValue(str(1.1))
        self.textbox2.SetValue(str(30))

        self.textbox_minP.SetValue(str(self.thresh))
        self.textbox_maxP.SetValue(str(1.2))
        self.textbox_lvlP.SetValue(str(30))


        
        self.text_pickfac=wx.StaticText(self, -1, 'CrossPeak cutoff:')
        self.text_savelist=wx.StaticText(self, -1, 'Save list:')

        self.textbox_pickfac = wx.TextCtrl(self,size=(50,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox_savelist = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER)
        self.textbox_pickfac.SetValue(str(3.0))
        self.textbox_savelist.SetValue(str('out/cross_man_save.out'))



        listy=[]
        for i in range(len(self.peak)):
            listy.append(self.peak[i][0])
        self.ComboBox1=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox1.SetSelection(0)


        #noes=[]
        #for i in range(len(self.peak)):
        #    noes.append(self.peak[i][0])
        #self.ComboBoxNoe=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=noes, style=wx.CB_READONLY)
        #self.ComboBoxNoe.SetSelection(0)

        #self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        #self.ComboBox2.SetSelection(0)





        self.Nbutton = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        #self.Nbutton2 = wx.Button(self, -1,"Next")
        #self.Bind(wx.EVT_BUTTON, self.on_N_button2, self.Nbutton2)

        #self.Pbutton2 = wx.Button(self, -1,"Previous")
        #self.Bind(wx.EVT_BUTTON, self.on_P_button2, self.Pbutton2)


        self.Upbutton = wx.Button(self, -1,"Up")
        self.Bind(wx.EVT_BUTTON, self.on_Up_button, self.Upbutton)

        self.Downbutton = wx.Button(self, -1,"Down")
        self.Bind(wx.EVT_BUTTON, self.on_Down_button, self.Downbutton)

        self.NOEbutton = wx.Button(self, -1,"NOE")
        self.Bind(wx.EVT_BUTTON, self.on_NOE_button, self.NOEbutton)



        #self.Upbutton2 = wx.Button(self, -1,"Up")
        #self.Bind(wx.EVT_BUTTON, self.on_Up_button2, self.Upbutton2)
        #
        #self.Downbutton2 = wx.Button(self, -1,"Down")
        #self.Bind(wx.EVT_BUTTON, self.on_Down_button2, self.Downbutton2)


        self.drawbutton = wx.Button(self, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)


        self.cb_grid = wx.CheckBox(self, -1,
            "Peak positions",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)

        self.cb_grid_auto = wx.CheckBox(self, -1,
            "Crosspeak labels",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid_auto, self.cb_grid_auto)
        self.cb_grid_auto.SetValue(1)


        self.cb_grid_select = wx.CheckBox(self, -1,
            "Select",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid_auto, self.cb_grid_auto)


        self.searchbutton = wx.Button(self, -1, "Search")
        self.Bind(wx.EVT_BUTTON, self.on_search_button, self.searchbutton)

        self.deletebutton = wx.Button(self, -1, "Delete")
        self.Bind(wx.EVT_BUTTON, self.on_delete_button, self.deletebutton)

        self.deselectbutton = wx.Button(self, -1, "Deselect")
        self.Bind(wx.EVT_BUTTON, self.on_deselect_button, self.deselectbutton)

        self.savebutton = wx.Button(self, -1, "Save List")
        self.Bind(wx.EVT_BUTTON, self.on_save_button, self.savebutton)


        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)

        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.ComboBox1, 0, border=3, flag=flags)
        self.hbox.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.Nbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.Upbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.Downbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.NOEbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.ComboBox2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Pbutton2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Nbutton2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Upbutton2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Downbutton2, 0, border=3, flag=flags)

        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid_auto, 0, border=3, flag=flags)
        self.hbox.Add(self.text_pickfac, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_pickfac, 0, border=3, flag=flags)


#        self.hbox.Add(self.cb_grid2, 0, border=3, flag=flags)
#        self.hbox.AddSpacer(30)
#        self.hbox.Add(self.slider_label, 0, flag=flags)
#        self.hbox.Add(self.slider_width, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.text_proj,0, border=3, flag=flags)
        self.hbox.Add(self.text1p, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_minP, 0, border=3, flag=flags)
        self.hbox.Add(self.text2p, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_maxP, 0, border=3, flag=flags)
        self.hbox.Add(self.text3p, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_lvlP, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)



        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.text_slice, 0, border=3, flag=flags)
        self.hbox.Add(self.text1, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox0, 0, border=3, flag=flags)
        self.hbox.Add(self.text2, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox1, 0, border=3, flag=flags)
        self.hbox.Add(self.text3, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox2, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.cb_grid_select, 0, border=3, flag=flags)
        self.hbox.Add(self.searchbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.deselectbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.deletebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.text_savelist, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_savelist, 0, border=3, flag=flags)
        self.hbox.Add(self.savebutton, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        #self.panel.SetSizer(self.vbox)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        input=self.readfile('out/slice2D/xy.out')#the CC plane
        self.xs_xy=[]
        self.ys_xy=[]
        self.zs_xy=[]
        self.Xs_xy=[]
        self.Ys_xy=[]
        self.Zs_xy=[]
        for i in range(len(input)):
            if(len(input[i])!=0):
                self.xs_xy.append(float(input[i][0]))
                self.ys_xy.append(float(input[i][1]))
                self.zs_xy.append(float(input[i][2]))
            else:
                self.Xs_xy.append(self.xs_xy)
                self.Ys_xy.append(self.ys_xy)
                self.Zs_xy.append(self.zs_xy)
                self.zs_xy=[]
                self.ys_xy=[]
                self.xs_xy=[]

        input=self.readfile('out/slice2D/yz.out')#the HC plane
        self.xs_yz=[]
        self.ys_yz=[]
        self.zs_yz=[]
        self.Xs_yz=[]
        self.Ys_yz=[]
        self.Zs_yz=[]
        for i in range(len(input)):
            if(len(input[i])!=0):
                self.xs_yz.append(float(input[i][0]))
                self.ys_yz.append(float(input[i][1]))
                self.zs_yz.append(float(input[i][2]))
            else:
                self.Xs_yz.append(self.xs_yz)
                self.Ys_yz.append(self.ys_yz)
                self.Zs_yz.append(self.zs_yz)
                self.zs_yz=[]
                self.ys_yz=[]
                self.xs_yz=[]


    #make an index
    def index(self,array):
        index=[]
        for i in range(len(array)):
            index.append((array[i][0]))
        return index

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=string.split(line)
            peak.append(linetosave)
        peakfile.close()
        return peak

    def findnear_index(self,test,array):
        inear=0
        itest=math.fabs(float(array[0])-float(test))
        for i in range(len(array)):
            if(math.fabs(float(array[i])-float(test))<itest):
                itest=math.fabs(float(array[i])-float(test))
                inear=i
        return int(inear)

    #get 2d strips from 3d data
    def GetSlice2d(self,spectrumfile):

        # Prefer the shared in-memory spectrum if it is already loaded.
        shared_data = getattr(self.tabOne, "data", None)
        shared_dic = getattr(self.tabOne, "dic", None)
        if shared_data is not None and shared_dic is not None and len(getattr(shared_data, "shape", ())) >= 3:
            dic = shared_dic
            self.data = shared_data
        else:
            dic,self.data = ng.pipe.read(spectrumfile)
        uc0 = ng.pipe.make_uc(dic,self.data,dim=0)
        uc1 = ng.pipe.make_uc(dic,self.data,dim=1)
        uc2 = ng.pipe.make_uc(dic,self.data,dim=2)
        Size=self.data.shape
        #print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
        #print "dimension 0 limits (ppm): ", uc0.ppm(0), uc0.ppm(Size[0]-1)  #carbon (max/min)
        #print "dimension 1 limits (ppm): ", uc1.ppm(0), uc1.ppm(Size[1]-1)  #carbon  (max/min)
        #print "dimension 2 limits (ppm): ", uc2.ppm(0), uc2.ppm(Size[2]-1)  #proton (max/min)

        self.index0=[]#make index of carbon chemical shifts for index 0
        for i in range((Size[0])):
            self.index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
        self.index1=[]#make index of carbon chemical shifts for index 1
        for i in range((Size[1])):
            self.index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
        self.index2=[]#make index of carbon chemical shifts for index 2
        for i in range((Size[2])):
            self.index2.append((uc2.ppm(0)-i*(-uc2.ppm(Size[2]-1)+uc2.ppm(0))/(Size[2]-1)))


    #get 2d strips from 3d data
    def ReSlice2d(self,inc,pkl,peak,width):

        #print out 2D slice for each peak correlation
        #print "Extracting slices from ",peak[pkl][0],"     proton:  ",peak[pkl][1], "ppm      carbon:  ",peak[pkl][2],"ppm of width ",width
        ptC=self.findnear_index(float(peak[pkl][2]),self.index1)#find the nearest point to desired chemical shift in carbon index
        ptC=ptC+inc
    
        ptH=self.findnear_index(float(peak[pkl][1]),self.index2)#find the nearest point to desired chemical shift in carbon index
        ptH_max=self.findnear_index(float(peak[pkl][1])+float(width)/2,self.index2)#find the nearest point to desired chemical shift in carbon index
        ptH_min=self.findnear_index(float(peak[pkl][1])-float(width)/2,self.index2)#find the nearest point to desired chemical shift in carbon index
        slice2d=self.data[:,ptC,ptH_max:ptH_min] #extract the relevant 2d slice
    
        Xs=[]
        Ys=[]
        Zs=[]
        for e in range(len(slice2d[0])):  #proton
            xs=[]
            ys=[]
            zs=[]
            for d in range(len(slice2d)): #carbon
                xs.append(self.index0[d])
                ys.append(self.index2[ptH_max+e])
                zs.append(slice2d[d][e])
            Xs.append(xs)
            Ys.append(ys)
            Zs.append(zs)

        #print 'Done!'
        return Xs,Ys,Zs



    #get 2d strips from 3d data
    def ReSlice2dNOE(self,i1,peak,width):

        ptC=self.findnear_index(float(peak[i1][2]),self.index1)#find the nearest point to desired chemical shift in carbon index
        ptH=self.findnear_index(float(peak[i1][1]),self.index2)#find the nearest point to desired chemical shift in carbon index

            
        ptC_min=self.findnear_index(float(peak[i1][2])-float(width)/2,self.index1)#find the nearest point to desired chemical shift in carbon index
        ptC_max=self.findnear_index(float(peak[i1][2])+float(width)/2,self.index1)#find the nearest point to desired chemical shift in carbon index


        slice2d=self.data[:,ptC_max:ptC_min,ptH] #extract the relevant 2d slice
    

        Xs=[]
        Ys=[]
        Zs=[]
        for e in range(len(slice2d[0])):  #proton
            xs=[]
            ys=[]
            zs=[]
            for d in range(len(slice2d)): #carbon
                xs.append(self.index0[d])
                ys.append(self.index1[ptC_max+e])
                zs.append(slice2d[d][e])
            Xs.append(xs)
            Ys.append(ys)
            Zs.append(zs)

        #print 'Done!'
        return Xs,Ys,Zs



    def draw_figure(self):
        """ Redraws the figure
        """
        colormap=cm.seismic


        CarbonWidth=2.0
        
        if(1==1):
            self.axes = self.fig.add_subplot(321)
            self.axes.clear()
            max_level=float(self.textbox1.GetValue())#set contour max level from box
            min_level=float(self.textbox0.GetValue())#set contour min level from box
            ctr_level=int(self.textbox2.GetValue())  #set the cnumber of contours from box                
            #levels=[]
            #for i in range(ctr_level):
            #    levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))
            levels=[] 
            levels.append(min_level)
            for i in range(ctr_level-1):
                levels.append(levels[i]*max_level)
            levels=numpy.array(levels)
            levels=numpy.concatenate((-1*levels[::-1],levels))
            self.axes.contour( self.Xs_xy, self.Ys_xy, self.Zs_xy,levels,cmap=colormap) #plot pdb network
            #find max and min y range
            #y_max=Ys[0][0]
            #y_min=Ys[0][(len(Ys[0]))-1]
            x_max=self.Xs_xy[0][0]
            x_min=self.Xs_xy[(len(self.Xs_xy))-1][0]

            self.axes.set_xlim(x_max,x_min)
            
            carbonline=self.peak[self.ComboBox1.GetSelection()][2]
            self.axes.set_ylim(float(carbonline)+float(CarbonWidth)/2.0,float(carbonline)-float(CarbonWidth)/2.0)
            self.axes.set_xlabel('Omega_C (ppm)')
            self.axes.set_ylabel('Omega_C (ppm)')

            ex=(x_min,x_max)
            ey=(carbonline,carbonline)
            self.axes.plot( ex, ey, c='r') #plot pdb network


            ###############################################################
            # Subplot 2 - the HC projection 
            
            self.axes = self.fig.add_subplot(322)
            self.axes.clear()
            #levels=[]
            #for i in range(ctr_level):
            #    levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))
            max_level=float(self.textbox_maxP.GetValue()) 
            min_level=float(self.textbox_minP.GetValue()) 
            ctr_level=int(self.textbox_lvlP.GetValue()) 
            levels=[] 
            levels.append(min_level)
            for i in range(ctr_level-1):
                levels.append(levels[i]*max_level)
            levels=numpy.array(levels)
            levels=numpy.concatenate((-1*levels[::-1],levels))
            self.axes.contour( self.Xs_yz, self.Ys_yz, self.Zs_yz,levels,cmap=colormap) #plot pdb network

            protonline=float(self.peak[self.ComboBox1.GetSelection()][1])
            carbonline=float(self.peak[self.ComboBox1.GetSelection()][2])
            label=self.peak[self.ComboBox1.GetSelection()][0]
            self.axes.text(protonline,carbonline,label)
            print(protonline,carbonline)

            self.axes.scatter(protonline,carbonline,c='r',s=50)

            #find max and min y range
            #y_max=Ys[0][0]
            #y_min=Ys[0][(len(Ys[0]))-1]
            #x_max=Xs[0][0]
            #x_min=Xs[(len(Xs))-1][0]

            ProtonWidth=1.0
            self.axes.set_xlim(float(protonline)+float(ProtonWidth)/2.0,float(protonline)-float(ProtonWidth)/2.0)
            self.axes.set_ylim(float(carbonline)+float(CarbonWidth)/2.0,float(carbonline)-float(CarbonWidth)/2.0)
            self.axes.set_xlabel('Omega_C (ppm)')
            self.axes.set_ylabel('Omega_C (ppm)')
        

        ############################################################
        #noe 1
        self.width=1.

        self.axes = self.fig.add_subplot(313)
        self.axes.clear()

        Xs,Ys,Zs=self.ReSlice2dNOE(self.ComboBox1.GetSelection(),self.peak,self.width)
        max_level=float(self.textbox1.GetValue())#set contour max level from box
        min_level=float(self.textbox0.GetValue())#set contour min level from box
        ctr_level=int(self.textbox2.GetValue())  #set the cnumber of contours from box
        y_max=Ys[0][0]
        y_min=Ys[(len(Ys))-1][0]
        x_max=Xs[0][0]
        x_min=Xs[0][(len(Xs[0]))-1]
        self.axes.set_xlim(x_max,x_min)

        #levels=[]                                #calculate the contour levels
        #for i in range(ctr_level):
        #    levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))
        levels=[] 
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels))

        #self.axes.set_title(self.ComboBox1.GetValue()+" to "+self.ComboBox2.GetValue(),fontsize=8)

        #plt.title('2D strip plot for '+self.ComboBox1.GetValue()+' at '+str(self.peak[self.ComboBox1.GetSelection()][2])+'ppm')
        self.axes.contour(Xs, Ys, Zs,levels,cmap=colormap) #plot pdb network
        #self.axes.scatter(float(self.peak[k][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),c='k',s=50)
        #self.axes.text(float(self.peak[self.ComboBox2.GetSelection()][2]),float(self.peak[self.ComboBox2.GetSelection()][1]),self.ComboBox2.GetValue())
        #self.axes.scatter(float(self.peak[self.ComboBox2.GetSelection()][2]),float(self.peak[self.ComboBox2.GetSelection()][1]),c='r',s=100)


        xl=[]
        xl.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
        xl.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
        hl=[]
        hl.append(x_min)
        hl.append(x_max)
        self.axes.plot(hl,xl,'green') #horizontal

        if(self.NOE==1):
            xl=[]
            xl.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
            xl.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
            hl=[]
            hl.append(x_min)
            hl.append(x_max)
            self.axes.plot(hl,xl,'green') #horizontal            


            self.conn_data=self.parent.tabOne.conn_data
            for i in range(len(self.conn_data)):
                if(self.conn_data[i][0]==self.ComboBox1.GetSelection()):
                    self.axes.scatter(float(self.conn_data[i][1]),float(self.peak[self.ComboBox1.GetSelection()][2]),c='b',s=100)
                    xl=[]
                    xl.append(float(self.conn_data[i][1]) )  
                    xl.append(float(self.conn_data[i][1]) ) 
                    hl=[]
                    hl.append(y_min)
                    hl.append(y_max)
                    self.axes.plot(xl,hl,'green') #horizontal            




        ##############################################################333
        #Subplot 3 - the slice
        self.width=0.1

        self.axes = self.fig.add_subplot(312)
        self.axes.clear()

        Xs,Ys,Zs=self.ReSlice2d(self.inc,self.ComboBox1.GetSelection(),self.peak,self.width)
        max_level=float(self.textbox1.GetValue())#set contour max level from box
        min_level=float(self.textbox0.GetValue())#set contour min level from box
        ctr_level=int(self.textbox2.GetValue())  #set the cnumber of contours from box
        #levels=[]                                #calculate the contour levels
        #for i in range(ctr_level):
        #    levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))
        levels=[] 
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels))

        #plt.title('2D strip plot for '+self.ComboBox1.GetValue()+' at '+str(self.peak[self.ComboBox1.GetSelection()][2])+'ppm')
        self.axes.contour(Xs, Ys, Zs,levels,cmap=colormap) #plot pdb network
        #self.axes.scatter(float(self.peak[k][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),c='k',s=50)
        self.axes.text(float(self.peak[self.ComboBox1.GetSelection()][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),self.ComboBox1.GetValue())
        self.axes.scatter(float(self.peak[self.ComboBox1.GetSelection()][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),c='r',s=100)

        y_max2=Ys[0][0]
        y_min2=Ys[(len(Ys))-1][0]
        x_max2=Xs[0][0]
        x_min2=Xs[0][(len(Xs[0]))-1]

        if(1==1): #plot the NOE line
            """
            y_min_lab=float(Ys[(len(Ys))-1][0])+0.005
            xl=[]
            xl.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
            xl.append((float(self.peak[self.ComboBox2.GetSelection()][2])))
            hl=[]
            hl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            hl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            self.axes.plot(xl,hl,'green') #horizontal
            yd=[] 
            yd.append(y_min2)
            yd.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            xd=[]
            xd.append((float(self.peak[self.ComboBox2.GetSelection()][2])))
            xd.append((float(self.peak[self.ComboBox2.GetSelection()][2])))
            self.axes.plot(xd,yd,'green') #vertical 1
            yd=[] 
            yd.append(y_min2)
            yd.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            xd=[]
            xd.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
            xd.append((float(self.peak[self.ComboBox1.GetSelection()][2])))
            self.axes.plot(xd,yd,'green') #vertical 2
            """

        if(self.cb_grid_auto.GetValue()==526):#if auto-detect has been run...
            for j in range(len(self.conn_data)):#search through the cross peaks
                if(self.index_data[self.conn_data[j][0]]==self.ComboBox1.GetValue()):#find references for the current plane...
                    for k in range(len(self.peak)):#look through the peak list...
                        if(self.peak[k][0]==self.index_data[self.conn_data[j][1]]):#to find the reference for the reciprocated cross peak...
                            self.axes.scatter(float(self.peak[k][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),c='k',s=50)

                            if(len(self.selection)==0):#if there is no selection, print all peak labels
                                self.axes.text(float(self.peak[k][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),self.index_data[self.conn_data[j][1]])
                            else:
                                for l in range(len(self.selection)):#if there is a selection, only show labels for these
                                    if(self.selection[l]==self.conn_data[j][1]):
                                        self.axes.text(float(self.peak[k][2]),float(self.peak[self.ComboBox1.GetSelection()][1]),self.index_data[self.conn_data[j][1]],color='r')
        y_max=Ys[0][0]
        y_min=Ys[(len(Ys))-1][0]
        x_max=Xs[0][0]
        x_min=Xs[0][(len(Xs[0]))-1]
    
        if(self.NOE==1):
            xl=[]
            xl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            xl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            hl=[]
            hl.append(x_min)
            hl.append(x_max)
            self.axes.plot(hl,xl,'green') #horizontal

            self.conn_data=self.parent.tabOne.conn_data
            for i in range(len(self.conn_data)):
                if(self.conn_data[i][0]==self.ComboBox1.GetSelection()):
                    self.axes.scatter(float(self.conn_data[i][1]),float(self.peak[self.ComboBox1.GetSelection()][1]),c='b',s=100)
                    xl=[]
                    xl.append(float(self.conn_data[i][1]) )  
                    xl.append(float(self.conn_data[i][1]) ) 
                    hl=[]
                    hl.append(y_min)
                    hl.append(y_max)
                    self.axes.plot(xl,hl,'green') #horizontal            



        if(self.cb_grid.IsChecked()):
            y_min_lab=float(Ys[(len(Ys))-1][0])+0.005
            xl=[]
            xl.append(x_max)
            xl.append(x_min)
            hl=[]
            hl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            hl.append((float(self.peak[self.ComboBox1.GetSelection()][1])))
            self.axes.plot(xl,hl,'black')
            yd=[]
            yd.append(y_max)
            yd.append(y_min)
            for i in range(len(self.peak)):
                xd=[]
                xd.append((float(self.peak[i][2])))
                xd.append((float(self.peak[i][2])))
                self.axes.plot(xd,yd,'black')
                #plt.text(xd[0],float(y_min_lab),self.peak[i][0],rotation=90,fontsize=8)
        if(self.ax_reset==1): #if we want to reset the axis
            self.axes.set_xlim(x_max,x_min)
            self.axes.set_ylim(y_min,y_max)
            self.xmin=x_min
            self.xmax=x_max
            self.ymin=y_min
            self.ymax=y_max
        else:#otherwise use the last saved values       
            self.axes.set_xlim(self.xmin,self.xmax)
            self.axes.set_ylim(self.ymin,self.ymax)
        #plt.xlabel('Omega_C (ppm)')
        self.ax_reset=0 #after a reset assume we have no further need to reset axes
        self.axes.set_ylabel('Omega_H (ppm)')

        self.xmin,self.xmax=self.axes.get_xlim()
        self.ymin,self.ymax=self.axes.get_ylim()




        self.canvas.draw()

    def on_cb_grid(self, event):
        self.draw_figure()

    def on_cb_grid_auto(self, event):
        self.draw_figure()



    def on_draw_button(self, event):
        self.inc=0
        self.ax_reset=1
        self.ax_reset2=1
        self.inc2=0
        self.draw_figure()


    def on_N_button(self, event):
        self.ax_reset=1
        self.inc=0
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button(self, event):
        self.ax_reset=1
        self.inc=0
        
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.selection=[]
        self.draw_figure()


    def on_Up_button(self, event):
        self.xmin,self.xmax=self.axes.get_xlim()
        self.ymin,self.ymax=self.axes.get_ylim()
        self.inc=self.inc+1
        self.selection=[]
        self.draw_figure()

    def on_Down_button(self, event):
        self.xmin,self.xmax=self.axes.get_xlim()
        self.ymin,self.ymax=self.axes.get_ylim()

        self.inc=self.inc-1
        self.selection=[]
        self.draw_figure()

    def on_NOE_button(self, event):
        bool=AssMan(self)


    def on_N_button2(self, event):
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()+1)
        self.selection=[]
        self.draw_figure()

    def on_P_button2(self, event):
        self.ax_reset2=1
        self.inc2=0
        self.ComboBox2.SetSelection(self.ComboBox2.GetSelection()-1)
        self.selection=[]
        self.draw_figure()


    #def on_Up_button2(self, event):
    #    self.xmin2,self.xmax2=plt.xlim()
    #    self.ymin2,self.ymax2=plt.ylim()
    #    self.inc2=self.inc2+1
    #    self.selection=[]
    #    self.draw_figure()

    #def on_Down_button2(self, event):
    #    self.xmin2,self.xmax2=plt.xlim()
    #    self.ymin2,self.ymax2=plt.ylim()
    #    self.inc2=self.inc2-1
    #    self.selection=[]
    #    self.draw_figure()




    def on_deselect_button(self, event):
        print('Clearing selection')
        self.selection=[]
        self.draw_figure()

    def on_delete_button(self, event):
        print('Removing selection')
        for i in range(len(self.selection)):
            for j in range(len(self.conn_data)):
                if(self.conn_data[j][0]==self.ComboBox1.GetSelection() and self.conn_data[j][1]==self.selection[i]):
                    print('Removing cross peak: ',self.index_data[self.conn_data[j][0]],'to',self.index_data[self.conn_data[j][1]])
                    self.conn_data.pop(j)
                    for k in range(len(self.conn_data)):
                        if(self.conn_data[k][0]==self.selection[i] and self.conn_data[k][1]==self.ComboBox1.GetSelection()):
                            print('Removing reciprocated cross peak: ',self.index_data[self.conn_data[k][0]],'to',self.index_data[self.conn_data[k][1]])
                            self.conn_data.pop(k)
                            break
                    break
        self.selection=[]
        self.draw_figure()




    def on_search_button(self, event):
        self.thresh=float(self.textbox0.GetValue())
        self.pick_fac=float(self.textbox_pickfac.GetValue())

        self.res=self.ComboBox1.GetSelection()
        self.conn_data=analslices1d_spec(self.res,self.selection,self.peak,self.index_data,self.conn_data,self.thresh,self.pick_fac)
        self.draw_figure()

    def on_save_button(self, event):
        self.outfile=self.textbox_savelist.GetValue()
        print()
        print('Saving ',len(self.conn_data),'entries in connectivity table to ',self.outfile)
        outy=open(self.outfile,'w')
        for i in range(len(self.conn_data)):
            outy.write('%s\t%s\n' % (self.index_data[self.conn_data[i][0]],self.index_data[self.conn_data[i][1]]))
        outy.close()






    #when search button is pressed make selection
    def on_pick(self, event):
        if(self.cb_grid_select.GetValue()==1):
            if(self.pick_cnt==0):
                self.selection=[]
                self.select=[]
                self.select.append((event.xdata,event.ydata))
                self.pick_cnt=1
                print('Click 1: ',event.xdata,event.ydata)

            else:
                self.selection=[]
                self.select.append((event.xdata,event.ydata))
                self.pick_cnt=0
                print('Click 2: ',event.xdata,event.ydata)
                for i in range(len(self.peak)):
                    if(float(self.peak[i][2])>float(self.select[0][0]) and float(self.peak[i][2])<float(self.select[1][0])):
                        self.selection.append(i)
                print('Peaks in this selection range:',ItoN1(self.selection,self.index_data))


                self.xmin,self.xmax=self.axes.get_xlim()
                self.ymin,self.ymax=self.axes.get_ylim()

                self.draw_figure()
                self.cb_grid_select.SetValue(0)




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


class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self




class AssMan(wx.App):
    def __init__(self,inherit):
        self.frame_AssManFrame=AssManFrame(None,10,'Assignment',inherit)
        self.frame_AssManFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class AssManFrame(wx.Frame):
#    title = 'AssBox'
    def __init__(self,parent, id, title,inherit):
        self._init_ctrls(parent,inherit)
        self.parent=inherit
        

    def _init_ctrls(self,prnt,parent):
        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=prnt,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title='ListBox Test ...')
        self.SetClientSize(wx.Size(900, 280))

        panel=wx.Panel(self,-1)

        corr=[]
        self.corrDict={}
        cnt=1
        inny=open('out/correlate','r')
        for line in inny.readlines():
            test=line.split()
            add=[]
            if(len(test)>0):
                corr.append(test)
                add.append(int(test[0]))
                add.append(float(test[1]))



                conf=1
                if(float(parent.tabOne.conn_data[cnt-1][2])>0.01):
                    conf+=1
                if(float(parent.tabOne.conn_data[cnt-1][2])>0.1):
                    conf+=2

                add.append(conf)
                add.append(float(parent.tabOne.conn_data[cnt-1][2]))


                for i in range(len(test)-2):
                    add.append(float(test[i+2]))

                corr.append(add)
                self.corrDict[cnt]=add
                cnt+=1
           

        self.lc=SortedListCtrl(panel,self.corrDict)

        self.lc.InsertColumn(0, 'Resonance')
        self.lc.InsertColumn(1, 'crossShift')
        self.lc.InsertColumn(2, 'Confidence')
        self.lc.InsertColumn(3, 'Kilter')
        self.lc.InsertColumn(4, 's/n')
        self.lc.InsertColumn(5, 'Shift')
        self.lc.SetColumnWidth(0, 140)
        self.lc.SetColumnWidth(1, 153)

        items=list(self.corrDict.items())
        for key,data in items:
            num_items = self.lc.GetItemCount()
            self.lc.InsertStringItem(num_items,str(data[0]))  #add assignment
            self.lc.SetStringItem(num_items, 0,str(parent.peak[data[0]][0])) #add atom        

            #self.lc.SetStringItem(num_items, 0,str(data[0]))  #add atom        
            #self.lc.SetStringItem(num_items, 1,str(data[1]))  #add atom            
            for i in range(len(data)-1):
                self.lc.SetStringItem(num_items, i+1,str(data[i+1]))  #add atom
            self.lc.SetItemData(num_items, key)




        self.Addbutton =  wx.Button(panel, 10, 'Show',(710,10))
        self.Removebutton= wx.Button(panel, 11, 'Remove',(710,60))
        self.Clearbutton = wx.Button(panel, 12, 'Clear',(710,110))
        self.Closebutton = wx.Button(panel, 13, 'Close',(710,160))
        self.Savebutton = wx.Button(panel, 14, 'Save',(710,210))

        #wx.StaticText(self, -1, 'Assman', (0,0))

        #self.pdbfile = infile
        self.textbox = wx.TextCtrl(
            panel,
            size=(150,-1),
            style=wx.TE_PROCESS_ENTER,pos=(690,240))
#        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
#        self.textbox.SetValue(' '.join(map(str, self.pdbfile)))
#self.textbox.SetValue(self.pdbfile)

        #self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnAdd,self.lc)

        self.Bind (wx.EVT_BUTTON, self.OnAdd, self.Addbutton)
        #self.Bind (wx.EVT_BUTTON, self.OnRemove, self.Removebutton)
        #self.Bind (wx.EVT_BUTTON, self.OnClear, self.Clearbutton)
        #self.Bind (wx.EVT_BUTTON, self.OnClose, self.Closebutton)
        #self.Bind (wx.EVT_BUTTON, self.OnSave, self.Savebutton)

        #self.vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)     
        hbox.Add(self.lc, 1, wx.EXPAND)

        vbox=wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.Addbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Removebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Clearbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Closebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Savebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.textbox, 0, wx.ALIGN_CENTER| wx.TOP)
        hbox.Add(vbox)
        panel.SetSizer(hbox)



        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        #hbox  = wx.BoxSizer(wx.HORIZONTAL)
        #hbox.Add(self.Addbutton, 1, wx.EXPAND)
        #hbox.Add(self.Removebutton, 1, wx.EXPAND)
        #hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        #hbox.Add(self.Closebutton, 1, wx.EXPAND)
        #hbox.Add(self.Savebutton, 1, wx.EXPAND)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        self.Centre()
        self.Show(True)





        #hbox.Add(vbox1, 1, wx.EXPAND)
        #hbox.Add(vbox2, 1, wx.EXPAND)
        #self.SetSizer(hbox)

        """
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        #TO ADD TO A LIST CONTROL
#        self.list.InsertColumn(0,"Data #1")
#        self.list.InsertColumn(1,"Data #2")
#        self.list.InsertColumn(2,"Data #3")

      # 0 will insert at the start of the list
#        pos = self.list.InsertStringItem(0,"hello")
#        # add values in the other columns on the same row
#        self.list.SetStringItem(pos,1,"world")
#        self.list.SetStringItem(pos,2,"!")


        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        self.vbox.Add(self.Addbutton,0,flag= wx.ALIGN_LEFT | wx.TOP)
#        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        self.hbox.Add(self.Pbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.Nbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox.Add(self.noisebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.AutoFitbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid_auto, 0, border=3, flag=flags)
#        self.hbox.AddSpacer(30)
#        self.hbox.Add(self.slider_label, 0, flag=flags)
#        self.hbox.Add(self.slider_width, 0, border=3, flag=flags)



        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.text1, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox0, 0, border=3, flag=flags)
        self.hbox.Add(self.text2, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox1, 0, border=3, flag=flags)
        self.hbox.Add(self.text3, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox2, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.textxn, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_xmin, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_xmax, 0, border=3, flag=flags)
        self.hbox.Add(self.textyn, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_ymin, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_ymax, 0, border=3, flag=flags)
        self.hbox.Add(self.textfac, 0, border=3, flag=flags)
        self.hbox.Add(self.textbox_fac, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        """

#        self.SetSizer(self.vbox)
#        self.vbox.Fit(self)

    def onItemSelected(self, event):
        """"""

        #currentItem = event.m_itemIndex
        #car = self.corrDict[currentItem]
        #print car

        #count = self.lc.GetItemCount()
        #self.sorted_artists = [self.list.GetItem(itemId=row, col=0).GetText() for row in xrange(count)]
        #print self.sorted_artists
        #print self.sorted_artists[currentItem]


    def AtoI(self,val):
        for i in range(len(self.parent.tabOne.peak)):
            if(val==self.parent.tabOne.peak[i][0]):
                return i

    def OnAdd(self, event):
        sele=self.lc.GetFirstSelected()
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(itemId=row, col=0).GetText() for row in range(count)][sele]
        col2 = [self.lc.GetItem(itemId=row, col=1).GetText() for row in range(count)][sele]
        self.parent.ComboBox1.SetSelection(self.AtoI(col1))
        self.parent.NOE=1

        self.parent.on_draw_button(True)


    def OnRemove(self, event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)

    def OnClose(self, event):
        self.Close()

    def OnClear(self, event):
        self.lc.DeleteAllItems()

    def OnSave(self, event):
        self.outfile=self.textbox.GetValue()
        print()
        print('Saving list to ',self.outfile)
        outy=open(self.outfile,'w')
        for itemIndex in range(self.lc.GetItemCount()):
#
#            print str(self.lc.GetItemText(self.lc.GetItem(itemIndex,0)))
#            print itemIndex,self.lc.GetItemData(self.lc.GetItem(itemIndex,0))

            item1 = self.lc.GetItem(itemIndex, 0) #GetItem(row, col)
            item2 = self.lc.GetItem(itemIndex, 1) #GetItem(row, col)
            #if(self.sum[0]==10):
            #    outy.write('%s\t%s\n' % (item1.GetText(),item2.GetText()))
            #    print item1.GetText(),item2.GetText()
            #else:
            #    outy.write('%s\t%s\n' % (item2.GetText(),item1.GetText()))
            #    print item2.GetText(),item1.GetText()
        print()
        #if(self.sum[0]!=10):
        #    outy.write('\n\n')
        outy.close()

#            print itemIndex,self.lc.GetItemText(itemIndex)

#        print self.lc.GetStringItem(0)
#        self.lc.DeleteAllItems()



 #           age = item.GetText()





#        self.button2 = wx.Button(id=wxID_FRAME1BUTTON2, label=u'Clear',
#              name='button2', parent=self, pos=wx.Point(104, 312),
#              size=wx.Size(87, 28), style=0)
#        self.button2.Bind(wx.EVT_BUTTON, self.OnButton2Button,
#              id=wxID_FRAME1BUTTON2)

#    def __init__(self, parent):


        return
    def OnListBox1Listbox(self, event):
        '''
        click list item and display the selected string in frame's title
        '''
#        selName = self.listBox1.GetStringSelection()
#        self.SetTitle(selName)
        return

    def OnButton2Button(self, event):
        '''
        click button to clear the listbox items
        '''
        self.listBox1.Clear()
