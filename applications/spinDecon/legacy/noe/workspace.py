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
import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
import nmrglue as ng
from matplotlib.figure import Figure

##########################################################################
# 2D plotting of NMR slices
#
class NOEFrame(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent,tabOne):
    #def __init__(self,uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile):
        wx.Panel.__init__(self, parent=parent)

        self.tabOne=tabOne
        self.state = getattr(tabOne, "state", getattr(parent, "state", None))
        
        #wx.Frame.__init__(self, None, -1, self.title)
        #copy in the previous variables
        #self.index_data=self.index(tabOne.peak)
        #self.thresh=(tabOne.noiseVal)
        #self.offset=copy.deepcopy(tabOne.offset)
        #self.offset=0.0
        #self.peak=(tabOne.peak)
        #self.conn_data=tabOne.conn_data
        #self.spectrumfile=tabOne.spectrumfile
        #get 2d strips from 3d data
        #self.GetSlice2d(self.spectrumfile)


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
        #self.canvas.mpl_connect('button_press_event', self.on_pick)

        #min max and lvls for slices
        self.text_slice=wx.StaticText(self, -1, 'Slices:',size=(80,-1))
        self.text1=wx.StaticText(self, -1, 'Min:')
        """
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

        self.textbox_minP.SetValue(str(self.thresh*10/2))
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

        self.ComboBox2=wx.ComboBox(self, -1, pos=(620, 180), size=(80, -1), choices=listy, style=wx.CB_READONLY)
        self.ComboBox2.SetSelection(0)





        self.Nbutton = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button, self.Nbutton)

        self.Pbutton = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button, self.Pbutton)

        self.Nbutton2 = wx.Button(self, -1,"Next")
        self.Bind(wx.EVT_BUTTON, self.on_N_button2, self.Nbutton2)

        self.Pbutton2 = wx.Button(self, -1,"Previous")
        self.Bind(wx.EVT_BUTTON, self.on_P_button2, self.Pbutton2)


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
        """

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)

        
        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.ComboBox1, 0, border=3, flag=flags)
        #self.hbox.Add(self.Pbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.Nbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.Upbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.Downbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.NOEbutton, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.ComboBox2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Pbutton2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Nbutton2, 0, border=3, flag=flags)
        #self.hbox.Add(self.Upbutton2, 0, border=3, flag=flags)
        ##self.hbox.Add(self.Downbutton2, 0, border=3, flag=flags)

        #self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        #self.hbox.Add(self.cb_grid_auto, 0, border=3, flag=flags)
        #self.hbox.Add(self.text_pickfac, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox_pickfac, 0, border=3, flag=flags)


#        self.hbox.Add(self.cb_grid2, 0, border=3, flag=flags)
#        self.hbox.AddSpacer(30)
#        self.hbox.Add(self.slider_label, 0, flag=flags)
#        self.hbox.Add(self.slider_width, 0, border=3, flag=flags)
        #self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.text_proj,0, border=3, flag=flags)
        #self.hbox.Add(self.text1p, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox_minP, 0, border=3, flag=flags)
        #self.hbox.Add(self.text2p, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox_maxP, 0, border=3, flag=flags)
        #self.hbox.Add(self.text3p, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox_lvlP, 0, border=3, flag=flags)
        #self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)



        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.text_slice, 0, border=3, flag=flags)
        #self.hbox.Add(self.text1, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox0, 0, border=3, flag=flags)
        #self.hbox.Add(self.text2, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox1, 0, border=3, flag=flags)
        #self.hbox.Add(self.text3, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox2, 0, border=3, flag=flags)
        #self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)


        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox.Add(self.cb_grid_select, 0, border=3, flag=flags)
        #self.hbox.Add(self.searchbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.deselectbutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.deletebutton, 0, border=3, flag=flags)
        #self.hbox.Add(self.text_savelist, 0, border=3, flag=flags)
        #self.hbox.Add(self.textbox_savelist, 0, border=3, flag=flags)
        #self.hbox.Add(self.savebutton, 0, border=3, flag=flags)
        #self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
       

        #self.panel.SetSizer(self.vbox)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        corr=[]
        inny=open('out/correlate','r')
        for line in inny.readlines():
            test=line.split()
            add=[]
            for i in range(len(test)):
                add.append(float(test[i]))
            if(len(test)>0):
                corr.append(numpy.array(add))
        self.corr=numpy.array(corr)

    def draw_figure(self):
        """ Redraws the figure
        """
        colormap=cm.seismic


        CarbonWidth=2.0
        
        self.axes = self.fig.add_subplot(211)
        self.axes.clear()
        
        maxy=numpy.max(numpy.fabs(self.corr[:,5]))
        miny=numpy.min(numpy.fabs(self.corr[:,5]))
        tot=30
        ref=numpy.linspace(0,tot-1.,tot)
        bins=miny*10**(numpy.log10(maxy/miny)*(ref/(tot-1.)))

        bins=miny+(maxy-miny)*ref/(tot-1.)

        n, bins, patches = self.axes.hist(numpy.fabs(self.corr[:,5]), bins=bins, facecolor='green')

        #self.axes.set_xscale('log')

        #max_level=float(self.textbox1.GetValue())#set contour max level from box
        #min_level=float(self.textbox0.GetValue())#set contour min level from box
        #ctr_level=int(self.textbox2.GetValue())  #set the cnumber of contours from box                
        #levels=[]
        #for i in range(ctr_level):
        #    levels.append(min_level+float(i)*(max_level-min_level)/(ctr_level-1))
        #levels=[] 
        #levels.append(min_level)
        #for i in range(ctr_level-1):
        #    levels.append(levels[i]*max_level)
        #levels=numpy.array(levels)
        #levels=numpy.concatenate((-1*levels[::-1],levels))
        #self.axes.contour( self.Ys_xy, self.Xs_xy, self.Zs_xy,levels,cmap=colormap) #plot pdb network
        #find max and min y range
        #y_max=Ys[0][0]
        #y_min=Ys[0][(len(Ys[0]))-1]
        #x_max=self.Xs_xy[0][0]
        #x_min=self.Xs_xy[(len(self.Xs_xy))-1][0]

        #self.axes.set_xlim(x_max,x_min)

        #carbonline=self.peak[self.ComboBox1.GetSelection()][2]
        #self.axes.set_ylim(float(carbonline)+float(CarbonWidth)/2.0,float(carbonline)-float(CarbonWidth)/2.0)
        #self.axes.set_xlabel('Omega_C (ppm)')
        #self.axes.set_ylabel('Omega_C (ppm)')

        #ex=(x_min,x_max)
        #ey=(carbonline,carbonline)
        #self.axes.plot( ex, ey, c='r') #plot pdb network
