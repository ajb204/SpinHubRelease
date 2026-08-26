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
import wx,string,os,numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
import copy
import sys
#import imp
import importlib
    
import nmrglue as ng
from spinDecon.domain.peaks import connEntry
from rdkit import Chem

from matplotlib.widgets import MultiCursor

from rdkit.Chem import Draw, AllChem, rdMolTransforms, PyMol
import wx.svg as svg
import wx.lib.scrolledpanel as scrolled
import time
# try:
#     import wx.lib.wxcairo
#     import cairo
#     haveCairo = True
# except ImportError:
#     haveCairo = False
# import rsvg
# import PIL2wx
############################################################################
# Frame for 1d slices
#
def PIL2wx (image):
    width, height = image.size
    return wx.BitmapFromBuffer(width, height, image.tobytes())

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def RunFrame(uc1min,uc1max,peak,noiseVal):
    app = wx.PySimpleApp()
    frame = SliceFrame(uc1min,uc1max,peak,noiseVal)
    app.MainLoop()

def rot_ar_x(radi):
    return  numpy.array([[1, 0, 0, 0],
                      [0, numpy.cos(radi), -numpy.sin(radi), 0],
                      [0, numpy.sin(radi), numpy.cos(radi), 0],
                     [0, 0, 0, 1]], dtype=numpy.double)

def rot_ar_y(radi):
    return  numpy.array([[numpy.cos(radi), 0, numpy.sin(radi), 0],
                      [0, 1, 0, 0],
                      [-numpy.sin(radi), 0, numpy.cos(radi), 0],
                     [0, 0, 0, 1]], dtype=numpy.double)

def rot_ar_z(radi):
    return  numpy.array([[numpy.cos(radi), -numpy.sin(radi), 0, 0],
                      [numpy.sin(radi), numpy.cos(radi), 0, 0],
                      [0, 0, 1, 0],
                     [0, 0, 0, 1]], dtype=numpy.double)
tforms = {0: rot_ar_x, 1: rot_ar_y, 2: rot_ar_z}

class STDFrame(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.dim=parent.tabOne.dim
        self.tabOne=tabOne
        self.sum=(0.,2.)
        self.peak=tabOne.peak
        self.thresh=tabOne.dmax
        self.offset=0
        self.DECON=0
        self.protein_read_in = False
        self.assigning=False
        self.conf=0
        self.pressed = False
        self.dragging = False
        

        # dmax=tabOne.uc0max
        # dmin=tabOne.uc0min
        self.create_projections()
        self.create_main_panel()
        self.draw_figure()
        #self.Show()
        self.Fit()

    def create_projections(self):
        self.raw_projection = numpy.sum(self.tabOne.data, axis=0)
        self.STD_projection = numpy.sum(self.tabOne.STD, axis=0)

    def GetConn(self,infile,sym='y'):
        if(os.path.exists(infile)==0):
            print('Cannot find file:',infile)
            return -1
        inny=open(infile)
        conn_data=[]
        for line in inny.readlines():
            test=line.split()
            conn_data.append(connEntry(test,sym=sym,peak=self.peak,dim=self.dim))
        inny.close()
        return conn_data


    def drawing_box(self):
        self.vbox2Lbl = wx.StaticBox(self,-1,'Drawing:')
        self.vbox2=wx.StaticBoxSizer(self.vbox2Lbl,wx.HORIZONTAL)

        self.drawbutton = wx.Button(self, -1, "Draw!", size=(-1,22))

        self.cb_grid = wx.CheckBox(self, -1,"Peaks",style=wx.ALIGN_RIGHT)
        self.cb_calc = wx.CheckBox(self, -1,"ShowCalc",style=wx.ALIGN_RIGHT)

        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_calc)

        self.vbox2.Add(self.cb_grid, border=10, flag=self.flags)
        self.vbox2.Add(self.cb_calc, border=10, flag=self.flags)
        self.vbox2.Add(self.drawbutton, border=10, flag=self.flags)

        self.vbox2.AddSpacer(10)

    def phasing_box(self):
        self.vbox3Lbl = wx.StaticBox(self,-1,'Phasing:')
        self.vbox3=wx.StaticBoxSizer(self.vbox3Lbl,wx.HORIZONTAL)

        self.p0_slider = wx.Slider(self, value = 0, minValue = -360, maxValue = 360,size=(150,-1),
        style = wx.SL_HORIZONTAL| wx.SL_AUTOTICKS)
        self.p0_slider.Bind(wx.EVT_SLIDER, self.on_p0)
        self.p1_slider = wx.Slider(self, value = 0, minValue = -360, maxValue = 360,size=(150,-1),
        style = wx.SL_HORIZONTAL| wx.SL_AUTOTICKS)
        self.p1_slider.Bind(wx.EVT_SLIDER, self.on_p0)

        self.p0Lab = wx.StaticText(self, label="P0:")
        self.p1Lab = wx.StaticText(self, label="P1:")
        self.p0 = wx.TextCtrl(self,size=(50,22) )
        self.p1 = wx.TextCtrl(self,size=(50,22) )
        self.p0.SetValue('0')
        self.p1.SetValue('0')
        # self.Bind(wx.EVT_CHECKBOX, self.on_cb_phase, self.cb_phase)
        self.vbox3.Add(self.p0Lab, border=10, flag=self.flags)
        self.vbox3.Add(self.p0, border=10, flag=self.flags)
        self.vbox3.Add(self.p0_slider, border=10, flag=wx.ALIGN_LEFT | wx.LEFT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox3.Add(self.p1Lab, border=10, flag=self.flags)
        self.vbox3.Add(self.p1, border=10, flag=self.flags)
        self.vbox3.Add(self.p1_slider, border=10, flag=wx.ALIGN_LEFT | wx.LEFT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox3.AddSpacer(10)


    def Projections_box(self):
        self.vbox4Lbl = wx.StaticBox(self,-1,'Projections:')
        self.vbox4=wx.StaticBoxSizer(self.vbox4Lbl,wx.HORIZONTAL)

        self.cb_1D = wx.CheckBox(self, -1,"1D",style=wx.ALIGN_RIGHT)
        self.cb_STD = wx.CheckBox(self, -1,"STD",style=wx.ALIGN_RIGHT)
        self.cb_1D.SetValue(1)
        self.cb_STD.SetValue(1)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_1D)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_STD)


        self.vbox4.Add(self.cb_1D, border=10, flag=self.flags)
        self.vbox4.Add(self.cb_STD, border=10, flag=self.flags)
        self.vbox4.AddSpacer(10)

    def slices_box(self):
        self.vboxSlicesLbl = wx.StaticBox(self,-1,'Slices:')
        self.vboxSlices=wx.StaticBoxSizer(self.vboxSlicesLbl,wx.HORIZONTAL)

        self.upSliceButton = wx.Button(self, size = (30,22),label="+",style=wx.ALIGN_RIGHT)
        self.downSliceButton = wx.Button(self,size= (30,22),label="-",style=wx.ALIGN_RIGHT)
        self.current_slice_box = wx.TextCtrl(self,size=(20,22))
        self.current_slice = 0
        self.current_slice_box.SetValue(str(self.current_slice))
        
        self.Bind(wx.EVT_BUTTON, self.on_up_slice, self.upSliceButton)
        self.Bind(wx.EVT_BUTTON, self.on_down_slice, self.downSliceButton)


        self.vboxSlices.Add(self.upSliceButton, border=10, flag=self.flags)
        self.vboxSlices.Add(self.downSliceButton, border=10, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL)
        self.vboxSlices.Add(self.current_slice_box, border=10, flag=self.flags)
        self.vboxSlices.AddSpacer(10)

    def on_up_slice(self, event):
        if self.current_slice == self.tabOne.data.shape[0]-1:
            return
        self.current_slice += 1
        self.current_slice_box.SetValue(str(self.current_slice+1))
        self.update_slices()

        

    def update_slices(self):

        self.data_1d.set_ydata(self.tabOne.data[self.current_slice,:])
        self.data_STD.set_ydata(self.tabOne.STD[self.current_slice,:])
        self.axes.set_ylim(numpy.min(self.tabOne.data[:,:])*1.1, numpy.max(self.tabOne.data[:,:])*1.1 )
        self.axesSTD.set_ylim(numpy.min(self.tabOne.STD[:,:])*1.1, numpy.max(self.tabOne.STD[:,:])*1.1 )
        self.canvas.draw_idle()
        
        
    def on_down_slice(self, event):
        if self.current_slice == 0:
            return
        self.current_slice -= 1
        self.current_slice_box.SetValue(str(self.current_slice+1))
        self.cb_1D.SetValue(0)
        self.cb_STD.SetValue(0)
        self.update_slices()

    def update_projections(self):
        self.data_1d.set_ydata(self.raw_projection)
        self.data_STD.set_ydata(self.STD_projection)
        self.axes.set_ylim(min(self.raw_projection)*1.1, max(self.raw_projection)*1.1 )
        self.axesSTD.set_ylim(min(self.STD_projection)*1.1, max(self.STD_projection)*1.1 )
        

    def on_cb_grid(self, event):
        if self.cb_1D.GetValue() == 1:
            self.update_projections()
        else:
            self.update_slices()
        if self.cb_STD.GetValue() == 1:
            self.data_STD.set_visible(True)
        else:
            self.data_STD.set_visible(False)

        self.canvas.draw_idle()

    def onGetDir(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a folder",         style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            splitPath = path.split(cwd)
            try:
                textBox.SetValue('.' + splitPath[1])
            except:
                textBox.SetValue(path)
            print("You chose the following file(s):")
            print(path)
            # os.chdir(path)
            # self.dirBox.SetValue(path)
            print("CWD: ",os.getcwd())

        dlg.Destroy()

    def read_protein(self, event):

        prot_file = self.prot_file.GetValue()
        if(os.path.exists(prot_file)==0):
            dlg = wx.MessageDialog(self, message="Protein file not found...", style = wx.ICON_ERROR| wx.OK)
            if dlg.ShowModal() == wx.ID_OK:
                return
        else:
            print("Protein File exists!")
            if '.ft2' in prot_file:
                path = os.path.split(prot_file)[0]
                print(path)
                if os.path.exists(os.path.join(path, "raw.ft2")) and os.path.exists(os.path.join(path, "sat.ft2")):
                    self.STD_prot_raw_path = os.path.join(path, "raw.ft2")
                    self.STD_prot_std_path = os.path.join(path, "sat.ft2")
                    self.protein_data, self.protein_STD, self.protein_index0 = self.tabOne.read_STD_spectra(self.STD_prot_raw_path, self.STD_prot_std_path)
            # self.protein_index0, self.protein_data, self.protein_STD, self.protein_times = self.parent.tabOne.read_STD(prot_file, 'raw')
            if (self.parent.tabOne.index0 == self.protein_times).all():
                print('Protein times and mixture times are the same')
                if self.protein_read_in==False:
                    self.protein_STD_plot,= self.axesSTD.plot(self.protein_index0, self.protein_STD, color='k', label='P STD')
                    self.protein_1D_plot,= self.axes.plot(self.protein_index0, self.protein_data, color='k', label='P 1D')
                    self.prot_remove_button.Enable()
                    self.prot_choose_decon_loc.Enable()
                    self.prot_decon_loc.Enable()
                else:
                    self.protein_STD_plot.set_xdata(self.protein_index0)
                    self.protein_STD_plot.set_ydata(self.protein_STD)
                    self.protein_1D_plot.set_xdata(self.protein_index0)
                    self.protein_1D_plot.set_ydata(self.protein_data)
                self.axes.legend()
                self.axesSTD.legend()
                self.canvas.draw_idle()
                self.protein_read_in = True
            else:
                print('Difference protein and mixture times...')

    def remove_protein(self, event):
        self.protein_STD_plot.remove()
        self.protein_1D_plot.remove()
        del(self.protein_index0)
        del(self.protein_data)
        del(self.protein_STD)
        # del(self.protein_times)
        self.prot_remove_button.Disable()
        self.prot_choose_decon_loc.Disable()
        self.prot_decon_loc.Disable()
        self.protein_read_in = False
        self.canvas.draw_idle()

    def peakshape_controls(self, parent):
        self.peakshape_vBoxLbl = wx.StaticBox(parent,-1,'Peak Shape Options:')
        self.peakshape_vBox=wx.StaticBoxSizer(self.peakshape_vBoxLbl,wx.HORIZONTAL)
        self.sig1Box = wx.TextCtrl(parent, size=(50, 22))
        self.voigt1Box = wx.TextCtrl(parent, size=(50, 22))
        self.lorentz1Box = wx.TextCtrl(parent, size=(50, 22))
        self.sigLab = wx.StaticText(parent, label="Gauss:")
        self.voigtLab = wx.StaticText(parent, label="Voigt:")
        self.lorLab = wx.StaticText(parent, label="Lor:")
        self.peak_fit_button = wx.Button(parent, -1, 'Peak fit')
        self.peak_fit_button.Bind(wx.EVT_BUTTON, self.OnButtonPeakFit)
        cnt=2
        self.peakshape_vBox.Add(self.sigLab, border=10, flag=self.flags)
        self.peakshape_vBox.Add(self.sig1Box, border=10, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL)
        self.peakshape_vBox.Add(self.voigtLab, border=10, flag=self.flags)
        self.peakshape_vBox.Add(self.voigt1Box, border=10, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL)
        self.peakshape_vBox.Add(self.lorLab, border=10, flag=self.flags)
        self.peakshape_vBox.Add(self.lorentz1Box, border=10, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL)
        self.peakshape_vBox.Add(self.peak_fit_button, border=10, flag=self.flags)
        self.sig1Box.SetValue(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile, 'sig1')));cnt += 1
        self.voigt1Box.SetValue(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'voigt1')));cnt += 1
        self.lorentz1Box.SetValue(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'lor1')));cnt += 1
        self.peakshape_vBox.AddSpacer(10)

        self.vbox_uSTA.Add(self.peakshape_vBox)

    def protein_controls(self, parent):
        self.protein_vBoxLbl = wx.StaticBox(parent,-1,'Protein Options:')
        self.protein_vBox=wx.StaticBoxSizer(self.protein_vBoxLbl,wx.VERTICAL)
        self.prot_label = wx.StaticText(parent, label="Protein Spectrum?")
        self.prot_file = wx.TextCtrl(parent,size=(200,22) )
        self.prot_file_open_button = wx.Button(parent, label="...", size=(40,22))
        self.prot_read_button = wx.Button(parent, label="Read", size=(50,22))
        self.prot_remove_button = wx.Button(parent, label="Remove", size=(75,22))

        self.prot_loc_label = wx.StaticText(parent, label="Protein Location:")
        self.prot_choose_decon_loc = wx.Button(parent, label="Choose Location", size=(125,22))
        self.prot_decon_loc = wx.TextCtrl(parent,size=(75,22) )


        self.prot_horizontal = wx.BoxSizer(wx.HORIZONTAL)
        self.prot_horizontal.Add(self.prot_label, border=10, flag=self.flags)
        self.prot_horizontal.Add(self.prot_file, border=10, flag=self.flags)
        self.prot_horizontal.Add(self.prot_file_open_button, border=10, flag=self.flags)
        self.prot_horizontal.Add(self.prot_read_button, border=10, flag=self.flags)
        self.prot_horizontal.Add(self.prot_remove_button, border=10, flag=self.flags)
        self.prot_horizontal.AddSpacer(10)

        self.protein_vBox.Add(self.prot_horizontal)

        self.prot_location = wx.BoxSizer(wx.HORIZONTAL)
        # self.prot_location.AddSpacer(20)
        self.prot_location.Add(self.prot_loc_label, border=10, flag=self.flags)
        self.prot_location.Add(self.prot_decon_loc, border=10, flag=self.flags)
        self.prot_location.Add(self.prot_choose_decon_loc, border=10, flag=self.flags)
        self.prot_choose_decon_loc.Disable()
        self.prot_decon_loc.Disable()
        self.protein_vBox.Add(self.prot_location)


        self.prot_file_open_button.Bind(wx.EVT_BUTTON, lambda evt: self.onGetDir(evt, self.prot_file))
        self.prot_read_button.Bind(wx.EVT_BUTTON, self.read_protein)
        self.prot_remove_button.Bind(wx.EVT_BUTTON, self.remove_protein)
        self.prot_choose_decon_loc.Bind(wx.EVT_BUTTON, self.choose_prot_decon_location)
        self.prot_remove_button.Disable()


        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile'))!='0'):
            self.prot_file.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile')))
        else:
            self.prot_file.SetValue('None')

        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_loc'))!='0'):
            self.prot_decon_loc.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_loc'))+' ppm')
        else:
            self.prot_decon_loc.SetValue('0 ppm')

        self.vbox_uSTA.Add(self.protein_vBox)

    def deconButtons(self, parent):
        self.vBox_goLbl = wx.StaticBox(parent,-1,'Go:')
        self.vBox_go=wx.StaticBoxSizer(self.vBox_goLbl,wx.HORIZONTAL)
        self.decon_button = wx.Button(parent, -1, 'Decon')
        self.decon_button.Bind(wx.EVT_BUTTON, self.OnButtonDecon)
        self.analyse_button = wx.Button(parent, -1, 'Analyse')
        self.analyse_button.Bind(wx.EVT_BUTTON, self.OnButtonAnalyse)
        self.vBox_go.Add(self.decon_button, border=10, flag=self.flags)
        self.vBox_go.Add(self.analyse_button, border=10, flag=self.flags)
        self.vBox_go.AddSpacer(10)

        self.vbox_uSTA.Add(self.vBox_go)

    def xlim_box(self, parent):
        self.xlim_min_lab = wx.StaticText(parent, label="Min (ppm)")
        self.xlim_min_box = wx.TextCtrl(parent,size=(60,22) )
        self.xlim_max_lab = wx.StaticText(parent, label="Max (ppm)")
        self.xlim_max_box = wx.TextCtrl(parent,size=(60,22) )
        self.thresh_label = wx.StaticText(parent, label="Noise threshold: ")
        self.thresh_box = wx.TextCtrl(parent, size=(40,22))
        self.thresh_box.SetValue(self.parent.tabOne.threshBox.GetValue())
        self.xlim_min_box.Bind(wx.EVT_KILL_FOCUS, self.rezoom)
        self.xlim_max_box.Bind(wx.EVT_KILL_FOCUS, self.rezoom)

        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_xlim_min'))!='0'):
            self.xlim_min_box.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile')))
        else:
            self.xlim_min_box.SetValue(str("%.2f" % min(self.parent.tabOne.index0)))
        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_xlim_max'))!='0'):
            self.xlim_max_box.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile')))
        else:
            self.xlim_max_box.SetValue(str("%.2f" % max(self.parent.tabOne.index0)))

        self.xlim_sizerLbl = wx.StaticBox(parent,-1,'xLim Options:')
        self.xlim_sizer=wx.StaticBoxSizer(self.xlim_sizerLbl,wx.HORIZONTAL)
        self.xlim_sizer.Add(self.xlim_min_lab, border=10, flag=self.flags)
        self.xlim_sizer.Add(self.xlim_min_box, border=10, flag=self.flags)
        self.xlim_sizer.Add(self.xlim_max_lab, border=10, flag=self.flags)
        self.xlim_sizer.Add(self.xlim_max_box, border=10, flag=self.flags)
        self.xlim_sizer.Add(self.thresh_label, border=10, flag=self.flags)
        self.xlim_sizer.Add(self.thresh_box, border=10, flag=self.flags)
        self.xlim_sizer.AddSpacer(10)

        self.vbox_uSTA.Add(self.xlim_sizer)

    def OnButtonNoise(self, event):
        self.noise_choosing=False
        self.noise_move = self.canvas.mpl_connect('motion_notify_event', self.noise_choose)
        self.noise_move_select = self.canvas.mpl_connect('button_press_event', self.set_noise_choose)
        self.noise1_line = self.axes.axvline(float(self.noise_min_box.GetValue().split(' ')[0]), color='r')
        self.noise1_line_STD = self.axesSTD.axvline(float(self.noise_min_box.GetValue().split(' ')[0]), color='r')
        self.canvas.draw_idle()

        return

    def set_noise_choose(self, event):
        if self.noise_choosing == False:
            self.noise2_line = self.axes.axvline(event.xdata, color='r')
            self.noise2_line_STD = self.axesSTD.axvline(event.xdata, color='r')
            self.noise1 = event.xdata
            self.noise_choosing = True
        else:
            self.noise_min_box.SetValue(str(min(self.noise1, event.xdata)))
            self.noise_max_box.SetValue(str(max(self.noise1, event.xdata)))
            self.noise1_line.remove()
            self.noise2_line.remove()
            self.noise1_line_STD.remove()
            self.noise2_line_STD.remove()
            self.canvas.mpl_disconnect(self.noise_move)
            self.canvas.mpl_disconnect(self.noise_move_select)
            self.canvas.draw_idle()

    def noise_choose(self, event):
        if self.noise_choosing == False:
            self.noise1_line.set_xdata(event.xdata)
            self.noise1_line_STD.set_xdata(event.xdata)
        else:
            self.noise2_line.set_xdata(event.xdata)
            self.noise2_line_STD.set_xdata(event.xdata)
        self.canvas.draw_idle()


    def noise_box(self, parent):
        self.noise_min_lab = wx.StaticText(parent, label="Min (ppm)")
        self.noise_min_box = wx.TextCtrl(parent,size=(60,22) )
        self.noise_max_lab = wx.StaticText(parent, label="Max (ppm)")
        self.noise_max_box = wx.TextCtrl(parent,size=(60,22) )
        self.thresh_label = wx.StaticText(parent, label="Noise threshold: ")
        self.thresh_box = wx.TextCtrl(parent, size=(40,22))
        self.thresh_box.SetValue(self.parent.tabOne.threshBox.GetValue())
        self.noise_min_box.Bind(wx.EVT_KILL_FOCUS, self.rezoom)
        self.noise_max_box.Bind(wx.EVT_KILL_FOCUS, self.rezoom)
        self.noise_button = wx.Button(parent, -1, 'Select Noise region')
        self.noise_button.Bind(wx.EVT_BUTTON, self.OnButtonNoise)

        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_noise_min'))!='0'):
            self.noise_min_box.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile')))
        else:
            self.noise_min_box.SetValue(str("%.2f" % min(self.parent.tabOne.index0)))
        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_noise_max'))!='0'):
            self.noise_max_box.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'uSTA_prot_infile')))
        else:
            self.noise_max_box.SetValue(str("%.2f" % max(self.parent.tabOne.index0)))

        self.noise_sizerLbl = wx.StaticBox(parent,-1,'Noise Options:')
        self.noise_sizer=wx.StaticBoxSizer(self.noise_sizerLbl,wx.VERTICAL)
        self.noise_first_line=wx.BoxSizer(wx.HORIZONTAL)
        self.noise_second_line=wx.BoxSizer(wx.HORIZONTAL)
        self.noise_first_line.Add(self.noise_min_lab, border=10, flag=self.flags)
        self.noise_first_line.Add(self.noise_min_box, border=10, flag=self.flags)
        self.noise_first_line.Add(self.noise_max_lab, border=10, flag=self.flags)
        self.noise_first_line.Add(self.noise_max_box, border=10, flag=self.flags)
        self.noise_first_line.Add(self.noise_button, border=10, flag=self.flags)
        self.noise_second_line.Add(self.thresh_label, border=10, flag=self.flags)
        self.noise_second_line.Add(self.thresh_box, border=10, flag=self.flags)
        self.noise_first_line.AddSpacer(10)
        self.noise_sizer.Add(self.noise_first_line)
        self.noise_sizer.Add(self.noise_second_line)
        self.vbox_uSTA.Add(self.noise_sizer)


    def uSTA_controls(self):
        self.vbox_uSTAlbl = wx.StaticBox(self,-1,'uSTA Options:')

        self.vbox_uSTA_main=wx.StaticBoxSizer(self.vbox_uSTAlbl,wx.VERTICAL)
        self.vbox_slider = scrolled.ScrolledPanel(self, -1, name="panel1", style=wx.VSCROLL)
        self.vbox_slider.SetAutoLayout(1)
        self.vbox_slider.SetupScrolling()
        self.vbox_uSTA = wx.BoxSizer(wx.VERTICAL)
        self.vbox_uSTA.AddSpacer(10)
        # words = "A Quick Brown Insane Fox Jumped Over the Fence and Ziplined to Cover".split()
        # self.spSizer = wx.BoxSizer(wx.VERTICAL)
        # for word in words:
            # text = wx.TextCtrl(self.vbox_slider, value=word)
            # self.vbox_uSTA.Add(text)


        self.protein_controls(self.vbox_slider)
        self.peakshape_controls(self.vbox_slider)
        self.xlim_box(self.vbox_slider)
        self.noise_box(self.vbox_slider)



        self.deconButtons(self.vbox_slider)
        self.smile_box(self.vbox_slider)
        self.vbox_slider.SetSizerAndFit(self.vbox_uSTA)
        self.vbox_uSTA_main.Add(self.vbox_slider,-1, flag=wx.GROW)
        self.slider_size = self.vbox_slider.Size

    def smile_box(self, parent):
        self.smile_boxLbl = wx.StaticBox(parent,-1,'Molecule Options:')
        self.smile_box=wx.StaticBoxSizer(self.smile_boxLbl,wx.VERTICAL)
        self.smile_first_line=wx.BoxSizer(wx.HORIZONTAL)

        self.smile_text_box_lbl = wx.StaticText(parent, label='SMILE:')
        self.smile_text_box = wx.TextCtrl(parent,size=(150,22) )
        self.smile_button = wx.Button(parent, -1, 'Draw' )
        self.smile_button.Bind(wx.EVT_BUTTON, self.draw_molecule_smile)
        self.smile_first_line.Add(self.smile_text_box_lbl, border=10, flag=self.flags)
        self.smile_first_line.Add(self.smile_text_box, border=10, flag=self.flags)
        self.smile_first_line.Add(self.smile_button, border=10, flag=self.flags)
        if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'SMILE_string'))!='0'):
            self.smile_text_box.SetValue(str(self.parent.tabOne.Parse(self.parent.deconParFile,'SMILE_string')))
        else:
            self.smile_text_box.SetValue(str("None"))

        self.smile_first_line.AddSpacer(10)
        self.smile_box.Add(self.smile_first_line)
        self.vbox_uSTA.Add(self.smile_box)

    def which_atom(self, event):
        print(event.xdata, event.ydata)

    def draw_molecule_smile(self,event):
        self.save()
        self.m = Chem.MolFromSmiles(self.smile_text_box.GetValue())
        self.m = Chem.AddHs(self.m)
        AllChem.EmbedMolecule(self.m, maxAttempts=5000)
        AllChem.MMFFOptimizeMolecule(self.m)
        self.assigned_atoms=[]

        # Draw.SetComicMode()
        # self.molsvg = Draw.MolToFile(self.m,'mol.svg', size=(500,500), drawOptionsbgColor=(0,0,0))
        self.assignment_labels=[]

        self.update_assignment_labels()



        # self.redraw_svg()
        self.update_highlights([])
        # self.atoms = self.read_svg(self.svg_txt)
        self.atoms_to_environments()

        x2,y2 = self.vbox_slider.Size

        x,y = self.vbox_slider.GetVirtualSize()
        self.vbox_slider.SetVirtualSize(x,y2+y)




        self.vbox_slider.Bind(wx.EVT_PAINT, self.OnPaint)
        self.vbox_slider.Bind(wx.EVT_LEFT_DOWN, self.OnClick_molecule)
        self.vbox_slider.Bind(wx.EVT_MOTION, self.OnMoleculeDrag)
        # for x2 in range(300):
        self.vbox_slider.Scroll(0,y+y2)
        self.OnPaint(event)
        self.canvas.mpl_connect('key_press_event', self.onkeypress)
        if self.DECON==0:
            self.OnButtonAnalyse(event)

    def OnMoleculeClick(self, event):
        x, y = event.GetPosition()
        self.x_drag = x
        self.y_drag = y

    def OnMoleculeDrag(self, event):
        x, y = event.GetPosition()

        if not event.Dragging():
            event.Skip()
            return
        event.Skip()
        #obj = event.GetEventObject()
        #sx, sy = obj.GetScreenPosition()
        #self.Move(sx+x,sy+y)
        # print("Dragging position", x, y)
        rdMolTransforms.TransformConformer(self.m.GetConformer(0), tforms[1](float(self.x_drag-x)/30))
        rdMolTransforms.TransformConformer(self.m.GetConformer(0), tforms[0](float(self.y_drag-y)/30))
        self.x_drag = x
        self.y_drag = y
        self.update_highlights([])
        self.Refresh()

    def onkeypress(self, event):
        # self.rotate_molecule(event.key)
        # self.conf+=1
        # d2d = Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
        # d2d.drawOptions().addAtomIndices = True
        # print(self.m.GetConformers())
        v = PyMol.MolViewer()
        v.ShowMol(self.m)
        # d2d.DrawMoleculeWithHighlights(self.m.GetConformer(self.conf),'',dict(full_highlights),{},{},{})
        #
        # d2d.FinishDrawing()


    def atoms_to_environments(self):
        environments = {}

        for at in self.atoms.keys():
            atom = self.m.GetAtomWithIdx(at)
            if atom.GetAtomicNum() == 1:
                environments[atom.GetIdx()]=[]
                # print(atom)
                for n1 in atom.GetNeighbors():
                    n2_count = 0
                    if n1.IsInRing():
                        environments[atom.GetIdx()].append(atom.GetIdx())

                    else:
                        for n2 in n1.GetNeighbors():
                            if n2.GetAtomicNum() == 1:
                                environments[atom.GetIdx()].append(n2.GetIdx())
                                n2_count +=1
        self.environments = environments





    def redraw_svg(self, full_highlights):
        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
        # d2d.SetBackgroundColour((1.0,1.0,1.0,0.0))
        d2d.drawOptions().clearBackground=False
        d2d.DrawMoleculeWithHighlights(self.m,'',dict(full_highlights),{},{},{})

        d2d.FinishDrawing()
        svg2 = str.encode(d2d.GetDrawingText().replace('svg:',''))

        self.img = svg.SVGimage.CreateFromBytes(svg2)

    def read_svg(self, string):
        reading = False
        atoms={}
        for line in string.split('\n'):
            # print(line)
            if reading==False:
                if "<path" in line and 'class=\'atom' in line:
                    fields = line.split('\'')
                    atom_num=int(fields[1].split('-')[1])
                    coord_x = fields[3].split(' ')[1]
                    coord_y = fields[3].split(' ')[2].rstrip()
                    atoms[atom_num] = [float(coord_x), float(coord_y)]
                    # print(atoms[atom_num])
                    reading=True
            if reading==True:
                # print(line)
                if '>' in line:
                    reading=False

                # print(line)
        return(atoms)

    def update_highlights(self, highlighted_atoms):
        full_highlights={}

        for x in range(len(self.assigned_atoms_colours)):
            # print(list(self.assigned_atoms_colours[x]))
            full_highlights[self.assigned_atoms[x]] = [tuple(self.assigned_atoms_colours[x])]

        full_highlights_color=copy.deepcopy(self.assigned_atoms_colours)
        for x in highlighted_atoms:
            # full_highlights[x] = [(213./255.,242./255.,227./255.)]
            full_highlights[x] = [(113./255.,242./255.,227./255.)]

        # print(full_highlights)

        # self.molsvg = Draw.MolToFile(self.m,'mol.svg', highlightAtoms=full_highlights, size=(500,500))
        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(500,500)
        d2d.drawOptions().addAtomIndices = True
        d2d.DrawMoleculeWithHighlights(self.m,'',dict(full_highlights),{},{},{})

        d2d.FinishDrawing()
        self.svg_txt = d2d.GetDrawingText().replace('svg:','')
        # exit()
        svg2 = str.encode(d2d.GetDrawingText().replace('svg:',''))
        self.atoms = self.read_svg(self.svg_txt)

        self.img = svg.SVGimage.CreateFromBytes(svg2)


        self.redraw_svg(full_highlights)

    def OnClick_molecule(self,event):
        x, y = event.GetPosition()
        self.x_drag = x
        self.y_drag = y
        x,y= event.GetLogicalPosition(self.dc)
        # print(x,y)
        min = 1e6
        num = -1
        for key, item in self.atoms.items():
            x0 = numpy.fabs(item[0]-x)
            y0 = numpy.fabs(item[1]-y)
            dist = (x0**2+y0**2)**0.5
            if dist < min:
                min = dist
                num = int(key)

        if min < 35.:
            self.env_selected=self.environments[num]

            # self.molsvg = Draw.MolToFile(self.m,'mol.svg', highlightAtoms=self.env_selected, size=(500,500))
            # self.redraw_svg()
            self.update_highlights(self.env_selected)
            x,y=self.vbox_slider.GetViewStart()
            self.vbox_slider.Scroll(x,y)
            self.assigning=True
            self.Refresh()
        else:
            self.update_highlights([])

            # self.molsvg = Draw.MolToFile(self.m,'mol.svg', size=(500,500))
            self.env_selected=[]
            # self.redraw_svg()
            self.assigning=False

            x,y=self.vbox_slider.GetViewStart()
            self.vbox_slider.Scroll(x,y)

            self.Refresh()

    def ppm_to_assign(self, ppm):
        if ppm in self.assignments.keys():
            assign = self.assignments[ppm]
            return assign
        else:
            return False

    def rotate_molecule(self, char):
        if char == 'q':
            rdMolTransforms.TransformConformer(self.m.GetConformer(0), tforms[0](2*numpy.pi/50))
        if char == 'e':
            rdMolTransforms.TransformConformer(self.m.GetConformer(0), tforms[1](2*numpy.pi/50))
        if char == 'w':
            rdMolTransforms.TransformConformer(self.m.GetConformer(0), tforms[2](2*numpy.pi/50))
        self.update_highlights([])
        self.Refresh()

    def onclick_assign(self, event):
        print(event.mouseevent.button)
        if event.mouseevent.dblclick and event.mouseevent.inaxes==self.axes and self.assigning:
            event.artist.set_linewidth(2.0)
            ppm = event.artist.get_xdata()[0]
            yval = event.artist.get_ydata()[-1]
            assignments = self.env_selected

            print(assignments, ppm, yval)
            self.canvas.draw_idle()
            self.on_assignment(assignments, ppm, yval)

        elif event.mouseevent.button==3:
            ppm = event.artist.get_xdata()[0]
            yval = event.artist.get_ydata()[-1]
            assignments = self.ppm_to_assign(ppm)[0].split(' ,')
            if assignments != -1:
                self.remove_assignment(assignments, ppm, yval)

    def update_assignment_labels(self):
        outfile = 'out/assigned.out'
        for x in self.assignment_labels:
            x.remove()
        self.assignment_labels=[]
        self.assignments={}
        self.assigned_atoms_colours=[]
        if os.path.exists(outfile):
            assignment_file = open(outfile, 'r')
            # assigned=[]
            self.assigned_atoms=[]
            # assigned_ppms=[]
            for line in assignment_file:
                fields = line.split('\t')
                assigned = fields[0]
                assigned_ppms = float(fields[1])
                self.assignments[float(fields[1])]=[fields[0], fields[2]]
                for x in assigned.split(', '):
                    self.assigned_atoms.append(int(x))

                self.assignment_labels.append(self.axes.text(assigned_ppms, 0, assigned,fontsize=8,rotation=90, horizontalalignment='center', verticalalignment='top', color=(115./255., 186./255., 155./255.)))

            self.assigned_atoms_colours = numpy.zeros((len(self.assigned_atoms), 3))
            self.assigned_atoms_colours[:]=(115./255., 186./255., 155./255.)
            # self.molsvg = Draw.MolToFile(self.m,'mol.svg', highlightAtoms=self.assigned_atoms, highlightAtomColors=self.assigned_atoms_colours, size=(500,500))


            self.canvas.draw_idle()

    def on_assignment(self, assignments, ppm, yval):
        # self.env_selected=[]
        outfile = 'out/assigned.out'
        if os.path.exists(outfile):
            assigned=[]
            assigned_ppms=[]
            assignment_file = open(outfile, 'r+')

            for line in assignment_file:
                fields = line.split()
                assigned.append(fields[0])
                assigned_ppms.append(fields[1])
            assi=str(assignments[0])
            for x in assignments[1:]:
                assi += ', '+str(x)
            assignment_file.write('%s\t%f\t%f\n' % (assi, ppm, yval))

            assignment_file.close()

            self.update_assignment_labels()


        else:
            assignment_file = open(outfile, 'w')
            assi=str(assignments[0])
            for x in assignments[1:]:
                assi += ', '+str(x)
            assignment_file.write('%s\t%f\t%f\n' % (assi, ppm, yval))
            assignment_file.close()

            # self.on_assignment()
        self.Refresh()

    def remove_assignment(self, assignments, ppm, yval):
        outfile = 'out/assigned.out'
        print('removing assignment')
        if os.path.exists(outfile):
            assigned=[]
            assigned_ppms=[]
            assignment_file = open(outfile, 'r')
            assignment_file2 = open(outfile+'2', 'w')
            print(assignments)
            assi=str(assignments[0])
            for x in assignments[1:]:
                assi += ', '+str(x)
            for line in assignment_file:
                fields = line.split('\t')
                print(assi, ppm, fields[0], fields[1])
                if assi != fields[0] or ppm != float(fields[1]):
                    # assigned.append(fields[0])
                    # assigned_ppms.append(fields[1])
                    assignment_file2.write('%s\t%s\t%s' % (fields[0], fields[1], fields[2]))



            assignment_file.close()
            assignment_file2.close()
            os.system('mv '+outfile+'2 '+outfile)
            # self.draw_molecule_smile(None)
            # self.assignment_labels=[]
            self.update_assignment_labels()
            self.update_highlights([])
            self.Refresh()
            # self.OnPaint(None)]



        else:
            print('error: File not found')


    def OnPaint(self, event):
        self.dc = wx.PaintDC(self.vbox_slider)
        # self.dc.SetBackground(wx.Brush('gray'))

        self.vbox_slider.DoPrepareDC(self.dc)
        # dc = wx.ClientDC(self)
        x,y = self.dc.GetDeviceOrigin().Get()

        self.dc.SetDeviceOrigin(x+0, y+self.slider_size.y)


        # dc.Clear()
        # print("Painting")

        dcdim = 500 #min(self.vbox_uSTA.Size.width, self.vbox_uSTA.Size.height)

        imgdim = min(self.img.width, self.img.height)
        scale = dcdim / imgdim
        width = int(self.img.width * scale)
        height = int(self.img.height * scale)

        ctx = wx.GraphicsContext.Create(self.dc)

        # bmp = self.img.ConvertToBitmap(scale=scale, width=width, height=height)
        # dc.DrawBitmap(bmp, 100, 500)
        # ctx = wx.GraphicsContext.Create(dc)
        self.img.RenderToGC(ctx, scale)
        # self.vbox_uSTA.Add(ctx)

    def OnButtonDecon(self,event):
        self.save()


        self.dim=1
        self.ncpus=4 #self.coreBox.GetValue()
        # self.thresh=0.0001 #float(self.threshBox.GetValue())
        self.sig1=(self.sig1Box.GetValue())
        self.voigt1=(self.voigt1Box.GetValue())
        self.lor1=(self.lorentz1Box.GetValue())
        self.fac= 1.4 #(self.facBox.GetValue())
        self.squash=1 #(self.squashBox.GetValue())
        self.maxiter=10000 #(self.maxiterBox.GetValue())

        print("Running binary")
        symmy=0    #symmetric mode?
        dec3d=0 #deconvolve in high dimensions?


        specstr=self.parent.tabOne.deconBin
        specstr+=' '+str(self.ncpus)
        specstr2='ncpus '+str(self.ncpus)
        specstr+=' '+"0"
        specstr2+='\tpeakList '+"0"
        specstr+=' '+str(self.dim)
        specstr2+='\tdim '+str(self.dim)
        specstr+=' '+str(1)
        specstr2+='\tuSTA '+str(1)
        specstr+=' '+self.parent.tabOne.STD_raw_path
        specstr2+='\tinfile '+self.parent.tabOne.STD_raw_path
        specstr+=' '+str(self.thresh)
        specstr2+='\tdmax '+str(float(self.thresh_box.GetValue())*self.tabOne.dmax)
        specstr+=' '+str(self.sig1)
        specstr2+='\tsig1 '+str(self.sig1)
        specstr+=' '+str(self.fac)
        specstr2+='\tfac '+str(self.fac)
        specstr+=' '+str(7)
        specstr2+='\trand '+str(7)

        specstr += ' ' + str(self.voigt1)
        specstr2+='\tvoigt1 '+str(self.voigt1)
        specstr += ' ' + str(self.lor1)
        specstr2+='\tlor1 '+str(self.lor1)
        if self.prot_file.GetValue() != 'None':
            path = os.path.split(self.prot_file.GetValue())[0]
            print(path)
            if os.path.exists(os.path.join(path, "raw.ft2")) and os.path.exists(os.path.join(path, "sat.ft2")):
                self.STD_prot_raw_path = os.path.join(path, "raw.ft2")
                self.STD_prot_std_path = os.path.join(path, "sat.ft2")
            if os.path.exists(self.STD_prot_raw_path):
                specstr += ' ' + self.STD_prot_raw_path
                specstr += " " + str(self.prot_decon_loc.GetValue()).split(' ')[0]
                specstr += ' 30' #windowB
                specstr += ' ' + str("False")
                specstr2+='\tbaseFile '+self.STD_prot_raw_path
                specstr += ' 30' #noiseval of std - deprecated
            else:
                dlg = wx.MessageDialog(self, message="Protein file not found...", style = wx.ICON_ERROR| wx.OK)
                if dlg.ShowModal() == wx.ID_OK:
                    return
        else:
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr2+='\tbaseFile False'


        print(specstr)
        print(specstr2)

        fields=specstr2.split('\t')
        print(fields)
        service = getattr(self.parent, "decon_service", None)
        if service is None:
            from spinDecon.project.decon_service import DeconService
            service = DeconService()

        service.write_init_lines(fields)

        sys.stdout.flush()

        service.launch(specstr)


        print('Done Deconning the raw 1D spectrum!')
        print()
        os.system('cp out/correlate.3 out/raw.3')
        os.system('cp out/correlate.base.3 out/raw_prot.3')


        specstr=self.parent.tabOne.deconBin
        specstr+=' '+str(self.ncpus)
        specstr2='ncpus '+str(self.ncpus)
        specstr+=' '+"0"
        specstr2+='\tpeakList '+"0"
        specstr+=' '+str(self.dim)
        specstr2+='\tdim '+str(self.dim)
        specstr+=' '+str(1)
        specstr2+='\tuSTA '+str(1)
        specstr+=' '+self.parent.tabOne.STD_std_path
        specstr2+='\tinfile '+self.parent.tabOne.STD_std_path
        self.threshSTD=float(numpy.max(self.tabOne.STD)*float(self.thresh_box.GetValue())*0.01)
        specstr+=' '+str(self.threshSTD)
        specstr2+='\tdmax '+str(float(self.thresh_box.GetValue())*self.tabOne.dmax)
        specstr+=' '+str(self.sig1)
        specstr2+='\tsig1 '+str(self.sig1)
        specstr+=' '+str(self.fac)
        specstr2+='\tfac '+str(self.fac)
        specstr+=' '+str(7)
        specstr2+='\trand '+str(7)
        #specstr+=' '+str(self.maxiter)
        specstr += ' ' + str(self.voigt1)
        specstr2+='\tvoigt1 '+str(self.voigt1)
        specstr += ' ' + str(self.lor1)
        specstr2+='\tlor1 '+str(self.lor1)
        if self.prot_file.GetValue() != 'None':
            if os.path.exists(self.STD_prot_std_path):
                specstr += ' ' + self.STD_prot_std_path
                specstr += " " + str(self.prot_decon_loc.GetValue()).split(' ')[0]
                specstr += ' 30' #windowB
                specstr += ' ' + str("out/raw.3")
                specstr2+='\tbaseFile '+self.STD_prot_std_path
                specstr += ' 30' #noiseval of std - deprecated
            else:
                dlg = wx.MessageDialog(self, message="Protein file not found...", style = wx.ICON_ERROR| wx.OK)
                if dlg.ShowModal() == wx.ID_OK:
                    return
        else:
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr += ' ' + str("False")
            specstr2+='\tbaseFile False'


        print(specstr)
        print(specstr2)

        fields=specstr2.split('\t')
        print(fields)
        service = getattr(self.parent, "decon_service", None)
        if service is None:
            from spinDecon.project.decon_service import DeconService
            service = DeconService()

        service.write_init_lines(fields)

        sys.stdout.flush()

        service.launch(specstr)
        print('Done Deconning the std 1D spectrum!')
        print()
        os.system('cp out/correlate.3 out/std.3')
        os.system('cp out/correlate.base.3 out/std_prot.3')
        self.OnButtonAnalyse(event)

    def calculate_transfer_efficiencies(self):
        self.TEs = {}
        for cn in self.raw_conn_data:
            for cn2 in self.std_conn_data:
                if cn.f1 == cn2.f1:
                    self.TEs[cn.f1] = float(cn2.s1/cn.s1)*100.0

    def OnButtonAnalyse(self, event):
        if os.path.exists('out/raw_prot.3') and os.path.exists('out/std_prot.3') and os.path.exists(self.parent.tabOne.STD_raw_path+'.decon') and os.path.exists(self.parent.tabOne.STD_std_path+'.decon'):
            for line in open('out/raw_prot.3'):
                self.raw_prot_decon_factor = float(line.split()[3])
            for line in open('out/std_prot.3'):
                self.std_prot_decon_factor = float(line.split()[3])

            raw_decon_path = self.parent.tabOne.STD_raw_path+'.decon'
            if getattr(self, "store", None) is not None and self.store.metadata.get("STD_raw_decon_path") == raw_decon_path:
                self.raw_dic_dec = self.store.metadata.get("STD_raw_dic_dec")
                self.raw_data_dec = self.store.metadata.get("STD_raw_data_dec")
            else:
                self.raw_dic_dec,self.raw_data_dec=ng.pipe.read(raw_decon_path)
                if getattr(self, "store", None) is not None:
                    self.store.metadata["STD_raw_decon_path"] = raw_decon_path
                    self.store.metadata["STD_raw_dic_dec"] = self.raw_dic_dec
                    self.store.metadata["STD_raw_data_dec"] = self.raw_data_dec
            self.conn_data=[]
            self.raw_conn_data = self.GetConn('out/raw.3',sym='n')
            if self.protein_read_in == False:
                self.read_protein(None)

            self.protein_1D_plot.set_ydata(self.protein_data*self.raw_prot_decon_factor)
            self.protein_1D_plot.set_label("P 1D x%.2f" %self.raw_prot_decon_factor)
            self.protein_STD_plot.set_label("P STD x%.2f" %self.std_prot_decon_factor)



            if self.DECON == 0:
                self.peaks_plotted = []
                self.raw_decon_line, = self.axes.plot(self.tabOne.index0, self.raw_data_dec, color='red', lw=0.5, label='Resimulated')
                for cn in self.raw_conn_data:
                        x=(cn.f1,cn.f1)
                        y=(0,cn.s1)
                        self.peaks_plotted.append(self.axes.plot(x,y,'orange', lw=0.5,picker=True, pickradius=5)[0])




            else:
                self.raw_decon_line.set_ydata(self.raw_data_dec)
                for peak in self.peaks_plotted:
                    peak.remove()
                self.peaks_plotted=[]
                for cn in self.raw_conn_data:
                        x=(cn.f1,cn.f1)
                        y=(0,cn.s1)
                        self.peaks_plotted.append(self.axes.plot(x,y,'orange', lw=0.5,picker=True, pickradius=5)[0])

            self.canvas.mpl_connect('pick_event', self.onclick_assign)




            std_decon_path = self.parent.tabOne.STD_std_path+'.decon'
            if getattr(self, "store", None) is not None and self.store.metadata.get("STD_std_decon_path") == std_decon_path:
                self.raw_dic_dec = self.store.metadata.get("STD_std_dic_dec")
                self.raw_data_dec = self.store.metadata.get("STD_std_data_dec")
            else:
                self.raw_dic_dec,self.raw_data_dec=ng.pipe.read(std_decon_path)
                if getattr(self, "store", None) is not None:
                    self.store.metadata["STD_std_decon_path"] = std_decon_path
                    self.store.metadata["STD_std_dic_dec"] = self.raw_dic_dec
                    self.store.metadata["STD_std_data_dec"] = self.raw_data_dec
            self.conn_data_STD=[]
            self.std_conn_data = self.GetConn('out/std.3',sym='n')

            self.calculate_transfer_efficiencies()

            if self.DECON == 0:
                self.peaks_plotted_STD = []
                self.raw_decon_line_STD, = self.axesSTD.plot(self.tabOne.index0, self.raw_data_dec, color='red', lw=0.5, label='Resimulated')
                # self.
                for cn in self.std_conn_data:
                        x=(cn.f1,cn.f1)
                        y=(0,cn.s1)
                        self.peaks_plotted_STD.append(self.axesSTD.plot(x,y,'orange', lw=0.5)[0])
                        self.axesSTD.text(cn.f1, 0, "%.2f" % self.TEs[cn.f1],fontsize=8,rotation=90, horizontalalignment='center', verticalalignment='top')
            else:
                self.raw_decon_line_STD.set_ydata(self.raw_data_dec)
                for peak in self.peaks_plotted_STD:
                    peak.remove()
                self.peaks_plotted_STD=[]
                for cn in self.std_conn_data:
                        x=(cn.f1,cn.f1)
                        y=(0,cn.s1)
                        self.peaks_plotted_STD.append(self.axesSTD.plot(x,y,'orange', lw=0.5)[0], horizontalalignment='center', verticalalignment='top')
            self.axes.legend()
            self.axesSTD.legend()

            self.canvas.draw_idle()
            self.DECON=1
        else:
            dlg = wx.MessageDialog(self, message="Rerun decon please - files not found...", style = wx.ICON_ERROR| wx.OK)
            if dlg.ShowModal() == wx.ID_OK:
                return

    def save(self):
        write={}
        write['sig1']=self.sig1Box.GetValue()
        write['voigt1'] = self.voigt1Box.GetValue()
        write['lor1'] = self.lorentz1Box.GetValue()
        write['uSTA_prot_infile'] = self.prot_file.GetValue()
        write['uSTA_prot_loc'] = self.prot_decon_loc.GetValue().split(' ')[0]
        write['thresh'] = self.thresh_box.GetValue()
        write['SMILE_string'] = self.smile_text_box.GetValue()

        dec=[]
        if(os.path.exists(self.parent.tabOne.deconParFile)):
            inny=open(self.parent.tabOne.deconParFile)
            for line in inny.readlines():
               dec.append(line)
            inny.close()

        outy=open(os.path.join(self.parent.tabOne.dirBox.GetValue(),self.parent.tabOne.deconParFile),'w')
        print('saving:',os.path.join(self.parent.tabOne.dirBox.GetValue(),self.parent.tabOne.deconParFile))
        for de in dec:
            test=de.split()
            if(len(test)>0):

                if(test[0] in list(write.keys())):
                    outy.write('%s = %s ' % (test[0],write[test[0]]))
                    for j in range(len(test)-3):
                        outy.write(' %s' % test[j+3])
                    outy.write('\n')
                    del write[test[0]]
                else:
                    outy.write(de)
            else:
                outy.write(de)

        for key,vals in list(write.items()):
            try:
                outy.write('%s = %s\n' % (key,vals[0]))
            except:
                print('problem with ',key)

        outy.close()



    def rezoom(self, event):
        print(self.axes.get_xlim())
        xmin = float(self.xlim_min_box.GetValue())
        xmax = float(self.xlim_max_box.GetValue())
        self.axes.set_xlim(xmax, xmin)
        self.canvas.draw_idle()
        # self.axis.set_xlim(xmax, xmin)

    def choose_prot_decon_location(self, event):
        self.prot_move = self.canvas.mpl_connect('motion_notify_event', self.protein_location)
        self.prot_move_select = self.canvas.mpl_connect('button_press_event', self.set_protein_location)
        self.protein_line = self.axes.axvline(float(self.prot_decon_loc.GetValue().split(' ')[0]))
        self.protein_line_STD = self.axesSTD.axvline(float(self.prot_decon_loc.GetValue().split(' ')[0]))
        self.canvas.draw_idle()

    def protein_location(self, event):
        if event.inaxes:
            self.protein_line.set_xdata(event.xdata)
            self.protein_line_STD.set_xdata(event.xdata)
            # self.protein_line.set_visible(True)
            self.canvas.draw_idle()

    def set_protein_location(self, event):
        if event.inaxes:
            self.prot_decon_loc.SetValue("%.2f ppm" % event.xdata)
            self.canvas.mpl_disconnect(self.prot_move)
            self.canvas.mpl_disconnect(self.prot_move_select)
            # self.protein_line.remove()
            # self.protein_line_STD.remove()
            # del(self.protein_line)
            # del(self.protein_line_STD)
            self.canvas.draw()


    def OnButtonPeakFit(self, event):
        from spinDecon.gui.workspaces import peak_fit
        peak_fit = imp.reload(peak_fit)
        fit_window = peak_fit.peakFitFrame(self.parent.tabOne)
        fit_window.Bind(wx.EVT_CLOSE, self.update_sigs)


    def update_sigs(self, event):
        print('updating')
        self.sig1Box.SetValue(self.parent.tabOne.sig1Box.GetValue())
        self.lorentz1Box.SetValue(self.parent.tabOne.lorentz1Box.GetValue())
        self.voigt1Box.SetValue(self.parent.tabOne.voigt1Box.GetValue())

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.axes = self.fig.add_subplot(211)
        self.axesSTD = self.fig.add_subplot(212, sharex=self.axes)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['top'].set_visible(False)
        self.axesSTD.spines['right'].set_visible(False)
        self.axesSTD.spines['top'].set_visible(False)
        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL

        # self.drawbutton = wx.Button(self, -1, "Draw!")
        # self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.toolbar = NavigationToolbar(self.canvas)
        self.hbox_top = wx.BoxSizer(wx.HORIZONTAL)
        self.uSTA_controls()
        self.hbox_top.AddSpacer(20)
        self.hbox_top.Add(self.vbox_uSTA_main, 2, wx.LEFT | wx.TOP | wx.GROW)
        self.hbox_top.AddSpacer(20)
        self.hbox_top.Add(self.canvas, 2, wx.LEFT | wx.TOP | wx.GROW)
        self.hbox_top.AddSpacer(20)



        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.AddSpacer(20)
        self.vbox.Add(self.hbox_top, 10, wx.GROW)
        self.vbox.AddSpacer(20)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(20)

        self.drawing_box()
        self.Projections_box()
        self.slices_box()
        self.phasing_box()

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vbox2)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vbox4)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vboxSlices)
        self.hbox.AddSpacer(5)
        self.hbox.Add(self.vbox3)




        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)


    def draw_figure(self):
        #try:
        self.hilb_data = 'a'
        self.draw_figureGO()
        #except:
        #    pass

    def hilb_T(self):
        hilb_data = ng.process.proc_base.ht(self.raw_projection, self.raw_projection.shape[0])
        return hilb_data
        

    def phaser(self, hilb_data, p0, p1):
        return ng.process.proc_base.ps(hilb_data, p0=p0, p1 =p1)

    def on_p0(self, event):
        self.phasing=True
        self.p0.SetValue(str(self.p0_slider.GetValue()))
        self.p1.SetValue(str(self.p1_slider.GetValue()))
        print('p0')
        # self.draw_figure()
        
        p0=float(self.p0.GetValue()) #/180.*numpy.pi
        p1=float(self.p1.GetValue()) #/180.*numpy.pi
        print(p0)

        if self.hilb_data == 'a':
            self.hilb_data = self.hilb_T()
        

        self.data_1d.set_ydata(self.phaser(self.hilb_data, p0, p1))
        self.canvas.draw_idle()



    def onFocus(self, event):
        print("uSTA has focus!")

    def on_press(self, event):
        self.coord1 = event.xdata, event.ydata
        self.pressed = True

    def on_release(self, event):
        if self.dragging == True:
            self.pressed = False
            self.dragging = False
            
            print('released')
            self.axes.set_xlim(max(self.coord2[0],self.coord1[0]),min(self.coord2[0],self.coord1[0]))
            self.coord1 = None
            self.coord2 = None
            self.canvas.draw_idle()

    def on_motion(self, event):
        if self.pressed:
            self.coord2 = event.xdata,event.ydata
            self.dragging = True
            # self.draw_figure()

    def draw_figureGO(self):
        """ Redraws the figure
        """
        self.axes.clear()
        print('drawing')
        #sele1=self.ComboBox1.GetSelection()
        self.thresh=float(self.tabOne.dmax*float(self.thresh_box.GetValue()))
        # self.threshSTD=float(numpy.max(self.tabOne.STD)*float(self.thresh_box.GetValue()))

        print('drawing')

        xs=self.tabOne.index0
        y2s=numpy.zeros_like(xs)
        y2s.fill(self.thresh)

        self.axes.set_xlabel(self.parent.tabOne.labb[0],fontsize=8)
        print('about to draw data')
        self.cursor = MultiCursor(self.canvas, (self.axes, self.axesSTD), color='r', lw=1)
       
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)


        self.data_1d, = self.axes.plot(xs,self.raw_projection,'darkblue',label='P+L 1D')
        self.data_STD, = self.axesSTD.plot(xs,self.STD_projection,'lightblue',label='P+L STD')

        self.axes.plot(xs,y2s,'g',label='Thresh', linewidth=0.5)
        if(self.tabOne.DECON==1):
            print(self.tabOne.DECON)
            # if(self.cb_calc.GetValue()==1):
            #     self.axes.plot(xs,self.datadec,'k',label='data', linewidth=0.5)
            # if self.cb_grid.GetValue()==1:
            #     for cn in self.parent.tabOne.conn_data:
            #             x=(cn.f1,cn.f1)
            #             y=(0,cn.s1)
            #             self.peaks_plot.append(self.axes.plot(x,y,'k', lw=0.5)[0])
                    #self.axes.text(cn.f3,-float(self.offset),cn.tag2,fontsize=8,rotation=90)

        self.xmin,self.xmax=self.axes.get_xlim()
        self.axes.set_xlim(self.xmax, self.xmin)
        self.ymin,self.ymax=self.axes.get_ylim()
        self.offset=-1*self.ymin/2
        self.fig.tight_layout()
        self.axes.legend()
        self.axesSTD.legend()
        self.canvas.draw()



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
        #box_points = event.artist.get_bbox().get_points()

        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        dlg = wx.MessageDialog(
            self,
            msg,
            "Click!",
            wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def on_scroll(self, event):
        # print('scrolling')
        # step = numpy.abs(event.step)
        # step = event.step
        step = numpy.sign(event.step)*min(19,numpy.abs(event.step))
        print(step)
        if event.inaxes==self.axes:
            self.ymin,self.ymax=self.axes.get_ylim()
            self.axes.set_ylim(self.ymin+(self.ymin*0.05*step), self.ymax+(self.ymax*0.05*step))
        elif event.inaxes==self.axesSTD:
            self.ymin,self.ymax=self.axesSTD.get_ylim()
            self.axesSTD.set_ylim(self.ymin+(self.ymin*0.05*step), self.ymax+(self.ymax*0.05*step))

        self.canvas.draw_idle()

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
