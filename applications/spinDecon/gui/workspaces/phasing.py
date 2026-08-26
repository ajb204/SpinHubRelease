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
###################################################################
# Phase nmr spectrum
###################################################################

import numpy,sys,os,wx,platform
import nmrglue as ng
from spinDecon.gui.dialogs.processing.process import path_escape
from spinDecon.processing.nmrpipe_scripts import MakeProj4D, MakeProj3D
# import matplotlib            #import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm
from matplotlib.widgets import RectangleSelector
import matplotlib.cm as cm
import matplotlib.colors as colors
# from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# from matplotlib.widgets import Slider
from matplotlib.lines import Line2D


import copy
from spinDecon.gui.context import context_for, project_for
from spinDecon.analysis.phasing_service import PhasingService
#import imp
import importlib

class Phasing(wx.Panel):

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.parent=parent
        self.app_context = context_for(tabOne, parent)
        self.phasing_service = (getattr(self.app_context, "phasing", None)
                                if self.app_context is not None else None) or PhasingService(tabOne)
        self.topology = self.phasing_service.topology
        self.spectral_dim_count = self.topology.spectral_dim_count
        self.physical_dim_count = self.topology.physical_dim_count
        self.dim = self.spectral_dim_count  # compatibility alias: spectral only
        self.labb = self.phasing_service.labels
        self.dirbox = self.phasing_service.working_directory
        self.state = project_for(tabOne, parent)
        self.sum=(0.,2.)
        self.peak = self.phasing_service.peaks
        self.pseudo_spectrum = self.phasing_service.pseudo_spectrum
        self.offset=0
        self.peaks_drawn=False
        self.peak_list = []
        self.peak_list_names = []
        self.plotting_atom_results = False

        dmin, dmax = self.phasing_service.axis_limits(0)

        self.p0_0=0.0
        self.p1_0=0.0
        self.p0_1=0.0
        self.p1_1=0.0
        self.p0_2=0.0
        self.p1_2=0.0



        self.create_main_panel()
        

        self.SetSizerAndFit(self.main_sizer)



    def create_main_panel(self):
        # create a panel existing of: 2/3 frames for each dimension of the 2D/3D spectrum showing the relevant projections allowing
        # the correct phasing to be accurately determined which can be subsequently inputted into the processing frame
        
        # Note: will need to include a slider for p0 and p1 for each dimension


        self.fig = Figure()
        self.fig.clear()
        self.canvas = FigCanvas(self, -1, self.fig)

       

        self.toolbar = NavigationToolbar(self.canvas)

        self.main_sizer=wx.BoxSizer(wx.VERTICAL)

        

        self.figures_sizer = wx.BoxSizer(wx.VERTICAL)
        self.FilePathToRead()
        self.desired1DSlices()

        self.axes1 = self.fig.add_subplot(121)
        self.axes2 = self.fig.add_subplot(122)
        
        
        self.figures_sizer.Add(self.fileSizer)
        self.figures_sizer.Add(self.desired_slices)


        self.buttonRead = wx.Button(self, label="ReadData",size=(100,44))
        self.buttonRead.Bind(wx.EVT_BUTTON, self.OnButtonRead)

        self.figures_sizer.Add(self.buttonRead)


        self.figures_sizer.Add(self.canvas, 1, wx.CENTER)
        self.figures_sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.main_sizer.Add(self.figures_sizer)

        self.SetSizerAndFit(self.main_sizer)

       
    
    def desired1DSlices(self):

        self.desired_slices_label = wx.StaticBox(self, label='Desired ppms for 1D slices')
        self.desired_slices=wx.StaticBoxSizer(self.desired_slices_label, wx.HORIZONTAL)
        self.dim1_slice_label = wx.StaticText(self, label="dim1 shift (ppm):")
        self.dim1_slice_box = wx.TextCtrl(self, size=(100,22), style=wx.TE_PROCESS_ENTER)
        self.dim2_slice_label = wx.StaticText(self, label="dim2 shift (ppm):")
        self.dim2_slice_box = wx.TextCtrl(self, size=(100,22), style=wx.TE_PROCESS_ENTER)

        self.desired_slices.Add(self.dim1_slice_label)
        self.desired_slices.Add(self.dim1_slice_box)
        self.desired_slices.Add(self.dim2_slice_label)
        self.desired_slices.Add(self.dim2_slice_box)



    def FilePathToRead(self):
        self.fileboxlabel = wx.StaticBox(self,-1,'Input 2D Data File:')
        self.fileSizer=wx.StaticBoxSizer(self.fileboxlabel,wx.HORIZONTAL)

        self.text=wx.StaticText(self, -1, 'File Path:')
        
        self.filebox = wx.TextCtrl(self, size=(200, 22), style=wx.TE_PROCESS_ENTER)

        self.fileSizer.Add(self.text)
        self.fileSizer.Add(self.filebox)

    
    def OnButtonRead(self, event):

        spectrumfile = self.dirbox + str(self.filebox.GetValue())
        # print(spectrumfile)
        # sys.exit()

        shared_data = None
        shared_dic = None
        if getattr(self, "store", None) is not None:
            phasing_cache = self.store.spectra.get("phasing", {})
            shared_data = phasing_cache.get("data")
            shared_dic = phasing_cache.get("dic")
            if shared_data is None:
                shared_data = getattr(self.store, "data", None)
            if shared_dic is None:
                shared_dic = getattr(self.store, "dic", None)
        if shared_data is not None and shared_dic is not None:
            dic = shared_dic
            self.twod_data = shared_data
        else:
            # COMPATIBILITY FALLBACK DISABLED:
            # dic,self.twod_data = ng.pipe.read(spectrumfile)
            raise RuntimeError("Phasing requires spectrum data in data_store")
        print(dic)
        print(self.twod_data)

        if getattr(self, "store", None) is not None:
            self.store.save_spectrum("phasing", dic=dic, data=self.twod_data, spectrumfile=spectrumfile, dim=2)
            if getattr(self.store, "data", None) is None and getattr(self.store, "dic", None) is None:
                self.store.save_spectrum("raw", dic=dic, data=self.twod_data, spectrumfile=spectrumfile, dim=2)

        # self.twod_data = numpy.fabs(self.twod_data)
        self.uc0 = ng.pipe.make_uc(dic, self.twod_data, dim=0)
        x0,x1=self.uc0.ppm_limits()
        self.uc0.ppms_scale=numpy.linspace(x0, x1, int(self.uc0._size))
        ppm_dim1 = self.uc0.ppm_scale()

        self.uc1 = ng.pipe.make_uc(dic, self.twod_data, dim=1)
        x0,x1=self.uc1.ppm_limits()
        self.uc1.ppms_scale=numpy.linspace(x0, x1, int(self.uc1._size))
        ppm_dim2 = self.uc1.ppm_scale()

        self.slice_dim1 = str(self.dim1_slice_box.GetValue()) +' ppm'
        self.slice_dim2 = str(self.dim2_slice_box.GetValue()) +' ppm'
        print('self.slice_dim1 = {}'.format(self.slice_dim1))
        print('self.slice_dim2 = {}'.format(self.slice_dim2))

        # find ppm list index of elements that are closest to the 

        dim1_slice = self.twod_data[self.uc0(self.slice_dim2),:]

        p0, p1 = ng.process.proc_autophase.manual_ps(data=dim1_slice)

        phased_dim1_slice = ng.proc_base.ps(data=dim1_slice, p0=p0, p1=p1)


        self.axes1.plot(ppm_dim2, phased_dim1_slice)
        self.axes1.set_title(dic["FDF1LABEL"])


        dim2_slice = self.twod_data[:,self.uc1(self.slice_dim1)]

        p0, p1 = ng.process.proc_autophase.manual_ps(data=dim2_slice)

        phased_dim2_slice = ng.proc_base.ps(data=dim2_slice, p0=p0, p1=p1)

        
        self.axes2.plot(ppm_dim1, phased_dim2_slice)
        self.axes2.set_title(dic["FDF2LABEL"])









    





   


   

 


   

       

    def onFocus(self, event):
        print("Pseudo has focus!")


   
