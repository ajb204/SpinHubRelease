# #!/usr/bin/python
# import wx,string,copy,math,numpy,os
# import matplotlib            #import matplotlib
# matplotlib.use('WXAgg')      #switch on the wxPython mode
# from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
# from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# import nmrglue as ng
# from matplotlib.figure import Figure
# from matplotlib.widgets import Slider
# # from ...uSTA.py.stdParse import std
# from file_reader import ReadBuild, ReadTrans, EvalConc
# from ...misc.errors import errorMessage
# # from .frameFeatures import drawing_box, contour_box

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


# class uSTA_sims_paramter():
#     def __init__(self, title, parent,  default, big = False):
#         self.label = wx.StaticText(parent, -1, title+':')
#         if big == False:
#             self.box = wx.TextCtrl(parent, size=(60,22), style=wx.TE_PROCESS_ENTER)
#         else:
#             self.box = wx.TextCtrl(parent, size=(150,22), style=wx.TE_PROCESS_ENTER)
#         self.box.SetValue(str(default))
#         self.hbox = wx.BoxSizer(wx.HORIZONTAL)
#         self.hbox.Add(self.label, 0, border=5, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)
#         self.hbox.Add(self.box, 0, border=5, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)
        

# class uSTA_sims_frame(wx.Panel):
#     """ The main frame of the application
#     """
#     title = 'Demo: wxPython with matplotlib'

#     def __init__(self,parent,tabOne):

#         wx.Panel.__init__(self,parent=parent)
#         self.thresh=tabOne.dmax*float(tabOne.threshBox.GetValue()) #get threshold from main tab
#         self.tabOne=tabOne    #get tabone panel from NMR tab
#         self.parent=parent    #get decon_tab main parent from notebook

#         self.create_main_panel()


#     def default_parameters(self):
#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'dfrq'))!='0'):
#             self.dfrq = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'dfrq')))
#         else:
#             self.dfrq = (600.)
        
#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'kd'))!='0'):
#             self.kd = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'kd')))
#         else:
#             self.kd = 1e-6

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'kex'))!='0'):
#             self.kex = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'kex')))
#         else:
#             self.kex = 5

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'tcE'))!='0'):
#             self.tcE = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'tcE')))
#         else:
#             self.tcE = 8e-8

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'tcG'))!='0'):
#             self.tcG = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'tcG')))
#         else:
#             self.tcG = 2.5e-9


#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISp'))!='0'):
#             self.rISp = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISp')))
#         else:
#             self.rISp = 1.5e-10

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISl'))!='0'):
#             self.rISl = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISl')))
#         else:
#             self.rISl =3e-10

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISmix'))!='0'):
#             self.rISmix = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rISmix')))
#         else:
#             self.rISmix =3e-10
        
#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'dw'))!='0'):
#             self.dw = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'dw')))
#         else:
#             self.dw = 1

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rbar'))!='0'):
#             self.rbar = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'rbar')))
#         else:
#             self.rbar = 0.318

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'pw'))!='0'):
#             self.pw = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'pw')))
#         else:
#             self.pw = 50000

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'w1'))!='0'):
#             self.w1 = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'w1')))
#         else:
#             self.w1 = 200

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'B1inhom'))!='0'):
#             self.B1inhom = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'B1inhom')))
#         else:
#             self.B1inhom = 0.1

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'fac'))!='0'):
#             self.fac = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'fac')))
#         else:
#             self.fac = 1.44

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'maxIterFit'))!='0'):
#             self.maxIterFit = (int(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'maxIterFit')))
#         else:
#             self.maxIterFit = 20
            
#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'maxIterBoot'))!='0'):
#             self.maxIterBoot = (int(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'maxIterBoot')))
#         else:
#             self.maxIterBoot = 20

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'bootGrid'))!='0'):
#             self.bootGrid = (int(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'bootGrid')))
#         else:
#             self.bootGrid = 20

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'bootFac'))!='0'):
#             self.bootFac = (int(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'bootFac')))
#         else:
#             self.bootFac = 10

#         if(str(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'protConc'))!='0'):
#             self.protein_conc = (float(self.parent.tabOne.ParseFlt(self.parent.deconParFile,'protConc')))
#         else:
#             self.protein_conc = 1e-5

        
        

#         if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'kd_indir'))!='0'):
#             self.infile = (str(self.parent.tabOne.Parse(self.parent.deconParFile,'kd_indir')))
#         else:
#             self.infile = './raw/kD/'
    
#     def on_fitting_box(self, event):
#         if self.fitter_box.IsChecked() or self.kex_error_box.IsChecked() or self.kd_error_box.IsChecked() or self.twoD_error_box.IsChecked() or self.jiggler_box.IsChecked():
#             self.fit_button.Enable()



#     def kinetic_box(self):
#         self.kinLbl = wx.StaticBox(self,-1,'uSTA Kinetic Initial Conditions:')
#         self.kinSizer=wx.StaticBoxSizer(self.kinLbl,wx.VERTICAL)
    
#         self.default_parameters()

#         self.kd_control = uSTA_sims_paramter('Kd', self, self.kd, big=True)
#         self.kex_control = uSTA_sims_paramter('kex', self, self.kex, big=True)
#         self.tcG_control = uSTA_sims_paramter('tcG', self, self.tcG, big=True)
#         self.tcE_control = uSTA_sims_paramter('tcE', self, self.tcE, big=True)
#         self.risp_control = uSTA_sims_paramter('rISp', self, self.rISp, big=True)
#         self.risl_control = uSTA_sims_paramter('rISl', self, self.rISl, big=True)
#         self.rismix_control = uSTA_sims_paramter('rISmix', self, self.rISmix, big=True)
#         self.dw_control = uSTA_sims_paramter('dw', self, self.dw, big=True)
#         self.rbar_control = uSTA_sims_paramter('rbar', self, self.rbar, big=True)
#         self.fac_control = uSTA_sims_paramter('fac', self, self.fac, big=True)
    
#         self.kinSizer.Add(self.kd_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.kex_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.tcG_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.tcE_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.risp_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.risl_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.rismix_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.dw_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.rbar_control.hbox, 0, border=5, flag=self.flags)
#         self.kinSizer.Add(self.fac_control.hbox, 0, border=5, flag=self.flags)

#         self.kinSizer.AddSpacer(5)

#     def spec_box(self):
#         self.specLbl = wx.StaticBox(self,-1,'uSTA Spectrometer Initial Conditions:')
#         self.specSizer=wx.StaticBoxSizer(self.specLbl,wx.VERTICAL)
    

#         self.dfrq_control = uSTA_sims_paramter('Spectrometer\nFrequency (Mhz)', self, self.dfrq)
#         self.w1_control = uSTA_sims_paramter('Saturation\nField (hz)', self, self.w1)
#         self.pw_control = uSTA_sims_paramter('Saturation\nPulse width (us)', self, self.pw)
#         self.B1Inhom_control = uSTA_sims_paramter('B1 Inhomogeneity (%)', self, self.B1inhom)
#         listy=['Gaussian','Seduce']

#         self.pulse_shape_label = wx.StaticText(self, -1, 'Pulse shape:')
#         self.pulse_shape=wx.ComboBox(self, -1,choices=listy, style=wx.CB_READONLY, size=(-1,22))
#         self.pulse_shape_box = wx.BoxSizer(wx.HORIZONTAL)
#         self.pulse_shape_box.Add(self.pulse_shape_label, 0, border=5, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)
#         self.pulse_shape_box.Add(self.pulse_shape, 0, border=5, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)

#         self.specSizer.Add(self.dfrq_control.hbox, 0, border=5, flag=self.flags)
#         self.specSizer.Add(self.w1_control.hbox, 0, border=5, flag=self.flags)
#         self.specSizer.Add(self.pw_control.hbox, 0, border=5, flag=self.flags)
#         self.specSizer.Add(self.B1Inhom_control.hbox, 0, border=5, flag=self.flags)
#         self.specSizer.Add(self.pulse_shape_box, 0, border=5, flag=self.flags)
        

#         self.specSizer.AddSpacer(5)


    
#     def fitting_box(self):
#         self.fitLbl = wx.StaticBox(self,-1,'uSTA Fitting Parameters:')
#         self.fitSizer=wx.StaticBoxSizer(self.fitLbl,wx.VERTICAL)
    

#         self.maxIterFit_control = uSTA_sims_paramter('Max Iterations\nFitting', self, self.maxIterFit)
#         self.maxIterBoot_control = uSTA_sims_paramter('Max Iterations\nErrors', self, self.maxIterBoot)
#         self.bootFac_control = uSTA_sims_paramter('Errors factor', self, self.bootFac)
#         self.bootGrid_control = uSTA_sims_paramter('Errors points', self, self.bootGrid)

        

#         self.fitSizer.Add(self.maxIterFit_control.hbox, 0, border=5, flag=self.flags)
#         self.fitSizer.Add(self.maxIterBoot_control.hbox, 0, border=5, flag=self.flags)
#         self.fitSizer.Add(self.bootFac_control.hbox, 0, border=5, flag=self.flags)
#         self.fitSizer.Add(self.bootGrid_control.hbox, 0, border=5, flag=self.flags)
        

#         self.fitSizer.AddSpacer(5)



#     def options_box(self):
#         self.optionsLbl = wx.StaticBox(self,-1,'uSTA Fitting:')
#         self.optionsSizer=wx.StaticBoxSizer(self.optionsLbl,wx.VERTICAL)
    

        

#         self.file_control = uSTA_sims_paramter('Folder', self, self.infile)
#         self.file_button = wx.Button(self, -1, '...', size=(40,22))
#         self.file_control.hbox.Add(self.file_button, 0, border=5, flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)


#         self.fitter_box = wx.CheckBox(self, -1, 'Fitting')
#         self.jiggler_box = wx.CheckBox(self, -1, 'Jiggler')
#         self.kex_error_box = wx.CheckBox(self, -1, 'kex Errors')
#         self.kd_error_box = wx.CheckBox(self, -1, 'Kd Errors')
#         self.twoD_error_box = wx.CheckBox(self, -1, 'kex/Kd 2D Errors')

#         self.fitter_box.Bind(wx.EVT_CHECKBOX, self.on_fitting_box)
#         self.jiggler_box.Bind(wx.EVT_CHECKBOX, self.on_fitting_box)
#         self.kex_error_box.Bind(wx.EVT_CHECKBOX, self.on_fitting_box)
#         self.kd_error_box.Bind(wx.EVT_CHECKBOX, self.on_fitting_box)
#         self.twoD_error_box.Bind(wx.EVT_CHECKBOX, self.on_fitting_box)

#         self.read_button = wx.Button(self, -1, 'Read')
#         self.sim_button = wx.Button(self, -1, 'Sim')
#         self.fit_button = wx.Button(self, -1, 'Fit!')
#         self.save_button = wx.Button(self, -1, 'Save')
    
#         self.optionsSizer.Add(self.file_control.hbox, 0, border=5, flag=self.flags)

#         self.optionsSizer.Add(self.fitter_box, 0, border=5, flag=self.flags)
#         self.optionsSizer.Add(self.jiggler_box, 0, border=5, flag=self.flags)
#         self.optionsSizer.Add(self.kex_error_box, 0, border=5, flag=self.flags)
#         self.optionsSizer.Add(self.kd_error_box, 0, border=5, flag=self.flags)
#         self.optionsSizer.Add(self.twoD_error_box, 0, border=5, flag=self.flags)
        
#         self.optionsSizer.Add(self.read_button, 0, border=5, flag=self.flags)
#         self.read_button.Bind(wx.EVT_BUTTON, self.onButtonRead)
#         self.optionsSizer.Add(self.sim_button, 0, border=5, flag=self.flags)
#         self.sim_button.Bind(wx.EVT_BUTTON, self.on_sim_button)
#         self.optionsSizer.Add(self.fit_button, 0, border=5, flag=self.flags)
#         self.fit_button.Bind(wx.EVT_BUTTON, self.on_fit_button)
#         self.optionsSizer.Add(self.save_button, 0, border=5, flag=self.flags)
#         self.save_button.Bind(wx.EVT_BUTTON, self.on_save_button)


#         self.optionsSizer.AddSpacer(5)

#     def on_fit_button(self, event):
#         self.update_params()
#         self.update_ligand_params()
#         self.plotting_errors = False
#         if self.fitter_box.IsChecked():
#             self.inst.FITflag = True
#         else:
#             self.inst.FITflag = False
#         if self.jiggler_box.IsChecked():
#             self.inst.JIGGLEflag = True
#         else:
#             self.inst.JIGGLEflag = False
        
#         if self.kex_error_box.IsChecked():
#             self.plotting_errors = True
#             self.inst.BOOTflag = True
#             self.inst.PARALLEL = True
#             self.inst.ncpus=8
#             self.inst.chiplot_kex = True
            

#         elif self.kd_error_box.IsChecked():
#             self.plotting_errors = True
#             self.inst.BOOTflag = True
#             self.inst.PARALLEL = True
#             self.inst.ncpus=8
#             self.inst.chiplot_kd = True
#         else:
#             self.inst.BOOTflag = False
#             self.inst.chiplot_kd = False
            
#         if self.twoD_error_box.IsChecked():
#             self.inst.BOOT2dflag = True
#         else:
#             self.inst.BOOT2dflag = False

#         import threading
#         anEVT_CALCULATED = wx.NewEventType()
#         EVT_CALCULATED = wx.PyEventBinder(anEVT_CALCULATED, 1)

#         self.Bind(EVT_CALCULATED, self.onFitting)
#         th = threading.Thread(target=self.inst.Go, args=(self, anEVT_CALCULATED))
#         th.start()
#         # self.inst.Go(self.axis)

#     def onFitting(self, event):  
#         ''' this is where your thread comes back '''
#         print('simming')
#         self.plot_sim()
#         if self.inst.running_errors == False:
#             self.save()

#     def on_save_button(self, event):
#         self.save()


#     def sample_box(self):
#         self.sampleLbl = wx.StaticBox(self,-1,'uSTA Sample parameters:')
#         self.sampleSizer=wx.StaticBoxSizer(self.sampleLbl,wx.VERTICAL)

#         self.protein_conc_control = uSTA_sims_paramter('Protein Conc', self, self.protein_conc)
#         listy = ['']

#         self.protein_label_label=wx.StaticText(self, -1, 'Protein Label:', size=(-1,22))
#         self.protein_label_control=wx.ComboBox(self, -1,choices=listy, style=wx.CB_READONLY, size=(80,22))
#         self.protein_label_control.Bind(wx.EVT_COMBOBOX, self.on_pro_ind_selection)

#         self.protein_label_box = wx.BoxSizer(wx.HORIZONTAL)
#         self.protein_label_box.Add(self.protein_label_label, 0, border=5, flag=wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)
#         self.protein_label_box.Add(self.protein_label_control, 0, border=5, flag= wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)


#         self.ligand_index_label=wx.StaticText(self, -1, 'Ligand Index:', size=(-1,22))
#         self.ligand_index_control=wx.ComboBox(self, -1,choices=listy, style=wx.CB_READONLY, size=(80,22))
#         self.ligand_index_control.Bind(wx.EVT_COMBOBOX, self.plot_raw)
        
#         self.ligand_index_box = wx.BoxSizer(wx.HORIZONTAL)
#         self.ligand_index_box.Add(self.ligand_index_label, 0, border=5, flag=wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)
#         self.ligand_index_box.Add(self.ligand_index_control, 0, border=5, flag= wx.BOTTOM | wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.TOP)


#         self.sampleSizer.Add(self.protein_label_box, 0, border=5, flag=self.flags)
#         self.sampleSizer.Add(self.protein_conc_control.hbox, 0, border=5, flag=self.flags)
#         self.protein_conc_control.box.Bind(wx.EVT_TEXT, self.on_pro_conc_change)
#         self.sampleSizer.Add(self.ligand_index_box, 0, border=5, flag=self.flags)

        
#         self.sampleSizer.AddSpacer(5)


#     def pre_read_disable(self):
#         self.protein_label_control.Disable()
#         self.protein_conc_control.box.Disable()
#         self.protein_conc_control.label.Disable()
#         self.ligand_index_control.Disable()
#         self.ligand_index_label.Disable()
#         self.protein_label_label.Disable()
#         self.fitter_box.Disable()
#         self.jiggler_box.Disable()
#         self.kex_error_box.Disable()
#         self.kd_error_box.Disable()
#         self.twoD_error_box.Disable()
#         self.fit_button.Disable()
#         self.sim_button.Disable()


#     def plot_raw(self, event):
#         self.axis.cla()
#         self.axis.set_xlabel('Mixing Time (s)')
#         self.axis.set_ylabel('uSTA Transfer Efficiency (%)')
#         concy = numpy.log10(self.inst.Conc)
#         concy = concy-numpy.min(concy)
#         self.raw_plots = []
#         for w in self.inst.w1:
#             for i2,key in enumerate(self.inst.concKeys):
#                 c = concy[i2]/numpy.max(concy)
#                 print(c)
#                 r = c
#                 g = 0.0
#                 b = 1-c
#                 for i,line in enumerate(self.inst.trans[key][0]):
#                     if(line.lab==self.ligand_index_control.GetValue()):
#                         self.raw_plots.append(self.axis.scatter(self.inst.pwslcest[i2], self.inst.rawBuild[w][self.ligand_index_control.GetValue()][i2], marker='x', color=(r,g,b), s=50, label = self.inst.concKeys[i2].split('_')[1]))
#         self.axis.legend()
#         self.canvas.draw()

#     def plot_sim(self):
#         for line in self.axis.lines:
#             line.remove()
#         concy = numpy.log10(self.inst.Conc)
#         concy = concy-numpy.min(concy)
#         self.sim_plots = []
#         for w in self.inst.w1:
#             for i2,key in enumerate(self.inst.concKeys):
#                 c = concy[i2]/numpy.max(concy)
#                 print(c)
#                 r = c
#                 g = 0.0
#                 b = 1-c
#                 for i,line in enumerate(self.inst.trans[key][0]):
#                     if(line.lab==self.ligand_index_control.GetValue()):
#                         self.sim_plots.append(self.axis.plot(self.inst.pwslcest[i2], self.inst.calcBuild_save[w][self.ligand_index_control.GetValue()][i2], color=(r,g,b))[0])
#         self.axis.legend()
#         self.update_fields()
#         ymin_1 = min(self.inst.calcBuild_save[w][self.ligand_index_control.GetValue()][i2])
#         ymin_2 = min(self.inst.rawBuild[w][self.ligand_index_control.GetValue()][i2])
#         ymin = min(ymin_1, ymin_2)
#         ymax_1 = max(self.inst.calcBuild_save[w][self.ligand_index_control.GetValue()][i2])
#         ymax_2 = max(self.inst.rawBuild[w][self.ligand_index_control.GetValue()][i2])
#         ymax = max(ymax_1, ymax_2)
#         self.axis.set_ylim(0, ymax)
#         self.axis_chi.cla()
#         self.axis_chi.plot(range(len(self.inst.chi2_array)), self.inst.chi2_array)
#         self.axis_chi.set_yscale('log')
#         self.axis_chi.set_xlabel('Iteration')
#         self.axis_chi.set_ylabel('Chi Squared')
#         self.axis_chi.text(0.6, 0.9, "chi2: {:.5e}".format(self.inst.chi2_array[-1]), transform=self.axis_chi.transAxes)
#         self.axis_chi.text(0.6, 0.8, "Best chi2: {:.5e}".format(min(self.inst.chi2_array)), transform=self.axis_chi.transAxes)

#         if self.plotting_errors == True:
#             self.axis_err.cla()
#             self.axis_err.scatter(self.inst.chi2_error_array_x, self.inst.chi2_error_array_y)

#         self.canvas.draw()

#     def update_fields(self):
#         self.kd_control.box.SetValue("{:.5e}".format(self.inst.Kd))
#         self.kex_control.box.SetValue("{:.5e}".format(self.inst.kex))
#         self.tcE_control.box.SetValue("{:.5e}".format(self.inst.tcE))
#         self.tcG_control.box.SetValue("{:.5e}".format(self.inst.tcG))
#         self.risl_control.box.SetValue("{:.5e}".format(self.inst.rISl[self.ligand_index_control.GetValue()]))
#         self.rismix_control.box.SetValue("{:.5e}".format(self.inst.rISmix[self.ligand_index_control.GetValue()]))
#         self.fac_control.box.SetValue("{:.5e}".format(self.inst.fac[self.ligand_index_control.GetValue()]))
#         self.dw_control.box.SetValue("{:.5e}".format(self.inst.dw))
#         self.rbar_control.box.SetValue("{:.5e}".format(self.inst.rbar))
        

#         self.risp_control.box.SetValue("{:.5e}".format(self.inst.rISp))


#     def setup_axis(self):
#         self.fig.clf()
#         self.axis = self.fig.add_subplot(121)
#         self.axis.set_ylabel('uSTA Transfer Efficiency')
#         self.axis.set_xlabel('Mixing Time (s)')
#         self.axis_chi = self.fig.add_subplot(222)
#         self.axis_err = self.fig.add_subplot(224)



#     def onButtonRead(self, event):
#         if(str(self.parent.tabOne.Parse(self.parent.deconParFile,'ligInd'))!='0'):
#             try:
#                 self.ligand_index_control.SetValue((str(self.parent.tabOne.Parse(self.parent.deconParFile,'ligInd'))))
#             except:
#                 errorMessage('ligand index Not found: reverting to default')
#         self.initial_read()
#         self.setup_axis()
#         self.plot_raw(event)
#         self.post_read_enable()
        
        


#     def extract_protein_names_ligand_concs(self):
#         self.protein_names = []
#         self.ligand_concs = {}
#         self.protein_label_control.Clear()
        
#         for x in os.listdir(self.infile):
#             print(x)
#             if 'data.Kd.build' in x:
#                 stub = x.split('.data.Kd.build')[0].split('raw.')[1]
#                 if os.path.isfile(self.infile + os.sep + 'raw.'+ stub + '.data.trans.kd.out'):
#                     try:
#                         protein_name = stub.split('_')[0]
#                         ligand_conc = stub.split('_')[1]
#                         print(stub)
#                         if protein_name not in self.protein_names:
#                             self.protein_names.append(protein_name)
#                             self.protein_label_control.Append(protein_name)
#                         if protein_name in self.ligand_concs.keys():
#                             self.ligand_concs[protein_name].append(ligand_conc)
#                         else:
#                             self.ligand_concs[protein_name] = [ligand_conc]
#                     except:
#                         pass

#     def extract_ligand_indices(self):
#         self.ligand_indices = []
#         self.ligand_index_control.Clear()
#         for x in os.listdir(self.infile):
#             if 'data.Kd.build' in x:
#                 stub = x.split('.data.Kd.build')[0].split('raw.')[1]
#                 if os.path.isfile(self.infile + os.sep + 'raw.'+ stub + '.data.trans.kd.out'):
#                     for line in open(self.infile + os.sep + 'raw.'+ stub + '.data.trans.kd.out', 'r'):
#                         fields = line.split()
#                         if len(fields) > 0:
#                             ligand_peak_name = fields[1]
#                             if ligand_peak_name not in self.ligand_indices:
#                                 self.ligand_indices.append(ligand_peak_name)
#                                 self.ligand_index_control.Append(ligand_peak_name)
#         # self.ligand_index_control.SetSelection(0)


#     def update_params(self):
#         self.inst.Kd= float(self.kd_control.box.GetValue())
#         # inst.fac['H-9']= self.fac
#         self.inst.kex = float(self.kex_control.box.GetValue())
#         self.inst.tcG = float(self.tcG_control.box.GetValue()) #1.0933760713396435e-08
#         self.inst.tcE = float(self.tcE_control.box.GetValue()) #2.6378290566204117e-05
#         self.inst.rISp = float(self.risp_control.box.GetValue()) #8.111206751650106e-09
#         self.inst.w1 = float(self.w1_control.box.GetValue()),
#         self.inst.dw= float(self.dw_control.box.GetValue())
#         self.inst.pw = float(self.pw_control.box.GetValue())
#         self.inst.B1inhom = float(self.B1Inhom_control.box.GetValue())
#         self.inst.dfrq = float(self.dfrq_control.box.GetValue())
#         self.inst.shape = self.pulse_shape.GetValue()
#         self.inst.rbar=float(self.rbar_control.box.GetValue())
#         self.inst.nsteps=40
#         self.inst.maxIterBOOT=int(self.maxIterBoot_control.box.GetValue())
#         self.inst.maxIterFIT=int(self.maxIterFit_control.box.GetValue())
#         self.inst.bootFac=int(self.bootFac_control.box.GetValue())
#         self.inst.bootGrid=int(self.bootGrid_control.box.GetValue())




#     def set_uncontrollable_params(self):
#         self.inst.ncpus=8
        
        
#         self.inst.dwflag=True
#         self.inst.rFree=False
#         self.inst.inhom='n' 
#         self.inst.pulse_flg = 'y'
#         if self.inst.shape == 'Gaussian':
#             print('gaussian')
#             self.inst.shape='GAUSSIAN'
            
#         self.inst.tcflag=True
#         self.inst.Rflag=True
#         self.inst.dw3=0.0   #proteinF (excitation on resonance)
#         self.inst.dw4=0.0 


#     def initial_read(self):
#         self.extract_protein_names_ligand_concs()
#         self.extract_ligand_indices()
        

#         # return
#         self.inst=std(self.infile)
        

#         for x in self.protein_names:
#             self.inst.proInd[x] = 0

#         self.inst.index = self.ligand_indices
#         self.inst.index = 'H-43',

#         self.update_params()
#         self.set_uncontrollable_params()

        
        

        
        
        
#         self.inst.Go_initial()

#         # self.tag, self.mol = ReadBuild(self.infile, self.protein_names) 
#         # print(self.tag, self.mol)
        
#     def post_read_enable(self):
#         self.protein_label_control.SetSelection(0)
#         self.protein_label_control.Enable()
#         self.protein_label_label.Enable()
#         self.protein_conc_control.label.Enable()
#         self.protein_conc_control.box.Enable()
#         self.ligand_index_control.Enable()
#         self.ligand_index_label.Enable()
#         self.fitter_box.Enable()
#         self.jiggler_box.Enable()
#         self.kex_error_box.Enable()
#         self.kd_error_box.Enable()
#         self.twoD_error_box.Enable()
#         self.sim_button.Enable()

    


#     def save(self):
#         write={}
#         write['kd'] = self.kd_control.box.GetValue()
#         write['kex'] = self.kex_control.box.GetValue()
#         write['tcG'] = self.tcG_control.box.GetValue()
#         write['tcE'] = self.tcE_control.box.GetValue()
#         write['rISp'] = self.risp_control.box.GetValue()
#         write['rISl']   = self.risl_control.box.GetValue()
#         write['rISmix'] = self.rismix_control.box.GetValue()
#         write['dw'] = self.dw_control.box.GetValue()
#         write['rbar'] = self.rbar_control.box.GetValue()
#         write['fac'] = self.fac_control.box.GetValue()
#         write['w1'] = self.w1_control.box.GetValue()
#         write['pw'] = self.pw_control.box.GetValue()
#         write['B1inhom'] = self.B1Inhom_control.box.GetValue()
#         write['dfrq'] = self.dfrq_control.box.GetValue()
#         write['shape'] = self.pulse_shape.GetValue()
#         write['maxIterFit'] = self.maxIterFit_control.box.GetValue()
#         write['maxIterBoot'] = self.maxIterBoot_control.box.GetValue()
#         write['bootFac'] = self.bootFac_control.box.GetValue()
#         write['bootGrid'] = self.bootGrid_control.box.GetValue()
#         write['protLabel'] = self.protein_label_control.GetValue()
#         write['protConc'] = self.protein_conc_control.box.GetValue()
#         write['ligInd'] = self.ligand_index_control.GetValue()
#         write['kd_indir'] = self.file_control.box.GetValue()




#         dec=[]
#         if(os.path.exists(self.parent.tabOne.deconParFile)):
#             inny=open(self.parent.tabOne.deconParFile)
#             for line in inny.readlines():
#                dec.append(line)
#             inny.close()

#         outy=open(os.path.join(self.parent.tabOne.dirBox.GetValue(),self.parent.tabOne.deconParFile),'w')
#         print('saving:',os.path.join(self.parent.tabOne.dirBox.GetValue(),self.parent.tabOne.deconParFile))
#         for de in dec:
#             test=de.split()
#             if(len(test)>0):

#                 if(test[0] in list(write.keys())):
#                     outy.write('%s = %s ' % (test[0],write[test[0]]))
#                     for j in range(len(test)-3):
#                         outy.write(' %s' % test[j+3])
#                     outy.write('\n')
#                     del write[test[0]]
#                 else:
#                     outy.write(de)
#             else:
#                 outy.write(de)

#         for key,vals in list(write.items()):
#             try:
#                 outy.write('%s = %s\n' % (key,vals[0]))
#             except:
#                 print('problem with ',key)

#         outy.close()

#     def on_pro_ind_selection(self, event):
#         self.protein_conc_control.box.SetValue(str(self.inst.proInd[self.protein_label_control.GetValue()]))


#     def on_pro_conc_change(self, event):
#         self.inst.proInd[self.protein_label_control.GetValue()] = float(self.protein_conc_control.box.GetValue())


#     def update_ligand_params(self):
#         self.inst.rISl[self.ligand_index_control.GetValue()]= float(self.risl_control.box.GetValue()) #2.1932532263494214e-10
#         self.inst.rISmix[self.ligand_index_control.GetValue()]= float(self.rismix_control.box.GetValue()) #1.2391092578874141e-08
#         self.inst.fac[self.ligand_index_control.GetValue()] = float(self.fac_control.box.GetValue())
#         # for x in self.protein_names:
#         #     self.inst.proInd[x] = float(self.protein_conc_control.box.GetValue())
#         self.inst.index = [self.ligand_index_control.GetValue()]

#     def on_sim_button(self, event):
#         self.plotting_errors = False

#         for x in self.sim_plots:
#             x.remove()
#         self.update_params()
#         self.update_ligand_params()
        
#         self.inst.Go()
        
#         self.plot_sim()

#     def create_main_panel(self):
#         """ Creates the main panel with all the controls on it:
#              * mpl canvas
#              * mpl navigation toolbar
#              * Control panel for interaction
#         """

#         ## Initialise matplotlib main figure
#         self.fig = Figure()
#         self.canvas = FigCanvas(self, -1, self.fig)
#         self.canvas.SetMinSize(wx.Size(400,400))
#         self.toolbar = NavigationToolbar(self.canvas)
#         self.sim_plots = []
#         self.raw_plots = []
#         ## Initialise matplotlib slider figure
#         # self.slider_fig = Figure()
#         # self.slider_canvas = FigCanvas(self, -1, self.slider_fig)
#         # self.slider_canvas.SetMinSize(wx.Size(-1,20))

#         ## Adding our control boxes
#         self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.LEFT | wx.ALIGN_RIGHT | wx.RIGHT
#         self.kinetic_box()
#         self.spec_box()
#         self.fitting_box()
#         self.sample_box()
#         self.options_box()
#         self.pre_read_disable()

#         ## Piece together the control boxes
#         self.hbox = wx.BoxSizer(wx.VERTICAL)
#         self.hbox.AddSpacer(5)
#         self.hbox.Add(self.kinSizer, flag = wx.GROW)
#         self.hbox.AddSpacer(5)
#         self.hbox.Add(self.specSizer, flag = wx.GROW)
        


#         self.matplotlib_vbox = wx.BoxSizer(wx.VERTICAL)

#         ## Main vertical sizer
#         self.vbox = wx.BoxSizer(wx.HORIZONTAL)
        
#         self.hbox2 = wx.BoxSizer(wx.VERTICAL)
#         self.hbox2.AddSpacer(5)
#         self.hbox2.Add(self.fitSizer, flag = wx.GROW)
#         self.hbox2.AddSpacer(5)
#         self.hbox2.Add(self.sampleSizer, flag = wx.GROW)
#         self.hbox2.AddSpacer(5)
#         self.hbox2.Add(self.optionsSizer, flag = wx.GROW)


#         self.vbox.AddSpacer(10)
#         self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT, border=0)
#         self.vbox.AddSpacer(10)
#         self.vbox.Add(self.hbox2, 0, flag = wx.ALIGN_LEFT, border=0)
#         self.vbox.AddSpacer(10)

        
#         self.matplotlib_vbox.Add(self.canvas, 5, wx.LEFT | wx.TOP | wx.GROW)
#         self.matplotlib_vbox.Add(self.toolbar, 0, wx.EXPAND)
#         self.vbox.Add(self.matplotlib_vbox, 5, flag = wx.ALIGN_LEFT | wx.GROW, border=0)

#         ## Define transpose axes:
        
#         self.resized = False
        

#         ## Final layout adjustments
#         self.fig.tight_layout()
#         self.SetSizer(self.vbox)
#         self.vbox.Fit(self)

    
    