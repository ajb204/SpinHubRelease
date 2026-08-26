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

import numpy,os,wx,glob,re
import nmrglue as ng

from spinDecon.project.parameter_store import read_structured_parameter_file, write_structured_parameter_file
from spinDecon.gui.dialogs import text_viewer as textEdit
from spinDecon.domain.pseudo_axis import PseudoAxisTable, PseudoAxisError, pseudo_axis_path, load_saved_column, save_selected_column
from spinDecon.gui.dialogs.pseudo_axis import show_pseudo_axis_table
from spinDecon.project.parameter_store import parse_value, update_parameter_file

import wx,copy,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.collections import LineCollection            
import matplotlib.cm as cm
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin


import wx.lib.mixins.listctrl  as  listmix


from scipy.optimize import leastsq

class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
    ''' TextEditMixin allows any column to be edited. '''
    
    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """Constructor"""
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)



class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty
    def GetListCtrl(self):
        return self


    
class DecayMan(wx.App):
    def __init__(self,inherit,pth='',auto_prepare=False):
        self.frame_ProcessFrame=DecayFrame(None,30,'Decay Analysis',inherit,pth=pth)
        #FGA added
        self.frame_ProcessFrame.Centre(direction=wx.BOTH)
        self.frame_ProcessFrame.Show(True)
        # Workflow Show analysis: regenerate curves, fit all peaks, then show R.
        if auto_prepare:
            wx.CallAfter(self.frame_ProcessFrame.prepare_workflow_analysis)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

class DecayFrame(wx.Frame):
    
    def __init__(self,parent,id,title,inherit,pth):
        #wx.Panel.__init__(self, parent=parent)
        self.parent=inherit
        self.pth=pth
        self.WXV=int(wx.__version__.split('.')[0])


        #FGA changed
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.Point(258, 184), size=wx.Size(800, 480),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'MAGMA results ...')
        #self.SetClientSize(wx.Size(900, 280))
        monitorWidth, monitorHeight = wx.GetDisplaySize()
        initial_w = min(900, max(700, int((monitorWidth - 120) * 0.75)))

        # Size the decay frame relative to the actual top-level main window, not
        # just the display.  ``inherit`` is sometimes a child panel/frame, so
        # resolve its top-level parent before deciding the opening height.
        work_h = wx.GetClientDisplayRect().height
        try:
            main_window = wx.GetTopLevelParent(inherit)
            main_h = main_window.GetSize().height if main_window else inherit.GetSize().height
        except Exception:
            main_h = work_h
        max_open_h = max(540, min(main_h, work_h))
        initial_h = min(760, max_open_h)

        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
              pos=wx.DefaultPosition, size=(initial_w, initial_h),
              style=wx.DEFAULT_FRAME_STYLE, title='Decay Analysis')
        self.SetMinSize((680, min(540, max_open_h)))
        self.SetMaxSize((-1, max_open_h))
        self.statusBar = self.CreateStatusBar(1)
        self.statusBar.SetStatusText('Ready')

        self.SetBackgroundColour('WHITE')
        self.panel=wx.Panel(self,-1)
        
        ########

        # Analysis products belong under WorkingDir/SpecPath.  Historically
        # this frame used <project>/raw for both acquisition input and analysis
        # output; keep acquisition data separate from the analysis workspace.
        self.raw = self._analysis_spec_dir(inherit, pth)
        self.acquisition_dir = self._acquisition_dir(inherit, pth)
        self.fuda_dir = inherit._fuda_dir()

        self.savefile=os.path.join(self.raw, 'decay', 'frame.save')
        if(os.path.exists(self.savefile)==0):
            self.PathExists((self.raw+'/decay',))
            outy=open(self.savefile,'w')
            outy.close()

            
        self.cpmgLocal={}
        self.plotMode = 'Peak'

        #self.lc=SortedListCtrl(panel,self.corrDict)

        self.datasets=AutoWidthListCtrl(self)
        #self.datasets.SetMinSize((650,300))
        self.datasets.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)
        self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        self.datasets.Bind(wx.EVT_LIST_COL_CLICK, self.OnButtonSort)
        self.datasets.Bind(wx.EVT_KEY_DOWN, self.OnDatasetKeyDown)
        
        
        self.buttonRefresh = wx.Button(self, label="Refresh")
        self.buttonRefresh.Bind(wx.EVT_BUTTON,self.OnButtonRefresh)
        

        #self.buttonAddDataset = wx.Button(self, label="Add Dataset")
        #self.buttonAddDataset.Bind(wx.EVT_BUTTON,self.OnButtonAddDataset)
        #self.buttonRemDataset = wx.Button(self, label="Remove Dataset")
        #self.buttonRemDataset.Bind(wx.EVT_BUTTON,self.OnButtonRemDataset)
        
        
        self.buttonSave = wx.Button(self, label="Save")
        self.buttonSave.Bind(wx.EVT_BUTTON,self.OnButtonSave)

        self.buttonLoad = wx.Button(self, label="Load")
        self.buttonLoad.Bind(wx.EVT_BUTTON,self.OnButtonLoad)

        #self.buttonWipe = wx.Button(self, label="Wipe")
        #self.buttonWipe.Bind(wx.EVT_BUTTON,self.OnButtonWipe)

        self.buttonGuess = wx.Button(self, label="FitAll")
        self.buttonGuess.Bind(wx.EVT_BUTTON,self.OnButtonGuess)
        

        self.buttonClose = wx.Button(self, label="Close")
        self.buttonClose.Bind(wx.EVT_BUTTON,self.OnButtonClose)

        #self.buttonRun = wx.Button(self, label="Run Catia")
        #self.buttonRun.Bind(wx.EVT_BUTTON,self.OnButtonRun)

        self.buttonPeakConvert = wx.Button(self, label="Create Decay curves")
        self.buttonPeakConvert.Bind(wx.EVT_BUTTON,self.OnButtonPeakConvert)





        #self.fieldboxtxt  = wx.StaticText(self, label="Field:")
        #self.fieldbox = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        #self.fieldbox.SetValue("750")

        #self.tempboxtxt  = wx.StaticText(self, label="Temperature:")
        #self.tempbox = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        #self.tempbox.SetValue("20")

        
        #self.nucboxtxt  = wx.StaticText(self, label="Nucleus:")
        #self.nucbox = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        #self.nucbox.SetValue("N15")



        self.TimeT2BoxTxt = wx.StaticText(self, label="Time Mult:")
        self.TimeT2Box = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.TimeT2Box.SetValue(self._load_time_mult())
        self.TimeT2Box.Bind(wx.EVT_TEXT_ENTER, self.OnTimeMultChanged)
        self.TimeT2Box.Bind(wx.EVT_KILL_FOCUS, self.OnTimeMultChanged)

        # The pseudo-axis is defined by conversion and stored in pseudo_axis.tsv.
        # Analysis chooses a named column; the choice is shared/persisted at
        # project level so Decay and CPMG reopen on the same experimental axis.
        self.pseudoAxisTxt = wx.StaticText(self, label="Pseudo-axis:")
        self.pseudoAxisCombo = wx.ComboBox(self, -1, choices=[], style=wx.CB_READONLY)
        self.pseudoAxisOpen = wx.Button(self, label="View")
        self.pseudoAxisCombo.Bind(wx.EVT_COMBOBOX, self.OnPseudoAxisColumn)
        self.pseudoAxisOpen.Bind(wx.EVT_BUTTON, self.OnOpenPseudoAxisTable)
        self._load_pseudo_axis_choices()

        

        
        #self.fitLocal = AutoWidthListCtrl(self)
        #self.fitLocal.SetMinSize((300,300))
        


        
        #self.setLocal.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClickSet)
        #self.parLocal.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClickPar)

        
        self.listbox=wx.BoxSizer(wx.HORIZONTAL)
        #setbox = wx.StaticBox(self,-1,'Set Parameters:')
        #setboxS=wx.StaticBoxSizer(setbox,wx.VERTICAL)
        #setboxS.Add(self.setLocal)
        
        #locbox = wx.StaticBox(self,-1,'Local Parameters:')
        #locboxS=wx.StaticBoxSizer(locbox,wx.VERTICAL)
        #locboxS.Add(self.parLocal)


        #fitbox = wx.StaticBox(self,-1,'Fitted Parameters:')
        #fitboxS=wx.StaticBoxSizer(fitbox,wx.VERTICAL)
        #fitboxS.Add(self.fitLocal)

        #self.listbox.AddSpacer(10)
        #self.listbox.Add(setboxS)
        #self.listbox.AddSpacer(10)
        #self.listbox.Add(locboxS)

        #self.listbox.AddSpacer(10)
        #self.listbox.Add(fitboxS)
        #self.listbox.AddSpacer(10)
        
        
        #self.parlocl.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)
        #self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        #self.DoLocal(True)


        

        
        # Responsive analysis layout.  The results list is the primary
        # navigation control and the plot receives all spare window space.
        # Keep the controls in a deliberately narrow panel.  Putting the
        # left-hand sizer directly in the frame allowed the best sizes of its
        # children (especially the action-button grid) to make the whole left
        # side wider than intended.
        self.optionsPanel = wx.Panel(self)
        self.optionsPanel.SetMinSize((255, -1))
        self.vboxOpts = wx.BoxSizer(wx.VERTICAL)

        results_static = wx.StaticBox(self.optionsPanel, label="Datasets / fits")
        results_box = wx.StaticBoxSizer(results_static, wx.VERTICAL)
        self.datasets.Reparent(results_static)
        self.datasets.SetMinSize((235, 260))
        results_box.Add(self.datasets, 1, wx.EXPAND | wx.ALL, 6)
        self.vboxOpts.Add(results_box, 1, wx.EXPAND | wx.BOTTOM, 8)

        settings_static = wx.StaticBox(self.optionsPanel, label="Decay settings")
        settings_box = wx.StaticBoxSizer(settings_static, wx.VERTICAL)
        self.TimeT2BoxTxt.Reparent(settings_static)
        self.TimeT2Box.Reparent(settings_static)
        self.pseudoAxisTxt.Reparent(settings_static)
        self.pseudoAxisCombo.Reparent(settings_static)
        self.pseudoAxisOpen.Reparent(settings_static)
        mode_row = wx.BoxSizer(wx.HORIZONTAL)
        mode_row.Add(wx.StaticText(settings_static, label="Plot:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        self.plotModeRadio = wx.RadioBox(settings_static, choices=['Peak', 'R', 'A0'], majorDimension=3, style=wx.RA_SPECIFY_COLS)
        self.plotModeRadio.SetSelection(0)
        self.plotModeRadio.Bind(wx.EVT_RADIOBOX, self.OnPlotMode)
        mode_row.Add(self.plotModeRadio, 1, wx.EXPAND)
        settings_box.Add(mode_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)
        axis_row = wx.BoxSizer(wx.HORIZONTAL)
        axis_row.Add(self.pseudoAxisTxt, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        axis_row.Add(self.pseudoAxisCombo, 1, wx.EXPAND | wx.RIGHT, 4)
        axis_row.Add(self.pseudoAxisOpen, 0)
        settings_box.Add(axis_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)
        time_row = wx.BoxSizer(wx.HORIZONTAL)
        time_row.Add(self.TimeT2BoxTxt, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        time_row.Add(self.TimeT2Box, 1, wx.EXPAND)
        settings_box.Add(time_row, 0, wx.EXPAND | wx.ALL, 6)
        self.vboxOpts.Add(settings_box, 0, wx.EXPAND | wx.BOTTOM, 8)

        actions_static = wx.StaticBox(self.optionsPanel, label="Actions")
        actions_box = wx.StaticBoxSizer(actions_static, wx.VERTICAL)
        # A two-column button grid made this panel's best width dominate the
        # horizontal sizer ("Create Decay curves" is quite a wide button).
        # A single column costs a little height but keeps the plot wide.
        action_grid = wx.GridSizer(rows=0, cols=1, vgap=4, hgap=0)
        for button in (self.buttonRefresh, self.buttonPeakConvert, self.buttonGuess,
                       self.buttonSave, self.buttonLoad, self.buttonClose):
            button.Reparent(actions_static)
            action_grid.Add(button, 0, wx.EXPAND)
        actions_box.Add(action_grid, 0, wx.EXPAND | wx.ALL, 6)
        self.vboxOpts.Add(actions_box, 0, wx.EXPAND)

        self.optionsPanel.SetSizer(self.vboxOpts)
        self.sizerMain = wx.BoxSizer(wx.VERTICAL)
        self.fullSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fullSizer.Add(self.optionsPanel, 0, wx.EXPAND | wx.ALL, 10)
        self.create_main_panel()
        self.sizerMain.Add(self.fullSizer, 1, wx.EXPAND)

        self.OnButtonRefresh(True)
        self.OnButtonLoad(True)

        #self.create_main_panel()
        #self.draw_figure()
        #self.canvas.draw()
        
        self.SetSizer(self.sizerMain)
        self.Layout()
        self.Centre()

        #FGA added
        #if self.ref=='subgraphMode':
        #    self.rbox.Disable()

        self.Show(True)
        # Do not Fit() the frame here: that collapses the matplotlib canvas to
        # its best size and defeats responsive resizing.
        self.Layout()
        #self.Thaw()




    def DecayCalc(self):
        self.cpmgF=self.A0*numpy.exp(-1.*self.cpmgX*self.kex)
    def PackDecay(self):
        x=[]
        x.append(self.A0)
        x.append(self.kex)
        return x
    def UnpackDecay(self,x):
        # Preserve the fitted amplitude sign so negative-intensity peaks can
        # be represented by the single-exponential model.  The decay rate
        # remains non-negative.
        self.A0=x[0]
        self.kex=numpy.fabs(x[1])
    def ChiDecay(self,x):
        self.UnpackDecay(x)
        self.DecayCalc()
        return self.cpmgF-self.cpmgY


    #def cpmg_Line(self):
    #    self.cpmgL=numpy.ones(len(self.cpmgX))*self.R0line
    #def PackLine(self):
    #    x=[]
    #    x.append(self.R0line)
    #    return x
    #def UnpackLine(self,x):
    #    x=numpy.fabs(x)
    #    cnt=0
    #    self.R0line=x[cnt];cnt+=1
    #def ChiLine(self,x):
    #    self.UnpackLine(x)
    #    self.cpmg_Line()
    #    return self.cpmgL-self.cpmgY

    #https://en.wikipedia.org/wiki/Gyromagnetic_ratio
    #units of 10^6 rad s-1 T-1
    Gamma={}
    Gamma['15N']=27.116
    Gamma['1H']=267.522
    Gamma['19F']=251.815
    Gamma['13C']=67.2828
    Gamma['2H']=41.065
    Gamma['31P']=108.291
    
    def _set_status(self, message):
        """Show concise analysis progress in the Decay window."""
        bar = getattr(self, 'statusBar', None)
        if bar is not None:
            bar.SetStatusText(str(message))
            wx.YieldIfNeeded()

    def _project_parameter_path(self):
        """Return the system parameter file path without relying on cwd."""
        owner = getattr(self, 'parent', None)
        seen = set()
        while owner is not None and id(owner) not in seen:
            seen.add(id(owner))
            tab_one = getattr(owner, 'tabOne', None)
            if tab_one is not None:
                owner = tab_one
                break
            if getattr(owner, 'deconParFile', None):
                break
            owner = getattr(owner, 'parent', None)
        parfile = getattr(owner, 'deconParFile', 'deconParFile') if owner is not None else 'deconParFile'
        if os.path.isabs(parfile):
            return parfile
        working_dir = '.'
        dir_box = getattr(owner, 'dirBox', None) if owner is not None else None
        if dir_box is not None:
            try:
                working_dir = str(dir_box.GetValue() or '.').strip() or '.'
            except Exception:
                pass
        return os.path.join(working_dir, parfile)

    def _load_time_mult(self):
        """Load the project-wide Decay time multiplier, retaining 0.04 for old projects."""
        path = self._project_parameter_path()
        value = parse_value(path, 'decayTimeMult', default='0.04')
        try:
            float(value)
        except (TypeError, ValueError):
            value = '0.04'
        return str(value)

    def _save_time_mult(self):
        text = self.TimeT2Box.GetValue().strip()
        try:
            value = float(text)
        except ValueError:
            self._set_status('Time Mult must be numeric')
            return False
        path = self._project_parameter_path()
        update_parameter_file(path, {'decayTimeMult': text}, source_path=path)
        return value

    def OnTimeMultChanged(self, event):
        value = self._save_time_mult()
        if value is not False:
            self._set_status('Time multiplier: %g' % value)
        if event is not None:
            event.Skip()

    def _load_pseudo_axis_choices(self):
        self.pseudo_axis_file = pseudo_axis_path(self.parent, fallback_spec_dir=getattr(self, 'raw', None))
        try:
            table = PseudoAxisTable.load(self.pseudo_axis_file)
            self.pseudo_axis_table = table
            self.pseudoAxisCombo.SetItems(table.data_columns)
            selected = table.default_column(load_saved_column(self.parent))
            if selected:
                self.pseudoAxisCombo.SetStringSelection(selected)
        except PseudoAxisError as exc:
            self.pseudo_axis_table = None
            self.pseudoAxisCombo.SetItems([])
            self._set_status(str(exc))

    def OnPseudoAxisColumn(self, event):
        column = self.pseudoAxisCombo.GetValue().strip()
        if column:
            save_selected_column(self.parent, column)
            self._set_status("Pseudo-axis: %s" % column)

    def OnOpenPseudoAxisTable(self, event):
        path = getattr(self, 'pseudo_axis_file', pseudo_axis_path(self.parent, fallback_spec_dir=getattr(self, 'raw', None)))
        if not os.path.exists(path):
            self._set_status('Pseudo-axis table not found')
            wx.MessageBox('Cannot find %s' % path, 'Pseudo-axis table', wx.OK | wx.ICON_WARNING)
            return
        try:
            table = getattr(self, 'pseudo_axis_table', None) or PseudoAxisTable.load(path)
            self.pseudo_axis_table = table
            show_pseudo_axis_table(self, table)
            self._set_status('Viewing pseudo-axis table')
        except PseudoAxisError as exc:
            self._set_status(str(exc))
            wx.MessageBox(str(exc), 'Pseudo-axis table', wx.OK | wx.ICON_WARNING)

    def _pseudo_axis_values(self):
        table = getattr(self, 'pseudo_axis_table', None)
        if table is None:
            self._load_pseudo_axis_choices()
            table = getattr(self, 'pseudo_axis_table', None)
        column = self.pseudoAxisCombo.GetValue().strip()
        if table is None or not column:
            raise PseudoAxisError('Choose a pseudo-axis column before analysis')
        values = table.numeric_values(column)
        save_selected_column(self.parent, column)
        return numpy.asarray(values, dtype=float)

    def GetActualField(self):
        count = self.setLocal.GetItemCount()
        if(count==0):
            self.field=750
            self.fieldLab='750'
            self.nuc='15N'
            self.temp=20
        else:
            col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
            col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])
            self.field=float(col2[col1=='sfrq'][0])
            self.fieldLab=str(int(self.field))
            self.nuc=col2[col1=='nucleus'][0]
            self.temp=col2[col1=='temperature'][0]

        if(self.nuc not in self.Gamma.keys()):
            print('We only have stored the following Gyros')
            print(self.Gamma)
            print(self.nuc,'not found')
        
        self.dfrq=self.field/self.Gamma['1H']*self.Gamma[self.nuc]


    def DoFit(self,pk):
        self.ReadFuda(pk) #set cpmgX, cpmgY for residue

        #self.R0line=numpy.average(self.cpmgY)
        #x0=leastsq(self.ChiLine,self.PackLine())

        #self.chi2Line=numpy.average(self.ChiLine(self.PackLine())**2.)

        #self.Time_T2=float(self.TimeT2Box.GetValue())

        #self.GetActualField() #set self.dfrq

        self.A0=self.cpmgY[numpy.argmax(numpy.abs(self.cpmgY))]
        self.kex=1./(numpy.average(self.cpmgX))
        self._set_status('Fitting %s...' % pk)
        x0, cov_x, infodict, mesg, ier = leastsq(self.ChiDecay, self.PackDecay(), full_output=True)
        self.UnpackDecay(x0)
        self.DecayCalc()
        self._set_status('Fit complete: %s' % pk)

        self.chi2Local=numpy.average(self.ChiDecay(self.PackDecay())**2./self.cpmgE**2.)
        # Parameter standard errors from the least-squares covariance matrix.
        # leastsq returns covariance before residual-variance scaling.
        self.A0_err = numpy.nan
        self.kex_err = numpy.nan
        if cov_x is not None:
            try:
                residual = numpy.asarray(infodict.get('fvec', []), dtype=float)
                dof = max(1, residual.size - len(x0))
                scale = float(numpy.dot(residual, residual)) / dof
                errs = numpy.sqrt(numpy.maximum(0.0, numpy.diag(cov_x) * scale))
                self.A0_err, self.kex_err = float(errs[0]), float(errs[1])
            except Exception:
                pass
        

        #print(sele)
        count = self.datasets.GetItemCount()
        col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])

        #self.RexScreen=float(self.RexScreenBox.GetValue())

        for i,c1 in enumerate(col1):
            if(c1==pk):
                Rex=self.kex #numpy.max(self.cpmgF)-numpy.min(self.cpmgF)
                self.datasets.SetItem(i,2,'%.2f' % self.kex)
                self.datasets.SetItem(i,3,'%.2f' % self.A0)
                self.datasets.SetItem(i,4,'%.2f' % (self.chi2Local**0.5/numpy.abs(self.A0)*100.) )
                #self.datasets.SetItem(i,3,'%.2f' % (1.-self.chi2Local/self.chi2Line))
                #self.datasets.SetItem(i,4,'%.2f' % self.chi2Local)


                
                #if(Rex<self.RexScreen or Rex!=Rex):
                #    self.datasets.SetItem(i,1,str(False))
                #else:
                #    self.datasets.SetItem(i,1,str(True))


                
        if(pk not in self.cpmgLocal.keys()):
            self.cpmgLocal[pk]={}
        self.cpmgLocal[pk]['A0']=self.A0
        self.cpmgLocal[pk]['chi2local']=self.chi2Local
        self.cpmgLocal[pk]['kex']=self.kex
        self.cpmgLocal[pk]['A0_err']=self.A0_err
        self.cpmgLocal[pk]['kex_err']=self.kex_err
        self.cpmgLocal[pk]['average_error_pct']=(self.chi2Local**0.5/numpy.abs(self.A0)*100.)

        
    def prepare_workflow_analysis(self):
        """Prepare Decay results when opened via Workflow -> Show analysis."""
        self.OnButtonPeakConvert(None)
        self.OnButtonGuess(None)
        selection = self.plotModeRadio.FindString('R')
        if selection != wx.NOT_FOUND:
            self.plotModeRadio.SetSelection(selection)
        self.plotMode = 'R'
        self._plot_current_mode()
        self.Raise()

    def OnButtonGuess(self,event):
        #clever way to do this is to fit flat line, fit to my equation,
        #do FTest to check if fitting to dispersion led to something good.
        #loop over peaks

        self.Time_T2=float(self.TimeT2Box.GetValue())

        #self.GetActualField() #set self.dfrq

        
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)]
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            self.DoFit(c1)
        self._sort_dataset_rows_naturally()
        self._plot_current_mode()


        
                
            
        #print('Current:',col1,col2)
        #print('Changing inclusion')
        #num_items = self.datasets.GetItemCount()

        
        #1. read CPMG curve.
        #2. fit to baldwin eqn.
        #3. fit to flat line
        #compare chi2 values


        

    def OnButtonWipe(self,event):
        self._set_status('Clearing previous analysis output...')
        os.system('rm '+self.raw+'/catia/OutPut/*')

        
    def DoLocal(self,event):
        

        #self.seqfil_cpmg=self.seqfilCombo.GetValue().split()[0]
        #self.basis_cpmg=self.basisCombo.GetValue().split()[0]

        #self.field=self.fieldbox.GetValue().split()[0]
        #self.GetActualField() #set self.dfrq field and fieldlab
        #self.Time_T2=float(self.TimeT2Box.GetValue())
        #xcar=118


        
        self.local=[]
        self.sett=[]

        self.sett.append(('sfrq',self.field,'set'))
        self.sett.append(('temperature',20,'set'))
        self.sett.append(('nucleus','15N','set'))

        self.sett.append(('xcar',xcar,'set'))
        self.sett.append(('pwx_cp',30E-6,'set'))
        self.sett.append(('time_T2',self.Time_T2,'set'))            
        
        self.local.append(('Omega',xcar,'',''))
        self.local.append(('DeltaO',1.0,'fit',''))


        if(self.basis_cpmg=='IphAph_13'):
            self.local.append(('JIS',-92,'',''))            
            self.local.append(('DeltaJ',0,'',''))            
            self.local.append(('R0iph_%s' % self.fieldLab,5.0,'fit',''))            
            self.local.append(('R1aph_%s' % self.fieldLab,2.0,''))            
            self.local.append(('R1iph_%s' % self.fieldLab,1.0,''))            
        elif(self.basis_cpmg=='TrATr_13'):
            self.local.append(('JIS',-92,''))            
            self.local.append(('DeltaJ',0,''))            
            self.local.append(('R0Tr_%s' % self.fieldLab,5.0,'fit'))
            self.local.append(('R0ATr_%s' % self.fieldLab,5.0,''))            
            self.local.append(('R1aph_%s' % self.fieldLab,2.0,''))            
            self.local.append(('R1iph_%s' % self.fieldLab,1.0,''))            
        elif(self.basis_cpmg=='Iph_7'):
            self.local.append(('R0iph_%s' % self.fieldLab,5.0,'fit'))            
            self.local.append(('R1iph_%s' % self.fieldLab,1.0,''))            


        if(self.seqfil_cpmg=='PE_CPMG'):
            self.sett.append(('couplednucleus','1H','set'))
            self.sett.append(('taub',0.002,'set'))
            self.sett.append(('time_equil',0.0,'set'))
        elif(self.seqfil_cpmg=='Trosy_CPMG'):
            self.sett.append(('couplednucleus','1H','set'))
            self.sett.append(('taub',0.002,'set'))
            self.sett.append(('time_equil',0.0,'set'))
        elif(self.seqfil_cpmg=='CW_CPMG'):
            self.sett.append(('time_equil',0.0,'set'))
            pass
        elif(self.seqfil_cpmg=='fast_cpmg'):
            pass
            
        self.SetLocal() #transfer pars to arrays

        
    def SetLocal(self):


        #self.setLocal.ClearAll()
        #self.setLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        #self.setLocal.InsertColumn(1, 'Set', width = 80,format=wx.LIST_FORMAT_CENTRE)

        #self.parLocal.ClearAll()        
        #self.parLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        #self.parLocal.InsertColumn(1, 'Initial', width = 80,format=wx.LIST_FORMAT_CENTRE)
        #self.parLocal.InsertColumn(2, 'Fit?', width = 80,format=wx.LIST_FORMAT_CENTRE)
        #self.parLocal.InsertColumn(3, 'File', width = 80,format=wx.LIST_FORMAT_CENTRE)

        #num_items = self.parLocal.GetItemCount()
        #cnt=0
        #for local in self.local:
        #    cnt+=1
        #    self.parLocal.InsertStringItem(num_items,str(cnt))
        #    self.parLocal.SetItem(num_items,0,local[0])
        #    self.parLocal.SetItem(num_items,1,str(local[1]))
        #    self.parLocal.SetItem(num_items,2,str(local[2]))
        #    if(len(local)==4):
        #        self.parLocal.SetItem(num_items,3,str(local[3]))
            
        #num_items = self.setLocal.GetItemCount()
        #cnt=0
        #for local in self.sett:
        #    cnt+=1
        #    self.setLocal.InsertStringItem(num_items,str(cnt))
        #    self.setLocal.SetItem(num_items,0,local[0])
        #    self.setLocal.SetItem(num_items,1,str(local[1]))
            
            
        return


        
    def PackOpts(self,a,b):
        hbox=wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(a)
        hbox.Add(b)
        self.vboxOpts.Add(hbox)

    def OnClickSet(self,event):

        sele=self.datasets.GetFirstSelected()
        #print(sele)
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)][sele]
        print('Current:',col1,col2)
        print('Changing inclusion')
        num_items = self.datasets.GetItemCount()


        pass

    def OnClickPar(self,event):
        pass
        
    def _peak_sort_key(self, name):
        """Natural peak order: largest integer in the name, then full name."""
        nums = [int(x) for x in re.findall(r'\d+', str(name))]
        return (max(nums) if nums else float('inf'), str(name).lower())

    def _selected_peak(self):
        row = self.datasets.GetFirstSelected()
        if row < 0 or row >= self.datasets.GetItemCount():
            return None
        return self.datasets.GetItem(row, 0).GetText()

    def _select_dataset_row(self, row):
        count = self.datasets.GetItemCount()
        if count == 0:
            return
        row = max(0, min(row, count - 1))
        current = self.datasets.GetFirstSelected()
        if current >= 0:
            self.datasets.SetItemState(current, 0, wx.LIST_STATE_SELECTED | wx.LIST_STATE_FOCUSED)
        self.datasets.SetItemState(row, wx.LIST_STATE_SELECTED | wx.LIST_STATE_FOCUSED,
                                   wx.LIST_STATE_SELECTED | wx.LIST_STATE_FOCUSED)
        self.datasets.EnsureVisible(row)

    def OnDatasetKeyDown(self, event):
        key = event.GetKeyCode()
        if key not in (wx.WXK_UP, wx.WXK_DOWN):
            event.Skip()
            return
        count = self.datasets.GetItemCount()
        if not count:
            return
        row = self.datasets.GetFirstSelected()
        if row < 0:
            row = 0
        else:
            row += -1 if key == wx.WXK_UP else 1
        self._select_dataset_row(row)
        # Selection event normally redraws; explicit redraw also makes keyboard
        # behaviour reliable across wx versions.
        self._plot_current_mode()

    def _sort_dataset_rows_naturally(self):
        selected = self._selected_peak()
        rows = []
        for row in range(self.datasets.GetItemCount()):
            rows.append([self.datasets.GetItem(row, col).GetText() for col in range(5)])
        rows.sort(key=lambda vals: self._peak_sort_key(vals[0]))
        self.datasets.DeleteAllItems()
        selected_row = -1
        for i, vals in enumerate(rows):
            self.datasets.InsertItem(i, vals[0])
            for col, value in enumerate(vals[1:], 1):
                self.datasets.SetItem(i, col, value)
            if vals[0] == selected:
                selected_row = i
        if selected_row >= 0:
            self._select_dataset_row(selected_row)

    def OnPlotMode(self, event):
        self.plotMode = self.plotModeRadio.GetStringSelection()
        self._plot_current_mode()

    def _ensure_summary_fits(self, parameter):
        key = 'kex' if parameter == 'R' else 'A0'
        errkey = key + '_err'
        peaks = [self.datasets.GetItem(i, 0).GetText() for i in range(self.datasets.GetItemCount())]
        missing = [p for p in peaks if p not in self.cpmgLocal or
                   key not in self.cpmgLocal[p] or errkey not in self.cpmgLocal[p]]
        if missing:
            self.OnButtonGuess(True)
        return all(p in self.cpmgLocal and key in self.cpmgLocal[p] for p in peaks)

    def _plot_parameter_summary(self, parameter):
        if not self._ensure_summary_fits(parameter):
            return
        key = 'kex' if parameter == 'R' else 'A0'
        errkey = key + '_err'
        peaks = sorted([self.datasets.GetItem(i, 0).GetText()
                        for i in range(self.datasets.GetItemCount())], key=self._peak_sort_key)
        values = [self.cpmgLocal[p][key] for p in peaks]
        errors = [self.cpmgLocal[p].get(errkey, numpy.nan) for p in peaks]
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        x = numpy.arange(len(peaks))
        self.ax.bar(x, values, yerr=errors, capsize=3)
        self.ax.set_xticks(x)
        self.ax.set_xticklabels(peaks, rotation=60, ha='right')
        self.ax.set_xlabel('Peak')
        self.ax.set_ylabel('R' if parameter == 'R' else 'A0')
        self.ax.set_title('%s by peak' % parameter)
        self.fig.tight_layout()
        self.canvas.draw()

    def _plot_peak_decay(self):
        pk = self._selected_peak()
        if not pk:
            return
        if self.ReadFuda(pk) is False:
            return
        if pk not in self.cpmgLocal:
            self.DoFit(pk)
        else:
            self.A0 = self.cpmgLocal[pk]['A0']
            self.kex = self.cpmgLocal[pk]['kex']
            self.DecayCalc()
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.errorbar(self.cpmgX, self.cpmgY / self.A0,
                         yerr=(self.cpmgE / abs(self.A0)), fmt='o', label='data')
        self.ax.plot(self.cpmgX, self.cpmgF / self.A0, label='fit')
        self.ax.legend(loc='upper right')
        self.ax.set_xlabel('t')
        self.ax.set_ylabel('I/I_0')
        self.ax.set_title(pk)
        self.fig.tight_layout()
        self.canvas.draw()

    @staticmethod
    def _format_value_error(value, error):
        """Format an estimate/error pair with a two-significant-figure error."""
        try:
            value, error = float(value), abs(float(error))
            if not numpy.isfinite(value) or not numpy.isfinite(error) or error == 0:
                return ('%.6g' % value, '%.2g' % error)
            exponent = int(numpy.floor(numpy.log10(error)))
            decimals = max(0, 1 - exponent)
            return (('%.*f' % (decimals, value)), ('%.*f' % (decimals, error)))
        except (TypeError, ValueError):
            return (str(value) if value != '' else '', str(error) if error != '' else '')

    def report_results_rows(self):
        """Return fitted parameters using the precision of their errors."""
        self.OnButtonGuess(None)
        rows = []
        peaks = sorted([self.datasets.GetItem(i, 0).GetText()
                        for i in range(self.datasets.GetItemCount())], key=self._peak_sort_key)
        for peak in peaks:
            result = self.cpmgLocal.get(peak, {})
            rval, rerr = self._format_value_error(result.get('kex', ''), result.get('kex_err', ''))
            aval, aerr = self._format_value_error(result.get('A0', ''), result.get('A0_err', ''))
            avg = result.get('average_error_pct', '')
            if avg == '' and result.get('chi2local', '') != '' and result.get('A0', ''):
                try: avg = float(result['chi2local'])**0.5 / abs(float(result['A0'])) * 100.
                except Exception: avg = ''
            try: avg = '%.2f' % float(avg)
            except (TypeError, ValueError): pass
            rows.append([peak, rval, rerr, aval, aerr, avg])
        return ['Peak', 'R', 'R error', 'A0', 'A0 error', 'Average error (%)'], rows

    def export_report_figures(self, report_dir):
        """Fit all peaks and export structured summary and square peak figures."""
        import os
        os.makedirs(os.fspath(report_dir), exist_ok=True)
        columns, rows = self.report_results_rows()
        summary_figures = []
        old_size = tuple(self.fig.get_size_inches())
        for parameter, filename in (('R', 'decay_R.pdf'), ('A0', 'decay_A0.pdf')):
            self._plot_parameter_summary(parameter)
            self.fig.savefig(os.path.join(os.fspath(report_dir), filename), bbox_inches='tight')
            summary_figures.append((filename, '%s by peak' % parameter))
        peak_figures = {}
        selected = self._selected_peak()
        try:
            self.fig.set_size_inches(4.0, 4.0, forward=True)
            for i in range(self.datasets.GetItemCount()):
                peak = self.datasets.GetItem(i, 0).GetText()
                self._select_dataset_row(i)
                self._plot_peak_decay()
                self.fig.set_size_inches(4.0, 4.0, forward=True)
                filename = 'decay_peak_%03d.pdf' % (i + 1)
                self.fig.savefig(os.path.join(os.fspath(report_dir), filename), bbox_inches='tight')
                peak_figures[peak] = filename
        finally:
            self.fig.set_size_inches(*old_size, forward=True)
        if selected:
            for i in range(self.datasets.GetItemCount()):
                if self.datasets.GetItem(i, 0).GetText() == selected:
                    self._select_dataset_row(i); break
        return {'columns': columns, 'rows': rows, 'summary_figures': summary_figures,
                'peak_figures': peak_figures}

    def _plot_current_mode(self):
        mode = getattr(self, 'plotMode', 'Peak')
        if mode == 'Peak':
            self._plot_peak_decay()
        else:
            self._plot_parameter_summary(mode)

    def OnClick(self,event): #when selecting a peak...plot what we have.

        sele=self.datasets.GetFirstSelected()
        #print(sele)
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)][sele]
        print('Current:',col1,col2)
        print('Changing inclusion')
        num_items = self.datasets.GetItemCount()

        if(col2=='True'):
            self.datasets.SetItem(sele,1,str(False))
        else:
            self.datasets.SetItem(sele,1,str(True))
            
        pass
        
    def OnButtonRefresh(self,event):

        self.datasets.ClearAll()

        self.datasets.InsertColumn(0, 'Peak', width = 45,format=wx.LIST_FORMAT_CENTRE) 
        self.datasets.InsertColumn(1, 'Include', width = 55,format=wx.LIST_FORMAT_CENTRE)
        self.datasets.InsertColumn(2, 'R', width = 55,format=wx.LIST_FORMAT_CENTRE)
        self.datasets.InsertColumn(3, 'A0', width = 55,format=wx.LIST_FORMAT_CENTRE)
        self.datasets.InsertColumn(4, 'AveError(%s)' % ('%',), width = 65,format=wx.LIST_FORMAT_CENTRE)
        #self.datasets.InsertColumn(4, 'Chi2local', width = 80,format=wx.LIST_FORMAT_CENTRE)

        files=os.listdir(self.fuda_dir)

        self.peaks=[]
        for file in files:
            if(file[-4:]=='.out'):
                self.peaks.append(file.split('.out')[0])

        self.peaks=sorted(self.peaks, key=self._peak_sort_key)
        cnt=0

        for peak in self.peaks:
            row = self.datasets.GetItemCount()
            cnt += 1
            self.datasets.InsertItem(row, peak)
            self.datasets.SetItem(row, 0, peak)
            self.datasets.SetItem(row, 1, str(True))
        self._sort_dataset_rows_naturally()
        if self.datasets.GetItemCount() and self.datasets.GetFirstSelected() < 0:
            self._select_dataset_row(0)
        pass

    def OnButtonSort(self,event):

        col=event.GetColumn()
        

        count = self.datasets.GetItemCount()
        col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.datasets.GetItem(row, 1).GetText() for row in range(count)])
        try:
            col3 = numpy.array([float(self.datasets.GetItem(row, 2).GetText()) for row in range(count)])
        except:
            col3 = numpy.array([self.datasets.GetItem(row, 2).GetText() for row in range(count)])
        try:
            col4 = numpy.array([float(self.datasets.GetItem(row, 3).GetText()) for row in range(count)])
        except:
            col4 = numpy.array([self.datasets.GetItem(row, 3).GetText() for row in range(count)])
        #print(numpy.argsort(col2))
        
        self.datasets.ClearAll()
        self.datasets.InsertColumn(0, 'Peak', width = 45,format=wx.LIST_FORMAT_CENTRE) 
        self.datasets.InsertColumn(1, 'Include', width = 55,format=wx.LIST_FORMAT_CENTRE)
        self.datasets.InsertColumn(2, 'Rex', width = 60,format=wx.LIST_FORMAT_CENTRE)
        self.datasets.InsertColumn(3, 'Chi2diff', width = 70,format=wx.LIST_FORMAT_CENTRE)

        
        if(col==0):
            s=numpy.argsort(col1)
        elif(col==1):
            s=numpy.flip(numpy.argsort(col2))
        elif(col==2):
            s=numpy.flip(numpy.argsort(col3))
        elif(col==3):
            s=numpy.flip(numpy.argsort(col4))
            
            
        for i,arg in enumerate(s):
            #print (arg,col1[arg],col2[arg],col3[arg],col4[arg])
            self.datasets.InsertStringItem(i,i)
            self.datasets.SetItem(i,0,str(col1[arg]))
            self.datasets.SetItem(i,1,str(col2[arg]))
            self.datasets.SetItem(i,2,str(col3[arg]))
            self.datasets.SetItem(i,3,str(col4[arg]))

            

        
        
    def OnButtonAddDataset(self,event):
        pass
    def OnButtonRemDataset(self,event):
        pass

    def PathExists(self,test):
        for t in test:
            if(os.path.exists(t)==False):
                os.system('mkdir '+t)

    
    def OnButtonRun(self,event):

        self.PathExists((self.raw+'/catia',self.raw+'/catia/dataset',self.raw+'/catia/OutPut'))


        self.datfile=[]
        self.globfile=[]
        self.locfile=[]

        self.KexFit=True
        self.PbFit=True
        self.Conv=1E-3
        self.MaxIter=100
        
        self.WriteCatiaDataset()
        self.WriteCatiaPar()
        self.WriteCatiaFile()

        
        if(os.uname()[0]=='Darwin'): #run catia
            rc = os.system('catia_Darwin_i386 < '+self.raw+'/catia/OneFieldFit.catia')
        else:
            rc = os.system('catia_Linux_x86_64 < '+self.raw+'/catia/OneFieldFit.catia')
        if rc == 0:
            mark_complete = getattr(self.parent, '_mark_pseudo_analysis_complete', None)
            if callable(mark_complete):
                mark_complete(model='decay', source='catia')


    def OnButtonPeakConvert(self,event):
        count = self.datasets.GetItemCount()
        peaks = [self.datasets.GetItem(row, 0).GetText() for row in range(count)]
        self._set_status('Creating decay curves for %d peaks...' % len(peaks))

        for index, pk in enumerate(peaks, 1):
            self._set_status('Creating decay curves: %d/%d (%s)' % (index, len(peaks), pk))
            cpmgfile=os.path.join(self.fuda_dir, pk + '.out.decay')
            if(os.path.exists(cpmgfile)):
                os.system('rm '+cpmgfile)
            self.ReadFuda(pk)
        self._set_status('Decay curves ready for %d peaks' % len(peaks))
        

    def WriteCatiaDataset(self):


        self.datfile.append(self.raw+'/catia/dataset.inp')
        outy=open(self.raw+'/catia/dataset.inp','w')

        self.GetActualField() #set dfrq and field, temp and fieldlab
        outy.write('ID=15N CPMG Trosy @ %s and %s\n' % (self.fieldLab,self.temp))
        
        count = self.setLocal.GetItemCount()
        col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])
        for c1,c2 in zip(col1,col2):
            outy.write(c1+' = '+c2+'\n')
        outy.write('minerror = (1.%;0.3/s)\n')
        outy.write('seqfil = %s\n' % self.seqfilCombo.GetValue())
        outy.write('basis = (%s)\n' % self.basisCombo.GetValue())

        outy.write('format = (0;1;2)\n')
        outy.write('DataDirectory = %s \n' % (self.fuda_dir + os.sep)   )
        outy.write('Data = (\n')
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)]
        for c1,c2 in zip(col1,col2):
            if(c2=='True'):
                outy.write('[%s;%s];\n' % (c1,c1+'.out.cpmg'))
        outy.write(')\n')
        outy.close()



    def WriteCatiaPar(self):
        #each dataset requires its own local parameters

        #sele=self.setLocal.GetFirstSelected()
        #print(sele)
        self.GetActualField() #set dfrq,field and fieldlab
        
        #count = self.setLocal.GetItemCount()
        #col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
        #col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])

        
        #sele=self.parLocal.GetFirstSelected()
        #print(sele)
        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])

        outy=open(self.raw+'/catia/ParamSet_'+self.fieldLab+'_'+self.temp+'.inp','w')
        self.locfile.append(self.raw+'/catia/ParamSet_'+self.fieldLab+'_'+self.temp+'.inp')
        outy.write('format = (')
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % c1)
        outy.write(')\n')
        outy.write('* = (')
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            if(i!=0):
                outy.write(';')
            outy.write('%s' % str(c2))        
        outy.write(')\n')
        outy.close()

        self.fix=[]
        for i,(c1,c3) in enumerate(zip(col1,col3)):
            if(c3=='fit'):
                pass
            else:
                self.fix.append(c1)
            
        self.globfile.append(self.raw+'/catia/ParamGlobal.inp')
        outy=open(self.raw+'/catia/ParamGlobal.inp','w')
        outy.write('kex=1000.\n')
        outy.write('pb=0.02\n')
        outy.close()

        
    def WriteCatiaFile(self):
        outy=open(self.raw+'/catia/OneFieldFit.catia','w')
        for i in range(len(self.datfile)):
            outy.write('ReadDataset(%s)  #Data summary file \n' % self.datfile[i])
        for i in range(len(self.locfile)):
            outy.write('ReadParam_Local(%s)      #Local parameters initial\n' % self.locfile[i])
        for i in range(len(self.globfile)):
            outy.write('ReadParam_Global(%s)   #Global parameters initial\n' % self.globfile[i])


        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])
        col4 = numpy.array([self.parLocal.GetItem(row, 3).GetText() for row in range(count)])


        for i,(c1,c4) in enumerate(zip(col1,col4)):
            if(len(c4)>0):
                fil=c4.split()[0]
                if(os.path.exists(fil)):
                    outy.write('ReadParam(%s;%s;0;1)\n' % (c1,c4))
            
        #outy.write('ReadParam(Omega;./raw/test.ft2.list;0;1)   #peak list\n')
        #if(os.path.exists('./raw/catia/DeltaOmega.inp')):
        #    outy.write('ReadParam(DeltaO;./raw/catia/DeltaOmega.inp;0;1)#delta omegas\n')

        outy.write('# Fix all the static parameters\n')

        
        for i in range(len(self.fix)):
            outy.write('FreeLocalParam(all;%s;false)\n' % self.fix[i])
            

        #    outy.write('SetGlobalParam(kex;250)\n')
        #    outy.write('SetGlobalParam(pb;0.05)\n')

        self.PathExists((self.raw+'/catia/OutPut',))
        outy.write('# Deal with global parameters\n')

        if(self.KexFit):
            outy.write('FreeGlobalParam(kex;true)\n')
        else:
            outy.write('FreeGlobalParam(kex;false)\n')
        if(self.PbFit):
            outy.write('FreeGlobalParam(pb;true)\n')
        else:
            outy.write('FreeGlobalParam(pb;false)\n')
        outy.write('# Minimize\n')
        outy.write('echo(\\n)\n')
        outy.write('#Minimize()\n')
        outy.write('Minimize(print=y;tol=%f;maxiter=%i)\n' % (self.Conv,self.MaxIter))
        outy.write('#\n')
        outy.write('#  // Print some files for plotting\n')
        outy.write('PrintParam('+self.raw+'/catia/OutPut/GlobalParam.fit;global)\n')
        outy.write('PrintParam('+self.raw+'/catia/OutPut/DeltaOmega.fit;DeltaO)\n')
        outy.write('PrintData('+self.raw+'/catia/OutPut/)\n')
        outy.write('echo(\n)\n')
        outy.write('ChiSq(all;all)\n')
        outy.write('exit(0)\n')
        outy.close()


        
    
    
    #if we have intensity plot and NOE tab selected,
    #then a click on the intensity graph selects an NOE.
    def on_pick(self, event):
        seleRbox=self.rbox.GetSelection()
        if(seleRbox==3):
            print('Intensity')
            print(event.xdata,event.ydata)
            x_min,x_max=self.ax.get_xlim()
            y_min,y_max=self.ax.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min        
            rad2=((self.inty[:,0]-event.xdata)/xdist)**2.+((self.inty[:,2]-event.ydata)/ydist)**2.

            argmin=numpy.argmin(rad2)
            print('argmin:',argmin)
            print(self.sheetbox.GetStringSelection())
            if(self.sheetbox.GetStringSelection()=='NOEs'):            

                sele=self.source.GetFirstSelected()
                self.source.Select(sele,on=0)
                self.source.Select(argmin,on=1)
                self.source.EnsureVisible(argmin) 
                self.OnTickPlot(True)
                #self.subgraph.SetSelection(argmin)            

            #self.xv = [self.source.GetItem(row, 7).GetText() for row in xrange(count)][sele]
            #self.yv = [self.source.GetItem(row, 3).GetText() for row in xrange(count)][sele]
                    #print xv,yv
            #self.ax.scatter(float(self.xv),float(self.yv),color='r',s=200)
    

    def onShiftXBtn(self,event):
        if 'shiftXFile' not in list(self.parent.pars.keys()):
        # if not self.parent.pars['shiftXFile']:
            print('Could not compare with calculated NMR shifts as the file is not there.')
            self.sheetbox.SetStringSelection('Assignments')
            self.OnTickFilt(True)
            return
        if 'ileShiftFile' not in list(self.parent.pars.keys()):
        # if not self.parent.pars['ileShiftFile']:
            print('Could not compare with calculated NMR shifts as the experimental NMR shift file is not there.')
            self.sheetbox.SetStringSelection('Assignments')
            self.OnTickFilt(True)
            return

        self.doShiftX()
        self.PopulateList(shift='y')

        """
        #populate list
        resultConf={}
        conf=0
        for key,vals in resDict.items():
            if(len(vals)==1):
                conf+=1
                resultConf[key]=vals
                del(resDict[key])
        cnt=0
        for i in range(self.colVal+1):
            if(os.path.exists(self.progress[i][1]+'/confident.res')): #is this mode done?
                new={}
                inny=open(self.progress[i][1]+'/confident.res')
                for line in inny.readlines():
                    test=line.split(':')
                    key=test[0].split()[0]
                    ass=test[1].split()[0]
                    if(key in resultConf.keys()):
                        if(ass==resultConf[key][0]):
                            new[key]=ass
                            del(resultConf[key])
                inny.close()

                for key,vals in new.items():
                    cnt+=1

                    num_items = self.source.GetItemCount()
                    self.source.InsertStringItem(num_items,str(cnt))
                    self.source.SetItem(num_items,0,str(cnt))
                    self.source.SetItem(num_items,1,key)
                    self.source.SetItem(num_items,2,self.progress[i][0])

                    #num_items = self.source.GetItemCount()
                    ##FGA changed- depreciated functions
                    #self.source.InsertItem(str(cnt))
                    #self.source.SetItem(0,str(cnt))
                    #self.source.SetItem(1,key)
                    #self.source.SetItem(2,self.progress[i][0])


                    stry=vals
                    self.source.SetItem(num_items,3,stry)

                    if(self.progress[i][0]=='subgraphMode' or self.progress[i][0]=='polishMode' or self.progress[i][0]=='nudgeMode' or self.progress[i][0]=='finalMode'):

                        color = (0,int(255), 0)
                        self.source.SetItemBackgroundColour(num_items,color)

        for key in sorted(resultConf,key=lambda k:len(resultConf[k])):
            val = resultConf[key][0]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(cnt))
            self.source.SetItem(num_items,1,key)
            self.source.SetItem(num_items,2,'compareShiftsFilter')
            self.source.SetItem(num_items,3,val)

        for key in sorted(resDict,key=lambda k:len(resDict[k])):
            vals=resDict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(cnt))
            self.source.SetItem(num_items,1,key)
            stry=''
            for val in vals:
                stry+=val+' '
            self.source.SetItem(num_items,3,stry)

        print 'Comparison with calculated shifts done.'
        """
    """
    def doShiftX(self):
        print('shiftXPostMode: Ruling out options from combined results if we can.')

        inst = self.parent.inst
        #Magma('input.magma',run='n')
        prior = str(self.colVal+1)
        if os.path.exists(inst.P.outdir+'/'+prior):
            shiftX = shiftXNMR(self.parent.pars['shiftXFile'].strip(), self.parent.pars['ileShiftFile'].strip(),
                                                self.parent.pars['chains'], self.parent.pars['residues'],
                                                inst.P.outdir+'/'+prior+'/combinedResults.res')
        else:
            shiftX = shiftXNMR(self.parent.pars['shiftXFile'].strip(), self.parent.pars['ileShiftFile'].strip(),
                                                self.parent.pars['chains'], self.parent.pars['residues'],
                                                inst.P.outdir+'/combinedResults.res')
        self.shiftDict=shiftX.Parse() #return assignment dictionary 







    def SetSubgraph(self):
        #self.subgraph=wx.ComboBox(self, -1, size=(80, -1), choices=listy, style=wx.CB_READONLY)
        #self.Bind(wx.EVT_COMBOBOX, self.OnTickPlot, self.subgraph)
        listy=[]
        listy.append('all')
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            for key,vals in list(self.parent.inst.subgraphRef.items()):
                listy.append(str(key+1))
        if(self.ref=='nudgeMode'):
            self.grps=self.DoNudge()
            for i in range(len(self.grps)):
                listy.append(str(i))

        self.subgraph.SetItems(listy)

        #FGA changed
        #if self.ref=='nudgeMode':
        #    self.subgraph.SetSelection(1)
        #else:
        #    self.subgraph.SetSelection(0)
        self.subgraph.SetSelection(0)


    def ParseMagma(self):
        inny=open(self.parent.inst.P.outdir+'/input.magma')
        check=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(test[0]=='subgraphMode' and test[1]=="on"):
                    check.append(0)
                if(test[0]=='polishMode' and test[1]=="on"):
                    check.append(1)
                if(test[0]=='nudgeMode' and test[1]=="on"):
                    check.append(2)
                if(test[0]=='finalMode' and test[1]=="on"):
                    check.append(3)
                if(test[0]=='longMode' and test[1]=="on"):
                    check.append(4)
        if(self.WXV==4):
            self.parent.protoBox.SetCheckedItems(check)
        else:
            self.parent.protoBox.SetChecked(check)


    def GetModes(self,curr='n'):
        if(curr=='n'):
            self.colVal=self.modes.GetCurrentSelection()

        self.ParseMagma()
        self.progress=[]

        cnt=0
        if(self.parent.protoBox.IsChecked(0)):
            cnt+=1
            self.progress.append(('subgraphMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(1)):
            cnt+=1
            self.progress.append(('polishMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(2)):
            cnt+=1
            self.progress.append(('nudgeMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(3) or self.parent.protoBox.IsChecked(0)==False):
            cnt+=1
            self.progress.append(('finalMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(4)):
            cnt+=1
            self.progress.append(('longMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))


        self.progress=numpy.array(self.progress)
        self.progress[-1][1]=self.parent.inst.P.outdir #adjust last folder

        cnt=0
        for i in range(len(self.progress)-1):
            if(os.path.exists(self.progress[i][1])):
                self.progress[i,3]='complete'
                cnt+=1
        if(cnt==len(self.progress)-1):
            if(os.path.exists(self.progress[-1,1]+'/combinedResults.res')):
                self.progress[-1,3]='complete'

        #print self.progress
        #print len(self.progress)
        tmp=[]
        tag=0
        for i in range(len(self.progress)):
            if(self.progress[i,3]=='incomplete'):
                tmp.append(self.progress[i,:])
                break
            else:
                tmp.append(self.progress[i,:])
        self.progress=numpy.array(tmp)
        print('after',len(self.progress))
        print(self.progress)
        self.modes.SetItems(self.progress[:,0])

        if(curr=='y'):
            self.colVal=len(self.progress)-1
            while(1==1):
                if(self.colVal>0):
                    if(self.progress[self.colVal-1,3]=='incomplete'):
                        self.colVal-=1
                    else:
                        break
                else:
                    break

        self.modes.SetSelection(self.colVal)        
        self.ref = self.progress[self.colVal][0]
        self.testdir=self.progress[self.colVal][1]
        if(os.path.exists(self.testdir)==0):
            self.testdir=self.parent.inst.P.outdir
            print('Calculation not finished.')
            print('Current mode is not saved')
             

    def ParseLogAll(self,var,full='n',split='y'):
        inny=open(self.testdir+'/log')
        listy=[]
        tag=0
        for line in inny.readlines():
            test=line.split(var)
            if(len(test)>1):
                tag=1
                if(full=='n'):
                    if(split=='y'):
                        return line.split(var)[1].split()
                    else:
                        return line.split(var)[1]
                else:
                    if(split=='y'):
                        listy.append(line.split(var)[1].split())
                    else:
                        listy.append(line.split(var)[1])
        if(tag==0):
            print('Could not find',var,'in')
            return 0
        else:
            return listy   

    def DoNudge(self):
        inny=open(self.testdir+'/log')
        iny=0
        grp=0
        tbl=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)!=0):
                if(test[0]=='INTERSECT:'):
                    iny=1
                if(test[0]=='Groupings:'):
                    grp=1
                if(iny==1 and grp==0 and test[0]!='INTERSECT:'):
                    tbl.append(test)
        inny.close()
        #st=''
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    st+='c'
        #outlat.write('\\begin{tabular}{c|%s}\n' % (st)) 
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    outlat.write('& %i' % (ii+1))                        
        #outlat.write('\\\\ \n')
        #outlat.write('\\hline\n')                        
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    outlat.write(' %i' % (ii+1))                        
        #    for jj in range(len(self.inst.subgraphRef.keys())):
        #        outlat.write('& %s' % (tbl[ii][jj+1].split(':')[1]))                        
        #    outlat.write('\\\\ \n')
        #outlat.write('\\end{tabular}\n') 
        #outlat.write('\n\n')                    
        #outlat.write('Grouping subgraphs with overlapped assignments:\n\n ')

        groups=self.ParseLogAll("Groupings:",split='n')
        test= re.findall(r'[0-9]+',groups)

        grps=[]
        cnt=0
        proggy=''
        for ii,t in enumerate(test):
            proggy+=str(int(t)+1)+' '
            if(ii%self.parent.inst.P.nudgy==self.parent.inst.P.nudgy-1):
                grps.append((cnt,proggy))
                proggy=''
                cnt+=1
        return grps


    def WhereAmI(self): #figure out where in the calc we are.
        #self.testdir=self.progress[self.colVal][1]
        print('working out calculation stage...')
        self.status=0

        #if save directory does not exist, calc is incomplete
        #if combinedResults.res does not exist, calc is incomplete
        for i in range(5):
            self.statusList[i].SetLabel("")

        self.statusList[0].SetLabel('mode '+self.progress[self.colVal,3])

        stt=2
        #check for existance of progress files
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            wait=[]
            win=[]
            for ref in list(self.parent.inst.subgraphRef.keys()):
                vf2file=self.testdir+'/vf2_%i.txt' % (int(ref)+1)
                mcsfile=self.testdir+'/mces_%i.txt' % (int(ref)+1)
                tag=0
                if(os.path.exists(vf2file)):
                    win.append(str(int(ref)+1)+'vf2' )
                    tag=1
                if(os.path.exists(mcsfile)):
                    win.append(str(int(ref)+1)+'mcs' )
                    tag=1
                if(tag==0):
                    wait.append(str(int(ref)+1))
                #else:

            strg=''
            strg='Finished: '
            for i in range(len(win)):
                strg+=str(win[i])+' '
                if(len(strg)>25):
                    self.statusList[stt].SetLabel(strg)
                    strg=''
                    stt+=1
            self.statusList[stt].SetLabel(strg);stt+=1

            if(len(wait)>0):
                strg=''
                strg='Remain: '
                for i in range(len(wait)):
                    strg+=str(wait[i])+' '
                self.statusList[stt].SetLabel(strg)

            if(len(wait)>0):
                print('currently running ',wait[0])
                curr=' on '+str(wait[0])
                if(len(wait)>1):

                    self.status=1
                    print('to go...')
                    for i,wa in enumerate(wait):
                        if(i!=0):
                            print(wa)
            else:
                curr=''
        elif(self.ref=='nudgeMode'):
            self.grps=self.DoNudge()
            stt=2
            strg='groups: '
            for i in range(len(self.grps)):
                strg+='('+self.grps[i][1]+') '
                if(len(strg)>25):
                    #print stt
                    self.statusList[stt].SetLabel(strg)
                    strg=''
                    stt+=1
            self.statusList[stt].SetLabel(strg)
        else:
            vf2file=self.testdir+'/vf2.txt' 
            mcsfile=self.testdir+'/mces.txt'
            tag=0
            if(os.path.exists(vf2file)):
                tag=1
            if(os.path.exists(mcsfile)):
                tag=1
            #if(tag==0):
            curr=''

        self.GetLastFile()
        #print 'fil:',self.fil        
        if(self.fil==''):
            return
            
        if(len(self.fil.split('community'))>1):
            self.status=2
            self.statusList[1].SetLabel('running communityMode %s' % curr)
        if(len(self.fil.split('_new'))>1):
            self.status=3
            self.statusList[1].SetLabel('running splitMode %s' % curr)
            return
        if(len(self.fil.split('vf2'))>1):
            self.status=4
            self.statusList[1].SetLabel('running vf2')
            return
        if(len(self.fil.split('.G'))>1):
            self.status=5
            self.statusList[1].SetLabel('running mces %s' % curr)
            return 
        if(len(self.fil.split('conv'))>1):
            self.status=6
            self.statusList[1].SetLabel('optimising mces %s' % curr)
            return

    def GetLastFile(self):
        list_of_files = glob.glob(self.testdir+'/*') # * means all if need specific format then *.csv
        adjusted= sorted(list_of_files,key=os.path.getctime)
        #print adjusted
        #latest_file = max(list_of_files, key=os.path.getctime)
        #print latest_file
        cnt=0
        while(1==1):
            ii=len(adjusted)-1-cnt
            test=adjusted[ii].split('/')[-1]
            #print test
            if(test!='mces_traj.png' and test!='mcesplt.gp' and test!='log'):
                self.fil=test
                return
            cnt+=1
            if(cnt==len(adjusted)):
                self.fil=''
                return

    def AtoI(self,val):
        for i in range(len(self.parent.parent.tabOne.peak)):
            if(val==self.parent.parent.tabOne.peak[i].name):
                return i

    def OnDoubleClick(self,event):
        self.OnTickPlot(True)
        print(self.sheetbox.GetStringSelection())
        if(self.sheetbox.GetStringSelection()=='NOEs'):
                sele=self.source.GetFirstSelected()
                print(sele)
                count = self.source.GetItemCount()
                col1 = [self.source.GetItem(row, 1).GetText() for row in range(count)][sele]
                col2 = [self.source.GetItem(row, 2).GetText() for row in range(count)][sele]
                print(col1,col2)

                t1=re.findall(r'[0-9]+',col1)[0]
                t2=re.findall(r'[0-9]+',col2)[0]
                inny=open('out/correlate.3')
                for line in inny.readlines():
                    test=line.split()
                    if(len(test)>0):
                        A=re.findall(r'[0-9]+',test[0])[0]
                        B=re.findall(r'[0-9]+',test[1])[0]
                        if(A==t1 and B==t2):
                            p1=test[0]
                            p2=test[1]
                inny.close()
                print(col1,col2,p1,p2,self.AtoI(p1),self.AtoI(p2))
                print(self.parent.parent.tabOne.peak[self.AtoI(p1)])
                self.parent.parent.tabFour.ComboBox1.SetSelection(self.AtoI(p1))
                self.parent.parent.tabFour.ComboBox2.SetSelection(self.AtoI(p2))
                self.parent.parent.tabFour.draw_figure()

    """

    def OnButtonClose(self,event):
        self.Close()

    def OnButtonLog(self,event):
        textEdit.MyFrame(self.testdir+'/log')


    def OnButtonLoad(self,event):
        self.sett = []
        self.local = []

        if os.path.exists(self.savefile) == False:
            print('Cannot find savefile:', self.savefile)
            return

        data = read_structured_parameter_file(self.savefile)

        for key, value in data.get('set', {}).items():
            self.sett.append((key, value))

        if 'Tdelay' in data:
            # Legacy per-window save files may contain Tdelay.  The project-level
            # decayTimeMult is authoritative when present.
            system_mult = parse_value(self._project_parameter_path(), 'decayTimeMult', default='')
            self.TimeT2Box.SetValue(str(system_mult or data['Tdelay']))

        for peak_id, value in data.get('peak', {}).items():
            count = self.datasets.GetItemCount()
            col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
            if peak_id not in col1:
                self.datasets.InsertItem(count, count)
                self.datasets.SetItem(count, 0, peak_id)
            count = self.datasets.GetItemCount()
            col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
            for i, c1 in enumerate(col1):
                if c1 == peak_id:
                    if value == 'True':
                        self.datasets.SetItem(i, 1, str(True))
                    elif value == 'False':
                        self.datasets.SetItem(i, 1, str(False))

        for par_name, par_vals in data.get('par', {}).items():
            val = []
            val.append(par_name)
            val.append(par_vals.get('1', ''))
            if par_vals.get('2', '') == 'fit':
                val.append('fit')
            else:
                val.append('')

            if par_vals.get('3', '') == 'None':
                val.append('')
            else:
                val.append(par_vals.get('3', ''))

            print('adding:', val)
            self.local.append(val)
            print(val)

        self._sort_dataset_rows_naturally()

        self.SetLocal() #transfer pars to arrays

    def OnButtonSave(self,event):

        write={}
        
        write['Tdelay']=self.TimeT2Box.GetValue()
        self._save_time_mult()
        #write['seqfil']=self.seqfilCombo.GetValue()
        #write['basis']=self.basisCombo.GetValue()
        #write['RexScreen']=self.RexScreenBox.GetValue()
        
        count = self.datasets.GetItemCount()
        col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.datasets.GetItem(row, 1).GetText() for row in range(count)])
        write['peak']={}
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            write['peak'][c1]=c2
        
        #count = self.setLocal.GetItemCount()
        #col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
        #col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])
        #write['set']={}
        #for i,(c1,c2) in enumerate(zip(col1,col2)):
        #    write['set'][c1]=c2
        

        #count = self.parLocal.GetItemCount()
        #col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        #col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        #col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])
        #col4 = numpy.array([self.parLocal.GetItem(row, 3).GetText() for row in range(count)])

        #write['par']={}
        #for i,(c1,c2,c3,c4) in enumerate(zip(col1,col2,col3,col4)):
        #    write['par'][c1]={}
        #    write['par'][c1][1]=c2

        #    if(c3==''):
        #        write['par'][c1][2]='fix'
        #    else:
        #        write['par'][c1][2]=c3

        #    if(c4==''):
        #        write['par'][c1][3]='None'
        #    else:
        #        write['par'][c1][3]=c4
        write_structured_parameter_file(self.savefile, write)
                
        pass
        
        """
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            defaultFile='results.res',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            self._set_status('Writing results: %s' % os.path.basename(path))
            outfileFile=path
            
            #if(os.path.exists(outfileFile)==0):
            #    outy=open(outfileFile,'w');outy.close()
            os.system('cp '+self.parent.inst.P.outdir+'/combinedResults.res '+outfileFile)

            #inst = self.parent.inst
            #Magma('input.magma',run='n')
            #prior = str(self.colVal+1)
            #if os.path.exists(inst.P.outdir+'/'+prior):

            #self.canvas.print_figure(path, dpi=self.dpi)
            #self.parent.parent.flash_status_message("Saved %s" % path)
        """

    def OnButtonShow(self,event):
        test=Visualise(self.parent.inst)
        
        tig=0
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1
        else:
            tig=1

        if(tig==1):
            test.monte_min()
            test.pymol_gen_noe()
        else:
            test.monte_min(sele=int(sele))
            test.pymol_gen_noe()

        #from os import popen
        self.parent.PymolExec('pyscript_noe.py')
        #self.calcy=popen('pymol report/pyscript_noe.py')


    def OnTickFilt(self,event):
        self.WhereAmI() #figure out where in the calc we are.
        #selebox=self.sheetbox.GetSelection()
        boxVal=self.sheetbox.GetStringSelection()
        #print selebox,self.sheetbox.GetStringSelection()
        self.source.ClearAll()
        if(boxVal=='Assignments'):        
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE) 

            self.source.InsertColumn(3, 'PDB',width=1000) 
            self.PopulateList()
            return
        elif(boxVal=='Filters'):        
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(3, 'PDB', width = 10000) 
            self.PopulateFilt()
            return
        elif(boxVal=='NOEs'):        
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(1, 'A', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(2, 'B', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(3, 'intensity', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(4, '#', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(5, 'tp', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(6, 'hits', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(7, 'dist', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(8, 'stdev', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.PopulateNOE()
            return
        elif boxVal == 'CompareShifts':
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE)
            #FGA changed
            # self.source.InsertColumn(3, 'PDB', width = 10)
            self.source.InsertColumn(3, 'PDB')
            self.onShiftXBtn()


    def OnRefresh(self,event):
        print('getting modes')
        self.GetModes()
        print('setting subgraph')
        self.SetSubgraph() #reset subgraph options
        self.parent.UpdateReport()
        self.sheetbox.SetSelection(0)
        self.OnTickPlotFilt(True)
        #FGA added
        #self.Layout()

    def OnTickPlotFilt(self,event):
        self.UpdateScores()
        self.WhereAmI()
        self.OnTickPlot(True)
        self.OnTickFilt(True)        

    def OnTickPlot(self,event):
        self.draw_figure()
        #FGA added
        #if self.ref=='subgraphMode':
        #    self.rbox.Disable()
        #else:
        #    self.rbox.Enable() 


    def _analysis_spec_dir(self, inherit, dataset_root=''):
        """Resolve analysis workspace as WorkingDir/SpecPath."""
        obj = inherit
        while obj is not None:
            state = getattr(obj, 'state', None)
            if state is not None and hasattr(state, 'spec_dir'):
                try:
                    return os.path.normpath(state.spec_dir())
                except Exception:
                    pass
            obj = getattr(obj, 'parent', None)

        # Compatibility fallback for callers without ProjectState.  A supplied
        # dataset root is the working directory; SpecPath defaults to ./spec.
        root = str(dataset_root or '.').strip() or '.'
        spec_path = './spec'
        obj = inherit
        while obj is not None:
            ctrl = getattr(obj, 'specPathBox', None)
            if ctrl is not None:
                try:
                    spec_path = str(ctrl.GetValue() or './spec').strip() or './spec'
                    break
                except Exception:
                    pass
            obj = getattr(obj, 'parent', None)
        if os.path.isabs(spec_path):
            return os.path.normpath(spec_path)
        return os.path.normpath(os.path.join(root, spec_path))

    def _acquisition_dir(self, inherit, dataset_root=''):
        """Resolve vendor/acquisition input independently of analysis output."""
        obj = inherit
        while obj is not None:
            state = getattr(obj, 'state', None)
            if state is not None and hasattr(state, 'raw_dir'):
                try:
                    return os.path.normpath(state.raw_dir())
                except Exception:
                    pass
            obj = getattr(obj, 'parent', None)
        root = str(dataset_root or '.').strip() or '.'
        return os.path.normpath(os.path.join(root, 'raw'))

    def create_main_panel(self):
        
        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        # FigureCanvas reports a large matplotlib best size.  In a horizontal
        # wx sizer that can leave no room for the canvas when the frame is
        # narrower than the combined best sizes.  Give wx a modest explicit
        # minimum; proportion=1 below then assigns all remaining width here.
        self.canvas.SetMinSize((320, 260))
        self.canvas.mpl_connect('button_press_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)


        plot_static = wx.StaticBox(self, label="Decay curve")
        plot_box = wx.StaticBoxSizer(plot_static, wx.VERTICAL)
        self.canvas.Reparent(plot_static)
        self.toolbar.Reparent(plot_static)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.TOP, 4)
        plot_box.Add(self.vbox, 1, wx.EXPAND | wx.ALL, 6)
        # Positive proportion is essential: the matplotlib canvas grows and
        # shrinks with the frame instead of retaining a fixed best size.
        self.fullSizer.Add(plot_box, 1, wx.EXPAND | wx.TOP | wx.RIGHT | wx.BOTTOM, 10)

    """
        
    #read in G2hist file for distance and intensity plots
    #combine results for subgraph/polishmodes
    def GetDistHist(self):
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                #print 'Reading all subgraphs'
                for ii in range(len(list(self.parent.inst.subgraphRef.keys()))):
                    testfile=self.testdir+'/G2hist_%i.hist' % (ii+1)
                    if(os.path.exists(testfile)):
                        histA,badNOEA,intyA=self.GetDistHistCore(testfile)
                        try:
                            hist[:,2]+=histA[:,2]
                        
                            for n in intyA:
                                inty.append(n)
                            for n in badNOEA:
                                    badNOE.append(n)
                        except:
                            hist=copy.deepcopy(histA)
                            badNOE=copy.deepcopy(badNOEA)
                            inty=copy.deepcopy(intyA)
                        #print testfile,len(inty),len(intyA)
            else:
                #for ii in range(len(self.parent.inst.subgraphRef.keys())):
                    print('Reading subgraph',sele)
                    testfile=self.testdir+'/G2hist_%s.hist' % (sele)
                    if(os.path.exists(testfile)):
                        histA,badNOEA,intyA=self.GetDistHistCore(testfile)
                        try:
                            hist[:,2]+=histA[:,2]
                        
                            for n in intyA:
                                inty.append(n)
                            for n in badNOEA:
                                    badNOE.append(n)
                        except:
                            hist=copy.deepcopy(histA)
                            badNOE=copy.deepcopy(badNOEA)
                            inty=copy.deepcopy(intyA)
                        #print testfile,len(inty),len(intyA)
        else:
            if(os.path.exists(self.testdir+'/G2hist.hist')):
                hist,badNOE,inty=self.GetDistHistCore(self.testdir+'/G2hist.hist')
            else:
                print('Cannot find',self.testdir+'/G2hist.hist')
        try:
            for i in range(len(hist)):
                if(hist[i,1]!=0):
                    hist[i,3]=hist[i,2]/hist[i,1]
                else:
                    hist[i,3]=0.
            return hist,badNOE,numpy.array(inty)
        except:
            return 0,0,0


    def GetDistHistCore(self,infile):
        tmp=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            tmp.append(test)
        cnt=0
        badNOE=[]
        hist=[]
        inty=[]
        for i in range(len(tmp)):
            try: #check if this line, and the next one are empty
                if(len(tmp[i])==0  and len(tmp[i+1])):
                    cnt+=1
            except:
                pass
                    
            if(cnt==1):
                try:
                    if(int(tmp[i][5])==0):
                        badNOE.append((tmp[i][1],tmp[i][2]))
                except:
                    pass

            if(cnt==2 and len(tmp[i])!=0):
                dist=float(tmp[i][0])
                pdb=float(tmp[i][1])
                soln=float(tmp[i][2])
                ratio=float(tmp[i][3])
                hist.append((dist,pdb,soln,ratio))

            if(cnt==1 and len(tmp[i])!=0):
                dist=float(tmp[i][9])
                disterr=float(tmp[i][10])
                iV=float(tmp[i][7])
                iVerr=float(tmp[i][8])
                if(disterr=='nan'):
                    disterr=0.
                if(disterr!=disterr):
                    disterr=0.

                #print dist,disterr,iV,iVerr
                inty.append((dist,disterr,iV,iVerr))

        hist=numpy.array(hist)
        return hist,badNOE,inty


    def GetDistHistNOE(self,infile):
        dat=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            dat.append(test)
        inny.close()
        cnt=0
        for i,tmp in enumerate(dat):
            try: #check if this line, and the next one are empty
                if(len(tmp)==0  and len(dat[i+1])):
                    cnt+=1
            except:
                pass
                    
            if(cnt==1 and len(tmp)>0):
                row=[]
                #row.append(int(tmp[0])) #index
                row.append(len(self.noe)) #index
                row.append(tmp[1]) #label1
                row.append(tmp[2]) #label2
                row.append(tmp[3]) #type
                row.append(int(tmp[4])) #weight
                row.append(int(tmp[5])) #number
                row.append(tmp[6]) #weighting
                row.append(tmp[7]) #intensity
                row.append(tmp[8]) #scr
                row.append('%.3f' % float(tmp[9])) #distance

                if(tmp[10]=='nan'):
                    row.append(0) #stdev
                else:
                    if(float(tmp[10])<0.1):
                        row.append(0) #stdev
                    else:
                        row.append('%.3f'% float(tmp[10])) #stdev


                self.noe.append(row)

    def GetConvFile(self,inconv1,inconv2):
        inny=open(inconv1)
        in1=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)==3):
                in1.append((float(test[0]),float(test[1]),float(test[2])))
        inny.close()

        inny=open(inconv2)
        in2=[]
        tost=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)==3):
                tost.append((float(test[0]),float(test[1]),float(test[2])))
            else:
                if(len(tost)>0):
                    in2.append(numpy.array(tost))
                tost=[]
        if(len(tost)>0):
            in2.append(numpy.array(tost))
        inny.close()
        return in1,in2
        

    def GetConv(self):
        analdir=self.progress[self.colVal][1]
        if(os.path.exists(analdir)==0):
            analdir=self.parent.inst.P.outdir
        #print analdir
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele!='all'):
                inconv1=analdir+'/conv.1_%i.out' % int(sele)
                inconv2=analdir+'/conv.2_%i.out' % int(sele)
                if(os.path.exists(inconv1)==1):
                    return self.GetConvFile(inconv1,inconv2)
        else:
                inconv1=analdir+'/conv.1.out' 
                inconv2=analdir+'/conv.2.out' 
                if(os.path.exists(inconv1)==1):
                    return self.GetConvFile(inconv1,inconv2)
        return 0,0
    """

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=line.split()
            peak.append(linetosave)
        peakfile.close()
        return peak
    def GetParVarianFlt(self,parfile,param):
        procpar=self.readfile(parfile)
        #print( procpar)
        tick=0
        args=[]
        for j in range(len(procpar)):
            if(procpar[j][0]==param):
                tick=1
                for k in range(int(procpar[j+1][0])):
                    try:
                        args.append(float(procpar[j+1][k+1]))
                    except:
                        args.append(float(procpar[j+1+k][0]))

        if(tick==0):
            sys.stdout.write('Could not find param %s in procpar\n' % (param))
            return False
        return numpy.array(args)


    
    def Get_Decay(self):
        """Return the conversion-defined pseudo coordinate for Decay analysis."""
        try:
            return self._pseudo_axis_values()
        except PseudoAxisError as exc:
            self._set_status(str(exc))
            return numpy.array([-1.0])


    #turn fuda files into  .cpmg files
    def MakeDecaycurve(self,pk,infile):

        import math
        self.typ='bruk'

        try:
            self.Time_T2=float(self.TimeT2Box.GetValue())
        except ValueError:
            self._set_status('Time Mult must be numeric')
            return
        self._save_time_mult()

        # pseudo_axis.tsv stores the conversion-defined coordinate.  Decay uses
        # physical times obtained by applying the user/project time multiplier.
        nucpmg=self.Get_Decay()
        if(nucpmg[0]==-1):
            self._set_status('No usable pseudo-axis values')
            return
            
        nucpmg=numpy.array(nucpmg, dtype=float) * self.Time_T2
        #nucpmgActual=nucpmg[nucpmg!=0]
        nucpmgActual=nucpmg

            
        #print ('Nonzero Xvalues:',nucpmgActual)
        #field=d[0]/TIME_T2

        
        # Get name of the peak
        Name=pk
        # Load data file
        inputfile = open(infile,'r')
        lines=inputfile.readlines()
        data_lines = [line for line in lines if line.strip() and not line.lstrip().startswith('#')]
        if len(data_lines) != len(nucpmg):
            inputfile.close()
            self._set_status('Pseudo-axis has %d values; fitted data contain %d planes' % (len(nucpmg), len(data_lines)))
            return
        Data=[]
        dublicate=[]
        j=0
        for line in lines:
            if ( len(line.split()) > 2 ):
                if line.split()[1]=='f01(ppm)':
                    Noffset=float(line.split()[2])
                if line.split()[1]=='f02(ppm)':
                    Hoffset=float(line.split()[2])
            if not ( line[0] == "#"):
                # Read the actual data.
                temp=line.split()
                if(self.typ=='bruk'):
                    if ( math.fabs(float(nucpmg[j]))<1e-6):
                        Ref=[0,float(temp[1]),float(temp[2])]
                    else:
                        Data.append([float(temp[0]),float(temp[1]),float(temp[2])])
                else:
                    if ( math.fabs(float(temp[0]))<1e-6):
                        Ref=[0,float(temp[1]),float(temp[2])]
                    else:
                        Data.append([float(temp[0]),float(temp[1]),float(temp[2])])
                j+=1    
                for d in range(len(Data)-1):
                    if ( math.fabs(float(temp[0])-Data[d][0]) < 1e-6 ):
                        dublicate.append(math.pow(float(temp[1])-Data[d][1],2.))
        #
        # Get dublicate data (Flemming's original spelling mistake)
        if(len(dublicate)==0):
            StdErr=0.3
        else:
            StdErr=0.
            for d in range(len(dublicate)):
                StdErr+=dublicate[d]
            StdErr=math.sqrt(StdErr)/len(dublicate)
        #
        # Calculate the R2 and field
        ofn=open(infile+".decay",'w')
        ofn.write("#%11s%15s%13s\n" % ('nu_cpmg(Hz)','R2(1/s)','Esd(R2)'))
        #print len(Data)
        #sys.exit(100)

        argy=numpy.argsort(nucpmgActual)
        
        self.cpmgX=[]
        self.cpmgY=[]
        self.cpmgE=[]

        for j in argy:
            #for j in range(len(Data)):
            d=Data[j]
            field=nucpmgActual[j]

            #R=math.log(math.fabs(Ref[1]/d[1]))/self.Time_T2
            #Esd=math.sqrt( math.pow(StdErr/Ref[1],2.)+math.pow(StdErr/d[1],2.))/self.Time_T2

            R=d[1]
            Esd=StdErr  #math.sqrt( math.pow(StdErr/Ref[1],2.)+math.pow(StdErr/d[1],2.))/self.Time_T2
            


            ofn.write(" %11.4e%15.6e%13.6e\n" % (field,R,Esd))
            
            self.cpmgX.append(field)
            self.cpmgY.append(R)
            self.cpmgE.append(Esd)

        ofn.close()
        inputfile.close()

        self.cpmgX=numpy.array(self.cpmgX)
        self.cpmgY=numpy.array(self.cpmgY)
        self.cpmgE=numpy.array(self.cpmgE)

        
        """
            if(j==0):
                resXmin=field
                resXmax=field
                resYmax=R
                resYmin=R
                resYmaxerr=Esd
                resYminerr=Esd


            if(field<xmin):
                xmin=field
            if(field>xmax):
                xmax=field
            if(R+Esd>ymax):
                if(R+Esd<100):
                    ymax=R+Esd
            if(R-Esd<ymin):
                if(R-Esd>-100):
                    ymin=R-Esd
                    
            if(field<resXmin):
                resXmin=field
                resYmin=R
                resYmaxerr=Esd
            if(field>resXmax):
                resXmax=field
                resYmax=R
                resYmaxerr=Esd
        #print resXmax,resYmin,resYmax,resYmin
        if( (resYmax+resYmaxerr)< (resYmin-resYminerr) ):
            print 'Dispersion!',val
            if(yminDisp>resYmin-resYminerr):
                yminDisp=resYmin-resYminerr
            if(ymaxDisp<resYmax-resYmaxerr):
                ymaxDisp=resYmax+resYmaxerr
        """

    def ReadDecaycurve(self,infile):
        if(os.path.exists(infile)==False):
            print ('no cpmg file:',infile)
            return
        inny=open(infile)
        cnt=0
        self.cpmgX=[]
        self.cpmgY=[]
        self.cpmgE=[]
        for line in inny.readlines():
            if(cnt!=0):
                test=line.split()
                self.cpmgX.append(float(test[0]))
                self.cpmgY.append(float(test[1]))
                self.cpmgE.append(float(test[2]))
            cnt+=1
        self.cpmgX=numpy.array(self.cpmgX)
        self.cpmgY=numpy.array(self.cpmgY)
        self.cpmgE=numpy.array(self.cpmgE)
        
    def ReadFuda(self,pk):
        infile=self.fuda_dir + os.sep+pk+'.out'
        if(os.path.exists(infile)==False):
            self._set_status('No fitted intensities for %s' % pk)
            return False
        cpmgfile=os.path.join(self.fuda_dir, pk + '.out.decay')
        if(os.path.exists(cpmgfile)==False):
            self.MakeDecaycurve(pk,infile)
            return True
        else:
            self.ReadDecaycurve(cpmgfile)
            return True

    """
    def ReadCatia(self,pk):

        self.fitLocal.ClearAll()
        self.fitLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.fitLocal.InsertColumn(1, 'Value', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.fitLocal.InsertColumn(2, 'Status', width = 80,format=wx.LIST_FORMAT_CENTRE)
        num_items = self.fitLocal.GetItemCount()

        infile=self.raw+'/catia/OutPut/'+pk+'.dat'
        if(os.path.exists(infile)==False):
            return False

        self.catiaX=[]
        self.catiaY=[]
        self.catiaE=[]
        self.catiaF=[]
        self.catiaL=''

        locFlag=0
        cnt=0
        inny=open(infile)
        for line in inny.readlines():
            if(len(line)>0):
                if(line[0]!='#'):
                    test=line.split()
                    if(len(test)==4):
                        self.catiaX.append(float(test[0]))
                        self.catiaY.append(float(test[1]))
                        self.catiaE.append(float(test[2]))
                        self.catiaF.append(float(test[3]))
                else:
                    test=line.split()
                    if(test[0]=='#DataSet:'):
                        self.catiaL=test[1]

                    if(test[0]=='#Atom:'):
                        locFlag=1
                        
                    if(locFlag==0 and len(test)==4):
                        self.fitLocal.InsertStringItem(num_items,str(cnt))
                        self.fitLocal.SetItem(num_items,0,test[1]) #id

                        try:
                            self.fitLocal.SetItem(num_items,1,'%.3f ' % (float(test[2]))) #A
                        except:
                            self.fitLocal.SetItem(num_items,1,'%s ' % ((test[2]))) #A


                        try:
                            self.fitLocal.SetItem(num_items,2,'%.3f ' % (float(test[3]))) #A
                        except:
                            self.fitLocal.SetItem(num_items,2,'%s ' % ((test[3]))) #A

                            
                        cnt+=1

                        

        self.catiaX=numpy.array(self.catiaX)
        self.catiaY=numpy.array(self.catiaY)
        self.catiaE=numpy.array(self.catiaE)
        self.catiaF=numpy.array(self.catiaF)


        self.chi2Global=numpy.average((self.catiaY-self.catiaF)**2.)

        

        self.globout=self.raw+'/catia/OutPut/GlobalParam.fit'
        if(os.path.exists(self.globout)==1):
            inny=open(self.globout)
            for line in inny.readlines():
                #test=line.split('=')
                #if(len(test)==2):
                #        self.fitLocal.InsertStringItem(num_items,str(cnt))
                #        self.fitLocal.SetItem(num_items,0,test[1]) #id
                    
                test=line.split()
                if(len(test)==3):
                    self.fitLocal.InsertStringItem(num_items,str(cnt))
                    self.fitLocal.SetItem(num_items,0,test[0]+' (global)') #id
                    try:
                        self.fitLocal.SetItem(num_items,1,'%.3f ' % (float(test[1]))) #A
                    except:
                        self.fitLocal.SetItem(num_items,1,'%s ' % ((test[1]))) #A


                    try:
                        self.fitLocal.SetItem(num_items,2,'%.3f ' % (float(test[2]))) #A
                    except:
                        self.fitLocal.SetItem(num_items,2,'%s ' % ((test[2]))) #A

                            
                    cnt+=1


        if(pk not in self.cpmgLocal.keys()):
            self.DoFit(pk)
                    
        for par in 'R0line','pb','kex','R0','dw':
            self.fitLocal.InsertStringItem(num_items,str(cnt))
            self.fitLocal.SetItem(num_items,0,par+' (local)') #id
            val=self.cpmgLocal[pk][par]
            self.fitLocal.SetItem(num_items,1,'%.3f ' % (val)) #A
            cnt+=1
            

        self.fitLocal.InsertStringItem(num_items,str(cnt))
        self.fitLocal.SetItem(num_items,0,'chi2Line') #id
        self.fitLocal.SetItem(num_items,1,'%.3f ' % (self.chi2Line)) #A
        cnt+=1

        self.fitLocal.InsertStringItem(num_items,str(cnt))
        self.fitLocal.SetItem(num_items,0,'chi2Local') #id
        self.fitLocal.SetItem(num_items,1,'%.3f ' % (self.chi2Local)) #A
        cnt+=1

        try:
            self.fitLocal.InsertStringItem(num_items,str(cnt))
            self.fitLocal.SetItem(num_items,0,'chi2Global') #id
            self.fitLocal.SetItem(num_items,1,'%.3f ' % (self.chi2Global)) #A
            cnt+=1
        except:
            pass

            
            
        return True
    """    
    def draw_figure(self,event):
        self._plot_current_mode()
        return

        sele=self.datasets.GetFirstSelected()
        #print(sele)
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)][sele]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)][sele]
        print('Selected:',col1,col2)

        if(self.ReadFuda(col1)==False):
            return



        self.DoFit(col1)
        
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.clear()

        self.ax.errorbar(self.cpmgX,self.cpmgY/self.A0,yerr=(self.cpmgE/abs(self.A0)),fmt='o')                



        #if(self.ReadCatia(col1)==True):
        #self.ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.catiaE,fmt='o',label='raw')                
        #    self.ax.plot(self.cpmgX,self.cpmgL,label='line')
        self.ax.plot(self.cpmgX,self.cpmgF/self.A0,label='local')
        #    self.ax.plot(self.catiaX,self.catiaF,label='global')
        #else:
        #self.ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.cpmgE,fmt='o',label='raw')                
        #self.ax.plot(self.cpmgX,self.cpmgL,label='line')
        #self.ax.plot(self.cpmgX,self.cpmgF,label='local')


        self.ax.legend(loc='upper right')
        self.ax.set_xlabel("t")
        self.ax.set_ylabel("I/I_0")
        self.canvas.draw()
        
        return
        
        """
        #plt.grid(True)
        #self.ax2=self.ax.twinx()
        #self.ax2.set_ylabel('Ratio',color='r')
        #    self.ax2.set_ylim([0.0,1.05])

        #self.ax.set_xlim([3.8,20])
        #self.ax2.set_xlim([3.8,20])

        l1=self.ax.plot(hist[:,0],hist[:,1],label='PDB')
        l2=self.ax.plot(hist[:,0],hist[:,2],label='Data')
        l3=self.ax2.plot(hist[:,0],hist[:,3],color='r',label='ratio = Data/PDB')

        self.ax.fill_between(hist[:,0],0,hist[:,1],facecolor='blue',alpha=0.1)
        self.ax2.fill_between(hist[:,0],0,hist[:,3],facecolor='red',alpha=0.1)
        self.ax.fill_between(hist[:,0],0,hist[:,2],facecolor='green',alpha=0.1)

        self.ax.set_xlabel('C-C Distance (A)')
        self.ax.set_ylabel('Counts')
        #ax.set_title(analdir)
    
            ymin,ymax=self.ax.get_ylim()
            self.ax.set_ylim([0,ymax])
            #print ymin,ymax
            
            #ax.legend((l1,l2,l3),('PDB','Data','Ratio'),'upper right')
            #ax.legend()
            #ax2.legend()
            
            import matplotlib.patches as mpatches
            red_patch = [mpatches.Patch(color='red', label='ratio')]
            blue_patch = [mpatches.Patch(color='blue', label='PDB')]
            green_patch = [mpatches.Patch(color='green', label='Data')]

            handles,labels=self.ax.get_legend_handles_labels()
            #handles=#handles+red_patch
            handles=blue_patch+green_patch+red_patch
            labels=labels+['ratio']
            self.ax.legend(handles,labels,loc=1)

            #fig.legend(handles=
            #ax.legend(('r','b','g'),("A","B","C"),loc=1)
            self.ax2.tick_params('y',colors='r')

            

        elif(seleRbox==3):
            hist,badNOE,inty=self.GetDistHist()
            try:
                if(hist==0):
                    print('No G2hist file')
                    return
            except:
                pass

            #print inty
            #print len(inty)

            self.inty=inty
            self.ax.errorbar(inty[:,0], inty[:,2], xerr=inty[:,1], yerr=inty[:,3],fmt='o')
            yvals=numpy.ones_like(inty[:,0])*0.5

            
            if(self.sheetbox.GetSelection()!=2):
                self.sheetbox.SetSelection(2)
                self.OnTickFilt(True)        
            if(self.sheetbox.GetStringSelection()=='NOEs'):
                #try:
                    sele=self.source.GetFirstSelected()
                    count = self.source.GetItemCount()
                    xv = [self.source.GetItem(row, 7).GetText() for row in range(count)][sele]
                    yv = [self.source.GetItem(row, 3).GetText() for row in range(count)][sele]

                    #print xv,yv
                    self.ax.scatter(float(xv),float(yv),color='r',s=200)

            self.ax.plot(inty[:,0],yvals)                

            self.ax.set_xlabel("distance(A)")
            self.ax.set_ylabel("CrosspeakIntensity")
            
            
            xmin,xmax=self.ax.get_xlim()
            ymin,ymax=self.ax.get_ylim()

            #print inty[:,0]
            print('xmax:',numpy.max(inty[:,0]))

            self.ax.set_xlim(0,numpy.max(inty[:,0])*1.05)
            self.ax.set_ylim(0,1.1)


        elif(seleRbox==0):
            self.OnTickFilt(True)
            try:
                prog,title=self.GetProg()
            except:
                self.fig.clear()
                print('No mces file to plot')
                return

            try:
                testmax=0
                for pro in prog:
                    test=numpy.max(pro[:,2])
                    if(test>testmax):
                        testmax=test
                testmin=testmax
                for pro in prog:
                    test=numpy.min(pro[:,2])
                    if(test<testmin):
                        testmin=test            
                xmax=0
                for pro in prog:
                    test=numpy.max(pro[:,0])
                    if(test>xmax):
                        xmax=test            
                ymax=0
                for pro in prog:
                    test=numpy.max(pro[:,1])
                    if(test>ymax):
                        ymax=test            
                if(xmax==0):
                    print('No plotable data')
                    return 
            except:
                print('Error in setting axes')
                return

            jet=cm.get_cmap('jet')

            lines=[]
            zed=[]
            for pro in prog:
                pro=numpy.array(pro)
                for i in range(len(pro)-1):
                    lines.append((list(pro[i:i+2,0]),list(pro[i:i+2,1])))
                    zed.append(pro[i+1,2])

            zed=numpy.array(zed)
            lines=[list(zip(x,y)) for x,y in lines]
            lines=LineCollection(lines,array=zed,cmap=jet)
            self.ax.add_collection(lines)
                
            #if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            self.ax.set_title("%s" % (title),fontsize=8)

            try:
                self.ax.set_xscale('log')
            except:
                pass
            self.ax.set_xlim(1,xmax*1.05)
            self.ax.set_ylim(-2,ymax+1)

            cb=self.fig.colorbar(lines)
            cb.set_label("NMR EdgeScore",fontsize=8)
            self.ax.set_xlabel("IterationNumber",fontsize=8)
            self.ax.set_ylabel("G1node",fontsize=8)

            self.fig.text(0.15,0.15,"current best score: %i" % testmax,fontsize=8)
            #plt.text(0.02, 0.5, textstr, fontsize=14, transform=plt.gcf().transFigure)

        elif(seleRbox==1):
            in1,in2=self.GetConv()
            if(in1==0):
                print('No optimisation data')
                return
            jet=cm.get_cmap('jet')
            #print in1,in2

            norm = matplotlib.colors.Normalize(vmin=0.0, vmax=in1[len(in1)-1][0])#normalise colour bar


            in1=numpy.array(in1)
            in2=numpy.array(in2)
            for i in range(len(in2)):
                ii=len(in2)-i-1
                self.ax.bar(in2[ii][:,1],in2[ii][:,2],color=(jet(norm(1+in2[ii][0,0]))),alpha=1)

            y = numpy.array([1,1+len(in1)])
            colors = cm.jet(y / float(max(y)))
            sm = cm.ScalarMappable(cmap=cm.jet, norm=matplotlib.colors.Normalize(vmin=1, vmax=(1+len(in1))))
            sm._A = []
            self.fig.colorbar(sm)


            self.ax2=self.ax.twinx()
            self.ax2.set_ylabel('Ratio',color='r')


            cax=self.ax2.plot(in1[:,1],in1[:,0],color='r',label='best',linewidth=1)


            if(self.ref=='subgraphMode' or self.ref=='polishMode'):
                self.ax.set_title("optimisation: %s" % (self.subgraph.GetStringSelection()),fontsize=8)
            else:
                self.ax.set_title("optimisation ",fontsize=8)


            self.ax.set_xlabel("Size of MCES",fontsize=8)
            self.ax.set_ylabel("Count",fontsize=8)
            self.ax2.set_ylabel('Largest MCES progress',color='r')

            inbig=in1[len(in1)-1,1]
            for i in range(len(in1)):
                if(in1[i,1]==inbig):
                    imin=i
                    break

            self.fig.text(0.15,0.75,"largest mces found: %i" % inbig,fontsize=8)
            self.fig.text(0.15,0.73, "discovered in iteration: %i" % imin,fontsize=8)

        """


    def ReadProg(self,infile):
        dat=[]
        dit=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)==0):
                if(len(dit)!=0):
                    dat.append(numpy.array(dit))
                dit=[]

            if(len(test)>0):
                x=float(test[0])
                y=float(test[2])
                c=float(test[4])
                dit.append((x,y,c))


        if(len(dat)==0):
            return numpy.array(dit)
        else:
            if(len(dit)!=0):
                dat.append(numpy.array(dit))
            return dat

    #get progress data
    def GetProg(self):

        dat=[]
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                for i in range(len(list(self.parent.inst.subgraphRef.keys()))):        
                    ii=len(list(self.parent.inst.subgraphRef.keys()))-i-1
                
                    #print 'looking for subgraph',ii+1,' and nproc',int(self.parent.processors.GetValue())
                    if(int(self.parent.processors.GetValue())==1):
                        testfile=self.testdir+'/mces_%i.txt.G' % ((ii+1))
                        if(os.path.exists(testfile)):
                            #print 'reading ',testfile
                            dat.append(self.ReadProg(testfile))
                            return dat,testfile
                    else:
                        testfile=self.testdir+'/mces_%i.txt.G.%i' % ((ii+1),0)
                        if(os.path.exists(testfile)): #if we are still running...
                            for j in range((int(self.parent.processors.GetValue()))):
                                #print 'reading ',testfile
                                testfile=self.testdir+'/mces_%i.txt.G.%i' % ((ii+1),j)
                                dat.append(self.ReadProg(testfile))
                            return dat,testfile
                        testfile=self.testdir+'/mces_%i.txt.G' % ((ii+1),)
                        if(os.path.exists(testfile)):                    
                            #print 'reading',testfile
                            dit=self.ReadProg(testfile)
                            try:
                                #print dit.shape
                                dat.append(dit)
                                return dat,testfile
                            except:
                                return dit,testfile
            else:
                #for i in range(len(self.parent.inst.subgraphRef.keys())):        
                    #print 'looking for subgraph',sele,' and nproc',int(self.parent.processors.GetValue())
                    if(int(self.parent.processors.GetValue())==1):
                        testfile=self.testdir+'/mces_%s.txt.G' % ((sele))
                        if(os.path.exists(testfile)):
                            #print 'reading ',testfile
                            dat.append(self.ReadProg(testfile))
                            return dat,testfile
                    else:
                        testfile=self.testdir+'/mces_%s.txt.G.%i' % ((sele),0)
                        if(os.path.exists(testfile)): #if we are still running...
                            for j in range((int(self.parent.processors.GetValue()))):
                                #print 'reading ',testfile
                                testfile=self.testdir+'/mces_%s.txt.G.%i' % ((sele),j)
                                dat.append(self.ReadProg(testfile))
                            return dat,testfile
                        testfile=self.testdir+'/mces_%s.txt.G' % ((sele),)
                        if(os.path.exists(testfile)):                    
                            print('reading',testfile)
                            dit=self.ReadProg(testfile)
                            if(len(dit)==int(self.parent.processors.GetValue())):
                                return dit,testfile
                            else:
                                dat.append(dit)
                                return dat,testfile



        else:
            if(int(self.parent.processors.GetValue())==1):
                testfile=self.testdir+'/mces.txt.G'
                if(os.path.exists(testfile)):
                    dat.append(self.ReadProg(testfile))
                return dat,testfile
            else:
                testfile=self.testdir+'/mces.txt.G.%i' % (0)
                if(os.path.exists(testfile)):
                    for j in range((int(self.parent.processors.GetValue()))):
                        testfile=self.testdir+'/mces.txt.G.%i' % (j)
                        if(os.path.exists(testfile)):
                            dat.append(self.ReadProg(testfile))
                    return dat,testfile
                testfile=self.testdir+'/mces.txt.G'
                if(os.path.exists(testfile)):                    
                    return self.ReadProg(testfile),testfile

        self.GetLastFile() #if the last file is in community mode...
        if(len(self.fil.split('community'))>1):
            testfile=self.testdir+'/'+self.fil
            return self.ReadProg(testfile),testfile
                    

    def PopulatelistNOE(self):
        result_dict={}
        self.noe=[]
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            
            if(sele=='all'):
                noe=[]
                for ii in range(len(list(self.parent.inst.subgraphRef.keys()))):
                    testfile=self.testdir+'/G2hist_'+str(ii+1)+'.hist'
                    if(os.path.exists(testfile)):
                        self.GetDistHistNOE(testfile)
            else:
                testfile=self.testdir+'/G2hist_'+str(sele)+'.hist'
                if(os.path.exists(testfile)):
                    self.GetDistHistNOE(testfile)

        else:
            testfile=self.testdir+'/G2hist.hist'
            if(os.path.exists(testfile)):
                self.GetDistHistNOE(testfile)

        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        cnt=0
        for i,noe in enumerate(self.noe):
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertStringItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(noe[0])) #id
            self.source.SetItem(num_items,1,str(noe[1])) #A
            self.source.SetItem(num_items,2,str(noe[2])) #B
            self.source.SetItem(num_items,3,str(noe[7])) #intensity
            self.source.SetItem(num_items,4,str(noe[4])) # number
            self.source.SetItem(num_items,5,str(noe[6])) #w/s
            self.source.SetItem(num_items,6,'%.2f' % (noe[5]*1./(1.*self.nsoln)) ) #hits
            self.source.SetItem(num_items,7,str(noe[9])) #dist
            self.source.SetItem(num_items,8,str(noe[10])) #stdev

            if(int(noe[5])==0):
                color = (int(255),0, 0)
                self.source.SetItemBackgroundColour(num_items,color)

            """
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(1, 'A', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(2, 'B', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(3, 'intensity', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(4, '#', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(5, 'tp', width = 50,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(6, 'hits', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(7, 'dist', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            self.source.InsertColumn(8, 'stdev', width = 80,format=wx.LIST_FORMAT_CENTRE) 
            """
            #stry=''
            #for val in vals:
            #    stry+=val+' '
            #self.source.SetItem(num_items,3,stry)

        return

    def PopulateFilt(self):
        result_dict={}
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()

            if(sele=='all'):
                for ii in range(len(list(self.parent.inst.subgraphRef.keys()))):
                    testfile=self.testdir+'/filter_'+str(ii+1)+'.res'
                    if(os.path.exists(testfile)):
                        analysis.AddFile(result_dict,testfile)
            else:
                testfile=self.testdir+'/filter_'+str(sele)+'.res'
                if(os.path.exists(testfile)):
                    analysis.AddFile(result_dict,testfile)
        else:
            testfile=self.testdir+'/filter.Full.res'
            if(os.path.exists(testfile)):
                analysis.AddFile(result_dict,testfile)

        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        cnt=0
        for key in sorted(result_dict,key=lambda k:len(result_dict[k])):
            vals=result_dict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertStringItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(cnt))
            self.source.SetItem(num_items,1,key)
            self.source.SetItem(num_items,2,'filter')
            stry=''
            for val in vals:
                stry+=val+' '
            self.source.SetItem(num_items,3,stry)

    def UpdateScores(self):
        tig=0  #set to 1 if all or not subgraph/polish modes
        if( self.ref=='subgraphMode' or self.ref=='polishMode' ):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1
        elif(self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1            
        else:
            tig=1

        ret=0 #return flag
        if(tig==1): #if not subgraph/polish nor sele='all'
            self.result_dict={}
            self.result_dict,nsoln,edges,edgestot,edgesG2=analysis.AnalAll(self.testdir)
        else:
            self.result_dict={}
            mcesfile=self.testdir+'/mces_%s.txt' % (int(sele))
            vf2file=self.testdir+'/vf2_%s.txt' % (int(sele))

            ret=0
            if(os.path.exists(mcesfile)):
                nsoln,edges,edgestot,edgesG2,g1node=analysis.AddFile(self.result_dict,mcesfile)
                mode='mces'
            elif(os.path.exists(vf2file)):
                mode='vf2 '
                nsoln,edges,edgestot,edgesG2,g1node=analysis.AddFile(self.result_dict,vf2file)
            else:
                print('Neither vf2 nor mces file exists for subgraph',sele)
                self.GetLastFile()
                nsoln=0
                edges=0
                edgestot=0
                edgesG2=0
                if(len(self.fil.split('vf2'))>1):
                    mode='vf2'
                else:
                    mode='mces'
                ret=1

        #print len(self.result_dict.keys())
        self.conf=0
        self.tot=len(list(self.result_dict.keys()))
        for key,vals in list(self.result_dict.items()):
            if(len(vals)==1):
                self.conf+=1

        self.nsoln=nsoln
        self.confLbl.SetLabel("Confident:     %i / %i" % (self.conf,self.tot))
        self.solnLbl.SetLabel("Solutions:     %i" % nsoln)
        self.nmreLbl.SetLabel("NMR EdgeScore: %i / %i" % (edges,edgestot))
        self.pdbeLbl.SetLabel("PDB EdgeScore: %i" % edgesG2)

        #get current result_dict
        if(len(list(self.result_dict.keys()))==0):
            ret=1

        return ret




    def PopulateList(self,shift='n'):
        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        if(self.UpdateScores()):
            return

        #populate list
        resultConf={}
        conf=0
        for key,vals in list(self.result_dict.items()):
            if(len(vals)==1):
                conf+=1
                resultConf[key]=vals
                del(self.result_dict[key])
        #print 'conf',conf
        cnt=0
        for i in range(self.colVal+1):
            if(os.path.exists(self.progress[i][1]+'/confident.res')): #is this mode done?
                #print 'testing:',self.progress[i][1]+'/confident.res'
                new={}
                inny=open(self.progress[i][1]+'/confident.res')
                for line in inny.readlines():
                    test=line.split(':')
                    key=test[0].split()[0]
                    ass=test[1].split()[0]
                    if(key in list(resultConf.keys())):
                        if(ass==resultConf[key][0]):
                            new[key]=ass
                            del(resultConf[key])
                        #else:
                            #print 'ODD? confident solution not in final results'
                            #print key,ass
                inny.close()

                for key,vals in list(new.items()):
                    cnt+=1

                    num_items = self.source.GetItemCount()
                    if(self.WXV==4):
                        self.source.InsertItem(num_items,str(cnt))
                        self.source.SetItem(num_items,0,str(cnt))
                        self.source.SetItem(num_items,1,key)
                        self.source.SetItem(num_items,2,self.progress[i][0])
                    else:
                        self.source.InsertStringItem(num_items,str(cnt))
                        self.source.SetItem(num_items,0,str(cnt))
                        self.source.SetItem(num_items,1,key)
                        self.source.SetItem(num_items,2,self.progress[i][0])



                    if(shift=='n'):
                        stry=vals
                    else:
                        stry=vals+"(%.2f)" % (self.shiftDict[key][vals]['sc'])
                    if(self.WXV==4):
                        self.source.SetItem(num_items,3,stry)
                    else:
                        self.source.SetItem(num_items,3,stry)

                    if(self.progress[i][0]=='subgraphMode' or self.progress[i][0]=='polishMode' or self.progress[i][0]=='nudgeMode' or self.progress[i][0]=='finalMode'):

                        #color = (int(255*percentage), 0, int(255*(1.-percentage)))
                        color = (0,int(255), 0)
                        #self.source.SetCellBackgroundColour(num_items, 2, color)
                        self.source.SetItemBackgroundColour(num_items,color)


                    #    color='green'
                    #else:
                    #    color='yellow'
                    #outlat.write('%s \\begin{footnotesize}\\colorbox{%s!30}{%s}\\end{footnotesize} ' % (vals,color,self.progress[i][0]))
                    #v1=(re.findall(r'[0-9]+',vals)[0])
                    #v2=(re.findall(r'[0-9]+',key)[0])
                    #if(v1!=v2):
                    #    outlat.write('\\begin{footnotesize}\\colorbox{red!30}{%s}\\end{footnotesize} ' % ('INCORRECT!'))
                    #outlat.write('\n\n')                    

        for key in sorted(self.result_dict,key=lambda k:len(self.result_dict[k])):
            vals=self.result_dict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            if(self.WXV==4):
                self.source.InsertItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,str(cnt))
                self.source.SetItem(num_items,1,key)
            else:
                self.source.InsertStringItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,str(cnt))
                self.source.SetItem(num_items,1,key)

 

            #self.source.SetItem(num_items,2,self.progress[i][0])
            stry=''
            for val in vals:
                if(shift=='n'):
                    stry+=val+' '
                else:
                    stry+=val+"(%.2f)" % (self.shiftDict[key][val]['sc'])
            if(self.WXV==4):
                self.source.SetItem(num_items,3,stry)
            else:
                self.source.SetItem(num_items,3,stry)
            
            #self.source.ForceRefresh()
        #print len(self.result_dict.keys())
