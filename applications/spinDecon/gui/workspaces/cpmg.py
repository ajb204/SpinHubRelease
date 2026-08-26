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
from spinDecon.domain.peaks import peak_sort_key
from spinDecon.gui.dialogs.pseudo_axis import show_pseudo_axis_table
from spinDecon.analysis.cpmg_service import (GAMMA, observe_frequency_mhz, build_r2eff, baldwin_r2eff,
    fit_local as fit_cpmg_local, fit_global as fit_cpmg_global, r2_infinity)

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


    
class CPMGMan(wx.App):
    def __init__(self,inherit,pth='',auto_prepare=False):
        self.frame_ProcessFrame=CPMGFrame(None,30,'CPMG Analysis',inherit,pth=pth)
        #FGA added
        self.frame_ProcessFrame.Centre(direction=wx.BOTH)
        self.frame_ProcessFrame.Show(True)
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

class CPMGAdvancedDialog(wx.Frame):
    """Advanced fitting and CATIA controls without exposing legacy frame widgets."""
    def __init__(self, owner):
        wx.Frame.__init__(self, owner, title='Advanced CPMG Analysis', size=(650, 620), style=wx.DEFAULT_FRAME_STYLE|wx.RESIZE_BORDER)
        self.owner=owner
        root=wx.BoxSizer(wx.VERTICAL)

        pars=wx.StaticBoxSizer(wx.StaticBox(self,label='Fit parameters'),wx.VERTICAL)
        self.paramList=wx.ListCtrl(self,style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        for i,(n,w) in enumerate((('Parameter',90),('Initial',85),('Units',70),('Fit',60),('Scope',90),('Meaning',200))): self.paramList.InsertColumn(i,n,width=w)
        # Keep the Initial column numeric-only.  Units are deliberately separate
        # so values can be parsed safely when constructing optimiser p0.
        rows=[('R0','from data','s-1','Yes','Local','baseline R2,eff'),('dw','1.0','ppm','Yes','Local','chemical-shift difference'),('kex','1000','s-1','Yes','Global*','exchange rate'),('pb','0.02','fraction','Yes','Global*','minor-state population')]
        for i,row in enumerate(rows):
            self.paramList.InsertItem(i,row[0])
            for j,val in enumerate(row[1:],1): self.paramList.SetItem(i,j,val)
        pars.Add(self.paramList,1,wx.EXPAND|wx.ALL,6)
        pars.Add(wx.StaticText(self,label='* kex and pb are shared in global fits; R0 and dw remain peak-specific.'),0,wx.LEFT|wx.RIGHT|wx.BOTTOM,6)
        root.Add(pars,1,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,8)

        fitbox=wx.StaticBoxSizer(wx.StaticBox(self,label='Internal fitter'),wx.HORIZONTAL)
        fit_included=wx.Button(self,label='Fit included residues')
        fit_included.SetToolTip('Globally fit every residue currently marked Include=True')
        fit_included.Bind(wx.EVT_BUTTON,self.OnFitIncluded)
        fitbox.Add(fit_included,1,wx.ALL,4)
        results=wx.Button(self,label='Global fit results...')
        results.SetToolTip('Show shared global parameters, per-residue parameters and goodness-of-fit metrics')
        results.Bind(wx.EVT_BUTTON,lambda e:self.owner.OnGlobalFitResults(e))
        root.Add(results,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,8)

        root.Add(fitbox,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,8)

        catia=wx.StaticBoxSizer(wx.StaticBox(self,label='CATIA'),wx.VERTICAL)
        catia_grid=wx.FlexGridSizer(2,4,6,8); catia_grid.AddGrowableCol(1,1); catia_grid.AddGrowableCol(3,1)
        catia_grid.Add(wx.StaticText(self,label='Sequence:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.sequence=wx.ComboBox(self,choices=['Trosy_CPMG','PE_CPMG','CW_CPMG'],style=wx.CB_READONLY); self.sequence.SetValue(owner.seqfilCombo.GetValue() or 'Trosy_CPMG'); catia_grid.Add(self.sequence,1,wx.EXPAND)
        catia_grid.Add(wx.StaticText(self,label='Basis:'),0,wx.ALIGN_CENTER_VERTICAL)
        self.basis=wx.ComboBox(self,choices=['IphAph_13','TrATr_13','Iph_7'],style=wx.CB_READONLY); self.basis.SetValue(owner.basisCombo.GetValue() or 'IphAph_13'); catia_grid.Add(self.basis,1,wx.EXPAND)
        catia.Add(catia_grid,0,wx.EXPAND|wx.ALL,6)
        catia_buttons=wx.BoxSizer(wx.HORIZONTAL)
        for label,handler in [('Export files',self.OnExportCatia),('Run CATIA',self.OnRunCatia),('Import selected',self.OnImportCatia)]:
            b=wx.Button(self,label=label); b.Bind(wx.EVT_BUTTON,handler); catia_buttons.Add(b,1,wx.ALL,4)
        catia.Add(catia_buttons,0,wx.EXPAND)
        persistence=wx.BoxSizer(wx.HORIZONTAL)
        for label,handler in [('Save CATIA settings',self.OnSaveCatiaSettings),('Load CATIA settings',self.OnLoadCatiaSettings)]:
            b=wx.Button(self,label=label); b.Bind(wx.EVT_BUTTON,handler); persistence.Add(b,1,wx.ALL,4)
        catia.Add(persistence,0,wx.EXPAND)
        root.Add(catia,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,8)
        close=wx.Button(self,wx.ID_CLOSE,label='Close'); close.Bind(wx.EVT_BUTTON,lambda e:self.Close()); root.Add(close,0,wx.ALIGN_RIGHT|wx.ALL,8)
        self.SetSizer(root); self.SetMinSize((580,500))

    def _sync_legacy_catia(self):
        o=self.owner; o.seqfilCombo.SetValue(self.sequence.GetValue()); o.basisCombo.SetValue(self.basis.GetValue()); o.DoLocal(None)

    def _selected_peak(self):
        return self.owner.get_selected_peak()

    def _global_initial_values(self):
        """Read numeric shared starting values from the parameter table.

        The Units column is intentionally separate, so strings such as
        '1000 s-1' never reach the optimiser.  Non-numeric cells simply fall
        back to the fitter defaults.
        """
        initial={}
        for row in range(self.paramList.GetItemCount()):
            name=self.paramList.GetItemText(row,0)
            if name not in ('kex','pb'):
                continue
            text=self.paramList.GetItemText(row,1).strip()
            try:
                initial[name]=float(text)
            except (TypeError,ValueError):
                pass
        return initial

    def _global_fit(self, peaks):
        o=self.owner; o.GetActualField(); curves={}
        for pk in peaks:
            if o.ReadFuda(pk): curves[pk]={'x':o.cpmgX.copy(),'y':o.cpmgY.copy(),'e':o.cpmgE.copy()}
        if not curves:
            wx.MessageBox('No readable CPMG curves were available for the global fit.','Global CPMG fit',wx.OK|wx.ICON_INFORMATION); return
        try: result=fit_cpmg_global(curves,float(o.TimeT2Box.GetValue()),o.dfrq,initial=self._global_initial_values())
        except Exception as exc:
            wx.MessageBox(str(exc),'Global CPMG fit',wx.OK|wx.ICON_ERROR); return
        if not result.get('success',False):
            wx.MessageBox(result.get('message','Global fit did not converge'),'Global CPMG fit',wx.OK|wx.ICON_WARNING); return
        # Keep local fits intact: they are the comparison curve.  Enrich the
        # global result independently with per-residue metrics used by the UI.
        for pk,r in result['peaks'].items():
            if pk not in o.cpmgLocal: o.DoFit(pk,update_screen=False)
            local=o.cpmgLocal.get(pk,{})
            r['valid']=True; r['success']=True
            r['Rex']=float(numpy.max(r['model'])-numpy.min(r['model']))
            r['R2inf']=r2_infinity(float(o.TimeT2Box.GetValue()),r['pb'],r['kex'],r['R0'],r['dw']*o.dfrq)
            line=local.get('chi2Line',numpy.nan); r['chi2Line']=line
            r['improvement']=1.-r['chi2local']/line if numpy.isfinite(line) and line>0 else numpy.nan
        o.cpmgGlobal=result
        o._refresh_global_fit_window()
        o._plot_current_mode()
        o.statusBar.SetStatusText('Global fit complete: kex=%s, pb=%s, chi2=%.3g'%(_fit_value_error_precision(result.get('kex'),result.get('kex_error'),' s-1'),_fit_value_error_precision(result.get('pb'),result.get('pb_error')),result['chi2']))
        o.OnGlobalFitResults(None)

    def OnFitIncluded(self,event):
        """Globally fit every residue currently marked include=True.

        This is the single advanced internal-fit action:
        kex and pb are shared across all included residues, while each residue
        retains its own R0/R2inf and delta-omega parameters.
        """
        peaks = self.owner._included_peaks()
        if not peaks:
            wx.MessageBox('No CPMG residues are marked Include=True.',
                          'Global CPMG fit', wx.OK|wx.ICON_INFORMATION)
            return
        self._global_fit(peaks)
    def OnSaveCatiaSettings(self,event):
        self._sync_legacy_catia()
        self.owner.SaveCatiaSettings()
        self.owner.statusBar.SetStatusText('CATIA settings saved')

    def OnLoadCatiaSettings(self,event):
        if not self.owner.LoadCatiaSettings():
            wx.MessageBox('No saved CATIA settings were found.','CATIA',wx.OK|wx.ICON_INFORMATION)
            return
        self.sequence.SetValue(self.owner.seqfilCombo.GetValue() or 'Trosy_CPMG')
        self.basis.SetValue(self.owner.basisCombo.GetValue() or 'IphAph_13')
        self.owner.statusBar.SetStatusText('CATIA settings loaded')

    def OnExportCatia(self,event):
        self._sync_legacy_catia(); o=self.owner; o.PathExists((o.raw+'/catia',o.raw+'/catia/dataset',o.raw+'/catia/OutPut')); o.datfile=[];o.globfile=[];o.locfile=[];o.KexFit=True;o.PbFit=True;o.Conv=1E-3;o.MaxIter=100; o.WriteCatiaDataset();o.WriteCatiaPar();o.WriteCatiaFile(); o.statusBar.SetStatusText('CATIA input files exported')
    def OnRunCatia(self,event): self._sync_legacy_catia(); self.owner.OnButtonRun(event)
    def OnImportCatia(self,event):
        pk=self._selected_peak()
        if pk and self.owner.ReadCatia(pk): self.owner.statusBar.SetStatusText('Imported CATIA results for '+pk)
        elif pk: wx.MessageBox('No CATIA result found for '+pk,'CATIA',wx.OK|wx.ICON_INFORMATION)

class CPMGPeakSelectionFrame(wx.Frame):
    """Detailed peak/fitting table shared by the simple and advanced CPMG views."""
    def __init__(self, owner):
        wx.Frame.__init__(self, owner, title='CPMG Peak Selection / Fit Details', size=(900, 560),
                          style=wx.DEFAULT_FRAME_STYLE|wx.RESIZE_BORDER)
        self.owner=owner
        root=wx.BoxSizer(wx.VERTICAL)
        note=wx.StaticText(self,label='Include reflects the Rex screen. Selection is synchronized with the main CPMG peak selector.')
        root.Add(note,0,wx.EXPAND|wx.ALL,8)
        owner.datasets.Reparent(self)
        owner.datasets.SetMinSize((650,320))
        root.Add(owner.datasets,1,wx.EXPAND|wx.LEFT|wx.RIGHT,8)
        buttons=wx.BoxSizer(wx.HORIZONTAL)
        refresh=wx.Button(self,label='Refresh'); refresh.Bind(wx.EVT_BUTTON,owner.OnButtonRefresh)
        fit=wx.Button(self,label='Analyse'); fit.Bind(wx.EVT_BUTTON,owner.OnAnalyse)
        close=wx.Button(self,label='Close'); close.Bind(wx.EVT_BUTTON,lambda e:self.Close())
        buttons.Add(refresh,0,wx.RIGHT,6); buttons.Add(fit,0,wx.RIGHT,6); buttons.AddStretchSpacer(); buttons.Add(close,0)
        root.Add(buttons,0,wx.EXPAND|wx.ALL,8)
        self.SetSizer(root); self.SetMinSize((700,420))
        self.Bind(wx.EVT_CLOSE,self.OnClose)

    def OnClose(self,event):
        # Keep the shared ListCtrl alive; park it on the main frame while this window is hidden.
        self.owner.datasets.Reparent(self.owner)
        self.owner.datasets.Hide()
        self.owner.peakSelectionFrame=None
        self.Destroy()

def _fit_value_error_precision(value, error, unit=''):
    """Format a fitted value/error to the precision justified by the 1-sigma error."""
    try:
        value=float(value); error=float(error)
        if not (numpy.isfinite(value) and numpy.isfinite(error)) or error <= 0:
            return ('%.4g' % value) + unit
        exponent=int(numpy.floor(numpy.log10(abs(error))))
        first=abs(error)/(10.0**exponent)
        sig=2 if first < 3.0 else 1
        decimals=max(0, -exponent + sig - 1)
        if decimals <= 6 and abs(value) < 1e7:
            fmt='%%.%df' % decimals
            return (fmt % value) + ' ± ' + (fmt % error) + unit
        digits=max(1, sig - exponent - 1)
        return ('%.*g' % (digits, value)) + ' ± ' + ('%.*g' % (sig, error)) + unit
    except (TypeError, ValueError, OverflowError):
        return '-'


class CPMGGlobalFitFrame(wx.Frame):
    """Modeless view of the latest global CPMG fit and its residue metrics."""
    def __init__(self, owner):
        wx.Frame.__init__(self, owner, title='CPMG Global Fit Results', size=(980,560),
                          style=wx.DEFAULT_FRAME_STYLE|wx.RESIZE_BORDER)
        self.owner=owner
        root=wx.BoxSizer(wx.VERTICAL)
        self.summary=wx.StaticText(self,label='No global fit has been run.')
        root.Add(self.summary,0,wx.EXPAND|wx.ALL,8)
        self.table=wx.ListCtrl(self,style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        headers=(('Residue',80),('R0 ± err',115),('Δω ± err (ppm)',125),('Rex',75),('R2inf',75),('χ²',75),('Fit gain',80),('Status',90))
        for i,(name,width) in enumerate(headers): self.table.InsertColumn(i,name,width=width,format=wx.LIST_FORMAT_CENTRE)
        self.table.Bind(wx.EVT_LIST_ITEM_SELECTED,self.OnSelect)
        self.table.Bind(wx.EVT_LIST_ITEM_ACTIVATED,self.OnSelect)
        root.Add(self.table,1,wx.EXPAND|wx.LEFT|wx.RIGHT,8)
        note=wx.StaticText(self,label='kex and pb are shared globally; R0 and Δω are fitted independently for each residue. Selecting a residue updates the main CPMG plot.')
        root.Add(note,0,wx.EXPAND|wx.ALL,8)
        close=wx.Button(self,label='Close'); close.Bind(wx.EVT_BUTTON,lambda e:self.Close()); root.Add(close,0,wx.ALIGN_RIGHT|wx.ALL,8)
        self.SetSizer(root); self.SetMinSize((760,420)); self.Bind(wx.EVT_CLOSE,self.OnClose)
        self.RefreshResults()

    def RefreshResults(self):
        result=getattr(self.owner,'cpmgGlobal',None); self.table.DeleteAllItems()
        if not result:
            self.summary.SetLabel('No global fit has been run. Use Advanced CPMG → Fit included residues.')
            return
        kex_text=_fit_value_error_precision(result.get('kex'),result.get('kex_error'),' s⁻¹')
        pb_text=_fit_value_error_precision(result.get('pb'),result.get('pb_error'))
        self.summary.SetLabel('Global:  kex = %s    pb = %s    overall χ² = %.4g    residues = %d    %s' %
                              (kex_text,pb_text,result.get('chi2',numpy.nan),len(result.get('peaks',{})),result.get('message','')))
        for row,pk in enumerate(sorted(result.get('peaks',{}),key=peak_sort_key)):
            r=result['peaks'][pk]; self.table.InsertItem(row,pk)
            vals=(self.owner._format_value_error(r,'R0','%.3f'),self.owner._format_value_error(r,'dw','%.4f'),
                  self.owner._format_fit_value(r.get('Rex'),'%.3f'),self.owner._format_fit_value(r.get('R2inf'),'%.3f'),
                  self.owner._format_fit_value(r.get('chi2local'),'%.4g'),self.owner._format_fit_value(r.get('improvement'),'%.3f'),
                  'OK' if r.get('valid',r.get('success',False)) else 'Fit failed')
            for col,val in enumerate(vals,1): self.table.SetItem(row,col,val)
        selected=self.owner.get_selected_peak()
        for row in range(self.table.GetItemCount()):
            if self.table.GetItem(row,0).GetText()==selected:
                self.table.Select(row); self.table.EnsureVisible(row); break

    def OnSelect(self,event):
        pk=self.table.GetItem(event.GetIndex(),0).GetText()
        if pk:
            self.owner.peakChoice.SetValue(pk); self.owner._refresh_peak_selector(pk); self.owner._sync_table_selection(pk)
            if self.owner.plotMode.GetStringSelection()=='Peak': self.owner._plot_peak(pk)

    def OnClose(self,event):
        self.owner.globalFitFrame=None; self.Destroy()

class CPMGFrame(wx.Frame):
    
    def __init__(self,parent,id,title,inherit,pth):
        #wx.Panel.__init__(self, parent=parent)
        print('PATH IS:',pth)



        self.parent=inherit
        self.pth=pth
        self.WXV=int(wx.__version__.split('.')[0])


        #FGA changed
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.Point(258, 184), size=wx.Size(800, 480),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'MAGMA results ...')
        #self.SetClientSize(wx.Size(900, 280))
        monitorWidth, monitorHeight = wx.GetDisplaySize()
        initial_w = min(1000, max(760, int((monitorWidth - 120) * 0.78)))
        work_h = wx.GetClientDisplayRect().height
        try:
            main_window = wx.GetTopLevelParent(inherit)
            main_h = main_window.GetSize().height if main_window else inherit.GetSize().height
        except Exception:
            main_h = work_h
        max_open_h = max(560, min(main_h, work_h))
        initial_h = min(780, max_open_h)
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
              pos=wx.DefaultPosition, size=(initial_w, initial_h),
              style=wx.DEFAULT_FRAME_STYLE, title='CPMG Analysis')
        self.SetMinSize((720, min(560, max_open_h)))
        self.SetMaxSize((-1, max_open_h))
        self.statusBar=self.CreateStatusBar(1)
        self.statusBar.SetStatusText('Ready')

        self.SetBackgroundColour('WHITE')
        # Historical code created an unsized child panel here.  Because it was
        # not owned by any sizer it could float over the real controls.
        self.panel=wx.Panel(self,-1)
        self.panel.Hide()
        
        ########

        self.raw=self._analysis_spec_dir(inherit,pth)
        self.acquisition_dir=self._acquisition_dir(inherit,pth)
        self.fuda_dir = inherit._fuda_dir()
            
        self.savefile=os.path.join(self.raw,'cpmg','frame.save')
        if(os.path.exists(self.savefile)==0):
            self.PathExists((os.path.dirname(self.savefile),))
            outy=open(self.savefile,'w')
            outy.close()

            
        self.cpmgLocal={}
        self.cpmgGlobal=None
        self.globalFitFrame=None

        #self.lc=SortedListCtrl(panel,self.corrDict)

        self.datasets=AutoWidthListCtrl(self)
        #self.datasets.SetMinSize((650,300))
        self.datasets.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)
        self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        self.datasets.Bind(wx.EVT_LIST_COL_CLICK, self.OnButtonSort)
        
        
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

        self.buttonWipe = wx.Button(self, label="Wipe")
        self.buttonWipe.Bind(wx.EVT_BUTTON,self.OnButtonWipe)

        self.buttonGuess = wx.Button(self, label="Analyse")
        self.buttonGuess.SetToolTip("Generate CPMG curves, fit every peak, and apply the Rex screen")
        self.buttonGuess.Bind(wx.EVT_BUTTON,self.OnAnalyse)
        

        self.buttonClose = wx.Button(self, label="Close")
        self.buttonClose.Bind(wx.EVT_BUTTON,self.OnButtonClose)

        self.buttonRun = wx.Button(self, label="Run Catia")
        self.buttonRun.Bind(wx.EVT_BUTTON,self.OnButtonRun)

        self.buttonPeakConvert = wx.Button(self, label="Create CPMG curves")
        self.buttonPeakConvert.Bind(wx.EVT_BUTTON,self.OnButtonPeakConvert)

        self.buttonAdvanced = wx.Button(self, label="Advanced CPMG...")
        self.buttonAdvanced.Bind(wx.EVT_BUTTON, self.OnAdvancedCPMG)
        self.buttonPeaks = wx.Button(self, label="Select peaks / details...")
        self.buttonPeaks.Bind(wx.EVT_BUTTON, self.OnPeakSelection)
        self.peakChoiceTxt=wx.StaticText(self,label="Peak:")
        self.peakChoice=wx.ComboBox(self,choices=[],style=wx.CB_READONLY)
        self.peakChoice.Bind(wx.EVT_COMBOBOX,self.OnPeakChoice)
        self.peakSpin=wx.SpinButton(self,style=wx.SP_VERTICAL)
        self.peakSpin.Bind(wx.EVT_SPIN_UP,self.OnPeakUp)
        self.peakSpin.Bind(wx.EVT_SPIN_DOWN,self.OnPeakDown)
        self.screenCheck=wx.CheckBox(self,label="Significant only")
        self.screenCheck.SetToolTip("Show only peaks whose Include flag is True after Rex screening")
        self.screenCheck.Bind(wx.EVT_CHECKBOX,self.OnScreenPeaks)
        self.peakSelectionFrame=None
        self.advancedFrame=None
        self.plotModeTxt = wx.StaticText(self, label="Plot:")
        self.plotMode = wx.RadioBox(self, choices=["Peak", "kex", "pb", "R2inf", "Δω", "Overlay"], majorDimension=2, style=wx.RA_SPECIFY_COLS)
        self.plotMode.Bind(wx.EVT_RADIOBOX, self.OnPlotMode)





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


        

        seqlist=[]
        seqlist.append('Trosy_CPMG')
        seqlist.append('PE_CPMG')
        seqlist.append('CW_CPMG')
        
        basislist=[]
        basislist.append('IphAph_13')
        basislist.append('TrATr_13')
        basislist.append('Iph_7')
        



        
        self.basistxt  = wx.StaticText(self, label="Basis:")
        self.seqfiltxt  = wx.StaticText(self, label="Seqfil:")
        self.basisCombo=wx.ComboBox(self, -1, size=(80, -1),choices=basislist , style=wx.CB_READONLY)
        self.seqfilCombo=wx.ComboBox(self, -1, size=(80, -1),choices=seqlist,  style=wx.CB_READONLY)


        self.RexScreenBoxTxt = wx.StaticText(self, label="RexScreen:")
        self.RexScreenBox = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.RexScreenBox.SetValue("2")
        self.RexScreenBox.SetToolTip("Minimum fitted Rex (s^-1) required for Include=True. Press Enter to re-analyse.")
        self.RexScreenBox.Bind(wx.EVT_TEXT_ENTER, self.OnRexScreenEnter)

        self.TimeT2BoxTxt = wx.StaticText(self, label="Time_T2:")
        self.TimeT2Box = wx.TextCtrl(self,size=(150,-1),style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.TimeT2Box.SetValue("0.04")

        # CPMG fit metadata is explicit and independent of legacy CATIA tables.
        self.fieldBoxTxt=wx.StaticText(self,label="1H field (MHz):")
        self.fieldBox=wx.TextCtrl(self,value="750")
        self.nucBoxTxt=wx.StaticText(self,label="Observed nucleus:")
        self.nucBox=wx.ComboBox(self,choices=sorted(GAMMA.keys()),style=wx.CB_DROPDOWN)
        self.nucBox.SetValue("15N")
        self.tempBoxTxt=wx.StaticText(self,label="Temperature (C):")
        self.tempBox=wx.TextCtrl(self,value="20")
        self.buttonGuessMeta=wx.Button(self,label="Guess metadata")
        self.buttonGuessMeta.Bind(wx.EVT_BUTTON,self.OnGuessMetadata)

        self.pseudoAxisTxt = wx.StaticText(self, label="Pseudo-axis:")
        self.pseudoAxisCombo = wx.ComboBox(self, -1, choices=[], style=wx.CB_READONLY)
        self.pseudoAxisOpen = wx.Button(self, label="View")
        self.pseudoAxisCombo.Bind(wx.EVT_COMBOBOX, self.OnPseudoAxisColumn)
        self.pseudoAxisOpen.Bind(wx.EVT_BUTTON, self.OnOpenPseudoAxisTable)
        self._load_pseudo_axis_choices()

        

        
        self.basisCombo.Bind(wx.EVT_COMBOBOX, self.DoLocal)
        self.seqfilCombo.Bind(wx.EVT_COMBOBOX, self.DoLocal)
        

        self.parLocal = EditableListCtrl(self,style=wx.LC_REPORT)
        self.parLocal.SetMinSize((300,300))

        self.setLocal = EditableListCtrl(self,style=wx.LC_REPORT)
        self.setLocal.SetMinSize((300,300))

        self.fitLocal = AutoWidthListCtrl(self)
        self.fitLocal.SetMinSize((300,300))
        


        
        #self.setLocal.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClickSet)
        #self.parLocal.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClickPar)

        
        self.listbox=wx.BoxSizer(wx.HORIZONTAL)
        
        setbox = wx.StaticBox(self,-1,'Set Parameters:')
        setboxS=wx.StaticBoxSizer(setbox,wx.VERTICAL)
        setboxS.Add(self.setLocal)
        
        locbox = wx.StaticBox(self,-1,'Local Parameters:')
        locboxS=wx.StaticBoxSizer(locbox,wx.VERTICAL)
        locboxS.Add(self.parLocal)


        fitbox = wx.StaticBox(self,-1,'Fitted Parameters:')
        fitboxS=wx.StaticBoxSizer(fitbox,wx.VERTICAL)
        fitboxS.Add(self.fitLocal)

        self.listbox.AddSpacer(10)
        self.listbox.Add(setboxS)
        self.listbox.AddSpacer(10)
        self.listbox.Add(locboxS)
        self.listbox.AddSpacer(10)
        self.listbox.Add(fitboxS)
        self.listbox.AddSpacer(10)
        
        
        #self.parlocl.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnClick)
        #self.datasets.Bind(wx.EVT_LIST_ITEM_SELECTED, self.draw_figure)
        self.DoLocal(True)


        

        
        # Responsive CPMG analysis layout.  Every visible control is owned by
        # a sizer; legacy CATIA widgets remain available to old methods but are
        # deliberately hidden so they cannot appear as untethered frame children.
        for legacy in (self.buttonWipe, self.buttonRun, self.basistxt, self.seqfiltxt,
                       self.basisCombo, self.seqfilCombo, self.parLocal,
                       self.setLocal, self.fitLocal, setbox, locbox, fitbox):
            legacy.Hide()

        self.optionsPanel = wx.Panel(self)
        self.optionsPanel.SetMinSize((285, -1))
        self.vboxOpts = wx.BoxSizer(wx.VERTICAL)

        # The detailed dataset/fitting table lives in the modeless Peak Selection window.
        self.datasets.Hide()

        settings_static = wx.StaticBox(self.optionsPanel, label="CPMG settings")
        settings_box = wx.StaticBoxSizer(settings_static, wx.VERTICAL)
        for ctrl in (self.pseudoAxisTxt, self.pseudoAxisCombo, self.pseudoAxisOpen,
                     self.TimeT2BoxTxt, self.TimeT2Box, self.RexScreenBoxTxt,
                     self.RexScreenBox, self.fieldBoxTxt, self.fieldBox,
                     self.nucBoxTxt, self.nucBox, self.tempBoxTxt, self.tempBox,
                     self.buttonGuessMeta):
            ctrl.Reparent(settings_static)

        for ctrl in (self.peakChoiceTxt,self.peakChoice,self.peakSpin,self.screenCheck): ctrl.Reparent(settings_static)
        peak_row=wx.BoxSizer(wx.HORIZONTAL)
        peak_row.Add(self.peakChoiceTxt,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,8)
        peak_row.Add(self.peakChoice,1,wx.EXPAND|wx.RIGHT,4)
        peak_row.Add(self.peakSpin,0,wx.EXPAND|wx.RIGHT,8)
        peak_row.Add(self.screenCheck,0,wx.ALIGN_CENTER_VERTICAL)
        settings_box.Add(peak_row,0,wx.EXPAND|wx.ALL,6)

        axis_row = wx.BoxSizer(wx.HORIZONTAL)
        axis_row.Add(self.pseudoAxisTxt, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        axis_row.Add(self.pseudoAxisCombo, 1, wx.EXPAND | wx.RIGHT, 4)
        axis_row.Add(self.pseudoAxisOpen, 0)
        settings_box.Add(axis_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)

        def add_setting_row(label, control):
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.AddStretchSpacer(1)
            control.SetMinSize((105, -1))
            row.Add(control, 0, wx.ALIGN_CENTER_VERTICAL)
            settings_box.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)

        add_setting_row(self.TimeT2BoxTxt, self.TimeT2Box)
        add_setting_row(self.RexScreenBoxTxt, self.RexScreenBox)
        add_setting_row(self.fieldBoxTxt, self.fieldBox)
        add_setting_row(self.nucBoxTxt, self.nucBox)
        add_setting_row(self.tempBoxTxt, self.tempBox)
        settings_box.Add(self.buttonGuessMeta, 0, wx.EXPAND | wx.ALL, 6)
        self.vboxOpts.Add(settings_box, 0, wx.EXPAND | wx.BOTTOM, 8)

        plot_static = wx.StaticBox(self.optionsPanel, label="Plot")
        plot_box = wx.StaticBoxSizer(plot_static, wx.VERTICAL)
        self.plotModeTxt.Reparent(plot_static); self.plotMode.Reparent(plot_static)
        self.plotModeTxt.Hide()  # the StaticBox label already says Plot
        plot_box.Add(self.plotMode, 0, wx.EXPAND | wx.ALL, 6)
        self.vboxOpts.Add(plot_box, 0, wx.EXPAND | wx.BOTTOM, 8)

        actions_static = wx.StaticBox(self.optionsPanel, label="Actions")
        actions_box = wx.StaticBoxSizer(actions_static, wx.VERTICAL)
        action_grid = wx.GridSizer(rows=0, cols=2, vgap=4, hgap=4)
        for button in (self.buttonRefresh, self.buttonPeakConvert, self.buttonGuess, self.buttonPeaks, self.buttonAdvanced,
                       self.buttonClose):
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
        
        #self.create_main_panel()
        #self.draw_figure()
        #self.canvas.draw()
        
        self.SetSizer(self.sizerMain)
        self.Layout()  

        self.OnButtonLoad(True)
        
        #FGA added
        #if self.ref=='subgraphMode':
        #    self.rbox.Disable()

        self.Centre()
        self.Show(True)
        # Do not Fit() here: matplotlib's best size would force the frame and
        # defeat responsive resizing.
        self.Layout()
        #self.Thaw()




    def cpmg_Baldwin(self):
        self.cpmgF=baldwin_r2eff(self.cpmgX,self.Time_T2,self.pb,self.kex,self.R0,self.dw)
    def PackCPMG(self):
        x=[]
        x.append(self.pb)
        x.append(self.kex)
        x.append(self.R0)
        x.append(self.dw)
        return x
    def UnpackCPMG(self,x):
        x=numpy.fabs(x)
        cnt=0
        self.pb=x[cnt];cnt+=1
        self.kex=x[cnt];cnt+=1
        self.R0=x[cnt];cnt+=1
        self.dw=x[cnt];cnt+=1
    def ChiCPMG(self,x):
        self.UnpackCPMG(x)
        self.cpmg_Baldwin()
        return self.cpmgF-self.cpmgY


    def cpmg_Line(self):
        self.cpmgL=numpy.ones(len(self.cpmgX))*self.R0line
    def PackLine(self):
        x=[]
        x.append(self.R0line)
        return x
    def UnpackLine(self,x):
        x=numpy.fabs(x)
        cnt=0
        self.R0line=x[cnt];cnt+=1
    def ChiLine(self,x):
        self.UnpackLine(x)
        self.cpmg_Line()
        return self.cpmgL-self.cpmgY

    #https://en.wikipedia.org/wiki/Gyromagnetic_ratio
    #units of 10^6 rad s-1 T-1
    Gamma={}
    Gamma['15N']=27.116
    Gamma['1H']=267.522
    Gamma['19F']=251.815
    Gamma['13C']=67.2828
    Gamma['2H']=41.065
    Gamma['31P']=108.291
    
    def GetActualField(self):
        """Use explicit CPMG controls rather than legacy CATIA parameter lists."""
        self.field=float(self.fieldBox.GetValue())
        self.fieldLab=str(int(round(self.field)))
        self.nuc=self.nucBox.GetValue().strip()
        self.temp=float(self.tempBox.GetValue())
        self.dfrq=observe_frequency_mhz(self.field,self.nuc)
        return self.dfrq

    def OnGuessMetadata(self,event):
        """Best-effort Varian procpar discovery; all guessed values remain editable."""
        guesses={}
        roots=[self.acquisition_dir,getattr(self.parent,'fidpath','')]
        for root in roots:
            if not root or not os.path.exists(root): continue
            proc=None
            for base,dirs,files in os.walk(root):
                if 'procpar' in files: proc=os.path.join(base,'procpar'); break
            if not proc: continue
            try:
                lines=open(proc,errors='ignore').read().splitlines()
                for i,line in enumerate(lines[:-1]):
                    bits=line.split(); vals=lines[i+1].split()
                    if not bits or len(vals)<2: continue
                    key=bits[0]; val=vals[1].strip('"')
                    if key=='sfrq': guesses['field']=val
                    elif key in ('temp','temperature'): guesses['temp']=val
                    elif key=='tn': guesses['nuc']={'N15':'15N','C13':'13C','H1':'1H','F19':'19F','P31':'31P'}.get(val.upper(),val)
            except Exception: pass
            break
        if 'field' in guesses: self.fieldBox.SetValue(str(guesses['field']))
        if 'temp' in guesses: self.tempBox.SetValue(str(guesses['temp']))
        if guesses.get('nuc') in GAMMA: self.nucBox.SetValue(guesses['nuc'])
        msg='Metadata guessed from acquisition data' if guesses else 'No acquisition metadata found; enter values manually'
        if hasattr(self,'statusBar'): self.statusBar.SetStatusText(msg)
        else: wx.MessageBox(msg,'CPMG metadata',wx.OK|wx.ICON_INFORMATION)


    def DoFit(self,pk, update_screen=True):
        if self.ReadFuda(pk) is False: return False
        try:
            self.Time_T2=float(self.TimeT2Box.GetValue()); self.GetActualField()
            result=fit_cpmg_local(self.cpmgX,self.cpmgY,self.Time_T2,self.dfrq,self.cpmgE)
        except Exception as exc:
            self._set_dataset_failure(pk, 'Fit failed')
            self.statusBar.SetStatusText('Fit failed for %s: %s' % (pk, exc)); return False
        if not result.get('valid', result.get('success', False)):
            self.cpmgLocal.pop(pk, None)
            self._set_dataset_failure(pk, 'Fit failed')
            self.statusBar.SetStatusText('Fit failed for %s: %s' % (pk, result.get('message', 'invalid fit')))
            return False
        self.cpmgLocal[pk]=result
        self.pb=result['pb']; self.kex=result['kex']; self.R0=result['R0']; self.dw=result['dw']*self.dfrq
        self.R0line=result['R0line']; self.chi2Line=result['chi2Line']; self.chi2Local=result['chi2local']
        self.cpmgF=result['model']; self.cpmgL=numpy.ones(len(self.cpmgX))*self.R0line
        if update_screen: self._update_dataset_result(pk,result)
        return True

    @staticmethod
    def _format_fit_value(value, fmt):
        try:
            value=float(value)
            return (fmt % value) if numpy.isfinite(value) else '-'
        except (TypeError, ValueError):
            return '-'

    def _format_value_error(self,result,key,fmt='%.3g'):
        value=self._format_fit_value(result.get(key),fmt)
        error=self._format_fit_value(result.get(key+'_error'),fmt)
        return value if error == '-' else value+' ± '+error

    def _set_dataset_failure(self, pk, status='Fit failed'):
        """Mark a peak as excluded and clear stale fitted parameters."""
        for row in range(self.datasets.GetItemCount()):
            if self.datasets.GetItem(row,0).GetText()==pk:
                vals=['False','-','-','-','-','-','-','-',status]
                for col,val in enumerate(vals,1): self.datasets.SetItem(row,col,val)
                self._refresh_peak_selector(self.peakChoice.GetValue())
                break

    def _update_dataset_result(self,pk,result):
        threshold=self._rex_threshold()
        valid=bool(result.get('valid', result.get('success', False)))
        rex=result.get('Rex',numpy.nan)
        include=bool(valid and numpy.isfinite(rex) and rex>=threshold)
        for row in range(self.datasets.GetItemCount()):
            if self.datasets.GetItem(row,0).GetText()==pk:
                vals=[str(include),
                      self._format_fit_value(rex,'%.2f'),
                      self._format_fit_value(result.get('R2inf'),'%.2f'),
                      self._format_value_error(result,'kex','%.3g'),
                      self._format_value_error(result,'pb','%.4f'),
                      self._format_value_error(result,'dw','%.3f'),
                      self._format_value_error(result,'R0','%.2f'),
                      self._format_fit_value(result.get('improvement'),'%.2f'),
                      'OK' if valid else 'Fit failed']
                for col,val in enumerate(vals,1): self.datasets.SetItem(row,col,val)
                self._refresh_peak_selector(self.peakChoice.GetValue())
                break

    def _included_peaks(self):
        return [self.datasets.GetItem(r,0).GetText() for r in range(self.datasets.GetItemCount()) if self.datasets.GetItem(r,1).GetText()=='True']

    def OnPlotMode(self,event): self._plot_current_mode()

    def _plot_current_mode(self):
        mode=self.plotMode.GetStringSelection()
        if mode=='Peak':
            pk=self.get_selected_peak()
            if pk: self._plot_peak(pk)
        elif mode in ('kex','pb','R2inf','Δω'): self._plot_parameter_summary('dw' if mode=='Δω' else mode)
        else: self._plot_included_curves()

    def _plot_peak(self,pk):
        if pk not in self.cpmgLocal and not self.DoFit(pk): return
        if self.ReadFuda(pk) is False: return
        r=self.cpmgLocal[pk]; self.GetActualField()
        self.fig.clear(); ax=self.fig.add_subplot(111)
        ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.cpmgE,fmt='o',label='data')
        xx=numpy.linspace(float(numpy.min(self.cpmgX)),float(numpy.max(self.cpmgX)),200)
        yy=baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),r['pb'],r['kex'],r['R0'],r['dw']*self.dfrq)
        ax.plot(xx,yy,label='Local fit')
        gr=(self.cpmgGlobal or {}).get('peaks',{}).get(pk)
        if gr is not None:
            gy=baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),gr['pb'],gr['kex'],gr['R0'],gr['dw']*self.dfrq)
            ax.plot(xx,gy,label='Global',color='tab:orange',linewidth=2.2)
        ax.axhline(r['R0line'],linestyle='--',label='flat')
        title='%s   local kex=%.1f s-1   pb=%.4f   Δω=%.3f ppm'%(pk,r['kex'],r['pb'],r['dw'])
        if gr is not None: title+='   | global kex=%.1f   pb=%.4f'%(gr['kex'],gr['pb'])
        ax.set_title(title); ax.set_xlabel('nu_CPMG (Hz)'); ax.set_ylabel('R2,eff (s-1)'); ax.legend(); self.canvas.draw()

    def _plot_parameter_summary(self,par):
        peaks=[]; vals=[]
        for pk in self._visible_peaks():
            if pk in self.cpmgLocal and numpy.isfinite(self.cpmgLocal[pk].get(par,numpy.nan)):
                peaks.append(pk); vals.append(self.cpmgLocal[pk][par])
        self.fig.clear(); ax=self.fig.add_subplot(111); ax.plot(numpy.arange(len(vals)),vals,'o',label='Local')
        global_result=self.cpmgGlobal or {}
        if par in ('kex','pb') and global_result and numpy.isfinite(global_result.get(par,numpy.nan)):
            ax.axhline(global_result[par],color='tab:orange',linewidth=2,label='Global')
            ax.legend()
        ax.set_xticks(numpy.arange(len(peaks))); ax.set_xticklabels(peaks,rotation=90)
        label='Δω' if par=='dw' else par
        units=' (s-1)' if par in ('kex','R2inf') else (' (ppm)' if par=='dw' else '')
        ax.set_xlabel('Peak'); ax.set_ylabel(label+units); ax.set_title('%s versus peak'%label)
        self.fig.tight_layout(); self.canvas.draw()

    def _plot_included_curves(self):
        self.fig.clear(); ax=self.fig.add_subplot(111); plotted=0; self.GetActualField()
        for pk in self._visible_peaks():
            if pk not in self.cpmgLocal and not self.DoFit(pk,update_screen=False): continue
            if self.ReadFuda(pk) is False: continue
            r=self.cpmgLocal[pk]; ax.errorbar(self.cpmgX,self.cpmgY,yerr=self.cpmgE,fmt='o',markersize=3,alpha=.7)
            xx=numpy.linspace(float(numpy.min(self.cpmgX)),float(numpy.max(self.cpmgX)),150)
            yy=baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),r['pb'],r['kex'],r['R0'],r['dw']*self.dfrq)
            ax.plot(xx,yy,label=pk); plotted+=1
        ax.set_xlabel('nu_CPMG (Hz)'); ax.set_ylabel('R2,eff (s-1)'); ax.set_title('Visible CPMG dispersions (%d)'%plotted)
        if plotted and plotted<=20: ax.legend(fontsize='small',ncol=2)
        self.canvas.draw()

    def _fit_global_for_report(self, peaks):
        """Run the internal global fit without opening any result windows."""
        self.GetActualField(); curves={}
        for pk in peaks:
            if self.ReadFuda(pk):
                curves[pk]={'x':self.cpmgX.copy(),'y':self.cpmgY.copy(),'e':self.cpmgE.copy()}
        if not curves:
            return None
        result=fit_cpmg_global(curves,float(self.TimeT2Box.GetValue()),self.dfrq)
        if result.get('success',False):
            self.cpmgGlobal=result
            return result
        return None

    def report_results_rows(self):
        """Return Rex-screened CPMG results, ranked by local Rex."""
        threshold=self._rex_threshold()
        peaks=[]
        for pk,r in self.cpmgLocal.items():
            rex=r.get('Rex',numpy.nan)
            if r.get('valid',r.get('success',False)) and numpy.isfinite(rex) and rex >= threshold:
                peaks.append(pk)
        peaks.sort(key=lambda pk:self.cpmgLocal[pk].get('Rex',-numpy.inf), reverse=True)
        columns=['Peak','Rank','Rex','Local R0','Local R0 error','Local R2inf','Local dw','Local dw error','Local kex','Local kex error','Local pb','Local pb error',
                 'Local chi2','Local gain','Global R0','Global R0 error','Global R2inf','Global dw','Global dw error','Global chi2','Global gain']
        rows=[]
        gp=(self.cpmgGlobal or {}).get('peaks',{})
        for rank,pk in enumerate(peaks,1):
            local=self.cpmgLocal[pk]; glob=gp.get(pk,{})
            f=lambda d,k,fmt: self._format_fit_value(d.get(k,numpy.nan),fmt)
            rows.append([pk,str(rank),f(local,'Rex','%.3g'),f(local,'R0','%.3g'),f(local,'R0_error','%.2g'),f(local,'R2inf','%.3g'),
                         f(local,'dw','%.4g'),f(local,'dw_error','%.2g'),f(local,'kex','%.4g'),f(local,'kex_error','%.2g'),f(local,'pb','%.5g'),f(local,'pb_error','%.2g'),f(local,'chi2local','%.3g'),
                         f(local,'improvement','%.3g'),f(glob,'R0','%.3g'),f(glob,'R0_error','%.2g'),f(glob,'R2inf','%.3g'),f(glob,'dw','%.4g'),f(glob,'dw_error','%.2g'),
                         f(glob,'chi2local','%.3g'),f(glob,'improvement','%.3g')])
        return columns,rows,peaks

    def export_report_figures(self, report_dir):
        """Analyse CPMG data and export ranked significant local/global results."""
        os.makedirs(os.fspath(report_dir),exist_ok=True)
        self.OnAnalyse(None)
        columns,rows,peaks=self.report_results_rows()
        if peaks:
            global_peaks=set((self.cpmgGlobal or {}).get('peaks',{}))
            if not set(peaks).issubset(global_peaks):
                self._fit_global_for_report(peaks)
            columns,rows,peaks=self.report_results_rows()
        summary_figures=[]; peak_figures={}; old_size=tuple(self.fig.get_size_inches())
        try:
            if peaks:
                self.fig.clear(); ax=self.fig.add_subplot(111)
                rex=[self.cpmgLocal[p]['Rex'] for p in peaks]
                ax.bar(numpy.arange(len(peaks)),rex)
                ax.axhline(self._rex_threshold(),linestyle='--',label='Rex screen')
                ax.set_xticks(numpy.arange(len(peaks))); ax.set_xticklabels(peaks,rotation=90)
                ax.set_ylabel('Rex (s-1)'); ax.set_title('Significant CPMG peaks ranked by Rex'); ax.legend()
                self.fig.tight_layout(); filename='cpmg_rex_ranking.pdf'
                self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                summary_figures.append((filename,'Significant peaks ranked by Rex'))

                # Compact parameter summaries for the selected (Rex-significant) set.
                global_result=self.cpmgGlobal or {}
                gp=global_result.get('peaks',{})
                x=numpy.arange(len(peaks))
                for par,label,units in (('kex','kex','s-1'),('pb','pb',''),('dw','deltaOmega','ppm')):
                    self.fig.clear(); ax=self.fig.add_subplot(111)
                    local_vals=[self.cpmgLocal[p].get(par,numpy.nan) for p in peaks]
                    local_err=[self.cpmgLocal[p].get(par+'_error',numpy.nan) for p in peaks]
                    ax.errorbar(x,local_vals,yerr=local_err,fmt='o',label='Local')
                    if par in ('kex','pb') and numpy.isfinite(global_result.get(par,numpy.nan)):
                        gv=global_result[par]; ge=global_result.get(par+'_error',numpy.nan)
                        ax.axhline(gv,linewidth=2,label='Global')
                        if numpy.isfinite(ge): ax.axhspan(gv-ge,gv+ge,alpha=.15)
                    elif par=='dw':
                        gvals=[gp.get(p,{}).get('dw',numpy.nan) for p in peaks]
                        gerr=[gp.get(p,{}).get('dw_error',numpy.nan) for p in peaks]
                        ax.errorbar(x,gvals,yerr=gerr,fmt='s',label='Global')
                    ax.set_xticks(x); ax.set_xticklabels(peaks,rotation=90)
                    ax.set_ylabel(label + ((' ('+units+')') if units else '')); ax.set_title(label+' for significant CPMG peaks'); ax.legend()
                    self.fig.tight_layout(); filename='cpmg_summary_%s.pdf'%par
                    self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                    summary_figures.append((filename,label+' for selected peaks'))

                # Overlay the selected relaxation-dispersion data and both fitted models.
                self.fig.clear(); ax=self.fig.add_subplot(111); self.GetActualField()
                for pk in peaks:
                    if self.ReadFuda(pk) is False: continue
                    r=self.cpmgLocal[pk]; xx=numpy.linspace(float(numpy.min(self.cpmgX)),float(numpy.max(self.cpmgX)),150)
                    ax.plot(self.cpmgX,self.cpmgY,'o',markersize=2.5,alpha=.45)
                    ax.plot(xx,baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),r['pb'],r['kex'],r['R0'],r['dw']*self.dfrq),linewidth=1,alpha=.65)
                    gr=gp.get(pk)
                    if gr is not None:
                        ax.plot(xx,baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),gr['pb'],gr['kex'],gr['R0'],gr['dw']*self.dfrq),linewidth=1.8,alpha=.85)
                ax.set_xlabel('nu_CPMG (Hz)'); ax.set_ylabel('R2,eff (s-1)')
                ax.set_title('Selected relaxation-dispersion curves: local and global fits')
                self.fig.tight_layout(); filename='cpmg_summary_dispersion.pdf'
                self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                summary_figures.append((filename,'Relaxation-dispersion curves for selected peaks'))

                # Normalised overlay: remove the fitted R2,infinity baseline so
                # exchange-dependent dispersion amplitudes can be compared directly.
                self.fig.clear(); ax=self.fig.add_subplot(111); self.GetActualField()
                for pk in peaks:
                    if self.ReadFuda(pk) is False: continue
                    r=self.cpmgLocal[pk]; xx=numpy.linspace(float(numpy.min(self.cpmgX)),float(numpy.max(self.cpmgX)),150)
                    lr2inf=r.get('R2inf',numpy.nan)
                    if numpy.isfinite(lr2inf):
                        ax.plot(self.cpmgX,numpy.asarray(self.cpmgY)-lr2inf,'o',markersize=2.5,alpha=.45)
                        ly=baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),r['pb'],r['kex'],r['R0'],r['dw']*self.dfrq)
                        ax.plot(xx,ly-lr2inf,linewidth=1,alpha=.65)
                    gr=gp.get(pk)
                    if gr is not None and numpy.isfinite(gr.get('R2inf',numpy.nan)):
                        gy=baldwin_r2eff(xx,float(self.TimeT2Box.GetValue()),gr['pb'],gr['kex'],gr['R0'],gr['dw']*self.dfrq)
                        ax.plot(xx,gy-gr['R2inf'],linewidth=1.8,alpha=.85)
                ax.axhline(0,linestyle='--',linewidth=.8)
                ax.set_xlabel('nu_CPMG (Hz)'); ax.set_ylabel('R2,eff - R2,infinity (s-1)')
                ax.set_title('Selected relaxation-dispersion curves normalised to R2,infinity')
                self.fig.tight_layout(); filename='cpmg_summary_dispersion_normalised.pdf'
                self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                summary_figures.append((filename,'Selected relaxation-dispersion curves after subtraction of fitted R2,infinity'))

                # Peak-specific R0 comparison from the local and global fits.
                self.fig.clear(); ax=self.fig.add_subplot(111)
                local_r0=[self.cpmgLocal[p].get('R0',numpy.nan) for p in peaks]
                local_r0e=[self.cpmgLocal[p].get('R0_error',numpy.nan) for p in peaks]
                global_r0=[gp.get(p,{}).get('R0',numpy.nan) for p in peaks]
                global_r0e=[gp.get(p,{}).get('R0_error',numpy.nan) for p in peaks]
                ax.errorbar(x,local_r0,yerr=local_r0e,fmt='o',label='Local')
                ax.errorbar(x,global_r0,yerr=global_r0e,fmt='s',label='Global')
                ax.set_xticks(x); ax.set_xticklabels(peaks,rotation=90)
                ax.set_xlabel('Residue'); ax.set_ylabel('R0 (s-1)')
                ax.set_title('R0 versus residue for significant CPMG peaks'); ax.legend()
                self.fig.tight_layout(); filename='cpmg_summary_r0.pdf'
                self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                summary_figures.append((filename,'Local and global R0 values for selected residues'))
            self.fig.set_size_inches(4.5,4.0,forward=True)
            for i,pk in enumerate(peaks):
                self._plot_peak(pk); self.fig.set_size_inches(4.5,4.0,forward=True)
                filename='cpmg_peak_%03d.pdf'%(i+1)
                self.fig.savefig(os.path.join(os.fspath(report_dir),filename),bbox_inches='tight')
                peak_figures[pk]=filename
        finally:
            self.fig.set_size_inches(*old_size,forward=True)
        global_result=self.cpmgGlobal or {}
        return {'columns':columns,'rows':rows,'summary_figures':summary_figures,'peak_figures':peak_figures,
                'screen':{'Rex_threshold':self._rex_threshold(),'n_total':len(getattr(self,'peaks',[])),'n_significant':len(peaks)},
                'global':{'kex':global_result.get('kex'),'kex_error':global_result.get('kex_error'),'pb':global_result.get('pb'),'pb_error':global_result.get('pb_error'),'chi2':global_result.get('chi2'),
                          'success':bool(global_result.get('success',False)),'n_peaks':len(global_result.get('peaks',{}))}}

    def OnAdvancedCPMG(self,event):
        if self.advancedFrame and self.advancedFrame:
            try: self.advancedFrame.Raise(); return
            except Exception: self.advancedFrame=None
        self.advancedFrame=CPMGAdvancedDialog(self)
        self.advancedFrame.Bind(wx.EVT_CLOSE,self._on_advanced_close)
        self.advancedFrame.Show()

    def _on_advanced_close(self,event):
        self.advancedFrame=None
        event.Skip()

    def OnGlobalFitResults(self,event):
        if self.globalFitFrame:
            try:
                self.globalFitFrame.RefreshResults(); self.globalFitFrame.Raise(); return
            except Exception: self.globalFitFrame=None
        self.globalFitFrame=CPMGGlobalFitFrame(self); self.globalFitFrame.Show()

    def _refresh_global_fit_window(self):
        if self.globalFitFrame:
            try: self.globalFitFrame.RefreshResults()
            except Exception: self.globalFitFrame=None

    def OnPeakSelection(self,event):
        if self.peakSelectionFrame:
            try: self.peakSelectionFrame.Raise(); return
            except Exception: self.peakSelectionFrame=None
        self.peakSelectionFrame=CPMGPeakSelectionFrame(self)
        self.datasets.Show()
        self.peakSelectionFrame.Show()
        self._sync_table_selection(self.get_selected_peak())

    def _visible_peaks(self):
        peaks=sorted(getattr(self,'peaks',[]), key=peak_sort_key)
        if self.screenCheck.GetValue():
            included=set(self._included_peaks())
            peaks=[p for p in peaks if p in included]
        return peaks

    def _refresh_peak_selector(self,keep=None):
        if keep is None: keep=self.peakChoice.GetValue()
        peaks=self._visible_peaks()
        self.peakChoice.SetItems(peaks)
        if keep in peaks: self.peakChoice.SetValue(keep)
        elif peaks: self.peakChoice.SetSelection(0)
        else: self.peakChoice.SetValue('')
        self.peakSpin.Enable(len(peaks)>1)

    def get_selected_peak(self):
        pk=self.peakChoice.GetValue()
        if pk: return pk
        i=self.datasets.GetFirstSelected()
        return self.datasets.GetItem(i,0).GetText() if i>=0 else None

    def get_selected_peaks(self):
        selected=[]; i=self.datasets.GetFirstSelected()
        while i!=-1:
            selected.append(self.datasets.GetItem(i,0).GetText()); i=self.datasets.GetNextSelected(i)
        return selected or ([self.get_selected_peak()] if self.get_selected_peak() else [])

    def _sync_table_selection(self,pk):
        if not pk: return
        for r in range(self.datasets.GetItemCount()):
            if self.datasets.GetItem(r,0).GetText()==pk:
                self.datasets.Select(r); self.datasets.Focus(r); self.datasets.EnsureVisible(r); break

    def OnPeakChoice(self,event):
        pk=self.get_selected_peak(); self._sync_table_selection(pk)
        if self.plotMode.GetStringSelection()=='Peak' and pk: self._plot_peak(pk)

    def _step_peak(self,delta):
        peaks=self._visible_peaks()
        if not peaks: return
        current=self.peakChoice.GetValue()
        try: i=peaks.index(current)
        except ValueError: i=0
        i=max(0,min(len(peaks)-1,i+delta)); self.peakChoice.SetValue(peaks[i]); self.OnPeakChoice(None)
    def OnPeakUp(self,event): self._step_peak(-1)
    def OnPeakDown(self,event): self._step_peak(1)
    def OnScreenPeaks(self,event):
        keep=self.peakChoice.GetValue(); self._refresh_peak_selector(keep); self._plot_current_mode()


    def prepare_workflow_analysis(self):
        """Regenerate curves and fit when opened from the pseudo3D workflow."""
        self.OnAnalyse(None)
        self.Raise()

    def _rex_threshold(self):
        try:
            value=float(self.RexScreenBox.GetValue())
        except (TypeError, ValueError):
            raise ValueError('RexScreen must be a numeric value')
        if not numpy.isfinite(value) or value < 0:
            raise ValueError('RexScreen must be a finite value >= 0')
        return value

    def OnRexScreenEnter(self,event):
        self.OnAnalyse(event)

    def OnAnalyse(self,event):
        """Generate curves, fit all peaks, apply Rex screening, and redraw."""
        try:
            threshold=self._rex_threshold()
            self.Time_T2=float(self.TimeT2Box.GetValue())
            self.GetActualField()
        except Exception as exc:
            wx.MessageBox(str(exc), 'CPMG analysis', wx.OK|wx.ICON_WARNING)
            self.statusBar.SetStatusText('Analysis not run: %s' % exc)
            return

        peaks=sorted(getattr(self,'peaks',[]), key=peak_sort_key)
        if not peaks:
            self.statusBar.SetStatusText('No CPMG peaks available to analyse')
            return
        self.buttonGuess.Enable(False)
        failures=[]
        try:
            for i,pk in enumerate(peaks,1):
                self.statusBar.SetStatusText('Analysing %d/%d: %s' % (i,len(peaks),pk))
                wx.YieldIfNeeded()
                infile=os.path.join(self.fuda_dir, pk+'.out')
                if not self.MakeCPMGcurve(pk,infile):
                    self.cpmgLocal.pop(pk, None)
                    self._set_dataset_failure(pk, 'No curve')
                    failures.append(pk)
                    continue
                if not self.DoFit(pk):
                    failures.append(pk)
            self._sync_peak_selection_window(self.peakChoice.GetValue())
            included=len(self._included_peaks())
            self._plot_current_mode()
            message='Analysis complete: %d/%d peaks have Rex >= %.3g s^-1' % (included,len(peaks),threshold)
            if failures:
                message += ' (%d failed)' % len(failures)
            self.statusBar.SetStatusText(message)
        finally:
            self.buttonGuess.Enable(True)

    def OnButtonGuess(self,event):
        """Backward-compatible alias for the Analyse action."""
        return self.OnAnalyse(event)

    def _legacy_OnButtonGuess(self,event):
        #clever way to do this is to fit flat line, fit to my equation,
        #do FTest to check if fitting to dispersion led to something good.
        #loop over peaks

        self.Time_T2=float(self.TimeT2Box.GetValue())

        self.GetActualField() #set self.dfrq

        
        count = self.datasets.GetItemCount()
        col1 = [self.datasets.GetItem(row, 0).GetText() for row in range(count)]
        col2 = [self.datasets.GetItem(row, 1).GetText() for row in range(count)]
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            self.DoFit(c1)


        
                
            
        #print('Current:',col1,col2)
        #print('Changing inclusion')
        #num_items = self.datasets.GetItemCount()

        
        #1. read CPMG curve.
        #2. fit to baldwin eqn.
        #3. fit to flat line
        #compare chi2 values


        

    def OnButtonWipe(self,event):
        print('Wiping output folder')
        os.system('rm '+self.raw+'/catia/OutPut/*')

        
    def DoLocal(self,event):
        

        self.seqfil_cpmg=self.seqfilCombo.GetValue().split()[0]
        self.basis_cpmg=self.basisCombo.GetValue().split()[0]

        #self.field=self.fieldbox.GetValue().split()[0]
        self.GetActualField() #set self.dfrq field and fieldlab
        self.Time_T2=float(self.TimeT2Box.GetValue())
        xcar=118


        
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


        self.setLocal.ClearAll()
        self.setLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.setLocal.InsertColumn(1, 'Set', width = 80,format=wx.LIST_FORMAT_CENTRE)

        self.parLocal.ClearAll()        
        self.parLocal.InsertColumn(0, 'Parameter', width = 50,format=wx.LIST_FORMAT_CENTRE) 
        self.parLocal.InsertColumn(1, 'Initial', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.parLocal.InsertColumn(2, 'Fit?', width = 80,format=wx.LIST_FORMAT_CENTRE)
        self.parLocal.InsertColumn(3, 'File', width = 80,format=wx.LIST_FORMAT_CENTRE)

        num_items = self.parLocal.GetItemCount()
        cnt=0
        for local in self.local:
            cnt+=1
            self.parLocal.InsertStringItem(num_items,str(cnt))
            self.parLocal.SetStringItem(num_items,0,local[0])
            self.parLocal.SetStringItem(num_items,1,str(local[1]))
            self.parLocal.SetStringItem(num_items,2,str(local[2]))
            if(len(local)==4):
                self.parLocal.SetStringItem(num_items,3,str(local[3]))
            
        num_items = self.setLocal.GetItemCount()
        cnt=0
        for local in self.sett:
            cnt+=1
            self.setLocal.InsertStringItem(num_items,str(cnt))
            self.setLocal.SetStringItem(num_items,0,local[0])
            self.setLocal.SetStringItem(num_items,1,str(local[1]))
            
            
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
        
    def OnClick(self,event): # activation selects; Include is analysis-owned/read-only.
        sele=self.datasets.GetFirstSelected()
        if sele < 0:
            return
        peak=self.datasets.GetItem(sele,0).GetText()
        self.peakChoice.SetValue(peak)
        self._refresh_peak_selector(peak)
        self._plot_current_mode()

    def _dataset_headers(self):
        return [('Peak',65),('Include',60),('Rex',65),('R2inf',70),('kex ± err',115),('pb ± err',105),('dw ± err',105),('R0 ± err',105),('Fit gain',70),('Status',85)]

    def _refresh_dataset_table(self, preserve_peak=None):
        if preserve_peak is None:
            preserve_peak=self.get_selected_peak()
        self.datasets.ClearAll()
        headers=self._dataset_headers()
        for i,(name,width) in enumerate(headers):
            self.datasets.InsertColumn(i,name,width=width,format=wx.LIST_FORMAT_CENTRE)
        for row,peak in enumerate(sorted(getattr(self,'peaks',[]), key=peak_sort_key)):
            self.datasets.InsertItem(row,peak)
            result=self.cpmgLocal.get(peak)
            valid=bool(result and result.get('valid', result.get('success', False)))
            if valid:
                rex=result.get('Rex',numpy.nan)
                try: include=bool(numpy.isfinite(rex) and rex >= self._rex_threshold())
                except ValueError: include=False
                vals=[str(include), self._format_fit_value(rex,'%.2f'),
                      self._format_fit_value(result.get('R2inf'),'%.2f'),
                      self._format_value_error(result,'kex','%.3g'),
                      self._format_value_error(result,'pb','%.4f'),
                      self._format_value_error(result,'dw','%.3f'),
                      self._format_value_error(result,'R0','%.2f'),
                      self._format_fit_value(result.get('improvement'),'%.2f'), 'OK']
            else:
                vals=['False','-','-','-','-','-','-','-','Not analysed']
            for col,val in enumerate(vals,1): self.datasets.SetItem(row,col,val)
        self._refresh_peak_selector(preserve_peak)
        self._sync_table_selection(preserve_peak)
        if self.peakSelectionFrame:
            try:
                self.peakSelectionFrame.Layout(); self.datasets.Refresh(); self.peakSelectionFrame.Refresh()
            except Exception:
                pass

    def _sync_peak_selection_window(self, preserve_peak=None):
        self._refresh_dataset_table(preserve_peak)

    def OnButtonRefresh(self,event):
        files=os.listdir(self.fuda_dir)
        self.peaks=sorted([file[:-4] for file in files if file.endswith('.out')], key=peak_sort_key)
        self._refresh_dataset_table(self.get_selected_peak())

    def OnButtonSort(self,event):
        col=event.GetColumn()
        headers=self._dataset_headers()
        selected=self.get_selected_peak()
        if not hasattr(self,'_dataset_sort_column'):
            self._dataset_sort_column=None; self._dataset_sort_ascending=True
        if self._dataset_sort_column==col:
            self._dataset_sort_ascending=not self._dataset_sort_ascending
        else:
            self._dataset_sort_column=col; self._dataset_sort_ascending=True
        rows=[[self.datasets.GetItem(r,c).GetText() for c in range(len(headers))] for r in range(self.datasets.GetItemCount())]
        def key(row):
            value=row[col]
            if col==0: return peak_sort_key(value)
            if col==1: return (value!='True', value.lower())
            try: return (0,float(value))
            except Exception: return (1,float('inf'))
        rows.sort(key=key, reverse=not self._dataset_sort_ascending)
        self.datasets.ClearAll()
        for i,(name,width) in enumerate(headers): self.datasets.InsertColumn(i,name,width=width,format=wx.LIST_FORMAT_CENTRE)
        for r,row in enumerate(rows):
            self.datasets.InsertItem(r,row[0])
            for c,val in enumerate(row[1:],1): self.datasets.SetItem(r,c,val)
        self._sync_table_selection(selected)

    def OnButtonAddDataset(self,event):
        pass
    def OnButtonRemDataset(self,event):
        pass

    def PathExists(self,test):
        for t in test:
            if(os.path.exists(t)==False):
                print('Creating:',t)
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
                mark_complete(model='cpmg', source='catia')


    def OnButtonPeakConvert(self,event):
        print('Making CPMG curves....')
        count = self.datasets.GetItemCount()
        peaks = [self.datasets.GetItem(row, 0).GetText() for row in range(count)]
        print(peaks)

        for pk in peaks:
            print('Making curve for '+pk)
            cpmgfile=os.path.join(self.fuda_dir, pk + '.out.cpmg')
            if(os.path.exists(cpmgfile)):
                os.system('rm '+cpmgfile)
            self.ReadFuda(pk)
        

    def WriteCatiaDataset(self):

        print('Writing Catia dataset file')
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

        print(self.temp,self.fieldLab)
        outy=open(self.raw+'/catia/ParamSet_'+self.fieldLab+'_'+str(self.temp)+'.inp','w')
        self.locfile.append(self.raw+'/catia/ParamSet_'+self.fieldLab+'_'+str(self.temp)+'.inp')
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
                    self.source.SetStringItem(num_items,0,str(cnt))
                    self.source.SetStringItem(num_items,1,key)
                    self.source.SetStringItem(num_items,2,self.progress[i][0])

                    #num_items = self.source.GetItemCount()
                    ##FGA changed- depreciated functions
                    #self.source.InsertItem(str(cnt))
                    #self.source.SetItem(0,str(cnt))
                    #self.source.SetItem(1,key)
                    #self.source.SetItem(2,self.progress[i][0])


                    stry=vals
                    self.source.SetStringItem(num_items,3,stry)

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


    def LoadCatiaSettings(self):
        """Load only CATIA-specific state from the shared CPMG save file."""
        if not os.path.exists(self.savefile):
            return False
        data = read_structured_parameter_file(self.savefile)
        if not data:
            return False
        if 'seqfil' in data:
            self.seqfilCombo.SetValue(data['seqfil'])
        if 'basis' in data:
            self.basisCombo.SetValue(data['basis'])
        self.sett = [(key, value) for key, value in data.get('set', {}).items()]
        self.local = []
        for par_name, par_vals in data.get('par', {}).items():
            self.local.append([
                par_name,
                par_vals.get('1', ''),
                'fit' if par_vals.get('2', '') == 'fit' else '',
                '' if par_vals.get('3', '') == 'None' else par_vals.get('3', ''),
            ])
        self.SetLocal()
        return any(key in data for key in ('seqfil','basis','set','par'))

    def SaveCatiaSettings(self):
        """Persist CATIA state while preserving CPMG include/Rex state."""
        data = read_structured_parameter_file(self.savefile) if os.path.exists(self.savefile) else {}
        data['seqfil'] = self.seqfilCombo.GetValue()
        data['basis'] = self.basisCombo.GetValue()
        data['set'] = {}
        for row in range(self.setLocal.GetItemCount()):
            data['set'][self.setLocal.GetItem(row,0).GetText()] = self.setLocal.GetItem(row,1).GetText()
        data['par'] = {}
        for row in range(self.parLocal.GetItemCount()):
            name = self.parLocal.GetItem(row,0).GetText()
            value = self.parLocal.GetItem(row,1).GetText()
            fit = self.parLocal.GetItem(row,2).GetText() or 'fix'
            bound = self.parLocal.GetItem(row,3).GetText() or 'None'
            data['par'][name] = {1:value, 2:fit, 3:bound}
        write_structured_parameter_file(self.savefile, data)

    def OnButtonLoad(self,event):
        self.sett = []
        self.local = []

        if os.path.exists(self.savefile) == False:
            print('Cannot find savefile:', self.savefile)
            return

        data = read_structured_parameter_file(self.savefile)

        for key, value in data.get('set', {}).items():
            self.sett.append((key, value))

        if 'seqfil' in data:
            self.seqfilCombo.SetValue(data['seqfil'])
        if 'basis' in data:
            self.basisCombo.SetValue(data['basis'])
        if 'RexScreen' in data:
            self.RexScreenBox.SetValue(data['RexScreen'])

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

        print('local', self.local)

        self.SetLocal() #transfer pars to arrays

    def OnButtonSave(self,event):

        write={}
        
        write['seqfil']=self.seqfilCombo.GetValue()
        write['basis']=self.basisCombo.GetValue()
        write['RexScreen']=self.RexScreenBox.GetValue()
        
        count = self.datasets.GetItemCount()
        col1 = numpy.array([self.datasets.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.datasets.GetItem(row, 1).GetText() for row in range(count)])
        write['peak']={}
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            write['peak'][c1]=c2
        
        count = self.setLocal.GetItemCount()
        col1 = numpy.array([self.setLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.setLocal.GetItem(row, 1).GetText() for row in range(count)])
        write['set']={}
        for i,(c1,c2) in enumerate(zip(col1,col2)):
            write['set'][c1]=c2
        

        count = self.parLocal.GetItemCount()
        col1 = numpy.array([self.parLocal.GetItem(row, 0).GetText() for row in range(count)])
        col2 = numpy.array([self.parLocal.GetItem(row, 1).GetText() for row in range(count)])
        col3 = numpy.array([self.parLocal.GetItem(row, 2).GetText() for row in range(count)])
        col4 = numpy.array([self.parLocal.GetItem(row, 3).GetText() for row in range(count)])

        write['par']={}
        for i,(c1,c2,c3,c4) in enumerate(zip(col1,col2,col3,col4)):
            write['par'][c1]={}


            write['par'][c1][1]=c2

            if(c3==''):
                write['par'][c1][2]='fix'
            else:
                write['par'][c1][2]=c3

            if(c4==''):
                write['par'][c1][3]='None'
            else:
                write['par'][c1][3]=c4
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
            
            print('Writing results to file: ',path)
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


    def _analysis_spec_dir(self,inherit,dataset_root=''):
        obj=inherit
        while obj is not None:
            state=getattr(obj,'state',None)
            if state is not None and hasattr(state,'spec_dir'):
                try: return os.path.normpath(state.spec_dir())
                except Exception: pass
            obj=getattr(obj,'parent',None)
        return os.path.normpath(os.path.join(str(dataset_root or '.'),'spec'))

    def _acquisition_dir(self,inherit,dataset_root=''):
        obj=inherit
        while obj is not None:
            state=getattr(obj,'state',None)
            if state is not None and hasattr(state,'raw_dir'):
                try: return os.path.normpath(state.raw_dir())
                except Exception: pass
            obj=getattr(obj,'parent',None)
        return os.path.normpath(os.path.join(str(dataset_root or '.'),'raw'))

    def create_main_panel(self):
        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize((360, 280))
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        self.toolbar = NavigationToolbar(self.canvas)

        plot_static = wx.StaticBox(self, label="CPMG dispersion")
        plot_box = wx.StaticBoxSizer(plot_static, wx.VERTICAL)
        self.canvas.Reparent(plot_static)
        self.toolbar.Reparent(plot_static)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.TOP, 4)
        plot_box.Add(self.vbox, 1, wx.EXPAND | wx.ALL, 6)
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
        print(param,':',args)
        return numpy.array(args)


    
    def _load_pseudo_axis_choices(self):
        self.pseudo_axis_file = pseudo_axis_path(self.parent)
        try:
            table = PseudoAxisTable.load(self.pseudo_axis_file)
            self.pseudo_axis_table = table
            self.pseudoAxisCombo.SetItems(table.data_columns)
            selected = table.default_column(load_saved_column(self.parent))
            if selected:
                self.pseudoAxisCombo.SetStringSelection(selected)
        except PseudoAxisError:
            self.pseudo_axis_table = None
            self.pseudoAxisCombo.SetItems([])

    def OnPseudoAxisColumn(self, event):
        column = self.pseudoAxisCombo.GetValue().strip()
        if column:
            save_selected_column(self.parent, column)

    def OnOpenPseudoAxisTable(self, event):
        path = getattr(self, 'pseudo_axis_file', pseudo_axis_path(self.parent))
        if not os.path.exists(path):
            wx.MessageBox('Cannot find %s' % path, 'Pseudo-axis table', wx.OK | wx.ICON_WARNING)
            return
        try:
            table = getattr(self, 'pseudo_axis_table', None) or PseudoAxisTable.load(path)
            self.pseudo_axis_table = table
            show_pseudo_axis_table(self, table)
        except PseudoAxisError as exc:
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

    def Get_nuCPMG(self):
        """Return the selected conversion-defined pseudo coordinate."""
        try:
            return self._pseudo_axis_values()
        except PseudoAxisError as exc:
            wx.MessageBox(str(exc), 'Pseudo-axis', wx.OK | wx.ICON_WARNING)
            return numpy.array([-1.0])


    #turn fuda files into  .cpmg files
    def MakeCPMGcurve(self,pk,infile):
        """Create an R2eff curve from fitted pseudo3D intensities.

        File parsing remains a GUI/I/O concern; the scientific conversion is
        implemented in :mod:`decon.analysis.cpmg_service`.
        """
        try:
            time_t2=float(self.TimeT2Box.GetValue())
            axis=self._pseudo_axis_values()
            rows=[]
            with open(infile,'r') as handle:
                for line in handle:
                    if not line.strip() or line.lstrip().startswith('#'): continue
                    bits=line.split()
                    if len(bits)>=3: rows.append((float(bits[1]),float(bits[2])))
            if len(rows)!=len(axis):
                raise ValueError('Pseudo-axis has %d values; fitted data contain %d planes' % (len(axis),len(rows)))
            intensities=numpy.asarray([r[0] for r in rows],dtype=float)
            errors=numpy.asarray([r[1] for r in rows],dtype=float)
            self.cpmgX,self.cpmgY,self.cpmgE=build_r2eff(axis,intensities,time_t2,errors=errors)
            outfile=infile+'.cpmg'
            with open(outfile,'w') as ofn:
                ofn.write('#%11s%15s%13s\n' % ('nu_cpmg(Hz)','R2(1/s)','Esd(R2)'))
                for x,y,e in zip(self.cpmgX,self.cpmgY,self.cpmgE):
                    ofn.write(' %11.4e%15.6e%13.6e\n' % (x,y,e))
            self.statusBar.SetStatusText('Created CPMG curve for %s' % pk)
            return True
        except (ValueError,PseudoAxisError,OSError) as exc:
            self.statusBar.SetStatusText(str(exc))
            wx.MessageBox(str(exc),'CPMG curve',wx.OK|wx.ICON_WARNING)
            return False

    def ReadCPMGcurve(self,infile):
        if(os.path.exists(infile)==False):
            print ('no cpmg file:',infile)
            return
        print ('Reading CPMG file:',infile)
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
        print('looking for',infile)
        if(os.path.exists(infile)==False):
            print('Cannot find fuda intensities')
            return False
        cpmgfile=os.path.join(self.fuda_dir, pk + '.out.cpmg')
        if(os.path.exists(cpmgfile)==False):
            print('Cannot find CPMG intensities - creating')
            result=self.MakeCPMGcurve(pk,infile)
            return result is not False
        else:
            self.ReadCPMGcurve(cpmgfile)
            return True


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
                        self.fitLocal.SetStringItem(num_items,0,test[1]) #id

                        try:
                            self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (float(test[2]))) #A
                        except:
                            self.fitLocal.SetStringItem(num_items,1,'%s ' % ((test[2]))) #A


                        try:
                            self.fitLocal.SetStringItem(num_items,2,'%.3f ' % (float(test[3]))) #A
                        except:
                            self.fitLocal.SetStringItem(num_items,2,'%s ' % ((test[3]))) #A

                            
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
                #        self.fitLocal.SetStringItem(num_items,0,test[1]) #id
                    
                test=line.split()
                if(len(test)==3):
                    self.fitLocal.InsertStringItem(num_items,str(cnt))
                    self.fitLocal.SetStringItem(num_items,0,test[0]+' (global)') #id
                    try:
                        self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (float(test[1]))) #A
                    except:
                        self.fitLocal.SetStringItem(num_items,1,'%s ' % ((test[1]))) #A


                    try:
                        self.fitLocal.SetStringItem(num_items,2,'%.3f ' % (float(test[2]))) #A
                    except:
                        self.fitLocal.SetStringItem(num_items,2,'%s ' % ((test[2]))) #A

                            
                    cnt+=1


        if(pk not in self.cpmgLocal.keys()):
            self.DoFit(pk)
                    
        for par in 'R0line','pb','kex','R0','dw':
            self.fitLocal.InsertStringItem(num_items,str(cnt))
            self.fitLocal.SetStringItem(num_items,0,par+' (local)') #id
            val=self.cpmgLocal[pk][par]
            self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (val)) #A
            cnt+=1
            

        self.fitLocal.InsertStringItem(num_items,str(cnt))
        self.fitLocal.SetStringItem(num_items,0,'chi2Line') #id
        self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Line)) #A
        cnt+=1

        self.fitLocal.InsertStringItem(num_items,str(cnt))
        self.fitLocal.SetStringItem(num_items,0,'chi2Local') #id
        self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Local)) #A
        cnt+=1

        try:
            self.fitLocal.InsertStringItem(num_items,str(cnt))
            self.fitLocal.SetStringItem(num_items,0,'chi2Global') #id
            self.fitLocal.SetStringItem(num_items,1,'%.3f ' % (self.chi2Global)) #A
            cnt+=1
        except:
            pass

            
            
        return True
        
    def draw_figure(self,event):
        try:
            pk_event=self.datasets.GetItem(event.GetIndex(),0).GetText()
            if pk_event:
                self.peakChoice.SetValue(pk_event)
        except Exception:
            pass
        if self.plotMode.GetStringSelection()=='Peak':
            pk=self.get_selected_peak()
            if pk: self._plot_peak(pk)
        else:
            self._plot_current_mode()
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
            self.source.SetStringItem(num_items,0,str(noe[0])) #id
            self.source.SetStringItem(num_items,1,str(noe[1])) #A
            self.source.SetStringItem(num_items,2,str(noe[2])) #B
            self.source.SetStringItem(num_items,3,str(noe[7])) #intensity
            self.source.SetStringItem(num_items,4,str(noe[4])) # number
            self.source.SetStringItem(num_items,5,str(noe[6])) #w/s
            self.source.SetStringItem(num_items,6,'%.2f' % (noe[5]*1./(1.*self.nsoln)) ) #hits
            self.source.SetStringItem(num_items,7,str(noe[9])) #dist
            self.source.SetStringItem(num_items,8,str(noe[10])) #stdev

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
            #self.source.SetStringItem(num_items,3,stry)

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
            self.source.SetStringItem(num_items,0,str(cnt))
            self.source.SetStringItem(num_items,1,key)
            self.source.SetStringItem(num_items,2,'filter')
            stry=''
            for val in vals:
                stry+=val+' '
            self.source.SetStringItem(num_items,3,stry)

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
                        self.source.SetStringItem(num_items,0,str(cnt))
                        self.source.SetStringItem(num_items,1,key)
                        self.source.SetStringItem(num_items,2,self.progress[i][0])



                    if(shift=='n'):
                        stry=vals
                    else:
                        stry=vals+"(%.2f)" % (self.shiftDict[key][vals]['sc'])
                    if(self.WXV==4):
                        self.source.SetItem(num_items,3,stry)
                    else:
                        self.source.SetStringItem(num_items,3,stry)

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
                self.source.SetStringItem(num_items,0,str(cnt))
                self.source.SetStringItem(num_items,1,key)

 

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
                self.source.SetStringItem(num_items,3,stry)
            
            #self.source.ForceRefresh()
        #print len(self.result_dict.keys())
