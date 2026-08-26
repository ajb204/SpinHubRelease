#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
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
import wx
from spinDecon.domain.dimensions.viewer_contract import topology_for, spectral_dim_count
from spinDecon.gui.context import context_for, project_for, data_for
import string
import copy
import math
import numpy
import os
import sys
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar
from spinDecon.gui.widgets.common import PersistentStateButton

class _ToolbarToggleState:
    """Small wx.CheckBox-compatible state holder for toolbar-owned toggles."""
    def __init__(self, value=False):
        self._value = bool(value)
        self._enabled = True
    def GetValue(self): return self._value
    def IsChecked(self): return self._value
    def SetValue(self, value): self._value = bool(value)
    def Enable(self, enabled=True): self._enabled = bool(enabled)
    def IsEnabled(self): return self._enabled

import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Ellipse
from matplotlib.widgets import MultiCursor
import re
from matplotlib.gridspec import GridSpec

import scipy.optimize as opt
from spinDecon.domain.peaks import peakEntry
from spinDecon.processing.nmrpipe_scripts import MakeProj4D, MakeProj3D
from spinDecon.gui.plotting.array_utils import ensure_xy_points, scatter_xy_points
from spinDecon.gui.plotting.display_utils import blit_artists
# from spinDecon.line_fitting.line_fitting import Unidec_line_fitting
import threading
##################################################################################################################

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def Gauss(x,x0,Gamma):
    sigma=Gamma
    return numpy.exp(-(x-x0)**2./(2*sigma**2.))*1./(numpy.sqrt(2*numpy.pi)*sigma)

def Lorentz(x,x0,Gamma):
    return (Gamma/2.)/((x-x0)**2.+(Gamma/2.)**2.)/(numpy.pi)

def PV(amp,xx,x0,Gamma,sigma,nu):
    # print(amp, xx, x0, Gamma, sigma, nu)
    yvals=nu*Gauss(xx,x0,Gamma)+(1-nu)*(Lorentz(xx,x0,sigma))
    return amp*(yvals/numpy.max(yvals))


class FileDrop(wx.FileDropTarget):

    def __init__(self, canvas,axis, levels, parent):

        wx.FileDropTarget.__init__(self)
        self.canvas = canvas
        self.axis = axis
        self.ucs= []
        self.data = []
        self.first_drop = True
        self.extra_plots = []
        self.levels = levels
        self.parent = parent
        self.color_list = ['orange', 'yellow', 'green', 'lightblue', 'blue', 'purple']
        # self.data_1d = data_1d

    def OnDropFiles(self, x, y, filenames):

        for name in filenames:
            if '.ft' in name:
                # try:
                    dic, data = ng.pipe.read(name)
                    if len(data.shape) == 2:
                        if self.first_drop:
                            msg = "Entering multiple plot mode: Please enter title of the first dataset"
                            dlg = wx.TextEntryDialog(None, msg)
                            res = dlg.ShowModal()
                            if res == wx.ID_CANCEL:
                                return False
                            # self.data.set_label(dlg.GetValue())
                            self.first_drop = False


                        uc0= ng.pipe.make_uc(dic,data, dim=0)
                        uc1= ng.pipe.make_uc(dic,data, dim=1)


                        self.data.append(data)


                        x0,x1=uc0.ppm_limits()
                        uc0.ppms_scale=numpy.linspace(x0, x1, int(uc0._size))
                        x0,x1=uc1.ppm_limits()
                        uc1.ppms_scale=numpy.linspace(x0, x1, int(uc1._size))

                        self.ucs.append([uc0,uc1])
                        msg = "Please enter title of this data!"
                        dlg = wx.TextEntryDialog(None, msg)
                        res = dlg.ShowModal()
                        if res == wx.ID_CANCEL:
                            self.canvas.draw_idle()
                            return False
                        
                        self.extra_plots.append(self.axis.contour(uc1.ppm_scale(), uc0.ppm_scale(), data, colors=self.color_list[len(self.extra_plots)], linewidths=0.5, levels=self.levels))
                        self.axis.legend()
                        # self.axis.draw_artist(self.extra_plots[-1])
                        msg = "Do you want to move the data?"
                        dlg = wx.MessageDialog(None, msg, style=wx.YES_NO)
                        res = dlg.ShowModal()
                        if res == wx.ID_YES:
                            self.parent.drag_extra_plot = True
                            
                        self.canvas.draw_idle()
                        self.parent.background_save(None)
                        

                    else:
                        msg = "This is not 1D data - currently more dims are not supported..."
                        dlg = wx.MessageDialog(None, msg)
                        dlg.ShowModal()

                        return False

            else:

                msg = "Can only deal with *.ft* files!"
                dlg = wx.MessageDialog(None, msg)
                dlg.ShowModal()

                return False



        return True



class Projection(wx.Panel):
    """ The main frame of the application
    """
    title = 'Demo: wxPython with matplotlib'

    def __init__(self,parent,tabOne):

        wx.Panel.__init__(self,parent=parent)


        
        # Project tabs can be constructed immediately after a project-file load.
        # In that path legacy projects may not yet have restored the cached dmax
        # value even though the spectrum itself has already been read.  Projection
        # only needs an intensity scale, so derive it from the authoritative loaded
        # data rather than assuming the cache is populated.
        self.app_context = context_for(tabOne, parent)
        self.projection_service = self.app_context.projection if self.app_context is not None else None
        if self.projection_service is None:
            from spinDecon.analysis.projection_service import ProjectionService
            self.projection_service = ProjectionService(tabOne)
        self.thresh = self.projection_service.intensity_threshold()
        self.topology = topology_for(tabOne)
        self.spectral_dim_count = self.topology.spectral_dim_count
        self.physical_dim_count = self.topology.physical_dim_count
        self.parent=parent    #get decon_tab main parent from notebook
        self.state = getattr(tabOne, "state", getattr(parent, "state", None))
        self.store = getattr(tabOne, "store", None)
        pass

        self.sum='yz','CH'
        #self.pdbfile=tabOne.textbox.GetValue() #get the pdbfile from tabOne
        #self.methlist=tabOne.methList          #get the methyl list from tabOne
        self.shiftXrun=0                      #initialise the shiftX run counter
        self.selected = ''
        self.peaks_drawn=False
        self.resized = False # the dirty flag
        # True-2D marginal slice viewer state (mirrors Full3D).
        self.horizontal = False
        self.vertical = False
        self.trackers_locked = False
        self.last_mouse_x = None
        self.last_mouse_y = None
        self._trace_blit_background = None
        # Transient X markers used by the Pseudo3D/Fitting peak table.  These
        # are intentionally not part of the normal peak overlay and can be
        # updated cheaply with Matplotlib blitting.
        self._fitting_selection_background = None
        self._fitting_selection_artists = []
        self._fitting_selection_names = []
        # Pseudo2D Projection is the 1D reference editor.  The authoritative
        # list remains store.peak_lists['full']; this is GUI/history state only.
        self.full_tool_mode = None
        self.full_selected_name = None
        self.full_undo_stack = []
        self.full_redo_stack = []
        self.create_main_panel()
        if not self._is_pseudo2d_projection_case():
            self.load_decon_data()


        self.draw_figure()

        self.ftol = 1e-10
        
        self.dragging_extra_plot = False
        self.drag_extra_plot = False
        self.first_open = False
    
        # self.canvas.draw()
        if self._is_true_2d_spectrum():
            # Match Full3D: use the canvas efficiently while retaining ppm labels.
            self.fig.subplots_adjust(left=0.10, right=0.94, bottom=0.11, top=0.94)
        else:
            self.fig.subplots_adjust(left=0.043, right=0.981, top=0.975, bottom=0.098, wspace=0.186, hspace=0.2)
        self.canvas.draw()
        if self._is_true_2d_spectrum():
            self._save_2d_trace_background()

        self.Show(True)
        self.Fit()
        # self.background_save(None)




    def _is_pseudo2d_projection_case(self):
        return (self.topology.spectral_dim_count == 1 and self.topology.has_pseudo_axis and self.topology.physical_dim_count == 2)

    def _pseudo2d_projection(self, decon=False):
        # Source-contract compatibility: pseudo2d_projection_decon

        # Historical equivalent for the common pseudo-axis-first layout: numpy.sum(data, axis=0)
        if decon:
            data = self.projection_service.pseudo2d_decon_projection()
            if data is None:
                return None, None
            raw = self.projection_service.pseudo2d_data(ensure_file=False)
            return numpy.asarray(raw['index']), numpy.asarray(data).squeeze()
        payload = self.projection_service.pseudo2d_data(ensure_file=False)
        if payload is None:
            return None, None
        return numpy.asarray(payload['index']), numpy.asarray(payload['data']).squeeze()

    def _pseudo2d_peak_overlay(self):
        """Return (ppm, intensity, label) from the authoritative Full 1D list."""
        store = getattr(self, 'store', None)
        payload = store.peak_lists.get('full', {}) if store is not None else {}
        peaks = payload.get('records') or payload.get('peaks') or []
        points = []
        for entry in peaks:
            try:
                if isinstance(entry, dict):
                    coords = entry.get('coordinates', ())
                    ppm = coords[0] if coords else entry.get('ppm')
                    intensity = entry.get('intensity')
                    label = entry.get('name', entry.get('number', ''))
                    if ppm is not None and intensity is not None:
                        points.append((float(ppm), float(intensity), str(label)))
                elif hasattr(entry, 'ppm') and hasattr(entry, 'intensity'):
                    points.append((float(entry.ppm), float(entry.intensity), str(getattr(entry, 'number', ''))))
            except (IndexError, KeyError, TypeError, ValueError):
                continue
        return points

    # ---- Pseudo2D Full 1D peak-list tools ---------------------------------
    def _make_modeless_window(self, title):
        frame = wx.Frame(self.GetTopLevelParent(), title=title,
                         style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(frame)
        return frame, panel

    def _show_full_tools(self):
        frame = self.fullToolsFrame
        owner = wx.GetTopLevelParent(self)
        if owner is not None:
            pos, size = owner.GetScreenPosition(), owner.GetSize()
            x, y = pos.x + size.width, pos.y
            display_index = wx.Display.GetFromWindow(owner)
            if display_index != wx.NOT_FOUND:
                area = wx.Display(display_index).GetClientArea(); tool_size = frame.GetSize()
                x = max(area.GetLeft(), min(x, area.GetRight() - tool_size.width + 1))
                y = max(area.GetTop(), min(y, area.GetBottom() - tool_size.height + 1))
            frame.SetPosition((x, y))
        if not frame.IsShown(): frame.Show()
        frame.Raise()

    def full_peak_tools_box(self):
        """Slice2D-style editor for the authoritative pseudo2D Full 1D list."""
        self.fullToolsFrame, panel = self._make_modeless_window('Tools')
        self.fullUndoButton = wx.Button(panel, label='Undo', size=(-1, 24))
        self.fullRedoButton = wx.Button(panel, label='Redo', size=(-1, 24))
        self.fullSelectButton = PersistentStateButton(panel, label='Select', size=(-1, 24))
        self.fullMoveButton = PersistentStateButton(panel, label='Move', size=(-1, 24))
        self.fullAddButton = PersistentStateButton(panel, label='Add', size=(-1, 24))
        self.fullMaxButton = wx.Button(panel, label='Maximise', size=(-1, 24))
        self.fullRemoveButton = wx.Button(panel, label='Remove', size=(-1, 24))
        for button, handler in ((self.fullUndoButton,self.on_full_undo),(self.fullRedoButton,self.on_full_redo),
                                (self.fullSelectButton,self.on_full_select),(self.fullMoveButton,self.on_full_move),
                                (self.fullAddButton,self.on_full_add),(self.fullMaxButton,self.on_full_maximise),
                                (self.fullRemoveButton,self.on_full_remove)):
            button.Bind(wx.EVT_BUTTON, handler)
        close = wx.Button(panel, label='Close', size=(-1, 24)); close.Bind(wx.EVT_BUTTON, self.on_full_tools_close_button)
        self.fullToolsFrame.Bind(wx.EVT_CLOSE, self.on_full_tools_close)
        root = wx.BoxSizer(wx.VERTICAL)
        # Keep the same command order requested for the pseudo2D reference editor.
        for widget in (self.fullUndoButton,self.fullRedoButton,self.fullSelectButton,self.fullMoveButton,
                       self.fullAddButton,self.fullMaxButton,self.fullRemoveButton,close):
            widget.SetMinSize((92,24)); root.Add(widget,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP,4)
        root.AddSpacer(4); panel.SetSizer(root); self.fullToolsFrame.SetClientSize(root.CalcMin())
        self._update_full_tool_controls()

    def _toolbar_tools(self, active):
        if active:
            self._show_full_tools(); self.toolbar.set_tools_active(True)
        else: self._hide_full_tools()

    def _hide_full_tools(self):
        if hasattr(self, 'fullToolsFrame'): self.fullToolsFrame.Hide()
        if hasattr(self, 'toolbar'): self.toolbar.set_tools_active(False)

    def on_full_tools_close_button(self, event=None): self._hide_full_tools()
    def on_full_tools_close(self, event): self._hide_full_tools(); event.Veto()

    def _full_payload(self):
        return self.projection_service.full_peak_payload()
    def _full_records(self): return copy.deepcopy(list(self._full_payload().get('records') or []))
    def _full_record(self, name=None):
        wanted = self.full_selected_name if name is None else name
        return next((r for r in self._full_records() if str(r.get('name','')) == str(wanted)), None)

    def _set_full_status(self, text):
        frame=wx.GetTopLevelParent(self); statusbar=getattr(frame,'statusbar',None) if frame else None
        if statusbar is not None: statusbar.SetStatusText(str(text))

    def _update_full_history_buttons(self):
        if hasattr(self,'fullUndoButton'):
            self.fullUndoButton.Enable(bool(self.full_undo_stack)); self.fullRedoButton.Enable(bool(self.full_redo_stack))

    def _update_full_tool_controls(self):
        if not hasattr(self,'fullSelectButton'): return
        selected=self._full_record() is not None; modal=self.full_tool_mode not in (None,'select')
        self.fullSelectButton.SetActive(self.full_tool_mode=='select' or (selected and self.full_tool_mode is None))
        self.fullMoveButton.SetActive(self.full_tool_mode=='move'); self.fullAddButton.SetActive(self.full_tool_mode=='add')
        self.fullSelectButton.Enable(not modal)
        self.fullMoveButton.Enable(((not modal) or self.full_tool_mode=='move') and selected)
        # Add is always available when another one-shot operation is not active.
        # Entering Add clears the current selection, matching the natural 1D
        # reference-editor workflow.
        self.fullAddButton.Enable((not modal) or self.full_tool_mode=='add')
        self.fullRemoveButton.Enable((not modal) and selected); self.fullMaxButton.Enable((not modal) and selected)
        for b in (self.fullSelectButton,self.fullMoveButton,self.fullAddButton): b.Refresh()
        self._update_full_history_buttons()

    def _cancel_projection_navigation(self):
        """Leave Matplotlib Pan/Zoom before using pseudo2D peak tools.

        NavigationToolbar2 keeps pan/zoom as a persistent mode.  The Full 1D
        editor deliberately ignores canvas clicks while either mode is active,
        so every editor command and Redraw must explicitly release that mode.
        """
        toolbar = getattr(self, 'toolbar', None)
        if toolbar is None:
            return
        mode = getattr(toolbar, 'mode', '')
        name = str(getattr(mode, 'name', '') or '').upper()
        text = str(mode or '').lower()
        try:
            if name == 'PAN' or 'pan/zoom' in text or text == 'pan':
                toolbar.pan()
            elif name == 'ZOOM' or 'zoom' in text:
                toolbar.zoom()
        except Exception:
            # Older Matplotlib/wx combinations can expose only the wx check
            # state.  Clear that state as a fallback; normal toolbar operation
            # remains untouched when no navigation mode is active.
            for label in ('Pan', 'Zoom'):
                tool_id = getattr(toolbar, 'wx_ids', {}).get(label)
                if tool_id is not None:
                    toolbar.ToggleTool(tool_id, False)
            try:
                toolbar.mode = ''
            except Exception:
                pass

    def _set_full_tool_mode(self, mode=None):
        if mode=='move' and self._full_record() is None:
            self._set_full_status('Move: select one Full peak first'); self._update_full_tool_controls(); return
        self.full_tool_mode=mode; self._update_full_tool_controls()
        if mode: self._set_full_status('%s: click the 1D projection' % mode.capitalize())

    def _clear_full_selection(self, redraw=True):
        self.full_selected_name=None; self._update_full_tool_controls()
        if redraw: self.draw_figure(keepaxes=True)

    def on_full_select(self,event=None):
        self._cancel_projection_navigation()
        if self.fullSelectButton.IsActive():
            self.full_tool_mode=None; self._clear_full_selection(True); self._set_full_status('Select: off')
        else:
            self._set_full_tool_mode('select'); self._set_full_status('Select: click the nearest Full 1D peak')
    def on_full_move(self,event=None):
        self._cancel_projection_navigation()
        self._set_full_tool_mode(None if self.fullMoveButton.IsActive() else 'move')
    def on_full_add(self,event=None):
        self._cancel_projection_navigation()
        if self.fullAddButton.IsActive():
            self._set_full_tool_mode(None)
            return
        # Adding is independent of the currently selected peak.
        self.full_selected_name = None
        self._set_full_tool_mode('add')

    def _full_snapshot(self): return copy.deepcopy(self._full_records())
    def _push_full_undo(self):
        self.full_undo_stack.append(self._full_snapshot()); self.full_redo_stack=[]; self._update_full_history_buttons()

    def _normalise_full_1d_records(self, records):
        rows=[]
        for i,r in enumerate(records):
            coords=tuple(float(v) for v in r.get('coordinates',())[:1])
            if not coords: continue
            r['coordinates']=coords; r['row_index']=i
            label = self.projection_service.pseudo2d_data(ensure_file=False).get('label','Direct')
            r['axis_values']={str(label):coords[0]}
            intensity=float(r.get('intensity') or 0.0)
            fields=list(r.get('fields') or [])
            while len(fields)<3: fields.append('0')
            fields[0]=str(r.get('name',fields[0])); fields[1]=str(coords[0]); fields[2]=str(intensity)
            r['fields']=fields[:3]; rows.append(list(r['fields']))
        return rows

    def _commit_full_records(self, records):
        records=copy.deepcopy(list(records)); rows=self._normalise_full_1d_records(records); old=self._full_payload()
        # Source-contract compatibility: self.tabOne.store.save_peak_list('full' and refresh_full_peak_list_viewers
        self.projection_service.save_full_peak_list(peaks=records, records=records, rows=rows,
                                                    dimension=1, source_path=old.get('source_path'))
        self.draw_figure(keepaxes=True); self._update_full_tool_controls()

    def on_full_undo(self,event=None):
        self._cancel_projection_navigation()
        if not self.full_undo_stack: return
        self.full_redo_stack.append(self._full_snapshot()); snap=self.full_undo_stack.pop(); self._commit_full_records(snap)
        if self._full_record() is None: self.full_selected_name=None
        self._update_full_tool_controls(); self._set_full_status('Undo')
    def on_full_redo(self,event=None):
        self._cancel_projection_navigation()
        if not self.full_redo_stack: return
        self.full_undo_stack.append(self._full_snapshot()); snap=self.full_redo_stack.pop(); self._commit_full_records(snap)
        if self._full_record() is None: self.full_selected_name=None
        self._update_full_tool_controls(); self._set_full_status('Redo')

    def _projection_point(self, ppm):
        x,y=self._pseudo2d_projection(False)
        if x is None or y is None or len(x)==0: return float(ppm),0.0,0
        i=int(numpy.nanargmin(numpy.abs(numpy.asarray(x,dtype=float)-float(ppm))))
        return float(x[i]),float(y[i]),i

    def _select_full_at(self, ppm):
        records=self._full_records()
        if not records: self._set_full_status('No Full 1D peaks'); return
        span=abs(self.axes1.get_xlim()[1]-self.axes1.get_xlim()[0]) or 1.0
        record=min(records,key=lambda r: abs(float(r.get('coordinates',(1e99,))[0])-float(ppm))/span)
        self.full_selected_name=str(record.get('name')); self.full_tool_mode='select'
        # Keep the Full 1D Peak List viewer synchronised with projection selection.
        self.projection_service.focus_full_peak_list_viewers(self.full_selected_name)
        self._set_full_status('Selected: %s' % self.full_selected_name); self._update_full_tool_controls(); self.draw_figure(keepaxes=True)

    def _move_full_at(self, ppm):
        if self._full_record() is None: return
        x,intensity,_=self._projection_point(ppm); self._push_full_undo(); records=self._full_records()
        for r in records:
            if str(r.get('name'))==str(self.full_selected_name): r['coordinates']=(x,); r['intensity']=intensity; break
        self.full_tool_mode=None; self._commit_full_records(records); self._set_full_status('Moved: %s' % self.full_selected_name)

    def _next_full_1d_name(self):
        names=[str(r.get('name','')) for r in self._full_records()]
        nums=[int(n) for n in names if n.isdigit()]
        if nums: return str(max(nums)+1)
        nums=[int(m.group(1)) for n in names for m in [re.match(r'^Peak_(\d+)$',n)] if m]
        return 'Peak_%d' % ((max(nums)+1) if nums else 1)

    def _add_full_at(self, ppm):
        x,intensity,_=self._projection_point(ppm); self._push_full_undo(); records=self._full_records(); name=self._next_full_1d_name()
        label=self.projection_service.pseudo2d_data(ensure_file=False).get('label','Direct')
        records.append({'name':name,'coordinates':(x,),'axis_values':{str(label):x},'intensity':intensity,
                        'row_index':len(records),'fields':[name,str(x),str(intensity)],'analysis':{}})
        self.full_selected_name=name; self.full_tool_mode=None; self._commit_full_records(records); self._set_full_status('Added: %s' % name)

    def on_full_remove(self,event=None):
        self._cancel_projection_navigation()
        if self._full_record() is None: self._set_full_status('Remove: select a Full peak first'); return
        name=self.full_selected_name; self._push_full_undo(); self.full_selected_name=None; self.full_tool_mode=None
        self._commit_full_records([r for r in self._full_records() if str(r.get('name'))!=str(name)]); self._set_full_status('Removed: %s' % name)

    def on_full_maximise(self,event=None):
        self._cancel_projection_navigation()
        record=self._full_record()
        if record is None: self._set_full_status('Maximise: select a Full peak first'); return
        x,y=self._pseudo2d_projection(False)
        if x is None or y is None or not len(x): return
        arr=numpy.asarray(y,dtype=float); xv=numpy.asarray(x,dtype=float); ppm=float(record['coordinates'][0])
        i=int(numpy.nanargmin(numpy.abs(xv-ppm)))
        # Slice2D-style local hill climb, reduced deliberately to one spectral dimension.
        for _ in range(30):
            lo,hi=max(0,i-1),min(len(arr),i+2); ni=lo+int(numpy.nanargmax(arr[lo:hi]))
            if ni==i: break
            i=ni
        self._push_full_undo(); records=self._full_records()
        for r in records:
            if str(r.get('name'))==str(self.full_selected_name): r['coordinates']=(float(xv[i]),); r['intensity']=float(arr[i]); break
        self._commit_full_records(records); self._set_full_status('Maximised: %s' % self.full_selected_name)

    def _handle_full_tool_click(self,event):
        if self.full_tool_mode is None or event.inaxes is not self.axes1 or event.xdata is None or getattr(event,'button',None)!=1: return False
        # Never consume clicks intended for Matplotlib pan/zoom.
        mode=getattr(getattr(self,'toolbar',None),'mode','')
        if (isinstance(mode,str) and mode) or (getattr(mode,'name',None) not in (None,'NONE')): return False
        if self.full_tool_mode=='select': self._select_full_at(event.xdata)
        elif self.full_tool_mode=='move': self._move_full_at(event.xdata)
        elif self.full_tool_mode=='add': self._add_full_at(event.xdata)
        else: return False
        return True

    def select_full_peak_from_list(self, peak_name, zoom=True, row_index=None):
        """Select a Full 1D peak from the list viewer and optionally zoom to it.

        Pseudo2D Full peaks have exactly one spectral coordinate.  The green
        marker is therefore an X-only selection marker.
        """
        records = self._full_records()
        record = None
        if row_index is not None:
            record = next((r for r in records if r.get('row_index') == row_index), None)
        if record is None:
            record = next((r for r in records
                           if str(r.get('name', '')) == str(peak_name)), None)
        if record is None:
            return False
        try:
            ppm = float(record.get('coordinates', ())[0])
        except (IndexError, TypeError, ValueError):
            return False
        self.full_selected_name = str(record.get('name', peak_name))
        self.full_tool_mode = None
        self.draw_figure(keepaxes=True)
        if zoom:
            x, _ = self._pseudo2d_projection(False)
            xv = numpy.asarray(x, dtype=float) if x is not None else numpy.asarray([])
            if xv.size:
                width = abs(float(numpy.nanmax(xv)) - float(numpy.nanmin(xv)))
                # Slice2D-like focused view: enough context to see the local peak
                # while remaining useful for neighbouring resonances.
                half = max(width * 0.025, abs(float(numpy.nanmedian(numpy.diff(xv)))) * 20.0 if xv.size > 1 else 0.05)
                self.axes1.set_xlim(ppm + half, ppm - half)
                # Autoscale Y from the raw points visible in the focused window.
                raw = numpy.asarray(self._pseudo2d_projection(False)[1], dtype=float)
                mask = (xv >= ppm-half) & (xv <= ppm+half)
                if raw.size == xv.size and numpy.any(mask):
                    vals = raw[mask]
                    lo, hi = float(numpy.nanmin(vals)), float(numpy.nanmax(vals))
                    pad = (hi-lo)*0.08 or max(abs(hi), 1.0)*0.08
                    self.axes1.set_ylim(lo-pad, hi+pad)
                self.canvas.draw_idle()
        self._update_full_tool_controls()
        try:
            notebook = self.GetParent()
            idx = notebook.FindPage(self)
            if idx != wx.NOT_FOUND:
                notebook.SetSelection(idx)
        except Exception:
            pass
        self.SetFocus()
        return True

    def _spectrum_view(self, decon=False, transpose='n'):
        """Return the centrally owned true-2D plotting view."""
        return self.projection_service.spectrum_view(decon=decon, transpose=transpose)

    @property
    def twod_data(self):
        view = self._spectrum_view(False)
        return None if view is None else view.get('ZZ')

    @property
    def twod_data_decon(self):
        view = self._spectrum_view(True)
        return None if view is None else view.get('ZZ')

    @property
    def uc0(self):
        view = self._spectrum_view(False)
        return None if view is None else view.get('uc0')

    @property
    def uc1(self):
        view = self._spectrum_view(False)
        return None if view is None else view.get('uc1')

    def _cached_decon_plane(self, plane):
        store = getattr(self, "store", None)
        if store is None:
            return None
        return store.projections.get(("decon_plane", plane))

    def _cached_peak_overlay(self):
        """Return authoritative Full Peak List entries for overlaying.

        Connectivity (legacy ``conn_data``) is deliberately excluded: it is a
        relationship layer, not an alternative peak collection.
        """
        payload = self.projection_service.full_peak_payload()
        if payload:
            if isinstance(payload, dict):
                peaks = payload.get("peak") or payload.get("peaks")
                if peaks:
                    return peaks
            elif payload:
                return payload
        store = getattr(self, "store", None)
        if store is not None:
            full = store.peak_lists.get("full", {})
            peaks = full.get("peak") or full.get("peaks")
            if peaks:
                return peaks
            peaks = getattr(store, "peak", None)
            if peaks:
                return peaks
        return []

    def _peak_points_for_overlay(self, entries, swap=False):
        """Convert peak/connection objects into x/y scatter points."""
        points = []
        labels = []
        for entry in entries or []:
            if hasattr(entry, "f1") and hasattr(entry, "f2"):
                x = float(entry.f2 if swap else entry.f1)
                y = float(entry.f1 if swap else entry.f2)
                # Prefer the raw correlate label (first column) when present.
                label = getattr(entry, "p1", "") or getattr(entry, "tag", "")
            elif hasattr(entry, "ppmI") and hasattr(entry, "ppmJ"):
                x = float(entry.ppmJ if swap else entry.ppmI)
                y = float(entry.ppmI if swap else entry.ppmJ)
                label = getattr(entry, "name", "")
            elif hasattr(entry, "x") and hasattr(entry, "y"):
                x = float(entry.y if swap else entry.x)
                y = float(entry.x if swap else entry.y)
                label = getattr(entry, "name", "")
            elif isinstance(entry, (tuple, list)) and len(entry) >= 2:
                x = float(entry[1] if swap else entry[0])
                y = float(entry[0] if swap else entry[1])
                label = ""
            else:
                continue
            points.append([x, y])
            labels.append(label)
        return ensure_xy_points(points), labels

    def _projection_view(self, l1, l2, transpose='n', decon=False):
        """Retrieve a plotting-ready projection view from the main controller."""
        return self.projection_service.projection_view(l1, l2, decon=decon, transpose=transpose)

    def _projected_peak_overlay(self, l1, l2, transpose='n'):
        """Retrieve GUI-ready projected peak markers from the controller.

        Projection.py deliberately does not interpret f1/f2/f3 coordinates;
        the main controller/DataStore owns that mapping alongside spectrum
        projection views.
        """
        return self.projection_service.projected_peak_overlay(l1, l2, transpose=transpose) or []

    @staticmethod
    def _draw_projected_peak_overlay(axes, peaks, labels=False, size=20):
        for peak in peaks or []:
            try:
                x = float(peak['x'])
                y = float(peak['y'])
            except (KeyError, TypeError, ValueError):
                continue
            label = str(peak.get('label', '') or '')
            color = str(peak.get('color', '#000000'))
            if labels and label:
                axes.text(x, y, label, fontsize=8, color=color)
            axes.scatter(x, y, marker='x', c=color, s=size, zorder=2, linewidth=1.0)

    _REAL_AXIS_LABELS = ('time_T2', 'ID', 'ncyc', 'ncyc_cp', 'gzlvl5', 'gzlvl1')

    def _is_3p_projection_case(self):
        """True for a physical 3D cube containing two spectral axes and one real axis."""
        # Canonical Stage-6 test: scientific topology, not legacy dim==3 or
        # label guessing, defines a physical pseudo-3D dataset.
        if not (self.topology.spectral_dim_count == 2 and
                self.topology.has_pseudo_axis and
                self.topology.physical_dim_count == 3):
            return False
        data = self.projection_service.data
        return getattr(data, 'ndim', len(getattr(data, 'shape', ()))) == 3

    def _is_true_2d_spectrum(self):
        """True only for two spectral axes with no pseudo/real axis."""
        return (self.topology.spectral_dim_count == 2 and
                not self.topology.has_pseudo_axis and
                self.topology.physical_dim_count == 2)

    def _spectral_labels_3p(self):
        """Return (X/direct, Y/indirect) labels for the sole 3p spectral projection."""
        labb = self.projection_service.labels[:self.topology.physical_dim_count]
        spectral = [str(labb[a.physical_index]) for a in self.topology.spectral_axes
                    if a.physical_index < len(labb)]
        if len(spectral) != self.topology.spectral_dim_count:
            return None
        # labb follows physical numpy axes.  The last spectral axis is the
        # direct/X dimension and the preceding spectral axis is indirect/Y.
        return spectral[-1], spectral[-2]

    def _projection_view_3p(self, decon=False):
        labels = self._spectral_labels_3p()
        if labels is None:
            return None
        return self._projection_view(labels[0], labels[1], transpose='n', decon=decon)

    def _projection_views_3d(self, decon=False):
        if self._is_3p_projection_case():
            view = self._projection_view_3p(decon=decon)
            return [view] if view is not None else None
        if len(self.projection_service.labels) < 3:
            return None
        specs = [
            (self.projection_service.labels[2], self.projection_service.labels[1], 'n'),
            (self.projection_service.labels[2], self.projection_service.labels[0], 'n'),
            (self.projection_service.labels[1], self.projection_service.labels[0], 'y'),
        ]
        views = [self._projection_view(a, b, transpose=t, decon=decon) for a, b, t in specs]
        return views if all(view is not None for view in views) else None

    def _projection_views_4d(self):
        if len(self.projection_service.labels) < 4:
            return None
        specs = [
            (self.projection_service.labels[2], self.projection_service.labels[3], 'n'),
            (self.projection_service.labels[0], self.projection_service.labels[1], 'y'),
            (self.projection_service.labels[1], self.projection_service.labels[3], 'n'),
            (self.projection_service.labels[0], self.projection_service.labels[2], 'n'),
        ]
        views = [self._projection_view(a, b, transpose=t, decon=False) for a, b, t in specs]
        return views if all(view is not None for view in views) else None

    @staticmethod
    def _unpack_projection_view(view):
        return view['XX'], view['YY'], view['ZZ'], tuple(view['labb'])

    def _ensure_3d_decon_views(self, force=False):
        """Check whether all calculated 3D projection views are available.

        ``force`` is retained for call-site compatibility but no longer causes
        GUI-side re-reading/recalculation.  The controller owns loading and
        view preparation, so this simply asks for the current stored views.
        """
        if not (self._is_3p_projection_case() or self.spectral_dim_count == 3):
            return False
        views = self._projection_views_3d(decon=True)
        self.no_decon = views is None
        return views is not None

    def OnSize(self,event):
        self.resized = True # set dirty
        if self.GetAutoLayout():
               self.Layout()

    def OnIdle(self,event):
        # print('Resized:', self.resized)
        if self.resized: 
            # take action if the dirty flag is set
            self.background_save(event)
            self.resized = False # reset the flag

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=line.split()
            peak.append(linetosave)
        peakfile.close()
        return peak

    def drawing_box(self):
        # Borderless toolbar controls, consistent with Slice2D.
        self.vbox2 = wx.BoxSizer(wx.HORIZONTAL)

        self.contourbutton = wx.Button(self, -1, "Contours", size=(-1,22))
        self.contourbutton.Hide()  # Replaced by the Matplotlib contour tool.
        # self.fitbutton = wx.Button(self, -1, "Fit!", size=(-1,22))
        self.cb_grid = wx.CheckBox(self, -1,"Peaks",style=wx.ALIGN_RIGHT)
        self.cb_grid.Hide()  # State-only; visible toggle is the Matplotlib Peaks tool.
        self.cb_calc = _ToolbarToggleState(False)
        if self._is_true_2d_spectrum():
            self.horizbutton = wx.ToggleButton(self, -1, "↔", size=(34, 22))
            self.vertbutton = wx.ToggleButton(self, -1, "↕", size=(34, 22))
            self.Bind(wx.EVT_TOGGLEBUTTON, self.on_horiz_button, self.horizbutton)
            self.Bind(wx.EVT_TOGGLEBUTTON, self.on_vert_button, self.vertbutton)

        self.Bind(wx.EVT_BUTTON, self.on_contour_button, self.contourbutton)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        # self.Bind(wx.EVT_BUTTON, self.peak_fit, self.fitbutton)

        # self.vbox2.Add(self.fitbutton, border=10, flag=self.flags)
        self.vbox2.AddSpacer(10)

    def contour_box(self):
        # Keep the contour controls owned by Projection so all existing 2D,
        # 3D and 4D plotting paths continue to read the same widget state.
        self.contourFrame = wx.Frame(self, title='Contours',
                                     style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        panel = wx.Panel(self.contourFrame)
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.text1 = wx.StaticText(panel, -1, 'Min:')
        self.text2 = wx.StaticText(panel, -1, 'Factor:')
        self.text3 = wx.StaticText(panel, -1, 'Number:')
        self.textbox0 = wx.TextCtrl(panel, size=(100, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox1 = wx.TextCtrl(panel, size=(50, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox2 = wx.TextCtrl(panel, size=(50, 22), style=wx.TE_PROCESS_ENTER)
        self.textbox0.SetValue(str(self.thresh))
        self.textbox1.SetValue(str(1.2))
        self.textbox2.SetValue(str(15))

        for ctrl in (self.textbox0, self.textbox1, self.textbox2):
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_draw_button)
        for widget in (self.text1, self.textbox0, self.text2, self.textbox1, self.text3, self.textbox2):
            sizer.Add(widget, 0, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)

        self.contourCloseButton = wx.Button(panel, -1, 'Close', size=(-1, 22))
        self.contourCloseButton.Bind(wx.EVT_BUTTON, self._hide_contour_frame)
        sizer.Add(self.contourCloseButton, 0, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)

        panel.SetSizerAndFit(sizer)
        frame_sizer = wx.BoxSizer(wx.VERTICAL)
        frame_sizer.Add(panel, 1, wx.EXPAND)
        self.contourFrame.SetSizerAndFit(frame_sizer)
        self.contourFrame.Bind(wx.EVT_CLOSE, self._hide_contour_frame)

    def _hide_contour_frame(self, event):
        self.contourFrame.Hide()
        if event is not None and hasattr(event, 'Veto'):
            event.Veto()

    def on_contour_button(self, event):
        if not self.contourFrame.IsShown():
            self.contourFrame.Show()
        self.contourFrame.Raise()
        self.textbox0.SetFocus()

    def shiftX_box(self):
        self.shiftxLbl = wx.StaticBox(self,-1,'ShiftX2:')
        self.shiftxSizer=wx.StaticBoxSizer(self.shiftxLbl,wx.HORIZONTAL)

        self.textboxPDB = wx.TextCtrl(self.shiftxLbl,size=(100,22),style=wx.TE_PROCESS_ENTER)
        self.textboxChains = wx.TextCtrl(self.shiftxLbl,size=(50,22),style=wx.TE_PROCESS_ENTER)
        self.textboxrH = wx.TextCtrl(self.shiftxLbl,size=(40,22),style=wx.TE_PROCESS_ENTER)
        self.textboxrC = wx.TextCtrl(self.shiftxLbl,size=(40,22),style=wx.TE_PROCESS_ENTER)

        self.shift1=wx.BoxSizer(wx.HORIZONTAL)
        self.shift2=wx.BoxSizer(wx.HORIZONTAL)
        self.shift3=wx.BoxSizer(wx.HORIZONTAL)

        self.textP=wx.StaticText(self.shiftxLbl,-1,'PDB: ')
        self.textC=wx.StaticText(self.shiftxLbl,-1,'Chain: ')
        self.textrHt=wx.StaticText(self.shiftxLbl,-1,'H:')
        self.textrCt=wx.StaticText(self.shiftxLbl,-1,'C:')

        self.textboxrH.SetValue('0')
        self.textboxrC.SetValue('0')

        self.shiftxButton = wx.Button(self.shiftxLbl, -1, "Run", size=(-1,22))
        self.openShiftxFileBtn = wx.Button(self.shiftxLbl, label="...", size=(40,22))
        self.cb_shiftx = wx.CheckBox(self.shiftxLbl, -1,"show",style=wx.ALIGN_RIGHT)

        self.Bind(wx.EVT_BUTTON, self.OnButtonShiftx, self.shiftxButton)
        self.openShiftxFileBtn.Bind(wx.EVT_BUTTON, lambda evt: self.onGetFile(evt, self.textboxPDB))
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_shiftx)
        border = 10
        self.flags2 = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL
        self.shift1.Add(self.textP, border=border, flag=self.flags)
        self.shift1.Add(self.textboxPDB, border=border, flag=self.flags2)
        self.shift1.Add(self.openShiftxFileBtn, border=border, flag=self.flags2|wx.LEFT)

        self.shift2.Add(self.textC, border=border, flag=self.flags)
        self.shift2.Add(self.textboxChains, border=border, flag=self.flags2)

        self.shift2.Add(self.textrHt, border=border, flag=self.flags)
        self.shift2.Add(self.textboxrH, border=border, flag=self.flags2)
        self.shift2.Add(self.textrCt, border=border, flag=self.flags)
        self.shift2.Add(self.textboxrC, border=border, flag=self.flags2)

        self.shift3.Add(self.shiftxButton, border=border, flag=self.flags)
        self.shift3.Add(self.cb_shiftx, border=border, flag=self.flags)

        self.shiftxSizer.Add(self.shift1)
        self.shiftxSizer.Add(self.shift2)
        self.shiftxSizer.Add(self.shift3)
        self.shiftxSizer.AddSpacer(5)
        self.shiftxSizer.Disable()


    def line_select_callback(self, eclick, erelease):
        """
        Callback for line selection.

        *eclick* and *erelease* are the press and release events.
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        # Plot coordinates are direct (X) / indirect (Y); stored data is [indirect, direct].
        x_a = int(self.uc1.f(str(eclick.xdata)+' ppm'))
        x_b = int(self.uc1.f(str(erelease.xdata)+' ppm'))
        y_a = int(self.uc0.f(str(eclick.ydata)+' ppm'))
        y_b = int(self.uc0.f(str(erelease.ydata)+' ppm'))
        print(x_a, x_b, y_a, y_b)

        y_1 = min(y_a, y_b)
        y_2 = max(y_a, y_b)
        x_1 = min(x_a, x_b)
        x_2 = max(x_a, x_b)

        count = 0

        lower_x = min(x1, x2)
        lower_y = min(y1, y2)
        higher_x = max(x1, x2)
        higher_y = max(y1, y2)

        bottom_right = (lower_x, lower_y)
        top_left = (higher_x, higher_y)
        

        for number, peak in enumerate(self.peak_list):
            x_peak, y_peak = peak[2], peak[3]
            if x_peak < top_left[0] and x_peak > bottom_right[0] and y_peak < top_left[1] and y_peak > bottom_right[1]:
                if count != 0:
                    print('This box contains more than one peak!')
                    return
                count+=1
                self.selected = self.peak_list_names[number]

        
        data_coords = (x_1,x_2,y_1,y_2)

        self.fuda_number = 0
        result = self.line_fitting.prelim_fuda_thread(data_coords, self.selected)
        self.draw_bar()
       
        self.axes_proj.remove()
        # self.axes_proj = self.fig.add_subplot(122, projection='3d')
        self.axes_proj = self.fig.add_subplot(self.gridspec[:2,1], projection='3d')



        self.line_fitting.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)
        
        # self.fuda_number += 1
   

        if self.peak_fitting_input == True:
            self.peak_fitted_input = True
            self.info_text.set_text("Thanks!  How does the fit look?")

    
    def draw_2d(self, keepaxes=False):
        levels=self.GetLevels()
        if keepaxes == True:
            print('Keeping')
            xlim = self.axes.get_xlim()
            ylim = self.axes.get_ylim()
        self.first_open = True

        self.Bind(wx.EVT_SIZE,self.OnSize)
        self.Bind(wx.EVT_IDLE,self.OnIdle)
        
        self.axes.clear()
        # A full contour redraw invalidates the cached background used by the
        # transient fitting-selection markers.
        self._fitting_selection_background = None
        self._fitting_selection_artists = []
        self.axes1D.clear()
        if hasattr(self, 'axesV'):
            self.axesV.clear()
        self.extra_lines = []

        colormap=cm.Blues
        colormap2=cm.Reds
        # colormap=cm.seismic

        # x = self.twod_data.shape[1] - 10
        # Stored 2D arrays are [indirect, direct], so direct is X and indirect is Y.
        positive = numpy.fabs(self.twod_data*(self.twod_data>0.))
        negative = numpy.fabs(self.twod_data*(self.twod_data<0.))
        self.axes.contour(self.uc1.ppm_scale(), self.uc0.ppm_scale(), positive, levels, colors='r', linewidths=0.5)
        if (negative > levels[0]).any():
            self.axes.contour(self.uc1.ppm_scale(), self.uc0.ppm_scale(), negative, levels, colors='b', linewidths=0.5)

        self.axes.set_xlabel(self.projection_service.labels[1] + " (ppm)", fontsize=8)
        self.axes.set_ylabel(self.projection_service.labels[0] + " (ppm)", fontsize=8)

        
        # self.axes.contour( self.uc0.ppm_scale(), self.uc1.ppm_scale(), numpy.fabs(self.twod_data.T*(self.twod_data.T>0.)),levels, cmap=cm.Reds, norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)),linewidths=0.5) #
        # self.axes.contour( self.uc0.ppm_scale(), self.uc1.ppm_scale(), numpy.fabs(self.twod_data.T*(self.twod_data.T>0.)),levels, colors = 'r',norm=colors.Normalize(vmin=0,vmax=numpy.max(levels)), linewidths=0.5) #
        # print(numpy.fabs(self.twod_data.T*(self.twod_data.T<0.)))
        # exit()
        # if (numpy.fabs(self.twod_data.T*(self.twod_data.T<0.))>levels[0]).any():
        #     self.axes.contour(self.uc0.ppm_scale(), self.uc1.ppm_scale(), numpy.fabs(self.twod_data.T*(self.twod_data.T<0.)),levels, colors = 'b', linewidths=0.5) #

        self.fig.canvas.mpl_connect('button_press_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_2d)

        self.dt = FileDrop(self.canvas, self.axes, levels, self)
        self.canvas.SetDropTarget(self.dt)


        if self.cb_calc.IsChecked():
            decon_view = self._spectrum_view(True)
            decon_data = None if decon_view is None else decon_view.get('ZZ')
            if decon_data is not None:
                colormap = cm.Blues
                self.axes.contour(
                    decon_view['x_axis'], decon_view['y_axis'], decon_data, levels,
                    cmap=colormap,
                    norm=colors.Normalize(vmin=-numpy.max(levels), vmax=numpy.max(levels)),
                    linewidths=0.5,
                )
                self.no_decon = False
            else:
                self.no_decon = True

        x_min,x_max=self.axes.get_xlim()
        y_min,y_max=self.axes.get_ylim()
        self.axes.set_ylim(y_max,y_min)
        self.axes.set_xlim(x_max,x_min)

        if keepaxes == True:
            self.axes.set_xlim(xlim[0], xlim[1])
            self.axes.set_ylim(ylim[0], ylim[1])
            
        # Full3D-style orthogonal marginal traces for true 2D data.
        if self._is_true_2d_spectrum():
            self.axes1D.set_ylabel('')
            self.axesV.set_xlabel('Intensity')
            x0 = float(self.uc1.ppms_scale[len(self.uc1.ppms_scale)//2])
            y0 = float(self.uc0.ppms_scale[len(self.uc0.ppms_scale)//2])
            if self.last_mouse_x is None: self.last_mouse_x = x0
            if self.last_mouse_y is None: self.last_mouse_y = y0
            self.line1, = self.axes1D.plot([], [], color='r', lw=0.6)
            self.line_decon, = self.axes1D.plot([], [], color='b', lw=0.6)
            self.line_v, = self.axesV.plot([], [], color='r', lw=0.6)
            self.line_v_decon, = self.axesV.plot([], [], color='b', lw=0.6)
            self.cross_h = self.axes.axhline(self.last_mouse_y, color='0.25', lw=0.5)
            self.cross_v = self.axes.axvline(self.last_mouse_x, color='0.25', lw=0.5)
            self._configure_2d_trace_visibility()
            self._update_2d_traces(self.last_mouse_x, self.last_mouse_y, redraw=False)
            self.fig.subplots_adjust(left=0.10, right=0.94, bottom=0.11, top=0.94)
            # Background is captured after the static contour is drawn below.
            self._trace_blit_background = None

        self.axes.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.axes.callbacks.connect('ylim_changed', self.on_ylims_change)
        #self.axes1D = self.fig.add_subplot(122)
        #self.info_text = self.axes.text(0.1,0.9,'ashdf', transform=self.axes.transAxes)
        # exit()

        # self.on_cb_grid(None)
        
        
        # self.canvas.draw()

    
    def show_fitting_peak_markers(self, peak_names):
        """Show the selected FUDA peak/group as transient X markers.

        ``peak_names`` are resolved against the authoritative reference peak
        list.  For true 2D data the reference peak coordinates already match
        the main Projection axes (X=direct, Y=indirect).  Blitting is used
        after the first update; if the backend/canvas cannot blit we fall back
        to a normal draw without affecting fitting behaviour.
        """
        if not self._is_true_2d_spectrum() or not hasattr(self, 'axes'):
            return
        names = [str(name) for name in (peak_names or [])]
        refs = {str(getattr(pk, 'name', '')): pk
                for pk in self.projection_service.reference_peaks()}
        peaks = [refs[name] for name in names if name in refs]

        try:
            # Restore the clean spectrum from the previous selection.
            if self._fitting_selection_background is not None:
                self.canvas.restore_region(self._fitting_selection_background)
            else:
                # Ensure the renderer contains the current contours before the
                # first background snapshot.
                self.canvas.draw()
                self._fitting_selection_background = self.canvas.copy_from_bbox(self.axes.bbox)

            for artist in self._fitting_selection_artists:
                try:
                    artist.remove()
                except (ValueError, RuntimeError):
                    pass
            self._fitting_selection_artists = []

            for pk in peaks:
                artist, = self.axes.plot([float(pk.x)], [float(pk.y)],
                                         marker='x', linestyle='None',
                                         markersize=9, markeredgewidth=1.5,
                                         zorder=20, animated=True)
                self._fitting_selection_artists.append(artist)
                self.axes.draw_artist(artist)
            self.canvas.blit(self.axes.bbox)
            self._fitting_selection_names = names
        except Exception:
            # Some wx/matplotlib combinations do not provide a usable blit
            # buffer.  Keep the feature functional with a conventional redraw.
            self._fitting_selection_background = None
            for artist in self._fitting_selection_artists:
                try:
                    artist.set_animated(False)
                except Exception:
                    pass
            self.canvas.draw_idle()
            self._fitting_selection_names = names


    def peak_fit(self, event):

        # self.uc0 = self.tabOne.uc1
        # self.uc1 = self.tabOne.uc2

        self.peak_list=[]
        self.peak_list_names=[]
        self.conn_data=self.projection_service.connections

        if(len(self.conn_data)>0):
            for cn in self.conn_data:
                loc1=float(cn.f1)
                loc2=float(cn.f2)
                index1 = self.uc0.f(str(loc1)+' ppm')
                index2 = self.uc1.f(str(loc2)+' ppm')
                lab=cn.p1
                self.peak_list.append([index1, index2, loc1, loc2, int(re.findall(r'[0-9]+', lab)[0])])
                self.peak_list_names.append(lab)



        if len(self.peak_list)==0:
                self.on_cb_grid(event)
        if len(self.peak_list)==0:
                print("No peaks read in!")
                return

        self.peak_list = numpy.array(self.peak_list)
        

        if self.peak_fitted_input == False:
            self.line_fitting = Unidec_line_fitting(self.projection_service.data, self.peak_list,self.peak_list_names, self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2, float(self.projection_service.threshold_fraction()),self.uc0, self.uc1)
            loading_result = self.line_fitting.load_results()
            if loading_result != False:
                print('Results already exist!')
            self.info_text.set_text('Please drag around an isolated peak!')
            self.fig.clf()
            self.gridspec = GridSpec(3,2, figure = self.fig)
            self.axes = self.fig.add_subplot(self.gridspec[:,0])
            self.axes_proj = self.fig.add_subplot(self.gridspec[:2,1], projection='3d')
            self.axes_bar = self.fig.add_subplot(self.gridspec[2:,1])
            self.axes_bar_amp = self.axes_bar.twinx()
            self.draw_bar()
            # self.line_fitting.bar_chart(self.axes_bar)
            self.draw_2d()
   
            self.selector = RectangleSelector(self.axes, self.line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # disable middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels',
                                           interactive=False)

            
            
            self.background_save(event)
            self.canvas.draw()
            self.peak_fitting_input = True
            return

        thread = threading.Thread(target=self.fuda_thread)
        thread.setDaemon(True)
        thread.start()

        self.canvas.draw()

    # def draw_bar(self):


        
       

        

    def fuda_thread(self):
        self.line_fitting.finding_overlaps()
        for peak in self.peak_list_names:
            self.line_fitting.fit_unoverlapped_peak(peak)
            wx.CallAfter(self.print_intensities)
            wx.CallAfter(self.draw_bar)


        self.overlap_thread()

    def draw_bar(self):
        self.axes_bar.cla()
        self.axes_bar_amp.cla()
        for key, item in self.line_fitting.intensities.items():
            if key == self.selected:
                intensity = item[0]
                fwhm_x = item[1]
                fwhm_y = item[2]
                amp = item[3]
                self.axes_bar.bar(int(re.findall(r'[0-9]+', key)[0]), intensity, color='r', edgecolor='k', width=0.5)
                self.axes_bar_amp.bar(int(re.findall(r'[0-9]+', key)[0])+0.5, amp, color='darkblue', edgecolor='k', width=0.5)
            else:
                intensity = item[0]
                fwhm_x = item[1]
                fwhm_y = item[2]
                amp = item[3]
                self.axes_bar.bar(int(re.findall(r'[0-9]+', key)[0]), intensity, color='gray', edgecolor='k', width=0.5)
                self.axes_bar_amp.bar(int(re.findall(r'[0-9]+', key)[0])+0.5, amp, color='lightblue', edgecolor='k', width=0.5)

        self.canvas.draw()

    def overlap_thread(self):
        for x in range(len(self.line_fitting.final_overlaps)):
            self.line_fitting.fit_overlapped_peaks(x)
            wx.CallAfter(self.print_intensities)
            wx.CallAfter(self.draw_bar)





    def print_intensities(self):
        outy = open('./out/fuda_outputs.txt', 'w')
        outy.write('Peak\tIntegral\tfwhm x (ppm)\tfwhm y (ppm)\tamplitude\n')
        for atom, intensity in self.line_fitting.intensities.items():
            outy.write('%s\t%f\t%f\t%f\t%f\n' % (atom, intensity[0], intensity[1], intensity[2], intensity[3]))
        outy.close()


    def determine_overlap_area(self, peaks):
        lowest_y = 1e10
        highest_y = -1
        lowest_x = 1e10
        highest_x = -1
        mask = numpy.ones_like(self.projection_service.data)
        peaks2=[]
        # i = 0
        for peak in peaks:

            x = peak[1]
            y = peak[0]
            print(x, y, self.fwhm_x, self.fwhm_y)
            x_1 = int(x-self.fwhm_x*self.overlap_distance)
            x_2 = int(x+self.fwhm_x*self.overlap_distance)
            y_1 = int(y-self.fwhm_y*self.overlap_distance)
            y_2 = int(y+self.fwhm_y*self.overlap_distance)
            print(y_1, y_2, x_1,x_2)
            if y_1 < lowest_y:
                lowest_y = y_1
            if x_1 < lowest_x:
                lowest_x = x_1
            if y_2 > highest_y:
                highest_y = y_2
            if x_2 > highest_x:
                highest_x = x_2
            mask[y_1:y_2, x_1:x_2] = 0

        for peak in peaks:
            x = peak[1]
            y = peak[0]
            peaks2.append([y-lowest_y, x-lowest_x])
            # break
        print('peaks2', peaks2, lowest_y, highest_y, lowest_x, highest_x, peak[1], peak[0])
        masked = numpy.ma.masked_array(self.projection_service.data, mask=mask, fill_value=0.0)[lowest_y:highest_y, lowest_x:highest_x]
        # masked = numpy.ma.set_fill_value(masked, 0)
        # print(masked.filled())
        # exit()
        return peaks2, masked.filled()

    def on_scroll(self, event):
        # print('scrolling')
        self.ymin,self.ymax=self.axes1D.get_ylim()
        self.axes1D.set_ylim(self.ymin+(self.ymin*0.05*event.step), self.ymax+(self.ymax*0.05*event.step))
        self.canvas.draw()

    # def OnSize(self, event):
    #     print('sized!')
        # self.background = self.canvas.copy_from_bbox(self.axes.bbox)



    def create_main_panel(self):
        self.no_decon=False
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fig = Figure(constrained_layout=False)
        self.fig.clear()
        pass
        if self._is_pseudo2d_projection_case():
            self.axes1 = self.fig.add_subplot(111)
            self.axes2 = self.axes3 = None
        elif self._is_true_2d_spectrum():
            self.axes = self.fig.add_subplot(111)
            self.axes1D = self.axes.twinx()
            self.axesV = self.axes.twiny()
            self.axes1D.patch.set_visible(False)
            self.axesV.patch.set_visible(False)
            self.axes1D.set_navigate(False)
            self.axesV.set_navigate(False)
            self.fig.subplots_adjust(left=0.10, right=0.94, bottom=0.11, top=0.94)
            self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)

            view = self._spectrum_view(False)
            if view is None or view.get('ZZ') is None:
                raise RuntimeError("Projection requires raw 2D spectrum view in data_store")
            self.peak_fitted_input=False
            self.peak_fitting_input=False
            self.fitted_unoverlapped=0
            self.fitted_overlapped=0
            self.unoverlapped=[]
            shape = self.projection_service.peak_shape_parameters(2)
            self.Gamma_x, self.Gamma_y = shape['lorentz']
            self.sigma_x, self.sigma_y = shape['sigma']
            self.nu1 = shape['voigt'][1]
            self.nu2 = shape['voigt'][1]

        elif self._is_3p_projection_case():
                # 3p has only one meaningful spectral projection: the sum over
                # the real axis.  Give that projection the full canvas.
            self.axes1 = self.fig.add_subplot(111)
            self.axes2 = None
            self.axes3 = None
            if self._projection_views_3d(decon=False) is None:
                raise RuntimeError("Projection requires cached raw 3p projection view in data_store")
            self.no_decon = self._projection_views_3d(decon=True) is None

        elif self.spectral_dim_count == 3:
            self.axes1 = self.fig.add_subplot(131)
            self.axes2 = self.fig.add_subplot(132)
            self.axes3 = self.fig.add_subplot(133)
            if self._projection_views_3d(decon=False) is None:
                raise RuntimeError("Projection requires cached raw projection views in data_store")
            self.no_decon = self._projection_views_3d(decon=True) is None

        elif(self.spectral_dim_count==4):
            from matplotlib.gridspec import GridSpec
            gs1=GridSpec(2,3)
            self.axes1 = self.fig.add_subplot(gs1[:,0])
            self.axes2 = self.fig.add_subplot(gs1[:,1])
            self.axes3 = self.fig.add_subplot(gs1[0,2])
            self.axes4 = self.fig.add_subplot(gs1[1,2])
            if self._projection_views_4d() is None:
                raise RuntimeError("Projection requires cached raw projection views in data_store")

        ## Initialise main matplotlib canvas
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(1,1))

        # Pseudo2D uses the Projection panel as its 1D reference editor.
        # Historically on_pick was only connected from the true-2D drawing
        # path below, so the pseudo2D Tools palette changed modes correctly
        # but clicks on the projection canvas never reached the editor.
        # Connect the pseudo2D canvas explicitly here, once, when it is built.
        self._pseudo2d_tool_click_cid = None
        if self._is_pseudo2d_projection_case():
            self._pseudo2d_tool_click_cid = self.canvas.mpl_connect(
                'button_press_event', self.on_pick)

        ## Initialise navigation toolbar
        self.toolbar = RedrawNavigationToolbar(self.canvas, self.redraw_view, peak_callback=self._toolbar_peaks, decon_callback=self._toolbar_decon, contour_callback=(None if self._is_pseudo2d_projection_case() else self._toolbar_contours), tools_callback=(self._toolbar_tools if self._is_pseudo2d_projection_case() else None), coordinates=False)
        self.toolbar.Realize()
        if self._is_pseudo2d_projection_case():
            self.full_peak_tools_box()

        ## Adding our control boxes
        self.flags = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.LEFT | wx.ALIGN_CENTER_VERTICAL
        self.drawing_box()
        self.contour_box()
        # self.shiftX_box()

        ## Put our controls inside Matplotlib's native wx.ToolBar.  Inserting
        ## them before Matplotlib's stretch/coordinate items makes the native
        ## wx.TB_BOTTOM border continue uninterrupted above every control.
        # Keep Matplotlib's native tools together, followed by our controls and
        # finally the live coordinate readout.
        if self._is_true_2d_spectrum():
            self.toolbar.AddSeparator()
            for widget in (self.horizbutton, self.vertbutton):
                widget.Reparent(self.toolbar)
                self.toolbar.AddControl(widget)
            self.toolbar.bind_control_status_help(self.horizbutton, 'Toggle horizontal projection trace')
            self.toolbar.bind_control_status_help(self.vertbutton, 'Toggle vertical projection trace')
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND) # Main matplotlib canvas
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 4)

        ## Make vbox the main sizer
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)


    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    """
    def GetData(self,infile):
        # print infile
        input=self.readfile(infile)
        # print input
        xs=[]
        ys=[]
        zs=[]
        Xs=[]
        Ys=[]
        Zs=[]
        for i in range(len(input)):
            if(len(input[i])!=0):
                xs.append(float(input[i][0]))
                ys.append(float(input[i][1]))
                zs.append(float(input[i][2]))
            else:
                Xs.append(xs)
                Ys.append(ys)
                Zs.append(zs)
                zs=[]
                ys=[]
                xs=[]
        if(len(xs)!=0):
            Xs.append(xs)
            Ys.append(ys)
            Zs.append(zs)
        return numpy.array(Xs),numpy.array(Ys),numpy.array(Zs)
    """
    def load_decon_data(self):
        """Update decon availability without copying projection data into GUI state."""
        if self._is_3p_projection_case() or self.spectral_dim_count == 3:
            self.no_decon = self._projection_views_3d(decon=True) is None
            return not self.no_decon
        if self.spectral_dim_count == 4:
            planes = [self._cached_decon_plane(name) for name in ('yz', 'xz', 'xy', 'za')]
            self.no_decon = not any(payload is not None for payload in planes)
            return not self.no_decon
        return not self.no_decon

    def background_save(self, event):
        # In true-2D mode the marginal traces are dynamic artists.  Do not
        # recreate a row trace here: that old behaviour made a horizontal
        # spectrum appear even when the slice toggle was off and also replaced
        # the red/blue dynamic artists with legacy black/C0 lines.
        if self._is_true_2d_spectrum() and hasattr(self, 'line1'):
            self._save_2d_trace_background()
            return

        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)

    def on_xlims_change(self, event_ax):
        if self._is_true_2d_spectrum() and self.first_open == False:

            self.background_save(event_ax)

    def on_ylims_change(self, event_ax):
        if self._is_true_2d_spectrum() and self.first_open == False:
            self.background_save(event_ax)


    def draw_figure(self,scale='y', keepaxes=False):
        """ Redraws the figure
        """
        levels=self.GetLevels()

        if self._is_pseudo2d_projection_case():
            # Toolbar overlays must not reset a user's ppm/intensity zoom.
            saved_xlim = self.axes1.get_xlim() if keepaxes and self.axes1.has_data() else None
            saved_ylim = self.axes1.get_ylim() if keepaxes and self.axes1.has_data() else None
            self.axes1.clear()
            x, projection = self._pseudo2d_projection()
            if x is None or projection is None:
                return
            self.axes1.plot(x, projection, color='r', lw=0.7)
            if self.cb_calc.IsChecked():
                xd, calculated = self._pseudo2d_projection(decon=True)
                self.no_decon = calculated is None
                if calculated is not None:
                    self.axes1.plot(xd, calculated, color='b', lw=0.7)
                else:
                    self.cb_calc.SetValue(False)
                    if hasattr(self, 'toolbar'):
                        self.toolbar.set_decon_active(False)
            else:
                _, calculated = self._pseudo2d_projection(decon=True)
                self.no_decon = calculated is None
            if self.cb_grid.IsChecked():
                ymax = float(numpy.max(projection)) if projection.size else 1.0
                ymin = float(numpy.min(projection)) if projection.size else 0.0
                yrange = ymax - ymin or max(abs(ymax), 1.0)
                for ppm, intensity, label in self._pseudo2d_peak_overlay():
                    # Match the established Slice1D convention: a peak is a stick
                    # from zero to its fitted/deconvolved intensity, not a full-height line.
                    selected = str(label) == str(self.full_selected_name)
                    self.axes1.vlines(ppm, 0.0, intensity, color=('green' if selected else 'k'),
                                      lw=(2.0 if selected else 0.8), alpha=0.9)
                    if label:
                        self.axes1.text(ppm, intensity, label, rotation=90, fontsize=7,
                                        ha='center', va='bottom' if intensity >= 0 else 'top',
                                        color=('green' if selected else 'k'))
            # Selection is independent of peak-overlay visibility: a Full-list
            # Show or Tools Select action always gets a full-height green cursor.
            selected_record = self._full_record()
            if selected_record is not None:
                try:
                    selected_ppm = float(selected_record.get('coordinates', ())[0])
                    self.axes1.axvline(selected_ppm, color='green', lw=1.5, alpha=0.95, zorder=20)
                except (IndexError, TypeError, ValueError):
                    pass
            label = self.projection_service.pseudo2d_data(ensure_file=False)['label']
            self.axes1.set_xlabel((label or 'Direct') + ' (ppm)', fontsize=8)
            self.axes1.set_ylabel('Projected intensity', fontsize=8)
            # NMR convention: chemical shift always decreases left-to-right.
            if saved_xlim is not None and saved_ylim is not None:
                lo, hi = sorted(saved_xlim)
                self.axes1.set_xlim(hi, lo)
                self.axes1.set_ylim(saved_ylim)
            elif len(x):
                self.axes1.set_xlim(float(numpy.nanmax(x)), float(numpy.nanmin(x)))

        elif self._is_true_2d_spectrum():

            self.draw_2d(keepaxes=keepaxes)
            self.Bind(wx.EVT_SIZE,self.OnSize)
        
            self.Bind(wx.EVT_IDLE,self.OnIdle)
            # self.peak_fitted_input = False


        elif self._is_3p_projection_case():
            self.axes1.clear()
            raw_views = self._projection_views_3d(decon=False)
            if raw_views is None:
                raise RuntimeError('Projection requires cached raw 3p projection view in data_store')
            XX1, YY1, ZZ1, axlab1 = self._unpack_projection_view(raw_views[0])
            colormap = cm.Reds
            colormap_neg = cm.Blues
            colormap_decon = cm.Greens
            self.axes1.contour(XX1, YY1, ZZ1, levels, cmap=colormap,
                               norm=colors.Normalize(vmin=-numpy.max(levels), vmax=numpy.max(levels)),
                               linewidths=0.5)
            self.axes1.contour(XX1, YY1, -ZZ1, levels, cmap=colormap_neg,
                               norm=colors.Normalize(vmin=-numpy.max(levels), vmax=numpy.max(levels)),
                               linewidths=0.5)
            decon_view = self._projection_view_3p(decon=True) if self.cb_calc.IsChecked() else None
            self.no_decon = decon_view is None
            if decon_view is not None:
                XXd, YYd, ZZd, _ = self._unpack_projection_view(decon_view)
                self.axes1.contour(XXd, YYd, ZZd, levels, cmap=colormap_decon,
                                   norm=colors.Normalize(vmin=-numpy.max(levels), vmax=numpy.max(levels)))
            x_min, x_max = self.axes1.get_xlim()
            y_min, y_max = self.axes1.get_ylim()
            self.axes1.set_ylim(y_max, y_min)
            self.axes1.set_xlim(x_max, x_min)
            if self.cb_grid.GetValue() == 1:
                labels = self._spectral_labels_3p()
                if labels is not None:
                    peaks = self._projected_peak_overlay(labels[0], labels[1], transpose='n')
                    self._draw_projected_peak_overlay(self.axes1, peaks, labels=True, size=50)
            self.axes1.set_xlabel(axlab1[0] + " (ppm)", fontsize=8)
            self.axes1.set_ylabel(axlab1[1] + " (ppm)", fontsize=8)

        elif(self.spectral_dim_count==3):
            self.axes1.clear()
            self.axes2.clear()
            self.axes3.clear()

            raw_views = self._projection_views_3d(decon=False)
            if raw_views is None:
                raise RuntimeError('Projection requires cached raw projection views in data_store')
            (XX1, YY1, ZZ1, axlab1), (XX2, YY2, ZZ2, axlab2), (XX3, YY3, ZZ3, axlab3) = [
                self._unpack_projection_view(view) for view in raw_views
            ]

            decon_views = self._projection_views_3d(decon=True) if self.cb_calc.IsChecked() else None
            self.no_decon = decon_views is None
            if decon_views is not None:
                (XX1_decon, YY1_decon, ZZ1_decon, _), (XX2_decon, YY2_decon, ZZ2_decon, _), (XX3_decon, YY3_decon, ZZ3_decon, _) = [
                    self._unpack_projection_view(view) for view in decon_views
                ]
            #plot1
            cnt=0 #for each combination of label, get nmrPipe projection

            colormap2=cm.Reds
            # colormap=cm.seismic
            colormap=cm.Reds
            colormap_decon=cm.Greens
            colormap_neg=cm.Blues



            self.axes1.contour( XX1, YY1, ZZ1,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #

            self.axes1.contour( XX1, YY1, -ZZ1,levels,cmap=colormap_neg,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #

            

            if(self.cb_calc.IsChecked() and not self.no_decon):
                self.axes1.contour( XX1_decon,YY1_decon,ZZ1_decon,levels,cmap=colormap_decon,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)))

            #flip axes (NMR)
            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            self.axes1.set_ylim(y_max,y_min)
            self.axes1.set_xlim(x_max,x_min)
            # if(self.cb_grid.GetValue()==1): #if plotting peaks
            #     for ii in range(len(self.projection_service.peaks)):
            #         loc1=self.projection_service.peaks[ii].ppmK
            #         loc2=self.projection_service.peaks[ii].ppmJ
            #         lab=self.projection_service.peaks[ii].name
            #         self.axes1.text(loc1,loc2,lab,fontsize=8)
            #         self.axes1.scatter(loc1,loc2,c='k',s=50,zorder=2,marker='x')
            if(self.cb_grid.GetValue()==1): # if plotting projected peaks
                peaks = self._projected_peak_overlay(
                    self.projection_service.labels[2], self.projection_service.labels[1], transpose='n'
                )
                self._draw_projected_peak_overlay(self.axes1, peaks, labels=True, size=50)
            # if(self.cb_grid.GetValue()==1): #if plotting peaks
            # if(self.cb_shiftx.IsChecked()):
                # for pk in self.shiftXlist:
                #     loc1=pk.ppmK
                #     loc2=pk.ppmJ
                #     self.axes1.text(loc1,loc2,pk.name,fontsize=8,color='r')
                #     self.axes1.scatter(loc1,loc2,c='r',s=50,zorder=2,marker='x')

            self.axes1.set_xlabel(axlab1[0]+" (ppm)", fontsize=8)
            self.axes1.set_ylabel(axlab1[1]+" (ppm)", fontsize=8)

            # colormap=cm.seismic
            self.axes2.contour( XX2, YY2, ZZ2,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #
            self.axes2.contour( XX2, YY2, -ZZ2,levels,cmap=colormap_neg,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #
            if(self.cb_calc.IsChecked() and not self.no_decon):
                self.axes2.contour( XX2_decon, YY2_decon, ZZ2_decon,levels,cmap=colormap_decon,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #

            #flip axes (NMR)
            x_min,x_max=self.axes2.get_xlim()
            y_min,y_max=self.axes2.get_ylim()
            self.axes2.set_ylim(y_max,y_min)
            self.axes2.set_xlim(x_max,x_min)
            if(self.cb_grid.GetValue()==1): # if plotting projected peaks
                peaks = self._projected_peak_overlay(
                    self.projection_service.labels[2], self.projection_service.labels[0], transpose='n'
                )
                self._draw_projected_peak_overlay(self.axes2, peaks, labels=False, size=20)
            self.axes2.set_xlabel(axlab2[0]+" (ppm)", fontsize=8)
            self.axes2.set_ylabel(axlab2[1]+" (ppm)", fontsize=8)

            # colormap=cm.seismic
            self.axes3.contour( XX3, YY3, ZZ3,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #
            self.axes3.contour( XX3, YY3, -ZZ3,levels,cmap=colormap_neg,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)), linewidths=0.5) #
            if(self.cb_calc.IsChecked() and not self.no_decon):
                self.axes3.contour( XX3_decon, YY3_decon, ZZ3_decon,levels,cmap=colormap_decon,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #

            #flip axes (NMR)
            x_min,x_max=self.axes3.get_xlim()
            y_min,y_max=self.axes3.get_ylim()
            self.axes3.set_ylim(y_max,y_min)
            self.axes3.set_xlim(x_max,x_min)
            if(self.cb_grid.GetValue()==1): # if plotting projected peaks
                peaks = self._projected_peak_overlay(
                    self.projection_service.labels[1], self.projection_service.labels[0], transpose='y'
                )
                self._draw_projected_peak_overlay(self.axes3, peaks, labels=False, size=20)
            self.axes3.set_xlabel(axlab3[0]+" (ppm)", fontsize=8)
            self.axes3.set_ylabel(axlab3[1]+" (ppm)", fontsize=8)

            # self.multi = MultiCursor(self.canvas, (self.axes1, self.axes2, self.axes3), horizOn = True, color='k', lw=0.5)




        elif(self.spectral_dim_count==4):

            self.axes1.clear()
            self.axes2.clear()
            self.axes3.clear()
            self.axes4.clear()
            raw_views = self._projection_views_4d()
            if raw_views is None:
                raise RuntimeError('Projection requires cached raw projection views in data_store')
            (XX1, YY1, ZZ1, axlab1), (XX2, YY2, ZZ2, axlab2), (XX3, YY3, ZZ3, axlab3), (XX4, YY4, ZZ4, axlab4) = [
                self._unpack_projection_view(view) for view in raw_views
            ]
            colormap=cm.Blues
            # colormap=cm.seismic
            colormap=cm.Reds
            #XX,YY,ZZ=self.GetData('out/'+'xy'+'.out')

            self.axes1.contour(XX1,YY1,ZZ1,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)))
            if(self.cb_calc.GetValue()==1): #if plotting peaks
                cached = self._cached_decon_plane("za")
                if cached is not None:
                    colormap=cm.Blues
                    Xs,Ys,Zs = cached["Xs"], cached["Ys"], cached["Zs"]
                    Zs=Zs/numpy.max(Zs)*numpy.max(ZZ1)
                    self.axes1.contour( Xs, Ys, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #
            if(self.cb_grid.GetValue()==1): #if plotting peaks
                for pk in self.projection_service.peaks:
                    self.axes1.text(pk.ppmI,pk.ppmJ,pk.name,fontsize=8)
                    self.axes1.scatter(pk.ppmI,pk.ppmJ,c='k',marker='x',s=50,zorder=2)
            #flip axes (NMR)
            x_min,x_max=self.axes1.get_xlim()
            y_min,y_max=self.axes1.get_ylim()
            self.axes1.set_ylim(y_max,y_min)
            self.axes1.set_xlim(x_max,x_min)
            self.axes1.set_xlabel(axlab1[0], fontsize=8)
            self.axes1.set_ylabel(axlab1[1], fontsize=8)

            if(self.cb_shiftx.IsChecked()):
                for pk in self.shiftXlist:
                    loc1=pk.ppmI
                    loc2=pk.ppmJ
                    self.axes1.text(loc1,loc2,pk.name,fontsize=8,color='r')
                    self.axes1.scatter(loc1,loc2,c='r',s=50,zorder=2,marker='x')


            colormap=cm.Blues
            colormap=cm.seismic
            colormap=cm.Reds

            #XX,YY,ZZ=self.GetData('out/'+'xy'+'.out')
            self.axes2.contour(XX2,YY2,ZZ2,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)))
            if(self.cb_calc.GetValue()==1): #if plotting peaks
                cached = self._cached_decon_plane("xy")
                if cached is not None:
                    colormap=cm.Blues
                    Xs,Ys,Zs = cached["Xs"], cached["Ys"], cached["Zs"]
                    Zs=Zs/numpy.max(Zs)*numpy.max(ZZ2)
                    self.axes2.contour( Xs, Ys, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #
            if(self.cb_grid.GetValue()==1): #if plotting peaks
                for pk in self.projection_service.peaks:
                    self.axes2.text(pk.ppmK,pk.ppmL,pk.name,fontsize=8)
                    self.axes2.scatter(pk.ppmK,pk.ppmL,c='k',marker='x',s=50,zorder=2)

            #flip axes (NMR)
            x_min,x_max=self.axes2.get_xlim()
            y_min,y_max=self.axes2.get_ylim()
            self.axes2.set_ylim(y_max,y_min)
            self.axes2.set_xlim(x_max,x_min)
            self.axes2.set_xlabel(axlab2[0], fontsize=8)
            self.axes2.set_ylabel(axlab2[1], fontsize=8)





            #self.axes = self.fig.add_subplot(130+cnt);self.axes.clear()
            # self.axes = self.fig.add_subplot(gs1[0,2])
            colormap=cm.Blues
            colormap=cm.seismic
            colormap=cm.Reds

            #XX,YY,ZZ=self.GetData('out/'+'yz'+'.out')
            self.axes3.contour(XX3,YY3,ZZ3,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)))
            if(self.cb_calc.GetValue()==1): #if plotting peaks
                cached = self._cached_decon_plane("yz")
                if cached is not None:
                    colormap=cm.Reds
                    colormap=cm.Blues
                    Xs,Ys,Zs = cached["Xs"], cached["Ys"], cached["Zs"]
                    Zs=Zs/numpy.max(Zs)*numpy.max(ZZ3)
                    self.axes3.contour( Xs, Ys, Zs,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #
            #flip axes (NMR)
            x_max,x_min=self.axes3.get_xlim()
            y_max,y_min=self.axes3.get_ylim()
            self.axes3.set_ylim(y_max,y_min)
            self.axes3.set_xlim(x_max,x_min)
            if(self.cb_grid.GetValue()==1): #if plotting peaks
                self.conn_data = self._cached_peak_overlay()
                for cn in self.conn_data:
                    # 3rd panel is the f2/f3 projection of correlate.3
                    self.axes3.scatter(float(cn.f2), float(cn.f3), c='k', marker='x', s=10, zorder=2)


            x_min,x_max=self.axes3.get_xlim()
            y_min,y_max=self.axes3.get_ylim()
            self.axes3.set_ylim(y_max,y_min)
            self.axes3.set_xlim(x_max,x_min)
            self.axes3.set_xlabel(axlab3[0], fontsize=8)
            self.axes3.set_ylabel(axlab3[1], fontsize=8)



            #cnt+=1
            #self.axes = self.fig.add_subplot(130+cnt);self.axes.clear()
            #self.axes = self.fig.add_subplot(gs1[1,2])
            colormap=cm.Blues
            colormap=cm.seismic
            colormap=cm.Reds

            #XX,YY,ZZ=self.GetData('out/'+'xz'+'.out')
            self.axes4.contour(XX4,YY4,ZZ4,levels,cmap=colormap,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels)))
            if(self.cb_calc.GetValue()==1): #if plotting peaks
                cached = self._cached_decon_plane("xz")
                if cached is not None:
                    colormap=cm.Reds
                    colormap=cm.Blues
                    Xs,Ys,Zs = cached["Xs"], cached["Ys"], cached["Zs"]
                    Zs=Zs/numpy.max(Zs)*numpy.max(ZZ4)
                    self.axes4.contour( Xs, Ys, Zs,levels,norm=colors.Normalize(vmin=-numpy.max(levels),vmax=numpy.max(levels))) #
            #flip axes (NMR)
            x_max,x_min=self.axes4.get_xlim()
            y_max,y_min=self.axes4.get_ylim()
            self.axes4.set_ylim(y_max,y_min)
            self.axes4.set_xlim(x_max,x_min)
            if(self.cb_grid.GetValue()==1): #if plotting peaks
                for cn in self.projection_service.connections:
                   # print cn.f1,cn.f2,cn.f3,cn.f4
                    self.axes4.scatter(cn.f2,cn.f4,c='k',marker='x',s=10,zorder=2)
            x_min,x_max=self.axes4.get_xlim()
            y_min,y_max=self.axes4.get_ylim()
            self.axes4.set_ylim(y_max,y_min)
            self.axes4.set_xlim(x_max,x_min)
            self.axes4.set_xlabel(axlab4[0], fontsize=8)
            self.axes4.set_ylabel(axlab4[1], fontsize=8)
        self.fig.tight_layout()

        self.canvas.draw()

    def on_key_2d(self, event):
        # print('blah')
        if event.key=='n':
            if self.selected != '':
                self.ftol = 1e-10
                # if self.fuda_number < len(self.projection_service.data) and self.fuda_number>-1:
                    # try:
                    #     for key in self.line_fitting.plotting_resim_data.keys():
                    #         if str(self.selected) == str(int(key)):
                    #             # print('fuda_fit', self.selected)
                    #             self.line_fitting.plot_fuda_fit(key, self.fuda_number, self.axes_proj, self.canvas)
                    #             self.fuda_number=0
                    #             self.selected = list(self.line_fitting.plotting_resim_data.keys())[numpy.argwhere(list(self.line_fitting.plotting_resim_data.keys()) == self.selected)+1]

                    # except:
                for la, key in enumerate(self.line_fitting.plotting_resim_data.keys()):
                    if str(self.selected) == str(key):
                        print('fuda_fit', self.selected)
                        self.selected = list(self.line_fitting.plotting_resim_data.keys())[la+1]
                        self.line_fitting.plot_fuda_fit(self.selected, self.fuda_number, self.axes_proj, self.canvas)
                        self.fuda_number=0
                        update_colors = []
                        for x in range(len(self.peaks_text)):
                            if self.peaks_text[x].get_text() == self.selected:
                                self.peaks_text[x].set_color('r')
                                self.axes.draw_artist(self.peaks_text[x])
                                update_colors.append('r')

                            elif self.peaks_text[x].get_text() in self.line_fitting.unoverlapped_names:
                                self.peaks_text[x].set_color('k')
                                self.axes.draw_artist(self.peaks_text[x])

                                update_colors.append('k')
                            else:
                                self.peaks_text[x].set_color('blue')
                                self.axes.draw_artist(self.peaks_text[x])

                                update_colors.append('blue')
                            
                        self.peaks_scatter[0].set_color(update_colors)
                        
                        blit_artists(self.canvas, self.axes, self.background, self.peaks_text + self.peaks_scatter)
                        self.background_save(None)
                        self.draw_bar()
                        # self.redraw_scatters()
                        # self.draw_2d()
                        break
        if event.key=='h':
            print('h clicked')
            if self.line1.get_visible():
                self.line1.set_visible(False)
                if(self.cb_calc.IsChecked()):
                    self.line_decon.set_visible(False)
            else:
                self.line1.set_visible(True)
                if(self.cb_calc.IsChecked()):
                    self.line_decon.set_visible(True)

        if event.key =='r' and self.selected != '':
           
            self.ftol = self.ftol*0.1
            overlapped = self.line_fitting.is_peak_overlapped(self.selected)
            if overlapped == None:
                return
            elif overlapped == False:  
                self.line_fitting.fit_unoverlapped_peak(self.selected, self.ftol)
                self.line_fitting.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)

            else:
                self.line_fitting.fit_overlapped_peaks(overlapped, self.ftol)
                self.line_fitting.plot_fuda_fit(self.selected, 0, self.axes_proj, self.canvas)
            self.print_intensities()
            self.print_lorentzian()
            self.draw_bar()


        if event.key=='enter' and self.peak_fitted_input == True:

            self.peak_fitting_input = False
            self.info_text.set_text('Thanks: now going to iterate through')
            self.peak_fit(event)
            self.canvas.draw()
            self.peak_fitted_input = False



    def print_lorentzian(self):
        outy = open('out/lorentzian_shape.out', 'a')
        nu2 = self.line_fitting.plotting_resim_data[self.selected][0][-1]
        nu1 = self.line_fitting.plotting_resim_data[self.selected][0][-2]
        sigma_y = self.line_fitting.plotting_resim_data[self.selected][0][-3]
        sigma_x = self.line_fitting.plotting_resim_data[self.selected][0][-4]
        Gamma_y = self.line_fitting.plotting_resim_data[self.selected][0][-5]
        Gamma_x = self.line_fitting.plotting_resim_data[self.selected][0][-6]
        amp = self.line_fitting.intensities[self.selected][3][0]
        outy.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (amp, Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2))

    def GetLevels(self):
        min_level=float(self.textbox0.GetValue())
        max_level=float(self.textbox1.GetValue())
        ctr_level=int(self.textbox2.GetValue())

        if(ctr_level==0):
            ctr_level=10
        if(max_level==0):
            max_level=1.2
        if(min_level==0):
            min_level=1E3

        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*max_level)
        levels=numpy.array(levels)
        # levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels

    # Projection file IO, axis construction, and transpose-view creation now
    # live in deconFrame/DataStore.  This GUI only retrieves plotting-ready
    # views through ``get_projection_view``.

    def OnButtonShiftx(self,event):
        if(os.path.exists(self.textboxPDB.GetValue())==0):
            print('pdb file not found. aborting.')
            return
        print('Running shiftx2')
        if(os.path.exists(self.textboxPDB.GetValue()+'.cs')==0):
            print('shiftx2.py -i'+self.textboxPDB.GetValue()+' -c '+self.textboxChains.GetValue())
            os.system('shiftx2.py -i'+self.textboxPDB.GetValue()+' -c '+self.textboxChains.GetValue())

        print('Parsing...')
        from spinDecon.analysis.shiftx_post_filter import shiftXNMR

        seleMet=self.parent.tabMagma.GetCombo(self.parent.tabMagma.methylBox,self.parent.tabMagma.methyls)

        shiftX = shiftXNMR(self.textboxPDB.GetValue()+'.cs','',self.textboxChains.GetValue(),seleMet,'')
        shiftX.SetMethylDict()
        shiftX.ReadShiftxFile()
        self.Make2Dlist(shiftX.shiftXDict) #make a 2D peak list
        self.draw_figure()

    def Make2Dlist(self,shiftXDict):

        refH=float(self.textboxrH.GetValue())
        refC=float(self.textboxrC.GetValue())

        self.shiftXlist=[]
        for key,vals in list(shiftXDict.items()):
            lab=str(key)+vals['type']
            for koi,vols in list(vals['name'].items()):
                if(koi[0]=='H'):
                    if('C'+koi[1:] in list(vals['name'].keys())):
                        test=lab+koi[1:],vals['name']['C'+koi[1:]]+refC,vols+refH #entry goes name,carbon,proton
                        self.shiftXlist.append(peakEntry(test))
        #now need to alias
        self.DoAlias()

    def DoAlias(self):
        if(self.parent.spectral_dim_count==2):
            for p in range(len(self.shiftXlist)):
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].y,0)
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].x,1)
        if(self.parent.spectral_dim_count==3):
            for p in range(len(self.shiftXlist)):
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].y,0)#C
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].y,1)#C
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].x,2)#H
        if(self.parent.spectral_dim_count==4):
            for p in range(len(self.shiftXlist)):
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].x,0) #i  #H
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].y,1) #j  #C
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].x,2) #k  #H
                self.projection_service.alias(self.shiftXlist[p],self.shiftXlist[p].y,3) #l  #C

    # def Multiple_cursor(self, )

    # def on_cb_grid(self, event):
    #     if(self.spectral_dim_count==3 and self.tabOne.DECON==0):
    #         print('No deconvolution data available')
    #         self.cb_calc.SetValue(0)

    # if(self.cb_grid.GetValue()==1): #if plotting peaks
        #     self.conn_data=self.projection_service.connections
        #     if(len(self.conn_data)>0):
        #         for cn in self.conn_data:
        #             loc1=float(cn.f1)
        #             loc2=float(cn.f2)
        #             lab=cn.p1
        #             if lab == self.selected:
        #                 self.axes.annotate(lab,
        #                 xy=(loc1, loc2), xycoords='data',
        #                 xytext=(3, 3), textcoords='offset points', fontsize=12, color='r')
        #                 self.axes.scatter(loc1,loc2,s=50,marker='x',zorder=2, color='r')
        #             else:
        #                 self.axes.annotate(lab,
        #                 xy=(loc1, loc2), xycoords='data',
        #                 xytext=(3, 3), textcoords='offset points', fontsize=12, color='k')
        #                 self.axes.scatter(loc1,loc2,s=50,marker='x',zorder=2, color='k')
    #     self.draw_figure(scale='n')

    def on_cb_grid(self, event):
        if self._is_true_2d_spectrum() and event.GetEventObject() == self.cb_grid:
            if(self.cb_grid.GetValue()==1): #if plotting peaks
                if self.peaks_drawn == False:
                    self.peak_list = []
                    self.peak_list_names = []
                self.peak_list_locs = []
                self.peaks_text = []
                self.peaks_scatter = []
                self.canvas.mpl_connect('pick_event', self.on_pick_peak)
                self.conn_data = self._cached_peak_overlay()
                self.peak_list_locs, peak_labels = self._peak_points_for_overlay(self.conn_data, swap=True)
                if len(self.peak_list_locs) > 0:
                    for (loc1, loc2), lab in zip(self.peak_list_locs, peak_labels):
                        if lab:
                            self.peaks_text.append(self.axes.text(loc1-0.01,loc2-0.01,lab,fontsize=12))
                self.peaks_scatter.append(scatter_xy_points(self.axes, self.peak_list_locs, c='k',s=50,zorder=2,marker='x', picker=True, pickradius=10))
                self.peaks_drawn = True

            else:
                self.peaks_drawn=False
                for x in self.peaks_scatter:
                    x.set_visible(False)
                    x.remove()
                for x in self.peaks_text:
                    x.set_visible(False)
                    x.remove()

            self.canvas.draw()
            self.background_save(event)
        if self._is_true_2d_spectrum() and event.GetEventObject() == self.cb_calc:
            self.no_decon = self.twod_data_decon is None
            self._configure_2d_trace_visibility()
            self._blit_2d_traces()
        if self._is_3p_projection_case() or self.spectral_dim_count == 3:
            self.draw_figure()


    

    def on_cb_grid_auto(self, event):
        self.draw_figure()

    def on_slider_width(self, event):
        self.draw_figure()

    def _toolbar_decon(self, active):
        self.cb_calc.SetValue(bool(active))
        if self._is_true_2d_spectrum():
            # The Deconvolution tool controls the main 2D contour overlay as
            # well as the horizontal/vertical traces.  Previously this path
            # only blitted the traces, leaving the main axes unchanged until
            # the user pressed Draw.  Redraw the 2D spectrum immediately and
            # preserve the current zoom limits so the toolbar behaves as a
            # true show/hide toggle.
            self.no_decon = self.twod_data_decon is None
            if active and self.no_decon:
                self.cb_calc.SetValue(False)
                self.toolbar.set_decon_active(False)
            self.draw_2d(keepaxes=True)
            self.canvas.draw_idle()
            self._configure_2d_trace_visibility()
            self._blit_2d_traces()
        elif self._is_pseudo2d_projection_case():
            self.draw_figure(keepaxes=True)
            self.canvas.draw_idle()
        else:
            self.draw_figure()

    def _toolbar_peaks(self, active):
        self.cb_grid.SetValue(bool(active))
        if self._is_pseudo2d_projection_case():
            self.draw_figure(keepaxes=True)
            self.canvas.draw_idle()
            return
        evt = wx.CommandEvent(wx.wxEVT_CHECKBOX, self.cb_grid.GetId())
        evt.SetEventObject(self.cb_grid)
        self.on_cb_grid(evt)

    def _toolbar_contours(self):
        self.on_contour_button(None)

    def redraw_view(self):
        # Redraw must also release persistent Matplotlib Pan/Zoom.  Otherwise
        # the pseudo2D Full tools continue to reject canvas clicks after zoom.
        if self._is_pseudo2d_projection_case():
            self._cancel_projection_navigation()
        # Redraw is the projection-window reset action: clear the transient
        # cross-view peak selection, then redraw the current projection state.
        try:
            self.projection_service.clear_peak_selection(redraw_full3d=True)
        except TypeError:
            self.projection_service.clear_peak_selection()
        self.draw_figure(keepaxes=False)

    def on_draw_button(self, event):
        self.redraw_view()

    def on_shiftx_button(self, event):
        self.shiftXrun=1

        chain=self.textbox_chain.GetValue()
        shiftx2.runShiftx2(self.pdbfile,self.methlist,chain)
        self.draw_figure()


    def on_AutoFit_button(self, event):
        self.thresh=float(self.textbox0.GetValue())
        self.plotty=analslices1d(self.projection_service.peaks,self.thresh)
        self.cb_grid_auto.SetValue(1)
        self.draw_figure()


    def on_N_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()+1)
        self.draw_figure()

    def on_P_button(self, event):
        self.ComboBox1.SetSelection(self.ComboBox1.GetSelection()-1)
        self.draw_figure()



    def on_pick(self, event):
        if self._is_pseudo2d_projection_case() and self._handle_full_tool_click(event):
            return
        if self.drag_extra_plot == True:
            self.dragging_extra_plot = True
            self.drag_extra_plot = False
            self.origin = (event.xdata, event.ydata)
            return
        if self._is_true_2d_spectrum() and getattr(event, 'button', None) == 1:
            x, y = self._mouse_main_coordinates_2d(event)
            if x is not None and y is not None:
                self.last_mouse_x, self.last_mouse_y = x, y
                self.trackers_locked = not self.trackers_locked
                self._update_2d_traces(x, y)


    def on_release(self, event):
        if self.dragging_extra_plot == True:
            self.dragging_extra_plot = False
            self.drag_extra_plot = False
            self.background_save(None)



            
    def on_move(self, event):
        if self.dragging_extra_plot == True:
            for line in self.dt.extra_plots[-1].collections:
                line.remove()
            del(self.dt.extra_plots[-1])
            if event.xdata is None or event.ydata is None:
                return
            dx = self.origin[0]-event.xdata
            dy = self.origin[1]-event.ydata
            self.dt.extra_plots.append(self.axes.contour(self.dt.ucs[-1][1].ppm_scale()-dx,self.dt.ucs[-1][0].ppm_scale()-dy, self.dt.data[-1], cmap=self.dt.color_list[len(self.dt.extra_plots)], levels=self.dt.levels, linewidths=0.5))
            self.canvas.draw_idle()
            return
        if self._is_true_2d_spectrum() and not self.trackers_locked:
            x, y = self._mouse_main_coordinates_2d(event)
            if x is not None and y is not None:
                self._update_2d_traces(x, y)

    def on_horiz_button(self, event):
        self.horizontal = bool(event.GetEventObject().GetValue())
        self._configure_2d_trace_visibility()
        self._blit_2d_traces()

    def on_vert_button(self, event):
        self.vertical = bool(event.GetEventObject().GetValue())
        self._configure_2d_trace_visibility()
        self._blit_2d_traces()

    def _configure_2d_trace_visibility(self):
        if not hasattr(self, 'line1') or not hasattr(self, 'line_v'):
            return
        show_calc = bool(self.cb_calc.IsChecked()) and not self.no_decon and self.twod_data_decon is not None
        self.axes1D.yaxis.set_visible(bool(self.horizontal))
        self.axesV.xaxis.set_visible(bool(self.vertical))
        self.line1.set_visible(bool(self.horizontal))
        self.line_decon.set_visible(bool(self.horizontal and show_calc))
        self.cross_h.set_visible(bool(self.horizontal))
        self.line_v.set_visible(bool(self.vertical))
        self.line_v_decon.set_visible(bool(self.vertical and show_calc))
        self.cross_v.set_visible(bool(self.vertical))

    def _update_2d_traces(self, xppm, yppm, redraw=True):
        if self.twod_data is None or not hasattr(self, 'line_v'):
            return
        xs = numpy.asarray(self.uc1.ppms_scale)
        ys = numpy.asarray(self.uc0.ppms_scale)
        ix = int(numpy.argmin(numpy.abs(xs - float(xppm))))
        iy = int(numpy.argmin(numpy.abs(ys - float(yppm))))
        self.last_mouse_x = float(xs[ix])
        self.last_mouse_y = float(ys[iy])
        # Stored arrays are [indirect, direct]: row -> horizontal, column -> vertical.
        self.line1.set_data(xs, self.twod_data[iy, :])
        self.line_v.set_data(self.twod_data[:, ix], ys)
        dec = self.twod_data_decon
        if dec is not None:
            self.line_decon.set_data(xs, dec[iy, :])
            self.line_v_decon.set_data(dec[:, ix], ys)
        self.cross_h.set_ydata([ys[iy], ys[iy]])
        self.cross_v.set_xdata([xs[ix], xs[ix]])
        lo, hi = float(numpy.nanmin(self.twod_data)), float(numpy.nanmax(self.twod_data))
        if lo == hi: hi = lo + 1.0
        self.axes1D.set_ylim(lo, hi)
        self.axesV.set_xlim(lo, hi)
        self._configure_2d_trace_visibility()
        if redraw:
            self._blit_2d_traces()

    def _mouse_main_coordinates_2d(self, event):
        # twinx/twiny axes can be reported as event.inaxes, so test the main
        # axes pixel bbox and transform through the contour axes explicitly.
        if not hasattr(self, 'axes') or event.x is None or event.y is None:
            return None, None
        if not self.axes.bbox.contains(event.x, event.y):
            return None, None
        try:
            x, y = self.axes.transData.inverted().transform((event.x, event.y))
            return float(x), float(y)
        except Exception:
            return None, None

    def _dynamic_2d_trace_artists(self):
        return [self.line1, self.line_decon, self.line_v, self.line_v_decon,
                self.cross_h, self.cross_v]

    def _save_2d_trace_background(self):
        if not hasattr(self, 'line1'):
            return
        dynamic = self._dynamic_2d_trace_artists()
        visible = [artist.get_visible() for artist in dynamic]
        for artist in dynamic:
            artist.set_visible(False)
        self.canvas.draw()
        self._trace_blit_background = self.canvas.copy_from_bbox(self.fig.bbox)
        self.background = self._trace_blit_background
        for artist, state in zip(dynamic, visible):
            artist.set_visible(state)
        self._blit_2d_traces()

    def _blit_2d_traces(self):
        if self._trace_blit_background is None or not hasattr(self, 'line1'):
            self.canvas.draw_idle()
            return
        try:
            self.canvas.restore_region(self._trace_blit_background)
            for artist in self._dynamic_2d_trace_artists():
                if artist.get_visible():
                    artist.axes.draw_artist(artist)
            self.canvas.blit(self.fig.bbox)
        except Exception:
            self._trace_blit_background = None
            self.canvas.draw_idle()

    def on_pick_peak(self, event):
        if event.mouseevent.inaxes==self.axes:
            ind = event.ind
            print('picked:', self.peaks_text[ind[0]].get_text())
            self.ftol = 1e-10
            self.selected = self.peaks_text[ind[0]].get_text()
            self.peaks_text[ind[0]].set_color('r')
            update_colors = []
            for x in range(len(self.peaks_text)):
                if x == ind[0]:
                    self.peaks_text[x].set_color('r')
                    update_colors.append('r')

                else:
                    self.peaks_text[x].set_color('k')
                    update_colors.append('k')

            event.artist.set_color(update_colors)
            self.fuda_number = 0

            for key in self.line_fitting.plotting_resim_data.keys():
                if self.selected == str(key):
                    self.line_fitting.plot_fuda_fit(str(self.selected), 0, self.axes_proj, self.canvas)
            blit_artists(self.canvas, self.axes, self.background, [event.artist, self.peaks_text[ind[0]]])
            self.draw_bar()
            self.background_save(None)




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

    def on_save_plot_file(self, file):
        self.canvas.print_figure(file+'.pdf')
        #self.flash_status_message("Saved to %s" % file+'.pdf')

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

    def onFocus(self, event):
        print("Projection has focus!")

    #FGA added
    def onGetFile(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="",
            wildcard="PDB file (*.pdb)|*.pdb|" , style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #print path
            #fu=self.dirBox.GetValue()
            #print fu
            #print path.split(fu)
            #splitPath = path.split(cwd)
            #textBox.SetValue('.' + splitPath[1])
            print("You chose the following file(s):")
            print(path)
            textBox.SetValue(path)
        dlg.Destroy()
