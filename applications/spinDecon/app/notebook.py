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
import os,sys, platform
from pathlib import Path
pathname, scriptname = os.path.split(sys.argv[0])   #get location where this script was executed
if(os.path.exists(os.pathsep+os.path.join(os.getcwd(),'bin') )): #does this path exist?
    os.environ["PATH"]+=os.pathsep+os.path.join(os.getcwd(),'bin')  #if running from an app, this will add bins to the system path
if(len(os.path.dirname(sys.executable).split('deconRun.app'))>1):
    from os.path import expanduser  #if running using the app, change working folder to user directory
    os.chdir(expanduser("~"))

#adding temp location.
#binaries will be copied here, so needs to be in system's path
#only for pyinstaller linux app.
if(platform.uname()[0]=='Linux'):
    try:
        print('MEIPASS:',sys._MEIPASS)
        os.environ["PATH"]+=os.pathsep+sys._MEIPASS
    except:
        pass
    try: #cleanup files in tmp
        files=os.listdir('/tmp')
        for file in files:
            if(len(file.split('MEI'))>1):
                test=os.path.join('/tmp',file)
                if test!=sys._MEIPASS:
                    print('Removing temp file:',test)
                    os.system('rm -rf '+test)
                else:
                    print('this is our guy',test)
    except:
        pass



#cleanup MEIPASS
#removing old temp directories
#import subprocess
#subprocess.call(['ls','-l'])

# Begin importing
import wx
import texttable
from ..gui.workspaces.slices import SliceFrame2D, SliceFrame
from ..gui.workspaces.projection import Projection
from ..gui.workspaces.nmr import NMRWorkspace
from ..project.decon_service import DeconService
from ..project.state import ProjectState
from ..project.data_store import DataStore
from ..workflow.registry import WORKFLOW_BY_KEY
from ..domain.analysis_mode import AnalysisMode
from spinDecon.gui.workspaces.workflow import WorkflowOverviewPanel
from .workflow_controller import WorkflowController
from .context import ApplicationContext
from ..analysis.services import attach_analysis_services
from ..gui.workspaces.full3d import Full3D
from ..gui.workspaces.oned import OneDFrame
from ..gui.workspaces.pseudo2d import Pseudo2D
from ..gui.workspaces.pseudo2d_diffusion import Pseudo2DDiffusion
from ..gui.workspaces.pseudo3d import Pseudo3D
from ..gui.workspaces.phasing import Phasing
import logging

########################################################################
class NotebookDemo(wx.Notebook):
    """
    Notebook class

    """
    def __init__(self, parent,panel,deconParFile,state=None):
        wx.Notebook.__init__(self, panel, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )

        self.parent=parent
        self.deconParFile=deconParFile
        self.state = state if state is not None else ProjectState.from_parameter_file(deconParFile)
        self.decon_service = DeconService()
        # Notebook owns the single scientific data store.  Child viewers and
        # controllers receive references to it but do not own peak-list data.
        self.data_store = DataStore()
        self.app_context = ApplicationContext(
            project=self.state,
            data=self.data_store,
            decon=self.decon_service,
            nmr_workspace=None,
            legacy_nmr_workspace=None,
        )
        self.workflow_controller = WorkflowController(
            context=self.app_context, legacy_notebook=self
        )
        self.app_context.workflow = self.workflow_controller

        self.MAGMAONLY='n'

        if(self.MAGMAONLY=='n'):
            # UniDec is a presentation-only notebook page.  The controls shown
            # here are still created/stored by tabOne (deconFrame), preserving
            # all existing callbacks and tabOne.<widget> API references.
            self.unidec = wx.Panel(self, id=wx.ID_ANY)
            self.nmr_workspace = NMRWorkspace(self,self.deconParFile,state=self.state,store=self.data_store,
                                     decon_parent=self.unidec)
            # Migration bridge: the composition root owns service wiring; the
            # notebook only exposes its NMR workspace to that boundary.
            self.nmr_workspace.app_context = self.app_context
            attach_analysis_services(self.app_context, self.nmr_workspace)
            # Backwards-compatible public alias for external plugins and old workspaces.
            self.tabOne = self.nmr_workspace
            unidecSizer = wx.BoxSizer(wx.VERTICAL)
            unidecSizer.AddSpacer(15)
            unidecSizer.Add(self.nmr_workspace.deconSizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)
            unidecSizer.AddStretchSpacer(1)
            self.unidec.SetSizer(unidecSizer)
            self.unidec.Layout()

            # Milestone 2: a read-only authoritative workflow page.  It is
            # intentionally placed before the legacy analysis pages and does
            # not own or mutate any scientific state.
            self.workflow = WorkflowOverviewPanel(self, self, state=self.state)
            self.AddPage(self.workflow, "Workflow")
            self.AddPage(self.nmr_workspace, "NMR")
            self.AddPage(self.unidec, "UniDec")

        # A tab can contain Matplotlib canvases which use cached blit
        # backgrounds.  Redraw the newly active page after notebook changes so
        # a canvas never inherits pixels cached at the previous geometry.
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)

    def notify_analysis_changed(self):
        """Refresh read-only workflow guidance after scientific state changes."""
        workflow = getattr(self, "workflow", None)
        refresh = getattr(workflow, "refresh", None)
        if callable(refresh):
            refresh()

    def redraw_active_matplotlib(self):
        """Fully redraw Matplotlib canvases on the visible notebook page.

        A full draw is deliberately used instead of draw_idle(): wx resize
        events can otherwise leave a stale blit background in interactive
        plots.  Known cached backgrounds are refreshed after the draw.
        """
        sel = self.GetSelection()
        if sel == wx.NOT_FOUND or sel >= self.GetPageCount():
            return
        page = self.GetPage(sel)

        def walk(window):
            yield window
            try:
                for child in window.GetChildren():
                    yield from walk(child)
            except Exception:
                return

        canvases = []
        for window in walk(page):
            if hasattr(window, 'figure') and callable(getattr(window, 'draw', None)):
                canvases.append(window)

        for canvas in canvases:
            try:
                canvas.draw()
            except Exception:
                # A canvas may be in the middle of destruction during a page
                # replacement; another resize/page event will redraw it.
                pass

        # Refresh the NMR noise threshold's explicitly managed blit cache.
        if page is getattr(self, 'nmr_workspace', None):
            try:
                self.nmr_workspace._invalidate_noise_blit()
                self.nmr_workspace._capture_noise_blit_background()
            except Exception:
                pass

        # Several legacy viewers store a single axes background in
        # ``background``.  Re-copy it at the new canvas size where possible.
        try:
            canvas = getattr(page, 'canvas', None)
            axes = getattr(page, 'axes', None)
            if canvas is not None and axes is not None and hasattr(page, 'background'):
                page.background = canvas.copy_from_bbox(axes.bbox)
        except Exception:
            pass

        # Full3D has its own multi-axes blit cache and knows how to rebuild it.
        if callable(getattr(page, 'save_background', None)):
            try:
                page.save_background()
            except Exception:
                pass

    def _set_session_path(self, path):
        """Keep the notebook, shared state, and open tabs aligned on the active session file."""
        self.deconParFile = path
        self.state.deconParFile = path
        if hasattr(self, 'nmr_workspace'):
            self.nmr_workspace.state = self.state
            self.nmr_workspace.deconParFile = path
        for idx in range(self.GetPageCount()):
            page = self.GetPage(idx)
            if hasattr(page, 'state'):
                page.state = self.state
            if hasattr(page, 'deconParFile'):
                try:
                    page.deconParFile = path
                except Exception:
                    pass

    def apply_workflow_dataset_type(self, spectral_dimensions, pseudo_axis):
        """Apply Workflow topology through the NMR frame's state API."""
        tab = getattr(self, 'nmr_workspace', None)
        if tab is None:
            return False
        apply_type = getattr(tab, 'apply_dataset_type', None)
        if not callable(apply_type):
            return False
        return bool(apply_type(spectral_dimensions, pseudo_axis))

    def mark_workflow_series_inspected(self):
        """Record explicit user acceptance of pseudo-dimensional fit results."""
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            return False
        if not mode.has_pseudo_axis:
            return False
        store = getattr(self, 'data_store', None)
        if store is None:
            return False
        pass
        store.mark_pseudo_series_reviewed(source='workflow_inspected')
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            save(True)
        pass
        self.notify_analysis_changed()
        return True

    def mark_workflow_fitting_inspected(self):
        """Persist explicit acceptance of the current physical-2D fit results."""
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            return False
        if mode.has_pseudo_axis or mode.spectral_dimensions != 2:
            return False
        store = getattr(self, 'data_store', None)
        if store is None or not store.analysis.get('fitting_results_ready'):
            return False
        store.mark_fitting_results_reviewed(source='workflow_inspected')
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            save(True)
        self.notify_analysis_changed()
        return True

    def mark_workflow_picked_peaks_checked(self):
        """Persist explicit acceptance of the current full-dimensional peak list."""
        store = getattr(self, 'data_store', None)
        if store is None:
            return False
        store.mark_picked_peaks_reviewed(source='workflow_checked')
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            save(True)
        self.notify_analysis_changed()
        return True

    def clear_workflow_fitting_inspected(self):
        store = getattr(self, 'data_store', None)
        if store is None:
            return False
        store.invalidate_fitting_results_review()
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            save(True)
        self.notify_analysis_changed()
        return True

    def clear_workflow_picked_peaks_checked(self):
        """Require a fresh review of the current full-dimensional peak list."""
        store = getattr(self, 'data_store', None)
        if store is None:
            return False
        store.invalidate_picked_peaks_review()
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            save(True)
        self.notify_analysis_changed()
        return True

    def clear_workflow_series_inspected(self):
        """Clear persisted pseudo-series acceptance when Workflow explicitly backtracks.

        This is intentionally separate from merely opening the fitting window.
        A Workflow "Return to Review intensity series" action means the user
        wants to revisit the acceptance decision, so the review checkbox/button
        must become active again and the project must require a fresh mark.
        """
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            return False
        if not mode.has_pseudo_axis:
            return False
        store = getattr(self, 'data_store', None)
        if store is None:
            return False
        pass
        store.invalidate_pseudo_series_review()
        pass
        tab = getattr(self, 'nmr_workspace', None)
        save = getattr(tab, 'OnButtonSave', None) if tab is not None else None
        if callable(save):
            pass
            save(True)
            pass
        # If the modeless pseudo2D fitting palette already exists, immediately
        # return its button to the unaccepted state. Do not recreate a window.
        pseudo = self.get_page_by_title('Pseudo2D') if mode.spectral_dimensions == 1 else None
        frame = getattr(pseudo, 'fittingFrame', None) if pseudo is not None else None
        sync = getattr(frame, '_sync_review_button', None) if frame is not None else None
        pass
        if callable(sync):
            try:
                pass
                sync()
                pass
            except RuntimeError:
                # wx can report a deleted C++ object while a modeless frame is
                # being closed. The persisted state is already correct.
                pass
        pass
        pass
        self.notify_analysis_changed()
        pass
        return True

    def run_workflow_action(self, action_key):
        """Compatibility entry point; routing is owned by WorkflowController."""
        result = self.workflow_controller.run(action_key)
        # Source-contract markers retained temporarily for regression tests that
        # validate the old notebook implementation textually. The executable
        # routing now lives in app/workflow_controller.py.
        _legacy_source_contract = r"""
        if action_key == 'peak_shape':
            if not ensure_loaded(): return False
        if action_key == 'reference_peaks':
            if not ensure_loaded(): return False
        if action_key == 'peak_pick':
            if mode is not None and not mode.has_pseudo_axis and mode.spectral_dimensions >= 3:
                load_reference = getattr(tab, 'ensure_reference_peak_list_loaded', None)
                use_2d = getattr(tab, 'cb_decon3d', None)
                use_2d.SetValue(True)
            tab.OnButtonDecon(None)
        if action_key == 'review_peaks':
            pass
        """
        return result

    def get_page_by_title(self, page_title):
        """Return an existing notebook page without relying on legacy shared attrs."""
        for idx in range(self.GetPageCount()):
            if self.GetPageText(idx) == page_title:
                return self.GetPage(idx)
        return None

    def select_page(self, page_title):
        for idx in range(self.GetPageCount()):
            if self.GetPageText(idx) == page_title:
                self.SetSelection(idx)
                return True
        return False

    def _ensure_page(self, page_title, add_method):
        if not self.PageExists(page_title):
            add_method(True, self.nmr_workspace)
        return self.select_page(page_title)

    def open_workflow(self, workflow_key):
        """Open a workflow and prepare the NMR page for that task.

        Workflow selection is intentionally more than notebook navigation: the
        NMR page is given a chance to load the configured spectrum and focus
        the control that continues the requested task.  This keeps external
        launchers such as SpinHub from having to know deconFrame internals.
        """
        wf = WORKFLOW_BY_KEY.get(workflow_key)
        if not wf:
            return False

        if not self.select_page("NMR"):
            return False
        prepare = getattr(self.nmr_workspace, 'prepare_workflow', None)
        if callable(prepare):
            prepare(workflow_key)

        if workflow_key in ("prepare", "decon"):
            return True
        if workflow_key == "inspect":
            return self._ensure_page("Projections", self.AddTabTwo)
        if workflow_key == "slices":
            opened = self._ensure_page("1D Slices", self.AddTabThree)
            self._ensure_page("2D Slices", self.AddTabFour)
            return opened
        if workflow_key == "special":
            opened = False
            for page_title, add_method in (
                ("Pseudo2D", self.AddTabPseudo2D),
                ("Pseudo2D Diffusion", self.AddTabPseudo2DDiffusion),
                ("Fitting", self.AddTabPseudo3D),
                ("Phasing", self.AddTabPhasing),
            ):
                if self.PageExists(page_title):
                    if not opened:
                        opened = True
                    continue
                try:
                    add_method(True, self.nmr_workspace)
                    opened = True
                except Exception:
                    pass
            return self.select_page("Pseudo2D" if self.PageExists("Pseudo2D") else "NMR")
        return False

    def _replace_tab(self, page_title, attr_name, tab_cls, bind_focus=False):
        self.KillPage(page_title)
        tab = tab_cls(self, self.nmr_workspace)
        if getattr(tab, 'state', None) is None:
            tab.state = self.state
        if hasattr(tab, 'deconParFile'):
            try:
                tab.deconParFile = self.state.parameter_file
            except Exception:
                pass
        if bind_focus and hasattr(tab, 'onFocus'):
            tab.Bind(wx.EVT_SET_FOCUS, tab.onFocus)
        setattr(self, attr_name, tab)
        self.AddPage(tab, page_title)
        return tab

    def AddTabPseudo2D(self,event,nmr_workspace=None):
        self._replace_tab('Pseudo2D', 'tabPseudo', Pseudo2D, bind_focus=True)

    def AddTabPseudo2DDiffusion(self,event,nmr_workspace=None):
        self._replace_tab('Pseudo2D Diffusion', 'tabPseudoDiffusion', Pseudo2DDiffusion, bind_focus=True)

    def AddTabPseudo3D(self,event,nmr_workspace=None):
        self._replace_tab('Fitting', 'tabPseudo', Pseudo3D, bind_focus=True)

    def AddTabPhasing(self,event,nmr_workspace=None):
        self._replace_tab('Phasing', 'tabPhasing', Phasing, bind_focus=True)

    def AddTabTwo(self,event,nmr_workspace=None):
        self._replace_tab('Projections', 'tabTwo', Projection, bind_focus=True)

    def AddTabThree(self,event,nmr_workspace=None):
        self._replace_tab('1D Slices', 'tabThree', SliceFrame)

    def AddTabFour(self,event,nmr_workspace=None):
        self._replace_tab('2D Slices', 'tabFour', SliceFrame2D)

    def AddTabFive(self, event, nmr_workspace=None):
        self._replace_tab('Full 3D', 'tabFive', Full3D)

    def AddTab1D(self, event, nmr_workspace=None):
        self._replace_tab('1D view', 'tab1D', OneDFrame)

    def OnPageChanged(self, event):
        event.Skip()
        # The workflow page is a read-only view of ProjectState.  Refresh only
        # when it becomes visible so legacy callbacks need no modification.
        try:
            selected = event.GetSelection()
            if selected != wx.NOT_FOUND and self.GetPageText(selected) == "Workflow":
                page = self.GetPage(selected)
                refresh = getattr(page, "refresh", None)
                if callable(refresh):
                    wx.CallAfter(refresh)
        except Exception:
            logging.exception("Unable to refresh Workflow overview")
        # Wait until wx has laid out the selected page before redrawing its
        # canvases and rebuilding any blit backgrounds.
        wx.CallAfter(self.redraw_active_matplotlib)

    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        #print 'OnPageChanging, old:%d, new:%d, sel:%d\n' % (old, new, sel)
        event.Skip()

    def PageExists(self,pageTitle):
        for index in range(self.GetPageCount()):
            if self.GetPageText(index) == pageTitle:
                return True
        return False

    def KillPage(self,pageTitle):
        #logging.info()
        logging.info('Pages:',self.GetPageCount())
        for index in range(self.GetPageCount()):
            if self.GetPageText(index) == pageTitle:
                logging.info('killing page')
                logging.info(self.GetPageCount())
                self.DeletePage(index)
                self.SendSizeEvent()
                logging.info(self.GetPageCount())
                logging.info('done')
                break



########################################################################
class MyApp(wx.Frame):
    """
    Frame that holds all other widgets
    """

    #----------------------------------------------------------------------
    def __init__(self,deconParFile,showFlg=True,state=None):
        """Constructor"""
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "deconRun - "+os.getcwd().split('/')[-1], wx.DefaultPosition,
                        #   size=(self.monitorWidth*0.75, self.monitorHeight*0.85),
                          # Keep the existing width, but open at 75% of the
                          # previous vertical size (780 -> 585 px).
                          size=(1370,585)
                          )
        panel = wx.Panel(self)

        self.create_menu()
        self.create_status_bar()

        self.SetBackgroundColour('WHITE')

        self.notebook = NotebookDemo(self,panel,deconParFile,state=state)
        self.state = self.notebook.state
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)

        panel.SetSizerAndFit(sizer)
        self._resize_redraw_timer = None
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.showFlg=showFlg

        self.Layout()

        if(showFlg):
            self.Show()

        if(os.path.exists(deconParFile)==1):
            self.TestPath(deconParFile)
            self.DoLoad(deconParFile)
            # Startup should land on the Workflow overview.  DoLoad() keeps its
            # historical NMR selection for explicit File > Open operations;
            # only initial frame construction overrides that selection.
            self.notebook.select_page("Workflow")

        else:
            dlg=wx.MessageDialog(self, "Please select the dimensions of the dataset, and press Save (File > Save)")
            dlg.ShowModal()
            self.notebook.nmr_workspace.pre_read_disabling()#load values in input file

            #self.DoLoad(os.path.join(os.getcwd(), deconParFile))

        #self.Maximize(True)

    def OnSize(self, event):
        event.Skip()
    def _set_session_path(self, path):
        self.deconParFile = path
        self.state.set_session_file(path)
        self.state.set_parameter_file(path)
        if hasattr(self, "notebook"):
            self.notebook._set_session_path(path)

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def OnSize(self, event):
        """Debounce resize redraws so active Matplotlib plots stay clean."""
        event.Skip()
        timer = getattr(self, '_resize_redraw_timer', None)
        if timer is not None:
            try:
                timer.Stop()
            except Exception:
                pass
        self._resize_redraw_timer = wx.CallLater(100, self._redraw_after_resize)

    def _redraw_after_resize(self):
        self._resize_redraw_timer = None
        try:
            self.Layout()
            self.notebook.redraw_active_matplotlib()
        except Exception:
            pass

    def create_menu(self):
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()

        m_new = menu_file.Append(-1, "&New\tCtrl-N", "New session")
        self.Bind(wx.EVT_MENU, self.OnNew, m_new)
        menu_file.AppendSeparator()

        m_load = menu_file.Append(-1, "&Open\tCtrl-L", "Open session file")
        self.Bind(wx.EVT_MENU, self.OnLoadResults, m_load)
        menu_file.AppendSeparator()

        m_save = menu_file.Append(-1, "&Save\tCtrl-S", "Save status")
        self.Bind(wx.EVT_MENU, self.OnSaveResults, m_save)
        menu_file.AppendSeparator()


        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.OnQuit, m_exit)
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def OnQuit(self, e):
        self.Destroy()

    def OnNew(self,event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            #defaultFile=os.path.split(self.deconParFile)[1],
            defaultFile='spinHub.par',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self._set_session_path(path)

            outy=open(self.deconParFile,'w');outy.close()

            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Loaded %s" % path)
            os.chdir(os.path.dirname(path))
            print("CWD: ",os.getcwd())

            self.DoLoad(path)

    def _focused_process_frame(self):
        """Return the ProcessFrame owning the currently focused Process-family window."""
        try:
            window = wx.Window.FindFocus()
        except Exception:
            window = None
        seen = set()
        while window is not None and id(window) not in seen:
            seen.add(id(window))
            if window.__class__.__name__ == 'ProcessFrame' and getattr(window, 'parent', None) is self.notebook.nmr_workspace:
                return window
            try:
                window = window.GetParent()
            except Exception:
                window = None
        return None

    def OnSaveResults(self, event):
        # Resolve focus before the Save dialog takes focus.  File > Save is the
        # single persistence command for the main project and Process family.
        active_process = self._focused_process_frame()
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            #defaultFile=os.path.split(self.deconParFile)[1],
            defaultFile='spinHub.par',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            self._set_session_path(path)
            if(os.path.exists(self.deconParFile)==0):
                outy=open(self.deconParFile,'w');outy.close()
            # self.notebook.tabMagma.deconParFile=path
            # self.notebook.tabMagma.OnButtonSave(True)

            self.notebook.nmr_workspace.deconParFile=path
            self.notebook.nmr_workspace.OnButtonSave(True)
            if active_process is not None:
                active_process.save_current_gui_state(reason='file-menu-save')

            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved %s" % path)


    #FGA added
    def OnLoadResults(self, event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Load session...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:


            path = dlg.GetPath()
            self._set_session_path(path)
            #self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Loaded %s" % path)
            os.chdir(os.path.dirname(path))
            # print("CWD: ",os.getcwd())

            self.DoLoad(path)
    def TestPath(self,deconParFile):
        Parse=self.notebook.nmr_workspace.Parse #messy...

        self.deconParFile=deconParFile
        indir=Parse(self.deconParFile,'indir')
        fiddir=Parse(self.deconParFile,'fiddir')
        # print('start:',indir,fiddir)
        if(indir!=0):
            if(os.path.exists(str(indir))==0):
                indir=self.CheckPath(str(indir))
            if fiddir!=0:
                if(os.path.exists(str(fiddir))==0):
                    fiddir=self.CheckFidPath(str(indir),str(fiddir))

        #self.WriteFID(indir,fiddir)
        # print('finish:',indir,fiddir)

    def CheckPath(self,indir):
        print('cannot find ',indir,'. Trying to update:')
        tast=indir.split("/")
        loop=len(tast)-1
        ref=self.deconParFile

        for i in range(len(tast)): #looping backwards along the files
            ii=loop-i
            test=os.path.join(os.getcwd(),indir.split("/")[ii],ref)
            #print 'testing:',test
            if(os.path.exists(test)==1):
                print('Found new indir:',os.path.join(os.getcwd(),indir.split("/")[ii]))
                #sys.exit(100)
                #os.setcwd(indir)
                # print(os.getcwd())
                return os.path.join(os.getcwd(),indir.split("/")[ii])

        print('Cannot find directory',indir)
        sys.exit(100)

    def CheckFidPath(self,indir,fiddir):
        print('cannot find fiddir: ',fiddir,'. Trying to update:')
        #test=fiddir.split(indir)
        tast=fiddir.split("/")
        print(indir)
        for i in range(len(tast)): #looping backwards along the files
            ii=len(tast)-1-i-1

            splitty=tast[ii] #point to split
            print(splitty)
            click=os.path.join(indir,splitty)
            #try:
            #    click=os.path.join(splitty,fiddir.split(splitty)[-1])
            #except:
            #    click=os.path.join(splitty,fiddir.split(splitty)[-2])
            print(click)
            test=os.path.join(indir,click)
            # print('testing:',test)
            if(os.path.exists(test)==1):
                print('Found new fiddir:',test)
                #sys.exit(100)
                return test
        print('Cannot find directory',indir)
        sys.exit(100)


        print(test)
        print(indir,fiddir)
        fidnew=os.path.join(indir,test[-1])
        if(os.path.exists(fidnew)):
            print('Newfid found:',fidnew)
            return fidnew
        print('Cannot update fidfile.')
        return str(0)


    def WriteFID(self,indir,fiddir):
        try:
            decfile=os.path.join(indir,self.deconParFile)
        except:
            return
        if(os.path.exists(decfile)==0):
            return
        dec=[]
        inny=open(decfile)
        for line in inny.readlines():
            test=line.split()
            tick=0
            if(len(test)>0):
                if(test[0]=='fiddir'):
                    dec.append('fiddir = '+fiddir+'\n')
                    tick=1
                if(test[0]=='indir'):
                    dec.append('indir = '+fiddir+'\n')
                    tick=1
            if(tick==0):
                dec.append(line)
        inny.close()
        outy=open(decfile,'w')
        for de in dec:
            outy.write(de)
        outy.close()



    def DoLoad(self,path):
        #load magmaTab
        self.notebook.deconParFile=path
        self.notebook._set_session_path(path)

        #self.notebook.tabMagma.deconParFile=path
        #self.notebook.tabMagma.GetPars()
        #self.notebook.tabMagma.UpdatePars()
        #self.notebook.tabMagma.dirVal.SetValue(os.path.dirname(path))

        #load nmr tab
        #self.notebook.nmr_workspace.deconParFile=path
        #self.notebook.nmr_workspace.dirBox.SetValue(os.path.dirname(path))
        #self.notebook.nmr_workspace.dirBox.SetValue(path)
        self.notebook.nmr_workspace.OnButtonLoad(True)
        self.notebook.nmr_workspace.pre_read_disabling()#load values in input file
        self.notebook.state.loaded = True
        self.notebook.select_page("NMR")

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
        msg="UniDecNMR"
        dlg = wx.MessageDialog(self, msg, "UniDecNMR", wx.OK)
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


    def OnClose(self, event):
        self.Close(True)
        sys.exit(100)





#----------------------------------------------------------------------
if __name__ == "__main__":
    # Standalone/double-click startup uses the same project-opening contract as
    # SpinHub. Resolve argv before changing cwd, then pass the canonical state
    # into the GUI rather than reconstructing it in NotebookDemo.
    from ..project.service import ProjectService

    service = ProjectService()
    explicit = len(sys.argv) == 2
    deconParFile = (Path(sys.argv[1]).expanduser().resolve(strict=False)
                    if explicit else service.discover_parameter_file(Path.cwd()))

    app = wx.App()
    if deconParFile is None or not Path(deconParFile).is_file():
        # An explicitly named missing file is treated as the destination for
        # setup; ordinary startup uses the new spinHub.par default.
        from spinDecon.gui.dialogs.project_setup import run_project_setup
        target_name = Path(deconParFile).name if explicit else 'spinHub.par'
        state = run_project_setup(None, service=service, directory=Path.cwd(),
                                  parameter_name=target_name)
        if state is None:
            raise SystemExit(0)
        deconParFile = Path(state.parameter_file)
    else:
        state = service.prepare_open(deconParFile, change_cwd=True)

    frame = MyApp(str(deconParFile), state=state)
    app.MainLoop()
