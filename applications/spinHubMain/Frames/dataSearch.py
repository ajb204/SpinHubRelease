#!/usr/bin/env python3
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
import wx,string,os,numpy,sys
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt          #plotting routines from matplotlib
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

from wx.lib.mixins.listctrl import ColumnSorterMixin

from pathlib import Path


############################################################################
# Frame for dataSearch
#

matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8


class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

    def Update(self,dicty):
        ColumnSorterMixin.__init__(self,len(list(dicty.keys())))
        self.itemDataMap = dicty
        #print(dicty[0])

    def CustColumnSorter(self, key1, key2):
        col = self._col
        print (key1,key2,col)
        ascending = self._colSortFlag[col]
        ascending=1
        item1 = self.itemDataMap[key1][col]
        item2 = self.itemDataMap[key2][col]

        self.num_cols=[0,2,3,]
        if col in self.num_cols:
            #just convert them to float, cmp do comparing float well
            item1 = float(item1)
            item2 = float(item2)

        cmpVal = cmp(item1, item2)

        # If the items are equal then pick something else to make the sort value unique
        if cmpVal == 0:
            cmpVal = cmp(*self.GetSecondarySortValues(col, key1, key2))

        if ascending:
            return cmpVal
        else:
            return -cmpVal

    def GetColumnSorter(self):
        return self.CustColumnSorter

    



class dataSearch(wx.Panel):
    """Responsive dataset browser with background recursive discovery."""
    # Keep the browser table focused on comparison/selection.  Rich metadata
    # belongs in the persistent inspector on the right.
    COLUMNS = (
        ('Dataset', 210), ('Pulse sequence', 165), ('Dim', 55),
        ('Source', 105), ('Status', 180),
    )

    def __init__(self, parent):
        import threading
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.SetClientSize(wx.Size(1150, 560))
        self.datasets = []
        self._row_datasets = []
        self._scan_thread = None
        self._scan_cancel = threading.Event()
        self._sort_column = 'Dataset'
        self._sort_descending = False
        from spinHubMain.core import BrowserLocation
        self.browser_location = BrowserLocation(Path.cwd())

        self.search = wx.SearchCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.search.ShowCancelButton(True)
        self.search.SetDescriptiveText('Search datasets, paths, sequences...')
        self.source_filter = wx.Choice(self, choices=['All', 'Decon project', 'Acquisition'])
        self.source_filter.SetSelection(0)
        self.status_filter = wx.Choice(self, choices=['All'])
        self.status_filter.SetSelection(0)
        self.search.Bind(wx.EVT_TEXT, self.OnFilterChanged)
        self.source_filter.Bind(wx.EVT_CHOICE, self.OnFilterChanged)
        self.status_filter.Bind(wx.EVT_CHOICE, self.OnFilterChanged)

        self.summary = wx.StaticText(self, label='No datasets scanned yet')
        summary_font = self.summary.GetFont()
        summary_font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.summary.SetFont(summary_font)

        # The browser location is explicit state.  Opening UniDecNMR may change
        # process CWD for legacy compatibility; that must not move SpinHub.
        self.location_text = wx.TextCtrl(self, value=str(self.browser_location.root), style=wx.TE_PROCESS_ENTER)
        self.location_text.Bind(wx.EVT_TEXT_ENTER, self.OnLocationEnter)
        self.back_button = wx.Button(self, label='Back')
        self.back_button.Bind(wx.EVT_BUTTON, self.OnBack)
        self.up_button = wx.Button(self, label='Up')
        self.up_button.Bind(wx.EVT_BUTTON, self.OnUp)
        self.browse_button = wx.Button(self, label='Browse...')
        self.browse_button.Bind(wx.EVT_BUTTON, self.OnBrowse)
        location_line = wx.BoxSizer(wx.HORIZONTAL)
        location_line.Add(wx.StaticText(self, label='Location:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        location_line.Add(self.location_text, 1, wx.RIGHT, 6)
        location_line.Add(self.back_button, 0, wx.RIGHT, 4)
        location_line.Add(self.up_button, 0, wx.RIGHT, 4)
        location_line.Add(self.browse_button, 0)
        self._update_location_controls()

        filters = wx.BoxSizer(wx.HORIZONTAL)
        filters.Add(wx.StaticText(self, label='Search:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        filters.Add(self.search, 1, wx.RIGHT, 10)
        filters.Add(wx.StaticText(self, label='Source:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        filters.Add(self.source_filter, 0, wx.RIGHT, 10)
        filters.Add(wx.StaticText(self, label='Status:'), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        filters.Add(self.status_filter, 0)

        # The three work areas live in splitter panes so users can tune the
        # workspace.  The inspector is deliberately protected with a useful
        # minimum width: it contains the next actions for the selected dataset.
        self.content_splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE | wx.SP_3D)
        self.work_splitter = wx.SplitterWindow(self.content_splitter, style=wx.SP_LIVE_UPDATE | wx.SP_3D)
        self.nav_panel = wx.Panel(self.content_splitter)
        self.table_panel = wx.Panel(self.work_splitter)
        self.inspector_panel = wx.Panel(self.work_splitter)

        self.tree = wx.TreeCtrl(self.nav_panel, style=wx.TR_HAS_BUTTONS | wx.TR_LINES_AT_ROOT | wx.TR_SINGLE)
        self.tree.SetMinSize(wx.Size(180, 220))
        self.tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnTreeSelect)

        self.lc = wx.ListCtrl(self.table_panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for idx, (label, width) in enumerate(self.COLUMNS):
            self.lc.InsertColumn(idx, label, width=width)
        self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.Select)
        self.lc.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnOpen)
        self.lc.Bind(wx.EVT_LIST_COL_CLICK, self.OnColumnClick)

        # Project dashboard: make the processing pipeline visible instead of
        # hiding resource state in a text dump.  The compact details box below
        # remains useful for spectrometer metadata and diagnostics.
        self.stage_panels = []
        self.stage_status = []
        self.stage_path = []
        self.stage_buttons = []
        dashboard = wx.BoxSizer(wx.VERTICAL)
        for title in ('1. Raw data', '2. Spectrum', '3. Peaks / analysis'):
            box = wx.StaticBoxSizer(wx.VERTICAL, self.inspector_panel, title)
            # wx requires controls owned by a StaticBoxSizer to be children of
            # the StaticBox itself (newer debug builds enforce this strictly).
            box_parent = box.GetStaticBox()
            status = wx.StaticText(box_parent, label='—')
            font = status.GetFont(); font.SetWeight(wx.FONTWEIGHT_BOLD); status.SetFont(font)
            path_label = wx.StaticText(box_parent, label='')
            path_label.Wrap(300)
            button = wx.Button(box_parent, label='—', size=(180, -1))
            button.SetMinSize(wx.Size(180, -1))
            button.Disable()
            button.Bind(wx.EVT_BUTTON, self.OnStageAction)
            box.Add(status, 0, wx.BOTTOM, 3); box.Add(path_label, 0, wx.EXPAND | wx.BOTTOM, 5); box.Add(button, 0, wx.EXPAND)
            dashboard.Add(box, 0, wx.EXPAND | wx.BOTTOM, 6)
            self.stage_status.append(status); self.stage_path.append(path_label); self.stage_buttons.append(button)

        self.details = wx.TextCtrl(self.inspector_panel, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP | wx.HSCROLL)
        self.details.SetMinSize(wx.Size(300, 110))
        self.open_button = wx.Button(self.inspector_panel, label='Open analysis')
        self.open_button.Disable(); self.open_button.Bind(wx.EVT_BUTTON, self.OnOpen)
        self.locate_raw_button = wx.Button(self, label='Locate raw data...')
        self.locate_raw_button.Disable(); self.locate_raw_button.Bind(wx.EVT_BUTTON, self.OnLocateRaw)
        self.locate_spectrum_button = wx.Button(self, label='Locate spectrum...')
        self.locate_spectrum_button.Disable(); self.locate_spectrum_button.Bind(wx.EVT_BUTTON, self.OnLocateSpectrum)
        self.refresh_button = wx.Button(self.inspector_panel, label='Refresh')
        self.refresh_button.Bind(wx.EVT_BUTTON, self.OnRefresh)
        self.cancel_button = wx.Button(self.inspector_panel, label='Cancel scan')
        self.cancel_button.Disable(); self.cancel_button.Bind(wx.EVT_BUTTON, self.OnCancelScan)
        self.progress = wx.Gauge(self, range=100, style=wx.GA_HORIZONTAL)
        self.progress.Hide()
        self.scan_status = wx.StaticText(self, label='')

        right = wx.BoxSizer(wx.VERTICAL)
        right.Add(wx.StaticText(self.inspector_panel, label='Dataset workflow'), 0, wx.BOTTOM, 6)
        right.Add(dashboard, 0, wx.EXPAND | wx.BOTTOM, 4)
        right.Add(wx.StaticText(self.inspector_panel, label='Details'), 0, wx.TOP | wx.BOTTOM, 4)
        right.Add(self.details, 1, wx.EXPAND)
        actions = wx.BoxSizer(wx.HORIZONTAL)
        actions.Add(self.open_button, 0, wx.RIGHT, 8)
        actions.Add(self.refresh_button, 0, wx.RIGHT, 8)
        actions.Add(self.cancel_button, 0)
        right.Add(actions, 0, wx.TOP, 8)

        main = wx.BoxSizer(wx.VERTICAL)
        # Keep the project summary in normal sizer flow and right-aligned.
        # On newer wxPython builds this also avoids it being painted over by
        # controls in the workflow pane during resize/layout.
        summary_line = wx.BoxSizer(wx.HORIZONTAL)
        summary_line.Add(self.locate_spectrum_button, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        summary_line.Add(self.locate_raw_button, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        summary_line.Add(self.summary, 0, wx.ALIGN_CENTER_VERTICAL)
        summary_line.AddStretchSpacer(1)
        main.Add(summary_line, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 8)
        main.Add(location_line, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)
        main.Add(filters, 0, wx.EXPAND | wx.ALL, 6)
        navigation = wx.BoxSizer(wx.VERTICAL)
        navigation.Add(wx.StaticText(self.nav_panel, label='Project navigator'), 0, wx.BOTTOM, 6)
        navigation.Add(self.tree, 1, wx.EXPAND)
        self.nav_panel.SetSizer(navigation)

        table_sizer = wx.BoxSizer(wx.VERTICAL)
        table_sizer.Add(self.lc, 1, wx.EXPAND)
        self.table_panel.SetSizer(table_sizer)

        self.inspector_panel.SetSizer(right)
        self.inspector_panel.SetMinSize(wx.Size(320, -1))

        # Nested splitters give independent navigator/table and table/inspector
        # resizing.  The right pane gets a minimum width and the table absorbs
        # most resizing, so the workflow remains visible on laptop screens.
        self.content_splitter.SetMinimumPaneSize(160)
        self.work_splitter.SetMinimumPaneSize(320)
        self.content_splitter.SetSashGravity(0.0)
        self.work_splitter.SetSashGravity(1.0)
        self.work_splitter.SplitVertically(self.table_panel, self.inspector_panel, -350)
        self.content_splitter.SplitVertically(self.nav_panel, self.work_splitter, 220)
        main.Add(self.content_splitter, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        progress_line = wx.BoxSizer(wx.HORIZONTAL)
        progress_line.Add(self.progress, 0, wx.RIGHT, 8)
        progress_line.Add(self.scan_status, 1, wx.ALIGN_CENTER_VERTICAL)
        main.Add(progress_line, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)
        self.SetSizer(main)
        self.StartScan()

    def Refresh(self):
        """Compatibility entry point: refresh means start a new background scan."""
        self.StartScan()

    def StartScan(self):
        import threading
        if self._scan_thread and self._scan_thread.is_alive():
            self._scan_cancel.set()
        self._scan_cancel = threading.Event()
        selected_key = self._selected_key()
        self._pending_selection_key = selected_key
        self.refresh_button.Disable(); self.cancel_button.Enable()
        self.progress.Show(); self.progress.Pulse()
        self.scan_status.SetLabel('Scanning recursively...')
        self.Layout()
        root = self.browser_location.root
        cancel = self._scan_cancel

        def worker():
            from spinHubMain.core import scan_tree
            def progress(p):
                wx.CallAfter(self._on_scan_progress, p, cancel)
            result = scan_tree(root, cancel_event=cancel, progress_callback=progress)
            wx.CallAfter(self._on_scan_complete, result, cancel)
        self._scan_thread = threading.Thread(target=worker, name='SpinHubScanner', daemon=True)
        self._scan_thread.start()

    def _on_scan_progress(self, progress, token):
        if token is not self._scan_cancel:
            return
        self.progress.Pulse()
        self.scan_status.SetLabel(
            f'Scanning... {progress.visited} entries; '
            f'{progress.acquisitions} acquisitions; {progress.projects} projects')

    def _on_scan_complete(self, result, token):
        if token is not self._scan_cancel:
            return
        from spinHubMain.core import resolve_datasets
        self.refresh_button.Enable(); self.cancel_button.Disable(); self.progress.Hide()
        if result.cancelled:
            self.scan_status.SetLabel('Scan cancelled.')
            self.Layout(); return
        self.datasets = resolve_datasets(list(result.acquisitions), list(result.projects))
        statuses = sorted({d.status_text for d in self.datasets})
        previous = self.status_filter.GetStringSelection() or 'All'
        self.status_filter.SetItems(['All'] + statuses)
        self.status_filter.SetStringSelection(previous if previous in ['All'] + statuses else 'All')
        from spinHubMain.core import browser_summary
        summary = browser_summary(self.datasets)
        self.summary.SetLabel(
            f"{summary['acquisitions']} acquisitions    {summary['projects']} projects    "
            f"{summary['ready']} ready    {summary['attention']} need attention")
        self.scan_status.SetLabel(
            f'{len(self.datasets)} dataset views from {len(result.acquisitions)} discovered acquisitions and '
            f'{len(result.projects)} Decon projects')
        self.ApplyFilters(self._pending_selection_key)
        self.Layout()

    def ApplyFilters(self, restore_key=None):
        from spinHubMain.core import filter_datasets, sort_datasets, row_for
        if restore_key is None:
            restore_key = self._selected_key()
        source = self.source_filter.GetStringSelection() or 'All'
        status = self.status_filter.GetStringSelection() or 'All'
        rows = filter_datasets(self.datasets, self.search.GetValue(), source, status)
        rows = sort_datasets(rows, self._sort_column, self._sort_descending)
        self._row_datasets = list(rows)
        self._populate_navigation(self._row_datasets)
        self.lc.DeleteAllItems(); restore = -1
        for idx, dataset in enumerate(self._row_datasets):
            row = row_for(dataset)
            vals = (row.name, row.sequence, row.dimension, row.source, row.status)
            self.lc.InsertItem(idx, vals[0])
            for col, val in enumerate(vals[1:], 1): self.lc.SetItem(idx, col, str(val))
            if restore_key and self._dataset_key(dataset) == restore_key: restore = idx
        if restore >= 0:
            self.lc.Select(restore); self.lc.Focus(restore)
        elif self._row_datasets:
            self.details.SetValue('Select a dataset to inspect its resources and available actions.')
            self.open_button.Disable()
            self.locate_raw_button.Disable(); self.locate_spectrum_button.Disable()
        else:
            self.details.SetValue('No datasets match the current filters.' if self.datasets else
                                  'No recognised NMR acquisitions or deconParFile projects found below this folder.')
            self.open_button.Disable()
            self.locate_raw_button.Disable(); self.locate_spectrum_button.Disable()

    def _populate_navigation(self, datasets):
        from spinHubMain.core import build_navigation
        self.tree.DeleteAllItems()
        root_model = build_navigation(list(datasets), self.browser_location.root)
        root = self.tree.AddRoot(root_model.label)

        def add(parent, model):
            item = self.tree.AppendItem(parent, model.label)
            if model.dataset is not None:
                self.tree.SetItemData(item, model.dataset)
            for child in model.children:
                add(item, child)
            if model.kind in ('group', 'project'):
                self.tree.Expand(item)
            return item

        for node in root_model.children:
            add(root, node)
        self.tree.Expand(root)

    def OnTreeSelect(self, event):
        item = event.GetItem()
        dataset = self.tree.GetItemData(item) if item.IsOk() else None
        if dataset is None:
            # Container nodes select their first actionable child when possible.
            child, cookie = self.tree.GetFirstChild(item) if item.IsOk() else (None, None)
            if child and child.IsOk():
                dataset = self.tree.GetItemData(child)
        if dataset is None:
            return
        key = self._dataset_key(dataset)
        for idx, candidate in enumerate(self._row_datasets):
            if self._dataset_key(candidate) == key:
                self.lc.Select(idx)
                self.lc.Focus(idx)
                self.lc.EnsureVisible(idx)
                break

    def _update_location_controls(self):
        if hasattr(self, 'location_text'):
            self.location_text.ChangeValue(str(self.browser_location.root))
            self.back_button.Enable(self.browser_location.can_back)
            self.up_button.Enable(self.browser_location.can_up)

    def _change_location(self, path):
        try:
            self.browser_location.go(path)
        except (OSError, ValueError) as exc:
            wx.MessageBox(str(exc), 'Could not open location', wx.OK | wx.ICON_ERROR, self)
            self._update_location_controls()
            return
        self._update_location_controls()
        self.StartScan()

    def OnLocationEnter(self, event):
        self._change_location(self.location_text.GetValue())

    def OnBrowse(self, event):
        with wx.DirDialog(self, 'Choose SpinHub local area', defaultPath=str(self.browser_location.root),
                          style=wx.DD_DIR_MUST_EXIST) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                self._change_location(dlg.GetPath())

    def OnUp(self, event):
        self.browser_location.up()
        self._update_location_controls()
        self.StartScan()

    def OnBack(self, event):
        self.browser_location.back()
        self._update_location_controls()
        self.StartScan()

    def OnFilterChanged(self, event): self.ApplyFilters()
    def OnRefresh(self, event): self.StartScan()
    def OnCancelScan(self, event):
        self._scan_cancel.set(); self.scan_status.SetLabel('Cancelling scan...')
    def OnColumnClick(self, event):
        column = self.COLUMNS[event.GetColumn()][0]
        if column == self._sort_column: self._sort_descending = not self._sort_descending
        else: self._sort_column, self._sort_descending = column, False
        self.ApplyFilters()

    @staticmethod
    def _dataset_key(dataset):
        if dataset.project: return ('project', str(dataset.project.parameter_file.resolve(strict=False)))
        if dataset.acquisition: return ('acquisition', str(dataset.acquisition.path.resolve(strict=False)))
        return None

    def _selected_key(self):
        idx = self.lc.GetFirstSelected() if hasattr(self, 'lc') else -1
        if 0 <= idx < len(self._row_datasets): return self._dataset_key(self._row_datasets[idx])
        return None

    def selected_dataset(self):
        idx = self.lc.GetFirstSelected()
        return self._row_datasets[idx] if 0 <= idx < len(self._row_datasets) else None

    def Select(self, event):
        from spinHubMain.core import detail_lines, primary_action_label
        dataset = self.selected_dataset()
        if dataset is None: return
        self.details.SetValue('\n'.join(detail_lines(dataset)))
        self.open_button.SetLabel(primary_action_label(dataset))
        c = dataset.capabilities
        self.open_button.Enable(bool(c and (c.can_open_project or c.can_create_project)))
        project = dataset.project
        valid_project = bool(project and project.valid)
        self.locate_raw_button.Enable(valid_project and project.resources.raw_state.name == 'MISSING')
        self.locate_spectrum_button.Enable(valid_project and project.resources.spectrum_state.name == 'MISSING')
        from spinHubMain.core import resource_cards
        for index, card in enumerate(resource_cards(dataset)):
            self.stage_status[index].SetLabel(card.state)
            self.stage_path[index].SetLabel(card.path or 'No resource configured')
            self.stage_path[index].Wrap(max(240, self.inspector_panel.GetClientSize().width - 30))
            button = self.stage_buttons[index]
            button.SetLabel(card.action); button.Enable(card.enabled)
            # wx.Button is a wx.Control, not wx.ItemContainer; recent wxPython
            # versions therefore do not expose SetClientData/GetClientData here.
            # Store this small piece of view state directly on the button.
            button._spinhub_action_kind = card.action_kind
        self.Layout()

    def _show_repair_error(self, title, exc):
        wx.MessageBox(str(exc), title, wx.OK | wx.ICON_ERROR, self)

    def OnLocateRaw(self, event):
        dataset = self.selected_dataset()
        if not dataset or not dataset.project or not dataset.project.valid:
            return
        project = dataset.project
        start = str(project.resources.raw_path.parent) if project.resources.raw_path else str(project.parameter_file.parent)
        with wx.DirDialog(self, 'Locate raw NMR acquisition', defaultPath=start, style=wx.DD_DIR_MUST_EXIST) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            try:
                from spinHubMain.core.project_service import ProjectService
                ProjectService().relink_raw(
                    project.parameter_file, dlg.GetPath(),
                    expected_old_path=project.resources.raw_path)
            except Exception as exc:
                self._show_repair_error('Could not relink raw data', exc); return
        self.StartScan()

    def OnLocateSpectrum(self, event):
        dataset = self.selected_dataset()
        if not dataset or not dataset.project or not dataset.project.valid:
            return
        project = dataset.project
        old = project.resources.spectrum_path
        start_dir = str(old.parent) if old else str(project.parameter_file.parent)
        default_file = old.name if old else ''
        with wx.FileDialog(self, 'Locate main spectrum', defaultDir=start_dir, defaultFile=default_file,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            try:
                from spinHubMain.core.project_service import ProjectService
                ProjectService().relink_spectrum(
                    project.parameter_file, dlg.GetPath(), expected_old_path=old)
            except Exception as exc:
                self._show_repair_error('Could not relink spectrum', exc); return
        self.StartScan()

    def OnStageAction(self, event):
        dataset = self.selected_dataset()
        if dataset is None: return
        action = getattr(event.GetEventObject(), '_spinhub_action_kind', None)
        if action == 'locate_raw': return self.OnLocateRaw(event)
        if action == 'locate_spectrum': return self.OnLocateSpectrum(event)
        from spinHubMain.core.project_service import ProjectService
        service = ProjectService()
        if action == 'create' and dataset.acquisition is not None:
            state = service.create_for_acquisition(dataset.acquisition.path)
            service.open_project(state.parameter_file, workflow='prepare')
            self.StartScan(); return
        if dataset.project is not None and action in ('prepare', 'decon', 'inspect'):
            service.open_dataset(dataset, workflow=action)

    def OnOpen(self, event):
        dataset = self.selected_dataset()
        if dataset is None: return
        from spinHubMain.core.project_service import ProjectService
        service = ProjectService()
        if dataset.project is not None:
            service.open_dataset(dataset); return
        if dataset.acquisition is not None and dataset.capabilities.can_create_project:
            raw = dataset.acquisition.path
            state = service.create_for_acquisition(raw)
            service.open_project(state.parameter_file, workflow='prepare')
            self.StartScan()

    def OnButtonSort(self, event): event.Skip()
