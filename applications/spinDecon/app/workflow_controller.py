"""Workflow action routing.

This module owns the application-level routing that used to live in
``NotebookDemo.run_workflow_action``.  During migration it deliberately
operates on the notebook host so established GUI/scientific callbacks remain
unchanged while the routing responsibility has a stable home.
"""
import os

import wx

from ..domain.analysis_mode import AnalysisMode


class WorkflowController:
    """Route workflow stage actions through the current application host."""

    def __init__(self, context=None, legacy_notebook=None):
        self.context = context
        self.legacy_notebook = legacy_notebook

    def run(self, action_key):
        host = self.legacy_notebook
        if host is None:
            raise RuntimeError("No workflow action host has been configured")
        return _run_workflow_action(host, action_key)


def _run_workflow_action(self, action_key):
    """Route Workflow actions through established GUI/scientific handlers.

    Modeless child windows launched from Workflow no longer force the
    notebook onto NMR.  Actions that genuinely need a notebook workspace
    (for example 2D Slices) still select that workspace.
    """
    tab = getattr(self, 'nmr_workspace', None)
    if tab is None:
        # Compatibility for third-party notebook hosts predating Stage 126.
        tab = getattr(self, 'tabOne', None)
    if tab is None:
        return False

    if action_key in ('extract_intensities', 'review_series', 'analyse_series'):
        if action_key == 'review_series':
            pass
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            return False
        debug = getattr(tab, '_workflow_debug', print)
        if not mode.has_pseudo_axis:
            debug('WORKFLOW ROUTE ABORT: project AnalysisMode says has_pseudo_axis=False')
            return False

        # A restored 2D+pseudo project may reach Workflow before the NMR
        # spectrum has been read into the live tab.  Do not ask
        # _is_pseudo3d_topology() to classify an unloaded dataset: that creates
        # a circular route where the guard prevents entry to the branch
        # that would load the spectrum.  Load the configured spectrum
        # first, then use its labels/shape as the scientific guard.
        if mode.spectral_dimensions == 2:
            infile = tab.infileBox.GetValue().strip() if hasattr(tab, 'infileBox') else ''
            spectrum = tab._resolve_input_path(infile) if infile else ''
            if not spectrum or not os.path.isfile(spectrum):
                wx.MessageBox(
                    'The processed spectrum file cannot be found. Process or select the spectrum first.',
                    'Workflow', wx.OK | wx.ICON_WARNING)
                return False
            current = getattr(tab, 'spectrumfile', None)
            if current != spectrum or getattr(tab, 'data', None) is None or not getattr(tab, 'READ', 0):
                debug('WORKFLOW ROUTE: loading spectrum before pseudo3D topology check: %r' % spectrum)
                tab.OnButtonRead(None)
            if not (getattr(tab, 'READ', 0) and getattr(tab, 'data', None) is not None):
                debug('WORKFLOW ROUTE ABORT: spectrum did not load before pseudo3D topology check')
                return False

        pseudo3d_topology = tab._is_pseudo3d_topology() if hasattr(tab, '_is_pseudo3d_topology') else None
        debug('WORKFLOW ROUTE action=%r mode.spectral_dimensions=%r mode.has_pseudo_axis=%r tab.dim=%r pseudo_axis=%r pseudo3d_topology=%r' %
              (action_key, mode.spectral_dimensions, mode.has_pseudo_axis, getattr(tab, 'dim', None),
               tab.state.pseudo_axis if getattr(tab, 'state', None) is not None else None, pseudo3d_topology))
        debug('WORKFLOW ROUTE metadata: labels=%r data_shape=%r specsize=%r' %
              (getattr(tab, 'labb', None), getattr(getattr(tab, 'data', None), 'shape', None), getattr(tab, 'specsize', None)))
        if mode.spectral_dimensions == 1:
            if action_key == 'extract_intensities':
                # Pseudo2D restrained extraction is the established Recon
                # path with Fit, Use 2D peaklist and phase fitting enabled.
                # Here the single authoritative spectral list is Full 1D.
                full_value = tab.fullPeakBox.GetValue().strip() if hasattr(tab, 'fullPeakBox') else ''
                full_path = tab._resolve_spec_file(full_value) if full_value else ''
                if not full_path or not os.path.isfile(full_path):
                    wx.MessageBox(
                        'A valid Full 1D peak list is required before extracting pseudo-axis intensities.',
                        'Workflow', wx.OK | wx.ICON_WARNING)
                    return False
                store = getattr(self, 'data_store', None)
                if store is not None:
                    store.invalidate_pseudo_series_review()
                tab.cb_decon3d.SetValue(True)
                tab.cb_decback.SetValue(True)
                tab.cb_fitphases.SetValue(True)
                save = getattr(tab, 'OnButtonSave', None)
                if callable(save):
                    save(True)
                tab.OnButtonRecon(None)
                return True
            if not self.PageExists('Pseudo2D'):
                self.AddTabPseudo2D(True, tab)
            self.select_page('Pseudo2D')
            if action_key == 'review_series':
                # Review belongs to the Pseudo2D workspace and opens its
                # existing modeless fitting palette.  Merely revisiting
                # the accepted fitting results must not revoke persisted
                # review evidence: that evidence is invalidated when the
                # underlying extraction/reconstruction is rerun above.
                self.notify_analysis_changed()
                pseudo = self.get_page_by_title('Pseudo2D')
                show = getattr(pseudo, 'show_fitting_window', None)
                pass
                if callable(show):
                    wx.CallAfter(show)
                    pass
            return True
        # Canonical 2D+pseudo is always two spectral dimensions plus one
        # real physical axis. Legacy physical counts are normalized on load.
        if mode.spectral_dimensions == 2 and (pseudo3d_topology is None or pseudo3d_topology):
            debug('WORKFLOW ROUTE: entering canonical pseudo3D branch')
            # Physical pseudo3D intensities depend on the *loaded* reference
            # 2D peak list and on one .out/.dat pair per reference peak.
            # Keep the workflow on the same handlers as the main NMR tab:
            # Load spectrum -> Load reference list -> inspect fit/ -> Recon
            # (Fit + Use 2D peak list) only when outputs are missing.
            if action_key == 'extract_intensities':
                # If Extract is visible/explicitly revisited, the user's
                # action means run extraction now, never merely accept old
                # fit files. A new extraction also requires a fresh review.
                store = getattr(self, 'data_store', None)
                if store is not None:
                    store.invalidate_pseudo_series_review()
                save = getattr(tab, 'OnButtonSave', None)
                if callable(save):
                    save(True)
                debug = getattr(tab, '_workflow_debug', print)
                debug('WORKFLOW extract_intensities ENTER: referenceBox=%r refs_in_store=%d' % (tab.referencePeakBox.GetValue() if hasattr(tab, 'referencePeakBox') else None, len(tab.get_reference_peaks() or [])))
                # Peak-shape/Fit Peaks normally loaded the main spectrum
                # already. Re-check here so a revisited workflow cannot run
                # Protocol3P against an unloaded/stale main dataset.
                infile = tab.infileBox.GetValue().strip() if hasattr(tab, 'infileBox') else ''
                spectrum = tab._resolve_input_path(infile) if infile else ''
                if not spectrum or not os.path.isfile(spectrum):
                    wx.MessageBox(
                        'The processed spectrum file cannot be found. Process or select the spectrum first.',
                        'Workflow', wx.OK | wx.ICON_WARNING)
                    return False
                current = getattr(tab, 'spectrumfile', None)
                if current != spectrum or getattr(tab, 'data', None) is None or not getattr(tab, 'READ', 0):
                    tab.OnButtonRead(None)
                if not (getattr(tab, 'READ', 0) and getattr(tab, 'data', None) is not None):
                    return False

                # The reference list has been decided by this stage, but it
                # may only exist as the filename in the main tab. Invoke the
                # same Load callback the user would press before examining
                # fit outputs.
                ref_value = tab.referencePeakBox.GetValue().strip() if hasattr(tab, 'referencePeakBox') else ''
                ref_path = tab._resolve_spec_file(ref_value) if ref_value else ''
                if not ref_path or not os.path.isfile(ref_path):
                    wx.MessageBox(
                        'Decide the reference peak list before extracting pseudo-axis intensities.',
                        'Workflow', wx.OK | wx.ICON_WARNING)
                    return False
                load_reference = getattr(tab, 'ensure_reference_peak_list_loaded', None)
                if callable(load_reference):
                    load_reference()
                else:
                    tab.OnButtonReadPeak(None)
                refs = list(tab.get_reference_peaks() or [])
                debug('WORKFLOW after reference-load request: refs_in_store=%d' % len(refs))
                if not refs:
                    wx.MessageBox(
                        'The selected reference peak list could not be loaded.',
                        'Workflow', wx.OK | wx.ICON_WARNING)
                    return False

                ensure = getattr(tab, 'ensure_pseudo3d_fit_results', None)
                if callable(ensure):
                    # Extraction owns only the existence check / Protocol3P
                    # launch.  Whether files already existed or Recon was
                    # started, the fitting UI belongs to the *next* workflow
                    # step (Inspect fitting results).
                    ensured = ensure(force_recompute=True)
                    debug('WORKFLOW ensure_pseudo3d_fit_results(force_recompute=True) returned %r' % (ensured,))
                return True
            # Pseudo3D/Fitting builds its rows from the shared Reference
            # 2D peak objects.  Ensure the configured list has actually
            # been read before constructing/opening that workspace.
            debug = getattr(tab, '_workflow_debug', print)
            debug('WORKFLOW %s ENTER: refs_in_store=%d referenceBox=%r' % (action_key, len(tab.get_reference_peaks() or []), tab.referencePeakBox.GetValue() if hasattr(tab, 'referencePeakBox') else None))
            load_reference = getattr(tab, 'ensure_reference_peak_list_loaded', None)
            if callable(load_reference) and not load_reference():
                wx.MessageBox(
                    'The reference 2D peak list must be loaded before opening fitting results.',
                    'Workflow', wx.OK | wx.ICON_WARNING)
                return False
            debug('WORKFLOW reference loaded for %s: refs_in_store=%d' % (action_key, len(tab.get_reference_peaks() or [])))
            if not self.PageExists('Fitting'):
                debug('WORKFLOW creating Fitting/Pseudo3D tab')
                self.AddTabPseudo3D(True, tab)
            self.select_page('Fitting')
            pseudo = self.get_page_by_title('Fitting')
            if action_key == 'review_series':
                # Re-entering Review means the currently displayed result
                # must be explicitly accepted again. Persist that
                # invalidation immediately so save/load mirrors the GUI.
                store = getattr(self, 'data_store', None)
                if store is not None:
                    store.invalidate_pseudo_series_review()
                save = getattr(tab, 'OnButtonSave', None)
                if callable(save):
                    save(True)
                self.notify_analysis_changed()
                show = getattr(pseudo, 'show_fitting_window', None)
                if callable(show):
                    # Let notebook selection/layout complete before raising
                    # the modeless fitting palette.
                    debug('WORKFLOW scheduling show_fitting_window via wx.CallAfter')
                    wx.CallAfter(show)
            elif action_key == 'analyse_series':
                # First visit: show the Pseudo3D Analysis palette so the
                # user can choose/confirm an analysis type.  Once a choice
                # has been persisted in the system file, Workflow becomes
                # a direct "Show analysis" action on subsequent runs.
                saved = str(getattr(tab, 'downstream_analysis', '') or '')
                available = getattr(pseudo, 'available_downstream_analyses', lambda: [])()
                if saved and saved in available:
                    open_saved = getattr(pseudo, 'open_saved_analysis', None)
                    if callable(open_saved):
                        wx.CallAfter(open_saved)
                else:
                    selector = getattr(pseudo, 'show_analysis_selector', None)
                    if callable(selector):
                        wx.CallAfter(selector)
            return True
        debug('WORKFLOW ROUTE ABORT: no pseudo workflow branch matched mode.spectral_dimensions=%r pseudo3d_topology=%r' %
              (mode.spectral_dimensions, tab._is_pseudo3d_topology() if hasattr(tab, '_is_pseudo3d_topology') else None))
        return False

    def ensure_loaded():
        infile = tab.infileBox.GetValue().strip() if hasattr(tab, 'infileBox') else ''
        spectrum = tab._resolve_input_path(infile) if infile else ''
        if not spectrum or not os.path.isfile(spectrum):
            wx.MessageBox('The processed spectrum file cannot be found. Process or select the spectrum first.', 'Workflow', wx.OK | wx.ICON_WARNING)
            return False
        current = getattr(tab, 'spectrumfile', None)
        if current != spectrum or getattr(tab, 'data', None) is None or not getattr(tab, 'READ', 0):
            tab.OnButtonRead(None)
        return bool(getattr(tab, 'READ', 0) and getattr(tab, 'data', None) is not None)

    if action_key == 'spectrum':
        # ProcMan is modeless.  Preserve Workflow as the selected notebook
        # page; processFrame already writes the resulting filename back to
        # infileBox/ProjectState when processing completes.
        tab.OnButtonProcess(None)
        return True
    if action_key == 'peak_shape':
        if not ensure_loaded(): return False
        tab.OnButtonPeakFit(None)
        return True
    if action_key == 'reference_peaks':
        if not ensure_loaded(): return False
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            mode = None
        if mode is not None and mode.has_pseudo_axis and mode.spectral_dimensions == 1:
            # Pseudo2D establishes its reference frequencies from the
            # projection workflow; Full 1D is the authoritative list.
            if not self.PageExists('Projections'):
                self.AddTabTwo(True, tab)
            return self.select_page('Projections')
        # PeakFrame is the established projection editor.  If a saved
        # reference list is configured but not yet in DataStore, load it
        # first so PeakFrame opens on the authoritative list.
        # Completing this stage must materialise the configured reference
        # list, even when Workflow already shows a tick because the file
        # exists.  This keeps all later stages dependent on the same loaded
        # authoritative list rather than merely on path evidence.
        load_reference = getattr(tab, 'ensure_workflow_reference_stage_loaded', None)
        if not callable(load_reference):
            load_reference = getattr(tab, 'ensure_reference_peak_list_loaded', None)
        if callable(load_reference):
            try: load_reference()
            except Exception: pass
        elif not tab.get_reference_peaks():
            ref = tab.referencePeakBox.GetValue().strip() if hasattr(tab, 'referencePeakBox') else ''
            if ref:
                try: tab.OnButtonReadPeak(None)
                except Exception: pass
        tab.OnButtonPeaky(None)
        return True
    if action_key == 'peak_pick':
        if not ensure_loaded(): return False
        # Physical 3D/4D picking is constrained by the reference list.
        # A file existing on disk is not enough: the deconvolution code
        # consumes the in-memory peak/index structures, so materialise it
        # before starting the job.
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            mode = None
        if mode is not None and not mode.has_pseudo_axis and mode.spectral_dimensions >= 3:
            load_reference = getattr(tab, 'ensure_reference_peak_list_loaded', None)
            if not callable(load_reference) or not load_reference():
                wx.MessageBox('The reference peak list must be loaded before picking peaks.',
                              'Workflow', wx.OK | wx.ICON_WARNING)
                return False
            # Workflow peak picking for a physical 3D/4D spectrum is
            # reference-constrained.  Keep the main-tab option in sync
            # with that workflow invariant rather than relying on the
            # user's previous checkbox state.
            use_2d = getattr(tab, 'cb_decon3d', None)
            if use_2d is not None:
                use_2d.SetValue(True)
        # Starting a new pick invalidates acceptance of the previous list
        # immediately; successful completion will also clear stale-shape
        # evidence and persist the newly generated list.
        self.clear_workflow_picked_peaks_checked()
        tab.OnButtonDecon(None)
        return True
    if action_key == 'review_peaks':
        # Review is a composite read operation.  Load its four project
        # products in dependency order before constructing either viewer:
        # spectrum -> reference peaks -> Full nD peaks -> deconvolution.
        # This is particularly important after a cold system-file load.
        prepare_review = getattr(tab, 'ensure_workflow_review_inputs_loaded', None)
        if callable(prepare_review):
            ok, message = prepare_review()
            if not ok:
                wx.MessageBox(message, 'Workflow', wx.OK | wx.ICON_WARNING)
                return False
        elif not ensure_loaded():
            return False
        tab.OnButtonFullPeakList(None)
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            mode = None
        if mode is not None and mode.spectral_dimensions == 2 and not mode.has_pseudo_axis:
            # Physical 2D review uses the actual 2D plane, not the 3D Slice2D viewer.
            if not self.PageExists('Projections'):
                self.AddTabTwo(True, tab)
            tab.OnButtonPeaky(None)
            return self.select_page('Projections')
        if not self.PageExists('2D Slices'):
            self.AddTabFour(True, tab)
        return self.select_page('2D Slices')
    if action_key == 'fit_spectrum':
        if not ensure_loaded(): return False
        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            mode = None
        if mode is not None and mode.spectral_dimensions == 2 and not mode.has_pseudo_axis:
            # Match the pseudo3D restrained-fitting workflow exactly: the
            # checked Full 2D list is the restraint, with both Fit and
            # Use 2D peaklist selected before Recon is launched.
            if not self.data_store.analysis.get('picked_peaks_reviewed'):
                wx.MessageBox('Review the picked peaks and mark them as checked before fitting.',
                              'Workflow', wx.OK | wx.ICON_WARNING)
                return False
            load_full = getattr(tab, 'ensure_full_peak_list_loaded', None)
            if callable(load_full) and not load_full():
                wx.MessageBox('The Full 2D peak list could not be loaded.',
                              'Workflow', wx.OK | wx.ICON_WARNING)
                return False
            self.data_store.invalidate_fitting_review()
            tab.cb_decon3d.SetValue(True)
            tab.cb_decback.SetValue(True)
            save = getattr(tab, 'OnButtonSave', None)
            if callable(save):
                save(True)
            tab.OnButtonRecon(None)
            return True
        tab.OnButtonPeakFit(None)
        return True
    if action_key == 'review_fitting':
        if not ensure_loaded(): return False
        if not self.data_store.analysis.get('fitting_results_ready'):
            wx.MessageBox('Fit the checked 2D peaks before reviewing fitting results.',
                          'Workflow', wx.OK | wx.ICON_WARNING)
            return False
        load_full = getattr(tab, 'ensure_full_peak_list_loaded', None)
        if callable(load_full) and not load_full():
            return False
        if not self.PageExists('Fitting'):
            self.AddTabPseudo3D(True, tab)
        self.select_page('Fitting')
        pseudo = self.get_page_by_title('Fitting')
        show = getattr(pseudo, 'show_fitting_window', None)
        if callable(show):
            wx.CallAfter(show)
        return True
    return False
