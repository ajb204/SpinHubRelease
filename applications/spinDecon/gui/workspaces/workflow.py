"""Authoritative workflow overview and dataset-topology selector.

Milestone 5 adds evidence-based status and next-step guidance for 1-4D spectral
workflows. Dataset changes and actions still delegate to legacy NMR callbacks;
scientific operations are not reimplemented here.
"""
import wx
import traceback


def _p2d_return_debug(message):
    return

from spinDecon.domain.analysis_mode import AnalysisMode
from spinDecon.workflow.model import StageRequirement, build_workflow_plan
from spinDecon.workflow.status import StageStatus, evaluate_workflow, recommended_action


class WorkflowOverviewPanel(wx.ScrolledWindow):
    """Guided workflow page layered over the existing NMR handlers."""

    def __init__(self, notebook, parent, state=None):
        super().__init__(parent)
        self.notebook = notebook
        self.state = state
        self._content = None
        self._updating_controls = False
        self._enter_action_key = None
        self.SetScrollRate(10, 10)
        # Handle Enter locally instead of installing a top-level default button.
        # wx default buttons remain active even while their notebook page is
        # hidden, which allowed Enter on NMR/Threshold to launch Workflow.
        self.Bind(wx.EVT_CHAR_HOOK, self._on_char_hook)
        self._build_shell()
        self.refresh()

    def _build_shell(self):
        root = wx.BoxSizer(wx.VERTICAL)

        title = wx.StaticText(self, label="Analysis workflow")
        title_font = title.GetFont()
        title_font.PointSize += 4
        title.SetFont(title_font.Bold())
        root.Add(title, 0, wx.ALL, 12)

        intro = wx.StaticText(
            self,
            label=(
                "This page summarises the recommended scientific workflow for the current "
                "dataset. Dataset topology can be set here; existing NMR, UniDec and specialist "
                "tabs continue to perform all analysis operations."
            ),
        )
        intro.Wrap(800)
        root.Add(intro, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)

        dataset_box = wx.StaticBoxSizer(wx.VERTICAL, self, "1. Define dataset")
        dataset_parent = dataset_box.GetStaticBox()
        dataset_help = wx.StaticText(
            dataset_parent,
            label=("Choose the dataset type here. This single decision defines the spectral dimensionality "
                   "and whether a real pseudo-axis is present, using the existing NMR controls underneath."),
        )
        dataset_help.Wrap(800)
        dataset_box.Add(dataset_help, 0, wx.ALL, 6)

        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(dataset_parent, label="Spectral dimensions:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.spectral_dim = wx.Choice(dataset_parent, choices=["1", "2", "3", "4"])
        row.Add(self.spectral_dim, 0, wx.RIGHT, 18)
        self.pseudo_axis = wx.CheckBox(dataset_parent, label="Contains a real pseudo-axis")
        row.Add(self.pseudo_axis, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 18)
        row.Add(wx.StaticText(dataset_parent, label="Physical dimensions:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.physical_dim = wx.StaticText(dataset_parent, label="-")
        row.Add(self.physical_dim, 0, wx.ALIGN_CENTER_VERTICAL)
        dataset_box.Add(row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)
        root.Add(dataset_box, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)

        self.spectral_dim.Bind(wx.EVT_CHOICE, self._on_dataset_type_changed)
        self.pseudo_axis.Bind(wx.EVT_CHECKBOX, self._on_dataset_type_changed)

        self._content = wx.BoxSizer(wx.VERTICAL)
        root.Add(self._content, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)

        note = wx.StaticText(
            self,
            label="Workflow actions use the existing NMR handlers. Use the detailed tabs whenever advanced settings need to be changed.",
        )
        note.Wrap(800)
        root.Add(note, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)

        self.SetSizer(root)


    def _on_char_hook(self, event):
        """Run the recommended action on Enter only while Workflow is visible."""
        keycode = event.GetKeyCode()
        if keycode not in (wx.WXK_RETURN, wx.WXK_NUMPAD_ENTER):
            event.Skip()
            return

        try:
            selected = self.notebook.GetSelection()
            workflow_is_active = (
                selected != wx.NOT_FOUND
                and self.notebook.GetPage(selected) is self
            )
        except Exception:
            workflow_is_active = False

        if not workflow_is_active or not self._enter_action_key:
            event.Skip()
            return

        self._on_workflow_action(event, self._enter_action_key)

    def _on_dataset_type_changed(self, event):
        """Delegate topology changes to NotebookDemo's compatibility boundary."""
        if self._updating_controls:
            return
        selection = self.spectral_dim.GetSelection()
        if selection == wx.NOT_FOUND:
            return
        spectral_dimensions = selection + 1
        pseudo = bool(self.pseudo_axis.GetValue())
        if pseudo and spectral_dimensions == 4:
            # Mirror the existing NMR restriction before invoking it, avoiding
            # a transient misleading workflow display.
            self._updating_controls = True
            self.pseudo_axis.SetValue(False)
            self._updating_controls = False
            pseudo = False
        apply_type = getattr(self.notebook, "apply_workflow_dataset_type", None)
        if callable(apply_type):
            apply_type(spectral_dimensions, pseudo)
        self.refresh()

    def _sync_dataset_controls(self):
        """Mirror canonical ProjectState topology without firing wx events."""
        tab = getattr(self.notebook, "tabOne", None)
        if tab is None:
            return
        try:
            state = tab.state
            dim = int(state.spectral_dimensions)
            pseudo = bool(state.pseudo_axis)
        except Exception:
            return
        self._updating_controls = True
        try:
            self.spectral_dim.SetSelection(dim - 1 if 1 <= dim <= 4 else wx.NOT_FOUND)
            self.pseudo_axis.SetValue(pseudo)
            self.physical_dim.SetLabel(str(dim + int(pseudo)) if 1 <= dim <= 4 else "-")
        finally:
            self._updating_controls = False

    def _on_workflow_action(self, event, action_key):
        """Delegate an action to the notebook compatibility boundary."""
        context = getattr(self.notebook, "app_context", None)
        controller = getattr(context, "workflow", None) if context is not None else None
        run_action = getattr(controller, "run", None)
        if not callable(run_action):
            # Compatibility fallback for external/older notebook hosts.
            run_action = getattr(self.notebook, "run_workflow_action", None)
        if callable(run_action):
            try:
                run_action(action_key)
            except Exception as exc:
                wx.MessageBox(
                    "The existing NMR action could not be opened or started.\n\n%s" % exc,
                    "Workflow",
                    wx.OK | wx.ICON_ERROR,
                )

    def _add_section_label(self, text):
        label = wx.StaticText(self, label=text)
        font = label.GetFont()
        font.PointSize += 1
        label.SetFont(font.Bold())
        self._content.Add(label, 0, wx.TOP | wx.BOTTOM, 6)
        return label

    def _add_stage_line(self, stage, state, marker, backtrack=False):
        """Add a compact stage row with a clear, coloured status symbol."""
        row = wx.BoxSizer(wx.HORIZONTAL)
        icon = wx.StaticText(self, label=marker)
        icon_font = icon.GetFont()
        icon_font.PointSize += 2
        icon.SetFont(icon_font.Bold())
        if state and state.status is StageStatus.COMPLETE:
            icon.SetForegroundColour(wx.Colour(34, 139, 74))
        elif marker in ("\u2715", "\u2717", "\u2718"):
            icon.SetForegroundColour(wx.Colour(196, 55, 55))

        line = wx.StaticText(self, label=stage.title)
        tooltip = state.detail if state and state.detail else ""
        if tooltip:
            icon.SetToolTip(tooltip)
            line.SetToolTip(tooltip)
        if backtrack:
            revisit_tip = tooltip + ("\n\n" if tooltip else "") + "Click to revisit this completed workflow step."
            for control in (icon, line):
                control.SetCursor(wx.Cursor(wx.CURSOR_HAND))
                control.SetToolTip(revisit_tip)
                control.Bind(
                    wx.EVT_LEFT_UP,
                    lambda event, current_stage=stage: self._show_backtrack_menu(event, current_stage),
                )
        row.Add(icon, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 7)
        row.Add(line, 0, wx.ALIGN_CENTER_VERTICAL)
        self._content.Add(row, 0, wx.LEFT | wx.BOTTOM, 4)

    def _show_backtrack_menu(self, event, stage):
        """Offer a non-destructive return to a completed scientific step.

        Backtracking deliberately does not delete downstream results. The user
        is returned to the existing scientific tool and can revise/re-run the
        step; the normal evidence callbacks then keep Workflow synchronised.
        """
        menu = wx.Menu()
        item = menu.Append(wx.ID_ANY, f"Return to {stage.title}...")
        self.Bind(
            wx.EVT_MENU,
            lambda evt, key=stage.key: self._on_backtrack_action(evt, key),
            item,
        )
        menu.AppendSeparator()
        info = menu.Append(wx.ID_ANY, "Existing later results will be kept")
        info.Enable(False)
        self.PopupMenu(menu)
        menu.Destroy()

    def _on_backtrack_action(self, event, action_key):
        """Return to a completed stage, with detailed diagnostics for pseudo2D review."""
        pass
        pass
        # Do not mutate/rebuild Workflow while wx is still dispatching the
        # popup-menu event.  On macOS that can invalidate native menu/window
        # objects underneath the current callback.  Defer the operation until
        # the menu event has completely unwound; keep tracing each boundary.
        wx.CallAfter(self._perform_backtrack_action, action_key)
        pass

    def _perform_backtrack_action(self, action_key):
        pass
        try:
            if action_key == "review_series":
                clear = getattr(self.notebook, "clear_workflow_series_inspected", None)
                pass
                if callable(clear):
                    clear_result = clear()
                    pass
            elif action_key == "review_fitting":
                clear = getattr(self.notebook, "clear_workflow_fitting_inspected", None)
                if callable(clear):
                    clear()
            elif action_key == "review_peaks":
                clear = getattr(self.notebook, "clear_workflow_picked_peaks_checked", None)
                if callable(clear):
                    clear()
            pass
            self._on_workflow_action(None, action_key)
            pass
        except Exception as exc:
            pass
            wx.MessageBox(
                "The workflow step could not be reopened.\n\n%s" % exc,
                "Workflow", wx.OK | wx.ICON_ERROR,
            )
        finally:
            pass

    def _add_next_step_card(self, stage, state, action_label):
        box = wx.StaticBoxSizer(wx.VERTICAL, self, "Next step")
        box_parent = box.GetStaticBox()

        heading = wx.StaticText(box_parent, label=stage.title)
        font = heading.GetFont()
        font.PointSize += 3
        heading.SetFont(font.Bold())
        box.Add(heading, 0, wx.ALL, 8)

        desc = wx.StaticText(box_parent, label=stage.description)
        desc.Wrap(760)
        box.Add(desc, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        if state and state.detail:
            detail = wx.StaticText(box_parent, label=state.detail)
            detail.Wrap(760)
            box.Add(detail, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        button = wx.Button(box_parent, label=action_label)
        self._enter_action_key = stage.key
        button.Bind(wx.EVT_BUTTON, lambda event, key=stage.key: self._on_workflow_action(event, key))
        box.Add(button, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)

        # Once a downstream pseudo3D analysis has been confirmed, the primary
        # action is deliberately just "Show analysis".  Keep an explicit
        # secondary route back to the selector so the user can revise the
        # persisted choice without clearing any scientific results.
        if stage.key == "analyse_series":
            tab_one = getattr(self.notebook, "tabOne", None)
            if str(getattr(tab_one, "downstream_analysis", "") or "").strip():
                change = wx.Button(box_parent, label="Change analysis type...")
                change.Bind(wx.EVT_BUTTON, self._on_change_analysis_type)
                box.Add(change, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)

        # Reviewing pseudo3D intensities is intentionally a two-step action:
        # first inspect the existing Fitting results, then explicitly accept
        # them. Opening the results alone is not completion evidence.
        if stage.key == "review_series":
            inspected = wx.Button(box_parent, label="Mark as inspected")
            inspected.Bind(wx.EVT_BUTTON, self._on_series_inspected)
            box.Add(inspected, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        elif stage.key == "review_fitting":
            inspected = wx.Button(box_parent, label="Mark as inspected")
            inspected.Bind(wx.EVT_BUTTON, self._on_fitting_inspected)
            box.Add(inspected, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        elif stage.key == "review_peaks":
            checked = wx.Button(box_parent, label="Mark as checked")
            checked.Bind(wx.EVT_BUTTON, self._on_picked_peaks_checked)
            box.Add(checked, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        self._content.Add(box, 0, wx.EXPAND | wx.BOTTOM, 14)


    def _add_terminal_analysis_card(self, stage, state, action_label):
        """Keep the terminal pseudo-series analysis visible after completion."""
        box = wx.StaticBoxSizer(wx.VERTICAL, self, "Workflow complete")
        box_parent = box.GetStaticBox()
        heading_row = wx.BoxSizer(wx.HORIZONTAL)
        tick = wx.StaticText(box_parent, label="\u2713")
        tick_font = tick.GetFont()
        tick_font.PointSize += 4
        tick.SetFont(tick_font.Bold())
        tick.SetForegroundColour(wx.Colour(34, 139, 74))
        heading = wx.StaticText(box_parent, label=stage.title)
        font = heading.GetFont()
        font.PointSize += 3
        heading.SetFont(font.Bold())
        heading_row.Add(tick, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 8)
        heading_row.Add(heading, 0, wx.ALIGN_CENTER_VERTICAL)
        box.Add(heading_row, 0, wx.ALL, 8)
        if state and state.detail:
            detail = wx.StaticText(box_parent, label=state.detail)
            detail.Wrap(760)
            box.Add(detail, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        button = wx.Button(box_parent, label=action_label)
        button.Bind(wx.EVT_BUTTON, lambda event, key=stage.key: self._on_workflow_action(event, key))
        box.Add(button, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)
        change = wx.Button(box_parent, label="Change analysis type...")
        change.Bind(wx.EVT_BUTTON, self._on_change_analysis_type)
        box.Add(change, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        self._content.Add(box, 0, wx.EXPAND | wx.BOTTOM, 14)


    def _on_change_analysis_type(self, event):
        """Open the Pseudo3D Analysis palette even when a choice is saved."""
        get_page = getattr(self.notebook, "get_page_by_title", None)
        pseudo = get_page("Fitting") if callable(get_page) else getattr(self.notebook, "tabPseudo", None)
        if pseudo is None:
            # Ensure the specialist workspace exists, but do not use the normal
            # analyse_series route because a saved choice would launch it.
            tab_one = getattr(self.notebook, "tabOne", None)
            add = getattr(self.notebook, "AddTabPseudo3D", None)
            if callable(add) and tab_one is not None:
                add(True, tab_one)
                pseudo = get_page("Fitting") if callable(get_page) else getattr(self.notebook, "tabPseudo", None)
        selector = getattr(pseudo, "show_analysis_selector", None) if pseudo is not None else None
        if callable(selector):
            wx.CallAfter(selector)

    def _on_picked_peaks_checked(self, event):
        """Explicitly accept the displayed full-dimensional peak list."""
        mark = getattr(self.notebook, "mark_workflow_picked_peaks_checked", None)
        if not callable(mark):
            return
        try:
            if mark():
                self.refresh()
        except Exception as exc:
            wx.MessageBox(
                "The picked peaks could not be marked as checked.\n\n%s" % exc,
                "Workflow", wx.OK | wx.ICON_ERROR,
            )

    def _on_fitting_inspected(self, event):
        """Explicitly accept the displayed physical-2D fitting results."""
        mark = getattr(self.notebook, "mark_workflow_fitting_inspected", None)
        if not callable(mark):
            return
        try:
            if mark():
                self.refresh()
        except Exception as exc:
            wx.MessageBox(
                "The fitting results could not be marked as inspected.\n\n%s" % exc,
                "Workflow", wx.OK | wx.ICON_ERROR,
            )

    def _on_series_inspected(self, event):
        """Explicitly accept the displayed pseudo-dimensional fit results."""
        mark = getattr(self.notebook, "mark_workflow_series_inspected", None)
        if not callable(mark):
            return
        try:
            if mark():
                self.refresh()
        except Exception as exc:
            wx.MessageBox(
                "The intensity series could not be marked as inspected.\n\n%s" % exc,
                "Workflow", wx.OK | wx.ICON_ERROR,
            )

    def _on_summarise_project(self, event):
        """Generate the NMR project-summary PDF and return to Workflow."""
        tab_one = getattr(self.notebook, "tabOne", None)
        generate = getattr(tab_one, "OnButtonSummariseProject", None)
        if not callable(generate):
            return
        try:
            generate(event)
        finally:
            # Report generation/viewer creation may manipulate focus; keep the
            # main GUI anchored on the Workflow page when it completes.
            select = getattr(self.notebook, "select_page", None)
            if callable(select):
                wx.CallAfter(select, "Workflow")

    def refresh(self):
        """Render one dominant recommendation over the evidence-based workflow."""
        if self._content is None:
            return
        self._sync_dataset_controls()
        self._enter_action_key = None
        self._content.Clear(delete_windows=True)

        try:
            mode = AnalysisMode.from_project_state(self.state)
        except (TypeError, ValueError):
            msg = wx.StaticText(
                self,
                label=(
                    "Dataset dimensionality has not been defined yet.\n\n"
                    "Choose the spectral dimensionality above. Once it is defined, this page "
                    "will show the recommended next scientific step."
                ),
            )
            msg.Wrap(800)
            self._content.Add(msg, 0, wx.ALL, 8)
            self.Layout()
            return

        plan = build_workflow_plan(mode)
        store = getattr(self.notebook, "data_store", None)
        stage_states = evaluate_workflow(plan, self.state, store, self.notebook)
        states_by_key = {item.key: item for item in stage_states}
        recommended = recommended_action(plan, stage_states)
        terminal_analysis_complete = bool(
            mode.has_pseudo_axis
            and mode.spectral_dimensions == 2
            and states_by_key.get("analyse_series") is not None
            and states_by_key["analyse_series"].status is StageStatus.COMPLETE
            and str(getattr(getattr(self.notebook, "tabOne", None), "downstream_analysis", "") or "").strip()
        )

        workflow_name = "Intensity-series workflow" if mode.has_pseudo_axis else "Peak-list workflow"
        dimensions = f"{mode.spectral_dimensions}D spectrum"
        if mode.has_pseudo_axis:
            dimensions += " + pseudo-axis"
        summary = wx.StaticText(self, label=f"{dimensions}  |  {workflow_name}\n{plan.objective}")
        summary.SetFont(summary.GetFont().Bold())
        summary.Wrap(800)
        self._content.Add(summary, 0, wx.BOTTOM, 14)

        action_labels = {
            "spectrum": "Prepare spectrum...",
            "peak_shape": "Determine peak shape...",
            "reference_peaks": "Establish reference peaks...",
            "peak_pick": "Pick peaks...",
            "fit_spectrum": "Fit spectrum...",
            "review_fitting": "Inspect fitting results...",
            "review_peaks": "Review picked peaks...",
            "extract_intensities": "Extract intensities...",
            "review_series": "Inspect fitting results...",
            "analyse_series": "Analyse intensity series...",
        }
        tab_one = getattr(self.notebook, "tabOne", None)
        if str(getattr(tab_one, "downstream_analysis", "") or "").strip():
            action_labels["analyse_series"] = "Show analysis..."

        completed = []
        later = []
        for stage in plan.stages:
            state = states_by_key.get(stage.key)
            if state is None:
                continue
            if state.status is StageStatus.COMPLETE:
                if not (terminal_analysis_complete and stage.key == "analyse_series"):
                    completed.append((stage, state))
            elif stage.key != recommended:
                later.append((stage, state))

        if completed:
            self._add_section_label("Completed")
            for stage, state in completed:
                self._add_stage_line(stage, state, "\u2713", backtrack=True)
            self._content.AddSpacer(8)

        if terminal_analysis_complete:
            stage = plan.stage("analyse_series")
            self._add_terminal_analysis_card(stage, states_by_key.get("analyse_series"), action_labels["analyse_series"])
        elif recommended:
            stage = plan.stage(recommended)
            self._add_next_step_card(stage, states_by_key.get(recommended), action_labels[recommended])
        else:
            done = wx.StaticBoxSizer(wx.VERTICAL, self, "Workflow")
            done_parent = done.GetStaticBox()
            heading = wx.StaticText(done_parent, label="No required next step")
            heading.SetFont(heading.GetFont().Bold())
            done.Add(heading, 0, wx.ALL, 8)
            detail = wx.StaticText(
                done_parent,
                label="The required workflow is complete, or only optional/advanced actions remain.",
            )
            detail.Wrap(760)
            done.Add(detail, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
            self._content.Add(done, 0, wx.EXPAND | wx.BOTTOM, 14)

        if later:
            self._add_section_label("Later / optional")
            markers = {
                StageStatus.READY: "\u25cb",
                StageStatus.BLOCKED: "\u00b7",
                StageStatus.OPTIONAL: "\u25cb",
                StageStatus.WARNING: "!",
            }
            for stage, state in later:
                marker = markers.get(state.status, "\u00b7")
                self._add_stage_line(stage, state, marker)
            self._content.AddSpacer(8)

        advanced = wx.StaticBoxSizer(wx.VERTICAL, self, "Advanced tools")
        advanced_parent = advanced.GetStaticBox()
        advanced_text = wx.StaticText(
            advanced_parent,
            label=(
                "Use the detailed NMR, fitting and specialist tabs for manual controls, "
                "diagnostics and optional refinement. Workflow remains synchronised with results saved there."
            ),
        )
        advanced_text.Wrap(760)
        advanced.Add(advanced_text, 0, wx.ALL, 8)
        summary_button = wx.Button(advanced_parent, label="Summarise project...")
        summary_button.SetToolTip("Generate the same project summary PDF available from the NMR tab.")
        summary_button.Bind(wx.EVT_BUTTON, self._on_summarise_project)
        advanced.Add(summary_button, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        self._content.Add(advanced, 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 8)

        self.Layout()
        self.FitInside()
        self.SendSizeEvent()
