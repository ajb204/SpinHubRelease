"""Read-only completion/status evaluation for the guided workflow."""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any
from .model import StageRequirement, WorkflowPlan

class StageStatus(str, Enum):
    BLOCKED="blocked"; READY="ready"; COMPLETE="complete"; OPTIONAL="optional"; WARNING="warning"

@dataclass(frozen=True)
class StageState:
    key: str
    status: StageStatus
    detail: str = ""



def _legacy_workspace(context):
    """Return the explicit legacy NMR bridge, accepting old test/host contexts."""
    if context is None:
        return None
    workspace = getattr(context, "nmr_workspace", None)
    if workspace is not None:
        return workspace
    workspace = getattr(context, "legacy_nmr_workspace", None)
    if workspace is not None:
        return workspace
    return getattr(context, "tabOne", None)

def _exists(path):
    try: return bool(path) and Path(path).is_file()
    except (OSError, TypeError, ValueError): return False

def _spectrum_file_exists(state):
    fn=getattr(state,"spectrum_path",None); return _exists(fn() if callable(fn) else "")

def _spectrum_loaded(store, context=None):
    spectra=getattr(store,"spectra",{}) or {}; raw=spectra.get("raw",{}) if isinstance(spectra,dict) else {}
    if raw and raw.get("data") is not None: return True
    if getattr(store,"data",None) is not None: return True
    tab=_legacy_workspace(context)
    return bool(getattr(tab,"READ",0) and getattr(tab,"data",None) is not None)

def _has_peak_list(state, store, key):
    lists=getattr(store,"peak_lists",{}) or {}; payload=lists.get(key,{}) if isinstance(lists,dict) else {}
    if payload and (payload.get("peaks") or payload.get("rows")): return True
    fn=getattr(state,"reference_peak_path" if key=="reference" else "full_peak_path",None)
    return _exists(fn() if callable(fn) else "")



def _pseudo2d_expected_peak_names(state, store, context=None):
    """Return Full 1D peak names, which are pseudo2D's fit/reference set."""
    names = []
    lists = getattr(store, "peak_lists", {}) or {}
    payload = lists.get("full", {}) if isinstance(lists, dict) else {}
    for peak in payload.get("peaks", []) or payload.get("records", []) or []:
        name = getattr(peak, "name", None)
        if name is None and isinstance(peak, dict):
            name = peak.get("name") or peak.get("Name")
        if name is not None and str(name).strip():
            names.append(str(name).strip())
    if not names:
        for row in payload.get("rows", []) or []:
            if row:
                names.append(str(row[0]).strip())
    if names:
        return names

    path = ""
    tab = _legacy_workspace(context)
    box = getattr(tab, "fullPeakBox", None) if tab is not None else None
    if box is not None:
        try:
            value = str(box.GetValue() or "").strip()
            resolver = getattr(tab, "_resolve_spec_file", None)
            path = resolver(value) if value and callable(resolver) else value
        except Exception:
            path = ""
    if not path:
        fn = getattr(state, "full_peak_path", None)
        path = fn() if callable(fn) else ""
    try:
        with open(path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 2:
                    continue
                try:
                    float(fields[1])
                except ValueError:
                    continue
                names.append(fields[0])
    except (OSError, TypeError):
        pass
    return names

def _physical2d_fit_files_complete(state, store, context=None):
    """True when every Full 2D peak has the .dat/.out pair used by Fitting."""
    names = []
    lists = getattr(store, "peak_lists", {}) or {}
    payload = lists.get("full", {}) if isinstance(lists, dict) else {}
    for peak in payload.get("peaks", []) or []:
        name = getattr(peak, "name", None)
        if name is not None and str(name).strip():
            names.append(str(name).strip())
    tab = _legacy_workspace(context)
    if not names:
        path = ""
        box = getattr(tab, "fullPeakBox", None) if tab is not None else None
        try:
            value = str(box.GetValue() or "").strip() if box is not None else ""
            resolver = getattr(tab, "_resolve_spec_file", None)
            path = resolver(value) if value and callable(resolver) else value
        except Exception:
            path = ""
        if not path:
            fn = getattr(state, "full_peak_path", None)
            path = fn() if callable(fn) else ""
        try:
            with open(path) as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) >= 3:
                        try:
                            float(fields[1]); float(fields[2])
                        except ValueError:
                            continue
                        names.append(fields[0])
        except (OSError, TypeError):
            pass
    if not names or tab is None:
        return False
    get_dir = getattr(tab, "get_fuda_dir", None)
    fit_dir = Path(get_dir()) if callable(get_dir) else None
    return bool(fit_dir) and all((fit_dir / (name + ".dat")).is_file() and
                                 (fit_dir / (name + ".out")).is_file() for name in names)

def _pseudo2d_fit_files_complete(state, store, context=None, pseudo_tab=None):
    """True only when every Full 1D peak has both fitting result files."""
    expected = _pseudo2d_expected_peak_names(state, store, context)
    if not expected:
        return False
    fit_dir = ""
    fitting = getattr(pseudo_tab, "fittingFrame", None) if pseudo_tab is not None else None
    if fitting is not None:
        fit_dir = getattr(fitting, "fit_dir", "") or ""
    tab = _legacy_workspace(context)
    if not fit_dir and tab is not None:
        getter = getattr(tab, "get_fuda_dir", None)
        if callable(getter):
            try: fit_dir = getter()
            except Exception: fit_dir = ""
    if not fit_dir or not Path(fit_dir).is_dir():
        return False
    return all((Path(fit_dir) / (name + ".dat")).is_file() and
               (Path(fit_dir) / (name + ".out")).is_file()
               for name in expected)

def _peak_shape_complete(store, context=None):
    metadata=getattr(store,"metadata",{}) or {}; analysis=getattr(store,"analysis",{}) or {}
    # The system file is authoritative across cold loads. Reading the spectrum
    # can rebuild/clear transient DataStore metadata, so also consult the
    # persisted flag restored onto the main NMR tab.
    tab=_legacy_workspace(context)
    return bool(metadata.get("peak_shape_determined") or analysis.get("peak_shape_determined")
                or getattr(tab,"peak_shape_fitted",False))

def evaluate_workflow(plan: WorkflowPlan, state: Any, store: Any=None, context: Any=None):
    """Evaluate status from observable evidence without mutating legacy state."""
    if store is None: store=object()
    spectrum_file=_spectrum_file_exists(state); loaded=_spectrum_loaded(store, context)
    # An explicitly configured input path is a persistence contract: if that
    # file is missing, an in-memory array must not make "Prepare spectrum"
    # complete. For transient/test projects with no input file configured,
    # loaded DataStore evidence is sufficient.
    configured_input = bool(str(getattr(state, "input_file", "") or "").strip())
    spectrum = spectrum_file or (loaded and not configured_input)
    reference=_has_peak_list(state,store,"reference"); full=_has_peak_list(state,store,"full"); shape=_peak_shape_complete(store, context)
    # Pseudo2D and physical 2D each have one authoritative spectral peak
    # collection: their Full list also supplies reference-peak evidence.
    # Higher-dimensional journeys retain an independent reference projection.
    if (plan.mode.has_pseudo_axis and plan.mode.spectral_dimensions == 1) or (not plan.mode.has_pseudo_axis and plan.mode.spectral_dimensions == 2):
        reference = full
        tab = _legacy_workspace(context)
        box = getattr(tab, "fullPeakBox", None) if (plan.mode.has_pseudo_axis and tab is not None) else None
        if box is not None:
            try:
                value = str(box.GetValue() or "").strip()
                resolver = getattr(tab, "_resolve_spec_file", None)
                path = resolver(value) if value and callable(resolver) else value
                reference = _exists(path)
            except Exception:
                reference = False
    analysis = getattr(store, "analysis", {}) or {}
    pseudo_tab = None
    if context is not None:
        get_page = getattr(context, "get_page_by_title", None)
        if callable(get_page):
            # Prefer the fitting workspace for 2D+pseudo; fall back to the
            # legacy Pseudo2D workspace for 1D+pseudo.
            pseudo_tab = get_page("Fitting") or get_page("Pseudo2D")
        else:
            pseudo_tab = getattr(context, "tabPseudo", None)
    # Pseudo-series evidence is deliberately observational.  Existing specialist
    # tabs remain authoritative; the workflow only recognises results they have
    # already produced.
    # A revised peak shape invalidates every previously extracted intensity
    # result. Treat the persisted stale flag as authoritative even when old
    # FUDA files or an already-open fitting tab still contain results.
    tab_one = _legacy_workspace(context)
    peak_pick_stale = bool(analysis.get("peak_pick_stale") or getattr(tab_one, "peak_pick_stale", False))
    picked_peaks_reviewed = bool(analysis.get("picked_peaks_reviewed"))
    fitting_results_ready = bool(analysis.get("fitting_results_ready"))
    if (not fitting_results_ready and not plan.mode.has_pseudo_axis
            and plan.mode.spectral_dimensions == 2 and context is not None):
        fitting_results_ready = _physical2d_fit_files_complete(state, store, context)
    fitting_results_reviewed = bool(analysis.get("fitting_results_reviewed"))
    pseudo_intensities_stale = bool(getattr(tab_one, "pseudo_intensities_stale", False))
    pseudo_intensities = False if pseudo_intensities_stale else bool(analysis.get("pseudo_intensities_ready"))
    if (not pseudo_intensities_stale and plan.mode.has_pseudo_axis
            and plan.mode.spectral_dimensions == 1 and context is not None):
        # In the live GUI, Pseudo2D completion is stricter than the persisted
        # ready flag: every Full 1D peak must have the .dat/.out pair consumed
        # by the fitting workspace. In headless/model evaluation there is no
        # fitting workspace to inspect, so DataStore's explicit evidence API
        # remains authoritative.
        pseudo_intensities = _pseudo2d_fit_files_complete(state, store, context, pseudo_tab)
    elif not pseudo_intensities_stale and not pseudo_intensities and pseudo_tab is not None:
        intensities = getattr(pseudo_tab, "intensities", None)
        scatters = getattr(pseudo_tab, "scatter_data", None)
        has_fuda = getattr(pseudo_tab, "_has_fuda_results", None)
        pseudo_intensities = bool(intensities) or bool(scatters)
        if not pseudo_intensities and callable(has_fuda):
            try: pseudo_intensities = bool(has_fuda())
            except Exception: pass
    pseudo_reviewed = bool(analysis.get("pseudo_series_reviewed"))
    pseudo_analysis = bool(analysis.get("pseudo_analysis_complete"))
    downstream_analysis = str(analysis.get("downstream_analysis", "") or "").strip()
    if not downstream_analysis and tab_one is not None:
        downstream_analysis = str(getattr(tab_one, "downstream_analysis", "") or "").strip()
    try:
        reference_required = plan.stage("reference_peaks").requirement is StageRequirement.REQUIRED
    except KeyError:
        reference_required = False
    out=[]
    for stage in plan.stages:
        if stage.key=="dataset": out.append(StageState(stage.key,StageStatus.COMPLETE,"Dataset topology is defined."))
        elif stage.key=="spectrum": out.append(StageState(stage.key,StageStatus.COMPLETE if spectrum else StageStatus.READY,"Spectrum file exists on disk." if spectrum else "Create a spectrum file to continue."))
        elif stage.key=="peak_shape":
            if not spectrum: out.append(StageState(stage.key,StageStatus.BLOCKED,"Prepare or select a processed spectrum first."))
            elif shape: out.append(StageState(stage.key,StageStatus.COMPLETE,"Peak-shape determination is recorded."))
            elif loaded: out.append(StageState(stage.key,StageStatus.READY,"Peak shape has not yet been explicitly confirmed."))
            else: out.append(StageState(stage.key,StageStatus.READY,"The processed spectrum will be loaded automatically before peak-shape fitting opens."))
        elif stage.key=="reference_peaks":
            if not spectrum: status,detail=StageStatus.BLOCKED,"Prepare the spectrum first."
            elif stage.requirement is StageRequirement.REQUIRED and plan.mode.has_pseudo_axis and not shape: status,detail=StageStatus.BLOCKED,"Determine the peak shape first."
            elif reference: status,detail=StageStatus.COMPLETE,("A Full 1D peak list is available." if plan.mode.has_pseudo_axis and plan.mode.spectral_dimensions == 1 else "A reference peak list is available.")
            elif stage.requirement is StageRequirement.OPTIONAL: status,detail=StageStatus.OPTIONAL,"Reference peaks are optional for this dataset."
            else: status,detail=StageStatus.READY,"Create or load reference peaks if this workflow uses them."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="peak_pick":
            if full and not peak_pick_stale: status,detail=StageStatus.COMPLETE,"A current full-dimensional peak list is available."
            elif not spectrum: status,detail=StageStatus.BLOCKED,"Prepare the spectrum first."
            elif reference_required and not reference: status,detail=StageStatus.BLOCKED,"Establish the required reference peaks first."
            elif peak_pick_stale: status,detail=StageStatus.READY,"Peak shape changed; run main peak picking again."
            else: status,detail=StageStatus.READY,"Main peak picking can be run."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="fit_spectrum":
            if fitting_results_ready:
                status, detail = StageStatus.COMPLETE, "Restrained 2D fitting results are available."
            elif not full:
                status, detail = StageStatus.BLOCKED, "Create the main peak list before spectrum fitting."
            elif plan.mode.spectral_dimensions == 2 and not plan.mode.has_pseudo_axis and not picked_peaks_reviewed:
                status, detail = StageStatus.BLOCKED, "Review the picked peaks and mark them as checked before fitting."
            elif stage.requirement is StageRequirement.OPTIONAL:
                status, detail = StageStatus.OPTIONAL, "Optional refinement using the selected peak-shape model."
            else:
                status, detail = StageStatus.READY, "Run Recon with Fit and Use 2D peak list selected."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="review_fitting":
            if fitting_results_reviewed and fitting_results_ready:
                status, detail = StageStatus.COMPLETE, "The current fitting results have been inspected."
            elif fitting_results_ready:
                status, detail = StageStatus.READY, "Open the Fitting tab and fitting window to inspect the fitted peaks."
            else:
                status, detail = StageStatus.BLOCKED, "Fit the checked peaks first."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="review_peaks":
            current_full = bool(full and not peak_pick_stale)
            if picked_peaks_reviewed and current_full:
                status, detail = StageStatus.COMPLETE, "The current picked peaks have been reviewed and accepted."
            elif current_full:
                status, detail = StageStatus.READY, "Open the full peak list and 2D slice viewer, then mark the picked peaks as checked."
            else:
                status, detail = StageStatus.BLOCKED, "Run main peak picking first."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="extract_intensities":
            if pseudo_intensities: status,detail=StageStatus.COMPLETE,"Fitted/extracted pseudo-axis intensity data are available."
            elif not spectrum: status,detail=StageStatus.BLOCKED,"Prepare the spectrum stack first."
            elif not shape: status,detail=StageStatus.BLOCKED,"Determine the peak shape first."
            elif not reference: status,detail=StageStatus.BLOCKED,"Select or load the reference peaks/frequencies first."
            else: status,detail=StageStatus.READY,"Open the existing pseudo-dimensional fitting workspace and extract intensities."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="review_series":
            if pseudo_reviewed: status,detail=StageStatus.COMPLETE,"The pseudo-axis intensity series has been reviewed."
            elif pseudo_intensities: status,detail=StageStatus.READY,"Inspect intensity versus pseudo-axis slice before downstream analysis."
            else: status,detail=StageStatus.BLOCKED,"Extract/fitted intensities are required first."
            out.append(StageState(stage.key,status,detail))
        elif stage.key=="analyse_series":
            # For the physical 2D+pseudo workflow, choosing/confirming the
            # downstream analysis type is the terminal workflow decision.
            # The specialist analysis window remains available for the actual
            # scientific fitting, but Workflow itself is complete at selection.
            analysis_selected = bool(plan.mode.spectral_dimensions == 2 and plan.mode.has_pseudo_axis and downstream_analysis)
            if pseudo_analysis or analysis_selected: status,detail=StageStatus.COMPLETE,"Analysis type '%s' is selected; the guided workflow is complete." % downstream_analysis if downstream_analysis else "Downstream pseudo-axis analysis results are recorded."
            elif pseudo_intensities: status,detail=StageStatus.READY,"Open the existing analysis tools for diffusion, relaxation or experiment-specific fitting."
            else: status,detail=StageStatus.BLOCKED,"Extract/fitted intensities are required first."
            out.append(StageState(stage.key,status,detail))
        else: out.append(StageState(stage.key,StageStatus.BLOCKED,"This workflow stage is not yet available."))
    return tuple(out)

def available_actions(states):
    """Return every action that can be started now, independent of display order."""
    return tuple(item.key for item in states if item.status is StageStatus.READY)

def recommended_action(plan, states):
    """Choose the best next action using explicit workflow recommendation ranks."""
    ready = set(available_actions(states))
    candidates = [stage for stage in plan.stages if stage.key in ready]
    if not candidates:
        return None
    return min(candidates, key=lambda stage: (stage.recommendation_rank, stage.key)).key

def next_action(states, plan=None):
    """Compatibility wrapper; pass ``plan`` for order-independent recommendation."""
    if plan is not None:
        return recommended_action(plan, states)
    for item in states:
        if item.status is StageStatus.READY:
            return item.key
    return None
