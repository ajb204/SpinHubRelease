#!/usr/bin/python
"""Dedicated projections window for indirect-dimension phasing."""
from __future__ import annotations

import glob
import re
import itertools
import os
import tempfile
from functools import partial
from pathlib import Path

import matplotlib
matplotlib.use('WXAgg')
import nmrglue as ng
import numpy as np
import wx
try:
    from scipy.signal import hilbert as scipy_hilbert
except Exception:
    scipy_hilbert = None
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from spinDecon.gui.plotting.toolbar import RedrawNavigationToolbar
from matplotlib.figure import Figure

from spinDecon.gui.plotting.display_utils import blit_artists
from spinDecon.project.parameter_store import parse_float, parse_int, parse_value


class ProjectionsFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(
            self,
            parent=parent,
            title='Projections',
            size=(720, 720),
            style=wx.DEFAULT_FRAME_STYLE,
        )
        self.process_parent = parent
        self.state = getattr(parent, 'state', None)
        self._initializing = True
        self._background_ready = False
        self._projection_cache: list[dict] = []
        self._phase_controls: list[dict] = []
        self._axis_backgrounds: dict[int, object] = {}
        self._redraw_timer = None
        self._debug_messages: list[str] = []
        self._projection_view_mode_name = 'matched'
        self.maxval = 0.0
        self._contour_baseline_ready = False
        self._syncing_contour_controls = False

        self.main_panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self.main_panel, -1, self.fig)
        self.toolbar = RedrawNavigationToolbar(
            self.canvas, self.redraw_view,
            contour_callback=self.OnContourButton,
            sliders_callback=self.OnSlidersButton,
            reprocess_callback=self.OnReprocess,
            coordinates=False,
        )
        self.toolbar.Realize()
        self.statusbar = self.CreateStatusBar()

        self._build_controls()
        self._build_layout()
        self._set_default_values()
        self._bind_events()

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_MOVE, self.OnMove)
        self.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.SetMinSize((620, 420))
        self._set_compact_frame_size()

        self.load_projections()
        self._initializing = False

    # ------------------------------------------------------------------
    # Small helpers
    def _project_dir(self) -> str:
        state = getattr(self.process_parent, 'state', None)
        return state.spec_dir() if state is not None else ''


    def _debug_widget_state(self, label: str) -> None:
        try:
            frame_sizer = self.GetSizer()
            panel_sizer = self.main_panel.GetSizer() if hasattr(self, 'main_panel') else None
            canvas_parent = self.canvas.GetParent() if hasattr(self, 'canvas') else None
            toolbar_parent = self.toolbar.GetParent() if hasattr(self, 'toolbar') else None
            self._debug_projection(
                f"{label}: frame_shown={self.IsShown()!r} main_panel_shown={self.main_panel.IsShown()!r} "
                f"canvas_shown={self.canvas.IsShown()!r} toolbar_shown={self.toolbar.IsShown()!r} "
                f"frame_size={tuple(self.GetSize())!r} panel_size={tuple(self.main_panel.GetSize())!r} "
                f"canvas_size={tuple(self.canvas.GetSize())!r} toolbar_size={tuple(self.toolbar.GetSize())!r} "
                f"frame_sizer={type(frame_sizer).__name__ if frame_sizer else None!r} "
                f"panel_sizer={type(panel_sizer).__name__ if panel_sizer else None!r} "
                f"canvas_parent={type(canvas_parent).__name__ if canvas_parent else None!r} "
                f"toolbar_parent={type(toolbar_parent).__name__ if toolbar_parent else None!r}"
            )
        except Exception as exc:
            self._debug_projection(f"{label}: widget state probe failed: {exc!r}")

    def _raw_dir(self) -> str:
        state = getattr(self.process_parent, 'state', None)
        return state.raw_dir() if state is not None else ''


    def _parameter_file_name(self) -> str:
        for obj in (self.process_parent, getattr(self.process_parent, 'parent', None)):
            if obj is not None and hasattr(obj, 'deconParFile'):
                try:
                    return str(obj.deconParFile).strip()
                except Exception:
                    continue
        return 'decon.par'

    def _parameter_file_candidates(self) -> list[str]:
        project_dir = self._project_dir()
        raw_dir = self._raw_dir()
        raw_name = self._parameter_file_name()
        candidates: list[str] = []

        for value in (raw_name,):
            if not value:
                continue
            candidates.append(value)
            for base in (project_dir, raw_dir):
                if base:
                    if not os.path.isabs(value):
                        candidates.append(os.path.join(base, value))
                    candidates.append(os.path.join(base, os.path.basename(value)))

        for base in (project_dir, raw_dir):
            if base:
                candidates.append(os.path.join(base, 'decon.par'))

        seen: set[str] = set()
        out: list[str] = []
        for cand in candidates:
            if not cand:
                continue
            norm = os.path.normpath(cand)
            if norm not in seen:
                seen.add(norm)
                out.append(norm)
        return out

    def _parameter_file_path(self) -> str:
        """Return the project parameter file from shared ProjectState."""
        state = getattr(self.process_parent, 'state', None)
        if state is not None:
            if state.parameter_file:
                return state.parameter_file
            return os.path.join(state.working_dir or '.', self._parameter_file_name())
        candidates = self._parameter_file_candidates()
        return candidates[0] if candidates else self._parameter_file_name()


    def _raw_dimension_value(self):
        for obj in (self.process_parent, getattr(self.process_parent, 'parent', None)):
            if obj is None:
                continue
            dim = getattr(obj, 'dim', None)
            if dim is not None:
                return dim
        return None

    def _spectral_dimension_count(self) -> int:
        dim = self._raw_dimension_value()
        try:
            return max(1, int(dim))
        except Exception:
            if isinstance(dim, str) and dim.endswith('p') and dim[:-1].isdigit():
                # Pseudo dimensions are named by the full spectral count, e.g.
                # 2p and 3p, so do not subtract 1 here.
                return max(1, int(dim[:-1]))
        for obj in (self.process_parent, getattr(self.process_parent, 'parent', None)):
            if obj is None:
                continue
            if hasattr(obj, '_spectral_dimension_count'):
                try:
                    return max(1, int(obj._spectral_dimension_count()))
                except Exception:
                    pass
        return 1

    def _parameter_value(self, key: str, default: str = '') -> str:
        try:
            value = parse_value(self._parameter_file_path(), key, default=default)
        except Exception:
            value = default
        if value is None:
            return default
        return str(value).strip()

    def _spectral_labels(self) -> list[str]:
        count = self._spectral_dimension_count()
        # Process owns the canonical live labels for the active session.  This
        # keeps Projection independent of whether Conversion has ever opened.
        getter = getattr(self.process_parent, 'get_spectral_labels', None)
        if callable(getter):
            try:
                labels = [str(x).strip() for x in getter() if str(x).strip()]
                if len(labels) >= count:
                    return labels[:count]
            except Exception:
                pass

        # Compatibility fallback for isolated/test callers without ProcessFrame.
        state = getattr(self, 'state', None)
        live = getattr(state, 'gui_settings', {}) if state is not None else {}
        labels = []
        for idx in range(1, count + 1):
            value = str(live.get(f'label{idx}', '')).replace(' ', '').strip()
            if not value:
                value = self._parameter_value(f'label{idx}', '')
            if value:
                labels.append(value)
        if len(labels) < count:
            labels.extend([f'Dim {i + 1}' for i in range(len(labels), count)])
        return labels[:count]

    def _processing_phase_defaults(self, label: str) -> dict[str, float]:
        savefile = self._parameter_file_path()

        # In 2D the projection window always represents the indirect dimension.
        # Its launch-time controls must therefore come directly from p0_1/p1_1.
        if self._spectral_dimension_count() == 2:
            p0 = parse_float(savefile, 'p0_1', 0.0)
            p1 = parse_float(savefile, 'p1_1', 0.0)
            self._debug_projection(
                f'_processing_phase_defaults 2D label={label!r} -> p0_1={p0!r} p1_1={p1!r} savefile={savefile!r}'
            )
            try:
                return {'p0': float(p0), 'p1': float(p1)}
            except Exception:
                return {'p0': 0.0, 'p1': 0.0}

        labels = self._spectral_labels()
        try:
            index = labels.index(label)
        except ValueError:
            # tolerate stale in-memory labels by comparing stripped values
            normalized = [str(x).strip() for x in labels]
            try:
                index = normalized.index(str(label).strip())
            except ValueError:
                index = -1

        if index <= 0:
            return {'p0': 0.0, 'p1': 0.0}

        suffix = str(index)
        p0 = parse_float(savefile, f'p0_{suffix}', 0.0)
        p1 = parse_float(savefile, f'p1_{suffix}', 0.0)
        try:
            return {'p0': float(p0), 'p1': float(p1)}
        except Exception:
            return {'p0': 0.0, 'p1': 0.0}

    def _raw_dimension_labels(self) -> tuple[list[str], str]:
        """Return raw labels from the authoritative Process-session store."""
        getter = getattr(self.process_parent, 'get_dimension_labels', None)
        if callable(getter):
            try:
                labels = [str(x).strip() for x in getter() if str(x).strip()]
                if labels:
                    return labels[:2], 'process_label_store'
            except Exception:
                pass
        # Compatibility fallback for isolated tests/legacy callers.
        labels = [
            str(self._parameter_value('label1', '')).replace(' ', '').strip(),
            str(self._parameter_value('label2', '')).replace(' ', '').strip(),
        ]
        labels = [x for x in labels if x]
        return (labels[:2] if labels else ['Direct', 'Indirect']), 'deconParFile' if labels else 'default'

    def _raw_spectrum_bundle(self):
        """Return the raw spectrum cached in the shared state, if available."""
        def _preferred_labels() -> list[str]:
            labels, source = self._raw_dimension_labels()
            self._debug_projection(f'_raw_spectrum_bundle preferred labels source={source!r} labels={labels!r}')
            return labels

        if self.state is not None:
            try:
                raw = getattr(self.state, 'spectra', {}).get('raw', {}) or {}
            except Exception:
                raw = {}
            dic = raw.get('dic')
            data = raw.get('data')
            if dic is not None and data is not None:
                dic = dict(dic)
                if not str(dic.get('FDF1LABEL', '')).strip() or not str(dic.get('FDF2LABEL', '')).strip():
                    try:
                        raw_path = str(raw.get('spectrumfile', '')).strip()
                    except Exception:
                        raw_path = ''
                    if raw_path and os.path.exists(raw_path):
                        try:
                            file_dic, file_data = ng.pipe.read(raw_path)
                            if file_dic is not None:
                                dic.update(dict(file_dic))
                                self._debug_projection(
                                    f'_raw_spectrum_bundle refreshed FDF metadata from file path={raw_path!r} f1={str(dic.get("FDF1LABEL", "")).strip()!r} f2={str(dic.get("FDF2LABEL", "")).strip()!r}'
                                )
                            if file_data is not None and np.asarray(data).shape != np.asarray(file_data).shape:
                                data = file_data
                        except Exception as exc:
                            self._debug_projection(f'_raw_spectrum_bundle file metadata refresh failed for {raw_path!r}: {exc!r}')
                return dic, np.asarray(data), _preferred_labels()

        for obj in (self.process_parent, getattr(self.process_parent, 'parent', None)):
            if obj is None:
                continue
            dic = getattr(obj, 'dic', None)
            data = getattr(obj, 'data', None)
            labels = getattr(obj, 'labb', None)
            if dic is not None and data is not None:
                if isinstance(labels, (list, tuple)) and labels:
                    labels = [str(x) for x in labels]
                else:
                    labels = _preferred_labels()
                return dic, np.asarray(data), labels

        return None, None, _preferred_labels()

    def _special_2d_axis_labels(self, source_dic) -> tuple[str, str, dict[str, object]]:
        """Return the raw 2D axis labels using FDF1/FDF2 metadata when available."""
        source_dic = dict(source_dic or {})
        f1_label = str(source_dic.get('FDF1LABEL', '')).strip()
        f2_label = str(source_dic.get('FDF2LABEL', '')).strip()
        f1_size = source_dic.get('FDF1FTSIZE', None)
        f2_size = source_dic.get('FDF2FTSIZE', None)
        direct_label = f2_label or f1_label or 'Direct'
        indirect_label = f1_label or f2_label or 'Indirect'
        info = {
            'f1_label': f1_label,
            'f2_label': f2_label,
            'f1_size': f1_size,
            'f2_size': f2_size,
            'direct_label': direct_label,
            'indirect_label': indirect_label,
        }
        self._debug_projection(f"_special_2d_axis_labels info={info!r}")
        return direct_label, indirect_label, info

    def _special_2d_projection_phase_context(
        self,
        source_dic,
        projected_1d: np.ndarray,
        *,
        direct_label: str,
        indirect_label: str,
        source_path: str | None = None,
    ) -> tuple[dict, np.ndarray | None, dict[str, object]]:
        """Normalize the header for the special raw 2D projection path."""
        projected_1d = np.asarray(projected_1d)
        source_dic = dict(source_dic or {})
        direct_label = str(direct_label).strip() or 'Direct'
        indirect_label = str(indirect_label).strip() or 'Indirect'

        def _as_int(value):
            try:
                return int(round(float(value)))
            except Exception:
                return None

        projected_size = int(projected_1d.shape[-1] if projected_1d.ndim else projected_1d.size)
        f1_size = _as_int(source_dic.get('FDF1FTSIZE'))
        f2_size = _as_int(source_dic.get('FDF2FTSIZE'))
        f1_label = str(source_dic.get('FDF1LABEL', '')).strip()
        f2_label = str(source_dic.get('FDF2LABEL', '')).strip()

        keep_dim = None
        # The 1D trace in the special 2D path is a projection retained on
        # the DIRECT dimension.  Match its size/label to the direct axis,
        # not the indirect axis (older 2D code had these reversed).
        if f2_label == direct_label and (f2_size == projected_size or f2_size is None):
            keep_dim = 2.0
        elif f1_label == direct_label and (f1_size == projected_size or f1_size is None):
            keep_dim = 1.0
        elif f2_size == projected_size and f1_size != projected_size:
            keep_dim = 2.0
        elif f1_size == projected_size and f2_size != projected_size:
            keep_dim = 1.0
        elif f1_size == projected_size and f2_size == projected_size:
            if f2_label == direct_label and f1_label != direct_label:
                keep_dim = 2.0
            elif f1_label == direct_label and f2_label != direct_label:
                keep_dim = 1.0
        if keep_dim is None:
            label_to_dim = {}
            for idx in range(1, 3):
                label = source_dic.get(f'FDF{idx}LABEL')
                if label is not None:
                    label_to_dim[str(label).strip()] = float(idx)
            keep_dim = label_to_dim.get(direct_label) or label_to_dim.get(indirect_label) or 2.0

        keep_dim = float(keep_dim)
        other_dim = 2.0 if keep_dim == 1.0 else 1.0

        special_dic = dict(source_dic)
        special_dic[f'FDF{int(keep_dim)}LABEL'] = direct_label
        special_dic[f'FDF{int(other_dim)}LABEL'] = indirect_label
        special_dic[f'FDF{int(keep_dim)}FTSIZE'] = projected_size

        self._debug_projection(
            f"_special_2d_projection_phase_context source_path={source_path!r} projected_size={projected_size!r} f1_size={f1_size!r} f2_size={f2_size!r} f1_label={f1_label!r} f2_label={f2_label!r} keep_dim={keep_dim!r} other_dim={other_dim!r} direct_label={direct_label!r} indirect_label={indirect_label!r}"
        )

        phase_dic = self._projection_1d_phase_dic(
            special_dic,
            projected_1d,
            direct_label=direct_label,
            indirect_label=indirect_label,
            source_path=source_path,
        )

        trace_scale = None
        try:
            trace_uc = ng.pipe.make_uc(phase_dic, projected_1d, dim=0)
            trace_scale = np.asarray(trace_uc.ppm_scale())
            self._debug_projection(
                f"_special_2d_projection_phase_context trace_scale built: size={trace_scale.size!r} first_last={(float(trace_scale[0]), float(trace_scale[-1])) if trace_scale.size else None!r}"
            )
        except Exception as exc:
            self._debug_projection(f"_special_2d_projection_phase_context trace_scale build failed: {exc!r}")

        info = {
            'projected_size': projected_size,
            'f1_size': f1_size,
            'f2_size': f2_size,
            'f1_label': f1_label,
            'f2_label': f2_label,
            'keep_dim': keep_dim,
            'other_dim': other_dim,
            'source_path': source_path,
            'trace_scale_size': int(trace_scale.size) if trace_scale is not None else None,
            'special_labels': {'direct_label': direct_label, 'indirect_label': indirect_label},
        }
        return phase_dic, trace_scale, info

    def _special_2d_projection_layout(
        self,
        source_dic,
        data: np.ndarray,
        *,
        direct_label: str,
        indirect_label: str,
        source_path: str | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
        """Use the raw 2D spectrum orientation directly for plotting."""
        source_dic = dict(source_dic or {})
        data = np.asarray(data)
        direct_label = str(direct_label).strip() or 'Direct'
        indirect_label = str(indirect_label).strip() or 'Indirect'

        def _as_int(value):
            try:
                return int(round(float(value)))
            except Exception:
                return None

        f1_label = str(source_dic.get('FDF1LABEL', '')).strip()
        f2_label = str(source_dic.get('FDF2LABEL', '')).strip()
        f1_size = _as_int(source_dic.get('FDF1FTSIZE'))
        f2_size = _as_int(source_dic.get('FDF2FTSIZE'))

        if f1_label and f2_label:
            # Use refreshed raw-spectrum metadata only for the displayed nucleus
            # names; do not remap the 2D data array itself.
            indirect_label = f1_label
            direct_label = f2_label

        if f1_label or f2_label:
            self._debug_projection(
                f"_special_2d_projection_layout raw-spectrum FDF mapping f1={f1_label!r}/{f1_size!r} f2={f2_label!r}/{f2_size!r} direct={direct_label!r} indirect={indirect_label!r}"
            )

        try:
            # Dedicated 3p projections are now written canonically as
            # (direct rows, indirect columns), so no corrective transpose is
            # required in the projection window.
            direct_uc = ng.pipe.make_uc(source_dic, data, dim=0)
            indirect_uc = ng.pipe.make_uc(source_dic, data, dim=1)
            direct_scale = np.asarray(direct_uc.ppm_scale())
            indirect_scale = np.asarray(indirect_uc.ppm_scale())
        except Exception as exc:
            self._debug_projection(f"_special_2d_projection_layout raw make_uc failed: {exc!r}")
            raise

        proj1d_raw = np.sum(data, axis=0)
        contour_expected = (int(direct_scale.size), int(indirect_scale.size))
        transpose_for_contour = False
        self._debug_projection(
            f"_special_2d_projection_layout raw layout: data_shape={data.shape!r} contour_expected={contour_expected!r} transpose_for_contour={transpose_for_contour!r} direct_scale.size={direct_scale.size!r} indirect_scale.size={indirect_scale.size!r} proj1d_raw.size={proj1d_raw.size!r} direct_label={direct_label!r} indirect_label={indirect_label!r}"
        )
        if data.shape != contour_expected:
            self._debug_projection(
                f"_special_2d_projection_layout WARNING canonical data shape does not match contour expectation: data.shape={data.shape!r} expected={contour_expected!r}"
            )

        layout_info = {
            'source_path': source_path,
            'f1_label': f1_label,
            'f2_label': f2_label,
            'f1_size': f1_size,
            'f2_size': f2_size,
            'direct_label': direct_label,
            'indirect_label': indirect_label,
            'direct_size': None,
            'indirect_size': None,
            'selected_transpose': False,
            'selected_shape': tuple(int(x) for x in data.shape),
            'selected_score': None,
            'selected_reasons': ['canonical_3p_projection_orientation'],
            'direct_scale_size': int(direct_scale.size),
            'indirect_scale_size': int(indirect_scale.size),
            'proj1d_raw_size': int(proj1d_raw.size),
            'plot_header_f1': (str(source_dic.get('FDF1LABEL', '')).strip(), _as_int(source_dic.get('FDF1FTSIZE'))),
            'plot_header_f2': (str(source_dic.get('FDF2LABEL', '')).strip(), _as_int(source_dic.get('FDF2FTSIZE'))),
            'raw_layout': True,
            'transpose_for_contour': transpose_for_contour,
        }
        self._debug_projection(f"_special_2d_projection_layout selected raw={layout_info!r}")
        return data, direct_scale, indirect_scale, layout_info

    def _direct_label(self) -> str:
        labels = self._spectral_labels()
        if labels:
            return labels[0]
        candidate = getattr(self.process_parent, 'labb', None)
        if candidate is None:
            candidate = getattr(getattr(self.process_parent, 'parent', None), 'labb', None)
        if isinstance(candidate, (list, tuple)) and candidate:
            return str(candidate[0])
        return 'Direct'

    def _projection_specs(self) -> list[dict[str, str]]:
        labels = self._spectral_labels()
        if len(labels) < 2:
            return []
        direct = self._direct_label()
        indirect_labels = [label for label in labels if label != direct]
        return [
            {'direct': direct, 'indirect': label, 'title': f'{direct}/{label}'}
            for label in indirect_labels
        ]

    def _is_special_2d_projection_case(self) -> bool:
        return self._spectral_dimension_count() == 2

    def _process_pipefile_name(self) -> str:
        for obj in (self.process_parent, getattr(self.process_parent, 'parent', None)):
            if obj is None:
                continue
            candidate = getattr(obj, '_process_pipefile', None)
            if callable(candidate):
                try:
                    value = str(candidate()).strip()
                    if value:
                        return value
                except Exception:
                    continue
        return ''

    def _raw_projection_source_candidates(self) -> list[str]:
        candidates: list[str] = []
        base = self._project_dir()
        pipefile = self._process_pipefile_name()

        # Prefer the newly processed spectrum.  The raw spectrum remains only
        # as a fallback for first-launch / pre-processing use.
        if base and pipefile:
            candidates.append(os.path.join(base, pipefile))
        if base:
            candidates.extend([
                os.path.join(base, 'test.ft2'),
                os.path.join(base, 'slice.ft1'),
                os.path.join(base, 'test.ft'),
            ])
        candidates.extend(['test.ft2', 'slice.ft1', 'test.ft'])

        try:
            raw_spectrumfile = str((self.state.spectra.get('raw', {}) or {}).get('spectrumfile', '')).strip()
        except Exception:
            raw_spectrumfile = ''
        if raw_spectrumfile:
            candidates.append(raw_spectrumfile)
            if base and not os.path.isabs(raw_spectrumfile):
                candidates.append(os.path.join(base, raw_spectrumfile))
        seen: set[str] = set()
        out: list[str] = []
        for cand in candidates:
            if not cand:
                continue
            norm = os.path.normpath(str(cand))
            if norm not in seen:
                seen.add(norm)
                out.append(norm)
        return out

    def _special_2d_projection_entry(self) -> tuple[str, str] | None:
        dic, _, labels = self._raw_spectrum_bundle()
        if not labels:
            labels = self._spectral_labels()
        direct, indirect, axis_info = self._special_2d_axis_labels(dic)
        # On a cold project open state.spectra['raw'] is not populated yet, so
        # there is no NMRPipe dictionary from which to obtain FDF1/FDF2 labels.
        # Use the persisted spectral labels in that case.  Previously this
        # produced the synthetic name Direct.Indirect, missed the already
        # existing H1.C13.dat-style projection, and only began working after
        # the main NMR tab happened to populate shared raw-spectrum state.
        if dic is None and len(labels) >= 2:
            direct, indirect = str(labels[0]), str(labels[1])
            axis_info = dict(axis_info)
            axis_info.update({'direct_label': direct, 'indirect_label': indirect,
                              'label_source': 'persisted spectral labels'})
        title = f'{direct}.{indirect}'
        self._debug_projection(f'_special_2d_projection_entry labels={labels!r} axis_info={axis_info!r} direct={direct!r} indirect={indirect!r} title={title!r}')
        # A dataset with two spectral dimensions plus a real/pseudo axis is
        # represented here as a special 2D projection case.  The dedicated
        # pseudo-3D projector writes the spectral plane in the GUI's native
        # direct.indirect naming/orientation (for example H1.C13.dat).
        projection_candidates = []
        for proj_dir in self._projection_search_dirs():
            projection_candidates.extend([
                os.path.join(proj_dir, f'{direct}.{indirect}.dat'),
                os.path.join(proj_dir, f'{indirect}.{direct}.dat'),
            ])

        candidates = projection_candidates + self._raw_projection_source_candidates()
        self._debug_projection(f'_special_2d_projection_entry candidates={candidates!r}')
        for path in candidates:
            try:
                if os.path.exists(path):
                    self._debug_projection(f'_special_2d_projection_entry selected path={path!r}')
                    return path, title
            except Exception as exc:
                self._debug_projection(f'_special_2d_projection_entry path check failed for {path!r}: {exc!r}')
        self._debug_projection('_special_2d_projection_entry no candidate file located')
        return None

    def _special_2d_projection_bundle(self) -> dict | None:
        # Always use the same source selected by _special_2d_projection_entry.
        # After re-processing this should be the newly generated test.ft2, not
        # the cached raw array retained in shared state.
        entry = self._special_2d_projection_entry()
        dic = data = None
        selected_path = entry[0] if entry else ''
        labels = self._spectral_labels()

        if selected_path:
            try:
                dic, data = ng.pipe.read(selected_path)
                self._debug_projection(
                    f'_special_2d_projection_bundle loaded selected file path={selected_path!r} data_shape={getattr(data, "shape", None)!r}'
                )
            except Exception as exc:
                self._debug_projection(
                    f'_special_2d_projection_bundle failed to read selected path={selected_path!r}: {exc!r}'
                )
                import traceback
                self._debug_projection(traceback.format_exc())
                dic = data = None

        if dic is None or data is None:
            dic, data, raw_labels = self._raw_spectrum_bundle()
            if not labels and raw_labels:
                labels = raw_labels
            if selected_path:
                self._debug_projection(
                    f'_special_2d_projection_bundle falling back to in-memory raw data after read failure; selected_path={selected_path!r}'
                )

        if not labels:
            labels = self._spectral_labels()
        direct, indirect, axis_info = self._special_2d_axis_labels(dic)
        title = f'{direct}.{indirect}'
        self._debug_projection(
            f'_special_2d_projection_bundle labels={labels!r} axis_info={axis_info!r} direct={direct!r} indirect={indirect!r} title={title!r} selected_path={selected_path!r}'
        )
        if dic is None or data is None:
            self._debug_projection('_special_2d_projection_bundle no usable spectrum data available')
            return None

        raw_path = selected_path
        if not raw_path:
            try:
                raw = getattr(self.state, 'spectra', {}).get('raw', {}) if self.state is not None else {}
                raw_path = str((raw or {}).get('spectrumfile', '')).strip()
            except Exception:
                raw_path = ''
        if not raw_path:
            raw_path = os.path.join(self._project_dir(), self._process_pipefile_name() or 'test.ft2') or 'in-memory-raw-spectrum'

        self._debug_projection(
            f'_special_2d_projection_bundle final source path={raw_path!r} data_shape={getattr(data, "shape", None)!r} FDF1={str(dic.get("FDF1LABEL", "")).strip()!r} FDF2={str(dic.get("FDF2LABEL", "")).strip()!r}'
        )
        return {'path': raw_path, 'title': title, 'dic': dic, 'data': np.asarray(data), 'direct': direct, 'indirect': indirect}

    def _projection_dir(self) -> str:
        base = self._project_dir()
        return os.path.join(base, 'projections') if base else 'projections'

    def _projection_search_dirs(self) -> list[str]:
        """Projection artefacts are always descendants of SpecPath."""
        base = self._project_dir()
        return [
            os.path.join(base, 'projections'),
            os.path.join(base, 'projection_decon'),
            os.path.join(base, 'projections1D'),
        ] if base else []

    def _find_projection_file(self, left: str, right: str) -> str | None:
        target_names = {f'{left}.{right}.dat', f'{right}.{left}.dat'}
        target_stems = {f'{left}.{right}', f'{right}.{left}'}
        self._debug_projection(f"_find_projection_file left={left!r} right={right!r} target_names={sorted(target_names)!r} target_stems={sorted(target_stems)!r}")

        for proj_dir in self._projection_search_dirs():
            try:
                if not os.path.isdir(proj_dir):
                    continue
            except Exception:
                continue

            for name in target_names:
                path = os.path.join(proj_dir, name)
                if os.path.exists(path):
                    return path

            try:
                dat_files = sorted(glob.glob(os.path.join(proj_dir, '*.dat')))
            except Exception:
                dat_files = []
            for path in dat_files:
                stem = Path(path).stem
                parts = stem.split('.')
                if stem in target_stems:
                    return path
                if left in parts and right in parts:
                    return path

        return None

    def _expected_projection_entries(self) -> list[tuple[str, str]]:
        labels = self._spectral_labels()
        direct = self._direct_label()
        self._debug_projection(f'_expected_projection_entries labels={labels!r} direct={direct!r}')
        if len(labels) < 1:
            return []

        if self._is_special_2d_projection_case():
            special = self._special_2d_projection_entry()
            return [special] if special else []

        entries: list[tuple[str, str]] = []
        for other in labels:
            if other == direct:
                continue
            path = self._find_projection_file(direct, other)
            if path:
                entries.append((path, f'{direct}.{other}'))

        if entries:
            return entries

        # Fallback for older projection folders: keep only files that include
        # the direct dimension name and are 2D projection files.
        for proj_dir in self._projection_search_dirs():
            try:
                dat_files = sorted(glob.glob(os.path.join(proj_dir, '*.dat')))
            except Exception:
                dat_files = []
            for path in dat_files:
                stem = Path(path).stem
                parts = stem.split('.')
                if direct in parts and len(parts) >= 2:
                    other = next((p for p in parts if p != direct), parts[-1])
                    entries.append((path, f'{direct}.{other}'))
                    if len(entries) >= max(1, len(labels) - 1):
                        return entries

        return entries

    def _cached_projection_entries(self) -> list[tuple[str, str, dict]]:
        entries: list[tuple[str, str, dict]] = []
        direct = self._direct_label()
        labels = self._spectral_labels()
        if self.state is None:
            return entries
        projections = getattr(self.state, 'projections', {}) or {}
        for key, payload in projections.items():
            try:
                data = np.asarray(payload.get('data'))
            except Exception:
                continue
            if data.ndim != 2:
                continue
            labb = payload.get('labb')
            if not isinstance(labb, (list, tuple)) or len(labb) < 2:
                continue
            a, b = str(labb[0]), str(labb[1])
            if direct not in (a, b):
                continue
            other = b if a == direct else a
            if other not in labels and labels:
                # tolerate partially cached metadata, but prefer a sensible label
                other = other if other else labels[-1]
            entries.append((str(payload.get('source', key)), f'{direct}.{other}', payload))
        return entries

    def _levels(self) -> np.ndarray:
        try:
            min_level = float(self.cmin.GetValue())
        except Exception:
            min_level = 1e5
        try:
            max_factor = float(self.cfac.GetValue())
        except Exception:
            max_factor = 1.3
        try:
            count = int(float(self.cnum.GetValue()))
        except Exception:
            count = 20

        if count <= 0:
            count = 10
        if max_factor == 0:
            max_factor = 1.2
        if min_level == 0:
            min_level = 1e3

        levels = [min_level]
        for _ in range(count - 1):
            levels.append(levels[-1] * max_factor)
        levels = np.asarray(levels, dtype=float)
        return np.concatenate((-levels[::-1], levels))

    def _set_text_ctrl_value(self, ctrl, value) -> None:
        self._syncing_contour_controls = True
        try:
            ctrl.SetValue(str(value))
        finally:
            self._syncing_contour_controls = False

    def _projection_maxval(self, entries, *, special_bundle=None, cached_entries=None) -> float:
        maxval = 0.0
        if special_bundle is not None:
            try:
                arr = np.asarray(special_bundle.get('data'))
                if arr.ndim == 2 and arr.size:
                    maxval = max(maxval, float(np.nanmax(np.abs(arr))))
            except Exception:
                pass
        if cached_entries:
            for _path, _title, payload in cached_entries:
                try:
                    arr = np.asarray((payload or {}).get('data'))
                    if arr.ndim == 2 and arr.size:
                        maxval = max(maxval, float(np.nanmax(np.abs(arr))))
                except Exception:
                    continue
            return maxval
        for path, _title in entries:
            if special_bundle is not None and path == special_bundle.get('path'):
                continue
            try:
                dic, data = ng.pipe.read(path)
                arr = np.asarray(data)
                if arr.ndim == 2 and arr.size:
                    maxval = max(maxval, float(np.nanmax(np.abs(arr))))
            except Exception:
                continue
        return maxval

    def _sync_contour_cmin_from_thresh(self) -> None:
        try:
            thresh = float(self.cthresh.GetValue())
        except Exception:
            thresh = 0.1
        if thresh <= 0:
            thresh = 0.1
        if self.maxval <= 0:
            return
        self._set_text_ctrl_value(self.cmin, self.maxval * thresh)

    def _sync_contour_thresh_from_cmin(self) -> None:
        try:
            cmin = float(self.cmin.GetValue())
        except Exception:
            return
        if cmin == 0 or self.maxval <= 0:
            return
        self._set_text_ctrl_value(self.cthresh, cmin / self.maxval)

    def _sync_contour_controls_from_baseline(self) -> None:
        self._sync_contour_cmin_from_thresh()

    def _safe_data_extent(self, arr: np.ndarray) -> tuple[float, float, float, float] | None:
        if arr.size == 0:
            return None
        try:
            return float(np.nanmin(arr)), float(np.nanmax(arr)), float(np.nanmin(arr)), float(np.nanmax(arr))
        except Exception:
            return None

    def _debug_projection(self, message: str) -> None:
        # Retain a small in-memory event history for on-screen diagnostics,
        # but do not create persistent debug files.
        self._debug_messages.append(message.rstrip())

    def _debug_summary(self) -> str:
        lines = [
            f'project_dir={self._project_dir()!r}',
            f'raw_dim={self._raw_dimension_value()!r}',
            f'spectral_count={self._spectral_dimension_count()}',
            f'labels={self._spectral_labels()!r}',
            f'direct_label={self._direct_label()!r}',
            f'projection_dirs={self._projection_search_dirs()!r}',
        ]
        for proj_dir in self._projection_search_dirs():
            try:
                files = sorted(Path(proj_dir).glob('*.dat'))
                lines.append(f'{proj_dir}: {[p.name for p in files]!r}')
            except Exception as exc:
                lines.append(f'{proj_dir}: ERROR {exc!r}')
        if self._debug_messages:
            lines.append('recent_events=')
            lines.extend(f'  {m}' for m in self._debug_messages[-40:])
        return '\n'.join(lines)

    # ------------------------------------------------------------------
    def _build_controls(self) -> None:
        # Keep the existing indirect-phasing controls and event logic, but host
        # them in a compact modeless companion frame rather than consuming plot
        # space in the Process Projections frame.
        self.phaseFrame = wx.Frame(
            self, title='Projection phasing',
            style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT,
        )
        self.phasePanel = wx.Panel(self.phaseFrame)
        self.phasePanelSizer = wx.BoxSizer(wx.VERTICAL)
        self.phasePanel.SetSizer(self.phasePanelSizer)
        phase_frame_sizer = wx.BoxSizer(wx.VERTICAL)
        phase_frame_sizer.Add(self.phasePanel, 1, wx.EXPAND | wx.ALL, 4)
        self.phaseFrame.SetSizer(phase_frame_sizer)
        self.phaseFrame.Bind(wx.EVT_CLOSE, self._hide_phase_frame)

        # Contour values remain attributes of this frame (and therefore keep
        # all existing dimensionality-specific plotting/persistence behaviour),
        # while their editors live in a modeless child window.
        self.contourFrame = wx.Frame(self, title='Contours', style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        contourPanel = wx.Panel(self.contourFrame)
        self.cminLab = wx.StaticText(contourPanel, label='Min:')
        self.cthreshLab = wx.StaticText(contourPanel, label='Thresh:')
        self.cfacLab = wx.StaticText(contourPanel, label='Fac:')
        self.cnumLab = wx.StaticText(contourPanel, label='Num:')
        enter_style = wx.TE_PROCESS_ENTER
        self.cmin = wx.TextCtrl(contourPanel, size=(100, 22), style=enter_style)
        self.cthresh = wx.TextCtrl(contourPanel, size=(70, 22), style=enter_style)
        self.cfac = wx.TextCtrl(contourPanel, size=(60, 22), style=enter_style)
        self.cnum = wx.TextCtrl(contourPanel, size=(60, 22), style=enter_style)
        contourSizer = wx.BoxSizer(wx.HORIZONTAL)
        for widget in (self.cthreshLab, self.cthresh, self.cminLab, self.cmin,
                       self.cfacLab, self.cfac, self.cnumLab, self.cnum):
            contourSizer.Add(widget, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 4)
        contourPanel.SetSizerAndFit(contourSizer)
        contourFrameSizer = wx.BoxSizer(wx.VERTICAL)
        contourFrameSizer.Add(contourPanel, 1, wx.EXPAND)
        self.contourFrame.SetSizerAndFit(contourFrameSizer)
        self.contourFrame.Bind(wx.EVT_CLOSE, self._hide_contour_frame)

        # Projection actions live in Matplotlib's native wx toolbar.  This is
        # the same arrangement used by the main Projection window: keeping the
        # controls inside wx.ToolBar preserves its continuous TB_BOTTOM rule.
        self.maskLab = wx.StaticText(self.main_panel, label='Selection:')
        self.maskStrengthMin = wx.StaticText(self.main_panel, label='0')
        self.maskStrengthMax = wx.StaticText(self.main_panel, label='1')
        self.maskStrengthValue = wx.StaticText(self.main_panel, label='0.50')
        self.maskStrengthSlider = wx.Slider(
            self.main_panel, value=50, minValue=0, maxValue=100,
            size=(130, 22), style=wx.SL_HORIZONTAL,
        )
        self.toolbar.AddSeparator()
        toolbar_controls = (
            self.maskLab, self.maskStrengthMin, self.maskStrengthSlider,
            self.maskStrengthMax, self.maskStrengthValue,
        )
        for widget in toolbar_controls:
            widget.Reparent(self.toolbar)
            self.toolbar.AddControl(widget)

        self.toolbar.bind_control_status_help(self.maskStrengthSlider, 'Selection strictness')

        # Restore Matplotlib's coordinate readout at the far right after our
        # controls, rather than allowing it to split the toolbar row.
        self.toolbar.AddStretchableSpace()
        self.toolbar._coordinates = True
        self.toolbar._label_text = wx.StaticText(self.toolbar, style=wx.ALIGN_LEFT)
        self.toolbar.AddControl(self.toolbar._label_text)
        self.toolbar.Realize()

        self._bind_live_contour_controls()

    def _bind_live_contour_controls(self) -> None:
        """Mirror persistent contour edits into shared live state.

        Projection phase sliders remain deliberately transient and are handled
        by the separate projection preview state until Re-process promotes them.
        """
        def changed(event):
            try:
                self.collect_updates()
            except Exception:
                pass
            try:
                event.Skip()
            except Exception:
                pass

        for ctrl in (self.cmin, self.cthresh, self.cfac, self.cnum):
            try:
                ctrl.Bind(wx.EVT_TEXT, changed)
            except Exception:
                pass

    def _set_compact_frame_size(self) -> None:
        """Size this frame compactly without exceeding the parent/main height."""
        try:
            parent_size = self.process_parent.GetSize()
            max_h = max(420, int(parent_size.height))
        except Exception:
            max_h = 720
        try:
            display = wx.Display(self.GetScreen()).GetClientArea()
            max_h = min(max_h, int(display.height))
        except Exception:
            pass
        self.SetSize((720, max_h))

    def _fit_phase_frame(self) -> None:
        if not hasattr(self, 'phaseFrame'):
            return
        try:
            self.phasePanel.Layout()
            self.phaseFrame.GetSizer().Fit(self.phaseFrame)
            width = max(620, int(self.GetSize().width))
            height = max(70, int(self.phaseFrame.GetSize().height))
            self.phaseFrame.SetSize((width, height))
            if self.phaseFrame.IsShown():
                self._position_phase_frame()
        except Exception:
            pass

    def _position_phase_frame(self) -> None:
        """Dock the phase-slider frame exactly beneath the Projections frame.

        The slider window is a visual extension of this frame: its left edge
        follows the Projections frame and its top edge is the first screen row
        immediately below the Projections frame.  Do not independently clamp
        or flip the companion window to the display, because doing so breaks
        that docking relationship.
        """
        try:
            rect = self.GetScreenRect()
            x = int(rect.x)
            y = int(rect.y + rect.height)
            self.phaseFrame.SetPosition((x, y))
        except Exception:
            pass

    def _hide_phase_frame(self, event):
        self.phaseFrame.Hide()
        try:
            self.toolbar.set_sliders_active(False)
        except Exception:
            pass
        if event is not None and event.CanVeto():
            event.Veto()

    def OnSlidersButton(self, active=None):
        if active is None:
            active = not self.phaseFrame.IsShown()
        active = bool(active)
        if active:
            self._fit_phase_frame()
            self._position_phase_frame()
            self.phaseFrame.Show(True)
            self.phaseFrame.Raise()
        else:
            self.phaseFrame.Hide()
        try:
            self.toolbar.set_sliders_active(active)
        except Exception:
            pass

    def _hide_contour_frame(self, event):
        self.contourFrame.Hide()
        if event is not None:
            event.Veto()

    def OnContourButton(self, event=None):
        if not self.contourFrame.IsShown():
            self.contourFrame.Show()
        self.contourFrame.Raise()
        self.cthresh.SetFocus()

    def OnContourEnter(self, event):
        # Synchronise paired threshold/min controls first, then use the same
        # redraw path as the existing projection implementation.
        source = event.GetEventObject()
        if source is self.cmin:
            self._sync_contour_thresh_from_cmin()
        elif source is self.cthresh:
            self._sync_contour_cmin_from_thresh()
        if self._redraw_timer is not None:
            try:
                self._redraw_timer.Stop()
            except Exception:
                pass
        self.refresh_contours(force_full=True)

    def _build_layout(self) -> None:
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.EXPAND)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 4)

        self.main_panel.SetSizer(self.vbox)
        self.main_panel.Layout()
        self.frame_sizer = wx.BoxSizer(wx.VERTICAL)
        self.frame_sizer.Add(self.main_panel, 1, wx.EXPAND)
        self.SetSizer(self.frame_sizer)
        self.Layout()
        try:
            self.canvas.Show()
            self.toolbar.Show(True)
            if hasattr(self.toolbar, 'update'):
                self.toolbar.update()
            if hasattr(self.toolbar, 'Enable'):
                self.toolbar.Enable(True)
        except Exception:
            pass
        wx.CallAfter(self._debug_widget_state, 'post-build-layout')

    def _set_default_values(self) -> None:
        savefile = self._parameter_file_path()
        # parse_* applies these defaults only when a key is absent (or
        # malformed).  Do not replace an explicitly saved zero with a GUI
        # default: the system save file is the source of truth.
        persisted_contours = {
            'cmin': str(parse_float(savefile, 'cmin', default=1e5)),
            'cthresh': str(parse_float(savefile, 'cthresh', default=0.1)),
            'cfac': str(parse_float(savefile, 'cfac', default=1.3)),
            'cnum': str(parse_int(savefile, 'cnum', default=20)),
        }
        if self.state is not None:
            self.state.seed_gui_settings(persisted_contours)
            contour_values = {
                key: str(self.state.resolved_gui_value(key, value))
                for key, value in persisted_contours.items()
            }
        else:
            contour_values = persisted_contours
        self.cmin.SetValue(contour_values['cmin'])
        self.cthresh.SetValue(contour_values['cthresh'])
        self.cfac.SetValue(contour_values['cfac'])
        self.cnum.SetValue(contour_values['cnum'])
        self._contour_baseline_ready = False

        mode = 'matched'
        if self.state is not None:
            try:
                mode = str(self.state.metadata.get('projection_view_mode', 'matched'))
            except Exception:
                mode = 'matched'
        mask_selection = 0.5
        if self.state is not None:
            try:
                mask_selection = float(self.state.metadata.get('projection_mask_selection', 0.5))
            except Exception:
                mask_selection = 0.5
        try:
            self._set_projection_view_mode(mode, redraw=False)
        except Exception:
            self._projection_view_mode_name = self._normalize_projection_view_mode(mode)
        try:
            self.maskStrengthSlider.SetValue(int(round(max(0.0, min(1.0, mask_selection)) * 100.0)))
        except Exception:
            pass
        self._update_mask_slider_widgets()

    def _bind_events(self) -> None:
        for ctrl in (self.cmin, self.cthresh, self.cfac, self.cnum):
            ctrl.Bind(wx.EVT_TEXT, self.OnContourTextChanged)
            ctrl.Bind(wx.EVT_TEXT_ENTER, self.OnContourEnter)
        self.maskStrengthSlider.Bind(wx.EVT_SLIDER, self.OnMaskStrictnessChanged)

    def _projection_tick_values(self, current: float, fine: bool) -> list[float]:
        current = float(current)
        offsets = [-10.0, -5.0, 0.0, 5.0, 10.0] if fine else [-180.0, -90.0, 0.0, 90.0, 180.0]
        return [current + offset for offset in offsets]

    def _format_projection_tick(self, value: float) -> str:
        try:
            value = float(value)
        except Exception:
            return ''
        if abs(value - round(value)) < 1e-9:
            return f'{int(round(value)):+d}°'
        return f'{value:+.2f}°'

    def _normalize_projection_view_mode(self, mode) -> str:
        mode_norm = str(mode).strip().lower()
        if mode_norm in {'proj', 'projection', 'sum'}:
            return 'Proj'
        if mode_norm in {'1d', 'seed'}:
            return '1D'
        if mode_norm in {'mask', 'matched', 'match'}:
            return 'matched'
        return 'matched'

    def _projection_view_mode(self) -> str:
        return self._normalize_projection_view_mode(getattr(self, '_projection_view_mode_name', 'matched'))

    def _projection_view_mode_branch(self) -> str:
        mode = self._projection_view_mode()
        return 'mask' if mode == 'matched' else mode

    def _set_projection_view_mode(self, mode, *, redraw: bool = True) -> None:
        normalized = self._normalize_projection_view_mode(mode)
        self._projection_view_mode_name = normalized
        try:
            if hasattr(self, 'modeLabel'):
                self.modeLabel.SetLabel(f'Mode: {normalized}')
        except Exception:
            pass
        if self.state is not None:
            try:
                self.state.metadata['projection_view_mode'] = normalized
            except Exception:
                pass
        if redraw and not self._initializing:
            self._full_redraw()

    def _mask_selection_value(self) -> float:
        slider = getattr(self, 'maskStrengthSlider', None)
        if slider is None:
            return 0.5
        try:
            return float(np.clip(float(slider.GetValue()) / 100.0, 0.0, 1.0))
        except Exception:
            return 0.5

    def _mask_selection_thresholds(self, selection: float) -> dict[str, float | int]:
        selection = float(np.clip(selection, 0.0, 1.0))
        if selection <= 0.5:
            alpha = selection / 0.5 if selection > 0 else 0.0
            peak_sigma = (1.0 - alpha) * 0.0 + alpha * 4.0
            run_sigma = (1.0 - alpha) * 0.0 + alpha * 2.5
            run_min_points = int(round((1.0 - alpha) * 1 + alpha * 3))
        else:
            alpha = (selection - 0.5) / 0.5
            peak_sigma = (1.0 - alpha) * 4.0 + alpha * 8.0
            run_sigma = (1.0 - alpha) * 2.5 + alpha * 4.0
            run_min_points = int(round((1.0 - alpha) * 3 + alpha * 8))
        return {
            'selection': selection,
            'peak_sigma': float(peak_sigma),
            'run_sigma': float(run_sigma),
            'run_min_points': max(1, int(run_min_points)),
        }

    def _transform_last_axis(self, data, transform):
        arr = np.asarray(data)
        if arr.size == 0:
            return arr
        if arr.ndim == 1:
            return transform(arr)
        try:
            transformed = transform(arr)
            if transformed.shape == arr.shape:
                return transformed
        except Exception:
            pass
        try:
            swapped = np.swapaxes(arr, -1, -1)
            transformed = transform(swapped)
            return np.swapaxes(transformed, -1, -1)
        except Exception:
            return transform(arr)

    def _transform_axis(self, data, axis: int, transform):
        arr = np.asarray(data)
        if arr.size == 0:
            return arr
        axis = int(axis) % arr.ndim
        if arr.ndim == 1 or axis == arr.ndim - 1:
            return transform(arr)
        swapped = np.swapaxes(arr, axis, -1)
        transformed = transform(swapped)
        return np.swapaxes(transformed, axis, -1)

    def _hilbert_trace_ng(self, trace):
        trace = np.asarray(trace)
        if trace.size == 0:
            return trace.astype(complex)
        if np.iscomplexobj(trace):
            return np.asarray(trace, dtype=complex)
        real_trace = np.asarray(trace, dtype=float)
        npts = int(real_trace.shape[-1])
        self._debug_projection(f"_hilbert_trace_ng explicit N={npts} input_shape={real_trace.shape!r} dtype={real_trace.dtype!r}")
        try:
            return np.asarray(ng.process.proc_base.ht(real_trace, N=npts))
        except Exception as exc:
            self._debug_projection(f"_hilbert_trace_ng explicit N failed: {exc!r}")
            return self._hilbert_projection_trace(real_trace)

    def _direct_dimension_seed_trace(self, proj2d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        proj2d = np.asarray(proj2d)
        self._debug_projection(f"_direct_dimension_seed_trace input: shape={proj2d.shape!r} dtype={proj2d.dtype!r} ndims={proj2d.ndim!r}")
        if proj2d.ndim != 2 or proj2d.size == 0:
            self._debug_projection('_direct_dimension_seed_trace early-exit: empty or non-2D input')
            empty = np.asarray([], dtype=float)
            return empty, empty.astype(complex), empty

        # The 2D projection plots direct dimension on axis 0 and indirect on axis 1.
        direct_npts = int(proj2d.shape[0])
        self._debug_projection(f'applying hilbert to direct axis=0 with explicit N={direct_npts}')
        direct_hilbert = self._transform_axis(
            proj2d,
            0,
            lambda arr, n=direct_npts: ng.process.proc_base.ht(arr, N=int(n)),
        )
        self._debug_projection(f"after hilbert: shape={getattr(direct_hilbert, 'shape', None)!r} dtype={getattr(direct_hilbert, 'dtype', None)!r}")
        self._debug_projection('applying inverse FFT to direct axis=0')
        direct_ifft = self._transform_axis(direct_hilbert, 0, ng.process.proc_base.ifft)
        self._debug_projection(f"after ifft: shape={getattr(direct_ifft, 'shape', None)!r} dtype={getattr(direct_ifft, 'dtype', None)!r}")
        display_real = np.asarray(np.real(direct_ifft), dtype=float)
        seed_complex = np.take(direct_ifft, 0, axis=0)
        seed_real = np.asarray(np.real(seed_complex), dtype=float)
        self._debug_projection(
            f"seed/display slices: seed_complex.shape={getattr(seed_complex, 'shape', None)!r} seed_real.shape={seed_real.shape!r} display_real.shape={display_real.shape!r} display_minmax={((float(np.nanmin(display_real)), float(np.nanmax(display_real))) if display_real.size else None)!r}"
        )
        return seed_real, np.asarray(seed_complex, dtype=complex), display_real

    def _longest_true_run(self, mask: np.ndarray) -> int:
        arr = np.asarray(mask, dtype=bool).ravel()
        longest = current = 0
        for flag in arr:
            if flag:
                current += 1
                longest = max(longest, current)
            else:
                current = 0
        return int(longest)

    def _masked_projection_trace(self, proj2d: np.ndarray, indirect_scale: np.ndarray, *, selection_strength: float = 0.5) -> tuple[np.ndarray, np.ndarray, dict]:
        proj2d = np.asarray(proj2d, dtype=float)
        indirect_scale = np.asarray(indirect_scale, dtype=float).ravel()
        self._debug_projection(
            f"_masked_projection_trace input: shape={proj2d.shape!r} dtype={proj2d.dtype!r} indirect_scale.size={indirect_scale.size!r}"
        )
        if proj2d.ndim != 2 or proj2d.size == 0:
            empty = np.asarray([], dtype=float)
            return empty, np.zeros(0, dtype=bool), {'accepted_rows': 0, 'total_rows': 0, 'fallback': 'empty'}

        if indirect_scale.size != proj2d.shape[1]:
            self._debug_projection(
                f"_masked_projection_trace indirect_scale mismatch: scale.size={indirect_scale.size!r} row_len={proj2d.shape[1]!r}; using index scale"
            )
            x = np.arange(proj2d.shape[1], dtype=float)
        else:
            x = indirect_scale

        thresholds = self._mask_selection_thresholds(selection_strength)
        peak_sigma = float(thresholds['peak_sigma'])
        run_sigma = float(thresholds['run_sigma'])
        run_min_points = int(thresholds['run_min_points'])
        self._debug_projection(
            f"_masked_projection_trace selection_strength={selection_strength!r} thresholds={{'peak_sigma': {peak_sigma!r}, 'run_sigma': {run_sigma!r}, 'run_min_points': {run_min_points!r}}}"
        )
        accepted = np.zeros(proj2d.shape[0], dtype=bool)
        debug_rows = []
        row_metrics: list[tuple[int, float, float, float, int, bool]] = []
        peak_only = selection_strength >= 0.999
        if peak_only:
            self._debug_projection('_masked_projection_trace peak-only mode active at selection_strength >= 1.0')
        for row_idx, row in enumerate(proj2d):
            finite = np.isfinite(row) & np.isfinite(x)
            finite_count = int(np.count_nonzero(finite))
            if finite_count < 4:
                debug_rows.append((row_idx, finite_count, None, None, None, None, False, 'too_few_points'))
                continue
            xf = x[finite]
            yf = np.asarray(row[finite], dtype=float)
            try:
                slope, intercept = np.polyfit(xf, yf, 1)
            except Exception as exc:
                debug_rows.append((row_idx, finite_count, None, None, None, None, False, f'polyfit_failed:{exc!r}'))
                continue
            trend = slope * xf + intercept
            resid = yf - trend
            resid_center = float(np.nanmedian(resid))
            centered = resid - resid_center
            mad = float(np.nanmedian(np.abs(centered)))
            sigma = 1.4826 * mad if np.isfinite(mad) else np.nan
            if not np.isfinite(sigma) or sigma <= 0:
                sigma = float(np.nanstd(centered))
            if not np.isfinite(sigma) or sigma <= 0:
                sigma = 1.0
            abs_centered = np.abs(centered)
            peak = float(np.nanmax(abs_centered)) if abs_centered.size else 0.0
            peak_z = peak / sigma if sigma else 0.0
            above = abs_centered >= (run_sigma * sigma)
            longest_run = self._longest_true_run(above)
            area_excess = float(np.nansum(np.clip(abs_centered - run_sigma * sigma, 0.0, None))) / sigma if sigma else 0.0
            if peak_only:
                row_metrics.append((row_idx, peak, peak_z, sigma, longest_run, finite_count))
            else:
                accepted_flag = bool(peak_z >= peak_sigma and longest_run >= run_min_points)
                if accepted_flag:
                    accepted[row_idx] = True
                debug_rows.append((row_idx, finite_count, float(slope), float(intercept), float(sigma), float(peak_z), accepted_flag, longest_run, area_excess))
                if row_idx < 5 or accepted_flag:
                    self._debug_projection(
                        f"mask row {row_idx}: finite={finite_count} slope={slope!r} intercept={intercept!r} sigma={sigma!r} peak_z={peak_z!r} longest_run={longest_run!r} area_excess={area_excess!r} accepted={accepted_flag!r}"
                    )

        if peak_only:
            if row_metrics:
                best = max(row_metrics, key=lambda item: (item[1], item[2], item[3], item[4], -item[0]))
                best_row = int(best[0])
                accepted[best_row] = True
                self._debug_projection(
                    f"_masked_projection_trace peak-only selection -> best_row={best_row!r} row_peak={best[1]!r} peak_z={best[2]!r} sigma={best[3]!r} longest_run={best[4]!r} finite_count={best[5]!r}"
                )
                debug_rows.append((best_row, best[5], None, None, best[3], best[2], True, best[4], best[1]))
            else:
                self._debug_projection('_masked_projection_trace peak-only mode found no eligible rows; returning zeros')

        if accepted.any():
            masked_proj = np.sum(proj2d[accepted], axis=0)
        else:
            masked_proj = np.zeros(proj2d.shape[1], dtype=float)
            self._debug_projection('_masked_projection_trace no slices passed the mask; returning zeros')

        stats = {
            'selection_strength': float(selection_strength),
            'accepted_rows': int(np.count_nonzero(accepted)),
            'total_rows': int(proj2d.shape[0]),
            'peak_sigma': peak_sigma,
            'run_sigma': run_sigma,
            'run_min_points': run_min_points,
            'debug_rows': debug_rows,
        }
        self._debug_projection(
            f"_masked_projection_trace summary: accepted_rows={stats['accepted_rows']!r} total_rows={stats['total_rows']!r} peak_sigma={peak_sigma!r} run_sigma={run_sigma!r} run_min_points={run_min_points!r}"
        )
        return np.asarray(masked_proj, dtype=float), accepted, stats

    def _projection_1d_phase_dic(
        self,
        source_dic,
        projected_1d: np.ndarray,
        *,
        direct_label: str,
        indirect_label: str,
        source_path: str | None = None,
    ) -> dict:
        projected_1d = np.asarray(projected_1d)
        source_dic = dict(source_dic or {})
        snapshot_keys = [
            'FDDIMCOUNT',
            'FDDIMORDER',
            'FDDIMORDER1',
            'FDDIMORDER2',
            'FDDIMORDER3',
            'FDDIMORDER4',
            'FDSIZE',
            'FDSPECNUM',
            'FDF1LABEL',
            'FDF2LABEL',
            'FDF1FTSIZE',
            'FDF2FTSIZE',
            'FDF1AQSIGN',
            'FDF2AQSIGN',
            'FDF1CAR',
            'FDF2CAR',
            'FDF1CENTER',
            'FDF2CENTER',
            'FD1DPHASE',
            'FD2DPHASE',
        ]
        self._debug_projection(
            f"_projection_1d_phase_dic start: source_path={source_path!r} direct_label={direct_label!r} indirect_label={indirect_label!r} projected_1d_shape={projected_1d.shape!r} dtype={projected_1d.dtype!r}"
        )
        self._debug_projection(f"_projection_1d_phase_dic source snapshot: {{k: source_dic.get(k) for k in snapshot_keys}}")

        def _coerce_order(value):
            if value is None:
                return []
            if isinstance(value, (list, tuple, np.ndarray)):
                out = []
                for item in value:
                    try:
                        out.append(float(item))
                    except Exception:
                        pass
                return out
            if isinstance(value, str):
                items = [part for part in re.split(r'[\s,]+', value.strip()) if part]
                out = []
                for item in items:
                    try:
                        out.append(float(item))
                    except Exception:
                        pass
                return out
            return []

        phase_dic = dict(source_dic)
        original_order = _coerce_order(phase_dic.get('FDDIMORDER'))
        if not original_order:
            original_order = []
            for idx in range(1, 5):
                try:
                    val = phase_dic.get(f'FDDIMORDER{idx}')
                    if val is not None:
                        original_order.append(float(val))
                except Exception:
                    pass
        original_order = [float(x) for x in original_order if float(x) != 4.0]
        if not original_order:
            original_order = [1.0, 2.0]
            self._debug_projection('_projection_1d_phase_dic missing FDDIMORDER; using fallback order [1.0, 2.0]')

        label_to_dim = {}
        for idx in range(1, 5):
            label = phase_dic.get(f'FDF{idx}LABEL')
            if label is not None:
                label_to_dim[str(label).strip()] = float(idx)

        keep_dim = label_to_dim.get(str(indirect_label).strip())
        if keep_dim is None:
            keep_dim = label_to_dim.get(str(direct_label).strip())
        if keep_dim is None and original_order:
            keep_dim = float(original_order[0])
            self._debug_projection(f"_projection_1d_phase_dic could not map labels to a dimension; defaulting keep_dim={keep_dim!r}")

        remaining = [float(x) for x in original_order if float(x) != float(keep_dim)]
        for candidate in (1.0, 2.0, 3.0):
            if candidate != float(keep_dim) and candidate not in remaining:
                remaining.append(candidate)
            if len(remaining) >= 2:
                break
        while len(remaining) < 2:
            remaining.append(4.0)

        projected_size = int(projected_1d.shape[-1] if projected_1d.ndim else projected_1d.size)
        old_size = phase_dic.get('FDSIZE')
        old_specnum = phase_dic.get('FDSPECNUM')
        old_order = phase_dic.get('FDDIMORDER')

        phase_dic['FDDIMCOUNT'] = '1'
        new_order = (float(keep_dim), float(remaining[0]), float(remaining[1]), 4.0)
        phase_dic['FDDIMORDER'] = new_order
        phase_dic['FDDIMORDER1'] = new_order[0]
        phase_dic['FDDIMORDER2'] = new_order[1]
        phase_dic['FDDIMORDER3'] = new_order[2]
        phase_dic['FDDIMORDER4'] = new_order[3]
        phase_dic['FDSIZE'] = projected_size
        phase_dic['FDSPECNUM'] = 0

        keep_key = f'FDF{int(float(keep_dim))}FTSIZE'
        if phase_dic.get(keep_key) is not None:
            phase_dic[keep_key] = projected_size

        if phase_dic.get(f'FDF{int(float(keep_dim))}LABEL') is None:
            phase_dic[f'FDF{int(float(keep_dim))}LABEL'] = str(indirect_label)

        updated_snapshot = {k: phase_dic.get(k) for k in snapshot_keys}
        self._debug_projection(
            f"_projection_1d_phase_dic updated: keep_dim={keep_dim!r} original_order={original_order!r} new_order={new_order!r} projected_size={projected_size!r} old_size={old_size!r} old_specnum={old_specnum!r} old_order={old_order!r} updated_snapshot={updated_snapshot!r}"
        )
        return phase_dic

    def _seed_overlay_line(self, ax2d, indirect_scale: np.ndarray, seed_trace: np.ndarray, direct_scale: np.ndarray):
        indirect_scale = np.asarray(indirect_scale)
        seed_trace = np.asarray(seed_trace, dtype=float)
        direct_scale = np.asarray(direct_scale, dtype=float)
        if indirect_scale.size == 0 or seed_trace.size == 0 or direct_scale.size == 0:
            return None
        if not np.isfinite(seed_trace).any():
            return None

        span = float(np.nanmax(direct_scale) - np.nanmin(direct_scale))
        if not np.isfinite(span) or span == 0.0:
            span = 1.0
        max_abs = float(np.nanmax(np.abs(seed_trace))) if np.isfinite(np.nanmax(np.abs(seed_trace))) else 0.0
        if max_abs == 0.0:
            max_abs = 1.0
        normalized = seed_trace / max_abs
        baseline = float(np.nanmax(direct_scale))
        y = baseline - normalized * (0.10 * span)
        try:
            return ax2d.plot(indirect_scale, y, color='0.2', linewidth=0.9, alpha=0.9)[0]
        except Exception:
            return None

    def _update_mask_slider_widgets(self) -> None:
        slider = getattr(self, 'maskStrengthSlider', None)
        label = getattr(self, 'maskStrengthValue', None)
        if slider is None or label is None:
            return
        try:
            value = self._mask_selection_value()
            label.SetLabel(f'{value:.2f}')
            label.Wrap(-1)
        except Exception:
            pass

    def _update_phase_slider_value_widgets(self, control: dict) -> None:
        for which in ('p0', 'p1'):
            slider = control.get(f'{which}_slider')
            value_label = control.get(f'{which}_value_label')
            mode_btn = control.get(f'{which}_mode_btn')
            if slider is None:
                continue
            current = self._projection_slider_value(slider, mode_btn)
            if value_label is not None:
                try:
                    if abs(current - round(current)) > 1e-9:
                        value_label.SetLabel(f'{current:+.2f}°')
                    else:
                        value_label.SetLabel(f'{int(round(current)):+d}°')
                except Exception:
                    pass

    def _update_phase_slider_tick_widgets(self, control: dict) -> None:
        for which in ('p0', 'p1'):
            slider = control.get(f'{which}_slider')
            tick_labels = control.get(f'{which}_tick_labels')
            mode_btn = control.get(f'{which}_mode_btn')
            if slider is None or not tick_labels:
                continue
            current = self._projection_slider_value(slider, mode_btn)
            fine = bool(mode_btn.GetValue()) if mode_btn is not None else False
            tick_values = self._projection_tick_values(current, fine)
            control[f'{which}_tick_values'] = tick_values
            for lab, value in zip(tick_labels, tick_values):
                try:
                    lab.SetLabel(self._format_projection_tick(value))
                except Exception:
                    pass

    def _projection_slider_scale(self, button=None):
        return 100

    def _projection_slider_value_for_mode(self, slider, fine=False):
        try:
            raw = float(slider.GetValue())
        except Exception:
            raw = 0.0
        scale = 100.0
        return raw / scale

    def _projection_slider_value(self, slider, button=None):
        try:
            raw = float(slider.GetValue())
        except Exception:
            raw = 0.0
        scale = self._projection_slider_scale(button)
        return raw / float(scale)

    def _set_projection_slider_value(self, slider, button, value):
        try:
            scale = self._projection_slider_scale(button)
            slider.SetValue(int(round(float(value) * float(scale))))
        except Exception:
            self._debug_projection('Failed to set projection slider value')

    def _set_projection_slider_mode(self, slider, button, fine=False, current_value=None):
        try:
            if slider is None or button is None:
                return
            if current_value is None:
                current_value = self._projection_slider_value_for_mode(slider, fine=not fine)
            cur = float(current_value)
            span = 10.0 if fine else 180.0
            scale = 100
            slider.SetRange(int(round((cur - span) * scale)), int(round((cur + span) * scale)))
            try:
                slider.SetLineSize(1)
                slider.SetPageSize(1 if fine else 10)
            except Exception:
                pass
            button.SetLabel('F' if fine else 'C')
            self._set_projection_slider_value(slider, button, cur)
        except Exception:
            self._debug_projection('Failed to update projection slider mode')

    def _update_phase_control_layout(self) -> None:
        for control in self._phase_controls:
            for key in ('p0_slider', 'p1_slider'):
                slider = control.get(key)
                if slider is None:
                    continue
                try:
                    slider.SetMinSize((-1, 28))
                except Exception:
                    pass
            for key in ('p0_mode_btn', 'p1_mode_btn'):
                button = control.get(key)
                if button is None:
                    continue
                try:
                    button.SetMinSize((26, 22))
                except Exception:
                    pass
            self._update_phase_slider_value_widgets(control)
        try:
            self.phasePanel.Layout()
            self.phasePanel.Layout()
            self._fit_phase_frame()
        except Exception:
            pass

    def _clear_phase_controls(self) -> None:
        try:
            self.phasePanel.DestroyChildren()
        except Exception:
            pass
        try:
            self.phasePanelSizer.Clear(False)
        except Exception:
            pass
        self._phase_controls = []

    def _projection_phase_state(self, label: str) -> dict[str, float]:
        baseline = self._processing_phase_defaults(label)
        if self.state is None:
            return {'p0': round(float(baseline['p0']), 2),
                    'p1': round(float(baseline['p1']), 2)}
        return self.state.projection_phase(
            label, p0=baseline['p0'], p1=baseline['p1']
        )

    def _update_projection_phase_state(self, label: str, *, p0: float | None = None, p1: float | None = None) -> None:
        if self.state is not None:
            self.state.update_projection_phase(label, p0=p0, p1=p1)

    def _clear_projection_phase_state(self) -> None:
        if self.state is not None:
            self.state.clear_projection_phase_preview()

    def _projection_phase_key(self, label: str, which: str) -> str | None:
        labels = self._spectral_labels()
        try:
            index = labels.index(label)
        except ValueError:
            normalized = [str(x).strip() for x in labels]
            try:
                index = normalized.index(str(label).strip())
            except ValueError:
                return None
        if index == 0:
            return which
        return f'{which}_{index}'

    def _format_projection_phase_value(self, value) -> str:
        try:
            return f'{float(value):.2f}'
        except Exception:
            return str(value)

    def _save_projection_phase_values(self) -> None:
        """Persist the projection sliders before an automatic re-process.

        In a normal 2D spectrum this window always represents the indirect
        dimension.  Its P0/P1 sliders therefore map explicitly to p0_1/p1_1.
        Do not infer the parameter-file index from nucleus-label ordering: the
        display label can differ from the label order in decon.par.
        """
        # Use the exact same canonical save-file location as ProcessFrame.
        # ProjectionsFrame is created with ProcessFrame as its parent.
        savefile = None
        try:
            parent_dir_box = None
            parent_par_name = getattr(self.process_parent, 'deconParFile', None)
            if parent_dir_box is not None and parent_par_name:
                savefile = os.path.join(parent_dir_box.GetValue(), str(parent_par_name))
        except Exception as exc:
            self._debug_projection(f'Could not resolve canonical ProcessFrame savefile: {exc!r}')

        if not savefile:
            savefile = self._parameter_file_path()

        updates: dict[str, str] = {}

        # 2D projection phasing is always the indirect-dimension phasing.
        if self._spectral_dimension_count() == 2:
            control = self._phase_controls[0] if self._phase_controls else None
            if control is not None:
                label = control.get('label') or ''
                p0_slider = control.get('p0_slider')
                p1_slider = control.get('p1_slider')
                p0_value = (
                    self._projection_slider_value(p0_slider, control.get('p0_mode_btn'))
                    if p0_slider is not None else None
                )
                p1_value = (
                    self._projection_slider_value(p1_slider, control.get('p1_mode_btn'))
                    if p1_slider is not None else None
                )

                if p0_value is not None:
                    updates['p0_1'] = self._format_projection_phase_value(p0_value)
                if p1_value is not None:
                    updates['p1_1'] = self._format_projection_phase_value(p1_value)
                if p0_value is not None or p1_value is not None:
                    self._update_projection_phase_state(label, p0=p0_value, p1=p1_value)

                self._debug_projection(
                    f'2D Re-process phase save: label={label!r} '
                    f'p0_1={updates.get("p0_1")!r} p1_1={updates.get("p1_1")!r} '
                    f'savefile={savefile!r}'
                )

        else:
            # Preserve the existing label-indexed behaviour for 3D/4D data.
            for control in self._phase_controls:
                label = control.get('label')
                if not label:
                    continue
                p0_slider = control.get('p0_slider')
                p1_slider = control.get('p1_slider')
                p0_value = (
                    self._projection_slider_value(p0_slider, control.get('p0_mode_btn'))
                    if p0_slider is not None else None
                )
                p1_value = (
                    self._projection_slider_value(p1_slider, control.get('p1_mode_btn'))
                    if p1_slider is not None else None
                )
                key0 = self._projection_phase_key(label, 'p0')
                key1 = self._projection_phase_key(label, 'p1')
                if key0 is not None and p0_value is not None:
                    updates[key0] = self._format_projection_phase_value(p0_value)
                if key1 is not None and p1_value is not None:
                    updates[key1] = self._format_projection_phase_value(p1_value)
                if p0_value is not None or p1_value is not None:
                    self._update_projection_phase_state(label, p0=p0_value, p1=p1_value)

        if not updates:
            self._debug_projection(
                f'Projection Re-process did not produce phase updates; savefile={savefile!r}'
            )
            return

        # Promote the accepted projection phase into the live processing
        # controls.  Slider movement itself remains temporary; only Re-process
        # calls this method.  The central Process save then persists one atomic
        # snapshot before script generation.
        processing = getattr(self.process_parent, 'processing_frame', None)
        for key, value in updates.items():
            ctrl = getattr(processing, key, None) if processing is not None else None
            if ctrl is not None:
                try:
                    ctrl.SetValue(str(value))
                except Exception:
                    pass
        if self.state is not None:
            self.state.promote_projection_phase(updates)
        self._promoted_phase_updates = dict(updates)

    def _refresh_projection_after_reprocess(self) -> None:
        self._debug_projection('_refresh_projection_after_reprocess start')
        self._clear_projection_phase_state()
        self._contour_baseline_ready = False

        # Rebuild controls from the system save file so the sliders reflect the
        # values that were actually persisted and then used for processing.
        savefile = self._parameter_file_path()
        if self._spectral_dimension_count() == 2:
            p0 = parse_float(savefile, 'p0_1', 0.0)
            p1 = parse_float(savefile, 'p1_1', 0.0)
            self._debug_projection(
                f'_refresh_projection_after_reprocess saved 2D phase values: p0_1={p0!r} p1_1={p1!r} savefile={savefile!r}'
            )

        try:
            self.load_projections()
            self._debug_projection('_refresh_projection_after_reprocess load_projections completed')
        except Exception as exc:
            self._debug_projection(f'Could not reload projections after re-process: {exc!r}')
            import traceback
            self._debug_projection(traceback.format_exc())
        try:
            self.main_panel.Layout()
        except Exception:
            pass
        self._debug_widget_state('post-reprocess-refresh')

    def _run_silent_reprocess(self) -> None:
        self._save_projection_phase_values()

        def _finish(*_args, **_kwargs):
            self._debug_projection('_run_silent_reprocess completion callback received')
            wx.CallAfter(self._refresh_projection_after_reprocess)

        runner = getattr(self.process_parent, '_run_processing_auto', None)
        if runner is None:
            raise AttributeError('Process frame has no _run_processing_auto method')
        runner(on_finish=_finish)

    def _on_projection_phase_scroll(self, event, label: str, which: str) -> None:
        try:
            slider = event.GetEventObject()
        except Exception:
            slider = None
        current_value = None
        mode_btn = None
        for control in self._phase_controls:
            if control.get('label') == label:
                mode_btn = control.get(f'{which}_mode_btn')
                if slider is not None:
                    current_value = self._projection_slider_value(slider, mode_btn)
                break
        if current_value is None:
            try:
                current_value = float(slider.GetValue()) if slider is not None else 0.0
            except Exception:
                current_value = 0.0
        if which == 'p0':
            self._update_projection_phase_state(label, p0=current_value)
        else:
            self._update_projection_phase_state(label, p1=current_value)

        for control in self._phase_controls:
            if control.get('label') == label:
                self._update_phase_slider_value_widgets(control)
                break
        self._update_projection_phase_plot(label)
        try:
            self.phasePanel.Layout()
        except Exception:
            pass
        event.Skip()

    def _on_projection_phase_mode_toggle(self, event, label: str, which: str) -> None:
        for control in self._phase_controls:
            if control.get('label') == label:
                button = control.get(f'{which}_mode_btn')
                slider = control.get(f'{which}_slider')
                if slider is not None and button is not None:
                    new_fine = bool(button.GetValue())
                    current = self._projection_slider_value_for_mode(slider, fine=not new_fine)
                    self._set_projection_slider_mode(slider, button, new_fine, current_value=current)
                    current_value = self._projection_slider_value(slider, button)
                    if which == 'p0':
                        self._update_projection_phase_state(label, p0=current_value)
                    else:
                        self._update_projection_phase_state(label, p1=current_value)
                    self._update_phase_slider_tick_widgets(control)
                    self._update_phase_slider_value_widgets(control)
                    self._update_projection_phase_plot(label)
                    try:
                        self.phasePanel.Layout()
                    except Exception:
                        pass
                break
        event.Skip()

    def _hilbert_projection_trace(self, projection_1d):
        projection_1d = np.asarray(projection_1d)
        self._debug_projection(f"_hilbert_projection_trace input: shape={projection_1d.shape!r} dtype={projection_1d.dtype!r} complex={np.iscomplexobj(projection_1d)}")
        if projection_1d.size == 0:
            self._debug_projection('_hilbert_projection_trace early-exit: empty input')
            return projection_1d.astype(complex)
        if np.iscomplexobj(projection_1d):
            self._debug_projection('_hilbert_projection_trace input already complex; returning copy')
            return np.asarray(projection_1d, dtype=complex)

        real_trace = np.asarray(projection_1d, dtype=float)
        if scipy_hilbert is not None:
            try:
                out = scipy_hilbert(real_trace)
                self._debug_projection(f'_hilbert_projection_trace used scipy.signal.hilbert -> shape={out.shape!r} dtype={out.dtype!r}')
                return out
            except Exception:
                pass

        # Fallback analytic signal construction for environments where SciPy is unavailable.
        n = real_trace.shape[0]
        spectrum = np.fft.fft(real_trace)
        h = np.zeros(n)
        if n % 2 == 0:
            h[0] = 1
            h[n // 2] = 1
            h[1:n // 2] = 2
        else:
            h[0] = 1
            h[1:(n + 1) // 2] = 2
        return np.fft.ifft(spectrum * h)

    def _phase_projection_trace(
        self,
        projection_1d,
        label: str,
        p0: int | float = 0,
        p1: int | float = 0,
        dic: dict | None = None,
    ):
        projection_1d = np.asarray(projection_1d)
        self._debug_projection(f"_phase_projection_trace input: label={label!r} shape={projection_1d.shape!r} dtype={projection_1d.dtype!r} p0={p0!r} p1={p1!r}")
        if projection_1d.size == 0:
            self._debug_projection('_phase_projection_trace early-exit: empty input')
            empty = projection_1d.astype(complex)
            return empty, empty

        baseline = self._processing_phase_defaults(label)
        complex_trace = np.asarray(projection_1d, dtype=complex)
        phase_dic = self._projection_1d_phase_dic(dic, projection_1d, direct_label=self._direct_label(), indirect_label=label, source_path='phase-input')

        pipe_proc = getattr(ng, 'pipe_proc', None)
        if pipe_proc is None:
            try:
                pipe_proc = ng.process.pipe_proc
            except Exception:
                pipe_proc = None

        if pipe_proc is not None:
            phase_snapshot = {k: phase_dic.get(k) for k in ['FDDIMCOUNT', 'FDDIMORDER', 'FDDIMORDER1', 'FDDIMORDER2', 'FDDIMORDER3', 'FDDIMORDER4', 'FDSIZE', 'FDSPECNUM', 'FDF1LABEL', 'FDF2LABEL', 'FDF1FTSIZE', 'FDF2FTSIZE']}
            self._debug_projection(f'_phase_projection_trace using pipe_proc.ps with baseline={baseline!r} phase_snapshot={phase_snapshot!r}')
            try:
                dic_out, unphased = pipe_proc.ps(
                    phase_dic,
                    complex_trace,
                    p0=float(baseline['p0']),
                    p1=float(baseline['p1']),
                    inv=True,
                    noup=True,
                )
                self._debug_projection(f"_phase_projection_trace pipe_proc.ps first pass returned dic_keys={sorted(dic_out.keys())[:25]!r} unphased.shape={getattr(unphased, 'shape', None)!r} dtype={getattr(unphased, 'dtype', None)!r}")
                dic_out2, phased = pipe_proc.ps(
                    dic_out,
                    unphased,
                    p0=float(p0),
                    p1=float(p1),
                    inv=False,
                    noup=True,
                )
                self._debug_projection(f"_phase_projection_trace pipe_proc.ps result: dic_keys={sorted(dic_out2.keys())[:25]!r} phased.shape={getattr(phased, 'shape', None)!r} phased.dtype={getattr(phased, 'dtype', None)!r}")
                return complex_trace, phased
            except Exception as exc:
                self._debug_projection(f'_phase_projection_trace pipe_proc.ps failed: {exc!r}')
                self._debug_projection(f"_phase_projection_trace failed phase_dic snapshot: {phase_snapshot!r}")

        # Fallback for environments where pipe_proc is unavailable.
        self._debug_projection('_phase_projection_trace using proc_base.ps fallback')
        unphased = ng.process.proc_base.ps(complex_trace, p0=float(baseline['p0']), p1=float(baseline['p1']), inv=True)
        phased = ng.process.proc_base.ps(unphased, p0=float(p0), p1=float(p1), inv=False)
        return complex_trace, phased

    def _save_projection_cache(self, cache: dict) -> None:
        if self.state is None:
            return

        payload = {
            'dic': cache.get('dic'),
            'data': cache.get('proj2d'),
            'labb': (cache.get('direct_label'), cache.get('indirect_label')),
            'source': cache.get('source'),
            'projection_mode': cache.get('projection_mode'),
            'projection_1d_raw': cache.get('proj1d_raw'),
            'projection_1d_hilbert': cache.get('proj1d_hilbert'),
            'projection_1d_phased': cache.get('proj1d_phased'),
            'phase_p0': cache.get('phase_p0'),
            'phase_p1': cache.get('phase_p1'),
            'direct_scale': cache.get('direct_scale'),
            'indirect_scale': cache.get('indirect_scale'),
            'projection_mask_selection': cache.get('proj1d_mask_selection'),
        }
        save_key = ('projection_phase', cache.get('direct_label'), cache.get('indirect_label'), cache.get('source'))

        try:
            if hasattr(self.state, 'save_projection'):
                self.state.save_projection(save_key, **payload)
            else:
                self.state.projections[save_key] = {k: v for k, v in payload.items() if v is not None}
        except Exception:
            pass

    def _update_projection_phase_plot(self, label: str) -> None:
        updated = False
        for cache in self._projection_cache:
            if cache.get('indirect_label') != label:
                continue

            line = cache.get('line1d')
            ax1d = cache.get('ax1d')
            proj1d_raw = cache.get('proj1d_raw')
            proj1d_hilbert = cache.get('proj1d_hilbert')
            indirect_scale = cache.get('indirect_scale')
            if line is None or ax1d is None or proj1d_raw is None or indirect_scale is None:
                continue

            if proj1d_hilbert is None:
                proj1d_hilbert = self._hilbert_trace_ng(proj1d_raw)

            phase_state = self._projection_phase_state(label)
            hilbert, phased = self._phase_projection_trace(proj1d_hilbert, label, phase_state['p0'], phase_state['p1'])
            y = np.real(phased)

            cache['proj1d_hilbert'] = hilbert
            cache['proj1d_phased'] = phased
            cache['phase_p0'] = phase_state['p0']
            cache['phase_p1'] = phase_state['p1']
            line.set_ydata(y)

            if y.size > 0:
                pmin = float(np.nanmin(y))
                pmax = float(np.nanmax(y))
                if pmin == pmax:
                    pad = abs(pmax) * 0.05 if pmax != 0 else 1.0
                else:
                    pad = max(abs(pmin), abs(pmax)) * 0.05
                ax1d.set_ylim(pmin - pad, pmax + pad)

            self._save_projection_cache(cache)
            updated = True

        if updated:
            try:
                self.canvas.draw_idle()
            except Exception:
                pass

    def _rebuild_phase_controls(self, entries: list[tuple[str, str]]) -> None:
        self._clear_phase_controls()
        if not entries:
            try:
                self.phaseFrame.Hide()
                self.toolbar.set_sliders_active(False)
            except Exception:
                pass
            return

        self._debug_projection(f'_rebuild_phase_controls entries={entries!r}')
        for idx, (_path, title) in enumerate(entries):
            indirect_label = title.split('.', 1)[-1] if title else f'Proj {idx + 1}'
            phase_state = self._projection_phase_state(indirect_label)
            self._debug_projection(
                f'_rebuild_phase_controls idx={idx!r} title={title!r} indirect_label={indirect_label!r} phase_state={phase_state!r}'
            )

            group = wx.Panel(self.phasePanel)
            row = wx.BoxSizer(wx.HORIZONTAL)
            header = wx.StaticText(group, label=indirect_label)
            header_font = header.GetFont()
            try:
                header_font = header_font.Bold()
            except Exception:
                pass
            header.SetFont(header_font)
            header.SetMinSize((55, -1))
            value_hint = wx.StaticText(group, label='')
            value_hint.Hide()
            row.Add(header, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)

            def build_slider(which: str, start_value: float):
                start_value = float(start_value)
                start_fine = abs(start_value - round(start_value)) > 1e-9
                span = 10.0 if start_fine else 180.0
                scale = 100
                label_ctrl = wx.StaticText(group, label=which.upper())
                slider = wx.Slider(
                    group, value=int(round(start_value * scale)),
                    minValue=int(round((start_value - span) * scale)),
                    maxValue=int(round((start_value + span) * scale)),
                    size=(190, 24), style=wx.SL_HORIZONTAL,
                )
                mode_btn = wx.ToggleButton(group, label='F' if start_fine else 'C', size=(26, 22))
                mode_btn.SetValue(start_fine)
                self._set_projection_slider_mode(slider, mode_btn, fine=start_fine, current_value=start_value)
                slider.Bind(wx.EVT_SLIDER, partial(self._on_projection_phase_scroll, label=indirect_label, which=which.lower()))
                mode_btn.Bind(wx.EVT_TOGGLEBUTTON, partial(self._on_projection_phase_mode_toggle, label=indirect_label, which=which.lower()))
                value_label = wx.StaticText(group, label='', size=(62, -1))
                row.Add(label_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
                row.Add(slider, 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
                row.Add(mode_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
                row.Add(value_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
                return slider, mode_btn, value_label, start_fine, start_value

            p0_slider, p0_mode_btn, p0_value_label, p0_fine, p0_start_value = build_slider('P0', phase_state['p0'])
            p1_slider, p1_mode_btn, p1_value_label, p1_fine, p1_start_value = build_slider('P1', phase_state['p1'])
            group.SetSizer(row)
            self.phasePanelSizer.Add(group, 0, wx.EXPAND | wx.BOTTOM, 2)
            self._phase_controls.append({
                'label': indirect_label, 'panel': group, 'header': header, 'header_value': value_hint,
                'p0_slider': p0_slider, 'p0_mode_btn': p0_mode_btn,
                'p1_slider': p1_slider, 'p1_mode_btn': p1_mode_btn,
                'p0_value_label': p0_value_label, 'p1_value_label': p1_value_label,
                # Tick widgets are intentionally omitted in the compact companion
                # window; range/mode/value semantics are unchanged.
                'p0_tick_labels': [], 'p1_tick_labels': [],
                'p0_tick_values': self._projection_tick_values(p0_start_value, p0_fine),
                'p1_tick_values': self._projection_tick_values(p1_start_value, p1_fine),
            })
            self._update_phase_slider_value_widgets(self._phase_controls[-1])

        self._update_phase_control_layout()
        self._fit_phase_frame()

    def _clear_projection_cache(self) -> None:
        self._projection_cache = []
        self._axis_backgrounds = {}
        self._background_ready = False

    # ------------------------------------------------------------------
    def _plot_projection_data(self, ax2d, ax1d, path: str, spec, dic, data) -> dict:
        direct_label = spec['direct']
        indirect_label = spec['indirect']
        view_mode_display = self._projection_view_mode()
        view_mode = self._projection_view_mode_branch()
        self._debug_projection(
            f"_plot_projection_data start path={path!r} spec={spec!r} view_mode_display={view_mode_display!r} canonical_view_mode={view_mode!r} direct_label={direct_label!r} indirect_label={indirect_label!r}"
        )

        self._debug_projection(
            f"loaded projection payload -> data_shape={getattr(data, 'shape', None)!r} dtype={getattr(data, 'dtype', None)!r} dic_keys={sorted(dic.keys())[:30]!r}"
        )
        data = np.asarray(data)
        if data.ndim != 2:
            raise ValueError(f'Expected 2D projection data, got {data.ndim}D from {path}')

        special_case = self._is_special_2d_projection_case()
        if special_case:
            self._debug_projection('_plot_projection_data special 2D raw-spectrum layout path active')
            proj2d, direct_scale, indirect_scale, layout_info = self._special_2d_projection_layout(
                dic,
                data,
                direct_label=direct_label,
                indirect_label=indirect_label,
                source_path=path,
            )
            self._debug_projection(
                f"special 2D layout resolved: layout_info={layout_info!r} direct_first_last={(float(direct_scale[0]), float(direct_scale[-1])) if direct_scale.size else None!r} indirect_first_last={(float(indirect_scale[0]), float(indirect_scale[-1])) if indirect_scale.size else None!r}"
            )
        else:
            # The saved projection files use the direct dimension as the second
            # plotted axis. Build the unit converters in that order so the contour
            # axes and the 1D trace follow the intended indirect (X) / direct (Y)
            # layout.
            direct_uc = ng.pipe.make_uc(dic, data, dim=1)
            indirect_uc = ng.pipe.make_uc(dic, data, dim=0)
            direct_scale = np.asarray(direct_uc.ppm_scale())
            indirect_scale = np.asarray(indirect_uc.ppm_scale())
            self._debug_projection(
                f"axis scales: direct_scale.size={direct_scale.size!r} indirect_scale.size={indirect_scale.size!r} direct_first_last={(float(direct_scale[0]), float(direct_scale[-1])) if direct_scale.size else None!r} indirect_first_last={(float(indirect_scale[0]), float(indirect_scale[-1])) if indirect_scale.size else None!r}"
            )

            expected_shape = (direct_scale.size, indirect_scale.size)
            self._debug_projection(f'shape check: data.shape={data.shape!r} expected_shape={expected_shape!r} transpose_shape={data.T.shape!r}')
            if data.shape == expected_shape:
                proj2d = data
                self._debug_projection('using data as-is for proj2d')
            elif data.T.shape == expected_shape:
                proj2d = data.T
                self._debug_projection('using transposed data for proj2d')
            else:
                self._debug_projection(f'shape mismatch path={path!r} data.shape={data.shape!r} expected_shape={expected_shape!r} transpose_shape={data.T.shape!r}')
                raise ValueError(
                    f'Projection array shape {data.shape} does not match expected '
                    f'{expected_shape} or transpose for {path}'
                )

        levels = self._levels()
        self._debug_projection(f'levels: count={levels.size} min={float(np.nanmin(levels)) if levels.size else None!r} max={float(np.nanmax(levels)) if levels.size else None!r}')
        mask_selection = self._mask_selection_value()
        self._debug_projection(f'mask selection slider value={mask_selection!r}')

        proj2d_display = proj2d
        mask_stats = None

        contour_x_scale = indirect_scale
        contour_y_scale = direct_scale
        contour_x_label = indirect_label
        contour_y_label = direct_label
        if special_case and isinstance(locals().get('layout_info'), dict) and layout_info.get('raw_layout'):
            if direct_scale.size == proj2d_display.shape[1]:
                contour_x_scale = direct_scale
                contour_x_label = direct_label
                contour_y_scale = indirect_scale
                contour_y_label = indirect_label
            elif indirect_scale.size == proj2d_display.shape[1]:
                contour_x_scale = indirect_scale
                contour_x_label = indirect_label
                contour_y_scale = direct_scale
                contour_y_label = direct_label
            self._debug_projection(
                f"special raw contour axis selection: proj2d_display.shape={proj2d_display.shape!r} contour_x_scale.size={getattr(contour_x_scale, 'size', None)!r} contour_y_scale.size={getattr(contour_y_scale, 'size', None)!r} contour_x_label={contour_x_label!r} contour_y_label={contour_y_label!r}"
            )

        if special_case and isinstance(locals().get('layout_info'), dict) and layout_info.get('raw_layout') and layout_info.get('transpose_for_contour'):
            trace_axis = 1 if contour_x_scale.size == proj2d.shape[1] else 0
            self._debug_projection(f'special raw 1D projection axis selection: trace_axis={trace_axis!r} (using the axis orthogonal to the contour x-axis so the indirect-dimension projection remains the bottom panel)')
        else:
            trace_axis = 0
        proj1d_proj_raw = np.sum(proj2d, axis=trace_axis)
        self._debug_projection(f"Proj sum result: proj1d_proj_raw.shape={proj1d_proj_raw.shape!r} dtype={proj1d_proj_raw.dtype!r} trace_axis={trace_axis!r} proj2d.shape={proj2d.shape!r} contour_x_scale.size={getattr(contour_x_scale, 'size', None)!r}")
        proj1d_raw = proj1d_proj_raw
        proj1d_seed_complex = None
        if view_mode == '1D':
            self._debug_projection('entering 1D seed path')
            try:
                proj1d_raw, proj1d_seed_complex, proj2d_display = self._direct_dimension_seed_trace(proj2d)
                self._debug_projection(f"1D seed result: proj1d_raw.shape={proj1d_raw.shape!r} proj1d_raw.dtype={proj1d_raw.dtype!r} seed_complex.shape={getattr(proj1d_seed_complex, 'shape', None)!r} seed_complex.dtype={getattr(proj1d_seed_complex, 'dtype', None)!r} proj2d_display.shape={getattr(proj2d_display, 'shape', None)!r} proj2d_display.dtype={getattr(proj2d_display, 'dtype', None)!r}")
            except Exception as exc:
                self._debug_projection(f'1D seed construction failed for {path!r}: {exc!r}')
                raise
        elif view_mode == 'mask':
            self._debug_projection('entering matched/mask projection path')
            try:
                mask_proj2d = proj2d
                mask_scale = contour_x_scale
                if special_case and isinstance(locals().get('layout_info'), dict) and layout_info.get('raw_layout'):
                    # Raw 2D spectra store the indirect axis in rows and the direct axis in columns.
                    # The matched/mask projection helper is row-oriented, so transpose only for the
                    # analysis step to keep the 1D result aligned with the indirect (15N) axis while
                    # leaving the displayed 2D contour untouched.
                    mask_proj2d = proj2d.T
                    mask_scale = direct_scale
                    self._debug_projection(
                        f"special raw mask analysis transpose applied: mask_proj2d.shape={getattr(mask_proj2d, 'shape', None)!r} mask_scale.size={getattr(mask_scale, 'size', None)!r} direct_scale.size={getattr(direct_scale, 'size', None)!r} indirect_scale.size={getattr(indirect_scale, 'size', None)!r}"
                    )
                proj1d_raw, accepted_rows, mask_stats = self._masked_projection_trace(mask_proj2d, mask_scale, selection_strength=mask_selection)
                proj1d_seed_complex = None
                self._debug_projection(
                    f"mask projection result: proj1d_raw.shape={proj1d_raw.shape!r} proj1d_raw.dtype={proj1d_raw.dtype!r} accepted_rows={int(np.count_nonzero(accepted_rows))!r}/{int(accepted_rows.size)!r} display_uses_raw_2d=True contour_x_scale.size={getattr(contour_x_scale, 'size', None)!r}"
                )
            except Exception as exc:
                self._debug_projection(f'mask projection failed for {path!r}: {exc!r}')
                raise
        self._debug_projection(
            f"2D display source: view_mode={view_mode!r} display_shape={getattr(proj2d_display, 'shape', None)!r} display_dtype={getattr(proj2d_display, 'dtype', None)!r} display_minmax={((float(np.nanmin(proj2d_display)), float(np.nanmax(proj2d_display))) if np.size(proj2d_display) else None)!r}"
        )
        if view_mode == 'mask' and mask_stats is not None:
            self._debug_projection(f"mask stats summary: accepted_rows={mask_stats.get('accepted_rows')!r} total_rows={mask_stats.get('total_rows')!r} peak_sigma={mask_stats.get('peak_sigma')!r} run_sigma={mask_stats.get('run_sigma')!r} run_min_points={mask_stats.get('run_min_points')!r}")
        contour_data = proj2d_display
        if special_case and isinstance(layout_info, dict) and layout_info.get('transpose_for_contour'):
            contour_data = np.asarray(proj2d_display).T
            self._debug_projection(
                f"contour transpose applied for special raw layout: contour_data.shape={getattr(contour_data, 'shape', None)!r} proj2d_display.shape={getattr(proj2d_display, 'shape', None)!r}"
            )

        pos = np.fabs(contour_data * (contour_data > 0.0))
        neg = np.fabs(contour_data * (contour_data < 0.0))
        self._debug_projection(
            f"contour input check: contour_data.shape={getattr(contour_data, 'shape', None)!r} pos.shape={pos.shape!r} neg.shape={neg.shape!r} x_size={indirect_scale.size!r} y_size={direct_scale.size!r} x_first_last={(float(indirect_scale[0]), float(indirect_scale[-1])) if indirect_scale.size else None!r} y_first_last={(float(direct_scale[0]), float(direct_scale[-1])) if direct_scale.size else None!r}"
        )

        if pos.size and np.nanmax(pos) >= levels[0]:
            try:
                ax2d.contour(indirect_scale, direct_scale, pos, levels, colors='r', linewidths=0.5)
            except Exception as exc:
                self._debug_projection(f"ax2d.contour positive failed: {exc!r}")
        if neg.size and np.nanmax(neg) >= levels[0]:
            try:
                ax2d.contour(indirect_scale, direct_scale, neg, levels, colors='b', linewidths=0.5)
            except Exception as exc:
                self._debug_projection(f"ax2d.contour negative failed: {exc!r}")
        self._debug_projection(f"before hilbert: proj1d_raw.size={proj1d_raw.size} minmax={((float(np.nanmin(proj1d_raw)), float(np.nanmax(proj1d_raw))) if proj1d_raw.size else None)!r}")
        if proj1d_proj_raw.size and proj1d_raw.size and proj1d_proj_raw.shape == proj1d_raw.shape:
            try:
                proj_abs_sum = float(np.nansum(np.abs(proj1d_proj_raw)))
                selected_abs_sum = float(np.nansum(np.abs(proj1d_raw)))
                delta_abs_sum = selected_abs_sum - proj_abs_sum
                delta_l2 = float(np.linalg.norm(np.asarray(proj1d_raw, dtype=float) - np.asarray(proj1d_proj_raw, dtype=float)))
                same_shape = proj1d_proj_raw.shape == proj1d_raw.shape
                same_values = bool(np.allclose(np.asarray(proj1d_raw, dtype=float), np.asarray(proj1d_proj_raw, dtype=float), rtol=1e-6, atol=1e-6, equal_nan=True))
                self._debug_projection(
                    f"1D validation fingerprints: proj_abs_sum={proj_abs_sum!r} selected_abs_sum={selected_abs_sum!r} delta_abs_sum={delta_abs_sum!r} delta_l2={delta_l2!r} same_shape={same_shape!r} same_values={same_values!r}"
                )
            except Exception as exc:
                self._debug_projection(f'1D validation fingerprinting failed: {exc!r}')

        proj1d_hilbert = self._hilbert_trace_ng(proj1d_raw)
        phase_state = self._projection_phase_state(indirect_label)
        self._debug_projection(f"phase state for {indirect_label!r}: p0={phase_state['p0']!r} p1={phase_state['p1']!r}")

        phase_dic = None
        trace_scale = indirect_scale
        special_projection_info = None
        if self._is_special_2d_projection_case():
            try:
                phase_dic, trace_scale, special_projection_info = self._special_2d_projection_phase_context(
                    dic,
                    proj1d_raw,
                    direct_label=direct_label,
                    indirect_label=indirect_label,
                    source_path=path,
                )
                self._debug_projection(f"special 2D projection context: {special_projection_info!r}")
            except Exception as exc:
                self._debug_projection(f"special 2D phase context build failed for {path!r}: {exc!r}")
                import traceback
                self._debug_projection(traceback.format_exc())
                phase_dic = None
                trace_scale = indirect_scale

        if phase_dic is None:
            phase_dic = self._projection_1d_phase_dic(dic, proj1d_raw, direct_label=direct_label, indirect_label=indirect_label, source_path=path)

        if trace_scale is None or getattr(trace_scale, 'size', 0) != proj1d_raw.size:
            if self._is_special_2d_projection_case():
                self._debug_projection(
                    f"special 2D trace_scale fallback: trace_scale.size={getattr(trace_scale, 'size', None)!r} proj1d_raw.size={proj1d_raw.size!r} contour_x_scale.size={getattr(contour_x_scale, 'size', None)!r}"
                )
            trace_scale = contour_x_scale

        proj1d_hilbert, proj1d_phased = self._phase_projection_trace(
            proj1d_hilbert, indirect_label, phase_state['p0'], phase_state['p1'], dic=phase_dic
        )
        self._debug_projection(f"phase result: hilbert.shape={getattr(proj1d_hilbert, 'shape', None)!r} phased.shape={getattr(proj1d_phased, 'shape', None)!r} phased.dtype={getattr(proj1d_phased, 'dtype', None)!r}")
        proj1d_display = np.real(proj1d_phased)

        line1d, = ax1d.plot(trace_scale, proj1d_display, color='k', linewidth=1.0)

        cache = {
            'ax2d': ax2d,
            'ax1d': ax1d,
            'line1d': line1d,
            'seed_overlay': None,
            'direct_label': direct_label,
            'indirect_label': indirect_label,
            'direct_scale': direct_scale,
            'indirect_scale': indirect_scale,
            'contour_x_scale': contour_x_scale,
            'contour_y_scale': contour_y_scale,
            'proj2d': proj2d,
            'proj1d_raw': proj1d_raw,
            'proj1d_seed_complex': proj1d_seed_complex,
            'proj1d_hilbert': proj1d_hilbert,
            'proj1d_phased': proj1d_phased,
            'proj1d_mask_stats': mask_stats,
            'proj1d_mask_selection': mask_selection,
            'projection_mode': view_mode_display,
            'phase_p0': phase_state['p0'],
            'phase_p1': phase_state['p1'],
            'source': path,
            'dic': dic,
            'phase_dic': phase_dic,
        }

        ax2d.set_title(f'{direct_label}/{indirect_label}', fontsize=8)
        ax2d.set_xlabel('')
        ax2d.set_ylabel(f'{direct_label} (ppm)', fontsize=8)
        ax2d.tick_params(labelsize=7)
        ax2d.tick_params(labelbottom=False)
        if view_mode == '1D':
            self._debug_projection('1D view mode: 2D panel is plotting the inverse-FFT display only; seed overlay intentionally disabled')
        elif view_mode == 'mask':
            self._debug_projection('matched view mode: 2D panel remains the raw 2D projection while the 1D trace uses accepted detrended slices only')
            self._debug_projection(f'matched view mode selection_strength={mask_selection!r}')

        ax1d.set_xlabel(f'{indirect_label} (ppm)', fontsize=8)
        ax1d.set_ylabel({'1D': 'Seed', 'mask': 'Matched'}.get(view_mode, 'Sum'), fontsize=8)
        ax1d.tick_params(labelsize=7)

        if indirect_scale.size > 1:
            ax2d.set_xlim(float(indirect_scale[0]), float(indirect_scale[-1]))
            ax1d.set_xlim(float(indirect_scale[0]), float(indirect_scale[-1]))
        if direct_scale.size > 1:
            ax2d.set_ylim(float(direct_scale[0]), float(direct_scale[-1]))

        if proj1d_display.size > 0:
            pmin = float(np.nanmin(proj1d_display))
            pmax = float(np.nanmax(proj1d_display))
            if pmin == pmax:
                pad = abs(pmax) * 0.05 if pmax != 0 else 1.0
            else:
                pad = max(abs(pmin), abs(pmax)) * 0.05
            ax1d.set_ylim(pmin - pad, pmax + pad)

        self._save_projection_cache(cache)
        return cache

    def _plot_projection_file(self, ax2d, ax1d, path: str, spec) -> dict:
        self._debug_projection(f'_plot_projection_file loading path={path!r}')
        dic, data = ng.pipe.read(path)
        return self._plot_projection_data(ax2d, ax1d, path, spec, dic, data)

    def _apply_projection_figure_margins(self, ncols: int) -> None:
        """Pack projection axes tightly while preserving title/label visibility."""
        if ncols <= 0:
            ncols = 1
        try:
            width, height = self.fig.get_size_inches()
        except Exception:
            width, height = 14.0, 9.0

        # Reserve enough space for axis labels/titles, but keep the plot area
        # as large as possible and adapt to the number of columns.
        left = min(0.08, 0.05 + 0.003 * max(0, ncols - 1))
        right = 0.992
        top = 0.96 if height >= 7.0 else 0.95
        bottom = 0.085 if height >= 7.0 else 0.10
        wspace = 0.18 if ncols <= 2 else 0.14
        hspace = 0.16

        try:
            self.fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom, wspace=wspace, hspace=hspace)
        except Exception:
            pass

    def redraw_view(self) -> None:
        self._full_redraw()

    def _full_redraw(self) -> None:
        self._debug_projection('_full_redraw start')
        self.fig.clf()
        self._clear_projection_cache()

        special_case = self._is_special_2d_projection_case()
        special_bundle = None
        cached_entries = None

        # A two-spectral-dimension projection must be self-loading.  On a cold
        # project open the main NMR tab has not necessarily populated
        # state.spectra['raw'] yet.  The old first pass found the projection
        # filename, but then plotted it through the generic path using
        # incomplete in-memory dimensional metadata.  Visiting the main tab
        # happened to initialise that metadata, so reopening this window then
        # worked.  Load the selected NMRPipe file/header here on every special
        # 2D pass instead; the window no longer depends on main-tab visit order.
        if special_case:
            special_bundle = self._special_2d_projection_bundle()
            if special_bundle is not None:
                entries = [(special_bundle['path'], special_bundle['title'])]
            else:
                entries = self._expected_projection_entries()
        else:
            entries = self._expected_projection_entries()
        self._debug_projection(f'_full_redraw special_case={special_case!r} expected entries count={len(entries)} entries={entries!r}')
        if not special_case and not entries and self.state is not None:
            cached = self._cached_projection_entries()
            if cached:
                entries = [(path, title) for path, title, _payload in cached]
                cached_entries = cached
                self._debug_projection('Using cached in-memory projections because files were not found')
                for _path, _title, payload in cached:
                    self._projection_cache.append({
                        'source': _path,
                        'title': _title,
                        'cached_projection': True,
                        'payload_keys': sorted(payload.keys()),
                    })
        if entries and not self._contour_baseline_ready:
            self.maxval = self._projection_maxval(entries, special_bundle=special_bundle, cached_entries=cached_entries)
            self._sync_contour_controls_from_baseline()
            self._contour_baseline_ready = True
        self._rebuild_phase_controls(entries)

        if not entries:
            debug_text = self._debug_summary()
            if special_case:
                self._debug_projection('2D special-case raw projection could not be resolved')
            else:
                self._debug_projection('No 2D projections found')
            self._debug_projection(debug_text)
            ax = self.fig.add_subplot(111)
            ax.axis('off')
            ax.text(0.5, 0.98, 'No 2D projections found', ha='center', va='top', transform=ax.transAxes, fontsize=12)
            ax.text(0.02, 0.92, debug_text, ha='left', va='top', transform=ax.transAxes, family='monospace', fontsize=8)
            self._apply_projection_figure_margins(1)
            self.canvas.draw()
            wx.CallAfter(self._debug_widget_state, 'post-redraw no-entries')
            return

        from matplotlib.gridspec import GridSpec
        ncols = len(entries)
        gs = GridSpec(2, ncols, figure=self.fig, height_ratios=[2.0, 1.0], hspace=0.10, wspace=0.20)

        for idx, (path, title) in enumerate(entries):
            ax2d = self.fig.add_subplot(gs[0, idx])
            ax1d = self.fig.add_subplot(gs[1, idx], sharex=ax2d)
            if special_case and special_bundle is not None and path == special_bundle['path']:
                spec = {
                    'direct': special_bundle['direct'],
                    'indirect': special_bundle['indirect'],
                    'title': special_bundle['title'],
                }
                try:
                    cache = self._plot_projection_data(ax2d, ax1d, path, spec, special_bundle['dic'], special_bundle['data'])
                except Exception as exc:
                    self._debug_projection(f'Special-case 2D plot error for {path!r}: {exc!r}')
                    import traceback
                    self._debug_projection(traceback.format_exc())
                    ax2d.axis('off')
                    ax1d.axis('off')
                    ax2d.text(0.5, 0.5, f"{title}\n{exc.__class__.__name__}", ha='center', va='center', transform=ax2d.transAxes)
                    cache = {'ax2d': ax2d, 'ax1d': ax1d, 'title': title, 'source': path}
            else:
                spec = {'direct': self._direct_label(), 'indirect': title.split('.', 1)[-1], 'title': title}
                try:
                    cache = self._plot_projection_file(ax2d, ax1d, path, spec)
                except Exception as exc:
                    self._debug_projection(f'Plot error for {path!r}: {exc!r}')
                    import traceback
                    self._debug_projection(traceback.format_exc())
                    ax2d.axis('off')
                    ax1d.axis('off')
                    ax2d.text(0.5, 0.5, f"{title}\n{exc.__class__.__name__}", ha='center', va='center', transform=ax2d.transAxes)
                    cache = {'ax2d': ax2d, 'ax1d': ax1d, 'title': title, 'source': path}
            self._projection_cache.append(cache)

        self._apply_projection_figure_margins(ncols)
        self.canvas.draw()
        wx.CallAfter(self._debug_widget_state, f'post-redraw ncols={ncols}')

    def load_projections(self) -> None:
        self._full_redraw()

    def refresh_contours(self, *, force_full: bool = False) -> None:
        self._full_redraw()

    def OnMaskStrictnessChanged(self, event):
        if self._initializing:
            event.Skip()
            return
        try:
            value = self._mask_selection_value()
            if self.state is not None:
                try:
                    self.state.metadata['projection_mask_selection'] = value
                except Exception:
                    pass
            self._update_mask_slider_widgets()
            mode = self._projection_view_mode()
            self._debug_projection(f'Mask strictness changed -> {value!r} (view_mode={mode!r})')
            self._full_redraw()
        except Exception as exc:
            self._debug_projection(f'Mask strictness refresh failed: {exc!r}')
        event.Skip()

    def set_projection_view_mode(self, mode: str) -> None:
        """Internal hook for future programmatic mode changes."""
        self._set_projection_view_mode(mode)

    def OnViewModeChanged(self, event):
        if not self._initializing:
            try:
                mode = self._projection_view_mode()
                if self.state is not None:
                    try:
                        self.state.metadata['projection_view_mode'] = mode
                    except Exception:
                        pass
                self._full_redraw()
            except Exception as exc:
                self._debug_projection(f'View mode refresh failed: {exc!r}')
        event.Skip()

    # ------------------------------------------------------------------
    def OnContourTextChanged(self, event):
        if self._initializing or self._syncing_contour_controls:
            event.Skip()
            return
        source = None
        try:
            source = event.GetEventObject()
        except Exception:
            source = None
        try:
            if source is self.cmin:
                self._sync_contour_thresh_from_cmin()
            elif source is self.cthresh:
                self._sync_contour_cmin_from_thresh()
        except Exception:
            pass
        if self._redraw_timer is not None:
            try:
                self._redraw_timer.Stop()
            except Exception:
                pass
        try:
            self._redraw_timer = wx.CallLater(200, self.refresh_contours)
        except Exception:
            self.refresh_contours()
        event.Skip()

    def OnMove(self, event):
        if not self._initializing and hasattr(self, 'phaseFrame') and self.phaseFrame.IsShown():
            wx.CallAfter(self._position_phase_frame)
        event.Skip()

    def OnSize(self, event):
        if not self._initializing:
            if hasattr(self, 'phaseFrame') and self.phaseFrame.IsShown():
                wx.CallAfter(self._fit_phase_frame)
            self._background_ready = False
            try:
                wx.CallAfter(self._update_phase_control_layout)
            except Exception:
                pass
            try:
                if self._projection_cache:
                    self._apply_projection_figure_margins(max(1, len(self._projection_cache)))
                    self.canvas.draw_idle()
            except Exception:
                pass
        event.Skip()

    # ------------------------------------------------------------------
    def OnReprocess(self, event):
        try:
            self._run_silent_reprocess()
        except Exception as exc:
            self._debug_projection(f'Re-process failed: {exc!r}')
            from spinDecon.gui.dialogs.errors import errorMessage
            errorMessage('Could not re-process projections.\n\n' + str(exc))
        if event is not None:
            event.Skip()

    def collect_updates(self, update_state=True):
        """Return persistent projection settings (temporary phase deltas excluded)."""
        updates = {
            'cmin': self.cmin.GetValue().strip(),
            'cthresh': self.cthresh.GetValue().strip(),
            'cfac': self.cfac.GetValue().strip(),
            'cnum': self.cnum.GetValue().strip(),
        }
        if self.state is not None and update_state:
            self.state.metadata.setdefault('projection_contours', {}).update(updates)
            self.state.update_gui_settings(updates)
        return updates

    def OnClose(self, event):
        try:
            self.collect_updates()
        except Exception:
            pass
        try:
            if hasattr(self, 'phaseFrame'):
                self.phaseFrame.Destroy()
        except Exception:
            pass
        try:
            if getattr(self.process_parent, 'projections_frame', None) is self:
                self.process_parent.projections_frame = None
        except Exception:
            pass
        self.Destroy()


__all__ = ['ProjectionsFrame']
