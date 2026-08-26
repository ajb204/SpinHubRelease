"""Single-project PDF summary generation for the decon GUI.

The GUI is the source of truth: figures are exported from the live matplotlib
views and settings are read from the current controls/shared ProjectState.
"""
from __future__ import annotations

import math
import os
import re
import csv
import json
import shutil
import subprocess
from pathlib import Path


def _latex_escape(value):
    text = str(value if value is not None else '')
    repl = {
        '\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
        '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
        '^': r'\textasciicircum{}',
    }
    return ''.join(repl.get(ch, ch) for ch in text)


def _ascii_detail(text):
    return (str(text).replace('μ', 'mu').replace('σ', 'sigma')
            .replace('×', 'x').replace('–', '-').replace('—', '-'))


def _ctrl_value(frame, name, default=''):
    ctrl = getattr(frame, name, None)
    if ctrl is None:
        return default
    try:
        if hasattr(ctrl, 'GetStringSelection'):
            value = ctrl.GetStringSelection()
            if value:
                return value
        return ctrl.GetValue()
    except Exception:
        try:
            return ctrl.GetSelection()
        except Exception:
            return default


def _count_rows(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as handle:
            return sum(1 for line in handle if line.strip())
    except Exception:
        return ''


def _write_kv_table(out, title, rows):
    rows = [(k, v) for k, v in rows if v not in ('', None)]
    if not rows:
        return
    # An empty title is intentional when the table sits directly beneath an
    # existing section heading (for example Noise).  Do not emit an empty
    # subsection: LaTeX still reserves vertical heading space for it.
    if title:
        out.write('\\subsection*{' + _latex_escape(title) + '}\n')
    # Data entries are deliberately smaller than report headings.  Fixed-width
    # paragraph columns also allow long project filenames to wrap inside their
    # minipage instead of running underneath the neighbouring figure.
    out.write('\\footnotesize\n')
    out.write('\\begin{tabular}{@{}>{\\bfseries}p{0.42\\linewidth}p{0.54\\linewidth}@{}}\n')
    for key, value in rows:
        out.write(_latex_escape(key) + ' & ' + _latex_escape(value) + '\\\\' + '\n')
    out.write('\\end{tabular}\n\\normalsize\n\n')


def _save_peak_views(frame, report_dir, warnings):
    """Render raw and reference-peak-overlay bore views without showing a window."""
    try:
        # PeakFrame copies the reference list during construction.  On a cold
        # project load the filename can be configured while the list itself is
        # not yet materialised, which produced an apparently valid but empty
        # ornamented figure.  Synchronise spectrum -> reference first.
        ensure_ref = getattr(frame, 'ensure_workflow_reference_stage_loaded', None)
        if callable(ensure_ref) and not ensure_ref():
            raise RuntimeError('the spectrum/reference peak list could not be loaded')
        from spinDecon.gui.workspaces.peaks import peakFrame
        pf = peakFrame(frame, showFlg=False)
        try:
            pf.canvas.print_figure(str(report_dir / 'peak.pdf'))
            # Preserve the historical ornamented peak view using the current
            # PeakFrame's own 'Peaks' overlay control and redraw logic.
            if hasattr(pf, 'cb_grid'):
                pf.cb_grid.SetValue(True)
                pf.draw_figure()
                pf.canvas.print_figure(str(report_dir / 'peak_ornament.pdf'))
        finally:
            pf.Destroy()
    except Exception as exc:
        warnings.append('Peak/bore figure unavailable: %s' % exc)


def _save_peak_shape(frame, report_dir, warnings):
    """Export the live Fit Peaks views and return reportable fit diagnostics."""
    info = []
    try:
        if hasattr(frame, 'OnButtonPeakFit'):
            pf = frame.OnButtonPeakFit(None, showFlg=False)
        else:
            from spinDecon.gui.workspaces.peaks import peakFitFrame
            pf = peakFitFrame(frame, showFlg=False)
        if pf is None:
            raise RuntimeError('Fit Peaks did not create a peak-fit window')
        try:
            # Automatic representative-peak detection runs on a worker thread.
            # Give it a bounded opportunity to finish so the report captures the
            # same overlay/histograms/statistics that the user sees in Fit Peaks.
            thread = getattr(pf, '_picker_thread', None)
            if thread is not None and thread.is_alive():
                import time
                try:
                    import wx
                except Exception:
                    wx = None
                deadline = time.time() + 15.0
                while thread.is_alive() and time.time() < deadline:
                    if wx is not None:
                        try: wx.YieldIfNeeded()
                        except Exception: pass
                    time.sleep(0.03)
                if wx is not None:
                    try: wx.YieldIfNeeded()
                    except Exception: pass

            result = getattr(pf, '_peak_search_result', None)
            representative = 0
            if result is not None and getattr(result, 'representative_indices', None) is not None:
                representative = len(result.representative_indices)
            if result is not None:
                info.append(('Candidate points', getattr(result, 'candidate_count', '')))
            maxima = getattr(pf, 'maxima', None)
            used = len(maxima) if maxima is not None else 0
            info.extend([
                ('Isolated peaks', representative),
                ('Peaks fit', used),
                ('Widths linked', 'yes' if bool(pf.linkWidths.GetValue()) else 'no'),
            ])
            pf.canvas.draw()
            pf.canvas.print_figure(str(report_dir / 'shape.pdf'))
            # The lower Fit Peaks canvas contains one calibrated width
            # histogram per spectral dimension.
            if hasattr(pf, 'canvas_widths') and hasattr(pf, 'fig_widths'):
                pf.canvas_widths.draw()
                pf.fig_widths.savefig(str(report_dir / 'shape_widths.pdf'))
        finally:
            pf.Destroy()
    except Exception as exc:
        warnings.append('Peak-shape figure unavailable: %s' % exc)
    return info

def _save_projection(frame, report_dir, warnings):
    """Save the Projection window twice: raw, then Peaks + ShowCalc.

    ``projection.pdf`` is drawn with both overlays disabled.  When a valid
    deconvolved spectrum is available, ``projection_decon.pdf`` is drawn from
    the *same Projection window* after enabling its Peaks and ShowCalc
    checkboxes, exactly as an interactive user would do before pressing Draw.
    """
    notebook = frame.parent
    try:
        # A project can have a valid .decon file on disk even if this particular
        # GUI session has not yet run Analyse.  Load it before constructing the
        # lazy Projection page so that the page sees the calculated data.
        spectrum = ''
        try:
            spectrum = frame._resolve_input_path(frame.infileBox.GetValue())
        except Exception:
            spectrum = _ctrl_value(frame, 'infileBox')
        if spectrum:
            spectrum = os.path.abspath(os.path.expanduser(str(spectrum)))
        try:
            decon_file = frame._active_deconvolution_path(spectrum)
        except Exception:
            decon_file = spectrum + '.decon' if spectrum else ''
        if decon_file and os.path.isfile(decon_file) and not bool(getattr(frame, 'DECON', 0)):
            loader = getattr(frame, '_load_decon_outputs', None)
            if loader is not None:
                # _load_decon_outputs accepts the product's unsuffixed source.
                loader(decon_file[:-6] if decon_file.endswith('.decon') else spectrum)

        notebook.open_workflow('inspect')
        proj = notebook.tabTwo
        old_calc = bool(proj.cb_calc.GetValue()) if hasattr(proj, 'cb_calc') else False
        old_peaks = bool(proj.cb_grid.GetValue()) if hasattr(proj, 'cb_grid') else False

        # Figure 1: raw Projection-window view.
        if hasattr(proj, 'cb_calc'):
            proj.cb_calc.SetValue(False)
        if hasattr(proj, 'cb_grid'):
            proj.cb_grid.SetValue(False)
        proj.draw_figure()
        proj.canvas.print_figure(str(report_dir / 'projection.pdf'))

        # Figure 2: reproduce the successful interactive operation literally:
        # select Peaks + ShowCalc first, then press Draw.  In particular, do
        # NOT inspect proj.no_decon before this draw.  Drawing the raw view
        # above intentionally has ShowCalc disabled and therefore sets
        # no_decon=True for 3D data; using that value as a precondition would
        # prevent the calculated view from ever being attempted.
        if hasattr(proj, 'cb_calc') and hasattr(proj, 'cb_grid'):
            proj.cb_grid.SetValue(True)   # Peaks
            proj.cb_calc.SetValue(True)   # ShowCalc
            proj.draw_figure()            # same action as the Draw button

            # draw_figure() has now asked the Projection window/DataStore for
            # the calculated projection views, so no_decon is meaningful only
            # at this point.  Save exactly what the Projection window drew.
            if not bool(getattr(proj, 'no_decon', False)):
                proj.canvas.print_figure(str(report_dir / 'projection_decon.pdf'))
            elif decon_file and os.path.isfile(decon_file):
                warnings.append('Calculated Projection view unavailable after Peaks + ShowCalc redraw: %s' % decon_file)

        # Restore the interactive Projection window exactly as we found it.
        if hasattr(proj, 'cb_calc'):
            proj.cb_calc.SetValue(old_calc)
        if hasattr(proj, 'cb_grid'):
            proj.cb_grid.SetValue(old_peaks)
        proj.draw_figure()
    except Exception as exc:
        warnings.append('Projection figure unavailable: %s' % exc)


def _parameter_file(frame):
    """Return the concrete project/system parameter file, if available."""
    state = getattr(frame, 'state', None)
    candidates = []
    if state is not None:
        candidates.append(getattr(state, 'parameter_file', ''))
    candidates.append(getattr(frame, 'deconParFile', ''))
    for raw in candidates:
        raw = str(raw or '').strip()
        if not raw:
            continue
        # Project/system parameter files belong to WorkingDir, not SpecPath.
        path = Path(raw).expanduser()
        if not path.is_absolute():
            working = getattr(state, 'working_dir', '') if state is not None else ''
            working = working or _ctrl_value(frame, 'dirBox') or os.getcwd()
            path = Path(working).expanduser() / path
        if path.is_file():
            return path
    return None


def _read_parameter_map(frame):
    """Read the saved system file without importing the legacy SettingsUnidec module."""
    path = _parameter_file(frame)
    if path is None:
        return {}, None
    values = {}
    try:
        for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
            fields = line.split()
            if not fields:
                continue
            key = fields[0]
            if len(fields) >= 3 and fields[1] == '=':
                values[key] = ' '.join(fields[2:])
            elif len(fields) >= 2:
                values[key] = ' '.join(fields[1:])
    except Exception:
        return {}, path
    return values, path


def _dimension_labels(frame, params=None, ndim=None):
    params = params or {}
    if ndim is None:
        try:
            ndim = int(getattr(frame, 'dim', 0) or getattr(getattr(frame, 'state', None), 'dimension', 0) or 1)
        except Exception:
            ndim = 1
    live = list(getattr(frame, 'labb', []) or [])
    labels = []
    for i in range(max(int(ndim), 1)):
        label = str(live[i] or '').strip() if i < len(live) else ''
        if not label:
            label = str(params.get('label%d' % i, '') or '').strip()
        labels.append(label or ('Dimension %d' % (i + 1)))
    return labels


def _processing_table(frame):
    params, _ = _read_parameter_map(frame)
    if not params:
        return [], []
    markers = ('window0', 'win2Val0', 'win3Val0', 'firstPoint0', 'flip0', 'p0_1', 'window1', 'firstPoint1')
    if not any(k in params for k in markers):
        return [], []
    try:
        ndim = int(str(params.get('dim', getattr(frame, 'dim', 1)))[0])
    except Exception:
        ndim = 1
    labels = list(reversed(_dimension_labels(frame, params, ndim)))
    window_names = {'0': 'GM', '1': 'SP', '2': 'EM'}
    rows = []
    specs = [('Phase p0','p0_','p0'), ('Phase p1','p1_','p1'), ('FT mode','flip',None),
             ('Linear prediction','lp',None), ('F1180','f180',None),
             ('Baseline correction','bl',None), ('Window','window',None), ('Window parameter 1','win2Val',None),
             ('Window parameter 2','win3Val',None), ('First point','firstPoint',None)]
    for title, key, direct in specs:
        vals=[]
        for i in range(ndim):
            if key == 'f180':
                v = params.get('f%d180' % i, '')
            elif i == 0 and direct and direct in params:
                v = params.get(direct, '')
            else:
                v = params.get(key + str(i), '')
            v = str(v or '').strip()
            if key == 'window' and v:
                v = window_names.get(v, v)
            vals.append(v)
        if any(vals):
            rows.append((title, vals))
    pol = str(params.get('pol', '') or '').strip()
    if pol:
        rows.append(('Polynomial order', [pol] + [''] * (ndim - 1)))

    # The Process window's Digital solvent suppress checkbox is a direct-
    # dimension-only setting.  processingFrame loads it from the saved `sol`
    # parameter in deconParFile, so report the same persisted value here.
    sol = str(params.get('sol', 'n') or 'n').strip().lower()
    sol = 'y' if sol == 'y' else 'n'
    rows.append(('Digital solvent removal', [sol] + [''] * (ndim - 1)))
    return labels, rows


def _write_matrix_table(out, title, columns, rows):
    """Write either a parameter matrix or an ordinary row table.

    Historical callers pass rows as ``(parameter, [dimension values...])`` and
    ``columns`` contains only the dimension headings.  Report adapters such as
    Pseudo3D and Decay instead return a conventional table where each row is a
    flat sequence matching ``columns``.  Supporting both shapes here keeps the
    existing acquisition/processing tables unchanged while allowing report
    data providers to expose the same rows used by their GUI controls.
    """
    if not rows or not columns:
        return
    if title:
        out.write('\\subsection*{' + _latex_escape(title) + '}\n')

    first = rows[0]
    legacy_matrix = (isinstance(first, (list, tuple)) and len(first) == 2 and
                     isinstance(first[1], (list, tuple)))
    if legacy_matrix:
        headings = ['Parameter'] + list(columns)
        table_rows = [[name] + list(values) for name, values in rows]
    else:
        headings = list(columns)
        table_rows = [list(row) for row in rows]

    out.write('\\footnotesize\n')
    out.write('\\begin{tabular}{@{}' + ('l' * len(headings)) + '@{}}\n\\toprule\n')
    out.write(' & '.join('\\textbf{' + _latex_escape(col) + '}' for col in headings))
    out.write('\\\\\n\\midrule\n')
    for row in table_rows:
        # Be defensive about partially populated adapters: pad short rows and
        # retain extra values rather than failing report generation.
        values = row + [''] * max(0, len(headings) - len(row))
        values = values[:len(headings)]
        out.write(' & '.join(_latex_escape(value) for value in values) + '\\\\\n')
    out.write('\\bottomrule\n\\end{tabular}\n\\normalsize\n\n')

def _write_long_table(out, columns, rows, font='\\scriptsize', compact=False):
    """Multipage table with repeated headers."""
    if not rows or not columns:
        return
    n = len(columns)
    if compact:
        out.write('\\begingroup' + font + '\\setlength{\\tabcolsep}{2.2pt}\\renewcommand{\\arraystretch}{0.92}\n')
    else:
        out.write(font + '\n')
    out.write('\\begin{longtable}{@{}' + ('l' * n) + '@{}}\n\\toprule\n')
    header = ' & '.join('\\textbf{' + _latex_escape(c) + '}' for c in columns) + '\\\\\n'
    out.write(header + '\\midrule\n\\endfirsthead\n\\toprule\n' + header + '\\midrule\n\\endhead\n')
    out.write('\\midrule\\multicolumn{' + str(n) + '}{r}{\\scriptsize Continued on next page}\\\\\n\\endfoot\n')
    out.write('\\bottomrule\n\\endlastfoot\n')
    for row in rows:
        vals = list(row) + [''] * max(0, n-len(row))
        out.write(' & '.join(_latex_escape(v) for v in vals[:n]) + '\\\\\n')
    out.write('\\end{longtable}\n')
    out.write('\\endgroup\n\n' if compact else '\\normalsize\n\n')


def _write_long_table_with_split(out, columns, rows, split_at, font='\\scriptsize', compact=False):
    """Multipage table with a vertical rule before ``split_at``."""
    if not rows or not columns:
        return
    n = len(columns)
    spec = []
    for i in range(n):
        if i == split_at:
            spec.append('|')
        spec.append('l')
    if compact:
        out.write('\\begingroup' + font + '\\setlength{\\tabcolsep}{2.2pt}\\renewcommand{\\arraystretch}{0.92}\n')
    else:
        out.write(font + '\n')
    out.write('\\begin{longtable}{@{}' + ''.join(spec) + '@{}}\n\\toprule\n')
    header = ' & '.join('\\textbf{' + _latex_escape(c) + '}' for c in columns) + '\\\\\n'
    out.write(header + '\\midrule\n\\endfirsthead\n\\toprule\n' + header + '\\midrule\n\\endhead\n')
    out.write('\\midrule\\multicolumn{' + str(n) + '}{r}{\\scriptsize Continued on next page}\\\\\n\\endfoot\n')
    out.write('\\bottomrule\n\\endlastfoot\n')
    for row in rows:
        vals = list(row) + [''] * max(0, n-len(row))
        out.write(' & '.join(_latex_escape(v) for v in vals[:n]) + '\\\\\n')
    out.write('\\end{longtable}\n')
    out.write('\\endgroup\n\n' if compact else '\\normalsize\n\n')


def _append_fitting_to_full_peaks(full_headers, full_rows, pseudo):
    """Append displayed fitting-window values to matching full-list peaks."""
    fit_cols = list(pseudo.get('columns') or [])
    fit_rows = [list(row) for row in (pseudo.get('rows') or [])]
    if pseudo.get('kind') != '2d' or not fit_cols or not fit_rows:
        return full_headers, full_rows, None
    fit_by_peak = {str(row[0]): row for row in fit_rows if row}
    # Peak is already present on the left side; Group and all fitted values are
    # useful review metadata, so append every fitting column except Peak.
    appended_cols = fit_cols[1:]
    joined = []
    for row in full_rows:
        name = str(row[0]) if row else ''
        fit = fit_by_peak.get(name)
        # Header and data are driven by the same fitting-column contract.
        # Pad/truncate legacy fitting rows so values can never drift into a
        # neighbouring labelled column.
        fit_values = list(fit[1:]) if fit else []
        fit_values = (fit_values + [''] * len(appended_cols))[:len(appended_cols)]
        joined.append(list(row) + fit_values)
    return list(full_headers) + appended_cols, joined, len(full_headers)


def _write_two_column_long_table(out, columns, rows):
    """Reference peaks in two blocks with a visible gutter and repeated headers."""
    if not rows or not columns:
        return
    half = (len(rows) + 1) // 2
    left, right = rows[:half], rows[half:]
    n = len(columns)
    combined = []
    blank = [''] * n
    for i in range(half):
        combined.append(list(left[i]) + (list(right[i]) if i < len(right) else blank))
    spec = '@{}' + ('l' * n) + '@{\\hspace{5mm}}' + ('l' * n) + '@{}'
    out.write('\\scriptsize\n\\begin{longtable}{' + spec + '}\n\\toprule\n')
    headers = list(columns) + list(columns)
    row_end = '\\\\' + '\n'
    header = ' & '.join('\\textbf{' + _latex_escape(c) + '}' for c in headers) + row_end
    out.write(header + '\\midrule\n\\endfirsthead\n\\toprule\n' + header + '\\midrule\n\\endhead\n')
    out.write('\\midrule\\multicolumn{' + str(2*n) + '}{r}{\\scriptsize Continued on next page}' + row_end + '\\endfoot\n')
    out.write('\\bottomrule\n\\endlastfoot\n')
    for row in combined:
        vals = list(row) + [''] * max(0, 2*n-len(row))
        out.write(' & '.join(_latex_escape(v) for v in vals[:2*n]) + row_end)
    out.write('\\end{longtable}\n\\normalsize\n\n')


def _noise_summary_rows(frame):
    """Return the concise statistical subset of the Noise Detail window."""
    stats = getattr(frame, 'noiseStats', None) or {}
    if not stats:
        return [('Status', 'No noise statistics available')]

    def number(key, fmt='.6g'):
        try:
            return format(float(stats.get(key)), fmt)
        except Exception:
            return ''

    rows = [
        ('Spectral points', '{:,}'.format(int(stats.get('points', 0)))),
        ('Points analysed', '{:,}'.format(int(stats.get('sampled_points', 0)))),
        ('Gaussian centre', number('centre')),
        ('Fitted noise sigma', number('noise_sigma')),
        ('MAD noise sigma', number('noise_mad_sigma')),
        ('Central population', number('core_fraction', '.2%')),
        ('Maximum intensity', number('max_intensity')),
        ('Maximum +S/N', number('max_snr', '.2f')),
    ]
    tails = stats.get('tail_counts', {}) or {}
    for level in (3, 4, 5):
        tail = tails.get(level, {}) or tails.get(str(level), {}) or {}
        try:
            pos = int(tail.get('positive', 0)); neg = int(tail.get('negative', 0))
            rows.append(('%s sigma (+ / -)' % level, '%s / %s' % (format(pos, ','), format(neg, ','))))
        except Exception:
            pass
    return rows


def _display_filename(value):
    """Return only the filename for a project file value."""
    text = str(value or '').strip()
    return Path(text).name if text else ''


def _relative_project_path(frame, value):
    """Return a project path relative to the current working directory."""
    text = str(value or '').strip()
    if not text:
        return ''
    state = getattr(frame, 'state', None)
    base = Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()
    try:
        path = Path(frame._project_path(text)).expanduser()
    except Exception:
        path = Path(text).expanduser()
    try:
        return os.path.relpath(str(path.resolve()), str(base.resolve()))
    except Exception:
        return text



def _resolve_project_file(frame, value):
    text = str(value or '').strip()
    if not text:
        return None
    try:
        path = Path(frame._project_path(text)).expanduser()
    except Exception:
        path = Path(text).expanduser()
    if not path.is_absolute():
        state = getattr(frame, 'state', None)
        base = Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()
        path = base / path
    return path


def _param_ci(params, *names, default=''):
    """Case-insensitive lookup for historical system-file parameter names."""
    wanted = {str(name).lower() for name in names}
    for key, value in (params or {}).items():
        if str(key).lower() in wanted:
            return value
    return default


def _resolve_nus_file(frame, params):
    """Resolve nusfil relative to fiddir, matching the conversion GUI semantics."""
    nus_value = str(_param_ci(params, 'nusfil', 'nusFil', default='') or '').strip()
    fid_value = str(_param_ci(params, 'fiddir', default='') or '').strip()
    state = getattr(frame, 'state', None)
    working = Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()

    if not nus_value:
        return None, nus_value, fid_value

    nus = Path(nus_value).expanduser()
    if nus.is_absolute():
        return nus, nus_value, fid_value

    # The saved NUS schedule is relative to the raw-FID directory, not to the
    # project directory.  This is also how ConversionFrame normalises it.
    fid = Path(fid_value).expanduser() if fid_value else Path('.')
    if not fid.is_absolute():
        fid = working / fid
    return fid / nus, nus_value, fid_value


def _nus_summary(frame, params, inst=None):
    """Summarise NUS sampling and return (report rows, diagnostic rows, maxima)."""
    rows = []
    debug = []
    path, nus_value, fid_value = _resolve_nus_file(frame, params)
    debug.append(('nusfil (system)', nus_value or '<not set>'))
    debug.append(('fiddir (system)', fid_value or '<not set>'))
    if path is None:
        debug.append(('NUS status', 'No nusfil value found'))
        return rows, debug, []

    debug.append(('Resolved NUS file', str(path)))
    debug.append(('NUS file exists', 'yes' if path.is_file() else 'no'))
    if not path.is_file():
        debug.append(('NUS status', 'Schedule could not be opened'))
        return rows, debug, []

    samples = []
    try:
        for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = []
            for field in line.split():
                try:
                    vals.append(int(float(field)))
                except Exception:
                    pass
            if vals:
                samples.append(vals)
    except Exception as exc:
        debug.append(('NUS status', 'Read failed: %s' % exc))
        return rows, debug, []
    if not samples:
        debug.append(('NUS status', 'Schedule contains no numeric samples'))
        return rows, debug, []

    width = max(len(row) for row in samples)
    maxima = []
    for col in range(width):
        vals = [row[col] for row in samples if col < len(row)]
        if vals:
            maxima.append(max(vals))

    # Schedule coordinates are zero based, hence max + 1 points per axis.
    sizes = [maximum + 1 for maximum in maxima]
    reconstructed = 1
    for size in sizes:
        reconstructed *= size

    rows = [('NUS schedule', _relative_project_path(frame, str(path))),
            ('Samples', str(len(samples)))]
    if sizes:
        rows.append(('Reconstructed', ' x '.join(str(x) for x in sizes)))
        rows.append(('Reconstructed size', str(reconstructed)))
    if reconstructed > 0:
        rows.append(('Sparsity', '{:.2f}%'.format(100.0 * len(samples) / float(reconstructed))))

    debug.append(('NUS numeric columns', str(len(maxima))))
    debug.append(('NUS maxima', ', '.join(str(x) for x in maxima)))

    # vpar is the authoritative source for acquisition sweep widths after the
    # Conversion window has configured it.  Include the direct dimension first
    # (sw and xN), followed by the NUS indirect dimensions (sw1/sw2/sw3).
    if inst is None:
        debug.append(('Sweep-width status', 'vpar unavailable'))
    else:
        debug.append(('vpar initialized', 'yes' if getattr(inst, 'initialized', False) else 'no'))
        labels = list(getattr(inst, 'labb', []) or [])
        raw_sw = getattr(inst, 'sw', None)
        raw_xn = getattr(inst, 'xN', getattr(inst, 'np2', None))
        debug.append(('sw', str(raw_sw) if raw_sw is not None else '<missing>'))
        debug.append(('xN', str(raw_xn) if raw_xn is not None else '<missing>'))
        try:
            direct_sw = float(raw_sw)
            direct_n = float(raw_xn)
            # Direct-dimension values are presented in the acquisition matrix
            # rather than repeated as key/value rows.
        except Exception:
            pass
        for col, maximum in enumerate(maxima):
            attr = 'sw%d' % (col + 1)
            raw_sw = getattr(inst, attr, None)
            debug.append((attr, str(raw_sw) if raw_sw is not None else '<missing>'))
            try:
                sw = float(raw_sw)
            except Exception:
                continue
            if sw <= 0:
                continue
            # Indirect-dimension values are presented in the acquisition matrix
            # rather than repeated as key/value rows.
    debug.append(('NUS status', 'Parsed %d samples' % len(samples)))
    return rows, debug, maxima

def _acquisition_table(inst, maxima):
    """Build the acquisition-dimension table from the configured vpar state.

    vpar is authoritative for the final vendor-resolved NMRPipe dimensions.
    The report shows both the canonical axis variable (xN/yN/zN/aN) and its
    numeric value.  If a required canonical value cannot be obtained, keep the
    dimension visible and return diagnostic information rather than silently
    dropping it from the table.
    """
    debug = []
    if inst is None:
        return [], [], [('XYZA size status', 'vpar instance unavailable')]

    labels = list(getattr(inst, 'labb', []) or [])
    try:
        ndim = int(str(getattr(inst, 'dim', len(labels) or 1))[0])
    except Exception:
        ndim = len(labels) or 1
    ndim = max(1, min(ndim, 4))

    try:
        canonical_sizes = inst.GetTimeDomainSizes()
    except Exception as exc:
        canonical_sizes = {axis: getattr(inst, axis + 'N', None)
                           for axis in ('x', 'y', 'z', 'a')}
        debug.append(('XYZA accessor status', 'GetTimeDomainSizes failed: %s' % exc))

    try:
        canonical_obs = inst.GetObservationFrequencies()
    except Exception as exc:
        canonical_obs = {axis: getattr(inst, axis + 'OBS', None)
                         for axis in ('x', 'y', 'z', 'a')}
        debug.append(('OBS accessor status', 'GetObservationFrequencies failed: %s' % exc))

    axes = ('x', 'y', 'z', 'a')
    sw_attrs = ('sw', 'sw1', 'sw2', 'sw3')
    columns = []
    axis_names = []
    size_values = []
    obs_values = []
    sweep = []
    sweep_ppm = []
    points = []
    times = []

    for i in range(ndim):
        axis = axes[i]
        axis_name = axis + 'N'
        columns.append(labels[i] if i < len(labels) and labels[i] else 'Dimension %d' % (i + 1))
        axis_names.append(axis_name)

        raw_n = canonical_sizes.get(axis)
        raw_obs = canonical_obs.get(axis)
        raw_sw = getattr(inst, sw_attrs[i], None)

        try:
            n = float(raw_n)
            if n <= 0:
                raise ValueError('non-positive value')
            n_text = str(int(n)) if n.is_integer() else '{:.6g}'.format(n)
            size_values.append(n_text)
            points.append(n_text)
        except Exception:
            n = None
            size_values.append('unavailable')
            points.append('unavailable')
            debug.append(('%s status' % axis_name,
                          'Could not resolve %s; canonical value=%r' % (axis_name, raw_n)))

        try:
            obs = float(raw_obs)
            if obs <= 0:
                raise ValueError('non-positive value')
            obs_values.append('{:.6g}'.format(obs))
        except Exception:
            obs = None
            obs_values.append('unavailable')
            debug.append(('%s OBS status' % axis_name,
                          'Could not resolve %sOBS; canonical value=%r' % (axis, raw_obs)))

        try:
            sw = float(raw_sw)
            if sw <= 0:
                raise ValueError('non-positive value')
            sweep.append('{:.6g}'.format(sw))
        except Exception:
            sw = None
            sweep.append('unavailable')
            debug.append(('%s sweep-width status' % axis_name,
                          'Could not resolve %s from vpar.%s; value=%r' %
                          (axis_name, sw_attrs[i], raw_sw)))

        if sw is not None and obs is not None:
            sweep_ppm.append('{:.6g}'.format(sw / obs))
        else:
            sweep_ppm.append('unavailable')

        if sw is not None and n is not None:
            if i > 0 and i - 1 < len(maxima or []):
                try:
                    coord = float(maxima[i - 1])
                except Exception:
                    coord = n
            else:
                coord = n
            times.append('{:.3f}'.format(1000.0 * coord / sw))
        else:
            times.append('unavailable')

    rows = [('OBS (MHz)', obs_values),
            ('Sweep width (Hz)', sweep),
            ('Sweep width (ppm)', sweep_ppm),
            ('Time points', points),
            ('Acquisition time (ms)', times)]
    return columns, rows, debug

def _format_fid_selection(proc, frame):
    """Return the current Process FID selection as ``selected of available``.

    The Process GUI exposes 2^(n-1) phase/FID choices for an n-D data set
    (2 for 2D, 4 for 3D, etc.).
    """
    try:
        selected = proc._current_fid_selection() if hasattr(proc, '_current_fid_selection') else getattr(proc, 'FIDsel', '')
    except Exception:
        selected = getattr(proc, 'FIDsel', '')
    try:
        ndim = int(getattr(frame, 'dim', 0) or getattr(getattr(frame, 'state', None), 'dimension', 0) or 1)
        total = 2 ** max(ndim - 1, 0)
    except Exception:
        total = None
    selected = str(selected).strip()
    return ('%s of %d' % (selected, total)) if selected and total else selected


def _load_pseudo_axis_tsv(frame):
    """Read the project-owned pseudo_axis.tsv verbatim enough for report display."""
    state = getattr(frame, 'state', None)
    candidates = []
    if state is not None:
        try: candidates.append(Path(state.spec_dir()) / 'pseudo_axis.tsv')
        except Exception: pass
        spec = str(getattr(state, 'spec_path', '') or '').strip()
        if spec: candidates.append(Path(spec) / 'pseudo_axis.tsv')
    for path in candidates:
        try:
            if not path.is_absolute():
                base = Path(getattr(state, 'working_dir', '') or os.getcwd())
                path = base / path
            if not path.is_file(): continue
            with path.open(newline='', encoding='utf-8', errors='replace') as handle:
                rows = list(csv.reader(handle, delimiter='\t'))
            if rows:
                return rows[0], rows[1:], path
        except Exception:
            continue
    return [], [], None


def _combine_report_rows(frame):
    """Return a compact table of raw sources recorded by Process -> Combine.

    Paths are displayed relative to the project working directory so reports do
    not expose or waste space on the absolute project prefix.
    """
    state = getattr(frame, 'state', None)
    raw = str(getattr(state, 'raw_path', '') or '').strip() if state is not None else ''
    if not raw: return []
    work = Path(getattr(state, 'working_dir', '') or os.getcwd()).expanduser().resolve()
    path = Path(raw).expanduser()
    if not path.is_absolute(): path = work / path
    manifest = path / 'combine_manifest.json'
    if not manifest.is_file(): return []
    try:
        payload = json.loads(manifest.read_text(encoding='utf-8'))
        rows = []
        for source in payload.get('sources') or []:
            directory = str(source.get('directory', '') or '')
            if directory:
                try:
                    dp = Path(directory).expanduser()
                    if dp.is_absolute(): directory = os.path.relpath(str(dp), str(work))
                except Exception:
                    pass
            row_start, row_end = source.get('row_start'), source.get('row_end')
            span = ('%s-%s' % (row_start, row_end)) if row_start is not None and row_end is not None else ''
            rows.append([str(source.get('experiment', '?')), directory, span])
        return rows
    except Exception:
        return []

def _process_report_data(frame, report_dir, warnings):
    """Use the hidden Process/Conversion GUI to collect metadata and figures."""
    result = {'rows': [], 'combine_rows': [], 'processing_rows': [], 'spectrum_rows': [], 'debug_rows': [], 'acquisition_columns': [], 'acquisition_rows': [], 'pseudo_axis_columns': [], 'pseudo_axis_rows': [], 'pseudo_axis_explanation': '', 'fid_path': ''}
    proc = None
    conv = None
    try:
        from spinDecon.gui.dialogs.processing.process import ProcessFrame
        proc = ProcessFrame(None, 20, 'Process', frame, showFlg=False)
        try:
            proc.UpdateLampLights()
        except Exception:
            pass
        fmt = ''
        try:
            fmt = proc.formatStatusValue.GetLabel().strip()
        except Exception:
            try:
                fmt = proc._format_status_text()
            except Exception:
                pass
        if fmt:
            result['rows'].append(('Spectrometer format', fmt))

        # Raw/FID and spectrum directories are project-owned state.  Process
        # consumes these values but no longer owns GUI controls for them.
        state = getattr(frame, 'state', None)
        if state is not None:
            fid_path = str(getattr(state, 'raw_path', '') or '').strip()
            spectrum_path = str(getattr(state, 'spec_path', '') or '').strip()
        else:
            fid_path = _ctrl_value(frame, 'outPathBox')
            spectrum_path = _ctrl_value(frame, 'specPathBox')
        if fid_path:
            fid = Path(fid_path).expanduser()
            if not fid.is_absolute():
                working = Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()
                fid = working / fid
            result['fid_path'] = str(fid.resolve())
        result['combine_rows'] = _combine_report_rows(frame)

        # Projection reduction is a Process-window setting.  Use the live
        # combobox value so the summary reflects the current GUI state; the
        # control itself is restored from the saved `projectionType` parameter
        # by ProcessFrame.set_default_values().
        try:
            projection_type = str(proc.projections.GetValue() or 'sum').strip().lower()
        except Exception:
            projection_type = 'sum'
        if projection_type not in ('sum', 'skyline'):
            projection_type = 'sum'
        result['processing_rows'].append(('Projection type', projection_type))

        # Script target is a processing setting, so report it alongside the
        # other processing controls immediately after Projection type.
        params, _ = _read_parameter_map(frame)
        target = str(_param_ci(params, 'ProcTarg', default='')).strip()
        if target.lower() in ('smile', 'mddnmr'):
            result['processing_rows'].append(('Script target', target.upper() if target.lower() == 'mddnmr' else 'SMILE'))

        # These are processing/display controls, not acquisition metadata.
        # Keep them immediately after Projection type in the Processing section.
        result['processing_rows'].extend([
            ('X range minimum (ppm)', _ctrl_value(proc, 'xminBox')),
            ('X range maximum (ppm)', _ctrl_value(proc, 'xmaxBox')),
            ('FID select', _format_fid_selection(proc, frame)),
        ])

        # Save exactly the two views displayed by the Process window.
        old_fid = bool(proc.cb_show_fid.GetValue()) if hasattr(proc, 'cb_show_fid') else False
        if hasattr(proc, 'cb_show_fid'):
            proc.cb_show_fid.SetValue(False)
        proc.draw_figure()
        proc.canvas.print_figure(str(report_dir / 'process_phased.pdf'))
        if hasattr(proc, 'cb_show_fid'):
            proc.cb_show_fid.SetValue(True)
            proc.draw_figure(reset_y=True)
            proc.canvas.print_figure(str(report_dir / 'process_fid.pdf'))
            proc.cb_show_fid.SetValue(old_fid)
            proc.draw_figure()

        # Export the phasing view from the Process -> Projections window.
        # Construct the same window used interactively by OnButtonProjections;
        # its constructor loads the current projection/phase state and performs
        # a full redraw, so the saved canvas is exactly the current GUI view.
        phase_frame = None
        try:
            from spinDecon.gui.dialogs.processing.projections import ProjectionsFrame
            phase_frame = ProjectionsFrame(proc)
            try:
                phase_frame._full_redraw()
            except Exception:
                # load_projections() is already called by the constructor; a
                # second redraw is only to make the report intent explicit.
                pass
            phase_frame.canvas.print_figure(str(report_dir / 'phasing.pdf'))
        except Exception as exc:
            warnings.append('Phasing figure unavailable: %s' % exc)
        finally:
            if phase_frame is not None:
                try: phase_frame.Destroy()
                except Exception: pass

        # Conversion's Show Script window obtains these fields from vpar.
        try:
            from spinDecon.gui.dialogs.processing.conversion import ConversionFrame
            conv = ConversionFrame(proc)
            params, _ = _read_parameter_map(frame)
            nus_path, nus_saved, fid_saved = _resolve_nus_file(frame, params)
            # ConversionFrame._build_vpar() passes the live nusFil control to
            # vpar.Setup().  Populate it from the saved system state when the
            # hidden Conversion window did not recover it itself.
            if nus_saved and hasattr(conv, 'nusFil') and not conv.nusFil.GetValue().strip():
                conv.nusFil.SetValue(nus_saved)
            inst = conv._build_vpar()
            # Resolve the native Varian/Bruker acquisition parameters without
            # writing or executing a conversion script.  vpar.Convert() now
            # publishes the final vendor-resolved xN/yN/zN/aN sizes for the
            # acquisition table while leaving the established conversion logic
            # untouched.
            inst.Convert()

            # For the canonical 2-spectral + 1-pseudo topology, report the
            # exact pseudo-axis table generated by vpar_decon.  pseudo_axis_info
            # is populated by the same vendor-specific detection used to write
            # pseudo_axis.tsv, so the report cannot drift from conversion.
            try:
                from spinDecon.domain.analysis_mode import AnalysisMode
                mode = AnalysisMode.from_project_state(frame.state)
                info = getattr(inst, 'pseudo_axis_info', None)
                if int(mode.spectral_dimensions) == 2 and bool(mode.has_pseudo_axis) and info and getattr(info, 'rows', None):
                    result['pseudo_axis_columns'] = ['spectrum'] + list(info.columns)
                    result['pseudo_axis_rows'] = [[str(i)] + [str(v) for v in row]
                                                  for i, row in enumerate(info.rows, start=1)]
                    groups = getattr(info, 'groups', []) or []
                    group_text = ', '.join('(' + ', '.join(map(str, group)) + ')' for group in groups)
                    if getattr(inst, 'tp', None) == 'bruk':
                        result['pseudo_axis_explanation'] = (
                            'The pseudo axis is reconstructed from the Bruker pulseprogram and its array-list files. '
                            'Parameters incremented within the same lo ... times loop are treated as synchronized '
                            'columns; independent or nested increment loops are combined as a Cartesian product. '
                            'The rows below are therefore the same ordered pseudo spectra written to pseudo_axis.tsv.'
                        )
                    elif getattr(inst, 'tp', None) == 'var':
                        result['pseudo_axis_explanation'] = (
                            'The pseudo axis is reconstructed from the Varian/VNMR procpar array expression. '
                            'Parameters grouped in parentheses are synchronized, independent array parameters are '
                            'combined as a Cartesian product, and phase/phase2/phase3 quadrature arrays are excluded. '
                            'The rows below are therefore the same ordered pseudo spectra written to pseudo_axis.tsv.'
                        )
                    else:
                        result['pseudo_axis_explanation'] = (
                            'The pseudo-axis rows are taken from the vendor-specific array detection used by vpar_decon '
                            'to generate pseudo_axis.tsv; synchronized parameters are zipped and independent arrays '
                            'are combined as a Cartesian product.'
                        )
                    if group_text:
                        result['pseudo_axis_explanation'] += ' Detected array group(s): ' + group_text + '.'
            except Exception as exc:
                result['debug_rows'].append(('Pseudo-axis table status', str(exc)))

            # Always include the project pseudo_axis.tsv for any pseudo-dimensional
            # topology, including one spectral + one pseudo (pseudo2D).  The
            # conversion-derived table above remains preferred when available.
            if not result['pseudo_axis_rows']:
                cols, rows, _axis_path = _load_pseudo_axis_tsv(frame)
                if rows:
                    result['pseudo_axis_columns'] = list(cols)
                    result['pseudo_axis_rows'] = [list(row) for row in rows]
                    result['pseudo_axis_explanation'] = (
                        'This is the pseudo_axis.tsv stored with the processed spectrum; '
                        'it is the authoritative mapping between pseudo-axis rows and raw experiments.'
                    )

            pulse = str(getattr(inst, 'seqfil', '') or '').strip()
            if pulse and pulse not in ('0', 'None'):
                result['rows'].append(('Pulse sequence', pulse))
            bfmt = str(getattr(inst, 'brukerFmt', '') or '').strip()
            if bfmt and bfmt not in ('0', 'None'):
                result['rows'].append(('Bruker format', bfmt))
            # Acquisition controls are exposed by vpar for both vendors:
            # Bruker NS / D1 and Varian nt / d1.
            scans = getattr(inst, 'ns', None)
            if scans not in (None, '', '0', 0):
                try:
                    scans_text = str(int(float(scans)))
                except Exception:
                    scans_text = str(scans).strip()
                if scans_text:
                    result['rows'].append(('Scans per transient', scans_text))
            recycle = getattr(inst, 'd1', None)
            if recycle not in (None, ''):
                try:
                    recycle_text = ('%.6g s' % float(recycle))
                except Exception:
                    recycle_text = str(recycle).strip()
                    if recycle_text:
                        recycle_text += ' s'
                if recycle_text:
                    result['rows'].append(('Recycle delay', recycle_text))
            try:
                measurement = str(inst.GetMeasurementTime() or '').strip()
            except Exception as exc:
                measurement = ''
                result['debug_rows'].append(('Measurement time status', str(exc)))
            if measurement:
                result['rows'].append(('Measurement time', measurement))
            else:
                result['debug_rows'].append(('Measurement time status', 'No usable spectrometer audit timing found'))
        except Exception as exc:
            warnings.append('Conversion metadata unavailable: %s' % exc)

        # params was loaded above for Processing metadata.
        nus_rows, nus_debug, nus_maxima = _nus_summary(frame, params, inst=inst if 'inst' in locals() else None)
        result['rows'].extend(nus_rows)
        result['debug_rows'].extend(nus_debug)
        acq_columns, acq_rows, acq_debug = _acquisition_table(inst if 'inst' in locals() else None, nus_maxima)
        result['acquisition_columns'] = acq_columns
        result['acquisition_rows'] = acq_rows
        result['debug_rows'].extend(acq_debug)
    except Exception as exc:
        warnings.append('Process-window summary unavailable: %s' % exc)
    finally:
        if conv is not None:
            try: conv.Destroy()
            except Exception: pass
        if proc is not None:
            try: proc.Destroy()
            except Exception: pass

    # Describe the actual loaded NMRPipe spectrum.
    try:
        import nmrglue as ng
        data = getattr(frame, 'data', None)
        dic = getattr(frame, 'dic', None)
        if data is not None and dic is not None:
            ndim = int(getattr(data, 'ndim', len(getattr(data, 'shape', []))))
            labels = _dimension_labels(frame, ndim=ndim)
            for i in range(ndim):
                uc = ng.pipe.make_uc(dic, data, dim=i)
                a, b = uc.ppm_limits()
                result['spectrum_rows'].append((labels[i] if i < len(labels) else 'Dimension %d' % (i + 1),
                                                [format(min(a,b), '.6g'), format(max(a,b), '.6g'), str(data.shape[i])]))
    except Exception as exc:
        warnings.append('Spectrum dimension metadata unavailable: %s' % exc)
    # Conversion metadata can fail independently of an otherwise valid project.
    # The saved TSV is still authoritative and should always be represented.
    if not result['pseudo_axis_rows']:
        cols, rows, _axis_path = _load_pseudo_axis_tsv(frame)
        if rows:
            result['pseudo_axis_columns'] = list(cols)
            result['pseudo_axis_rows'] = [list(row) for row in rows]
            result['pseudo_axis_explanation'] = (
                'This is the pseudo_axis.tsv stored with the processed spectrum; '
                'it is the authoritative mapping between pseudo-axis rows and raw experiments.'
            )
    return result


def _collect_rows(frame):
    state = getattr(frame, 'state', None)
    project = []
    if state is not None:
        # Keep this section concise and source file names from the current
        # project state wherever that state owns the value.  The full nD peak
        # list is currently GUI-owned, so use its live control value.
        # Spectrum-associated controls are canonically relative to SpecPath.
        # Report that canonical value rather than re-prefixing it with SpecPath.
        spectrum_name = _relative_project_path(frame, getattr(state, 'input_file', ''))
        reference_name = _relative_project_path(frame, getattr(state, 'reference_peak_file', getattr(state, 'peak_file', '')))
        full_name = _relative_project_path(frame, _ctrl_value(frame, 'fullPeakBox') or getattr(state, 'full_peak_file', ''))
        working_dir = str(Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser().resolve())
        session_name = _relative_project_path(frame, getattr(state, 'session_file', ''))
        spec_value = str(getattr(state, 'spec_path', '') or '').strip()
        spec_path = os.path.relpath(str(Path(state.spec_dir()).expanduser().resolve()), working_dir) if spec_value else ''
        project.extend([
            ('Session file', session_name),
            ('Working dir', working_dir),
            ('SpecPath', spec_path),
            ('Spectrum', spectrum_name),
            ('Reference peaklist', reference_name),
            ('Full peaklist', full_name),
            ('Dimensions', state.dimension),
            ('Pseudo axis', 'yes' if state.pseudo_axis else 'no'),
        ])
        for key, value in (getattr(state, 'metadata', {}) or {}).items():
            if value not in ('', None) and key not in (
                    'parameter_file', 'input_file', 'peak_file', 'session_file',
                    'full_peak_file', 'working_dir', 'spec_path', 'fidsel', 'FIDsel', 'projection_view_mode',
                    'projection_phasing', 'peakframe_projection_key', 'peakframe_peak_count'):
                display_key = str(key).replace('_', ' ').strip().lower()
                if display_key in ('fidsel', 'fid sel', 'fid selection', 'projection view mode', 'projection phasing', 'peakframe projection key', 'peakframe peak count'):
                    continue
                project.append((str(key).replace('_', ' ').title(), value))
    try:
        fit_value = 'yes' if bool(frame.cb_decback.GetValue()) else 'no'
    except Exception:
        fit_value = 'no'
    decon = [('CPUs', _ctrl_value(frame, 'coreBox')), ('Factor', _ctrl_value(frame, 'facBox')),
             ('Convergence', _ctrl_value(frame, 'convBox')), ('Maximum iterations', _ctrl_value(frame, 'maxiterBox')),
             ('Fit', fit_value),
             ('Impose symmetry', 'yes' if (state is not None and state.sym_mode) else 'no'),
             ('Use 2D peak list', 'yes' if (state is not None and state.decon_bore) else 'no')]
    processing = _processing_table(frame)
    try:
        ndim = int(getattr(frame, 'dim', 0) or getattr(state, 'dimension', 0) or 1)
    except Exception:
        ndim = 1
    labels = _dimension_labels(frame, ndim=ndim)
    shape_rows=[]
    for title, prefix in [('Sigma','sig'), ('Lorentz','lorentz'), ('Voigt','voigt')]:
        vals=[_ctrl_value(frame, '%s%dBox' % (prefix, i+1)) for i in range(ndim)]
        if any(v not in ('', None) for v in vals):
            shape_rows.append((title, vals))
    return project, processing, decon, (labels, shape_rows)



def _read_processing_scripts(frame):
    """Return the current NMRPipe scripts from the canonical SpecPath.

    SpecPath is rooted at WorkingDir and owned by ProjectState/NMR.  Do not
    route the directory itself through ``_project_path`` because that helper
    resolves *files inside* SpecPath.
    """
    state = getattr(frame, 'state', None)
    if state is not None:
        try:
            base = Path(state.spec_dir()).expanduser()
        except Exception:
            base = None
    else:
        base = None
    if base is None:
        working = Path(_ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()
        spec = Path(_ctrl_value(frame, 'specPathBox') or './spec').expanduser()
        base = spec if spec.is_absolute() else working / spec

    scripts = []
    for label, filename in (('Conversion script', 'fid.test.com'),
                            ('Processing script', 'nmrproc.test.com')):
        path = base / filename
        try:
            text = path.read_text(encoding='utf-8', errors='replace')
        except Exception:
            continue
        scripts.append((label, filename, text))
    return scripts


def _write_script_verbatim(out, title, filename, text):
    """Write a shell/NMRPipe script literally in tiny LaTeX verbatim text."""
    # LaTeX commands need a single backslash in the generated .tex file,
    # whereas line breaks here must be real newline characters (not the two
    # literal characters ``\\n``).
    out.write('\\subsection*{' + _latex_escape(title) + ': ' + _latex_escape(filename) + '}\n')
    out.write('\\begingroup\\tiny\n')
    out.write('\\begin{verbatim}\n')
    # A literal end-verbatim token in a shell comment would terminate the
    # environment. It is not expected in generated NMRPipe scripts, but keep
    # the historical defensive guard.
    text = str(text).replace('\\end{verbatim}', '\\end{verbatim} ')
    out.write(text)
    if text and not text.endswith('\n'):
        out.write('\n')
    out.write('\\end{verbatim}\n')
    out.write('\\endgroup\n\n')

def _workflow_report_data(frame, warnings):
    """Evaluate the same workflow model used by the Workflow page."""
    try:
        from spinDecon.domain.analysis_mode import AnalysisMode
        from spinDecon.workflow.model import build_workflow_plan
        from spinDecon.workflow.status import evaluate_workflow, StageStatus
        mode = AnalysisMode.from_project_state(frame.state)
        plan = build_workflow_plan(mode)
        notebook = getattr(frame, 'parent', None)
        store = getattr(frame, 'store', None) or getattr(notebook, 'data_store', None)
        states = evaluate_workflow(plan, frame.state, store, notebook)
        by_key = {item.key: item for item in states}

        # In the physical 2D + pseudo workflow, ``analyse_series`` is the
        # terminal guided-workflow decision: selecting an analysis type marks
        # the project workflow complete.  The interactive Workflow page keeps
        # the review stage independently actionable so that revisiting it can
        # invalidate an old inspection.  A project summary, however, is a
        # snapshot of the completed workflow.  Do not report the prerequisite
        # Review Intensity Series stage as failed when the terminal analysis
        # stage is already complete.
        terminal_2d_pseudo = bool(
            mode.spectral_dimensions == 2
            and mode.has_pseudo_axis
            and 'analyse_series' in by_key
            and by_key['analyse_series'].status is StageStatus.COMPLETE
        )

        rows = []
        for stage in plan.stages:
            item = by_key[stage.key]
            complete = item.status is StageStatus.COMPLETE
            status = item.status.value
            detail = item.detail
            if terminal_2d_pseudo and stage.key == 'review_series' and not complete:
                complete = True
                status = StageStatus.COMPLETE.value
                detail = 'The intensity-series review is complete as part of the completed 2D + pseudo workflow.'
            rows.append((stage.title, complete, status, detail))

        return {
            'pseudo': bool(mode.has_pseudo_axis),
            'rows': rows,
        }
    except Exception as exc:
        warnings.append('Workflow summary unavailable: %s' % exc)
        return {'pseudo': False, 'rows': []}


def _pseudo_report_data(frame, report_dir, warnings, workflow, progress_callback=None):
    """Collect fitting tables/figures through the specialist GUI APIs.

    The physical-2D fitting workflow deliberately reuses the Pseudo3D/Fitting
    workspace as a single-plane fit.  Treat that workspace as the report
    authority too, rather than maintaining a second plotting/parser path.
    """
    result = {'columns': [], 'rows': [], 'units': [], 'figures': [], 'analysis': None, 'kind': 'pseudo3d'}
    try:
        from spinDecon.domain.analysis_mode import AnalysisMode
        mode = AnalysisMode.from_project_state(frame.state)
        if int(mode.spectral_dimensions) == 1 and bool(mode.has_pseudo_axis):
            result['kind'] = 'pseudo2d'
        elif int(mode.spectral_dimensions) == 2 and not bool(mode.has_pseudo_axis):
            result['kind'] = '2d'
    except Exception:
        pass
    is_physical_2d = result.get('kind') == '2d'
    notebook = getattr(frame, 'parent', None)
    getter = getattr(notebook, 'get_page_by_title', None)
    existing_fitting = getter('Fitting') if callable(getter) else getattr(notebook, 'tabPseudo', None)
    fitting_ready = bool(getattr(getattr(frame, 'store', None), 'analysis', {}).get('fitting_results_ready'))
    if not workflow.get('pseudo') and not (is_physical_2d and (fitting_ready or existing_fitting is not None)):
        return result
    pseudo = None
    try:
        getter = getattr(notebook, 'get_page_by_title', None)
        if callable(getter):
            pseudo = getter('Fitting') or getter('Pseudo2D')
        if pseudo is None and is_physical_2d:
            add = getattr(notebook, 'AddTabPseudo3D', None)
            if callable(add):
                add(True, frame)
                if callable(getter):
                    pseudo = getter('Fitting')
        if pseudo is None:
            pseudo = getattr(notebook, 'tabPseudo', None)
        if pseudo is None:
            raise RuntimeError('pseudo-dimensional fitting workspace is not available')
        # Refresh the actual SpinUniDec Fitting window silently, then read its
        # displayed rows/groups.  Project-memory grouping is not authoritative.
        result['columns'], result['rows'], result['units'] = pseudo.fitting_window_report_data()
        if callable(progress_callback):
            if result.get('kind') == 'pseudo2d':
                progress_callback('Saving pseudo2D group and slice-fit figures')
            elif result.get('kind') == '2d':
                progress_callback('Saving 2D group fitting figures')
            else:
                progress_callback('Saving pseudo3D reference and fitting figures')
        result['figures'] = pseudo.export_fitting_report_figures(report_dir, result['units'])
    except Exception as exc:
        warnings.append('%s fitting report unavailable: %s' % (result['kind'], exc))
        return result

    # pseudo2D currently has no downstream Analysis mode.  Keep the report
    # structure ready for it, but do not manufacture an analysis section.
    if result.get('kind') in ('pseudo2d', '2d'):
        return result

    if callable(progress_callback):
        progress_callback('Collecting pseudo3D analysis results')
    selected = ''
    try:
        selected = pseudo.selected_downstream_analysis()
    except Exception:
        pass
    if selected == 'Decay':
        decay = None
        try:
            from spinDecon.gui.workspaces.decay import DecayFrame
            decay = DecayFrame(None, 30, 'Decay Analysis', pseudo, pth='')
            decay.Hide()
            decay.OnButtonPeakConvert(None)
            result['analysis'] = {'name': 'Decay', **decay.export_report_figures(report_dir)}
        except Exception as exc:
            warnings.append('Decay analysis report unavailable: %s' % exc)
        finally:
            if decay is not None:
                try: decay.Destroy()
                except Exception: pass
    elif selected == 'CPMG':
        cpmg = None
        try:
            from spinDecon.gui.workspaces.cpmg import CPMGFrame
            cpmg = CPMGFrame(None, 30, 'CPMG Analysis', pseudo, pth='')
            cpmg.Hide()
            result['analysis'] = {'name': 'CPMG', **cpmg.export_report_figures(report_dir)}
        except Exception as exc:
            warnings.append('CPMG analysis report unavailable: %s' % exc)
        finally:
            if cpmg is not None:
                try: cpmg.Destroy()
                except Exception: pass
    elif selected:
        warnings.append('%s analysis is selected; detailed report export is not yet provided by that analysis frame.' % selected)
    return result


def _write_workflow_table(out, workflow):
    rows = workflow.get('rows', [])
    if not rows:
        return
    out.write('\\section*{Workflow status}\n')
    out.write('\\begin{tabular}{@{}p{0.07\\linewidth}p{0.30\\linewidth}p{0.57\\linewidth}@{}}\n')
    out.write('\\toprule Status & Stage & Detail\\\\ \\midrule\n')
    for title, complete, status, detail in rows:
        mark = (r'\textcolor{green!55!black}{\Large $\boldsymbol{\checkmark}$}' if complete else r'\textcolor{red!70!black}{\Large $\boldsymbol{\times}$}')
        out.write(mark + ' & ' + _latex_escape(title) + ' & ' + _latex_escape(detail) + r'\\' + '\n')
    out.write('\\bottomrule\\end{tabular}\n\n')


def _fit_value_error_text(value, error):
    """Format a value/error to the decimal place supported by its 1-sigma error."""
    try:
        value=float(value)
        if error in (None, '', '-'):
            return str(value)
        error=float(error)
        if not (math.isfinite(value) and math.isfinite(error)) or error <= 0:
            return '%.4g' % value
        exponent=int(math.floor(math.log10(abs(error))))
        first=abs(error)/(10.0**exponent)
        sig=2 if first < 3.0 else 1
        decimals=max(0, -exponent + sig - 1)
        if decimals <= 6 and abs(value) < 1e7:
            fmt='%%.%df' % decimals
            return (fmt % value) + r' $\pm$ ' + (fmt % error)
        return ('%.4g' % value) + r' $\pm$ ' + ('%.*g' % (sig,error))
    except (TypeError, ValueError, OverflowError):
        return '-'


def _write_analysis_summary(out, pseudo):
    analysis = pseudo.get('analysis') or {}
    figures = analysis.get('summary_figures', [])
    screen = analysis.get('screen') or {}
    global_fit = analysis.get('global') or {}
    if not figures and not screen and not global_fit:
        return
    out.write('\\section*{Analysis summary}\n')
    if analysis.get('name') == 'CPMG':
        out.write('\\textbf{CPMG Rex screen:} threshold %s s$^{-1}$; %s of %s peaks significant.\\\\\n' %
                  (_latex_escape(screen.get('Rex_threshold','-')), _latex_escape(screen.get('n_significant','-')),
                   _latex_escape(screen.get('n_total','-'))))
        if global_fit.get('success'):
            kex=_fit_value_error_text(global_fit.get('kex'),global_fit.get('kex_error'))
            pb=_fit_value_error_text(global_fit.get('pb'),global_fit.get('pb_error'))
            out.write('\\textbf{Global fit:} $k_{ex}$=%s s$^{-1}$, $p_b$=%s, $\\chi^2$=%s (%s peaks).\\par\\medskip\n' %
                      (kex, pb, _latex_escape(global_fit.get('chi2','-')), _latex_escape(global_fit.get('n_peaks','-'))))
        elif screen.get('n_significant'):
            out.write('\\textit{No successful global CPMG fit was available.}\\par\\medskip\n')
    for index, (filename, title) in enumerate(figures):
        if index and index % 2 == 0:
            out.write('\\par\\medskip\n')
        out.write('\\begin{minipage}[t]{0.49\\linewidth}\\centering\n')
        out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(filename) + '}\\\\\n')
        out.write('\\small ' + _latex_escape(title) + '\n\\end{minipage}')
        if index % 2 == 0 and index + 1 < len(figures):
            out.write('\\hfill\n')
    out.write('\\par\\medskip\n\n')
    # Keep the numerical CPMG summary with the analysis, after its figures.
    if analysis.get('name') == 'CPMG':
        _write_cpmg_parameter_tables(out, analysis, section_level='subsection')


def _write_cpmg_parameter_tables(out, analysis, section_level='section'):
    """Write CPMG results separately from the SpinUniDec peak-fit table."""
    columns=list(analysis.get('columns') or [])
    rows=[list(r) for r in analysis.get('rows') or []]
    if not columns or not rows:
        return
    idx={name:i for i,name in enumerate(columns)}
    def val(row,name):
        i=idx.get(name); return row[i] if i is not None and i < len(row) else ''
    def pm(row,name,err):
        v=val(row,name); e=val(row,err)
        return ('%s +/- %s' % (v,e)) if e not in ('', '-', None) else v

    out.write('\\%s*{CPMG parameter summary}\n' % section_level)
    out.write('\\noindent\\textit{Significant peaks are ranked by the local-fit Rex screening value. Errors are approximate 1$\\sigma$ standard errors from the least-squares Jacobian covariance.}\\par\\smallskip\n')
    local_cols=['Rank','Peak','Rex (s$^{-1}$)','$R_0$ (s$^{-1}$)','$\\Delta\\omega$ (ppm)','$k_{ex}$ (s$^{-1}$)','$p_b$','Fit gain']
    local_rows=[]
    for row in rows:
        local_rows.append([val(row,'Rank'),val(row,'Peak'),val(row,'Rex'),pm(row,'Local R0','Local R0 error'),
                           pm(row,'Local dw','Local dw error'),pm(row,'Local kex','Local kex error'),
                           pm(row,'Local pb','Local pb error'),val(row,'Local gain')])
    _write_long_table_with_split(out,local_cols,local_rows,len(local_cols),font='\\scriptsize',compact=True)

    heading='subsubsection' if section_level == 'subsection' else 'subsection'
    out.write('\\%s*{Global-fit peak parameters}\n' % heading)
    out.write('\\noindent\\scriptsize Shared $k_{ex}$ and $p_b$ (with their errors) are reported above; only peak-specific global parameters are repeated here.\\par\\smallskip\n')
    global_cols=['Rank','Peak','$R_0$ (s$^{-1}$)','$\\Delta\\omega$ (ppm)','$R_{2,\\infty}$ (s$^{-1}$)','Local $\\chi^2$','Global $\\chi^2$','Fit gain']
    global_rows=[]
    for row in rows:
        global_rows.append([val(row,'Rank'),val(row,'Peak'),pm(row,'Global R0','Global R0 error'),
                            pm(row,'Global dw','Global dw error'),val(row,'Global R2inf'),val(row,'Local chi2'),val(row,'Global chi2'),val(row,'Global gain')])
    _write_long_table_with_split(out,global_cols,global_rows,len(global_cols),font='\\scriptsize',compact=True)

def _joined_pseudo_rows(pseudo):
    """Join SpinUniDec fitting and downstream-analysis rows by peak name."""
    fit_cols = list(pseudo.get('columns') or [])
    fit_rows = [list(row) for row in (pseudo.get('rows') or [])]
    analysis = pseudo.get('analysis') or {}
    ana_cols = list(analysis.get('columns') or [])
    ana_rows = [list(row) for row in (analysis.get('rows') or [])]
    ana_by_peak = {str(row[0]): row for row in ana_rows if row}
    # Peak is the join key and is not repeated on the analysis side.
    joined_cols = fit_cols + ana_cols[1:]
    joined_rows = []
    for fit in fit_rows:
        peak = str(fit[0]) if fit else ''
        ana = ana_by_peak.get(peak, [])
        joined_rows.append(fit + (ana[1:] if ana else [''] * max(0, len(ana_cols) - 1)))
    return joined_cols, joined_rows, len(fit_cols)


def _write_joined_results_table(out, columns, rows, split_at):
    """Compact, page-breakable master fitting/analysis table."""
    if not columns or not rows:
        return
    # Keep the reduced presentation, but use longtable so arbitrarily long
    # fitting result sets can continue over page boundaries with a repeated header.
    _write_long_table_with_split(out, columns, rows, split_at, font='\\tiny', compact=True)


def _write_compact_peak_table(out, peak, columns, joined_by_peak):
    """Render one small-footprint peak summary card."""
    index = {name: i for i, name in enumerate(columns)}
    row = joined_by_peak.get(str(peak), [])
    if not row:
        return

    def get(name):
        i = index.get(name)
        return row[i] if i is not None and i < len(row) else ''

    out.write('\\textbf{' + _latex_escape(peak) + '}')
    perr = get('%err')
    if perr:
        out.write('\\hfill {\\scriptsize peak fit ' + _latex_escape(perr) + '\\%}')
    out.write('\\\\[-1mm]\n')
    out.write('\\begin{tabular}{@{}lrr@{}}\n')
    out.write(' & \\textbf{Dim 1} & \\textbf{Dim 2}\\\\\n')
    out.write('Position & ' + _latex_escape(get('f01(ppm)')) + ' & ' + _latex_escape(get('f02(ppm)')) + r'\\' + '\n')
    out.write('Width & ' + _latex_escape(get('w1(Hz)')) + ' & ' + _latex_escape(get('w2(Hz)')) + r'\\' + '\n')
    out.write('Shape & ' + _latex_escape(get('g1')) + ' & ' + _latex_escape(get('g2')) + r'\\' + '\n')
    out.write('\\end{tabular}\\\\[-0.5mm]\n')

    # Deliberately separate the SpinUniDec peak-shape result from the
    # downstream analysis result inside every peak card.
    out.write('\\noindent\\rule{\\linewidth}{0.35pt}\\vspace{0.6mm}\n')
    analysis_names = [c for c in columns if c in ('R', 'R error', 'A0', 'A0 error', 'Average error (%)')]
    if analysis_names:
        out.write('\\begin{tabular}{@{}lr@{}}\n')
        for name in analysis_names:
            label = {'R error':'R err.', 'A0 error':'A0 err.', 'Average error (%)':'Avg. err. (\\%)'}.get(name, name)
            out.write(label + ' & ' + _latex_escape(get(name)) + r'\\' + '\n')
        out.write('\\end{tabular}\n')


def _write_peak_summary_minipage(out, peak, columns, joined_by_peak, width='0.31'):
    """Write one peak card as a minipage suitable for a three-column array."""
    out.write('\\begin{minipage}[t]{' + width + '\\linewidth}\\vspace{0pt}\\scriptsize\n')
    _write_compact_peak_table(out, peak, columns, joined_by_peak)
    out.write('\\end{minipage}')

def _write_pseudo2d_report(out, pseudo):
    if not pseudo.get('rows') and not pseudo.get('figures'):
        return
    out.write('\\section*{Pseudo2D fitting results}\n')
    columns = list(pseudo.get('columns') or [])
    rows = [list(row) for row in pseudo.get('rows') or []]
    by_peak = {str(row[0]): row for row in rows if row}
    index = {name: i for i, name in enumerate(columns)}

    def value(row, name):
        i = index.get(name); return row[i] if i is not None and i < len(row) else ''

    for _filename, unit in pseudo.get('figures', []):
        peaks = unit.get('peaks', [])
        label = ('Group %s: %s' % (unit.get('group'), ', '.join(peaks))) if unit.get('group') is not None else (peaks[0] if peaks else 'Peak')
        out.write('\\subsection*{' + _latex_escape(label) + '}\n')
        # Compact textual summary derived from the exact Fitting-window rows.
        out.write('\\begingroup\\scriptsize\n')
        for peak in peaks:
            row = by_peak.get(str(peak), [])
            if not row: continue
            summary = [('%err', value(row, '%err')), ('Position (ppm)', value(row, 'f01(ppm)')),
                       ('Width (Hz)', value(row, 'w1(Hz)')), ('Shape', value(row, 'g1')),
                       ('Phase (deg)', value(row, 'Phase (deg)')), ('wD/wA', value(row, 'wD/wA'))]
            text = ', '.join('%s: %s' % (k, v) for k, v in summary if v not in ('', None))
            out.write('\\textbf{' + _latex_escape(peak) + '}: ' + _latex_escape(text) + r'\\' + '\n')
        out.write('\\endgroup\\smallskip\n')

        # One contour/3D snapshot represents the whole overlap group.
        overviews = unit.get('overview_figures', [])
        if overviews:
            _peak, overview = overviews[0]
            out.write('\\begin{center}\n')
            out.write('\\includegraphics[width=0.88\\linewidth]{./pdf/' + _latex_escape(overview) + '}\\\\[-1mm]\n')
            out.write('\\scriptsize Group overview -- pseudo axis fixed to exp\\_no\n')
            out.write('\\end{center}\n')

        # One group fit per pseudo slice.  Each plot marks every fitted member
        # and uses the complete x range present in the saved fitting data.
        slice_files = unit.get('slice_figures', [])
        if slice_files:
            out.write('\\noindent\\textit{Group slice fits}\\par\\smallskip\n')
            out.write('\\begin{center}\n')
            for i, slice_file in enumerate(slice_files):
                out.write('\\begin{minipage}[t]{0.19\\linewidth}\\centering\n')
                out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(slice_file) + '}\n')
                out.write('\\end{minipage}')
                if i + 1 < len(slice_files): out.write('\\hfill\n')
                if (i + 1) % 5 == 0 and i + 1 < len(slice_files): out.write('\\par\\smallskip\n')
            out.write('\n\\end{center}\n')

    if rows:
        out.write('\\section*{Fitted parameters and metrics}\n')
        _write_joined_results_table(out, columns, rows, len(columns))
        _write_pseudo2d_fitting_model(out)


def _write_pseudo2d_fitting_model(out):
    """Document the restrained pseudo2D model implemented by SpinUniDec.

    Keep this description beside the fitted-parameter table: it describes the
    actual Protocol2PFit/PeakFraction1D model rather than a generic peak fit.
    """
    out.write('\\subsection*{Peak fitting model}\n')
    out.write('\\small\n')
    out.write(
        'Pseudo2D fitting is a restrained, group-wise fit. Peaks whose extraction '
        'windows overlap are fitted simultaneously, while different pseudo-axis '
        'slices share the same resonance shape parameters. For each slice $s$, '
        'the calculated spectrum in an overlap group is a linear sum of the fitted '
        'resonances,\\\\\n'
    )
    out.write('\\[\n')
    out.write('I_s(x)=\\sum_p A_{p,s}\\left[P_A(x-\\delta_p)+r_p P_D(x-\\delta_p)\\right],\n')
    out.write('\\]\n')
    out.write(
        'where $A_{p,s}$ is the absorptive intensity of peak $p$ in slice $s$, '
        '$\\delta_p$ is its restrained chemical shift, and $r_p$ is the signed '
        'dispersive-to-absorptive amplitude ratio. The amplitudes are obtained by '
        'linear least squares independently for every pseudo slice, so overlapping '
        'peaks are resolved simultaneously rather than fitted one at a time.\n'
    )
    out.write('\\par\\smallskip\n')
    out.write(
        'The absorptive line shape is the unit-height pseudo-Voigt used by SpinUniDec,\n'
    )
    out.write('\\[\n')
    out.write('P_A(\\Delta)= (1-g)\\exp\\!\\left[-\\frac{\\Delta^2}{2\\sigma^2}\\right] + g\\frac{(\\Gamma/2)^2}{\\Delta^2+(\\Gamma/2)^2},\n')
    out.write('\\qquad \\sigma=w_A/2.355,\n')
    out.write('\\]\n')
    out.write(
        'with mixing fraction $g$ ($g=0$ Gaussian, $g=1$ Lorentzian). The fitted '
        'absorptive width scales the input Gaussian and Lorentzian widths together. '
        'When phase/distortion fitting is enabled, a quadrature partner $P_D$ is '
        'added: its Gaussian term is evaluated from Dawson\'s integral and its '
        'Lorentzian term is $\\left(\\Gamma_D/2\\right)\\Delta/'
        '[\\Delta^2+(\\Gamma_D/2)^2]$. The absorptive and dispersive components '
        'may have different fitted widths, but share the same pseudo-Voigt mixing '
        'fraction $g$; one signed ratio $r_p$ is shared by all pseudo slices for a '
        'given resonance.\n'
    )
    out.write('\\par\\smallskip\n')
    out.write(
        'The nonlinear shape parameters are refined by bounded coordinate searches '
        'to reduce the summed squared residual over all slices. Peak centres remain '
        'strongly restrained to the supplied Full 1D peak positions, and the fitting '
        'radius defines only the data/overlap region---it is not used as a linewidth. '
        'For fixed centres, widths, mixing and dispersive ratios, slice amplitudes are '
        're-solved exactly by linear least squares.\n'
    )
    out.write('\\normalsize\n')


def _write_pseudo_report(out, pseudo):
    if pseudo.get('kind') == 'pseudo2d':
        _write_pseudo2d_report(out, pseudo)
        return
    if not pseudo.get('rows') and not pseudo.get('figures') and not pseudo.get('analysis'):
        return
    heading = '2D fitting results' if pseudo.get('kind') == '2d' else 'Pseudo3D fitting results'
    out.write('\\section*{' + heading + '}\n')
    columns, joined_rows, split_at = _joined_pseudo_rows(pseudo)
    joined_by_peak = {str(row[0]): row for row in joined_rows if row}

    # Each SpinUniDec unit is a self-contained result: concise numerical
    # summary on the left, canonical two-pane fit on the right, then its
    # downstream per-peak fits directly underneath.
    analysis = pseudo.get('analysis') or {}
    peak_figures = analysis.get('peak_figures', {})
    for filename, unit in pseudo.get('figures', []):
        label = ('Group %s: %s' % (unit['group'], ', '.join(unit['peaks']))) if unit.get('group') is not None else unit['peaks'][0]
        out.write('\\subsection*{' + _latex_escape(label) + '}\n')
        peaks = [peak for peak in unit.get('peaks', []) if str(peak) in joined_by_peak]
        # Peak summaries form a compact three-column array.  When the final
        # row contains exactly one card, use the otherwise empty two columns
        # for the Pseudo3D plot; in all other cases put the plot on its own
        # line so neither the cards nor the figure are cramped.
        # A single-peak fitting unit keeps the original compact two-column
        # composition: data card on the left and Pseudo3D fit immediately to
        # its right.  The three-card grid rules are only for overlap groups.
        if len(peaks) == 1:
            out.write('\\noindent\n')
            _write_peak_summary_minipage(out, peaks[0], columns, joined_by_peak, width='0.30')
            out.write('\\hfill\n')
            out.write('\\begin{minipage}[t]{0.67\\linewidth}\\vspace{0pt}\\centering\n')
            out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(filename) + '}\n')
            out.write('\\end{minipage}\\par\\smallskip\n')
            plot_written = True
            full_rows = remainder = cursor = 0
        else:
            full_rows = len(peaks) // 3
            remainder = len(peaks) % 3
            cursor = 0
            plot_written = False
        for _ in range(full_rows):
            out.write('\\noindent\n')
            for col in range(3):
                _write_peak_summary_minipage(out, peaks[cursor], columns, joined_by_peak)
                cursor += 1
                if col < 2:
                    out.write('\\hfill\n')
            out.write('\\par\\smallskip\n')

        if remainder == 1:
            out.write('\\noindent\n')
            _write_peak_summary_minipage(out, peaks[cursor], columns, joined_by_peak)
            out.write('\\hfill\n')
            out.write('\\begin{minipage}[t]{0.65\\linewidth}\\centering\n')
            out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(filename) + '}\n')
            out.write('\\end{minipage}\\par\\smallskip\n')
            plot_written = True
        elif remainder:
            out.write('\\noindent\n')
            for col in range(remainder):
                _write_peak_summary_minipage(out, peaks[cursor], columns, joined_by_peak)
                cursor += 1
                if col + 1 < remainder:
                    out.write('\\hfill\n')
            out.write('\\par\\smallskip\n')

        if not plot_written:
            out.write('\\begin{center}\n')
            out.write('\\includegraphics[width=0.82\\linewidth]{./pdf/' + _latex_escape(filename) + '}\n')
            out.write('\\end{center}\\smallskip\n')
        slice_figures = unit.get('slice_figures', []) if pseudo.get('kind') != '2d' else []
        if slice_figures:
            out.write('\\noindent\\textit{Fit across pseudo-axis slices}\\par\\smallskip\n')
            out.write('\\begin{center}\n')
            for i, slice_file in enumerate(slice_figures):
                out.write('\\begin{minipage}[t]{0.19\\linewidth}\\centering\n')
                out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(slice_file) + '}\n')
                out.write('\\end{minipage}')
                if i + 1 < len(slice_figures): out.write('\\hfill\n')
                if (i + 1) % 5 == 0 and i + 1 < len(slice_figures): out.write('\\par\\smallskip\n')
            out.write('\n\\end{center}\n')

        fits = [(peak, peak_figures.get(peak)) for peak in unit['peaks'] if peak_figures.get(peak)]
        if fits:
            out.write('\\begin{center}\n')
            for i, (peak, decay_file) in enumerate(fits):
                out.write('\\begin{minipage}[t]{0.19\\linewidth}\\centering\n')
                out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(decay_file) + '}\\\\[-1mm]\n')
                out.write('\\scriptsize ' + _latex_escape(peak) + '\n\\end{minipage}')
                if i + 1 < len(fits): out.write('\\hfill\n')
                if (i + 1) % 5 == 0 and i + 1 < len(fits): out.write('\\par\\smallskip\n')
            out.write('\n\\end{center}\n')

    if joined_rows:
        if (pseudo.get('analysis') or {}).get('name') == 'CPMG':
            out.write('\\section*{Peak fitting results}\n')
            fit_columns=list(pseudo.get('columns') or [])
            fit_rows=[list(r) for r in pseudo.get('rows') or []]
            _write_joined_results_table(out, fit_columns, fit_rows, len(fit_columns))
        else:
            out.write('\\section*{Fitted parameters and metrics}\n')
            _write_joined_results_table(out, columns, joined_rows, split_at)


def _peak_tables(frame):
    """Return report-ready reference/full peak lists for multidimensional data.

    Reference peaks remain a physical-3D audit feature.  The full peak list is
    dimension-independent so physical 2D reports can append fitting results.
    """
    try:
        dim = int(getattr(frame, 'dim', 0) or 0)
        if dim < 2:
            return None
    except Exception:
        return None
    refs = list(frame.get_reference_peaks() or []) if dim == 3 else []
    ref_headers = list(frame.get_reference_peak_headers() or [])
    ref_rows = []
    for i, pk in enumerate(refs):
        match = re.findall(r'[0-9]+', str(pk.name))
        ref_rows.append([match[0] if match else str(i + 1), str(pk.name), str(pk.x), str(pk.y)])

    payload = frame.get_full_peak_payload() or {}
    raw_rows = payload.get('rows') or payload.get('peaks') or []
    width = max([len(row) for row in raw_rows], default=0)
    full_headers = list(frame.get_full_peak_headers(row_width=width or 5) or [])
    # Match the Full Peak List Show window's 3D generated-name presentation.
    def split_name(value):
        left, sep, right = str(value).rpartition('_')
        digits = re.findall(r'\d+', right)
        return (left, max(enumerate(digits), key=lambda item: (len(item[1]), -item[0]))[1]) if sep and left and digits else None
    split_names = dim == 3 and bool(raw_rows) and all(row and split_name(row[0]) is not None for row in raw_rows)
    if split_names and full_headers:
        full_headers = ['nResID', 'Number'] + full_headers[1:]
    intensity_col = next((i for i, h in enumerate(full_headers) if 'intensity' in str(h).lower()), None)
    if intensity_col is not None:
        full_headers[intensity_col] = 'SNR'
    ppm_cols = [i for i, h in enumerate(full_headers) if 'ppm' in str(h).lower()]
    sigma = frame.get_noise_sigma() if hasattr(frame, 'get_noise_sigma') else None
    full_rows = []
    import math
    for row in raw_rows:
        raw_values = ([split_name(row[0])[0], split_name(row[0])[1]] + list(row[1:])) if split_names else list(row)
        values = [str(v) for v in raw_values]
        for ppm_col in ppm_cols:
            if ppm_col < len(raw_values):
                try: values[ppm_col] = '%.3f' % (math.trunc(float(raw_values[ppm_col]) * 1000.0) / 1000.0)
                except Exception: pass
        if intensity_col is not None and intensity_col < len(values):
            if sigma is None:
                values[intensity_col] = 'N/A'
            else:
                try: values[intensity_col] = '%.4g' % (float(raw_values[intensity_col]) / sigma)
                except Exception: pass
        full_rows.append(values)
    return ref_headers, ref_rows, full_headers, full_rows


def _save_three_d_reference_slices(frame, report_dir, warnings):
    """Export one compact Slice2D review figure for every physical-3D reference peak."""
    try:
        if int(getattr(frame, 'dim', 0) or 0) != 3:
            return []
        # The 1D lower pane obtains calculated traces from pkSlice1Ddec and
        # peak markers from the authoritative Full list.  Both are populated
        # only after the complete Review dependency chain has been loaded.
        ensure_review = getattr(frame, 'ensure_workflow_review_inputs_loaded', None)
        if callable(ensure_review):
            ok, message = ensure_review()
            if not ok:
                raise RuntimeError(message or '3D review inputs could not be loaded')
        notebook = frame.parent
        # Report export must not navigate through open_workflow('slices').  That
        # public workflow entry point may reload the main spectrum as part of
        # prepare_workflow(), which can close transient Fit Peaks windows while
        # their background picker still has wx.CallAfter callbacks queued.  The
        # report only needs the already-loaded Slice2D panel, so create/select
        # that page directly without changing project/workflow state.
        if not notebook.PageExists('2D Slices'):
            notebook.AddTabFour(True, frame)
        view = notebook.get_page_by_title('2D Slices') or getattr(notebook, 'tabFour', None)
        if view is None:
            raise RuntimeError('The 2D Slices view could not be created.')
        refs = list(frame.get_reference_peaks() or [])
        old = (view.ComboBox1.GetSelection(), view.ComboBox2.GetSelection(),
               bool(view.cb_decon.GetValue()), bool(view.cb_1d.GetValue()),
               bool(view.cb_grid_auto.GetValue()))
        files = []
        try:
            # Use the same toggle paths as the Slice2D toolbar so the
            # rendered state is identical to an interactive Decon + 1D + Show
            # peaks view.
            view._toolbar_decon(True)
            view._toolbar_1d(True)
            view._toolbar_peaks(True)
            for i, pk in enumerate(refs):
                view.ComboBox1.SetSelection(i)
                # Right is irrelevant in this report view; keeping it on Left
                # also prevents it from selecting a different reference plane.
                if i < view.ComboBox2.GetCount(): view.ComboBox2.SetSelection(i)
                view.ax_reset1 = view.ax_reset2 = 1
                view.draw_figure()
                # Projection guides are animated for interactive blitting.  For
                # the static report make only the green (Left/current) guides
                # ordinary artists and suppress all blue (Right) guides.
                artists = list(getattr(view, '_projection_artists', []) or [])
                for j, artist in enumerate(artists):
                    artist.set_animated(False)
                    artist.set_visible(j in (0, 2, 4))
                view.canvas.draw()
                filename = 'slice2d_reference_%04d.pdf' % (i + 1)
                view.canvas.print_figure(str(report_dir / filename))
                files.append((str(pk.name), filename))
                for artist in artists:
                    artist.set_visible(True); artist.set_animated(True)
        finally:
            a, b, dec, one_d, peaks = old
            if a >= 0: view.ComboBox1.SetSelection(a)
            if b >= 0: view.ComboBox2.SetSelection(b)
            view._toolbar_decon(dec)
            view._toolbar_1d(one_d)
            view._toolbar_peaks(peaks)
            view.ax_reset1 = view.ax_reset2 = 1
            view.draw_figure()
        return files
    except Exception as exc:
        warnings.append('3D reference Slice2D figures unavailable: %s' % exc)
        return []

def _write_report(frame, tex_path, report_dir, warnings, process_info=None, workflow=None, pseudo=None, reference_slices=None, peak_shape_info=None):
    project, processing, decon, peak_shape_params = _collect_rows(frame)
    workflow = workflow or {'rows': [], 'pseudo': False}
    pseudo = pseudo or {}
    reference_slices = reference_slices or []
    process_info = process_info or {'rows': [], 'processing_rows': [], 'spectrum_rows': [], 'debug_rows': [], 'acquisition_columns': [], 'acquisition_rows': [], 'pseudo_axis_columns': [], 'pseudo_axis_rows': [], 'pseudo_axis_explanation': ''}
    peak_path = frame._project_path(_ctrl_value(frame, 'referencePeakBox'))
    full_path = frame._project_path(_ctrl_value(frame, 'fullPeakBox'))
    decon_path = (getattr(frame, 'spectrumfile', '') or frame._project_path(_ctrl_value(frame, 'infileBox'))) + '.decon'
    # Available-result status is part of the project summary rather than a
    # separate late report section.
    is_pseudo2d = pseudo.get('kind') == 'pseudo2d'
    # Match the main-window deconvolution lamp: it is lit exactly when the
    # active deconvolution product exists, rather than when a guessed suffix does.
    try:
        spectrum_for_lamp = frame._project_path(_ctrl_value(frame, 'infileBox'))
        decon_lamp_path = frame._active_deconvolution_path(spectrum_for_lamp)
    except Exception:
        decon_lamp_path = decon_path
    availability = []
    if not is_pseudo2d:
        availability.append(('Reference peaks', _count_rows(peak_path) if peak_path and os.path.isfile(peak_path) else 'not available'))
    availability.extend([
        ('Full peaks', _count_rows(full_path) if full_path and os.path.isfile(full_path) else 'not available'),
        ('Deconvolved?', 'available' if decon_lamp_path and os.path.isfile(decon_lamp_path) else 'not available'),
    ])
    project.extend(availability)
    if is_pseudo2d:
        project[:] = [row for row in project if row[0] not in ('Reference peaklist', 'Reference peaks')]
    with open(tex_path, 'w', encoding='utf-8') as out:
        out.write(r"""\documentclass[a4paper,11pt]{article}
\usepackage[a4paper,margin=18mm]{geometry}
\usepackage{graphicx,booktabs,array,longtable,fancyhdr,xcolor,amssymb,amsmath}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\pagestyle{fancy}
\fancyhf{}
\lhead{UniDec NMR Project Summary}
\rhead{\thepage}
\setlength{\parindent}{0pt}
\begin{document}
""")
        _write_workflow_table(out, workflow)
        _write_analysis_summary(out, pseudo)
        out.write('\\section*{Project Summary}\n')
        working_row = next((row for row in project if row[0] == 'Working dir'), None)
        project = [row for row in project if row[0] != 'Working dir']
        if working_row:
            out.write('\\noindent\\scriptsize\\textbf{Working dir:} \\texttt{' + _latex_escape(working_row[1]) + '}\\normalsize\\par\\smallskip\n')

        # The compact 1D spectrum view belongs beside the project metadata.
        # projection.pdf is the Projection-window figure (2D projections) and
        # is deliberately shown later in the Spectrum section.
        spectrum = report_dir / 'spectrum.pdf'
        out.write('\\begin{minipage}[t]{0.46\\linewidth}\n\\vspace{0pt}\n')
        _write_kv_table(out, '', project)
        out.write('\\end{minipage}\\hfill\\begin{minipage}[t]{0.52\\linewidth}\n\\vspace{0pt}\n')
        if spectrum.exists():
            out.write('\\includegraphics[width=\\linewidth]{./pdf/spectrum.pdf}\n')
        else:
            out.write('Spectrum projection unavailable.\n')
        out.write('\\end{minipage}\n\n')

        # Acquisition metadata is presented as one compact two-column block:
        # the general acquisition values on the left and the per-dimension
        # acquisition table on the right.  Both are part of the Acquisition
        # section, so there is deliberately no separate "Acquisition dimensions"
        # subsection/heading.
        if process_info.get('rows') or process_info.get('acquisition_rows') or process_info.get('pseudo_axis_rows') or process_info.get('debug_rows'):
            out.write('\\section*{Acquisition}\n')
            out.write('\\begin{minipage}[t]{0.48\\linewidth}\n\\vspace{0pt}\n')
            if process_info.get('rows'):
                _write_kv_table(out, '', process_info['rows'])
            out.write('\\end{minipage}\\hfill\\begin{minipage}[t]{0.50\\linewidth}\n\\vspace{0pt}\n')
            if process_info.get('acquisition_rows'):
                _write_matrix_table(out, '', process_info.get('acquisition_columns', []), process_info['acquisition_rows'])
            out.write('\\end{minipage}\n\n')
            if process_info.get('fid_path'):
                out.write('\\noindent\\scriptsize\\textbf{FID path:} \\texttt{' + _latex_escape(process_info['fid_path']) + '}\\normalsize\\par\\smallskip\n')
            if process_info.get('combine_rows'):
                out.write('\\noindent\\textbf{Combined raw-data sources}\\par\\smallskip\n')
                _write_matrix_table(out, '', ['Experiment', 'Source (relative to working dir)', 'Rows'], process_info['combine_rows'])
            if process_info.get('pseudo_axis_rows'):
                out.write('\\subsection*{Pseudo-axis sampling}\n')
                explanation = process_info.get('pseudo_axis_explanation', '')
                if explanation:
                    out.write(_latex_escape(explanation) + '\n\n')
                _write_matrix_table(out, '', process_info.get('pseudo_axis_columns', []),
                                    process_info.get('pseudo_axis_rows', []))
            # NUS/XYZA diagnostics are intentionally retained in process_info for
            # debugging, but are not rendered in the normal project report.

        proc_columns, proc_rows = processing
        if proc_rows or process_info.get('processing_rows') or (report_dir / 'process_phased.pdf').exists() or (report_dir / 'process_fid.pdf').exists():
            out.write('\\section*{Processing}\n')
            if process_info.get('processing_rows'):
                _write_kv_table(out, '', process_info['processing_rows'])
            out.write('\\begin{minipage}[t]{0.48\\linewidth}\n\\vspace{0pt}\n')
            if proc_rows:
                _write_matrix_table(out, '', proc_columns, proc_rows)
            out.write('\\end{minipage}\\hfill')
            out.write('\\begin{minipage}[t]{0.24\\linewidth}\n\\vspace{0pt}\n')
            if (report_dir / 'process_phased.pdf').exists():
                out.write('\\includegraphics[width=\\linewidth]{./pdf/process_phased.pdf}\n')
            out.write('\\end{minipage}\\hfill')
            out.write('\\begin{minipage}[t]{0.24\\linewidth}\n\\vspace{0pt}\n')
            if (report_dir / 'process_fid.pdf').exists():
                out.write('\\includegraphics[width=\\linewidth]{./pdf/process_fid.pdf}\n')
            out.write('\\end{minipage}\n\n')

        # The Process -> Projections window is the interactive phasing view.
        # Show its current matplotlib canvas immediately before the scripts.
        if pseudo.get('kind') != 'pseudo2d' and (report_dir / 'phasing.pdf').exists():
            out.write('\\subsection*{Indirect Phasing}\n')
            out.write('\\begin{center}\n')
            out.write('\\includegraphics[width=0.9\\linewidth]{./pdf/phasing.pdf}\n')
            out.write('\\end{center}\n\n')

        # Record the exact current NMRPipe scripts immediately before the
        # spectrum-dimension summary.  Verbatim mode is required because shell
        # syntax contains many characters that have special meaning in LaTeX.
        for script_title, script_name, script_text in _read_processing_scripts(frame):
            _write_script_verbatim(out, script_title, script_name, script_text)

        if process_info.get('spectrum_rows'):
            out.write('\\section*{Spectrum dimensions}\n')
            _write_matrix_table(out, '', ['Minimum ppm', 'Maximum ppm', 'Points'], process_info['spectrum_rows'])

        # Noise comes before UniDec.  Put the concise statistics immediately
        # into the left-hand minipage (without a redundant Statistics heading)
        # and keep the live Noise-window plot on the right.  Threshold belongs
        # here because it is interpreted relative to the measured noise.
        out.write('\\section*{Noise}\n')
        noise_rows = [('Threshold', _ctrl_value(frame, 'threshBox'))] + _noise_summary_rows(frame)
        out.write('\\begin{minipage}[t]{0.46\\linewidth}\n\\vspace{0pt}\n')
        _write_kv_table(out, '', noise_rows)
        out.write('\\end{minipage}\\hfill\\begin{minipage}[t]{0.52\\linewidth}\n\\vspace{0pt}\n')
        if (report_dir/'noise.pdf').exists():
            out.write('\\includegraphics[width=\\linewidth]{./pdf/noise.pdf}\n')
        else:
            out.write('Noise plot unavailable.\n')
        out.write('\\end{minipage}\n\n')

        _write_kv_table(out, 'UniDec', decon)

        peak_tables = _peak_tables(frame)

        # Peak-shape parameters and the Fit Peaks matplotlib figure are one
        # logical report item, so keep them together in facing minipages.
        shape_columns, shape_rows = peak_shape_params
        if pseudo.get('kind') != 'pseudo2d' and (shape_rows or (report_dir / 'shape.pdf').exists()):
            out.write('\\section*{Peak shape parameters}\n')
            out.write('\\begin{minipage}[t]{0.46\\linewidth}\n\\vspace{0pt}\n')
            if peak_shape_info:
                _write_kv_table(out, '', peak_shape_info)
            if shape_rows:
                _write_matrix_table(out, '', shape_columns, shape_rows)
            else:
                out.write('Peak shape parameters unavailable.\n')
            out.write('\\end{minipage}\\hfill\\begin{minipage}[t]{0.52\\linewidth}\n\\vspace{0pt}\n')
            if (report_dir / 'shape.pdf').exists():
                out.write('\\includegraphics[width=\\linewidth]{./pdf/shape.pdf}\n')
                if (report_dir / 'shape_widths.pdf').exists():
                    out.write('\\includegraphics[width=\\linewidth]{./pdf/shape_widths.pdf}\n')
            else:
                out.write('Peak shape fit unavailable.\n')
            out.write('\\end{minipage}\n\n')

        # Bore plane precedes the Projection-window spectrum section.
        if pseudo.get('kind') != 'pseudo2d' and (report_dir/'peak.pdf').exists():
            out.write('\\section*{Bore plane}\n')
            if (report_dir/'peak_ornament.pdf').exists():
                out.write('\\begin{minipage}[t]{0.49\\linewidth}\\includegraphics[width=\\linewidth]{./pdf/peak.pdf}\\end{minipage}\\hfill\n')
                out.write('\\begin{minipage}[t]{0.49\\linewidth}\\includegraphics[width=\\linewidth]{./pdf/peak_ornament.pdf}\\end{minipage}\n')
            else:
                out.write('\\begin{center}\\includegraphics[width=0.82\\linewidth]{./pdf/peak.pdf}\\end{center}\n')

        if (report_dir/'projection.pdf').exists():
            out.write('\\section*{Spectrum}\n')
            out.write('\\begin{center}\\includegraphics[width=0.9\\linewidth]{./pdf/projection.pdf}\\end{center}\n')
            if (report_dir/'projection_decon.pdf').exists():
                out.write('\\begin{center}\\includegraphics[width=0.9\\linewidth]{./pdf/projection_decon.pdf}\\end{center}\n')
        _write_pseudo_report(out, pseudo)

        if warnings:
            out.write('\\section*{Report notes}\n\\begin{itemize}\n')
            for warning in warnings:
                out.write('\\item ' + _latex_escape(warning) + '\n')
            out.write('\\end{itemize}\n')

        # Physical-3D review material is intentionally last: it is detailed
        # audit material rather than part of the compact project narrative.
        if peak_tables is not None:
            ref_headers, ref_rows, full_headers, full_rows = peak_tables
            if ref_rows:
                out.write('\\section*{Reference peak list}\n')
                _write_two_column_long_table(out, ref_headers, ref_rows)
            if full_rows:
                out.write('\\section*{Full peak list}\n')
                full_headers, full_rows, fit_split = _append_fitting_to_full_peaks(full_headers, full_rows, pseudo)
                if fit_split is not None:
                    _write_long_table_with_split(out, full_headers, full_rows, fit_split, font='\\tiny', compact=True)
                else:
                    _write_long_table(out, full_headers, full_rows, font='\\tiny', compact=True)
            if reference_slices:
                out.write('\\section*{3D reference-peak review}\n')
                out.write('\\footnotesize Each figure reproduces the Slice2D review view with Decon, Peaks and 1D enabled. Only the current Left reference indicator is shown on the projections.\\normalsize\n\n')
                for peak_name, filename in reference_slices:
                    out.write('\\begin{minipage}[t]{0.49\\linewidth}\\centering\n')
                    out.write('\\includegraphics[width=\\linewidth]{./pdf/' + _latex_escape(filename) + '}\\\\[-1mm]\n')
                    out.write('\\scriptsize ' + _latex_escape(peak_name) + '\n\\end{minipage}\\hfill\n')
                out.write('\\par\\smallskip\n')
        out.write('\\end{document}\n')


def project_summary_stages(frame):
    """Return the coarse report-build stages for the current dataset topology.

    The list intentionally follows real report work rather than a synthetic
    percentage.  It is consumed by the standard magnet/spin progress control.
    """
    spectral = int(getattr(getattr(frame, 'state', None), 'dimension', 0) or getattr(frame, 'dim', 1) or 1)
    pseudo_axis = bool(getattr(getattr(frame, 'state', None), 'pseudo_axis', False))
    stages = [
        'Saving spectrum and noise views',
        'Saving projection views',
        'Collecting acquisition and processing details',
        'Collecting workflow status',
    ]
    if pseudo_axis:
        if spectral == 1:
            stages += ['Collecting pseudo2D fitting results', 'Saving pseudo2D group and slice-fit figures']
        else:
            stages += ['Collecting pseudo3D fitting results', 'Saving pseudo3D reference and fitting figures',
                       'Collecting pseudo3D analysis results']
    else:
        if spectral == 2 and bool(getattr(getattr(frame, 'store', None), 'analysis', {}).get('fitting_results_ready')):
            stages += ['Collecting 2D fitting results', 'Saving 2D group fitting figures']
        # Reference-peak material is part of the ordinary spectral reports,
        # irrespective of whether the spectrum is 1D, 2D, 3D or 4D.
        stages += ['Saving reference-peak views']
        if spectral >= 2:
            stages += ['Saving peak-shape and dimensional views']
    stages += ['Writing LaTeX report', 'Compiling PDF', 'Finalising summary']
    return stages

def generate_project_summary(frame, output_pdf=None, progress_callback=None):
    """Generate a PDF summary from the current deconFrame state.

    Returns ``(pdf_path, warnings)``.  Optional sections fail softly and are
    recorded in the report notes; LaTeX/compiler failures are fatal.
    """
    state = getattr(frame, 'state', None)
    base = Path(getattr(state, 'working_dir', '') or _ctrl_value(frame, 'dirBox') or os.getcwd()).expanduser()
    if not base.is_absolute():
        base = (Path.cwd() / base).resolve()
    base.mkdir(parents=True, exist_ok=True)
    report_dir = base / 'pdf'
    report_dir.mkdir(parents=True, exist_ok=True)
    warnings = []

    def progress(label):
        if callable(progress_callback):
            progress_callback(label)

    # Every live report figure depends on the configured spectrum being
    # materialised.  A freshly loaded system/project file can have valid paths
    # while the GUI canvases still contain their blank startup state.
    ensure_spectrum = getattr(frame, 'ensure_workflow_spectrum_loaded', None)
    if callable(ensure_spectrum):
        try:
            if not ensure_spectrum():
                warnings.append('Configured spectrum could not be loaded before report figures were created.')
        except Exception as exc:
            warnings.append('Spectrum synchronisation before report failed: %s' % exc)
    try:
        plot_noise = getattr(frame, 'plot_noise_histogram', None)
        if callable(plot_noise) and getattr(frame, 'data', None) is not None:
            plot_noise(frame.data)
    except Exception as exc:
        warnings.append('Noise plot refresh unavailable: %s' % exc)

    # Report figures must be generated from materialised project data, not
    # merely from filenames restored by the system file.  This is deliberately
    # done before any figure/window is constructed so projection, noise, bore
    # and peak-shape views all see the same loaded spectrum.
    ensure_spectrum = getattr(frame, 'ensure_workflow_spectrum_loaded', None)
    if callable(ensure_spectrum):
        try:
            if not ensure_spectrum():
                warnings.append('Report spectrum could not be loaded before figure export.')
        except Exception as exc:
            warnings.append('Report spectrum load failed: %s' % exc)
    try:
        topology = frame._active_topology()
        ref_value = _ctrl_value(frame, 'referencePeakBox').strip()
        if ref_value and int(getattr(topology, 'spectral_dim_count', 1) or 1) >= 2:
            ensure_ref = getattr(frame, 'ensure_workflow_reference_stage_loaded', None)
            if callable(ensure_ref) and not ensure_ref():
                warnings.append('Reference peak list could not be loaded before figure export.')
    except Exception as exc:
        warnings.append('Reference peak-list preflight failed: %s' % exc)

    # Export the exact live matplotlib figures used by the NMR GUI.
    progress('Saving spectrum and noise views')
    try:
        frame.canvas.print_figure(str(report_dir / 'spectrum.pdf'))
    except Exception as exc:
        warnings.append('Spectrum figure unavailable: %s' % exc)
    try:
        frame.noiseFig.savefig(str(report_dir / 'noise.pdf'))
    except Exception as exc:
        warnings.append('Noise figure unavailable: %s' % exc)

    progress('Saving projection views')
    _save_projection(frame, report_dir, warnings)
    progress('Collecting acquisition and processing details')
    process_info = _process_report_data(frame, report_dir, warnings)
    progress('Collecting workflow status')
    workflow = _workflow_report_data(frame, warnings)

    # Pseudo-dimensional export is deliberately reported as separate stages:
    # these paths silently open fitting/analysis windows and can dominate the
    # report build time.  The exporter itself remains the single authority for
    # the generated data.
    try:
        from spinDecon.domain.analysis_mode import AnalysisMode
        _summary_mode = AnalysisMode.from_project_state(frame.state)
    except Exception:
        _summary_mode = None
    if _summary_mode is not None and _summary_mode.has_pseudo_axis:
        if _summary_mode.spectral_dimensions == 1:
            progress('Collecting pseudo2D fitting results')
        else:
            progress('Collecting pseudo3D fitting results')
    elif (_summary_mode is not None and int(_summary_mode.spectral_dimensions) == 2
          and bool(getattr(getattr(frame, 'store', None), 'analysis', {}).get('fitting_results_ready'))):
        progress('Collecting 2D fitting results')
    pseudo = _pseudo_report_data(frame, report_dir, warnings, workflow, progress_callback=progress_callback)
    if pseudo.get('kind') != 'pseudo2d' and not getattr(_summary_mode, 'has_pseudo_axis', False):
        progress('Saving reference-peak views')
    _has_reference_peaks = bool(getattr(frame, 'PEAK', 0) or (getattr(frame, 'store', None) and getattr(frame.store, 'peak_lists', {}).get('reference')))
    if pseudo.get('kind') != 'pseudo2d' and _has_reference_peaks:
        _save_peak_views(frame, report_dir, warnings)
    if not getattr(_summary_mode, 'has_pseudo_axis', False) and int(getattr(_summary_mode, 'spectral_dimensions', 1) or 1) >= 2:
        progress('Saving peak-shape and dimensional views')
    peak_shape_info = []
    if pseudo.get('kind') != 'pseudo2d' and _summary_mode is not None and int(getattr(_summary_mode, 'spectral_dimensions', 1) or 1) >= 1:
        peak_shape_info = _save_peak_shape(frame, report_dir, warnings)
    reference_slices = []
    if _has_reference_peaks and _summary_mode is not None and not _summary_mode.has_pseudo_axis and int(_summary_mode.spectral_dimensions or 0) == 3:
        # The 3D report consumes the same four products as Review picked peaks.
        # Materialise them in the established dependency order before tables or
        # Slice2D figures inspect DataStore.  This is a no-op when already loaded.
        ensure_review = getattr(frame, 'ensure_workflow_review_inputs_loaded', None)
        if callable(ensure_review):
            try:
                ok, message = ensure_review()
                if not ok:
                    warnings.append('3D report review inputs incomplete: %s' % message)
            except Exception as exc:
                warnings.append('3D report could not synchronise review inputs: %s' % exc)
        progress('Saving 3D reference-peak review figures')
        reference_slices = _save_three_d_reference_slices(frame, report_dir, warnings)

    progress('Writing LaTeX report')
    tex_path = report_dir / 'summary.tex'
    _write_report(frame, tex_path, report_dir, warnings, process_info=process_info, workflow=workflow, pseudo=pseudo, reference_slices=reference_slices, peak_shape_info=peak_shape_info)

    progress('Compiling PDF')
    pdflatex = shutil.which('pdflatex')
    if not pdflatex:
        raise RuntimeError('pdflatex was not found on PATH. summary.tex was created at %s' % tex_path)
    # Compile from the project working directory.  Figure references in the
    # generated LaTeX are therefore explicitly rooted at ./pdf/.  Keep all
    # LaTeX build artefacts in the same pdf directory.
    tex_arg = str(Path('pdf') / tex_path.name)
    proc = subprocess.run(
        [pdflatex, '-interaction=nonstopmode', '-halt-on-error',
         '-output-directory=pdf', tex_arg],
        cwd=str(base), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding='utf-8', errors='replace',
    )
    if proc.returncode != 0:
        log_path = report_dir / 'pdflatex.log'
        log_path.write_text(proc.stdout, encoding='utf-8', errors='replace')
        raise RuntimeError('pdflatex failed; see %s' % log_path)

    progress('Finalising summary')
    built_pdf = report_dir / 'summary.pdf'
    final_pdf = Path(output_pdf).expanduser() if output_pdf else base / 'summary.pdf'
    if not final_pdf.is_absolute():
        final_pdf = (base / final_pdf).resolve()
    if built_pdf.resolve() != final_pdf.resolve():
        shutil.copy2(built_pdf, final_pdf)

    # Report generation opens several lazy workflow pages internally.  Return
    # the main application to the NMR tab so generating a summary never leaves
    # the user on a report-helper/viewer tab.
    try:
        notebook = getattr(frame, 'parent', None)
        if notebook is not None and hasattr(notebook, 'select_page'):
            notebook.select_page('NMR')
        elif notebook is not None and hasattr(notebook, 'SetSelection'):
            # NMR is the first notebook page now that the obsolete Start tab
            # has been removed.
            notebook.SetSelection(0)
    except Exception:
        pass
    return str(final_pdf), warnings
