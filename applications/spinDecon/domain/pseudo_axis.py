"""Shared access to conversion-generated pseudo_axis.tsv tables."""
from __future__ import annotations

import csv
import os


class PseudoAxisError(ValueError):
    pass


def _project_owner(obj):
    """Return the object exposing project state/decon parameter helpers."""
    seen = set()
    while obj is not None and id(obj) not in seen:
        seen.add(id(obj))
        tab_one = getattr(obj, 'tabOne', None)
        if tab_one is not None:
            return tab_one
        if getattr(obj, 'state', None) is not None:
            return obj
        obj = getattr(obj, 'parent', None)
    return None


def pseudo_axis_path(obj, fallback_spec_dir=None):
    owner = _project_owner(obj)
    state = getattr(owner, 'state', None) if owner is not None else None
    if state is not None and hasattr(state, 'spec_dir'):
        try:
            return os.path.join(os.path.normpath(state.spec_dir()), 'pseudo_axis.tsv')
        except Exception:
            pass
    if fallback_spec_dir:
        return os.path.join(os.path.normpath(fallback_spec_dir), 'pseudo_axis.tsv')
    return os.path.join('spec', 'pseudo_axis.tsv')


class PseudoAxisTable:
    def __init__(self, path, headers, rows):
        self.path = path
        self.headers = headers
        self.rows = rows

    @classmethod
    def load(cls, path):
        if not os.path.exists(path):
            raise PseudoAxisError('Pseudo-axis table not found: %s' % path)
        with open(path, 'r', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            headers = list(reader.fieldnames or [])
            if not headers:
                raise PseudoAxisError('Pseudo-axis table has no header')
            rows = [row for row in reader if any(str(v or '').strip() for v in row.values())]
        if not rows:
            raise PseudoAxisError('Pseudo-axis table contains no data rows')
        return cls(path, headers, rows)

    @property
    def data_columns(self):
        return [h for h in self.headers if h.strip().lower() != 'spectrum']

    def numeric_values(self, column):
        if column not in self.headers:
            raise PseudoAxisError("Pseudo-axis column '%s' was not found" % column)
        values = []
        for row in self.rows:
            text = str(row.get(column, '')).strip()
            try:
                values.append(float(text))
            except (TypeError, ValueError):
                raise PseudoAxisError("Pseudo-axis column '%s' contains a non-numeric value: %s" % (column, text))
        return values

    def varying_columns(self):
        varying = []
        for column in self.data_columns:
            try:
                vals = self.numeric_values(column)
            except PseudoAxisError:
                continue
            if len(set(vals)) > 1:
                varying.append(column)
        return varying

    def default_column(self, saved_name=''):
        if saved_name in self.data_columns:
            return saved_name
        columns = self.data_columns
        if len(columns) == 1:
            return columns[0]
        varying = self.varying_columns()
        if len(varying) == 1:
            return varying[0]
        return columns[0] if columns else ''


def load_saved_column(obj, key='pseudoAxisAnalysisColumn'):
    owner = _project_owner(obj)
    if owner is None:
        return ''
    parser = getattr(owner, 'Parse', None)
    parfile = getattr(owner, 'deconParFile', None)
    if parser is None or not parfile:
        return ''
    try:
        value = parser(parfile, key, default='')
    except TypeError:
        try:
            value = parser(parfile, key)
        except Exception:
            return ''
    except Exception:
        return ''
    return str(value or '').strip()


def save_selected_column(obj, column, key='pseudoAxisAnalysisColumn'):
    owner = _project_owner(obj)
    if owner is None or not column:
        return False
    parfile = getattr(owner, 'deconParFile', None)
    if not parfile:
        return False
    working_dir = '.'
    dir_box = getattr(owner, 'dirBox', None)
    if dir_box is not None:
        try:
            working_dir = str(dir_box.GetValue() or '.').strip() or '.'
        except Exception:
            pass
    savefile = parfile if os.path.isabs(parfile) else os.path.join(working_dir, parfile)
    from spinDecon.project.parameter_store import update_parameter_file
    update_parameter_file(savefile, {key: column}, source_path=parfile)
    return True
