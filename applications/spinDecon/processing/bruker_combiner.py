"""Combine numbered 1D Bruker experiments into a decon pseudo-2D dataset.

The destination is the directory containing the numbered experiments.  The
source experiments are never modified.  The generated files are intentionally
only the Bruker subset required by decon/bruk2pipe, not a TopSpin experiment.
"""
from __future__ import annotations

import json
import itertools
import math
import os
import re
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


STRUCTURAL_KEYS = (
    'TD', 'DTYPA', 'BYTORDA', 'AQ_mod', 'SW_h', 'SFO1', 'BF1', 'O1', 'NUC1',
    'DECIM', 'DSPFVS', 'GRPDLY',
)
# Parameters which are bookkeeping rather than useful experimental coordinates.
IGNORE_VARYING = {
    'DATE', 'DATE_START', 'OWNER', 'ORIGIN', 'TITLE', 'NS', 'DS', 'RG',
    'AQSEQ', 'TD', 'DTYPA', 'BYTORDA', 'AQ_mod', 'SW_h', 'SFO1', 'BF1',
    'O1', 'NUC1', 'DECIM', 'DSPFVS', 'GRPDLY',
}


@dataclass
class BrukerExperimentInfo:
    number: int
    path: Path
    raw_file: Path
    raw_kind: str
    raw_bytes: int
    record_bytes: int
    row_count: int
    parameters: dict[str, Any] = field(default_factory=dict)
    pseudo_columns: list[str] = field(default_factory=list)
    pseudo_rows: list[list[str]] = field(default_factory=list)


@dataclass
class CombinationInspection:
    experiments: list[BrukerExperimentInfo]
    varying_parameters: dict[str, list[str]]
    numeric_varying_parameters: dict[str, list[str]]
    errors: list[str]
    warnings: list[str]


def _normalise_value(value: str) -> str:
    value = str(value).strip()
    if value.startswith('<') and value.endswith('>'):
        value = value[1:-1].strip()
    return value


def read_jcamp(path: Path) -> dict[str, Any]:
    """Read Bruker JCAMP scalars and arrays sufficiently for comparison."""
    result: dict[str, Any] = {}
    if not path.is_file():
        return result
    lines = path.read_text(errors='replace').splitlines()
    i = 0
    header = re.compile(r'^##\$([^=]+)=\s*(.*)$')
    while i < len(lines):
        m = header.match(lines[i])
        if not m:
            i += 1
            continue
        key, rest = m.group(1).strip(), m.group(2).strip()
        if rest.startswith('('):
            values = []
            i += 1
            while i < len(lines) and not lines[i].startswith('##'):
                values.extend(lines[i].strip().split())
                i += 1
            result[key] = tuple(_normalise_value(v) for v in values)
            continue
        result[key] = _normalise_value(rest)
        i += 1
    # Bruker stores common pulse-program coordinates as indexed arrays.
    # Expose their elements using the same familiar names used by TopSpin/decon
    # (D20 == D[20], P1 == P[1], CNST3 == CNST[3]) while retaining the
    # original tuple for provenance.
    for array_name in ('D', 'P', 'CNST'):
        values = result.get(array_name)
        if isinstance(values, tuple):
            for index, value in enumerate(values):
                result[f'{array_name}{index}'] = value
    return result


def _scalar_text(value: Any) -> str | None:
    if isinstance(value, tuple):
        return None
    if value is None:
        return None
    return str(value).strip()


def _as_float(value: Any) -> float | None:
    value = _scalar_text(value)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _same_parameter(a: Any, b: Any, key: str) -> bool:
    if a == b:
        return True
    # Numeric acquisition values may differ only in formatting.
    fa, fb = _as_float(a), _as_float(b)
    if fa is not None and fb is not None:
        tol = 1e-7 * max(1.0, abs(fa), abs(fb))
        return abs(fa - fb) <= tol
    return False


def bruker_record_bytes(parameters: dict[str, Any]) -> int:
    td = int(float(parameters.get('TD', 0) or 0))
    dtypa = int(float(parameters.get('DTYPA', 0) or 0))
    if td <= 0:
        raise ValueError('TD is missing or invalid')
    word_bytes = 8 if dtypa == 2 else 4
    return int(math.ceil(float(td * word_bytes) / 1024.0) * 1024)


def _source_pseudo_axis(path: Path, params: dict[str, Any]) -> tuple[list[str], list[list[str]]]:
    """Lightweight pseudo-axis parser used during raw-data discovery.

    Kept local so the combiner does not import the full conversion module (and
    its nmrglue dependency) merely to inspect source folders.
    """
    pp = path / 'pulseprogram'
    if not pp.is_file():
        return [], []
    lines = [x.split(';', 1)[0].strip() for x in pp.read_text(errors='replace').splitlines()]
    aliases = {'vd': 'VDLIST', 'vc': 'VCLIST', 'vp': 'VPLIST', 'va': 'VALIST'}
    define_re = re.compile(r'define\s+list<[^>]+>\s+([A-Za-z_]\w*)\s*=\s*<\$([A-Za-z0-9_]+)>', re.I)
    for line in lines:
        m = define_re.search(line)
        if m:
            aliases[m.group(1)] = m.group(2).upper()
    inc = {}
    for idx, line in enumerate(lines):
        found = [a for a in aliases if re.search(r'\b%s\.(?:inc|dec)\b' % re.escape(a), line, re.I)]
        for a in ('vd', 'vc', 'vp', 'va'):
            if re.search(r'\b[di]%s\b' % a, line, re.I):
                found.append(a)
        if found:
            inc[idx] = list(dict.fromkeys(found))
    groups, covered = [], set()
    for end, line in enumerate(lines):
        m = re.search(r'\blo\s+to\s+(\d+)\s+times\s+([^\s]+)', line, re.I)
        if not m:
            continue
        label = m.group(1); start = None
        for j in range(end - 1, -1, -1):
            if re.match(r'^%s(?:\s|$)' % re.escape(label), lines[j]):
                start = j; break
        if start is None:
            continue
        group = []
        for idx in sorted(inc):
            if start <= idx < end:
                group.extend(inc[idx]); covered.add(idx)
        group = list(dict.fromkeys(group))
        if group and group not in groups:
            groups.append(group)
    for idx, group in inc.items():
        if idx not in covered:
            groups.append(group)
    clean_groups, seen = [], set()
    for group in groups:
        clean = [a for a in group if a not in seen]
        if clean:
            clean_groups.append(clean); seen.update(clean)
    values_by_name, valid_groups = {}, []
    for group in clean_groups:
        valid = []
        for alias in group:
            parameter = aliases[alias]
            configured = _scalar_text(params.get(parameter))
            candidates = [path / parameter.lower()]
            if configured:
                configured = _normalise_value(configured)
                candidates += [path / configured, path / Path(configured).name]
            lp = next((x for x in candidates if x.is_file()), None)
            if lp is None:
                continue
            vals = []
            unit = None
            for raw in lp.read_text(errors='replace').splitlines():
                token = raw.strip().split()
                if not token or raw.lstrip().startswith((';', '#')):
                    continue
                token = token[0]
                try: float(token); numeric = True
                except ValueError: numeric = False
                if not vals and unit is None and not numeric:
                    unit = token; continue
                vals.append(token)
            if vals:
                values_by_name[alias] = vals; valid.append(alias)
        if valid:
            valid_groups.append(valid)
    if not valid_groups:
        return [], []
    group_rows = []
    for group in valid_groups:
        lengths = [len(values_by_name[a]) for a in group]
        if len(set(lengths)) != 1:
            raise ValueError(f'{path.name}: synchronized pseudo-axis lists have different lengths')
        group_rows.append([list(x) for x in zip(*(values_by_name[a] for a in group))])
    rows = []
    for combo in itertools.product(*group_rows):
        rows.append([v for part in combo for v in part])
    return [a for group in valid_groups for a in group], rows


def inspect_experiment(path: Path, number: int | None = None) -> BrukerExperimentInfo:
    path = Path(path)
    acqus = path / 'acqus'
    if not acqus.is_file():
        raise ValueError(f'{path.name}: no acqus file')
    raw = path / 'fid'
    kind = 'fid'
    if not raw.is_file():
        raw = path / 'ser'
        kind = 'ser'
    if not raw.is_file():
        raise ValueError(f'{path.name}: no fid or ser file')
    params = read_jcamp(acqus)
    record = bruker_record_bytes(params)
    size = raw.stat().st_size
    if size % record:
        raise ValueError(f'{path.name}: raw size {size} is not a multiple of Bruker record size {record}')
    if number is None:
        number = int(path.name)
    row_count = size // record
    pseudo_columns: list[str] = []
    pseudo_rows: list[list[str]] = []
    if kind == 'ser':
        # Reuse decon's existing Bruker pulseprogram/list parser.  A source SER
        # is itself pseudo-2D, so its inner pseudo coordinate must be retained
        # when several experiments are flattened into the combined SER.
        pseudo_columns, pseudo_rows = _source_pseudo_axis(path, params)
    return BrukerExperimentInfo(int(number), path, raw, kind, size, record, row_count,
                                params, pseudo_columns, pseudo_rows)


def discover_numbered_experiments(parent: str | Path, start: int | None = None,
                                    finish: int | None = None) -> list[BrukerExperimentInfo]:
    """Discover numeric child folders, using an inclusive start/finish range."""
    parent = Path(parent)
    found = []
    if not parent.is_dir():
        return found
    for child in parent.iterdir():
        if not child.is_dir() or not child.name.isdigit():
            continue
        number = int(child.name)
        if start is not None and number < int(start):
            continue
        if finish is not None and number > int(finish):
            continue
        try:
            found.append(inspect_experiment(child, number))
        except ValueError:
            # Discovery deliberately ignores numeric folders which are not raw data.
            continue
    return sorted(found, key=lambda x: x.number)


def inspect_combination(experiments: list[BrukerExperimentInfo]) -> CombinationInspection:
    errors, warnings = [], []
    if not experiments:
        return CombinationInspection([], {}, {}, ['No Bruker experiments selected.'], [])
    for e in experiments:
        if e.raw_kind == 'fid' and e.row_count != 1:
            errors.append(f'Experiment {e.number}: fid contains {e.row_count} direct records; expected one.')
        if e.raw_kind == 'ser':
            if not e.pseudo_rows:
                errors.append(f'Experiment {e.number}: ser contains {e.row_count} rows but its pseudo axis could not be detected.')
            elif len(e.pseudo_rows) != e.row_count:
                errors.append(f'Experiment {e.number}: ser contains {e.row_count} raw rows but the detected pseudo axis has {len(e.pseudo_rows)} rows.')

    ref = experiments[0]
    for exp in experiments[1:]:
        if exp.record_bytes != ref.record_bytes:
            errors.append(f'Experiment {exp.number}: record size {exp.record_bytes} differs from {ref.number}: {ref.record_bytes}.')
        for key in STRUCTURAL_KEYS:
            if not _same_parameter(ref.parameters.get(key), exp.parameters.get(key), key):
                errors.append(f'Experiment {exp.number}: structural parameter {key} differs '
                              f'({ref.parameters.get(key)!r} -> {exp.parameters.get(key)!r}).')

    common_keys = set(experiments[0].parameters)
    for exp in experiments[1:]:
        common_keys &= set(exp.parameters)
    varying: dict[str, list[str]] = {}
    numeric: dict[str, list[str]] = {}
    for key in sorted(common_keys):
        vals = [e.parameters.get(key) for e in experiments]
        if all(_same_parameter(vals[0], v, key) for v in vals[1:]):
            continue
        if key in IGNORE_VARYING:
            continue
        # Keep scalar changes for the generated pseudo-axis lists. Arrays are
        # still preserved in the manifest's per-experiment metadata summary.
        texts = [_scalar_text(v) for v in vals]
        if any(v is None for v in texts):
            continue
        varying[key] = [str(v) for v in texts]
        if all(_as_float(v) is not None for v in texts):
            numeric[key] = [str(v) for v in texts]

    if not numeric:
        warnings.append('No numeric acquisition parameter varies across the selection; experiment number will be the pseudo axis.')
    return CombinationInspection(experiments, varying, numeric, errors, warnings)


def _replace_or_append_scalar(text: str, key: str, value: str) -> str:
    pattern = re.compile(r'(?m)^##\$%s=.*$' % re.escape(key))
    line = f'##${key}= {value}'
    if pattern.search(text):
        return pattern.sub(line, text, count=1)
    if text and not text.endswith('\n'):
        text += '\n'
    return text + line + '\n'


def _safe_alias(key: str, used: set[str]) -> str:
    base = re.sub(r'[^A-Za-z0-9_]', '_', key).lower().strip('_') or 'value'
    if base[0].isdigit():
        base = 'p_' + base
    alias = base[:24]
    n = 2
    while alias in used:
        suffix = f'_{n}'
        alias = base[:24-len(suffix)] + suffix
        n += 1
    used.add(alias)
    return alias


def _flatten_axis(experiments: list[BrukerExperimentInfo],
                  inspection: CombinationInspection) -> list[tuple[str, str, list[str]]]:
    """Build synchronized list columns, one value for every raw row in output SER."""
    total_rows = sum(e.row_count for e in experiments)
    columns: list[tuple[str, str, list[str]]] = [
        ('expno', 'DECON_EXPLIST', [str(e.number) for e in experiments for _ in range(e.row_count)])
    ]
    used = {'expno'}
    for key, values in inspection.numeric_varying_parameters.items():
        alias = _safe_alias(key, used)
        expanded = [value for e, value in zip(experiments, values) for _ in range(e.row_count)]
        columns.append((alias, 'DECON_' + re.sub(r'[^A-Za-z0-9_]', '_', key).upper(), expanded))

    # Preserve pseudo coordinates already present inside source SER files.
    # All generated columns are synchronized and therefore describe the exact
    # flattened row order of the concatenated binary SER.
    inner_names: list[str] = []
    for e in experiments:
        for name in e.pseudo_columns:
            if name not in inner_names:
                inner_names.append(name)
    for name in inner_names:
        alias = _safe_alias('src_' + name, used)
        vals: list[str] = []
        for e in experiments:
            if e.raw_kind == 'fid':
                vals.append('0')
                continue
            if name not in e.pseudo_columns:
                raise ValueError(f'Experiment {e.number}: source pseudo axis lacks column {name!r}.')
            idx = e.pseudo_columns.index(name)
            vals.extend(str(row[idx]) for row in e.pseudo_rows)
        if len(vals) != total_rows:
            raise ValueError(f'Internal pseudo-axis length error for {name}: {len(vals)} != {total_rows}.')
        columns.append((alias, 'DECON_SRC_' + re.sub(r'[^A-Za-z0-9_]', '_', name).upper(), vals))
    return columns


def _write_synthetic_metadata(parent: Path, experiments: list[BrukerExperimentInfo],
                              inspection: CombinationInspection) -> dict[str, str]:
    ref = experiments[0]
    acqus_text = (ref.path / 'acqus').read_text(errors='replace')
    columns = _flatten_axis(experiments, inspection)
    total_rows = sum(e.row_count for e in experiments)

    generated_lists: dict[str, str] = {}
    for alias, parameter, values in columns:
        filename = 'decon_' + alias + '.list'
        (parent / filename).write_text('\n'.join(values) + '\n')
        generated_lists[parameter] = filename
        acqus_text = _replace_or_append_scalar(acqus_text, parameter, f'<{filename}>')
    (parent / 'acqus').write_text(acqus_text)

    # The generated indirect dimension is a real pseudo dimension with one
    # point per stored direct-dimension record in the concatenated SER.
    source_acqu2s = ref.path / 'acqu2s'
    acqu2s_text = source_acqu2s.read_text(errors='replace') if source_acqu2s.is_file() else acqus_text
    acqu2s_text = _replace_or_append_scalar(acqu2s_text, 'TD', str(total_rows))
    acqu2s_text = _replace_or_append_scalar(acqu2s_text, 'FnMODE', '1')
    (parent / 'acqu2s').write_text(acqu2s_text)

    definitions = [f'define list<delay> {alias}= <${parameter}>' for alias, parameter, _ in columns]
    increments = ' '.join(f'{alias}.inc' for alias, _, _ in columns)
    pulseprogram = [
        '; Synthetic pulseprogram generated by decon Bruker combiner.',
        '; All lists are synchronized and map one-to-one to combined SER rows.',
        *definitions,
        '1 ze',
        f'2 {increments}',
        f'lo to 2 times {total_rows}',
        'exit',
    ]
    (parent / 'pulseprogram').write_text('\n'.join(pulseprogram) + '\n')
    return generated_lists


def combine_bruker_experiments(parent: str | Path, experiments: list[BrukerExperimentInfo],
                               overwrite: bool = False) -> Path:
    """Populate *parent* with ser/acqus/acqu2s/list metadata for decon."""
    parent = Path(parent)
    inspection = inspect_combination(experiments)
    if inspection.errors:
        raise ValueError('\n'.join(inspection.errors))
    targets = [parent / name for name in ('ser', 'acqus', 'acqu2s', 'pulseprogram', 'combine_manifest.json')]
    if not overwrite and any(p.exists() for p in targets):
        raise FileExistsError('Combined raw-data files already exist in the fid path.')

    parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_name = tempfile.mkstemp(prefix='.decon_ser_', dir=str(parent))
    os.close(tmp_fd)
    tmp = Path(tmp_name)
    try:
        with tmp.open('wb') as out:
            for exp in experiments:
                with exp.raw_file.open('rb') as src:
                    shutil.copyfileobj(src, out, length=1024 * 1024)
        expected = sum(e.raw_bytes for e in experiments)
        if tmp.stat().st_size != expected:
            raise IOError(f'Combined SER has {tmp.stat().st_size} bytes; expected {expected}.')
        os.replace(tmp, parent / 'ser')
        generated_lists = _write_synthetic_metadata(parent, experiments, inspection)

        manifest = {
            'format': 'decon-bruker-combination-v1',
            'destination': str(parent.resolve()),
            'record_bytes': experiments[0].record_bytes,
            'rows': sum(e.row_count for e in experiments),
            'varying_parameters': inspection.varying_parameters,
            'pseudo_axis_numeric_parameters': inspection.numeric_varying_parameters,
            'generated_lists': generated_lists,
            'warnings': inspection.warnings,
            'sources': [
                {
                    'row_start': 1 + sum(x.row_count for x in experiments[:i]),
                    'row_end': sum(x.row_count for x in experiments[:i + 1]),
                    'rows': e.row_count,
                    'experiment': e.number,
                    'directory': str(e.path.resolve()),
                    'raw_file': e.raw_file.name,
                    'raw_bytes': e.raw_bytes,
                    'source_pseudo_columns': e.pseudo_columns,
                    'source_pseudo_rows': e.pseudo_rows,
                    'parameters': {k: e.parameters.get(k) for k in sorted(set(STRUCTURAL_KEYS) | set(inspection.varying_parameters))},
                }
                for i, e in enumerate(experiments)
            ],
        }
        (parent / 'combine_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    finally:
        if tmp.exists():
            tmp.unlink()
    return parent
