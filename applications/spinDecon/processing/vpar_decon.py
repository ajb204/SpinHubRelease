#!/usr/bin/python
#####################################################
# Functions to work out an fid.com conversion script
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
import os,sys,numpy
import re
from datetime import datetime
from tabnanny import check
import shutil
import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import nmrglue as ng



@dataclass
class PseudoAxisInfo:
    """Description of arrayed acquisition parameters forming a pseudo axis."""
    size: int = 1
    columns: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)
    parameters: list[dict[str, Any]] = field(default_factory=list)
    groups: list[list[str]] = field(default_factory=list)


def _read_bruker_list(path: Path) -> tuple[list[str], str | None]:
    """Read a TopSpin list file, preserving values as text and an optional unit header."""
    values, unit = [], None
    try:
        lines = path.read_text(errors='replace').splitlines()
    except OSError:
        return values, unit
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith((';', '#')):
            continue
        token = line.split()[0]
        # Power lists commonly start with dB; other list files can similarly
        # contain a non-numeric unit/type marker.  It is metadata, not a point.
        try:
            float(token)
            numeric = True
        except ValueError:
            numeric = False
        if not values and unit is None and not numeric:
            unit = token
            continue
        values.append(token)
    return values, unit


def _cartesian_group_rows(groups, values_by_name):
    """Zip parameters within a group; take the Cartesian product between groups."""
    group_rows = []
    for group in groups:
        lengths = [len(values_by_name[name]) for name in group]
        if not lengths or min(lengths) == 0:
            continue
        if len(set(lengths)) != 1:
            raise ValueError('Pseudo-axis parameters incremented in the same loop have different lengths: %s' %
                             ', '.join('%s=%d' % (n, len(values_by_name[n])) for n in group))
        group_rows.append([list(x) for x in zip(*(values_by_name[name] for name in group))])
    if not group_rows:
        return []
    rows = []
    for combination in itertools.product(*group_rows):
        row = []
        for part in combination:
            row.extend(part)
        rows.append(row)
    return rows


def detect_bruker_pseudo_axis(path) -> PseudoAxisInfo:
    """Detect Bruker pseudo-axis lists from pulseprogram increment/loop syntax.

    Parameters incremented in the same ``lo ... times`` loop are synchronized
    (columns in one array).  Parameters incremented in independent/nested loops
    form Cartesian dimensions, so their lengths multiply.
    """
    base = Path(path)
    pp = base / 'pulseprogram'
    if not pp.is_file():
        return PseudoAxisInfo()
    raw_lines = pp.read_text(errors='replace').splitlines()
    # Strip TopSpin ';' comments before syntax matching.
    lines = [line.split(';', 1)[0].strip() for line in raw_lines]

    alias_to_parameter = {'vd': 'VDLIST', 'vc': 'VCLIST', 'vp': 'VPLIST', 'va': 'VALIST'}
    # e.g. define list<power> protCW= <$VALIST>
    define_re = re.compile(r'define\s+list<[^>]+>\s+([A-Za-z_]\w*)\s*=\s*<\$([A-Za-z0-9_]+)>', re.I)
    for line in lines:
        m = define_re.search(line)
        if m:
            alias_to_parameter[m.group(1)] = m.group(2).upper()

    increments = {}
    for idx, line in enumerate(lines):
        found = []
        for alias in alias_to_parameter:
            if re.search(r'\b%s\.(?:inc|dec)\b' % re.escape(alias), line, re.I):
                found.append(alias)
        # Built-in list pointer operations: ivd/dvd, ivc/dvc, etc.
        for alias in ('vd', 'vc', 'vp', 'va'):
            if re.search(r'\b[di]%s\b' % alias, line, re.I):
                found.append(alias)
        if found:
            increments[idx] = list(dict.fromkeys(found))

    # Associate increments with the innermost Bruker loop body that contains them.
    loop_groups = []
    covered = set()
    loop_re = re.compile(r'\blo\s+to\s+(\d+)\s+times\s+([^\s]+)', re.I)
    for end, line in enumerate(lines):
        m = loop_re.search(line)
        if not m:
            continue
        label = m.group(1)
        start = None
        for j in range(end - 1, -1, -1):
            if re.match(r'^%s(?:\s|$)' % re.escape(label), lines[j]):
                start = j
                break
        if start is None:
            continue
        aliases = []
        for idx in sorted(increments):
            if start <= idx < end:
                aliases.extend(increments[idx])
                covered.add(idx)
        aliases = list(dict.fromkeys(aliases))
        if aliases and aliases not in loop_groups:
            loop_groups.append(aliases)
    for idx, aliases in increments.items():
        if idx not in covered:
            loop_groups.append(aliases)

    # Remove duplicate aliases from outer loops: an increment belongs to its
    # innermost loop.  Identical synchronized groups collapse naturally.
    groups, seen = [], set()
    for group in loop_groups:
        clean = [a for a in group if a not in seen]
        if clean:
            groups.append(clean)
            seen.update(clean)

    values_by_name, params = {}, []
    valid_groups = []
    for group in groups:
        valid = []
        for alias in group:
            parameter = alias_to_parameter[alias]
            configured = _bruker_scalar(base / 'acqus', parameter)
            candidates = [base / parameter.lower()]
            if configured:
                candidates.extend((base / configured, base / Path(configured).name))
            list_path = next((x for x in candidates if x.is_file()), None)
            if list_path is None:
                continue
            values, unit = _read_bruker_list(list_path)
            if not values:
                continue
            values_by_name[alias] = values
            params.append({'name': alias, 'parameter': parameter, 'configured_file': configured,
                           'file': str(list_path), 'unit': unit, 'values': values})
            valid.append(alias)
        if valid:
            valid_groups.append(valid)
    rows = _cartesian_group_rows(valid_groups, values_by_name) if valid_groups else []
    columns = [name for group in valid_groups for name in group]
    return PseudoAxisInfo(size=len(rows) if rows else 1, columns=columns, rows=rows,
                          parameters=params, groups=valid_groups)


def _parse_varian_array_expression(expression: str) -> list[list[str]]:
    """Parse VNMR array syntax, including synchronized groups such as (a,b),c."""
    expression = expression.strip().strip('"')
    groups, token, depth = [], '', 0
    for ch in expression:
        if ch == '(':
            depth += 1
            if depth == 1:
                continue
        elif ch == ')':
            depth -= 1
            if depth == 0:
                if token.strip():
                    groups.append([x.strip() for x in token.split(',') if x.strip()])
                token = ''
                continue
        if ch == ',' and depth == 0:
            if token.strip():
                groups.append([token.strip()])
            token = ''
        else:
            token += ch
    if token.strip():
        groups.append([token.strip()])
    return groups


def detect_varian_pseudo_axis(path) -> PseudoAxisInfo:
    """Detect VNMR pseudo-axis arrays directly from procpar's ``array`` parameter."""
    base = Path(path)
    procpar = base / 'procpar'
    array_values = _varian_values(procpar, 'array')
    if not array_values:
        return PseudoAxisInfo()
    groups = _parse_varian_array_expression(array_values[0])
    # phase/phase2/phase3 encode spectral quadrature, not pseudo parameters.
    groups = [[name for name in group if not re.fullmatch(r'phase\d*', name, re.I)] for group in groups]
    groups = [group for group in groups if group]
    values_by_name, params, valid_groups = {}, [], []
    for group in groups:
        valid = []
        for name in group:
            values = _varian_values(procpar, name)
            if not values:
                continue
            values_by_name[name] = values
            params.append({'name': name, 'parameter': name, 'file': str(procpar), 'unit': None, 'values': values})
            valid.append(name)
        if valid:
            valid_groups.append(valid)
    rows = _cartesian_group_rows(valid_groups, values_by_name) if valid_groups else []
    columns = [name for group in valid_groups for name in group]
    return PseudoAxisInfo(size=len(rows) if rows else 1, columns=columns, rows=rows,
                          parameters=params, groups=valid_groups)


@dataclass(frozen=True)
class AcquisitionInfo:
    """Read-only description of a spectrometer acquisition directory.

    ``vendor`` deliberately uses the established vpar codes (``bruk`` and
    ``var``), so existing Decon callers do not need a translation layer.
    Metadata is best-effort: malformed or old parameter files still identify
    as acquisitions when their structural markers are valid.
    """
    path: Path
    vendor: str
    fid_file: Path
    compressed: bool = False
    dimension: int | None = None
    sequence: str | None = None
    nuclei: tuple[str, ...] = ()
    observation_frequency_mhz: float | None = None
    temperature_k: float | None = None
    acquired_at: str | None = None
    acquisition_time: str | None = None
    pseudo_axis: bool = False
    pseudo_axis_size: int | None = None
    pseudo_axis_columns: tuple[str, ...] = ()
    dimension_confidence: str = 'best-effort'
    dimension_reason: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


def detect_spectrometer_type(path='.'):
    """Return the established vpar vendor code, without side effects.

    The markers mirror the formats already understood by :class:`vpar`.
    ``acqus`` is the normal modern Bruker marker; ``acqu``/``acqu2`` retain
    compatibility with legacy data handled by the original code.
    """
    base = Path(path)
    if any((base / name).is_file() for name in ('acqu', 'acqu2', 'acqus')):
        return 'bruk'
    if (base / 'procpar').is_file():
        return 'var'
    return None


def _bruker_scalar(path: Path, key: str) -> str | None:
    """Read a scalar JCAMP-DX value from a Bruker acquisition file."""
    try:
        prefix = f'##${key}='
        for line in path.read_text(errors='replace').splitlines():
            if line.startswith(prefix):
                value = line[len(prefix):].strip()
                if value.startswith('<') and value.endswith('>'):
                    value = value[1:-1]
                return value
    except OSError:
        pass
    return None


def _varian_values(path: Path, key: str) -> list[str]:
    """Read values using the same procpar record layout as vpar.GetParVarian."""
    try:
        rows = [line.split() for line in path.read_text(errors='replace').splitlines()]
    except OSError:
        return []
    for i, row in enumerate(rows[:-1]):
        if row and row[0] == key:
            values = rows[i + 1]
            if not values:
                return []
            try:
                count = int(values[0])
            except (ValueError, TypeError):
                return []
            result = values[1:1 + count]
            # Long procpar arrays can continue on following physical lines.
            j = i + 2
            while len(result) < count and j < len(rows):
                result.extend(rows[j][:count - len(result)])
                j += 1
            return [v.strip('"') for v in result]
    return []


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalise_nucleus(value: str | None) -> str | None:
    if not value:
        return None
    value = value.strip().strip('<>').strip('"')
    # Keep the spectrometer spelling in metadata, but expose the convention
    # historically used by Decon (1H, 13C, 15N) where possible.
    m = re.fullmatch(r'([A-Za-z]+)(\d+)', value)
    if m:
        return f'{m.group(2)}{m.group(1)}'
    return value


def _format_elapsed(seconds) -> str | None:
    try:
        minutes = max(0, int(round(float(seconds) / 60.0)))
    except (TypeError, ValueError):
        return None
    days, rem = divmod(minutes, 1440); hours, mins = divmod(rem, 60)
    parts = []
    if days: parts.append(f"{days} day" + ("s" if days != 1 else ""))
    if hours: parts.append(f"{hours} hour" + ("s" if hours != 1 else ""))
    if mins or not parts: parts.append(f"{mins} minute" + ("s" if mins != 1 else ""))
    return ', '.join(parts)

def _actual_acquisition_time(base: Path, vendor: str) -> str | None:
    if vendor == 'bruk':
        audit = base / 'audita.txt'
        if audit.is_file():
            start, end = vpar._bruker_completed_acquisition(str(audit))
            if start is not None and end is not None:
                return _format_elapsed((end-start).total_seconds())
    else:
        log = base / 'log'
        if log.is_file():
            stamps=[]
            for line in log.read_text(errors='replace').splitlines():
                dt=vpar._parse_audit_datetime(line) or vpar._parse_audit_datetime(line.split(': ',1)[0])
                if dt is not None: stamps.append(dt)
            if len(stamps) >= 2: return _format_elapsed((max(stamps)-min(stamps)).total_seconds())
    return None

def _bruker_indirect_points(f: Path) -> int | None:
    # NusTD is the most useful indication for NUS acquisitions; otherwise TD.
    raw = _bruker_scalar(f, 'NusTD') or _bruker_scalar(f, 'TD')
    val = _as_float(raw)
    return max(0, int(round(val))) if val is not None else None

def _inspect_bruker(base: Path) -> dict[str, Any]:
    acqus = base / 'acqus'
    primary = acqus if acqus.is_file() else next((base/n for n in ('acqu','acqu2') if (base/n).is_file()), acqus)
    nuclei=[]; indirect=[]
    for index, filename in enumerate(('acqus','acqu2s','acqu3s','acqu4s'), start=1):
        f=base/filename
        if not f.is_file(): continue
        nuc=_normalise_nucleus(_bruker_scalar(f,'NUC1'))
        if nuc: nuclei.append(nuc)
        if index > 1: indirect.append((index, _bruker_indirect_points(f)))
    try: pseudo=detect_bruker_pseudo_axis(base)
    except Exception: pseudo=PseudoAxisInfo()
    pseudo_size=pseudo.size if pseudo.size and pseudo.size > 1 else None
    active=[(i,n) for i,n in indirect if n is not None and n > 1]
    # If an indirect acquisition count matches a detected array length, treat that
    # physical dimension as pseudo rather than as another spectral dimension.
    spectral_active=list(active)
    if pseudo_size:
        matches=[x for x in spectral_active if x[1] in (pseudo_size, pseudo_size*2)]
        if matches: spectral_active.remove(matches[-1])
    dimension=1+len(spectral_active)
    if dimension == 1 and ((base/'ser').is_file() or (base/'ser.gz').is_file()) and not indirect:
        dimension=2
    reason='indirect TD/NusTD increments > 1'
    if pseudo_size: reason += f'; detected pseudo axis with {pseudo_size} points'
    return {'dimension':dimension,'sequence':_bruker_scalar(primary,'PULPROG'),'nuclei':tuple(nuclei),
        'observation_frequency_mhz':_as_float(_bruker_scalar(primary,'SFO1')),'temperature_k':_as_float(_bruker_scalar(primary,'TE')),
        'acquired_at':_bruker_scalar(primary,'DATE'),'acquisition_time':_actual_acquisition_time(base,'bruk'),
        'pseudo_axis':bool(pseudo_size),'pseudo_axis_size':pseudo_size,'pseudo_axis_columns':tuple(pseudo.columns),
        'dimension_reason':reason,'metadata':{'parameter_file':str(primary),'indirect_points':indirect}}


def _inspect_varian(base: Path) -> dict[str, Any]:
    procpar = base / 'procpar'
    def first(key):
        vals = _varian_values(procpar, key)
        return vals[0] if vals else None
    dimension = 1
    for index, key in ((2, 'ni'), (3, 'ni2'), (4, 'ni3')):
        val = _as_float(first(key))
        if val is not None and val > 1:
            dimension = index
    if (base / 'ser').is_file() or (base / 'ser.gz').is_file():
        dimension = max(dimension, 2)
    nuclei = tuple(n for n in (
        _normalise_nucleus(first('tn')),
        _normalise_nucleus(first('dn')),
        _normalise_nucleus(first('dn2')),
        _normalise_nucleus(first('dn3')),
    ) if n)
    temp_c = _as_float(first('temp'))
    try: pseudo = detect_varian_pseudo_axis(base)
    except Exception: pseudo = PseudoAxisInfo()
    pseudo_size = pseudo.size if pseudo.size and pseudo.size > 1 else None
    return {
        'dimension': dimension,
        'sequence': first('seqfil'),
        'nuclei': nuclei,
        'observation_frequency_mhz': _as_float(first('sfrq')),
        'temperature_k': temp_c + 273.15 if temp_c is not None else None,
        'acquired_at': first('time_complete'),
        'acquisition_time': _actual_acquisition_time(base, 'var'),
        'pseudo_axis': bool(pseudo_size),
        'pseudo_axis_size': pseudo_size,
        'pseudo_axis_columns': tuple(pseudo.columns),
        'dimension_reason': 'ni/ni2/ni3 values > 1',
        'metadata': {'parameter_file': str(procpar)},
    }


def inspect_acquisition(path='.'):
    """Inspect an acquisition without modifying it or changing cwd.

    Recognition and metadata extraction are intentionally tolerant. A broken
    optional metadata field must not make an otherwise valid acquisition
    disappear from SpinHub.
    """
    base = Path(path)
    fid_file = next((base / name for name in ('fid', 'ser', 'fid.gz', 'ser.gz')
                     if (base / name).is_file()), None)
    if fid_file is None:
        return None
    vendor = detect_spectrometer_type(base)
    if vendor is None:
        return None
    details = _inspect_bruker(base) if vendor == 'bruk' else _inspect_varian(base)
    return AcquisitionInfo(
        path=base,
        vendor=vendor,
        fid_file=fid_file,
        compressed=fid_file.suffix == '.gz',
        **details,
    )


def find_child_acquisitions(path='.', max_results=12):
    """Return immediate child acquisition summaries when *path* is one level too high."""
    base=Path(path)
    if not base.is_dir(): return []
    found=[]
    try: children=sorted((p for p in base.iterdir() if p.is_dir()), key=lambda p:p.name.lower())
    except OSError: return []
    for child in children:
        info=inspect_acquisition(child)
        if info is not None:
            found.append(info)
            if len(found) >= max_results: break
    return found


class vpar():
    def __init__(self, *args, **kwargs):
        self._reset_runtime_state()
        if args or kwargs:
            self.Setup(*args, **kwargs)

    def _reset_runtime_state(self):
        self.parent = None
        self.outdir = ''
        self.FidPath = ''
        self.AMX = None
        self.brukerFmt = None
        self.dim = None
        self.labb = []
        self.rk = []
        self.nuslist = ''
        self.o1p = 'auto'
        self.abort = 0
        self.tp = None
        self.parfile = ''
        self.seqfil = ''
        self._pipe_rows = []
        self.phasing = False
        self.spa = 10
        self.axis = ('x', 'y', 'z', 'a')
        self.initialized = False
        self.setup_error = None
        # Canonical NMRPipe time-domain sizes.  These are populated only after
        # the vendor-specific conversion path has resolved acquisition order,
        # NUS reconstruction sizes and pseudo-dimensional handling.
        self.xN = None
        self.yN = None
        self.zN = None
        self.aN = None
        # Canonical NMRPipe observation frequencies (MHz), published after
        # vendor-specific frequency and acquisition-order handling.
        self.xOBS = None
        self.yOBS = None
        self.zOBS = None
        self.aOBS = None
        self.pseudo_axis_info = PseudoAxisInfo()

    def Setup(self,parent,outdir,dim,labb,rk,nuslist='',o1p='auto',AMX=None,process=True,FidPath=''):
        self._reset_runtime_state()
        self.parent=parent
        self.outdir=outdir  #folder to dump outputs.
        if(FidPath==''):  #set path to raw data.
            self.FidPath=outdir
        else:
            self.FidPath=FidPath
        #print outdir
        #print os.path.exists(os.path.join(outdir,'ser'))
        self.AMX=AMX
        self.brukerFmt = None
        #possible names of raw fid file
        testfiles=[]
        testfiles.append('fid')
        testfiles.append('ser')
        testfiles.append('fid.gz')
        testfiles.append('ser.gz')
        pass
        tick=0
        for test in testfiles:
            pass
            if(os.path.exists(os.path.join(self.FidPath,test))==True):
               tick=1
               infid=os.path.join(self.FidPath,test)
               pass
               break
        if(tick==0):
            self.setup_error = 'missing_raw_data'
            return self

        if(infid[-3:]=='.gz'): #fid is zipped up...
            #unzipping.
            os.system('gunzip '+infid)




        if(os.path.exists(os.path.join(self.FidPath,'procpar'))==0 and os.path.exists(os.path.join(self.FidPath,'acqus'))==0):
            self.setup_error = 'missing_source_metadata'
            return self
        self.dim=dim
        self.labb=labb
        self.rk=rk
        self.nuslist=nuslist
        self.o1p=o1p

        self.abort=0
        # Default nucleus attributes so partial Bruker metadata does not
        # trigger attribute errors in later conversion/reference steps.
        self.NUC1 = None
        self.NUC1_raw = None
        self.NUC2 = None
        self.NUC2_raw = None
        self.NUC3 = None
        self.NUC3_raw = None
        self.NUC4 = None
        self.NUC4_raw = None
        self.n1 = None
        self.n2 = None
        self.n3 = None
        self.n4 = None
        # Spectrometer detection is owned by ProcessFrame because FidPath is
        # fixed for that window. Reuse its cached result; retain a fallback for
        # non-GUI/legacy callers that construct vpar directly.
        cached_tp = getattr(parent, 'tp', None)
        if cached_tp in ('bruk', 'var', 'omega'):
            self.tp = cached_tp
            # GetSpectrometerType() normally initialises vendor-specific input
            # paths as a side effect.  When the GUI supplies a cached vendor we
            # skip that function, so initialise the Varian procpar path here.
            # Without this, GetParVarian() attempts to open self.parfile == ''
            # when Show Script/Guess builds a fresh conversion script.
            if self.tp == 'var':
                self.parfile = os.path.join(self.FidPath, 'procpar')
        else:
            self.GetSpectrometerType(path=self.FidPath)
        self.GetSequence()
        self.initialized = True
        return self

    ##############################################################
    #compute and run fid.test.com
    def Process(self):
        script_path = self.BuildConversionScript()
        if script_path == -1:
            return

        os.system('csh ' + script_path)
        pass
        if self._slice_dim_count() >= 1:
            self.ProcessSlice()



    def _slice_dim_count(self):
        if type(self.dim) == str:
            try:
                return int(str(self.dim).split('p')[0])
            except Exception:
                return 0
        try:
            return int(self.dim)
        except Exception:
            return 0

    def _needs_slice_conversion(self):
        return self._slice_dim_count() >= 2

    def _dim_at_least(self, count):
        """Compare legacy numeric/``2p``/``3p`` dimensions safely.

        ``self.dim`` is still a legacy wire value in this backend.  Python 3
        cannot order strings such as ``'2p'`` against integers, so every
        dimensional threshold test must pass through the normalized physical
        dimension count.
        """
        return self._slice_dim_count() >= int(count)

    def _dim_greater_than(self, count):
        return self._slice_dim_count() > int(count)

    def BuildConversionScript(self):
        """Generate fid.test.com without executing it."""
        self.Convert()
        if self.abort == 1:
            pass
            return -1
        if self.outdir:
            os.makedirs(self.outdir, exist_ok=True)
        if type(self.dim) == str and self.dim in ('2p', '3p'):
            self._write_pseudo_axis_table()
        self.PipeParse()
        return os.path.join(self.outdir, 'fid.test.com')

    def _write_pseudo_axis_table(self, filename='pseudo_axis.tsv'):
        """Write the conversion-defined pseudo axis as TSV and CSV.

        The TSV remains the authoritative legacy file.  A CSV twin is emitted
        from the same rows so downstream tools which expect comma-separated
        metadata see exactly the same spectrum-to-source mapping.
        """
        info = self.pseudo_axis_info
        if not info.rows or not self.outdir:
            return None
        path = Path(self.outdir) / filename
        with path.open('w', encoding='utf-8') as out:
            out.write('spectrum\t' + '\t'.join(info.columns) + '\n')
            for number, row in enumerate(info.rows, start=1):
                out.write(str(number) + '\t' + '\t'.join(map(str, row)) + '\n')

        import csv
        csv_path = path.with_suffix('.csv')
        with csv_path.open('w', encoding='utf-8', newline='') as out:
            writer = csv.writer(out)
            writer.writerow(['spectrum'] + list(info.columns))
            for number, row in enumerate(info.rows, start=1):
                writer.writerow([number] + list(row))
        return str(path)

    def _infer_bruker_amx(self):
        """Infer the Bruker digitizer family using acquisition metadata.

        NMRPipe's bruk2pipe defaults to AMX, with DMX reserved for explicit
        legacy DMX data. We follow the same rule here and only select DMX when
        the acquisition instrument string clearly indicates a DMX-family console.
        """
        if isinstance(self.AMX, bool):
            return self.AMX

        try:
            instr = str(self.GetParBruk('acqus', ('', 'INSTRUM'))[0])
        except Exception:
            instr = ''
        instr = instr.strip().strip('<>').upper()

        # Legacy DMX-family instruments should explicitly request DMX.
        if any(tag in instr for tag in (' DMX', 'DMX', ' DRX', 'DRX', ' DPX', 'DPX')):
            return False

        # Everything else defaults to AMX, matching bruk2pipe's default.
        return True

    def _bruker_format_flag(self):
        return 'AMX' if self._infer_bruker_amx() else 'DMX'

    def _set_bruker_format(self):
        amxtext = self._bruker_format_flag()
        self.AMX = (amxtext == 'AMX')
        self.brukerFmt = amxtext
        return amxtext

    def _write_bruker_bad_line(self, flag, amxtext, decim, dspfvs, grpdly, dtext='', decim_as_int=False):
        decim_text = '%i' % decim if decim_as_int else '%f' % decim
        extra = (' %s' % dtext) if dtext else ''
        self.outy.write(' -bad 0.0 -%s -%s -decim %s -dspfvs %i -grpdly %f%s' % (flag, amxtext, decim_text, dspfvs, grpdly, extra))
        self.EndPipeLine(self.outy)

    def _pipe_parse_varian_3p_slice(self):
        """Write the 2D first-plane preview conversion for Varian 3p data."""
        self.phasing = False
        self.spa = 10
        self.axis = ('x', 'y', 'z', 'a')
        self._pipe_rows = []

        outpath = os.path.join(self.outdir, 'fid.test.slice.com')
        self.outy = open(outpath, 'w')
        self.outy.write('#!/bin/csh\n\n')
        self.outy.write('set ft4trec=%s/slice.fid\n' % self.outdir)
        self.outy.write('if( -e %s/slice.fid) rm -rf %s/slice.fid\n' %
                        (self.outdir, self.outdir))

        # Match the main Varian 3p conversion: RelaxFix reads raw/fid and puts
        # its processed output in spec/fid.final; var2pipe reads that file.
        relax_out = os.path.join(self.outdir, 'fid.final')
        relax_in = os.path.join(self.FidPath, 'fid')
        self.outy.write('RelaxFix.out {} {} {} 0 {} {}\n\n'.format(
            self.np, self.ni, self.nz, relax_out, relax_in))
        self.outy.write('var2pipe -in %s \\\n' % relax_out)
        self.outy.write(' -noaswap ')
        if self.acqORD != '':
            self.outy.write(' -aqORD %i \\\n' % self.acqORD)
        self.outy.write('\\\n')

        # Keep the direct dimension unchanged and retain exactly one complex
        # point in the indirect spectral dimension: yN=2, yT=1.
        # Select the actual spectral indirect axis (ni or ni2), just as the
        # main 3p conversion does.
        if self.ni > 1:
            ysw, yobs, ycar = self.sw1, self.dfrq, self.f1ppm
        elif self.ni2 > 1:
            ysw, yobs, ycar = self.sw2, self.dfrq, self.f1ppm
        else:
            ysw, yobs, ycar = self.sw1, self.dfrq, self.f1ppm

        N = (self.np, 2)
        T = (self.np2, 1)
        M = ('Complex', 'Complex')
        sw = (self.sw, ysw)
        O = (self.sfrq, yobs)
        C = (self.waterppm, ycar)
        labels = tuple(self.labb[:2])

        # AddPipe normally derives the number of rows from self.dim.  This
        # object is still a 3p dataset, but the preview script deliberately
        # describes only X and Y.  Temporarily present it as a 2D conversion
        # while writing these rows; otherwise AddPipe asks for vals[2] and
        # raises IndexError before the buffered script is flushed/closed.
        saved_dim = self.dim
        try:
            self.dim = 2
            self.AddPipe(self.outy, self.axis, 'N', N, self.spa)
            self.AddPipe(self.outy, self.axis, 'T', T, self.spa)
            self.AddPipe(self.outy, self.axis, 'MODE', M, self.spa)
            self.AddPipe(self.outy, self.axis, 'SW', sw, self.spa)
            self.AddPipe(self.outy, self.axis, 'OBS', O, self.spa)
            self.AddPipe(self.outy, self.axis, 'CAR', C, self.spa)
            self.AddPipe(self.outy, self.axis, 'LAB', labels, self.spa)
        finally:
            self.dim = saved_dim
        self.outy.write(' -ndim  %s -aq2D  %s \\\n' %
                        (str(2).ljust(self.spa), 'States'.ljust(self.spa)))
        self.outy.write('  -out $ft4trec -verb -ov\n')
        self.outy.close()
        return outpath

    def PipeParseSlice(self):
        if not self._needs_slice_conversion():
            return

        # A Varian 3p acquisition contains two spectral dimensions followed by
        # one real pseudo axis.  The preview must be converted from the raw
        # acquisition in exactly the same way as the main conversion, but with
        # only the first point of the indirect spectral dimension retained.
        # Do not describe the pseudo axis to var2pipe at all: this is a 2D
        # preview, so there must be no zN/zT/zMODE/zSW/zOBS/zCAR/zLAB fields.
        if self.tp == 'var' and self.dim == '3p':
            return self._pipe_parse_varian_3p_slice()

        self.phasing = False
        self.spa = 10
        self.axis = 'x', 'y', 'z', 'a'
        self._pipe_rows = []

        outpath = os.path.join(self.outdir, 'fid.test.slice.com')
        self.outy = open(outpath, 'w')
        self.outy.write('#!/bin/csh\n\n')

        if self.tp == 'bruk':
            self.outy.write('cp %s/acqus ./\n' % self.FidPath)

        self.outy.write('set ft4trec=%s/slice.fid\n' % self.outdir)
        self.outy.write('if( -e %s/slice.fid) rm -rf %s/slice.fid\n' % (self.outdir, self.outdir))

        if self.tp == 'var':
            infile = os.path.join(self.FidPath, 'fid')
            self.outy.write('var2pipe -in %s \\\n' % infile)
            self.outy.write(' -noaswap ')
            if self.acqORD != '':
                self.outy.write(' -aqORD %i \\\n' % (self.acqORD))
            self.outy.write('\\\n')
        elif self.tp == 'bruk':
            if self.dim == 1 or self.dim == '1p':
                infile = os.path.join(self.FidPath, 'fid')
            else:
                infile = os.path.join(self.FidPath, 'ser')
            self.outy.write('bruk2pipe -in %s  \\\n' % infile)
            GRPDLY = getattr(self, 'GRPDLY', float(self.GetParBruk('acqus', ('', 'GRPDLY'))[0]))
            DSPFVS = getattr(self, 'DSPFVS', float(self.GetParBruk('acqus', ('', 'DSPFVS'))[0]))
            DECIM = getattr(self, 'DECIM', float(self.GetParBruk('acqus', ('', 'DECIM'))[0]))

            BYTORDA = getattr(self, 'BYTORDA', int(self.GetParBruk('acqus', ('', 'BYTORDA'))[0]))
            if BYTORDA == 1:
                flag = 'noaswap'
            else:
                flag = 'aswap'

            DTYPA = getattr(self, 'DTYPA', int(self.GetParBruk('acqus', ('', 'DTYPA'))[0]))
            dtext = getattr(self, 'dtext', '')
            if DTYPA == 2 and dtext == '':
                dtext = '-ws 8 -noi2f'

            AMXtext = self._set_bruker_format()

            self._write_bruker_bad_line(flag, AMXtext, DECIM, DSPFVS, GRPDLY, dtext=dtext)
        else:
            infile = os.path.join(self.FidPath, 'fid')
            self.outy.write('var2pipe -in %s \\\n' % infile)
            self.outy.write(' -noaswap \\\n')

        M = []
        if self.tp == 'bruk':
            M.append('DQD')
        else:
            M.append('Complex')

        self.modDict = {}
        self.modDict['0'] = 'Complex'
        self.modDict['1'] = 'QF'
        self.modDict['2'] = 'QSEQ'
        self.modDict['3'] = 'TPPI'
        self.modDict['4'] = 'States'
        self.modDict['5'] = 'States-TPPI'
        self.modDict['6'] = 'Echo-Antiecho'

        for i in range(len(self.rk)):
            if self.rk[i] == 0:
                if self.tp == 'bruk':
                    M.append(self.modDict[self.mode[i]])
                else:
                    M.append('Complex')
            else:
                M.append('Rance-Kay')
        if type(self.dim) == str:
            # Keep the final pseudo axis explicitly real; do not allow extra
            # rk entries to turn it into a spectral quadrature axis.
            pseudo_ndim = int(self.dim.split('p')[0])
            M = M[:pseudo_ndim - 1] + ['Real']

        nd = self._slice_dim_count()
        if nd == 2:
            N = (self.np, 2)
            T = (self.np2, 1)
            sw = (self.sw, self.sw1)
            O = (self.sfrq, self.frq1)
            C = (self.waterppm, self.f1ppm)
        elif nd == 3:
            N = (self.np, 2, 2)
            T = (self.np2, 1, 1)
            sw = (self.sw, self.sw1, self.sw2)
            O = (self.sfrq, self.frq1, self.frq2)
            C = (self.waterppm, self.f1ppm, self.f2ppm)
        elif nd == 4:
            N = (self.np, 2, 2, 2)
            T = (self.np2, 1, 1, 1)
            sw = (self.sw, self.sw1, self.sw2, self.sw3)
            O = (self.sfrq, self.frq1, self.frq2, self.frq3)
            C = (self.waterppm, self.f1ppm, self.f2ppm, self.f3ppm)
        else:
            self.outy.close()
            return

        if type(self.dim) == str:
            if self.dim == '2p':
                N = (self.np, 2)
                T = (self.np2, 1)
                sw = (self.sw, self.sw)
                O = (self.sfrq, self.sfrq)
                C = (self.waterppm, self.waterppm)
            elif self.dim == '3p':
                N = (self.np, 2, 2)
                T = (self.np2, 1, 1)
                sw = (self.sw, self.sw, self.sw)
                O = (self.sfrq, self.sfrq, self.sfrq)
                C = (self.waterppm, self.waterppm, self.waterppm)

        if type(self.dim) == str:
            if self.tp == 'bruk' and self.aqseq == '312' and len(M) >= 3:
                N = (N[0], N[2], N[1])
                T = (T[0], T[2], T[1])
                sw = (sw[0], sw[2], sw[1])
                O = (O[0], O[2], O[1])
                C = (C[0], C[2], C[1])
                M = (M[0], M[2], M[1])
                self.labb = (self.labb[0], self.labb[2], self.labb[1])

        self.AddPipe(self.outy, self.axis, 'N', N, self.spa)
        self.AddPipe(self.outy, self.axis, 'T', T, self.spa)
        self.AddPipe(self.outy, self.axis, 'MODE', M, self.spa)
        self.AddPipe(self.outy, self.axis, 'SW', sw, self.spa)
        self.AddPipe(self.outy, self.axis, 'OBS', O, self.spa)
        self.AddPipe(self.outy, self.axis, 'CAR', C, self.spa)
        self.AddPipe(self.outy, self.axis, 'LAB', self.labb, self.spa)

        if nd in (2, 3):
            self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(nd).ljust(self.spa), 'States'.ljust(self.spa)))
            self.outy.write('  -out $ft4trec -verb -ov\n')
        elif nd == 4:
            self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(nd).ljust(self.spa), 'States'.ljust(self.spa)))
            self.outy.write('| pipe2xyz -x -out $ft4trec -verb -ov -to 0\n')
        else:
            self.outy.write('  -out $ft4trec -verb -ov\n')

        self.outy.close()
        return

    def _prepare_1d_slice_input(self):
        src = os.path.join(self.outdir, 'test.fid')
        src_gz = src + '.gz'
        dst = os.path.join(self.outdir, 'slice.fid')
        dst_gz = dst + '.gz'

        if os.path.exists(dst):
            try:
                os.remove(dst)
            except Exception:
                pass
        if os.path.exists(dst_gz):
            try:
                os.remove(dst_gz)
            except Exception:
                pass

        if os.path.exists(src):
            shutil.copy2(src, dst)
            pass
            return dst

        if os.path.exists(src_gz):
            import gzip
            with gzip.open(src_gz, 'rb') as fin, open(dst, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
            pass
            return dst

        raise FileNotFoundError(f'Could not find 1D conversion output {src!r} or {src_gz!r}')

    def _prepare_bruker_pseudo_y_slice_input(self):
        """Prepare records 1/NZ+1 and a 2D bruk2pipe conversion script.

        The normal ``fid.test.com`` is authoritative for Bruker conversion
        parameters and for the AQSEQ-resolved physical axes.  In the special
        X-pseudoY-spectralZ layout, the spectral indirect axis is physical Z
        in the raw SER.  We therefore copy raw records 1 and NZ+1 byte-for-byte
        and build a 2D conversion using X parameters from the main script and
        the main script's Z parameters renamed to Y.  This keeps bruk2pipe's
        digital-filter correction while removing the pseudo dimension.
        """
        if self.tp != 'bruk' or self.dim != '3p':
            raise ValueError('Bruker pseudo-Y slice extraction requires a Bruker 3p dataset')
        if getattr(self, 'pseudo_acq_axis', None) != 'y':
            raise ValueError('Bruker pseudo-Y slice extraction called for a non-Y pseudo axis')

        import math
        import re
        import shlex

        ser = os.path.join(self.FidPath, 'ser')
        if not os.path.isfile(ser):
            raise FileNotFoundError('Bruker SER file not found: %s' % ser)
        os.makedirs(self.outdir, exist_ok=True)
        raw_slice = os.path.join(self.outdir, 'slice.ser')
        dst = os.path.join(self.outdir, 'slice.fid')
        for path in (raw_slice, dst):
            if os.path.exists(path):
                os.remove(path)

        # Bruker pads each direct-dimension record to a 1024-byte boundary.
        # Preserve those records exactly: bruk2pipe, not Python, must perform
        # the vendor-specific decoding and digital-filter correction.
        td_words = int(float(self.GetParBruk('acqus', ('', 'TD'))[0]))
        dtypa = int(getattr(self, 'DTYPA', int(self.GetParBruk('acqus', ('', 'DTYPA'))[0])))
        word_bytes = 8 if dtypa == 2 else 4
        record_bytes = int(math.ceil(float(td_words * word_bytes) / 1024.0) * 1024)
        total_bytes = os.path.getsize(ser)
        if record_bytes <= 0 or total_bytes % record_bytes:
            raise ValueError('Invalid Bruker SER record geometry: size=%d record=%d TD=%d'
                             % (total_bytes, record_bytes, td_words))
        nrecords = total_bytes // record_bytes
        nz = int(self.nz)
        second = nz  # zero based => one-based record NZ+1
        if nz < 1 or second >= nrecords:
            raise ValueError('Cannot extract Bruker records 1 and %d from %d SER records'
                             % (nz + 1, nrecords))

        with open(ser, 'rb') as src, open(raw_slice, 'wb') as out:
            first_block = src.read(record_bytes)
            if len(first_block) != record_bytes:
                raise IOError('Short read while extracting Bruker record 1')
            out.write(first_block)
            src.seek(second * record_bytes)
            second_block = src.read(record_bytes)
            if len(second_block) != record_bytes:
                raise IOError('Short read while extracting Bruker record %d' % (nz + 1))
            out.write(second_block)

        # Do not reconstruct Bruker conversion/header parameters here.  The
        # already-successful main conversion script has resolved AQSEQ and is
        # the single source of truth.  Before its TP/ZTP/TP rearrangement the
        # pseudo axis is physical Y and the spectral N15 axis is physical Z,
        # so copy X verbatim and rename the main Z fields to Y.
        main_script = os.path.join(self.outdir, 'fid.test.com')
        if not os.path.isfile(main_script):
            raise FileNotFoundError('Main conversion script not found: %s' % main_script)
        text = open(main_script, 'r', encoding='utf-8').read()
        # Join genuine csh continuation lines before tokenising.  Also tolerate
        # older generated scripts containing extra whitespace after '\\'.
        text = re.sub(r'\\[ \t]*\r?\n', ' ', text)
        bruk_line = None
        for line in text.splitlines():
            if 'bruk2pipe ' in line:
                bruk_line = line[line.index('bruk2pipe '):]
                break
        if bruk_line is None:
            raise ValueError('Could not find bruk2pipe command in %s' % main_script)
        # Ignore downstream nmrPipe/pipe2xyz stages; only bruk2pipe arguments
        # describe the raw vendor data.
        bruk_line = bruk_line.split('|', 1)[0].strip()
        tokens = shlex.split(bruk_line)

        def option_value(name):
            try:
                i = tokens.index(name)
            except ValueError:
                raise ValueError('Main conversion is missing required option %s' % name)
            if i + 1 >= len(tokens):
                raise ValueError('Main conversion has no value for option %s' % name)
            return tokens[i + 1]

        # General bruk2pipe options which control byte interpretation and the
        # Bruker digital filter.  Preserve their values exactly as Guess wrote
        # them.  Boolean switches are copied separately.
        general_pairs = []
        for opt in ('-bad', '-decim', '-dspfvs', '-grpdly', '-ws'):
            if opt in tokens:
                general_pairs.extend((opt, option_value(opt)))
        general_switches = [opt for opt in ('-aswap', '-noaswap', '-AMX', '-DMX', '-noi2f')
                            if opt in tokens]

        x_keys = ('N', 'T', 'MODE', 'SW', 'OBS', 'CAR', 'LAB')
        xvals = {key: option_value('-x' + key) for key in x_keys}
        zvals = {key: option_value('-z' + key) for key in x_keys}

        # The artificial indirect dimension consists of two raw spectral
        # increments.  Its quadrature/header metadata otherwise comes directly
        # from physical Z of the correct full conversion (e.g. N15), not from
        # the physical-Y pseudo axis (e.g. ncyc).
        zvals['N'] = '2'
        zvals['T'] = '1'

        script = os.path.join(self.outdir, 'fid.test.slice.com')
        lines = []
        lines.append('bruk2pipe -in %s' % raw_slice)
        general = general_switches + general_pairs
        if general:
            lines.append('  ' + ' '.join(general))
        lines.append('  ' + ' '.join('-%s%s %s' % ('x', key, xvals[key]) for key in x_keys))
        lines.append('  ' + ' '.join('-%s%s %s' % ('y', key, zvals[key]) for key in x_keys))
        lines.append('  -ndim 2 -aq2D States -out $ft4trec -verb -ov')

        with open(script, 'w', encoding='utf-8') as f:
            f.write('#!/bin/csh\n\n')
            f.write('set ft4trec=%s/slice.fid\n' % self.outdir)
            f.write('if( -e $ft4trec) rm -f $ft4trec\n\n')
            # A csh continuation backslash MUST be the final character before
            # the newline.  Writing each physical line separately prevents the
            # previous "\\ -bad" malformed command and argument shifting.
            for line in lines[:-1]:
                f.write(line + ' \\\n')
            f.write(lines[-1] + '\n')
        return script

    def _pseudo2d_phase_slice_mode(self):
        """Return the Process-window pseudo-2D phasing trace policy."""
        parent = getattr(self, 'parent', None)
        value = getattr(parent, 'phaseSliceModeValue', None)
        if value is None and parent is not None and hasattr(parent, 'phaseSliceMode'):
            try:
                value = parent.phaseSliceMode.GetStringSelection()
            except Exception:
                value = None
        if value is None and parent is not None:
            state = getattr(parent, 'state', None)
            if state is not None:
                value = getattr(state, 'metadata', {}).get('phaseSliceMode')
        value = str(value or 'First').strip().lower()
        return 'Summed' if value == 'summed' else 'First'

    def _prepare_pseudo2d_slice_input(self):
        """Create a 1D phasing FID from a converted pseudo-2D stack.

        ``test.fid`` is the normal converted 2D time-domain stack (pseudo x
        direct).  First selects row zero; Summed adds every pseudo row in the
        time domain.  The resulting NMRPipe file is deliberately flattened to
        one dimension before the existing 1D preview processing is run.
        """
        import gzip
        import tempfile
        import numpy as np
        import nmrglue as ng

        src = os.path.join(self.outdir, 'test.fid')
        src_gz = src + '.gz'
        dst = os.path.join(self.outdir, 'slice.fid')
        source = src
        temporary_source = None
        if not os.path.exists(source):
            if not os.path.exists(src_gz):
                raise FileNotFoundError('Could not find pseudo-2D conversion output %r or %r' % (src, src_gz))
            fd, temporary_source = tempfile.mkstemp(suffix='.fid')
            os.close(fd)
            with gzip.open(src_gz, 'rb') as fin, open(temporary_source, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
            source = temporary_source
        try:
            dic, data = ng.pipe.read(source)
            data = np.asarray(data)
            if data.ndim < 2:
                trace = np.asarray(data)
            elif self._pseudo2d_phase_slice_mode() == 'Summed':
                trace = np.sum(data, axis=tuple(range(data.ndim - 1)), dtype=data.dtype)
            else:
                trace = data.reshape((-1, data.shape[-1]))[0].copy()
            trace = np.asarray(trace)
            dic = dict(dic)
            dic['FDDIMCOUNT'] = 1.0
            dic['FDSIZE'] = float(trace.shape[-1])
            dic['FDSPECNUM'] = 1.0
            dic['FDFILECOUNT'] = 1.0
            dic['FDPIPEFLAG'] = 0.0
            if 'FDPIPECOUNT' in dic:
                dic['FDPIPECOUNT'] = 0.0
            ng.pipe.write(dst, dic, trace, overwrite=True)
        finally:
            if temporary_source is not None and os.path.exists(temporary_source):
                os.remove(temporary_source)
        return dst

    def _prepare_pseudo3d_slice_input(self):
        """Write the first converted pseudo plane as a standalone 2D file.

        This remains the normal 3p path when the pseudo dimension is already
        the outer physical axis.  Bruker pseudo-Y acquisitions are handled by
        :meth:`_prepare_bruker_pseudo_y_slice_input` instead.
        """
        src = os.path.join(self.outdir, 'fids', 'test001.fid')
        src_gz = src + '.gz'
        dst = os.path.join(self.outdir, 'slice.fid')

        if os.path.exists(dst):
            try:
                os.remove(dst)
            except Exception:
                pass

        source = src
        temporary_source = None
        if not os.path.exists(source):
            if not os.path.exists(src_gz):
                raise FileNotFoundError(
                    f'Could not find first pseudo-axis conversion plane {src!r} or {src_gz!r}'
                )
            import gzip
            import tempfile
            fd, temporary_source = tempfile.mkstemp(suffix='.fid')
            os.close(fd)
            with gzip.open(src_gz, 'rb') as fin, open(temporary_source, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
            source = temporary_source

        try:
            import nmrglue as ng
            if hasattr(ng.pipe, 'read_2D'):
                dic, data = ng.pipe.read_2D(source)
            else:
                from nmrglue.fileio import pipe as ng_pipe
                dic, data = ng_pipe.read_2D(source)
            dic['FDDIMCOUNT'] = 2.0
            dic['FDFILECOUNT'] = 1.0
            dic['FDPIPEFLAG'] = 1.0
            if 'FDPIPECOUNT' in dic:
                dic['FDPIPECOUNT'] = 0.0
            ng.pipe.write(dst, dic, data, overwrite=True)
        finally:
            if temporary_source is not None and os.path.exists(temporary_source):
                try:
                    os.remove(temporary_source)
                except Exception:
                    pass
        return dst

    def _build_slice_processing_script(self, parent):
        script_builder = None
        if parent is not None and hasattr(parent, 'nmrPipe'):
            script_builder = parent.nmrPipe
        if script_builder is None:
            from spinDecon.processing.nmrpipe_scripts import nmrPipe as _nmrPipe
            script_builder = _nmrPipe(parent)
            if parent is not None:
                parent.nmrPipe = script_builder
        script_path = self.outdir + '/nmrproc.1D.com'
        script_builder.make_proc_script_1d_slice(parent, script_path)
        return script_path

    def ProcessSlice(self, output_frame=None, on_finish=None):
        """Prepare and run preview/slice processing without blocking wx.

        When called by ConversionFrame, all subprocess output is appended to the
        same Conversion Output window used by the primary conversion command.
        """
        from spinDecon.gui.dialogs.shell_output import run_command_with_output

        parent = getattr(self, 'parent', None)

        def done(*_args):
            if on_finish is not None:
                on_finish()

        def run_processing(*_args):
            if output_frame is not None:
                output_frame.start_step('Process preview spectrum')
            try:
                script_path = self._build_slice_processing_script(parent)
            except Exception as exc:
                pass
                if output_frame is not None:
                    output_frame.append_text('\nCould not prepare 1D preview processing: %s\n' % exc)
                    output_frame.set_status('Failed')
                done()
                return
            run_command_with_output(
                ['csh', script_path], parent=parent, title='Conversion Output',
                output_frame=output_frame, on_finish=done, final=False,
                label='Process 1D preview slice')

        if self._slice_dim_count() <= 0:
            done()
            return

        try:
            if self._slice_dim_count() == 1:
                self._prepare_1d_slice_input()
                run_processing()
                return

            # Varian 3p data must build the phasing preview directly from the
            # vendor acquisition.  PipeParseSlice() deliberately converts only
            # X plus the first complex Y point (two raw traces), omitting the
            # real pseudo axis, and writes the result as spec/slice.fid.
            # Keep this ahead of the generic 3p handling below: that path uses
            # an already-converted pseudo plane and is appropriate for Bruker,
            # but it bypasses the Varian slice extraction script.
            if self.tp == 'var' and self.dim == '3p':
                if output_frame is not None:
                    output_frame.start_step('Extract Varian phasing traces')
                self.PipeParseSlice()
                slice_script = self.outdir + '/fid.test.slice.com'
                run_command_with_output(
                    ['csh', slice_script], parent=parent, title='Conversion Output',
                    output_frame=output_frame, on_finish=run_processing, final=False,
                    label='Extract Varian phasing traces')
                return

            # Two spectral dimensions plus a pseudo axis (3p): the pseudo axis
            # is real, not spectral.  The main conversion has already emitted
            # one 2D file per pseudo point.  For preview/phasing use only the
            # first converted plane; do not build a second indirect projection
            # and, importantly, do not ask for non-existent sw2/frq2 metadata.
            # Varian 3p is handled above from the vendor binary.  This branch
            # retains the Bruker pseudo-plane handling, including pseudo-in-Y.
            if self.dim == '3p':
                # Bruker can physically acquire the real pseudo dimension in Y.
                # In that layout adjacent SER traces are different pseudo
                # values, so the phasing pair is trace 1 and trace NZ+1 rather
                # than traces 1 and 2.  Extract that pair directly from the
                # vendor SER with nmrglue and write a true 2D slice.fid.
                if self.tp == 'bruk' and getattr(self, 'pseudo_acq_axis', None) == 'y':
                    if output_frame is not None:
                        output_frame.start_step('Extract spectral phasing traces')
                    slice_script = self._prepare_bruker_pseudo_y_slice_input()
                    if output_frame is not None:
                        output_frame.start_step('Convert selected Bruker traces')
                    def _slice_converted(*_args):
                        if not os.path.isfile(os.path.join(self.outdir, 'slice.fid')):
                            if output_frame is not None:
                                output_frame.append_text('\nbruk2pipe did not create slice.fid\n')
                                output_frame.set_status('Failed')
                            done()
                            return
                        run_processing()
                    run_command_with_output(
                        ['csh', slice_script], parent=parent, title='Conversion Output',
                        output_frame=output_frame, on_finish=_slice_converted, final=False,
                        label='Convert selected Bruker phasing traces')
                    return

                # The main conversion stores the pseudo dimension as the outer
                # file axis, so test001.fid is exactly the first XY spectral
                # plane we want for phasing.  A plane file copied out of an
                # NMRPipe 3D file family still carries the parent 3D header,
                # however.  If passed straight to xyz2pipe it therefore tries
                # to find all pseudo planes again.  Copy the first plane and
                # explicitly flatten only its header to a standalone 2D
                # NMRPipe file before invoking the existing 1D phasing path.
                if output_frame is not None:
                    output_frame.start_step('Prepare first pseudo-axis plane')
                self._prepare_pseudo3d_slice_input()
                run_processing()
                return

            # One spectral dimension plus a real pseudo axis: phasing can
            # use either the first converted time-domain row or the coherent
            # time-domain sum of every row in the pseudo stack.
            if self.dim == '2p':
                if output_frame is not None:
                    output_frame.start_step('Prepare pseudo-2D phasing FID (%s)' % self._pseudo2d_phase_slice_mode())
                self._prepare_pseudo2d_slice_input()
                run_processing()
                return

            # Ordinary multidimensional spectral data use the existing slice
            # conversion logic.
            if output_frame is not None:
                output_frame.start_step('Prepare preview slice')
            self.PipeParseSlice()
            slice_script = self.outdir + '/fid.test.slice.com'
            run_command_with_output(
                ['csh', slice_script], parent=parent, title='Conversion Output',
                output_frame=output_frame, on_finish=run_processing, final=False,
                label='Extract conversion preview slice')
        except Exception as exc:
            pass
            if output_frame is not None:
                output_frame.append_text('\nCould not prepare preview slice: %s\n' % exc)
                output_frame.set_status('Failed')
            done()

    #work on the conversion values.
    def Convert(self):
        if(self.tp=='var'):
            self.ConvertVarian()
        elif(self.tp=='bruk'):
            self.ConvertBruker()
        elif(self.tp=='omeg'):
            self.ConvertOmega()
        # Publish a vendor-independent xyza view only after the selected vendor
        # path has finished.  This does not alter np/ni/ni2/ni3 or any script
        # calculation; it only exposes the already-resolved conversion sizes.
        self._publish_time_domain_sizes()
        self._publish_observation_frequencies()

    # Based on Patrik Lundstrom 011126
    #take water and sfrq, calc ppms of C and N
    def shift(self,dfrq,nuc='C13'):
        if(nuc=='H1'):
            CONV=1.
        if(nuc=='C13'):
            CONV=0.251449530
        if(nuc=='N15'):
            CONV=0.101329118
        if(nuc=='P31'):
            CONV=0.4048064954
        if(nuc=='F19'):
            CONV=0.9412866605363297

        sfrq0  = self.sfrq / (1.0 + self.waterppmTOF*1e-6);
        dfrq0 = sfrq0*CONV;
        ppm = (dfrq-dfrq0)/dfrq0*1e6;
        return ppm


    #figure out which nuclei are begin analysed.
    def GetNuc(self,lab):
        """
        if(nuc=='H1'):
            CONV=1.
        if(nuc=='C13'):
            CONV=0.251449530
        if(nuc=='N15'):
            CONV=0.101329118
        if(nuc=='P31'):
            CONV=0.4048064954
        if(nuc=='F19'):
            CONV=0.9412866605363297
        """
        pass

        if(lab[0]=='C'):
            if(self.tn=='H1'):
                if(self.dim==2 and self.ni>1):
                    if(numpy.fabs(self.dfrq/self.sfrq-0.25)<2E-2):
                        return self.dfrq,self.dn
                elif(self.dim==2 and self.ni==1):
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.25)<2E-2):
                        return self.dfrq2,self.dn2
                else:
                    if(numpy.fabs(self.dfrq/self.sfrq-0.25)<2E-2):
                        return self.dfrq,self.dn
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.25)<2E-2):
                        return self.dfrq2,self.dn2
            if(self.tn=='F19'):
                if(self.dim==2 and self.ni>1):
                    if(numpy.fabs(self.dfrq/self.sfrq-0.26)<1E-1):
                        return self.dfrq,self.dn
                elif(self.dim==2 and self.ni==1):
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.26)<1E-1):
                        return self.dfrq2,self.dn2
                else:
                    if(numpy.fabs(self.dfrq/self.sfrq-0.26)<1E-1):
                        return self.dfrq,self.dn
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.26)<1E-1):
                        return self.dfrq2,self.dn2
                
        if(lab[0]=='N'):
            if(self.tn=='H1'):
                if(self.dim==2 and self.ni>1):
                    if(numpy.fabs(self.dfrq/self.sfrq-0.1)<2E-2):
                        
                        return self.dfrq, self.dn
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.1)<2E-2):
                        return self.dfrq2,self.dn2
                elif(self.dim==2 and self.ni==1):
                    
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.1)<2E-2):
                        return self.dfrq2,self.dn2
                else:
                    if(numpy.fabs(self.dfrq/self.sfrq-0.1)<2E-2):
                        return self.dfrq,self.dn
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.1)<2E-2):
                        return self.dfrq2,self.dn2
            if(self.tn=='F19'):
                if(self.dim==2 and self.ni>1):
                    if(numpy.fabs(self.dfrq/self.sfrq-0.26)<1E-2):
                        return self.dfrq,self.dn
                elif(self.dim==2 and self.ni==1):
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.26)<1E-2):
                        return self.dfrq2,self.dn2
                else:
                    if(numpy.fabs(self.dfrq/self.sfrq-0.26)<1E-2):
                        return self.dfrq,self.dn
                    if(numpy.fabs(self.dfrq2/self.sfrq-0.26)<1E-2):
                        return self.dfrq2,self.dn2
        if(lab[0]=='H'):
            if(numpy.fabs(self.sfrq/self.sfrq-1)<1E-2):
                return self.sfrq,self.tn
        if(lab[0]=='F'):
            if(numpy.fabs(self.sfrq/self.sfrq-1)<1E-2):
                return self.sfrq,self.tn
        pass
        self.abort=1
        return -1

    #calculte fid.test.com pars for bruker
    def ConvertBruker(self):

        self.ns=self.GetParBruk('acqus',('','NS'))[0]
        # Acquisition metadata used by the project summary.  Bruker D1 is
        # the first recycle-delay element in the D array; GetParBruk already
        # supports indexed array parameters such as D1.
        try:
            self.d1 = float(self.GetParBruk('acqus', ('', 'D1'))[0])
        except Exception:
            self.d1 = None
        #N=self.np,self.ni*2,self.ni2*2,self.ni3*2
        #T=self.np/2,self.ni,self.ni2,self.ni3



        # The user-selected dimensionality is authoritative; use TD from the
        # acquisition parameters rather than trying to infer sizes from a data read.
        self.np = int(self.GetParBruk('acqus',('','TD'))[0])
        if(self.np <= 0):
            self.np = 1
        self.np2 = int(self.np/2)
        self.xN = self.np2

        self.BYTORDA = int(self.GetParBruk('acqus',('','BYTORDA'))[0])
        self.bytor = 'noaswap' if self.BYTORDA == 1 else 'aswap'
        self.DTYPA = int(self.GetParBruk('acqus',('','DTYPA'))[0])
        self.DECIM = float(self.GetParBruk('acqus',('','DECIM'))[0])
        self.DSPFVS = float(self.GetParBruk('acqus',('','DSPFVS'))[0])
        self.GRPDLY = float(self.GetParBruk('acqus',('','GRPDLY'))[0])
        self.dtext = ''
        if(self.DTYPA==2):
            self.dtext='-ws 8 -noi2f'

        if(self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
            if(os.path.exists(os.path.join(self.FidPath,'acqu2s'))==0):
                pass
                aq2=False
            else:
                aq2=True
        else:
            aq2=True


        if(self.dim=='3p' or self._dim_at_least(3)):
            if(os.path.exists(os.path.join(self.FidPath,'acqu3s'))==0):
                pass
                aq3=False
            else:
                aq3=True
        else:
            aq3=True

        if(self.nuslist==''):
            if(aq3): #get from acquX or L345
                if(self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
                    self.ni=int(self.GetParBruk('acqu2s',('','TD'))[0])/2
                if(self.dim=='3p' or self._dim_at_least(3)):
                    self.ni2=int(self.GetParBruk('acqu3s',('','TD'))[0])/2
                if(type(self.dim)!=str and self._dim_at_least(4)):
                    self.ni3=int(self.GetParBruk('acqu4s',('','TD'))[0])/2
            else: #to handle Marius Clores' 4D
                self.ni=int(self.GetParBruk('acqus',('','L3'))[0])
                if(self._dim_at_least(2)):
                    self.ni2=int(self.GetParBruk('acqus',('','L4'))[0])
                if(type(self.dim)!=str and self.dim==4):
                    self.ni3=int(self.GetParBruk('acqus',('','L5'))[0])
                if(type(self.dim)!=str and self.dim==4):
                    self.ni3=int(self.GetParBruk('acqu4s',('','TD'))[0])/2
        else:
            if(self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
                self.ni=int(self.GetParBruk('acqu2s',('','NusTD'))[0])/2
            if(self.dim=='3p' or self._dim_at_least(3)):
                self.ni2=int(self.GetParBruk('acqu3s',('','NusTD'))[0])/2
            if(self._dim_at_least(4)):
                self.ni3=int(self.GetParBruk('acqu4s',('','NusTD'))[0])/2


        if type(self.dim) == str and self.dim in ('2p', '3p'):
            self.pseudo_axis_info = detect_bruker_pseudo_axis(self.FidPath)
            self.nz = self.pseudo_axis_info.size
            self.pseudo_acq_axis = None
            self.pseudo_logical_axis = None
            if self.dim == '3p':
                # A Bruker pseudo-3D can acquire the real pseudo dimension as
                # either physical Y (acqu2s) or physical Z (acqu3s).
                td2 = int(float(self.GetParBruk('acqu2s', ('', 'TD'))[0])) if aq2 else 0
                td3 = int(float(self.GetParBruk('acqu3s', ('', 'TD'))[0])) if aq3 else 0
                y_match = (td2 == self.nz)
                z_match = (td3 == self.nz)
                if y_match and not z_match:
                    self.pseudo_logical_axis = 'y'
                elif z_match:
                    self.pseudo_logical_axis = 'z'
                elif y_match:
                    self.pseudo_logical_axis = 'y'
                else:
                    raise ValueError(
                        'Bruker pseudo-axis size %d does not match acqu2s TD=%d or acqu3s TD=%d'
                        % (self.nz, td2, td3))

        if(self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
            self.yMode=self.GetParBruk('acqu2s',('','FnMODE'))[0]
        if(self.dim=='3p' or self._dim_at_least(3)):
            if(aq3):
                self.zMode=self.GetParBruk('acqu3s',('','FnMODE'))[0]
                if(self.dim==4):
                    self.aMode=self.GetParBruk('acqu4s',('','FnMODE'))[0]
            else: #deal with Marius Clore's 4D
                self.zMode=str(0)
                if(self.dim==4):
                    self.aMode=str(0)

        self.aqseq=self.GetAcqseq()

        # For pseudo-3D Bruker data, acqu2s/acqu3s identify the TopSpin
        # logical F2/F1 dimensions, but AQSEQ determines their physical order
        # in the SER stream.  aqseq 312 exchanges the two indirect dimensions
        # for bruk2pipe: a pseudo list attached to acqu3s (F1) is therefore
        # physical Y, while one attached to acqu2s (F2) is physical Z.
        if self.dim == '3p' and self.pseudo_logical_axis in ('y', 'z'):
            # Map TopSpin logical acqu2s/acqu3s axes to their physical order
            # in the SER stream exactly once.  Keep pseudo_logical_axis intact:
            # it is needed later to select the correct acquisition metadata.
            self.pseudo_acq_axis = self.pseudo_logical_axis
            if self.aqseq == '312':
                self.pseudo_acq_axis = 'y' if self.pseudo_logical_axis == 'z' else 'z'
        pass
        if(type(self.dim)!=str and self.dim==3):
            if(self.aqseq=='312'):
                pass
                self.ni2,self.ni=self.ni,self.ni2
                self.yMode,self.zMode=self.zMode,self.yMode
            else: #self.aqseq='321'
                pass

        self.mode=[]
        if (self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
            self.mode.append(self.yMode)
        if(self.dim=='3p' or self._dim_at_least(3)):
            self.mode.append(self.zMode)
        if(type(self.dim)!=str and self._dim_at_least(4)):
            self.mode.append(self.aMode)


        #self.ni=int(self.GetParBruk('acqus',('','L3'))[0])
        #self.ni2=int(self.GetParBruk('acqus',('','L4'))[0])
        #if(self.dim==4):
        #    self.ni3=int(self.GetParBruk('acqus',('','L5'))[0])

        #sw=self.sw,self.sw1,self.sw2,self.sw3
        #try:
        #    self.sw=float(self.GetParBruk('acqus',('','SW_h'))[0])
        #    self.sw1=1/(2*float(self.GetParBruk('acqus',('','IN0'))[0]))
        #    self.sw2=1/(2*float(self.GetParBruk('acqus',('','IN8'))[0]))
        #    if(self.dim==4):
        #        self.sw3=1/(2*float(self.GetParBruk('acqus',('','IN19'))[0]))
        #except:
        self.sw=float(self.GetParBruk('acqus',('','SW_h'))[0])
        self.SFO1=float(self.GetParBruk('acqus',('','SFO1'))[0])  #typically H
        self.BF1=float(self.GetParBruk('acqus',('','BF1'))[0])  #typically H
        self.O1=float(self.GetParBruk('acqus',('','O1'))[0])  #typically H
        NUC1=self.GetParBruk('acqus',('','NUC1'))[0].replace('<', '').replace('>', '').strip()
        #print(f'[vpar] acqus NUC1 raw={NUC1!r}', flush=True)
        self.NUC1_raw = NUC1
        self.NUC1=NUC1[-1:]+NUC1[:-1]
        self.n1 = self.NUC1
        self.frq=self.SFO1

        if(self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
            if(aq2):
                self.sw1=(float(self.GetParBruk('acqu2s',('','SW_h'))[0]))
                

                self.SFO2=float(self.GetParBruk('acqu2s',('','SFO1'))[0])  #typically C
                self.BF2=float(self.GetParBruk('acqu2s',('','BF1'))[0])  #typically C
                self.O2=float(self.GetParBruk('acqu2s',('','O1'))[0])  #typically C

                if self.sw1 == 0.:
                    self.sw1 = (float(self.GetParBruk('acqu2s',('','SW'))[0]))*self.SFO2 ## Old bruker par files don't store SW_h
                NUC2=self.GetParBruk('acqu2s',('','NUC1'))[0].replace('<', '').replace('>', '').strip()
                #print(f'[vpar] acqu2s NUC1 raw={NUC2!r}', flush=True)
                self.NUC2_raw = NUC2
                self.NUC2=NUC2[-1:]+NUC2[:-1]
                self.n2 = self.NUC2
                self.frq1=self.SFO2
            else:
                pass
                return
        if(self.dim=='3p' or self._dim_at_least(3)):
            if(aq3):
                self.sw2=(float(self.GetParBruk('acqu3s',('','SW_h'))[0]))
                self.SFO3=float(self.GetParBruk('acqu3s',('','SFO1'))[0]) #typically N
                self.BF3=float(self.GetParBruk('acqu3s',('','BF1'))[0]) #typically N
                self.O3=float(self.GetParBruk('acqu3s',('','O1'))[0]) #typically N
                if self.sw2 == 0.:
                    self.sw2 = (float(self.GetParBruk('acqu3s',('','SW'))[0]))*self.SFO3  ## Old bruker par files don't store SW_h
                NUC3=self.GetParBruk('acqu3s',('','NUC1'))[0].replace('<', '').replace('>', '').strip()
                #print(f'[vpar] acqu3s NUC1 raw={NUC3!r}', flush=True)
                self.NUC3_raw = NUC3
                self.NUC3=NUC3[-1:]+NUC3[:-1]
                self.n3 = self.NUC3
                self.frq2=self.SFO3
                if(self.dim==4):
                    self.sw3=(float(self.GetParBruk('acqu4s',('','SW_h'))[0]))
                    self.SFO4=float(self.GetParBruk('acqu4s',('','SFO1'))[0]) #typically N
                    self.BF4=float(self.GetParBruk('acqu4s',('','BF1'))[0]) #typically N
                    self.O4=float(self.GetParBruk('acqu4s',('','O1'))[0]) #typically N
                    NUC4=self.GetParBruk('acqu4s',('','NUC1'))[0].replace('<', '').replace('>', '').strip()
                    #print(f'[vpar] acqu4s NUC1 raw={NUC4!r}', flush=True)
                    self.NUC4_raw = NUC4
                    self.NUC4=NUC4[-1:]+NUC4[:-1]
                    self.n4 = self.NUC4
                    self.frq3=self.SFO4
            else: #deal with Marius Clore's 4D
                self.sw1=1/(2*float (self.GetParBruk('acqus',('','IN0'))[0]))
                self.frq1,self.f1ppm=self.GetShiftBruk(self.labb[1])
                if(self._dim_at_least(3)):
                    self.sw2=1/(2*float(self.GetParBruk('acqus',('','IN8'))[0]))
                    self.frq2,self.f2ppm=self.GetShiftBruk(self.labb[2])
                    if(self.labb[2][0]==self.labb[0][0]): #if the labelled nucleus for Z and X are the same...
                        self.SFO3=float(self.GetParBruk('acqus',('','SFO1'))[0]) #typically N
                        self.BF3=float(self.GetParBruk('acqus',('','BF1'))[0]) #typically N
                        self.O3=float(self.GetParBruk('acqus',('','O1'))[0]) #typically N
                    elif(self.labb[2][0]==self.labb[1][0]):  #if the labelled nucleus for Z and Y are the same...
                        self.SFO3=float(self.GetParBruk('acqu2s',('','SFO1'))[0]) #typically N
                        self.BF3=float(self.GetParBruk('acqu2s',('','BF1'))[0]) #typically N
                        self.O3=float(self.GetParBruk('acqu2s',('','O1'))[0]) #typically N

                    if(self.labb[2][0]=='C'):
                        self.NUC3='C13'
                    elif(self.labb[2][0]=='H'):
                        self.NUC3='H1'
                    elif(self.labb[2][0]=='N'):
                        self.NUC3='N15'
                    if(self.dim==4):
                        self.sw3=1/(2*float(self.GetParBruk('acqus',('','IN19'))[0]))
                        self.frq3,self.f3ppm=self.GetShiftBruk(self.labb[3])

                        if(self.labb[3][0]==self.labb[0][0]):
                            self.SFO4=float(self.GetParBruk('acqus',('','SFO1'))[0]) #typically N
                            self.BF4=float(self.GetParBruk('acqus',('','BF1'))[0]) #typically N
                            self.O4=float(self.GetParBruk('acqus',('','O1'))[0]) #typically N
                        elif(self.labb[3][0]==self.labb[1][0]):
                            self.SFO4=float(self.GetParBruk('acqu2s',('','SFO1'))[0]) #typically N
                            self.BF4=float(self.GetParBruk('acqu2s',('','BF1'))[0]) #typically N
                            self.O4=float(self.GetParBruk('acqu2s',('','O1'))[0]) #typically N

                        if(self.labb[3][0]=='C'):
                            self.NUC4='C13'
                        elif(self.labb[3][0]=='H'):
                            self.NUC4='H1'
                        elif(self.labb[3][0]=='N'):
                            self.NUC4='N15'

        if(type(self.dim)!=str and self.dim==3):
            pass
        #print self.SFO1,self.SFO2,self.SFO3
        #print self.O1,self.O2,self.O3

        #O=self.sfrq,self.frq1,self.frq2,self.frq3
        #C=self.waterppm,self.f1ppm,self.f2ppm,self.f3ppm
        # The direct dimension observation frequency must come from the
        # Bruker base frequency in acqus (BF1). This is the value expected
        # by nmrPipe as the direct-dimension OBS/xOBS setting.
        self.sfrq = float(self.BF1)
        self.waterppm = self.O1 / self.BF1
        self.waterppmTOF = self.waterppm  ##CHARLIE ADDED NOT SURE
        #self.frq1,self.f1ppm=self.GetShiftBruk(self.labb[1])
        #self.frq2,self.f2ppm=self.GetShiftBruk(self.labb[2])
        #if(self.dim==4):
        #    self.frq3,self.f3ppm=self.GetShiftBruk(self.labb[3])



        if(self.o1p=='Water'): #re-reference assuming carrier is water
            #print self.NUC1,self.NUC2,self.NUC3
            self.temp=float(self.GetParBruk('acqus',('','TE'))[0])
            try:
                self.waterppm=self.WaterPPM()
                pass
            except:
                pass
                self.waterppm=4.7
            self.waterppmTOF=self.waterppm

            if (self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
                self.f1ppm=self.shift(self.frq1,nuc=self.NUC2)
            if(self.dim=='3p' or self._dim_at_least(3)):
                self.f2ppm=self.shift(self.frq2,nuc=self.NUC3)
                pass
                if(self.dim==4):
                    self.f3ppm=self.shift(self.frq3,nuc=self.NUC4)

        elif(type(self.o1p)!=str): #if not water or auto, then it's a ppm value set manually by the user
            #reference against the value maunally set
            #print self.NUC1,self.NUC2,self.NUC3
            self.waterppm=self.o1p
            self.waterppmTOF=self.waterppm #self.shift will take THIS number
            if self._dim_at_least(2):
                self.f1ppm=self.shift(self.frq1,nuc=self.NUC2)
            if(self._dim_at_least(3)):
                self.f2ppm=self.shift(self.frq2,nuc=self.NUC3)
                pass
                if(self.dim==4):
                    self.f3ppm=self.shift(self.frq3,nuc=self.NUC4)



        #if(self.o1p=='Auto'):
        else:
            if self.dim in ('2p', '3p') or (not isinstance(self.dim, str) and self._dim_at_least(2)):
                self.f1ppm=self.O2/self.BF2
            if not isinstance(self.dim, str) and self._dim_at_least(3):
                self.f2ppm=self.O3/self.BF3
                if(type(self.dim)!=str and self.dim==4): #flemming's data implies we can't trust o1 in acqu4s
                    if(self.NUC4==self.NUC1):
                        self.f3ppm=self.waterppm
                        self.frq3=self.sfrq
                    elif(self.NUC4==self.NUC2):
                        self.f3ppm=self.f1ppm
                        self.frq3=self.frq1
                    elif(self.NUC4==self.NUC3):
                        self.f3ppm=self.f2ppm
                        self.frq3=self.frq2
                    else:
                        self.f3ppm=self.O4/self.BF4

        
    

            #print self.O1,self.BF1,self.SFO1
            #print self.O4,self.BF4,self.SFO4

            #print self.waterppm
            #sys.exit(100)

 

        

        """
        #if(self.tn!='H1'):
        #    print 'direct nucleus isnt proton: you will need to write your own conversion script.'
        #    self.abort=1
        #    return
        print self.labb,self.labb[0]
        self.frq1,self.n1=self.GetNuc(self.labb[1])
        self.frq2,self.n2=self.GetNuc(self.labb[2])


        if(self.frq1==-1 or self.frq2==-1):
            return

        if(self.dim==4):
            self.frq3,self.n3=self.GetNuc(self.labb[3])
            if(self.frq3==-1):
                return

        self.f1ppm=self.shift(self.frq1,self.n1)
        self.f2ppm=self.shift(self.frq2,self.n2)
        if(self.dim==4):
            self.f3ppm=self.shift(self.frq3,self.n3)

        #npAdj=BrukFidAdjust(np)
        """

    #calculate reference frequencies, bruker.
    def GetShiftBruk(self,lab):
        pass
        lab0 = lab[0] if lab else ''
        nuc1 = getattr(self, 'NUC1', None)
        nuc2 = getattr(self, 'NUC2', None)
        nuc3 = getattr(self, 'NUC3', None)
        if(nuc1 and lab0==nuc1[0]):
            return self.BF1,self.O1/self.BF1
        if(nuc2 and lab0==nuc2[0]):
            return self.SFO2,self.O2/self.BF2
        if(nuc3 and lab0==nuc3[0]):
            return self.SFO3,self.O3/self.BF3
        pass
        return -1,-1






    #calculate fid.test.com pars for varian
    def ConvertVarian(self):
        # Keep acquisition metadata on vpar so callers do not need to parse
        # procpar independently.
        try:
            self.ns = self.GetParVarian(('', 'nt'))[0]
        except Exception:
            self.ns = None
        try:
            self.d1 = float(self.GetParVarian(('', 'd1'))[0])
        except Exception:
            self.d1 = None
        self.ni=int(self.GetParVarian(('','ni'))[0])
        self.np=int(self.GetParVarian(('','np'))[0])
        self.np2=int(self.GetParVarian(('','np'))[0])/2
        self.xN = self.np2
        
        pass
        array=self.GetParVarian(('','array'))[0].split('"')[1].split(',')

        

        #nz=len(vpar.GetParVarian('./procpar','n',('','ncyc_cp')))
        if(type(self.dim)!=str):
            if(self.dim==3):
                self.ni2=int(self.GetParVarian(('','ni2'))[0])
                pass
            if(self.dim==4):
                self.ni2=int(self.GetParVarian(('','ni2'))[0])
                self.ni3=int(self.GetParVarian(('','ni3'))[0])
                pass
        elif self.dim in ('2p', '3p'):
            self.pseudo_axis_info = detect_varian_pseudo_axis(self.FidPath)
            self.nz = self.pseudo_axis_info.size


                    
        
        if(type(self.dim)!=str):
            if self.dim ==1:
                self.acqORD=''
            if(self.dim==2):
                self.acqORD=''
            if(self.dim==3):
                if(array[0]=='phase'): #phase,phase2
                    self.acqORD=1
                else: #phase2,phase
                    self.acqORD=0
            elif(self.dim==4):
                if(array[0]=='phase' and array[1]=='phase2'): #phase,phase2,phase3
                    self.acqORD=1
                elif(array[0]=='phase3' and array[1]=='phase2'): #phase3,phase2,phase
                    self.acqORD=0
        else:
            self.acqORD=''

        #-aqORD  aqOrd  [0] Acquisition Order Code:
        #0 = 2D d2,phase
        #0 = 3D d3,d2,phase2,phase
        #0 = 4D d4,d3,d2,phase3,phase2,phase
        #1 = 3D d3,d2,phase,phase2
        #1 = 4D d4,d3,d2,phase,phase2,phase3
        #2 = 4D d3,d2,d4,phase3,phase2,phase  #NOT COVERED




        self.sw=float(self.GetParVarian(('','sw'))[0])
        self.sw1=float(self.GetParVarian(('','sw1'))[0])
        if(type(self.dim)!=str):
            if(self._dim_greater_than(2)):
                self.sw2=float(self.GetParVarian(('','sw2'))[0])
            if(self.dim==2 and self.ni == 1):
                self.sw2=float(self.GetParVarian(('','sw2'))[0])
            if self.dim==2 and self.ni == 1:
                self.ni2 = int(self.GetParVarian(('', 'ni2'))[0])
            if(self.dim==4):
                self.sw3=float(self.GetParVarian(('','sw3'))[0])
        #elif(self.dim=='3p'):
        #    self.sw2=float(self.GetParVarian(('','sw2'))[0])

            

        self.sfrq=float(self.GetParVarian(('','sfrq'))[0])
        self.dfrq=float(self.GetParVarian(('','dfrq'))[0])
        self.dfrq2=float(self.GetParVarian(('','dfrq2'))[0])
        self.f1180_flg=self.GetParVarian(('','f1180'))[0]
        self.f2180_flg=self.GetParVarian(('','f2180'))[0]

        self.tn=self.GetParVarian(('','tn'))[0].split('"')[1]
        self.dn=self.GetParVarian(('','dn'))[0].split('"')[1]
        self.dn2=self.GetParVarian(('','dn2'))[0].split('"')[1]
        self.n1=self.tn
        self.n2=self.dn
        self.n3=self.dn2
        

        if(self.o1p=='Water' or self.o1p=='Auto'):
            self.temp=float(self.GetParVarian(('','temp'))[0])+273.19
            pass
            try:
                self.waterppm=self.WaterPPM()
                self.waterppmTOF=self.waterppm
                pass
            except:
                pass
                self.waterppm=4.7
                self.waterppmTOF=self.waterppm
            try:
                if(self.tn!='F19'):
                    pass
                    self.tof=float(self.GetParVarian(('','tof'))[0])
                    self.tof_me=float(self.GetParVarian(('','tof_me'))[0])
                    self.waterppm+=-(self.tof-self.tof_me)/self.sfrq
            except:
                pass
                pass
        else:
            self.waterppm=self.o1p
            self.waterppmTOF=self.waterppm

        if(type(self.dim)!=str):
            if(self._dim_greater_than(2)):
                self.dn2=self.GetParVarian(('','dn2'))[0].split('"')[1]
                pass
            if(self.dim==2 and self.ni == 1):
                self.dn2=self.GetParVarian(('','dn2'))[0].split('"')[1]
                pass
        elif(self.dim=='3p'):
            self.dn2=self.GetParVarian(('','dn2'))[0].split('"')[1]
            pass

        if(self.tn!='H1' and self.tn!='F19'):
            pass
            pass
            sys.exit()

        if(type(self.dim)!=str):
            if self._dim_greater_than(1):
                pass
                self.frq1,self.n1=self.GetNuc(self.labb[1])
            if(self.dim==3):
                self.frq2,self.n2=self.GetNuc(self.labb[2])
            if self._dim_greater_than(1):
                if(self.frq1==-1):
                    return
            if(self._dim_at_least(3)):
                if(self.frq2==-1):
                    return
            if(self.dim==4):
                self.frq3,self.n3=self.GetNuc(self.labb[3])
                if(self.frq3==-1):
                    return
            if self._dim_greater_than(1):
                self.f1ppm=self.shift(self.frq1,self.n1)
            if(self._dim_at_least(3)):
                self.f2ppm=self.shift(self.frq2,self.n2)
            if(self.dim==4):
                self.f3ppm=self.shift(self.frq3,self.n3)
        
        elif(self.dim=='3p'):
            self.frq1,self.n1=self.GetNuc(self.labb[1])
            if(self.frq1==-1):
                    return
            self.f1ppm=self.shift(self.frq1,self.n1)



        
        pass

        # if direct detect nucleus is F19, need to calculate central ppm for F19 rather than water
        if(self.tn=='F19'):
            self.tof=float(self.GetParVarian(('','tof'))[0])
            self.F19centre=(1.771789E-3)*self.tof - 111.26377
            
            # change self.waterppm to self.F19centre as not observing on proton
            self.waterppm=self.F19centre
            self.waterppmTOF=self.F19centre

            # enable second dimension ppm to be changed to correct value if 2D/pseudo3D/3D
            if(self.n1=='C13'):
                self.dof=float(self.GetParVarian(('','dof'))[0])
                self.f1ppm=(6.6301732E-3)*self.dof + 94.77335
                pass
            

        pass
        




    #parse an NUS schedule
    def GetNUSsamp(self):
        inny=open(self.FidPath+'/'+self.nuslist)
        samp=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                row=[]
                for i in range(len(test)):
                    row.append(int(test[i]))
                samp.append(row)
        samp=numpy.array(samp)

        maxy=numpy.max(samp,axis=0)
        pass
        pass
        self.ni=maxy[0]+1

        uni=self.ni*1. #uniform sampling size
        try:
            if self._dim_greater_than(2):
                self.ni2=maxy[1]+1
                uni*=self.ni2
            if(self.dim==4):
                self.ni3=maxy[2]+1
                uni*=self.ni3
        except:
            pass
            pass

        self.schedule=samp
        self.samp=len(samp)
        self.nusdim=len(samp[0])
        self.comp=float(self.samp)/float(uni)*100. #compression factor for sampling


    #parse a file to return a list
    def ParseList(self,infile):
        inny=open(infile)
        vdlist=[]
        for line in inny.readlines():
            test=line.split()
            vdlist.append(test[0])
        return vdlist
    


    #take internal variables and write fid.test.com
    def PipeParse(self):
        self.phasing=False

        self.spa=10               #spacing in output file
        self.axis='x','y','z','a' #names of axes
        self._pipe_rows = []

        if self.outdir:
            os.makedirs(self.outdir, exist_ok=True)
        if self.outdir:
            os.makedirs(self.outdir, exist_ok=True)
        self.outy=open(self.outdir+'/fid.test.com','w')
        self.outy.write('#!/bin/csh\n\n')


        if(self.tp=='var'):
            self.seqfil = self.GetParVarian('seqfil')[0]
            pass
             

        if(type(self.dim)==str): # pseudo-dimensional: 2p is pseudo2D, 3p is pseudo3D
            if(self.tp=='bruk'):
                info = self.pseudo_axis_info or detect_bruker_pseudo_axis(self.FidPath)
                self.pseudo_axis_info = info
                self.nz = info.size
            else:
                if(self.dim=='2p'):
                    z_variable = self.labb[1]
                    #self.nz = len(self.GetParVarian(z_variable))  
                    
                if(self.dim == '3p'):
                    # Varian pseudo3D data needs RelaxFix before var2pipe.
                    # Pseudo2D (2p) is already in the form expected by var2pipe
                    # and must be converted directly from the acquisition fid.
                    # Raw acquisition data remains in FidPath, while RelaxFix
                    # output belongs in outdir (the project's spec path).
                    relax_out = os.path.join(self.outdir, 'fid.final')
                    relax_in = os.path.join(self.FidPath, 'fid')
                    self.outy.write('RelaxFix.out {} {} {} 0 {} {}\n\n'.format(
                        self.np, self.ni, self.nz, relax_out, relax_in))

                #    z_variable= self.labb[2]


                #if(self.seqfil=='"HtoC_CH3_exchange_600_DC_dfh_v2_forAB"'):


        if(self.nuslist!=''):
            self.GetNUSsamp()
            if(self.tp=='var'):

                self.outy.write('nusExpand.tcl -mode varian -sampleCount %i -off 0 \\\n' % (self.samp))
                self.outy.write('-in %s/fid -out %s/fid_full -sample %s/%s -procpar %s/procpar\n' % (self.FidPath,self.outdir,self.FidPath,self.nuslist,self.FidPath))
            elif(self.tp=='bruk'):
                #
                self.outy.write('cp %s/acqus ./\n' % self.FidPath) #dodgy script. needs acqus
                self.outy.write('nusExpand.tcl -mode bruker -sampleCount %i -off 0 \\\n' % (self.samp))
                self.outy.write('-in %s/ser -out %s/ser_full -sample %s/%s\n' % (self.FidPath,self.outdir,self.FidPath,self.nuslist))
                #self.outy.write('rm ./acqus\n')
        # if(self.dim==1):
        #     self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/test.fid'))
        #     self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/test.fid',self.outdir+'/test.fid'))
        if(type(self.dim)!=str):
            if(self.dim==2 or self.dim==1):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/test.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/test.fid',self.outdir+'/test.fid'))
            if(self.dim==3):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/fids/test%03d.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/fids/test001.fid',self.outdir+'/fids'))
            if(self.dim==4):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/fids/test%03d%03d.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/fids/test001001.fid',self.outdir+'/fids'))
        if(self.dim=='2p'):
            self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/test.fid'))
            self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/test.fid',self.outdir+'/test.fid'))
        if(self.dim=='3p'):
            self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/fids/test%03d.fid'))
            self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/fids/test001.fid',self.outdir+'/fids'))


        if(self.tp=='omega'):
            pass
            #self.outy.write('bin2pipe -in %s -ge -neg \\\n' % (infile))
        elif(self.tp=='bruk'):
            if (self.dim=='2p' or self.dim=='3p' or self._dim_at_least(2)):
                infile=os.path.join(self.FidPath,'ser')
            else:
                infile=os.path.join(self.FidPath,'fid')
            if(self.nuslist!=''):
                infile=os.path.join(self.outdir,'ser_full')
            self.outy.write('bruk2pipe -in %s  \
' % (infile))
            GRPDLY=getattr(self,'GRPDLY',float(self.GetParBruk('acqus',('','GRPDLY'))[0]))
            DSPFVS=getattr(self,'DSPFVS',float(self.GetParBruk('acqus',('','DSPFVS'))[0]))
            DECIM=getattr(self,'DECIM',float(self.GetParBruk('acqus',('','DECIM'))[0]))

            BYTORDA=getattr(self,'BYTORDA',int(self.GetParBruk('acqus',('','BYTORDA'))[0]))
            if(BYTORDA==1):
                flag='noaswap'
            else:
                flag='aswap'

            DTYPA=getattr(self,'DTYPA',int(self.GetParBruk('acqus',('','DTYPA'))[0]))
            dtext=getattr(self,'dtext','')
            if(DTYPA==2 and dtext==''):
                dtext='-ws 8 -noi2f'

            AMXtext = self._set_bruker_format()

            #if wanting to remove digital filter, might need to set bad point threshold (normally 0)
            # -bad ' + str(self.nmrdata.bad_point_threshold) + ' -ext '

            self._write_bruker_bad_line(flag, AMXtext, DECIM, DSPFVS, GRPDLY, dtext=dtext)
           
        elif(self.tp=='var'):
            infile=os.path.join(self.FidPath,'fid')
            if(self.nuslist!=''):
                infile=os.path.join(self.outdir,'fid_full')
            #if(self.seqfil=='"HtoC_CH3_exchange_600_DC_dfh_v2_forAB"'):
            #    infile='fid.final'
            # Only pseudo3D passes through RelaxFix.  Varian pseudo2D must
            # feed the raw acquisition fid directly to var2pipe.
            if(self.dim == '3p'):
                infile=os.path.join(self.outdir,'fid.final')
            
            self.outy.write('var2pipe -in %s \\\n' % (infile))

            self.outy.write(' -noaswap ')
            if(self.acqORD!=''):
                self.outy.write(' -aqORD %i \\\n' % (self.acqORD))
            self.outy.write('\\\n')

        M=[]
        if(self.tp=='bruk'):
            M.append('DQD')
        else:
            M.append('Complex')

        self.modDict={}
        self.modDict['0']='Complex'
        self.modDict['1']='QF'
        self.modDict['2']='QSEQ'
        self.modDict['3']='TPPI'
        self.modDict['4']='States'
        self.modDict['5']='States-TPPI'
        self.modDict['6']='Echo-Antiecho'

        for i in range(len(self.rk)):
            if(self.rk[i]==0):
                if(self.tp=='bruk'):
                    pass
                    pass
                    M.append(self.modDict[self.mode[i]])
                else:
                    M.append('Complex')
            else:
                M.append('Rance-Kay')
        if(type(self.dim)==str):
            # A pseudo dimension is a real-valued axis, not another spectral
            # quadrature dimension.  rk can contain GUI entries for every
            # displayed axis, so simply appending Real can shift the pseudo
            # mode onto aMODE (e.g. 3p -> zMODE Complex, aMODE Real).  Trim
            # the spectral modes to the number of true spectral axes and put
            # Real explicitly on the final pseudo axis.
            pseudo_ndim = int(self.dim.split('p')[0])
            M = M[:pseudo_ndim - 1] + ['Real']
                
            
        if(type(self.dim)!=str):
            if(self.dim ==1):
                N=self.np,
                T=self.np/2,
                #M='Complex','Complex','Complex'
                sw=self.sw,
                O=self.sfrq,
                C=self.waterppm,

            if(self.dim==2):
                if self.ni > 1:
                    N=self.np,self.ni*2
                    T=self.np2,self.ni
                    #M='Complex','Complex','Complex'
                    sw=self.sw,self.sw1
                    O=self.sfrq,self.frq1
                    C=self.waterppm,self.f1ppm
                elif self.ni2 >1:
                    N=self.np,self.ni2*2
                    T=self.np2,self.ni2
                    #M='Complex','Complex','Complex'
                    sw=self.sw,self.sw2
                    O=self.sfrq,self.frq1
                    C=self.waterppm,self.f1ppm

            if(self.dim==3):
                N=self.np,self.ni*2,self.ni2*2
                T=self.np2,self.ni,self.ni2
                #M='Complex','Complex','Complex'
                sw=self.sw,self.sw1,self.sw2
                O=self.sfrq,self.frq1,self.frq2
                C=self.waterppm,self.f1ppm,self.f2ppm

            elif(self.dim==4):
                N=self.np,self.ni*2,self.ni2*2,self.ni3*2
                T=self.np2,self.ni,self.ni2,self.ni3
                #M='Complex','Complex','Complex','Complex'
                sw=self.sw,self.sw1,self.sw2,self.sw3
                O=self.sfrq,self.frq1,self.frq2,self.frq3
                C=self.waterppm,self.f1ppm,self.f2ppm,self.f3ppm
        elif(self.dim=='2p'):
            N=self.np,self.nz
            T=self.np2,self.nz
            #M='Complex','Complex','Complex'
            sw=self.sw,self.sw
            O=self.sfrq,self.sfrq
            C=self.waterppm,self.waterppm
        
        elif(self.dim=='3p'):
            if self.tp == 'bruk':
                # Build the two indirect axes from their *logical* Bruker
                # parameter files first, then map them into physical SER Y/Z
                # order once via AQSEQ.  This keeps size, T, mode, SW/OBS/CAR
                # and label attached to the same dimension as a unit.
                td2 = int(float(self.GetParBruk('acqu2s', ('', 'TD'))[0]))
                td3 = int(float(self.GetParBruk('acqu3s', ('', 'TD'))[0]))
                logical = {
                    'y': {
                        'N': td2, 'T': td2 / 2, 'MODE': self.modDict.get(str(self.yMode), 'Complex'),
                        'SW': self.sw1, 'OBS': self.frq1, 'CAR': self.f1ppm,
                        'LAB': self.labb[1],
                    },
                    'z': {
                        'N': td3, 'T': td3 / 2, 'MODE': self.modDict.get(str(self.zMode), 'Complex'),
                        'SW': self.sw2, 'OBS': self.frq2, 'CAR': self.f2ppm,
                        'LAB': self.labb[2],
                    },
                }

                pseudo_logical = getattr(self, 'pseudo_logical_axis', None)
                if pseudo_logical not in ('y', 'z'):
                    raise ValueError('Bruker pseudo-3D logical pseudo axis was not determined')
                spectral_logical = 'z' if pseudo_logical == 'y' else 'y'

                # The pseudo axis is real and has one stored value per list
                # entry.  Its spectral header fields are placeholders only;
                # use the real spectral indirect axis as the harmless source
                # while pseudo_axis.tsv remains authoritative for its values.
                spectral = logical[spectral_logical]
                pseudo = dict(logical[pseudo_logical])
                pseudo.update({
                    'N': self.nz, 'T': self.nz, 'MODE': 'Real',
                    'SW': spectral['SW'], 'OBS': spectral['OBS'], 'CAR': spectral['CAR'],
                })
                # Preserve the GUI pseudo label (normally ncyc) and spectral
                # label, independent of AQSEQ.  labb is logical X,Y,Z here.
                pseudo['LAB'] = self.labb[1 if pseudo_logical == 'y' else 2]
                spectral['LAB'] = self.labb[2 if pseudo_logical == 'y' else 1]

                # AQSEQ 312 exchanges logical Y/Z in the physical SER stream.
                # Do not perform any second generic AQSEQ swap below.
                if self.aqseq == '312':
                    physical_y = logical['z']
                    physical_z = logical['y']
                    physical_y = pseudo if pseudo_logical == 'z' else spectral
                    physical_z = pseudo if pseudo_logical == 'y' else spectral
                else:
                    physical_y = pseudo if pseudo_logical == 'y' else spectral
                    physical_z = pseudo if pseudo_logical == 'z' else spectral

                N = (self.np, physical_y['N'], physical_z['N'])
                T = (self.np2, physical_y['T'], physical_z['T'])
                M = (M[0], physical_y['MODE'], physical_z['MODE'])
                sw = (self.sw, physical_y['SW'], physical_z['SW'])
                O = (self.sfrq, physical_y['OBS'], physical_z['OBS'])
                C = (self.waterppm, physical_y['CAR'], physical_z['CAR'])
                self.labb = (self.labb[0], physical_y['LAB'], physical_z['LAB'])
            elif self.ni > 1:
                N=self.np,self.ni*2, self.nz
                T=self.np2,self.ni, self.nz
                sw=self.sw,self.sw1, self.sw1
                yobs = self.dfrq
                O=self.sfrq,yobs,yobs
                C=self.waterppm,self.f1ppm,self.f1ppm
            elif self.ni2 >1:
                N=self.np,self.ni2*2, self.nz
                T=self.np2,self.ni2, self.nz
                sw=self.sw,self.sw2, self.sw2
                yobs = self.dfrq
                O=self.sfrq,yobs,yobs
                C=self.waterppm,self.f1ppm,self.f1ppm

        # T is the final time-domain size tuple that will be written to the
        # NMRPipe conversion script.  At this point all vendor-specific logic
        # has already been applied (Varian procpar interpretation, Bruker
        # TD/NusTD handling, acquisition-order swaps and pseudo dimensions).
        # Store exactly these final values as the canonical xN/yN/zN/aN state
        # so other consumers never need to duplicate the vendor-specific rules.
        self._set_time_domain_sizes(T)
        self._set_observation_frequencies(O)

        self.AddPipe(self.outy,self.axis,'N',N,self.spa)
        self.AddPipe(self.outy,self.axis,'T',T,self.spa)
        self.AddPipe(self.outy,self.axis,'MODE',M,self.spa)
        self.AddPipe(self.outy,self.axis,'SW',sw,self.spa)
        self.AddPipe(self.outy,self.axis,'OBS',O,self.spa)

        
        self.AddPipe(self.outy,self.axis,'CAR',C,self.spa)
        # else:
        #     if(type(self.dim)!=str):
        #         self.AddPipe(self.outy,self.axis,'CAR',numpy.zeros(self.dim),self.spa)
        #     else:
        #         dim = int(float(self.dim.split('p')[0]))
        #         self.AddPipe(self.outy,self.axis,'CAR',numpy.zeros(dim),self.spa)
        self.AddPipe(self.outy,self.axis,'LAB',self.labb,self.spa)

        if(self.dim==1):
            self.outy.write(' -ndim  %s \\\n' % (str(self.dim).ljust(self.spa)))
            self.outy.write('  -out $ft4trec -verb -ov\n')
        elif(type(self.dim)!=str):
            self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(self.dim).ljust(self.spa),'States'.ljust(self.spa)))
            if(self.dim==2):
                self.outy.write('  -out $ft4trec -verb -ov\n')
            if(self.dim==3):
                self.outy.write('  -out $ft4trec -verb -ov\n')
            if(self.dim==4):
                self.outy.write('| pipe2xyz -x -out $ft4trec -verb -ov -to 0\n')
        else:

            if(self.tp=='bruk'):
                self.outy.write('  -ndim               3  -aq2D         Complex                         \\\n')
                if getattr(self, 'pseudo_acq_axis', 'z') == 'y':
                    # TP-ZTP-TP swaps Y and Z while preserving X.  This is
                    # needed only when the pseudo dimension was acquired as Y.
                    self.outy.write('| nmrPipe -fn TP  -exch -noord -nohdr \\\n')
                    self.outy.write('| nmrPipe -fn ZTP -exch -noord \\\n')
                    self.outy.write('| nmrPipe -fn TP  -exch -noord -nohdr \\\n')
                # Acquired-Z pseudo axes are already outermost.
                self.outy.write('| pipe2xyz -x -out $ft4trec -ov \n')


            else:
                self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (self.dim.split('p')[0].ljust(self.spa),'States'.ljust(self.spa)))
                self.outy.write('  -out $ft4trec -verb -ov\n')


        """
        if(self.nuslist!=''):
            if(self.tp=='var'):
                infile='fid_full'
            elif(self.tp=='bruk'):
                infile='ser_full'
            else:
                infile='a'
            #self.outy.write('rm %s/%s\n\n' % (self.outdir,infile))
        """

        # if(self.dim!=1 and self.dim!='2p'):
        #     self.phasing=True
        #     self.MakePhasing1D()
        #     self.phasing=False


        self.outy.close()
        #self.outy.write('| nmrPipe -ov -verb -out test.fid\n') #spit into a giant fid

        return


    ##############################################################################
    #extract a 1D for indirect phasing.
    #experimental at this stage.
    def MakePhasing1D(self):
        self.outy.write('\n\n')
        self.outy.write('# Making Phasing 1D\n')

        if(type(self.dim)!=str):
            if(self.dim==2 or self.dim==1):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/Phasing1D.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/Phasing1D.fid',self.outdir+'/Phasing1D.fid'))
            if(self.dim==3):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/Phasing1D.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/Phasing1D.fid',self.outdir+'/Phasing1D.fid'))
            if(self.dim==4):
                self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/Phasing1D.fid'))
                self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/Phasing1D.fid',self.outdir+'/Phasing1D.fid'))
        if(self.dim=='2p'):
            self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/Phasing1D.fid'))
            self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/Phasing1D.fid',self.outdir+'/Phasing1D.fid'))
        if(self.dim=='3p'):
            self.outy.write('%s\n' % ('set ft4trec='+self.outdir+'/Phasing1D.fid'))
            self.outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/Phasing1D.fid',self.outdir+'/Phasing1D.fid'))


        if(self.tp=='omega'):
            pass
            #self.outy.write('bin2pipe -in %s -ge -neg \\\n' % (infile))
        elif(self.tp=='bruk'):
            if self._dim_at_least(2):
                infile='ser'
            else:
                infile='fid'
            if(self.nuslist!=''):
                infile='ser_full'
            self.outy.write('bruk2pipe -in %s/%s  \\\n' % (self.outdir,infile))
            GRPDLY=float(self.GetParBruk('acqus',('','GRPDLY'))[0])
            DSPFVS=int(self.GetParBruk('acqus',('','DSPFVS'))[0])
            DECIM=float(self.GetParBruk('acqus',('','DECIM'))[0])

            BYTORDA=int(self.GetParBruk('acqus',('','BYTORDA'))[0])
            if(BYTORDA==1):
                flag='noaswap'
            else:
                flag='aswap'

            AMXtext = self._set_bruker_format()

            self._write_bruker_bad_line(flag, AMXtext, DECIM, DSPFVS, GRPDLY, decim_as_int=True)
        elif(self.tp=='var'):
            infile='fid'
            if(self.nuslist!=''):
                infile='fid_full'
            if(self.seqfil=='"HtoC_CH3_exchange_600_DC_dfh_v2_forAB"'):
                infile='fid.final'
            
            self.outy.write('var2pipe -in %s/%s \\\n' % (self.outdir,infile))

            self.outy.write(' -noaswap ')
            if(self.acqORD!=''):
                self.outy.write(' -aqORD %i \\\n' % (self.acqORD))
            self.outy.write('\\\n')

        M=[]
        if(self.tp=='bruk'):
            M.append('DQD')
        else:
            M.append('Complex')

        self.modDict={}
        self.modDict['0']='Complex'
        self.modDict['1']='QF'
        self.modDict['2']='QSEQ'
        self.modDict['3']='TPPI'
        self.modDict['4']='States'
        self.modDict['5']='States-TPPI'
        self.modDict['6']='Echo-Antiecho'

        for i in range(len(self.rk)):
            if(self.rk[i]==0):
                if(self.tp=='bruk'):
                    pass
                    pass
                    M.append(self.modDict[self.mode[i]])
                else:
                    M.append('Complex')
            else:
                M.append('Rance-Kay')
        if(type(self.dim)==str):
            # A pseudo dimension is a real-valued axis, not another spectral
            # quadrature dimension.  rk can contain GUI entries for every
            # displayed axis, so simply appending Real can shift the pseudo
            # mode onto aMODE (e.g. 3p -> zMODE Complex, aMODE Real).  Trim
            # the spectral modes to the number of true spectral axes and put
            # Real explicitly on the final pseudo axis.
            pseudo_ndim = int(self.dim.split('p')[0])
            M = M[:pseudo_ndim - 1] + ['Real']
                
        M = M[0]

        if(type(self.dim)!=str):
            if(self.dim ==1):
                N=self.np
                T=self.np/2
                #M='Complex','Complex','Complex'
                sw=self.sw
                O=self.sfrq
                C=self.waterppm

            if(self.dim==2):
                N=self.np
                T=self.np/2
                #M='Complex','Complex','Complex'
                sw=self.sw
                O=self.sfrq
                C=self.waterppm


            if(self.dim==3):
                N=self.np
                T=self.np/2
                #M='Complex','Complex','Complex'
                M=M[0]
                sw=self.sw
                O=self.sfrq
                C=self.waterppm

            elif(self.dim==4):
                N=self.np
                T=self.np/2
                #M='Complex','Complex','Complex','Complex'
                sw=self.sw
                O=self.sfrq
                C=self.waterppm
        elif(self.dim=='2p'):
            N=self.np
            T=self.np/2
            #M='Complex','Complex','Complex'
            sw=self.sw
            O=self.sfrq
            C=self.waterppm
        
        elif(self.dim=='3p'):
            N=self.np
            T=self.np/2
            #M='Complex','Complex','Complex'
            sw=self.sw
            O=self.sfrq
            C=self.waterppm

        #print([T])
        #sys.exit(100)
        self.AddPipe(self.outy,self.axis,'N',[N],self.spa)
        self.AddPipe(self.outy,self.axis,'T',[T],self.spa)
        self.AddPipe(self.outy,self.axis,'MODE',[M],self.spa)
        self.AddPipe(self.outy,self.axis,'SW',[sw],self.spa)
        self.AddPipe(self.outy,self.axis,'OBS',[O],self.spa)

        
        self.AddPipe(self.outy,self.axis,'CAR',[C],self.spa)
        # else:
        #     if(type(self.dim)!=str):
        #         self.AddPipe(self.outy,self.axis,'CAR',numpy.zeros(self.dim),self.spa)
        #     else:
        #         dim = int(float(self.dim.split('p')[0]))
        #         self.AddPipe(self.outy,self.axis,'CAR',numpy.zeros(dim),self.spa)
        self.AddPipe(self.outy,self.axis,'LAB',self.labb,self.spa)

        if(type(self.dim)!=str):
            self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(self.dim).ljust(self.spa),'States'.ljust(self.spa)))
            if(self.dim==2 or self.dim==1):
                self.outy.write('  -out $ft4trec -verb -ov\n')
            if(self.dim==3):
                self.outy.write('  -out $ft4trec -verb -ov\n')
            if(self.dim==4):
                self.outy.write('| pipe2xyz -x -out $ft4trec -verb -ov -to 0\n')
        else:
            self.outy.write(' -ndim  %s -aq2D  %s \\\n' % (self.dim.split('p')[0].ljust(self.spa),'States'.ljust(self.spa)))
            self.outy.write('  -out $ft4trec -verb -ov\n')






    ############################################################################
    # functions to help with writing files
    # add an entry to fid.test.com
    def _publish_time_domain_sizes(self):
        """Publish final conversion sizes as xN/yN/zN/aN.

        Varian and Bruker deliberately keep their existing native parsing paths.
        This method runs afterwards and maps their resolved state to the axis
        sizes used by NMRPipe.  PipeParse later overwrites these from its final T
        tuple, providing an additional guarantee that generated scripts and the
        public xyza state agree.
        """
        try:
            dim = self.dim
            if type(dim) != str:
                if dim == 1:
                    values = (self.np2,)
                elif dim == 2:
                    indirect = self.ni if getattr(self, 'ni', 0) > 1 else self.ni2
                    values = (self.np2, indirect)
                elif dim == 3:
                    values = (self.np2, self.ni, self.ni2)
                elif dim == 4:
                    values = (self.np2, self.ni, self.ni2, self.ni3)
                else:
                    return
            elif dim == '2p':
                values = (self.np2, self.nz)
            elif dim == '3p':
                indirect = self.ni if getattr(self, 'ni', 0) > 1 else self.ni2
                values = (self.np2, indirect, self.nz)
                if self.tp == 'bruk' and getattr(self, 'aqseq', '') == '312':
                    values = (values[0], values[2], values[1])
            else:
                return
            self._set_time_domain_sizes(values)
        except (AttributeError, TypeError, ValueError):
            # Preserve historical conversion behaviour if a partial/legacy data
            # set does not provide enough information to publish every size.
            pass

    def _set_time_domain_sizes(self, values):
        """Store final NMRPipe time-domain sizes as xN/yN/zN/aN.

        ``values`` is the final T tuple used by the conversion script.  Keeping
        this assignment here deliberately avoids changing the established
        Varian and Bruker conversion calculations themselves.
        """
        attrs = ('xN', 'yN', 'zN', 'aN')
        for attr in attrs:
            setattr(self, attr, None)
        for attr, value in zip(attrs, values):
            try:
                number = float(value)
                value = int(number) if number.is_integer() else number
            except (TypeError, ValueError):
                pass
            setattr(self, attr, value)

    def GetTimeDomainSizes(self):
        """Return canonical conversion sizes using the x/y/z/a axis names.

        Values correspond exactly to the final ``-xT/-yT/-zT/-aT`` sizes
        written by the current conversion path.  Missing dimensions are omitted.
        """
        result = {}
        for axis in ('x', 'y', 'z', 'a'):
            value = getattr(self, axis + 'N', None)
            if value not in (None, ''):
                result[axis] = value
        return result

    def _publish_observation_frequencies(self):
        """Publish final conversion OBS values as xOBS/yOBS/zOBS/aOBS.

        This mirrors _publish_time_domain_sizes(): existing Varian and Bruker
        frequency calculations remain untouched.  We only expose the resolved
        frequencies using NMRPipe axis names, including the established Bruker
        312 acquisition-order swap and pseudo-dimensional conventions.
        """
        try:
            dim = self.dim
            if type(dim) != str:
                count = int(dim)
                values = (self.sfrq, getattr(self, 'frq1', None),
                          getattr(self, 'frq2', None), getattr(self, 'frq3', None))[:count]
            elif dim == '2p':
                values = (self.sfrq, self.sfrq)
            elif dim == '3p':
                values = (self.sfrq, self.sfrq, self.sfrq)
            else:
                return
            if self.tp == 'bruk' and getattr(self, 'aqseq', '') == '312' and len(values) >= 3:
                values = (values[0], values[2], values[1]) + tuple(values[3:])
            self._set_observation_frequencies(values)
        except (AttributeError, TypeError, ValueError):
            pass

    def _set_observation_frequencies(self, values):
        """Store final NMRPipe observation frequencies as x/y/z/a OBS (MHz)."""
        attrs = ('xOBS', 'yOBS', 'zOBS', 'aOBS')
        for attr in attrs:
            setattr(self, attr, None)
        for attr, value in zip(attrs, values):
            try:
                value = float(value)
            except (TypeError, ValueError):
                pass
            setattr(self, attr, value)

    def GetObservationFrequencies(self):
        """Return canonical NMRPipe observation frequencies in MHz by axis."""
        result = {}
        for axis in ('x', 'y', 'z', 'a'):
            value = getattr(self, axis + 'OBS', None)
            if value not in (None, ''):
                result[axis] = value
        return result

    def AddPipeLine(self,outy,lab,par,val,spa):
        outy.write(' -%s%s %s' % (lab, par.ljust(5), str(val).ljust(spa)))

    #end a line.
    def EndPipeLine(self,outy):
        outy.write(' ' + chr(92) + '\n')

    def _flushPipeRows(self, outy, spa):
        rows = getattr(self, '_pipe_rows', None)
        if not rows:
            return
        for row in rows:
            if not row:
                continue
            line = ''
            for lab, par, val in row:
                line += ' -%s%s %s' % (lab, par.ljust(5), str(val).ljust(spa))
            outy.write(line.rstrip() + ' ' + chr(92) + '\n')
        self._pipe_rows = []

    #add a complete line. to fid.test.com
    def AddPipe(self,outy,axis,par,vals,spa):
        if not hasattr(self, '_pipe_rows') or self._pipe_rows is None:
            self._pipe_rows = []

        if(self.phasing==False):
            if(type(self.dim)!=str):
                count = self.dim
            elif(self.dim=='2p'):
                count = 2
            elif(self.dim=='3p'):
                count = 3
            else:
                count = 1
        else:
            count = 1

        for i in range(count):
            if len(self._pipe_rows) <= i:
                self._pipe_rows.append([])
            value = vals[i] if count > 1 else vals[0]
            self._pipe_rows[i].append((axis[i], par, value))

        if par == 'LAB':
            self._flushPipeRows(outy, spa)

    ############################################################################
    #return water chemical shift in range 0-100oC (self.temp is in K)
    def WaterPPM(self):
        return 5.060 - 0.0122*(self.temp-273.19) + (2.11E-5)*(self.temp-273.19)**2.
        #return 7.83 - (self.temp/96.9)

    ##########################################################################
    #What type of spectrometer was used?
    def GetSpectrometerType(self,path='./'):
        # Keep the legacy method/API, but delegate identification to the
        # side-effect-free common inspector used by SpinHub.
        detected = detect_spectrometer_type(path)
        self.tp = detected if detected is not None else 'omega'
        if self.tp == 'var':
            self.parfile = os.path.join(path, 'procpar')

            
    #Get acqisition order.
    def GetAcqseq(self):
        """Return Bruker pulse-program acquisition order (for example 312).

        The pulseprogram belongs to the raw acquisition directory (FidPath).
        Older code looked only in ``outdir``; in the GUI that can be the
        conversion/output directory, silently losing AQSEQ and therefore the
        Y/Z ordering of pseudo-3D data.
        """
        candidates = []
        for base in (getattr(self, 'FidPath', None), getattr(self, 'outdir', None)):
            if not base:
                continue
            for name in ('pulseprogram', 'pulseprogram.precomp'):
                path = os.path.join(base, name)
                if path not in candidates:
                    candidates.append(path)
        for path in candidates:
            if not os.path.isfile(path):
                continue
            with open(path, errors='replace') as inny:
                for line in inny:
                    test = line.split()
                    if len(test) > 1 and test[0].lower() == 'aqseq':
                        return test[1].strip()
        return -1

    #get pulse sequence
    def GetSequence(self):
        if(self.tp=='var'):
            seqfil=self.GetParVarian(('','seqfil'))[0].split('"')[1]
        elif(self.tp=='bruk'):
            seqfil=self.GetParBruk('acqus',('','PULPROG',))[0].replace('<', '').replace('>', '').strip()
        elif(self.tp=='omeg'):
            parfile=self.GetOmegaParFile()
            test=self.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
            seqfil=test[len(test)-1]
        self.seqfil=seqfil
        pass

    #parse a file to a list.
    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=line.split()
            peak.append(linetosave)
        peakfile.close()
        return peak


    #analyse either acqu and acqu2
    def GetParBruk(self,infile,argv,verb='n',):
        verb='y'
        args=[]
        full_path=self.FidPath+'/'+infile
        #print(f'[GetParBruk] checking {full_path}: exists={os.path.exists(full_path)}', flush=True)
        if not os.path.exists(full_path):
            pass
            return (0,)
        procpar=self.readfile(full_path)
        for i in range(len(argv)-1):
            param=argv[i+1]
            tick=0
            for j in range(len(procpar)):
                #print procpar[j]
                test=procpar[j][0].split('##$')

                if(len(test)>1):
                    test2=test[1].split('=')[0]
                    if(test2==param):
                        if(verb=='y'):
                            pass
                        args.append(procpar[j][1])
                        tick=1
                else:
                    #we have a line of zeros
                    #is the previous line what we're after?
                    test=procpar[j-1][0].split('##$')
                    if(len(test)>1):
                        test2=test[1].split('=')[0]
                        for i in range(100):
                            parT=test2+str(i)
                            if(parT==param):
                                if(len(param.split(test2))>1):
                                    if(verb=='y'):
                                        pass
                                    #parameters are in rows in j,j+1,j+2...
                                    go=0
                                    cnt=0
                                    while(go==0):
                                        if(i<len(procpar[j+cnt])):
                                            val=procpar[j+cnt][i]
                                            go=1
                                        else:
                                            i-=len(procpar[j+cnt])
                                            cnt+=1

                                    args.append(val)
                                    if(verb=='y'):
                                        pass
                                    tick=1
                #sys.exit(100)
            if(tick==0):
                if(verb=='y'):
                    pass
                return 0,
            else:
                return args

    def GetParVarian(self,argv,verb='n'):
        verb='y'
        args=[]
        procpar=self.readfile(self.parfile)
        if(argv=='time_T2' or argv=='ncyc' or argv=='gzlvl1' or argv=='ncyc_cp' or argv=='time_T1' or argv=='select_flg' or argv=='seqfil'):
            for j in range(len(procpar)):
                if(procpar[j][0]==argv):
                    if(verb=='y'):
                        pass
                    if(int(procpar[j+1][0])>1):
                        pass
                    pass
                    for k in range(int(procpar[j+1][0])):
                        try:
                            if(verb=='y'):
                                pass
                            args.append(procpar[j+1][k+1])
                        except:
                            if(verb=='y'):
                                pass
                            args.append(procpar[j+1+k][0])
                    if(verb=='y'):
                        pass
                    tick=1
        else:
            for i in range(len(argv)-1):
                param=argv[i+1]
                tick=0
                for j in range(len(procpar)):
                    if(procpar[j][0]==param):
                        if(verb=='y'):
                            pass
                        if(int(procpar[j+1][0])>1):
                            pass
                        pass
                        for k in range(int(procpar[j+1][0])):
                            try:
                                if(verb=='y'):
                                    pass
                                args.append(procpar[j+1][k+1])
                            except:
                                if(verb=='y'):
                                    pass
                                args.append(procpar[j+1+k][0])
                        if(verb=='y'):
                            pass
                        tick=1
        if(tick==0):
            if(verb=='y'):
                pass
            return 'fail'
        else:
            return args


    def GetParOmega(self,infile,verb,argv):
        args=[]
        procpar=self.readfile(infile)
        for i in range(len(argv)-1):
            param=argv[i+1]
            tick=0
            for j in range(len(procpar)):
                if(procpar[j][0]==param):
                    if(verb=='y'):
                        pass
                    if(len(procpar[j])>2):
                        for k in range(len(procpar[j])-2):
                            if(verb=='y'):
                                pass
                            args.append(procpar[j][k+1])
                    else:
                        for k in range(len(procpar[j])-1):
                            if(verb=='y'):
                                pass
                            args.append(procpar[j][k+1])


                    if(verb=='y'):
                        pass
                    tick=1
        if(tick==0):
            if(verb=='y'):
                pass
            return 'fail'
        else:
            return args


    def GetOmegaVal(self,infile,param):
        if(self.GetParOmega(infile,'n',('',param))!='Fail'):
            return self.GetParOmega(infile,'n',('',param))


    def GetBrukVal(self,infile,param):
        test=self.GetParBruk(infile,'n',('',param))
        if(test!='Fail'):
            return test

    #################################################################
    #taken from datasearch.py
    def GetTime(self):
        tim=''
        if(self.tp=='var'):

            t=self.GetParVarian((self.outdir,'time_complete'))
            
            if(len(t)!=0 and t!='fail'):

                tim=self.GetParVarian((self.outdir,'time_complete'))[0].split('"')[1]
                tim=tim[:4]+'/'+tim[4:6]+'/'+tim[6:8]+' '+tim[9:11]+':'+tim[11:13]+':'+tim[13:15]
            else:
                test=os.path.join(self.outdir,'log')
                if(os.path.exists(test)):
                    inny=open(test)
                    for line in inny.readlines():
                        test=line.split(': ')
                        if(len(test)>1):
                            tim=test[0]
            
                
            #tim=yr+' '+tim.split('T')[1]
            
        elif(self.tp=='bruk'):
            tim=self.GetTimeBruk(self.outdir+'/acqus')
            tim=tim.replace('-','/').split('.')[0]
        return tim
        
    def GetTimeBruk(self,infile):
        procpar=self.readfile(infile)
        for j in range(len(procpar)):        
            #print (procpar[j])
            if(procpar[j][0]=='$$'):
                return procpar[j][1]+' '+procpar[j][2]
        return 'False'

    
    @staticmethod
    def _parse_audit_datetime(text):
        """Parse date/time spellings used in spectrometer audit logs."""
        value = str(text).strip().strip('"').strip()
        value = re.sub(r'\s+', ' ', value)
        try:
            return datetime.fromisoformat(value.replace('Z', '+00:00'))
        except Exception:
            pass
        for fmt in (
                '%Y-%m-%d %H:%M:%S', '%Y/%m/%d %H:%M:%S',
                '%Y-%m-%d %H:%M', '%Y/%m/%d %H:%M',
                '%d-%b-%Y %H:%M:%S', '%d-%b-%Y %H:%M',
                '%a %b %d %H:%M:%S %Y', '%b %d %H:%M:%S %Y'):
            try:
                return datetime.strptime(value, fmt)
            except ValueError:
                pass
        match = re.search(r'(\d{4}[-/]\d{1,2}[-/]\d{1,2}[ T]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?)', value)
        if match:
            try:
                return datetime.fromisoformat(match.group(1).replace('/', '-'))
            except Exception:
                pass
        return None

    @staticmethod
    def FormatMeasurementTime(seconds):
        """Return elapsed time as days, hours and minutes."""
        try:
            total_minutes = max(0, int(round(float(seconds) / 60.0)))
        except Exception:
            return ''
        days, rem = divmod(total_minutes, 24 * 60)
        hours, minutes = divmod(rem, 60)
        parts = []
        if days:
            parts.append('%d %s' % (days, 'day' if days == 1 else 'days'))
        if hours:
            parts.append('%d %s' % (hours, 'hour' if hours == 1 else 'hours'))
        if minutes or not parts:
            parts.append('%d %s' % (minutes, 'minute' if minutes == 1 else 'minutes'))
        return ', '.join(parts)

    @staticmethod
    def _parse_bruker_audit_timestamp(value):
        """Parse a TopSpin audita.txt acquisition timestamp, including UTC offset."""
        value = str(value).strip().rstrip(',').strip()
        for fmt in (
                '%Y-%m-%d %H:%M:%S.%f %z',
                '%Y-%m-%d %H:%M:%S %z'):
            try:
                return datetime.strptime(value, fmt)
            except ValueError:
                pass
        return None

    @classmethod
    def _bruker_completed_acquisition(cls, audit_file):
        """Return (start, completion) for the last completed TopSpin acquisition.

        TopSpin audita.txt records begin ``( NUMBER, ...`` and may span many
        physical lines.  Acquisition attempts that were stopped contain
        ``started at`` plus ``terminated by command``; a successfully completed
        acquisition contains both ``started at`` and ``completed at``.  Search
        records from newest to oldest and use the last completed acquisition.
        """
        try:
            with open(audit_file, 'r', errors='replace') as handle:
                text = handle.read()
        except OSError:
            return None, None

        starts = list(re.finditer(r'(?m)^\s*\(\s*(\d+)\s*,', text))
        if not starts:
            return None, None

        records = []
        for i, match in enumerate(starts):
            stop = starts[i + 1].start() if i + 1 < len(starts) else len(text)
            records.append((int(match.group(1)), text[match.start():stop]))

        stamp = r'(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}(?:\.\d+)?\s+[+-]\d{4})'
        for _number, block in reversed(records):
            started = re.search(r'(?i)started\s+at\s+' + stamp, block)
            completed = re.search(r'(?i)completed\s+at\s+' + stamp, block)
            if not (started and completed):
                continue
            start_dt = cls._parse_bruker_audit_timestamp(started.group(1))
            end_dt = cls._parse_bruker_audit_timestamp(completed.group(1))
            if start_dt is not None and end_dt is not None and end_dt >= start_dt:
                return start_dt, end_dt
        return None, None

    def GetMeasurementTime(self):
        """Return actual elapsed acquisition time from the spectrometer audit log."""
        raw_dir = self.FidPath or self.outdir
        if self.tp == 'var':
            logfile = os.path.join(raw_dir, 'log')
            if not os.path.isfile(logfile):
                return ''
            stamps = []
            with open(logfile, 'r', errors='replace') as handle:
                for line in handle:
                    dt = self._parse_audit_datetime(line)
                    if dt is None:
                        dt = self._parse_audit_datetime(line.split(': ', 1)[0])
                    if dt is not None:
                        stamps.append(dt)
            if len(stamps) >= 2:
                return self.FormatMeasurementTime((max(stamps) - min(stamps)).total_seconds())
            return ''

        if self.tp == 'bruk':
            # audita.txt lives beside acqus.  FidPath is normally that directory,
            # but retain outdir as a compatibility fallback for older projects.
            candidates = []
            for directory in (raw_dir, self.outdir):
                if directory:
                    candidate = os.path.join(directory, 'audita.txt')
                    if candidate not in candidates:
                        candidates.append(candidate)
            audit = next((x for x in candidates if os.path.isfile(x)), None)
            if audit is None:
                return ''
            start_dt, end_dt = self._bruker_completed_acquisition(audit)
            if start_dt is None or end_dt is None:
                return ''
            return self.FormatMeasurementTime((end_dt - start_dt).total_seconds())
        return ''

    def GetTemp(self):
        if(self.tp=='var'):
            #print(self.GetParVarian((self.outdir,'temp')))[0]
            self.temp=float(self.GetParVarian((self.outdir,'temp'))[0])+273.19 #.split('"')[1]
        elif(self.tp=='bruk'):
            self.temp=self.GetParBruk('acqus',('','TE',))[0]

    def GetP1(self):
        if(self.tp=='var'):
            self.pw=self.GetParVarian((self.outdir,'pw'))[0] #.split('"')[1]
        elif(self.tp=='bruk'):
            self.pw=self.GetParBruk('acqus',('','P1',))[0]

            
    def GetSfrq(self):
        if(self.tp=='var'):
            self.sfrq=float(self.GetParVarian((self.outdir,'sfrq'))[0])  
        elif(self.tp=='bruk'):
            self.sfrq=float(self.GetParBruk('acqus',('','SFO1',))[0])

    def Write(self):

        self.Convert() #running conversions.

        spec={}
        spec['bruk']='bruker'
        spec['var']='varian'
        self.GetTemp()
        self.GetP1()

        
        self.labels={} #initialise labels.

        self.labels['Time completed']='%s ' % self.GetTime()
        self.labels['Spectrometer OS']='%s' % spec[self.tp]
        self.labels['Pulse sequence']='\\verb+%s+' % self.seqfil
        self.labels['Temperature']='%s K ' % self.temp
        self.labels['90$^o$ pulse time:'] = '%s $\\mu$s' % self.pw
        self.labels['Data path']='%s' % self.outdir
        self.labels['Reference mode']='%s' % self.o1p

        self.texTab=[]

        

        if(self.tp=='bruk'):
            nt=int(self.GetParBruk('acqus',('','NS',))[0]) #[0].replace('<', '').replace('>', '').strip()
            d1=float(self.GetParBruk('acqus',('','D1',))[0])
        elif(self.tp=='var'):
            nt=int(self.GetParVarian((self.outdir,'nt'))[0])
            d1=float(self.GetParVarian((self.outdir,'d1'))[0])

            
        self.labels['NS']='%i' % nt
        self.labels['D1']='%.2f s' % d1

        #print(self.labb,self.n1,self.n2)
        #sys.exit(100)
        at=self.np2/self.sw
        #self.labels[cnt]='Direct: %s aq: %.2f ms np: %i sw: %.2f Hz swp: %.2f ppm car: %.2f ppm' % (self.n1,self.np2/self.sw*1000,self.np2,self.sw,self.sw/self.sfrq,self.waterppm);cnt+=1
        self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.labb[0],self.sfrq,at*1000,self.np2,self.sw,self.sw/self.sfrq,self.waterppm))

        dim=1
        incr=1
        self.nuslist='nuslist'            
        if(os.path.exists(os.path.join(self.outdir,self.nuslist))):
            self.GetNUSsamp()
            self.labels['NUS']='Schedule length: %i Compression: %.2f%s' % (self.samp,self.comp,'\\%')

            dim+=self.nusdim
                    
            #self.labels[cnt]='Indirect: %s at: %.2f ms ni: %i sw: %.2f Hz car: %.2f ppm' % (self.n2,self.ni/self.sw1*1000,self.ni,self.sw1,self.f1ppm);cnt+=1
            self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.labb[1],self.frq1,self.ni/self.sw1*1000,self.ni,self.sw1,self.sw1/self.frq1,self.f1ppm))
            if(self.nusdim>=2):
                #self.labels[cnt]='Indirect2: %s at: %.2f ms ni: %i sw: %.2f Hz car: %.2f ppm' % (self.n3,self.ni2/self.sw2*1000,self.ni2,self.sw2,self.f2ppm);cnt+=1
                self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.labb[2],self.frq2,self.ni2/self.sw2*1000,self.ni2,self.sw2,self.sw2/self.frq2,self.f2ppm))                    
            incr*=self.samp
            #not yet handling 4Ds.
                    
        else:
            #uniformly sampled data.
            if(self.ni>1):
                dim+=1
                incr*=self.ni
                #self.labels[cnt]='Indirect: %s at: %.2f ms ni: %i sw: %.2f Hz car: %.2f ppm' % (self.n2,ni/sw1*1000,ni,sw1,self.f1ppm);cnt+=1
                self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.labb[1],self.frq1,self.ni/self.sw1*1000,self.ni,self.sw1,self.sw1/self.frq1,self.f1ppm))

            try:
                self.ni2
            except:
                self.ni2=1
                pass
            
            if(self.ni2>1):
                dim+=1
                incr*=self.ni2
                #self.labels[cnt]='Indirect2: %s at: %.2f ms ni: %i sw: %.2f Hz car: %.2f ppm' % (self.n3,ni2/sw2*1000,ni2,sw2,self.f2ppm);cnt+=1
                self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.labb[2],self.frq2,self.ni2/self.sw2*1000,self.ni2,self.sw2,self.sw2/self.frq2,self.f2ppm))                                                                            
        pseudo=False
        nz=1
        if(self.tp=='bruk'):

            if(os.path.exists('vdlist')):
                pseudo=True
                self.labels['Pseudo']='PseudoAxis found';cnt+=1
                pass

        elif(self.tp=='var'):
            arrs=self.GetParVarian((self.outdir,'array'))

            pass
            if(len(arrs)!=0 and arrs[0]!='""'):
            
                arrs=arrs[0].split('"')[1].split(',')
                arrStr=''
                for i,arr  in enumerate(arrs):
                    if(i!=0):
                        arrStr+=','
                    a=self.GetParVarian((self.outdir,arr))            
                    nz*=len(a)
                    arrStr+=arr+'('+str(len(a))+')'
                    if(arr!='phase' and arr!='phase2'):
                        pseudo=True
                
                self.labels['Array']=arrStr

            
            """
            arrs=self.GetParVarian((self.outdir,'array'))
            
            print('arrays:',arrs)
            nz=1
            if(len(arrs)!=0 and arrs[0]!='""'):
            
            arrs=arrs[0].split('"')[1].split(',')
            arrStr=''
            for i,arr  in enumerate(arrs):
            if(i!=0):
            arrStr+=','
            a=self.GetParVarian((self.outdir,arr))            
            nz*=len(a)
            arrStr+=arr+'('+str(len(a))+')'
            if(arr!='phase' and arr!='phase2'):
            pseudo=True
            
            self.labels[cnt]='array:'+arrStr;cnt+=1
            """
            #self.labels[7]='fids: %i' % (nz);cnt+=1
            #print(nz,nt,ni,ni2,(d1+at))
        self.labels['Measurement time']='%s' % self.FormatTime(nz*nt*incr*(d1+at))

        if(pseudo):
            dim=str(dim)+'p'

        if(self.dim==dim):
            self.labels['Dimensions']='%s' % str(dim)
        else:
            self.labels['Dimensions set']= '%s' % str(self.dim)
            self.labels['Dimensions expected']='%s' % str(dim)

                
        #elif(self.tp=='var'):

        #print(self.GetParVarian((self.outdir,'nt'))[0])

        #sfrq=float(self.GetParVarian((self.outdir,'sfrq'))[0])
        #sw=float(self.GetParVarian((self.outdir,'sw'))[0])
        #at=float(self.GetParVarian((self.outdir,'at'))[0])

        #np=int(self.GetParVarian((self.outdir,'np'))[0])

            
        
        #self.labels[cnt]='sfrq: %.0f MHz' % sfrq;cnt+=1
        #self.labels['NT']='NT: %i' % nt

        #at=self.np2/self.sw
        #self.labels[cnt]='Direct: %s aq: %.2f ms np: %i sw: %.2f Hz swp: %.2f ppm car: %.2f ppm' % (self.n1,self.np2/self.sw*1000,self.np2,self.sw,self.sw/self.sfrq,self.waterppm);cnt+=1
        #self.labels[cnt]='Direct: at: %.2f ms np: %i sw: %.2f Hz swp: %f ppm' % (at*1000,np,sw,sw/sfrq);cnt+=1
        #    self.texTab.append('%s & %.2f & %.2f & %i & %.2f & %.2f & %.2f \n' % (self.n1,self.sfrq,at*1000,self.np2,self.sw,self.sw/self.sfrq,self.waterppm))
        """

            #self.labels[cnt]='at: %f s' % at;cnt+=1
            #self.labels[cnt]='NP: %i' % np;cnt+=1
            #self.labels[cnt]='sw:'+'%f ppm' % (sw/sfrq);cnt+=1

            nis=self.GetParVarian((self.outdir,'ni'))
            if(len(nis)>0):
                ni=int(self.GetParVarian((self.outdir,'ni'))[0])
                if(ni>1):
                    #dfrq=int(self.GetParVarian((self.outdir,'dfrq'))[0])
                    sw1=float(self.GetParVarian((self.outdir,'sw1'))[0])
                    self.labels[cnt]='Indirect: at: %.2f ms ni: %i sw: %.2f Hz' % (ni/sw1*1000,ni,sw1);cnt+=1
                else:
                    ni=1
            else:
                ni=1

            ni2s=self.GetParVarian((self.outdir,'ni2'))
            if(len(ni2s)>0):
                ni2=int(self.GetParVarian((self.outdir,'ni2'))[0])
                if(ni2>1):
                    #dfrq=int(self.GetParVarian((self.outdir,'dfrq'))[0])
                    sw2=float(self.GetParVarian((self.outdir,'sw2'))[0])
                    self.labels[cnt]='Indirect2: at: %.2f ms ni: %i sw: %.2f Hz' % (ni2/sw2*1000,ni2,sw2);cnt+=1
                else:
                    ni2=1
            else:
                ni2=1
        """
            

        #self.labels[7]='fids: %i' % (nz);cnt+=1
        """
        self.labels[cnt]='time: %s' % self.FormatTime(nz*nt*ni*(d1+at));cnt+=1


            if(ni==1 and ni2==1):
                dim=1
            elif(ni2>1 and ni>1):
                dim=3
            elif(ni2>2 or ni>1):
                dim=2
            
            if(pseudo):
                dim=str(dim)+'p'

            if(self.dim==dim):
                self.labels[cnt]='Dimensions: %s' % str(dim);cnt+=1
            else:
                self.labels[cnt]='Dimensions set: %s' % str(self.dim);cnt+=1
                self.labels[cnt]='Dimensions expect: %s' % str(dim);cnt+=1


        else:
            print('Spectrometer type not recognised')
        """

            
    def FormatTime(self,sec):
        if(sec<60):
            return '%i s' % (sec)
        else:
            return '%.2f min' % (sec/60.)
        return
        

#parse a file to a list.
def readfile(infile):
    peak=[]
    peakfile=open(infile,'r')
    for line in peakfile.readlines():
        linetosave=line.split()
        peak.append(linetosave)
    peakfile.close()
    return peak

        
def GetParBrukFile(filey):
    pass
    if not os.path.exists(filey):
        pass
        return []
    inny=open(filey)
    pars=[]
    for line in inny.readlines():
        test=line.split()
        if(len(test)>0):
            pars.append(test)
    pass
    return pars


#analyse either acqu and acqu2
def GetParBruk(infile,argv,verb='n',):
    args=[]
    full_path=infile
    pass
    if not os.path.exists(full_path):
        pass
        return (0,)
    procpar=readfile(full_path)
    for i in range(len(argv)-1):
        param=argv[i+1]
        tick=0
        for j in range(len(procpar)):
            #print procpar[j]
            test=procpar[j][0].split('##$')

            if(len(test)>1):
                test2=test[1].split('=')[0]
                if(test2==param):
                    if(verb=='y'):
                        pass
                    args.append(procpar[j][1])
                    tick=1
            else:
                #we have a line of zeros
                #is the previous line what we're after?
                test=procpar[j-1][0].split('##$')
                if(len(test)>1):
                    test2=test[1].split('=')[0]
                    for i in range(100):
                        parT=test2+str(i)
                        if(parT==param):
                            if(len(param.split(test2))>1):
                                if(verb=='y'):
                                    pass
                                #parameters are in rows in j,j+1,j+2...
                                go=0
                                cnt=0
                                while(go==0):
                                    if(i<len(procpar[j+cnt])):
                                        val=procpar[j+cnt][i]
                                        go=1
                                    else:
                                        i-=len(procpar[j+cnt])
                                        cnt+=1

                                args.append(val)
                                if(verb=='y'):
                                    pass
                                tick=1
            #sys.exit(100)
        if(tick==0):
            if(verb=='y'):
                pass
            return 0,
        else:
            return args



import numpy as np
import nmrglue as ng
from scipy.optimize import leastsq
class estNoise():
    def __init__(self,binFac=1000,fit=False,write=False,folder=''):
        self.binFac=binFac
        self.folder=folder
        self.read()
        if(fit):

            self.binNo=self.binFac
            bins=np.linspace(-2*self.dstd,2*self.dstd,self.binNo)
            self.hist,self.edges=np.histogram(self.fat,bins=bins)
            self.edges=(self.edges[:-1]+self.edges[1:])*0.5
            
            argh=np.argmax(self.hist)                
            self.noise=np.fabs(self.edges[argh])



            #print('noise         :',self.noise)
            self.x=np.array(self.edges)
            self.ydat=np.array(self.hist)
            self.fit()
        if(write):
            self.write()
            
    def Parse(self,par):
        inny=open('deconParFile')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0 and test[0]==par):
                return test[2]
        return None
            
    def read(self):
        dic,dat=ng.pipe.read(self.Parse('infile'))
        self.dmax=np.max(dat)
        self.dstd=np.std(dat)
        #fat=np.ndarray.flatten(dat)
        self.fat=np.ravel(dat)
        mask=(self.fat>-2*self.dstd)*(self.fat<2*self.dstd)
        fot=self.fat[mask]
        self.noise=np.std(fot)
        self.thresh=float(self.Parse('thresh'))
        self.noiseStd=self.noise

        pass
        pass
        pass
        pass
        pass
                    
        
    def guess(self):
        argy=np.argmax(self.ydat)
        self.x0=self.x[argy]
        self.A=self.ydat[argy]
        ytest=self.ydat[argy:]
        for i in range(len(ytest)):
            if(ytest[i]<self.A/2.):
                self.sig=(self.x[argy+i]-self.x[argy])
                break
        #print('guess:',self.A,self.x0,self.sig)
    def pack(self):
        x=[]
        x.append(self.A)
        x.append(self.x0)
        x.append(self.sig)
        return x
    def unpack(self,x):
        cnt=0
        self.A=x[cnt];cnt+=1
        self.x0=x[cnt];cnt+=1
        self.sig=x[cnt];cnt+=1
    def ycalc(self):
        self.y=self.A*np.exp(-(self.x-self.x0)**2./(2.*self.sig**2.))
    def chi(self,x):
        self.unpack(x)
        self.ycalc()
        return self.ydat-self.y
    def fit(self):
        self.guess()
        x0=leastsq(self.chi,self.pack())
        #print (x0)
        self.noise=self.sig
        pass
        pass
        pass
    def write(self):
        self.outfile=os.path.join(self.folder,'noise.out')
        out=open(self.outfile,'w')
        for i in range(len(self.x)):
            out.write('%e\t%e\t%e\n' % (self.x[i],self.ydat[i],self.y[i]))
        out.close()
        self.makeplot()

    def makeplot(self):
        pass
        self.figfile=self.outfile.replace('.out','.pdf')
        gnu=open('gnu.gp','w')
        gnu.write('set term pdf enh color solid\n')
        gnu.write('set output \'%s\'\n' % self.figfile)
        gnu.write('set size square\n')
        import math
        self.ord=math.floor(math.log10(np.max(self.x)))

        gnu.write('set xlabel \'x10^{%i}\'\n' % self.ord)
        gnu.write('set ylabel \'Count\'\n')
        #gnu.write('set xtics rotate\n')
        #gnu.write('set format x "%.0sx10^{%T}"\n')
        gnu.write('set label sprintf(\'MaxSignal  : %.2e\') at graph 0.1,0.9\n'  % self.dmax)
        gnu.write('set label sprintf(\'Noise(std) : %.2e\') at graph 0.1,0.85\n' % self.noiseStd)
        gnu.write('set label sprintf(\'S/N(std)   : %.2f\') at graph 0.1,0.8\n'  % (self.dmax/self.noiseStd))
        gnu.write('set label sprintf(\'Noise(fit) : %.2e\') at graph 0.1,0.75\n' % self.noise)
        gnu.write('set label sprintf(\'S/N(fit)   : %.2f\') at graph 0.1,0.7\n'  % (self.dmax/self.noise)) 
        gnu.write('plot \'%s\' u ($1/%i):2 ti \'histogram\',\'\' u ($1/%i):3 ti \'fit\' w li\n' % (self.outfile,10**self.ord,10**self.ord))
        gnu.close()
        import os
        os.system('gnuplot gnu.gp')
        os.system('rm gnu.gp')
        os.system('rm noise.out')




        
