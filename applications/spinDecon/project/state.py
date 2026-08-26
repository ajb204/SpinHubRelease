"""Shared project state and canonical path resolution for the decon GUI."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from spinDecon.domain.topology import DatasetTopology


@dataclass
class ProjectState:
    session_file: str = ""
    working_dir: str = ""
    parameter_file: str = ""
    raw_path: str = "./raw"
    spec_path: str = "./spec"
    input_file: str = ""
    reference_peak_file: str = ""
    full_peak_file: str = ""
    # Canonical contract: ``dimension`` is retained as a compatibility field
    # name, but from Stage 2 onward it ALWAYS means the number of SPECTRAL
    # dimensions.  The pseudo axis, when present, is one additional physical
    # data dimension.  New code should prefer ``spectral_dimensions``.
    dimension: int = 0
    pseudo_axis: bool = False
    fid_selection: int = 1
    sym_mode: bool = False
    decon_bore: bool = False
    loaded: bool = False
    dirty: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)
    # Live, unsaved GUI parameter values. Disk is a persistence boundary, not
    # the message bus between windows.
    gui_settings: dict[str, Any] = field(default_factory=dict)
    # Projection phasing is deliberately transient until Re-process promotes it.
    # Keeping it outside gui_settings prevents an ordinary Save from committing
    # a slider preview as a processing phase.
    projection_phase_preview: dict[str, dict[str, float]] = field(default_factory=dict)

    def update_gui_settings(self, updates: dict[str, Any]) -> None:
        """Merge user edits into the authoritative live settings."""
        self.gui_settings.update({str(k): v for k, v in updates.items()})
        self.dirty = True

    def seed_gui_settings(self, updates: dict[str, Any]) -> None:
        """Hydrate missing live settings without overwriting newer GUI edits.

        This is used after a frame loads persisted values.  Once a key exists in
        ``gui_settings`` the live value wins over disk until an explicit project
        reload creates/reseeds state.
        """
        for key, value in updates.items():
            self.gui_settings.setdefault(str(key), value)

    def hydrate_gui_settings(self, updates: dict[str, Any], *, overwrite: bool = False) -> None:
        """Hydrate persisted session values without marking them as user edits.

        ``overwrite=True`` is reserved for a new Process-session boundary, where
        the saved system file is authoritative before any Process child can edit
        shared state.  Normal child-window hydration should continue to use
        :meth:`seed_gui_settings` so newer live edits always win.
        """
        for key, value in updates.items():
            key = str(key)
            if overwrite or key not in self.gui_settings:
                self.gui_settings[key] = value

    def resolved_gui_value(self, key: str, fallback: Any = None) -> Any:
        """Return the live value when present, otherwise *fallback*."""
        return self.gui_settings[key] if key in self.gui_settings else fallback

    def gui_value(self, key: str, default: Any = None) -> Any:
        return self.gui_settings.get(key, default)

    def projection_phase(self, label: str, *, p0: float = 0.0, p1: float = 0.0) -> dict[str, float]:
        """Return the transient projection phase for *label*.

        Projection slider values are preview state, not persistent processing
        parameters.  They only enter ``gui_settings`` through explicit
        promotion by Re-process.
        """
        entry = self.projection_phase_preview.get(str(label), {})
        try:
            out0 = round(float(entry.get('p0', p0)), 2)
        except Exception:
            out0 = round(float(p0), 2)
        try:
            out1 = round(float(entry.get('p1', p1)), 2)
        except Exception:
            out1 = round(float(p1), 2)
        return {'p0': out0, 'p1': out1}

    def update_projection_phase(self, label: str, *, p0: Any = None, p1: Any = None) -> None:
        """Update projection preview state without marking the project dirty."""
        entry = self.projection_phase_preview.setdefault(str(label), {})
        if p0 is not None:
            entry['p0'] = round(float(p0), 2)
        if p1 is not None:
            entry['p1'] = round(float(p1), 2)
        # Mirror for older code that still inspects metadata.
        self.metadata['projection_phasing'] = self.projection_phase_preview

    def clear_projection_phase_preview(self) -> None:
        self.projection_phase_preview.clear()
        self.metadata['projection_phasing'] = self.projection_phase_preview

    def promote_projection_phase(self, updates: dict[str, Any]) -> dict[str, Any]:
        """Promote accepted projection phases into authoritative live settings."""
        promoted = {str(k): v for k, v in updates.items()}
        self.update_gui_settings(promoted)
        return promoted


    @property
    def spectral_dimensions(self) -> int:
        """Canonical number of spectral (frequency/chemical-shift) axes."""
        return int(self.dimension or 0)

    @spectral_dimensions.setter
    def spectral_dimensions(self, value: int) -> None:
        self.dimension = int(value)

    @property
    def physical_dimensions(self) -> int:
        """Canonical physical array dimensionality implied by project state."""
        if self.spectral_dimensions < 1:
            return 0
        return self.spectral_dimensions + int(bool(self.pseudo_axis))

    def topology(self, **kwargs) -> DatasetTopology:
        """Return the canonical topology represented by this project state."""
        if self.spectral_dimensions < 1:
            raise ValueError("project spectral dimensionality is not set")
        # Loaded projects retain physical-axis identity in metadata so every
        # consumer sees the same axis ordering.  Callers may still override
        # these values explicitly (useful for new/unloaded projects and tests).
        if self.pseudo_axis and "pseudo_physical_index" not in kwargs:
            kwargs["pseudo_physical_index"] = int(self.metadata.get("pseudo_physical_index", 0))
        if "spectral_labels" not in kwargs:
            labels = self.metadata.get("spectral_axis_labels")
            if labels and len(labels) == self.spectral_dimensions:
                kwargs["spectral_labels"] = tuple(labels)
        if self.pseudo_axis and "pseudo_label" not in kwargs:
            kwargs["pseudo_label"] = str(self.metadata.get("pseudo_axis_label", ""))
        return DatasetTopology.from_counts(
            self.spectral_dimensions, bool(self.pseudo_axis), **kwargs
        )

    def canonicalize_loaded_dimensions(
        self,
        physical_ndim: int,
        axis_labels: Iterable[str] = (),
        *,
        real_axis_labels: Iterable[str] = (),
    ) -> bool:
        """Normalize an ambiguous legacy loaded project at the load boundary.

        Historical projects could persist ``dim`` as the physical NMRPipe
        dimension for pseudo data (notably dim=3,pseudo=1 for 2D+pseudo).
        Canonical state always stores the spectral count.  We only rewrite an
        ambiguous value when the loaded array has exactly ``dim`` physical
        axes and exactly one supplied axis label is recognised as real/pseudo.

        Returns True when a legacy physical count was migrated.
        """
        physical = int(physical_ndim)
        dim = self.spectral_dimensions
        pseudo = bool(self.pseudo_axis)
        labels = tuple(str(x) for x in axis_labels)[:physical]
        real = frozenset(str(x).strip().lower() for x in real_axis_labels)
        real_indices = [i for i, label in enumerate(labels) if label.strip().lower() in real]
        real_count = len(real_indices)

        pass

        migrated = False
        legacy_unresolved = self.metadata.get('dimension_contract') == 'legacy_unresolved'
        # A state read from a legacy parameter file explicitly records that its
        # dimension contract is unresolved.  For pseudo data, the historical
        # value was the physical count, so a declared dim equal to ndarray.ndim
        # can be migrated without relying on a particular pseudo-axis label.
        # For ad-hoc/in-memory states we retain the stricter label evidence to
        # avoid guessing an ambiguous scientific topology.
        if (pseudo and dim == physical and physical >= 2 and
                (legacy_unresolved or real_count == 1)):
            self.dimension = physical - 1
            migrated = True
            self.metadata['legacy_dimension_migrated_from'] = dim
            self.metadata['dimension_contract'] = 'spectral'
        elif dim:
            expected = dim + int(pseudo)
            if physical != expected:
                raise ValueError(
                    f"loaded data has {physical} physical dimensions but project "
                    f"topology requires {expected} ({dim} spectral + {int(pseudo)} pseudo)"
                )
            self.metadata['dimension_contract'] = 'spectral'

        # Persist physical-axis identity at the one authoritative load boundary.
        # This prevents viewers/processors from falling back to an assumed
        # pseudo-axis position after the project has been canonicalised.
        if pseudo and real_count == 1:
            pseudo_index = real_indices[0]
            self.metadata['pseudo_physical_index'] = pseudo_index
            self.metadata['pseudo_axis_label'] = labels[pseudo_index]
            self.metadata['spectral_axis_labels'] = tuple(
                label for i, label in enumerate(labels) if i != pseudo_index
            )
        elif not pseudo and labels:
            self.metadata['spectral_axis_labels'] = tuple(labels)
        pass
        return migrated

    def _spec_relative(self, value: str | Path) -> str:
        """Store spectrum-associated files as paths relative to SpecPath.

        Subdirectories (for example ``projections/1H.15N.2D.list``) are
        significant and must never be collapsed to a basename.
        """
        value = str(value or '').strip()
        if not value:
            return ''
        p = Path(value)
        if p.is_absolute():
            try:
                p = p.resolve().relative_to(Path(self.spec_dir()).resolve())
            except (ValueError, OSError):
                return str(p)  # compatibility: resolver will keep it absolute
        norm = Path(str(p).replace('\\', '/'))
        if any(part == '..' for part in norm.parts):
            raise ValueError('Spectrum-associated paths may not escape SpecPath')

        # Older parameter files (and some legacy GUI paths) may already include
        # the project-relative SpecPath prefix, e.g. ``./spec/test.ft3``.  The
        # canonical value stored in the NMR text boxes is *relative to*
        # SpecPath, so remove that prefix exactly once while preserving genuine
        # descendants such as ``projections/1H.15N.2D.list``.
        spec = Path(str(self.spec_path or './spec').replace('\\', '/'))
        if not spec.is_absolute():
            spec_parts = tuple(part for part in spec.parts if part not in ('', '.'))
            norm_parts = tuple(part for part in norm.parts if part not in ('', '.'))
            if spec_parts and norm_parts[:len(spec_parts)] == spec_parts:
                norm = Path(*norm_parts[len(spec_parts):])

        return norm.as_posix().lstrip('./')

    def _project_resolve(self, value: str | Path) -> Path:
        p = Path(str(value or '').strip() or '.')
        if p.is_absolute():
            return p
        base = Path(self.working_dir or '.')
        return base / p

    def working_path(self) -> str:
        return str(Path(self.working_dir or '.'))

    def raw_dir(self) -> str:
        return str(self._project_resolve(self.raw_path or './raw'))

    def spec_dir(self) -> str:
        return str(self._project_resolve(self.spec_path or './spec'))

    def resolve_spec_file(self, value: str | Path) -> str:
        value = str(value or '').strip()
        if not value:
            return ''
        p = Path(value)
        # Absolute paths are accepted only as an operational compatibility
        # boundary; persisted/displayed values remain names under SpecPath.
        if p.is_absolute():
            return str(p)
        # Canonicalise here as well as at persistence/UI boundaries.  This
        # prevents a legacy value such as ``./spec/test.ft3`` from resolving
        # to ``<WorkingDir>/spec/spec/test.ft3``.
        rel = self._spec_relative(value)
        return str(Path(self.spec_dir()) / rel)

    def spectrum_path(self) -> str:
        return self.resolve_spec_file(self.input_file)

    def reference_peak_path(self) -> str:
        return self.resolve_spec_file(self.reference_peak_file)

    def full_peak_path(self) -> str:
        return self.resolve_spec_file(self.full_peak_file)

    def housekeeping_paths(self) -> dict[str, str]:
        """Return the canonical runtime paths mirrored by UniDecNMR housekeeping.

        ProjectState remains authoritative.  Housekeeping/UI code may display
        the persisted project-relative values, but runtime consumers that need
        concrete locations should use this mapping rather than cwd.
        """
        return {
            "workingDirectory": self.working_path(),
            "FIDpath": self.raw_dir(),
            "specPath": self.spec_dir(),
            "infile": self.spectrum_path(),
        }

    # Compatibility alias for legacy code/parameter naming.
    @property
    def peak_file(self) -> str:
        return self.reference_peak_file

    @peak_file.setter
    def peak_file(self, value: str | Path) -> None:
        self.reference_peak_file = self._spec_relative(value)

    def set_session_file(self, path: str | Path) -> None:
        self.session_file = str(path)
        if not self.working_dir:
            self.working_dir = str(Path(path).parent)

    def set_parameter_file(self, path: str | Path) -> None:
        self.parameter_file = str(path)
        self.set_session_file(path)


    @classmethod
    def from_parameter_file(cls, path: str | Path) -> "ProjectState":
        """Load the path-bearing core of a deconParFile.

        Resource existence is deliberately not validated here. A project with
        an offline raw acquisition or missing spectrum remains loadable and can
        be inspected by SpinHub. Relative paths are interpreted relative to the
        directory containing the deconParFile via ``working_dir``.
        """
        from .parameter_store import parse_int, parse_value

        state = cls()
        state.set_parameter_file(Path(path))

        def clean(value, default=''):
            value = str(value if value is not None else default).strip()
            return default if value in ('', '0') else value

        state.raw_path = clean(parse_value(path, 'fiddir', './raw'), './raw')
        state.spec_path = clean(parse_value(path, 'specPath', './spec'), './spec')
        state.dimension = parse_int(path, 'dim', parse_int(path, 'dimension', 0))
        state.pseudo_axis = bool(parse_int(path, 'pseudo', 0))
        # Missing legacy path fields default to the processing pipeline's
        # dimensional output name.  These values are stored relative to
        # SpecPath; e.g. 2 spectral + pseudo -> test.ft2.
        spectral_dim = max(1, min(4, int(state.dimension or 2)))
        default_spectrum = 'test.ft' if spectral_dim == 1 else 'test.ft%d' % spectral_dim
        state.input_file = state._spec_relative(clean(parse_value(path, 'infile', default_spectrum), default_spectrum))
        state.reference_peak_file = state._spec_relative(clean(parse_value(path, 'peakfile', default_spectrum + '.list'), default_spectrum + '.list'))
        state.full_peak_file = state._spec_relative(clean(parse_value(path, 'fullPeakFile', '')))
        state.metadata['dimension_contract'] = 'legacy_unresolved' if state.pseudo_axis else 'spectral'
        state.loaded = True
        state.dirty = False
        return state

    def sync_from_values(self, *, working_dir=None, raw_path=None, input_file=None,
                         peak_file=None, reference_peak_file=None, full_peak_file=None,
                         spec_path=None, dimension=None, spectral_dimensions=None,
                         pseudo_axis=None, fid_selection=None, sym_mode=None,
                         decon_bore=None) -> None:
        if working_dir is not None: self.working_dir = str(working_dir)
        if raw_path is not None: self.raw_path = str(raw_path)
        if spec_path is not None: self.spec_path = str(spec_path)
        if input_file is not None: self.input_file = self._spec_relative(input_file)
        ref = reference_peak_file if reference_peak_file is not None else peak_file
        if ref is not None: self.reference_peak_file = self._spec_relative(ref)
        if full_peak_file is not None: self.full_peak_file = self._spec_relative(full_peak_file)
        if spectral_dimensions is not None:
            self.spectral_dimensions = int(spectral_dimensions)
        elif dimension is not None:
            # Compatibility input; ``dimension`` is canonical spectral count.
            self.spectral_dimensions = int(dimension)
        if pseudo_axis is not None: self.pseudo_axis = bool(pseudo_axis)
        if fid_selection is not None: self.fid_selection = int(fid_selection)
        if sym_mode is not None: self.sym_mode = bool(sym_mode)
        if decon_bore is not None: self.decon_bore = bool(decon_bore)
        self.dirty = True

    @property
    def deconParFile(self): return self.parameter_file
    @deconParFile.setter
    def deconParFile(self, path): self.set_parameter_file(path)
    @property
    def parameter_path(self): return self.parameter_file
    @parameter_path.setter
    def parameter_path(self, path): self.set_parameter_file(path)

    def as_dict(self) -> dict[str, Any]:
        return {k: getattr(self, k) for k in (
            'session_file','working_dir','parameter_file','raw_path','spec_path',
            'input_file','reference_peak_file','full_peak_file','dimension',
            'pseudo_axis','fid_selection','sym_mode','decon_bore','loaded','dirty'
        )} | {
            'spectral_dimensions': self.spectral_dimensions,
            'physical_dimensions': self.physical_dimensions,
            'peak_file': self.reference_peak_file,
            'metadata': dict(self.metadata),
        }


def coerce_project_state(state=None, *, session_file='', parameter_file='', spec_path=''):
    if state is None:
        state = ProjectState(session_file=session_file, parameter_file=parameter_file)
    if session_file and not state.session_file: state.session_file = session_file
    if parameter_file: state.parameter_file = parameter_file
    if spec_path: state.spec_path = spec_path
    return state
