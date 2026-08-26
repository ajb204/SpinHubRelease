"""Immutable plain-value state used by NMRPipe script generation.

This module is the single boundary where mutable wx controls are sampled.
Everything downstream receives ordinary immutable Python values.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType




# Processing controls are part of the processing/script contract, not ProcessFrame.
# Keeping the schema beside ProcessingScriptState prevents the top-level process
# window from having to know which wx controls the NMRPipe renderer consumes.
PROCESSING_CONTROL_NAMES = (
    'cb_baseLin', 'cb_basepol', 'cb_basePol', 'cb_baseSol', 'maxIterBox', 'mddMethodBox', 'mddAlgorithmBox', 'mddIterBox', 'mddVEBox',
    'p0', 'p1', 'cb_ft0', 'windowBox0', 'firstPoint0', 'win2Val0', 'win3Val0',
    'windowOp1', 'windowOp2', 'firstPointFactor', 'f1180Lab', 'lpLab', 'polyLab',
    'p0Lab', 'p1Lab', 'flipLab', 'windowLab', 'baseLab', 'lab0', 'lab1', 'lab2', 'lab3',
    'cb_f1180', 'cb_lp1', 'cb_basepol1', 'p0_1', 'p1_1', 'cb_ft1', 'windowBox1', 'firstPoint1', 'win2Val1', 'win3Val1',
    'cb_f2180', 'cb_lp2', 'cb_basepol2', 'p0_2', 'p1_2', 'cb_ft2', 'windowBox2', 'firstPoint2', 'win2Val2', 'win3Val2',
    'cb_f3180', 'cb_lp3', 'cb_basepol3', 'p0_3', 'p1_3', 'cb_ft3', 'windowBox3', 'firstPoint3', 'win2Val3', 'win3Val3',
)

@dataclass(frozen=True)
class ProcessingScriptState:
    """Plain-value snapshot of processing controls used for script rendering.

    This is the explicit data boundary between wx widgets/shared GUI state and
    NMRPipe generation.  Script builders never need to inspect live widgets to
    decide which value is current.
    """
    values: object

    @classmethod
    def capture_current(cls, processing, project_state=None):
        """Capture the complete script parameter contract from current GUI state.

        An open ProcessingFrame wins for every control it owns. Missing controls
        fall back to ProjectState.gui_settings. This is the only place where the
        script layer crosses from mutable GUI state into an immutable snapshot.
        """
        live = getattr(project_state, 'gui_settings', {}) if project_state is not None else {}
        return cls.capture(processing, live, PROCESSING_CONTROL_NAMES)

    @classmethod
    def capture(cls, processing, live, control_names=PROCESSING_CONTROL_NAMES):
        values = {}
        live = live or {}
        for name in control_names:
            control = getattr(processing, name, None) if processing is not None else None
            if control is not None:
                values[name] = cls._read_control_value(name, control)
            elif name in live:
                value = live[name]
                values[name] = value
        return cls(MappingProxyType(values))

    @staticmethod
    def _read_control_value(name, control):
        """Extract the script-semantic value from a wx-like control.

        wx.TextCtrl.GetSelection() returns the selected character range as a
        ``(start, end)`` tuple.  The old generic probe tried GetSelection before
        GetValue, so text processing parameters (notably the 3D phase/window
        fields) were captured as tuples and later failed at ``float(...)``.

        Only FT combo boxes need their numeric selection index. Window combo
        boxes need their displayed string (GM/SP/EM), checkboxes need booleans,
        and ordinary text controls need their text value.
        """
        if name.startswith('cb_ft'):
            method = getattr(control, 'GetSelection', None)
            if callable(method):
                try:
                    return method()
                except Exception:
                    pass

        method = getattr(control, 'IsChecked', None)
        if callable(method):
            try:
                return bool(method())
            except Exception:
                pass

        method = getattr(control, 'GetValue', None)
        if callable(method):
            try:
                return method()
            except Exception:
                pass

        method = getattr(control, 'GetLabel', None)
        if callable(method):
            try:
                return method()
            except Exception:
                pass
        return control

    def control(self, name, default=None):
        """Compatibility name for a plain value; never returns a wx object."""
        return self.values.get(name, default)

    def value(self, name, default=None):
        return self.values.get(name, default)

    def checked(self, name, default=False):
        return bool(self.values.get(name, default))

    def selection(self, name, default=0):
        value = self.values.get(name, default)
        try:
            return int(value)
        except (TypeError, ValueError):
            return int(default)
