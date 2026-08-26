"""Conservative status-bar hover helper for secondary plotting panels.

The helper is deliberately non-invasive: it binds only leaf interactive controls
(or explicitly paired StaticText labels), always propagates wx events, and never
assumes a particular parent/frame hierarchy.
"""
import wx

# Keep this list to stable wx core controls used by these frames.  Do not bind
# panels, StaticBoxes, canvases, toolbars, sizers, or other containers.
_INTERACTIVE = tuple(
    cls for cls in (
        getattr(wx, 'TextCtrl', None), getattr(wx, 'Button', None),
        getattr(wx, 'ToggleButton', None), getattr(wx, 'SpinButton', None),
        getattr(wx, 'Slider', None), getattr(wx, 'CheckBox', None),
        getattr(wx, 'ComboBox', None), getattr(wx, 'Choice', None),
        getattr(wx, 'RadioButton', None), getattr(wx, 'SpinCtrl', None),
        getattr(wx, 'ListCtrl', None),
    ) if cls is not None
)
_STATIC_TEXT = getattr(wx, 'StaticText', None)


def _find_status_bar(owner):
    """Find an existing status bar without depending on a specific frame type."""
    current = owner
    seen = set()
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        for attr in ('statusbar', 'status_bar', 'statusBar'):
            try:
                bar = getattr(current, attr, None)
                if bar is not None and hasattr(bar, 'SetStatusText'):
                    return current, bar
            except Exception:
                pass
        try:
            getter = getattr(current, 'GetStatusBar', None)
            bar = getter() if callable(getter) else None
            if bar is not None and hasattr(bar, 'SetStatusText'):
                return current, bar
        except Exception:
            pass
        try:
            current = current.GetParent()
        except Exception:
            current = None
    return None, None


def bind_status_help(owner, widget, text, allow_label=False):
    """Bind status help to a leaf control, never altering its normal events."""
    if widget is None or not text:
        return
    try:
        if not isinstance(widget, wx.Window):
            return
        is_label = bool(allow_label and _STATIC_TEXT is not None and isinstance(widget, _STATIC_TEXT))
        if not isinstance(widget, _INTERACTIVE) and not is_label:
            return
    except Exception:
        return

    if getattr(widget, '_secondary_status_help_bound', False):
        widget._secondary_status_help_text = str(text)
        return
    widget._secondary_status_help_bound = True
    widget._secondary_status_help_text = str(text)

    def enter(evt):
        try:
            frame, bar = _find_status_bar(owner)
            if bar is not None:
                # Save the pre-hover message only when starting a new hover.
                # Widget identity (not message text) prevents paired labels and
                # controls with identical help from restoring each other badly.
                active = getattr(frame, '_secondary_hover_widget', None)
                if active is None:
                    try:
                        frame._secondary_hover_previous = bar.GetStatusText()
                    except Exception:
                        frame._secondary_hover_previous = ''
                frame._secondary_hover_widget = widget
                frame._secondary_hover_text = widget._secondary_status_help_text
                bar.SetStatusText(widget._secondary_status_help_text)
        except Exception:
            # Hover help must never be capable of breaking a scientific GUI.
            pass
        finally:
            evt.Skip()

    def leave(evt):
        try:
            frame, bar = _find_status_bar(owner)
            if bar is not None and getattr(frame, '_secondary_hover_widget', None) is widget:
                try:
                    if bar.GetStatusText() == widget._secondary_status_help_text:
                        bar.SetStatusText(getattr(frame, '_secondary_hover_previous', '') or '')
                finally:
                    frame._secondary_hover_widget = None
                    frame._secondary_hover_text = None
                    frame._secondary_hover_previous = ''
        except Exception:
            pass
        finally:
            evt.Skip()

    try:
        widget.Bind(wx.EVT_ENTER_WINDOW, enter)
        widget.Bind(wx.EVT_LEAVE_WINDOW, leave)
    except Exception:
        # A control that cannot accept these bindings simply gets no hover help.
        pass


def bind_map(owner, mapping):
    """Bind a sequence of (widget, text[, allow_label]) entries defensively."""
    for item in mapping:
        try:
            if len(item) == 2:
                widget, text = item
                allow_label = False
            else:
                widget, text, allow_label = item
            bind_status_help(owner, widget, text, allow_label=allow_label)
        except Exception:
            # A missing/unsupported optional control must not stop frame setup.
            continue
