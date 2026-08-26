"""Reusable wxPython controls shared by decon GUI frames."""
import wx


class PersistentStateButton(wx.Control):
    """Owner-drawn button with a pressed state that persists across windows."""
    def __init__(self, parent, id=wx.ID_ANY, label='', size=wx.DefaultSize):
        super().__init__(parent, id=id, size=size, style=wx.BORDER_NONE)
        self._label = label
        self._active = False
        self._mouse_down = False
        self._active_colour = wx.Colour(80, 145, 220)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetMinSize(size)
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_LEFT_DOWN, self._on_left_down)
        self.Bind(wx.EVT_LEFT_UP, self._on_left_up)
        self.Bind(wx.EVT_MOUSE_CAPTURE_LOST, self._on_capture_lost)
        self.Bind(wx.EVT_KEY_DOWN, self._on_key_down)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda event: None)

    def SetActive(self, active):
        active = bool(active)
        if active != self._active:
            self._active = active
            self.Refresh()

    def IsActive(self):
        return self._active

    def SetLabel(self, label):
        self._label = label
        self.Refresh()

    def GetLabel(self):
        return self._label

    def _on_left_down(self, event):
        if not self.IsEnabled():
            return
        self._mouse_down = True
        self.SetFocus()
        if not self.HasCapture():
            self.CaptureMouse()
        self.Refresh()

    def _on_left_up(self, event):
        if self.HasCapture():
            self.ReleaseMouse()
        was_down = self._mouse_down
        self._mouse_down = False
        self.Refresh()
        if was_down and self.IsEnabled() and self.GetClientRect().Contains(event.GetPosition()):
            command = wx.CommandEvent(wx.EVT_BUTTON.typeId, self.GetId())
            command.SetEventObject(self)
            self.GetEventHandler().ProcessEvent(command)

    def _on_capture_lost(self, event):
        self._mouse_down = False
        self.Refresh()

    def _on_key_down(self, event):
        if self.IsEnabled() and event.GetKeyCode() in (wx.WXK_SPACE, wx.WXK_RETURN, wx.WXK_NUMPAD_ENTER):
            command = wx.CommandEvent(wx.EVT_BUTTON.typeId, self.GetId())
            command.SetEventObject(self)
            self.GetEventHandler().ProcessEvent(command)
        else:
            event.Skip()

    def _on_paint(self, event):
        dc = wx.AutoBufferedPaintDC(self)
        rect = self.GetClientRect()
        renderer = wx.RendererNative.Get()

        # Let wx/the platform draw the ordinary (unlatched) button.  This is
        # deliberately different from hand-emulating a native button: themes,
        # DPI settings and platform-specific button geometry are then inherited
        # automatically from the rest of the application.
        if not self._active:
            flags = 0
            if self._mouse_down:
                flags |= wx.CONTROL_PRESSED
            if not self.IsEnabled():
                flags |= wx.CONTROL_DISABLED
            if self.HasFocus():
                flags |= wx.CONTROL_FOCUSED
            renderer.DrawPushButton(self, dc, rect, flags)
            text = (wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT)
                    if not self.IsEnabled()
                    else wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT))
            pressed = self._mouse_down
        else:
            # The latched state keeps the same system-button vocabulary, but
            # reverses the edge perspective and uses a coloured face so that
            # the persistent state remains obvious even after focus moves to a
            # different window.
            face = self._active_colour
            if not self.IsEnabled():
                face = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
            dc.SetPen(wx.Pen(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNSHADOW)))
            dc.SetBrush(wx.Brush(face))
            dc.DrawRectangle(rect)

            if rect.width > 3 and rect.height > 3:
                hi = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNHIGHLIGHT)
                shadow = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNSHADOW)
                dark = wx.SystemSettings.GetColour(wx.SYS_COLOUR_3DDKSHADOW)
                # Recessed top/left, highlighted bottom/right.
                dc.SetPen(wx.Pen(dark))
                dc.DrawLine(0, 0, rect.width - 1, 0)
                dc.DrawLine(0, 0, 0, rect.height - 1)
                dc.SetPen(wx.Pen(shadow))
                dc.DrawLine(1, 1, rect.width - 2, 1)
                dc.DrawLine(1, 1, 1, rect.height - 2)
                dc.SetPen(wx.Pen(hi))
                dc.DrawLine(rect.width - 2, 1, rect.width - 2, rect.height - 2)
                dc.DrawLine(1, rect.height - 2, rect.width - 2, rect.height - 2)

            text = (wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT)
                    if not self.IsEnabled() else wx.WHITE)
            pressed = True

        dc.SetFont(self.GetFont())
        dc.SetTextForeground(text)
        tw, th = dc.GetTextExtent(self._label)
        x = max(0, (rect.width - tw) // 2)
        y = max(0, (rect.height - th) // 2)
        if pressed:
            x += 1
            y += 1
        dc.DrawText(self._label, x, y)
