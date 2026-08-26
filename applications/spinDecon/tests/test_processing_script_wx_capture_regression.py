from spinDecon.processing.script_context import ProcessingScriptState


class WxTextLike:
    def __init__(self, value):
        self._value = value
    def GetSelection(self):
        # wx.TextCtrl exposes this too: it is a character-range tuple.
        return (1, 3)
    def GetValue(self):
        return self._value


class WxComboLike(WxTextLike):
    def __init__(self, value, selection):
        super().__init__(value)
        self._selection = selection
    def GetSelection(self):
        return self._selection


class Processing3DLike:
    p0_1 = WxTextLike('37.5')
    p1_1 = WxTextLike('12.0')
    win2Val1 = WxTextLike('20')
    maxIterBox = WxTextLike('64')
    windowBox1 = WxComboLike('SP', 1)
    cb_ft1 = WxComboLike('Alt', 2)


def test_text_controls_capture_value_not_text_selection_tuple():
    state = ProcessingScriptState.capture_current(Processing3DLike())
    assert state.value('p0_1') == '37.5'
    assert state.value('p1_1') == '12.0'
    assert state.value('win2Val1') == '20'
    assert state.value('maxIterBox') == '64'
    assert float(state.value('p0_1')) == 37.5


def test_window_combo_captures_display_value_but_ft_combo_captures_index():
    state = ProcessingScriptState.capture_current(Processing3DLike())
    assert state.value('windowBox1') == 'SP'
    assert state.value('cb_ft1') == 2
