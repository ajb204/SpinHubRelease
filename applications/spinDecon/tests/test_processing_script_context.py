from types import MappingProxyType

from spinDecon.processing.script_context import ProcessingScriptState


class Control:
    def __init__(self, value): self.value = value
    def GetValue(self): return self.value


class Check:
    def __init__(self, value): self.value = value
    def IsChecked(self): return self.value


class Choice:
    def __init__(self, value): self.value = value
    def GetSelection(self): return self.value


class Processing:
    def __init__(self):
        self.p0_1 = Control('51')
        self.cb_lp1 = Check(True)
        self.cb_ft1 = Choice(3)


def snapshot(processing=None, live=None, names=('p0_1', 'cb_lp1', 'cb_ft1')):
    return ProcessingScriptState.capture(processing, live or {}, names)


def test_state_freezes_current_widget_value_at_capture():
    processing = Processing()
    state = snapshot(processing, {'p0_1': '37'})
    processing.p0_1.value = '99'
    assert state.value('p0_1') == '51'


def test_state_uses_shared_live_state_when_widget_is_absent():
    state = snapshot(None, {'p0_1': '37'})
    assert state.value('p0_1') == '37'


def test_state_contains_plain_values_not_control_adapters():
    state = snapshot(Processing())
    assert state.value('p0_1') == '51'
    assert state.checked('cb_lp1') is True
    assert state.selection('cb_ft1') == 3
    assert not hasattr(state.control('p0_1'), 'GetValue')
    assert not hasattr(state.control('cb_lp1'), 'IsChecked')


def test_state_mapping_is_immutable():
    state = snapshot(Processing())
    assert isinstance(state.values, MappingProxyType)
    try:
        state.values['p0_1'] = '0'
    except TypeError:
        pass
    else:
        raise AssertionError('script state unexpectedly allowed mutation')
